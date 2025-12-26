// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2007-2008 Kristian Duske
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
// world.c -- world query functions
#include "quakedef.h"

// entities never clip against themselves, or their owner
// line of sight checks trace->crosscontent, but bullets don't
static hull_t box_hull;
static mclipnode_t box_clipnodes[6]; //johnfitz -- was dclipnode_t
static mplane_t box_planes[6];
static areanode_t sv_areanodes[AREA_NODES];
static s32 sv_numareanodes;

s32 SV_HullPointContents(hull_t *hull, s32 num, vec3_t p);

//Set up the planes and clipnodes so that the six floats of a bounding box
//can just be stored out and get a proper hull_t structure.
void SV_InitBoxHull()
{
	box_hull.clipnodes = box_clipnodes;
	box_hull.planes = box_planes;
	box_hull.firstclipnode = 0;
	box_hull.lastclipnode = 5;
	for(s32 i = 0; i < 6; i++) {
		box_clipnodes[i].planenum = i;
		s32 side = i&1;
		box_clipnodes[i].children[side] = CONTENTS_EMPTY;
		if(i != 5) box_clipnodes[i].children[side^1] = i + 1;
		else box_clipnodes[i].children[side^1] = CONTENTS_SOLID;
		box_planes[i].type = i>>1;
		box_planes[i].normal[i>>1] = 1;
	}

}

// To keep everything totally uniform, bounding boxes are turned into small
// BSP trees instead of being compared directly.
hull_t *SV_HullForBox(vec3_t mins, vec3_t maxs)
{
	box_planes[0].dist = maxs[0];
	box_planes[1].dist = mins[0];
	box_planes[2].dist = maxs[1];
	box_planes[3].dist = mins[1];
	box_planes[4].dist = maxs[2];
	box_planes[5].dist = mins[2];
	return &box_hull;
}

// Returns a hull that can be used for testing or clipping an object of
// mins/maxs size.
// Offset is filled in to contain the adjustment that must be added to the
// testing object's origin to get a point to use with the returned hull.
hull_t *SV_HullForEntity(edict_t *ent, vec3_t mins, vec3_t maxs, vec3_t offset)
{
	model_t *model;
	vec3_t size;
	vec3_t hullmins, hullmaxs;
	hull_t *hull;
	// decide which clipping hull to use, based on the size
	if(ent->v.solid == SOLID_BSP) { // explicit hulls in the BSP model
		if(ent->v.movetype != MOVETYPE_PUSH)
		   Host_Error("SOLID_BSP without MOVETYPE_PUSH(%s at %f %f %f)",
			PR_GetString(ent->v.classname), ent->v.origin[0],
			ent->v.origin[1], ent->v.origin[2]);
		model = sv.models[(s32)ent->v.modelindex];
		if(!model || model->type != mod_brush)
		   Host_Error("SOLID_BSP with a non bsp model(%s at %f %f %f)",
			PR_GetString(ent->v.classname), ent->v.origin[0],
			ent->v.origin[1], ent->v.origin[2]);

		VectorSubtract(maxs, mins, size);
		if(size[0] < 3) hull = &model->hulls[0];
		else if(size[0] <= 32) hull = &model->hulls[1];
		else hull = &model->hulls[2];
		// calculate an offset value to center the origin
		VectorSubtract(hull->clip_mins, mins, offset);
		VectorAdd(offset, ent->v.origin, offset);
	} else { // create a temp hull from bounding box sizes
		VectorSubtract(ent->v.mins, maxs, hullmins);
		VectorSubtract(ent->v.maxs, mins, hullmaxs);
		hull = SV_HullForBox(hullmins, hullmaxs);
		VectorCopy(ent->v.origin, offset);
	}
	return hull;
}

areanode_t *SV_CreateAreaNode(s32 depth, vec3_t mins, vec3_t maxs)
{
	vec3_t size;
	vec3_t mins1, maxs1, mins2, maxs2;
	areanode_t *anode = &sv_areanodes[sv_numareanodes];
	sv_numareanodes++;
	ClearLink(&anode->trigger_edicts);
	ClearLink(&anode->solid_edicts);
	if(depth == AREA_DEPTH) {
		anode->axis = -1;
		anode->children[0] = anode->children[1] = NULL;
		return anode;
	}
	VectorSubtract(maxs, mins, size);
	if(size[0] > size[1]) anode->axis = 0;
	else anode->axis = 1;
	anode->dist = 0.5 * (maxs[anode->axis] + mins[anode->axis]);
	VectorCopy(mins, mins1);
	VectorCopy(mins, mins2);
	VectorCopy(maxs, maxs1);
	VectorCopy(maxs, maxs2);
	maxs1[anode->axis] = mins2[anode->axis] = anode->dist;
	anode->children[0] = SV_CreateAreaNode(depth+1, mins2, maxs2);
	anode->children[1] = SV_CreateAreaNode(depth+1, mins1, maxs1);
	return anode;
}

void SV_ClearWorld()
{
	SV_InitBoxHull();
	memset(sv_areanodes, 0, sizeof(sv_areanodes));
	sv_numareanodes = 0;
	SV_CreateAreaNode(0, sv.worldmodel->mins, sv.worldmodel->maxs);
}

void SV_UnlinkEdict(edict_t *ent)
{
	if(!ent->area.prev) return; // not linked in anywhere
	RemoveLink(&ent->area);
	ent->area.prev = ent->area.next = NULL;
}

// Spike -- just builds a list of entities within the area, rather than walking
// them and risking the list getting corrupt.
static void SV_AreaTriggerEdicts(edict_t *ent, areanode_t *node, edict_t **list,
		s32 *listcount, const s32 listspace )
{
	link_t *next;
	s32 i = 0;
	// touch linked edicts
	for(link_t*l=node->trigger_edicts.next;l!=&node->trigger_edicts;l=next){ 
		next = l->next;
		edict_t *touch = EDICT_FROM_AREA(l);
		if(touch == ent) continue;
		if(!touch->v.touch || touch->v.solid != SOLID_TRIGGER) continue;
		if(ent->v.absmin[0] > touch->v.absmax[0]
		|| ent->v.absmin[1] > touch->v.absmax[1]
		|| ent->v.absmin[2] > touch->v.absmax[2]
		|| ent->v.absmax[0] < touch->v.absmin[0]
		|| ent->v.absmax[1] < touch->v.absmin[1]
		|| ent->v.absmax[2] < touch->v.absmin[2] ) continue;
		if(*listcount == listspace) return; // should never happen
		list[*listcount] = touch;
		(*listcount)++;
		++i;
	}
	if(node->axis == -1) return; // recurse down both sides
	if(ent->v.absmax[node->axis] > node->dist)
	   SV_AreaTriggerEdicts(ent,node->children[0],list,listcount,listspace);
	if(ent->v.absmin[node->axis] < node->dist)
	   SV_AreaTriggerEdicts(ent,node->children[1],list,listcount,listspace);
}

// ericw -- copy the touching edicts to an array(on the hunk) so we can avoid
// iteating the trigger_edicts linked list while calling PR_ExecuteProgram
// which could potentially corrupt the list while it's being iterated.
// Based on code from Spike.
void SV_TouchLinks(edict_t *ent)
{
	s32 old_self, old_other;
	s32 mark = Hunk_LowMark();
	edict_t **list = (edict_t **)Hunk_Alloc(sv.num_edicts*sizeof(edict_t*));
	s32 listcount = 0;
	SV_AreaTriggerEdicts(ent, sv_areanodes, list, &listcount,sv.num_edicts);
	for(s32 i = 0; i < listcount; i++) {
		edict_t *touch = list[i];
		// re-validate in case of PR_ExecuteProgram having side effects
		// that make edicts later in the list no longer touch
		if(touch == ent) continue;
		if(!touch->v.touch || touch->v.solid != SOLID_TRIGGER) continue;
		if(ent->v.absmin[0] > touch->v.absmax[0]
		|| ent->v.absmin[1] > touch->v.absmax[1]
		|| ent->v.absmin[2] > touch->v.absmax[2]
		|| ent->v.absmax[0] < touch->v.absmin[0]
		|| ent->v.absmax[1] < touch->v.absmin[1]
		|| ent->v.absmax[2] < touch->v.absmin[2] ) continue;
		old_self = pr_global_struct->self;
		old_other = pr_global_struct->other;
		pr_global_struct->self = EDICT_TO_PROG(touch);
		pr_global_struct->other = EDICT_TO_PROG(ent);
		pr_global_struct->time = sv.time;
		PR_ExecuteProgram(touch->v.touch);
		pr_global_struct->self = old_self;
		pr_global_struct->other = old_other;
	}
	Hunk_FreeToLowMark(mark); // free hunk-allocated edicts array
}

void SV_FindTouchedLeafs(edict_t *ent, mnode_t *node)
{
	if(node->contents == CONTENTS_SOLID) return;
	if(ent->num_leafs == MAX_ENT_LEAFS) return;
	if(node->contents < 0) { // add an efrag if the node is a leaf
		mleaf_t *leaf = (mleaf_t *)node;
		s32 leafnum = leaf - sv.worldmodel->leafs - 1;
		ent->leafnums[ent->num_leafs] = leafnum;
		ent->num_leafs++;
		return;
	}
	mplane_t *splitplane = node->plane; // NODE_MIXED
	s32 sides = BoxOnPlaneSide(ent->v.absmin, ent->v.absmax, splitplane);
	// recurse down the contacted sides
	if(sides & 1) SV_FindTouchedLeafs(ent, node->children[0]);
	if(sides & 2) SV_FindTouchedLeafs(ent, node->children[1]);
}

void SV_LinkEdict(edict_t *ent, bool touch_triggers)
{
	if(ent->area.prev) SV_UnlinkEdict(ent); // unlink from old position
	if(ent == sv.edicts) return; // don't add the world
	if(ent->free) return;
	VectorAdd(ent->v.origin, ent->v.mins, ent->v.absmin); // set the abs box
	VectorAdd(ent->v.origin, ent->v.maxs, ent->v.absmax);

	if((s32)ent->v.flags & FL_ITEM) {
		ent->v.absmin[0] -= 15;	// to make items easier to pick up and
		ent->v.absmin[1] -= 15; // allow them to be grabbed off of
		ent->v.absmax[0] += 15; // shelves, the abs sizes are expanded
		ent->v.absmax[1] += 15;
	} else {
		ent->v.absmin[0] -= 1; // because movement is clipped an epsilon
		ent->v.absmin[1] -= 1; // away from an actual edge, we must
		ent->v.absmin[2] -= 1; // fully check even when bounding boxes
		ent->v.absmax[0] += 1; // don't quite touch
		ent->v.absmax[1] += 1;
		ent->v.absmax[2] += 1;
	}
	ent->num_leafs = 0; // link to PVS leafs
	if(ent->v.modelindex) SV_FindTouchedLeafs(ent, sv.worldmodel->nodes);
	if(ent->v.solid == SOLID_NOT) return;
	areanode_t *node = sv_areanodes;
	while(1) { // find the first node that the ent's box crosses
		if(node->axis == -1) break;
		if(ent->v.absmin[node->axis] > node->dist)
			node = node->children[0];
		else if(ent->v.absmax[node->axis] < node->dist)
			node = node->children[1];
		else break; // crosses the node
	}
	if(ent->v.solid == SOLID_TRIGGER) // link it in
		InsertLinkBefore(&ent->area, &node->trigger_edicts);
	else InsertLinkBefore(&ent->area, &node->solid_edicts);
	//if touch_triggers, touch all entities at this node and decend for more
	if(touch_triggers) SV_TouchLinks(ent);
}

s32 SV_HullPointContents(hull_t *hull, s32 num, vec3_t p)
{
	f32 d;
	while(num >= 0) {
		if(num < hull->firstclipnode || num > hull->lastclipnode)
			Sys_Error("SV_HullPointContents: bad node number");
		mclipnode_t *node = hull->clipnodes + num;
		mplane_t *plane = hull->planes + node->planenum;
		if(plane->type < 3) d = p[plane->type] - plane->dist;
		else d=DotProductF64(plane->normal, p)-plane->dist;
		if(d < 0) num = node->children[1];
		else num = node->children[0];
	}
	return num;
}

s32 SV_PointContents(vec3_t p)
{
	s32 cont = SV_HullPointContents(&sv.worldmodel->hulls[0], 0, p);
	if(cont <= CONTENTS_CURRENT_0 && cont >= CONTENTS_CURRENT_DOWN)
		cont = CONTENTS_WATER;
	return cont;
}

s32 SV_TruePointContents(vec3_t p)
{ return SV_HullPointContents(&sv.worldmodel->hulls[0], 0, p); }

edict_t *SV_TestEntityPosition(edict_t *ent)
{
	trace_t trace = SV_Move(ent->v.origin, ent->v.mins, ent->v.maxs,
			ent->v.origin, 0, ent);
	if(trace.startsolid) return sv.edicts;
	return NULL;
}

bool SV_RecursiveHullCheck(hull_t *hull, s32 num, f32 p1f, f32 p2f,
				vec3_t p1, vec3_t p2, trace_t *trace) {
	if(num < 0) { // check for empty
		if(num != CONTENTS_SOLID) { trace->allsolid = 0;
					if(num==CONTENTS_EMPTY)trace->inopen=1;
					else trace->inwater = 1;
		} else trace->startsolid = 1; return 1;} // empty
	if(num < hull->firstclipnode || num > hull->lastclipnode)
		Sys_Error("SV_RecursiveHullCheck: bad node number");
	mclipnode_t *node = hull->clipnodes + num; // find the point distances
	mplane_t *plane = hull->planes + node->planenum;
	f32 t1 = plane->type < 3 ? p1[plane->type] - plane->dist :
		DotProductF64(plane->normal, p1)-plane->dist;
	f32 t2 = plane->type < 3 ? p2[plane->type] - plane->dist :
		DotProductF64(plane->normal, p2)-plane->dist;
	if(t1 >= 0 && t2 >= 0) return SV_RecursiveHullCheck(
			hull, node->children[0], p1f, p2f, p1, p2, trace);
	if(t1 < 0 && t2 < 0) return SV_RecursiveHullCheck(
			hull, node->children[1], p1f, p2f, p1, p2, trace);
	// put the crosspoint DIST_EPSILON pixels on the near side
	f32 frac = t1<0 ? (t1+DIST_EPSILON)/(t1-t2) : (t1-DIST_EPSILON)/(t1-t2);
	frac = CLAMP(0, frac, 1);
	f32 midf = p1f + (p2f - p1f)*frac;
	vec3_t mid;
	for(s32 i = 0; i < 3; i++) mid[i] = p1[i] + frac*(p2[i] - p1[i]);
	s32 side = (t1 < 0);
	if(!SV_RecursiveHullCheck( // move up to the node
		hull, node->children[side], p1f, midf, p1, mid, trace))return 0;
	if(SV_HullPointContents(hull, node->children[side^1], mid) != -2)
		return SV_RecursiveHullCheck(hull, node->children[side^1],
				midf, p2f, mid, p2, trace); // go past the node
	if(trace->allsolid) return 0; // never got out of the solid area
	// the other side of the node is solid, this is the impact point
	if(!side) {  VectorCopy(plane->normal, trace->plane.normal);
		     trace->plane.dist = plane->dist;
	} else { VectorSubtract(vec3_origin, plane->normal,trace->plane.normal);
		trace->plane.dist = -plane->dist; }
	while(SV_HullPointContents(hull, hull->firstclipnode, mid) == -2) {
		frac -= 0.1; // shouldn't really happen, but does occasionally
		if(frac < 0) { trace->fraction = midf;
				VectorCopy(mid, trace->endpos);
				Con_DPrintf("backup past 0\n");
				return 0; }
		midf = p1f + (p2f - p1f)*frac;
		for(s32 i = 0; i < 3; i++)mid[i]=p1[i] + frac*(p2[i] - p1[i]); }
	trace->fraction = midf;
	VectorCopy(mid, trace->endpos);
	return 0;
}


// Handles selection or creation of a clipping hull, and offseting (and
// eventually rotation) of the end points
trace_t SV_ClipMoveToEntity(edict_t *ent, vec3_t start, vec3_t mins, 
		vec3_t maxs, vec3_t end) {
	vec3_t offset;
	vec3_t start_l, end_l;
	trace_t trace;
	memset(&trace, 0, sizeof(trace_t)); // fill in a default trace
	trace.fraction = 1;
	trace.allsolid = 1;
	VectorCopy(end, trace.endpos);
	hull_t *h=SV_HullForEntity(ent, mins, maxs, offset); //get clipping hull
	VectorSubtract(start, offset, start_l);
	VectorSubtract(end, offset, end_l);
	// trace a line through the apropriate clipping hull
	SV_RecursiveHullCheck(h, h->firstclipnode, 0, 1, start_l, end_l,&trace);
	if(trace.fraction != 1) // fix trace up by the offset
		VectorAdd(trace.endpos, offset, trace.endpos);
	if(trace.fraction < 1 || trace.startsolid ) // did we clip the move?
		trace.ent = ent;
	return trace;
}


void SV_ClipToLinks( areanode_t *node, moveclip_t *clip )
{ // Mins and maxs enclose the entire area swept by the move
	link_t *next; // touch linked edicts
	for(link_t *l=node->solid_edicts.next; l!=&node->solid_edicts; l=next){
		next = l->next;
		edict_t *touch = EDICT_FROM_AREA(l);
		if(touch->v.solid == SOLID_NOT) continue;
		if(touch == clip->passedict) continue;
		if(touch->v.solid == SOLID_TRIGGER)
			Sys_Error("Trigger in clipping list");
		if(clip->type == MOVE_NOMONSTERS && touch->v.solid != SOLID_BSP)
			continue;
		if(clip->boxmins[0] > touch->v.absmax[0]
		|| clip->boxmins[1] > touch->v.absmax[1]
		|| clip->boxmins[2] > touch->v.absmax[2]
		|| clip->boxmaxs[0] < touch->v.absmin[0]
		|| clip->boxmaxs[1] < touch->v.absmin[1]
		|| clip->boxmaxs[2] < touch->v.absmin[2] ) continue;
		if(clip->passedict&&clip->passedict->v.size[0]&&!touch->v.size[0])
			continue; // points never interact
		// might intersect, so do an exact clip
		if(clip->trace.allsolid) return;
		if(clip->passedict) {
			if(PROG_TO_EDICT(touch->v.owner) == clip->passedict)
				continue; // don't clip against own missiles
			if(PROG_TO_EDICT(clip->passedict->v.owner) == touch)
				continue; // don't clip against owner
		}
		trace_t trace = (s32)touch->v.flags & FL_MONSTER ?
       SV_ClipMoveToEntity(touch,clip->start,clip->mins2,clip->maxs2,clip->end):
       SV_ClipMoveToEntity(touch,clip->start,clip->mins,clip->maxs,clip->end);
		if(trace.allsolid || trace.startsolid ||
				trace.fraction < clip->trace.fraction) {
			trace.ent = touch;
			if(clip->trace.startsolid) {
				clip->trace = trace;
				clip->trace.startsolid = 1;
			} else clip->trace = trace;
		} else if(trace.startsolid) clip->trace.startsolid = 1;
	}
	if(node->axis == -1) return; // recurse down both sides
	if(clip->boxmaxs[node->axis] > node->dist)
		SV_ClipToLinks( node->children[0], clip );
	if(clip->boxmins[node->axis] < node->dist)
		SV_ClipToLinks( node->children[1], clip );
}

void SV_MoveBounds(vec3_t start, vec3_t mins, vec3_t maxs, vec3_t end,
		vec3_t boxmins, vec3_t boxmaxs) {
	for(s32 i = 0; i < 3; i++) {
		if(end[i] > start[i]) {
			boxmins[i] = start[i] + mins[i] - 1;
			boxmaxs[i] = end[i] + maxs[i] + 1;
		} else {
			boxmins[i] = end[i] + mins[i] - 1;
			boxmaxs[i] = start[i] + maxs[i] + 1;
		}
	}
}

trace_t SV_Move(vec3_t start, vec3_t mins, vec3_t maxs, vec3_t end,
		s32 type, edict_t *passedict) {
	moveclip_t clip;
	memset(&clip, 0, sizeof( moveclip_t ));
	// clip to world
	clip.trace = SV_ClipMoveToEntity(sv.edicts, start, mins, maxs, end);
	clip.start = start;
	clip.end = end;
	clip.mins = mins;
	clip.maxs = maxs;
	clip.type = type;
	clip.passedict = passedict;
	if(type == MOVE_MISSILE) {
		for(s32 i = 0; i < 3; i++) {
			clip.mins2[i] = -15;
			clip.maxs2[i] = 15;
		}
	} else {
		VectorCopy(mins, clip.mins2);
		VectorCopy(maxs, clip.maxs2);
	}
	// create the bounding box of the entire move
	SV_MoveBounds(start, clip.mins2, clip.maxs2, end,
			clip.boxmins, clip.boxmaxs);
	// clip to entities
	SV_ClipToLinks( sv_areanodes, &clip );
	return clip.trace;
}
