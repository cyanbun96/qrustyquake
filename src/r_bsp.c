// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

static mvertex_t pbverts[MAX_BMODEL_VERTS]; // brush faces clipped to the BSP
static bedge_t pbedges[MAX_BMODEL_EDGES];
static s32 numbverts, numbedges;
static f32 entity_rotation[3][3];
// The world faces the main pass left for the alpha pass, in BSP order
typedef struct { msurface_t *surf; s32 key, clipflags; } alphaface_t;
static alphaface_t *r_alphafaces;
static s32 r_numalphafaces;
static bool r_alphafaces_lost; // ran out of memory, traverse instead

void R_RotateBy(f32 m[3][3], vec3_t vec)
{
	vec3_t tvec;
	VectorCopy(vec, tvec);
	vec[0] = DotProduct(m[0], tvec);
	vec[1] = DotProduct(m[1], tvec);
	vec[2] = DotProduct(m[2], tvec);
}

void R_EntityRotate(vec3_t vec)
{ R_RotateBy(entity_rotation, vec); }

// The turn of e's model, which R_RotateBy applies
void R_EntityRotation(const entity_t *e, f32 m[3][3])
{
	f32 s[3], c[3];
	for(s32 i = 0; i < 3; i++){
		f32 angle = e->angles[i] * M_PI * 2 / 360;
		s[i] = sin(angle);
		c[i] = cos(angle);
	}
	f32 yaw[3][3] = {{c[YAW], s[YAW], 0}, {-s[YAW], c[YAW], 0}, {0, 0, 1}};
	f32 pitch[3][3] = {{c[PITCH], 0, -s[PITCH]}, {0, 1, 0},
		{s[PITCH], 0, c[PITCH]}};
	f32 roll[3][3] = {{1, 0, 0}, {0, c[ROLL], s[ROLL]},
		{0, -s[ROLL], c[ROLL]}};
	f32 temp[3][3];
	R_ConcatRotations(pitch, yaw, temp);
	R_ConcatRotations(roll, temp, m);
}

void R_RotateBmodel()
{
	R_EntityRotation(currententity, entity_rotation);
	R_EntityRotate(modelorg); // rotate modelorg and transformation matrix
	R_EntityRotate(vpn);
	R_EntityRotate(vright);
	R_EntityRotate(vup);
	R_TransformFrustum();
}

static void R_BEdgeBounds(const bedge_t *e, f32 *bounds)
{ // of the vertices of a list of edges, mins then maxs
	for(s32 k = 0; k < 3; k++){
		bounds[k] = FLT_MAX;
		bounds[3 + k] = -FLT_MAX;
	}
	for(; e; e = e->pnext)
	for(s32 j = 0; j < 2; j++)
	for(s32 k = 0; k < 3; k++){
		f32 x = e->v[j]->position[k];
		if(x < bounds[k]) bounds[k] = x;
		if(x > bounds[3 + k]) bounds[3 + k] = x;
	}
}

static bool R_BEdgesInPlane(const bedge_t *e, const mplane_t *p)
{ // all the vertices of a closed list of edges, within float error
	for(; e; e = e->pnext){
		f32 d = DotProduct(e->v[0]->position, p->normal) - p->dist;
		if(d > 0.01f || d < -0.01f)
			return 0;
	}
	return 1;
}

// Puts a new edge from v0 to v1 on a side list
static void R_NewBEdge(bedge_t **side, mvertex_t *v0, mvertex_t *v1)
{
	bedge_t *e = &pbedges[numbedges++];
	e->pnext = *side;
	*side = e;
	e->v[0] = v0;
	e->v[1] = v1;
}

// Closes both sides along the clip plane like R_RecursiveClipBPoly does for a
// convex face, which crosses it twice. Float error can zigzag an edge that
// runs along the plane across it more often, and joining only two of the
// crossings left both sides open: a face with no edge to close it on a row got
// a span to the right of the screen in the alpha pass. The crossings are the
// last n clipped vertices. In order along the line the face meets the plane
// on, each joins the nearest one of the other kind.
static Q_NOINLINE bool R_JoinBCrossings(bedge_t **sides, s32 n,
		mnode_t *pnode, msurface_t *psurf)
{
	typedef struct { mvertex_t *v; f32 along; bool enter; } bcross_t;
	bcross_t c[64];
	if(n > (s32)Q_COUNTOF(c))
		return 0;
	vec3_t normal, dir; // the node's plane in model space, as clipped
	const f32 *pn = pnode->plane->normal;
	for(s32 k = 0; k < 3; k++)
		normal[k] = DotProduct(entity_rotation[k], pn);
	CrossProduct(normal, psurf->plane->normal, dir);
	for(s32 i = 0; i < n; i++){
		bcross_t x = { .v = &pbverts[numbverts - n + i] };
		x.along = DotProduct(x.v->position, dir);
		// the front side goes on from where it comes in
		bedge_t *e = sides[0];
		while(e && e->v[0] != x.v)
			e = e->pnext;
		x.enter = e != NULL;
		s32 j = i;
		for(; j > 0 && c[j - 1].along > x.along; j--)
			c[j] = c[j - 1];
		c[j] = x;
	}
	s32 open = 0; // c[0..open) wait for one of the other kind
	for(s32 i = 0; i < n; i++){
		bcross_t x = c[i];
		if(!open || c[open - 1].enter == x.enter){
			c[open++] = x;
			continue;
		}
		mvertex_t *o = c[--open].v;
		mvertex_t *in = x.enter ? x.v : o, *out = x.enter ? o : x.v;
		R_NewBEdge(&sides[0], out, in);
		R_NewBEdge(&sides[1], in, out);
	}
	return 1;
}

static void R_RecursiveClipBPoly(bedge_t *pedges, mnode_t *pnode,
		msurface_t *psurf, const f32 *bounds)
{
	bedge_t *psideedges[2] = {NULL, NULL}, *pnextedge;
	bool split = 0; // anything clipped
	// transform the BSP plane into model space
	mplane_t *splitplane = pnode->plane;
	mplane_t tplane;
	tplane.dist=splitplane->dist-DotProduct(r_entorigin,splitplane->normal);
	tplane.normal[0] = DotProduct(entity_rotation[0], splitplane->normal);
	tplane.normal[1] = DotProduct(entity_rotation[1], splitplane->normal);
	tplane.normal[2] = DotProduct(entity_rotation[2], splitplane->normal);
	// Most of the time the plane misses the fragment's bounds, and the
	// loop below would only put all its edges on one side, in reverse.
	// The margin is for the rounding in the vertices' own distances.
	f32 mid = -tplane.dist, half = 0;
	for(s32 k = 0; k < 3; k++){
		f32 center = (bounds[k] + bounds[3 + k]) * 0.5f;
		f32 extent = (bounds[3 + k] - bounds[k]) * 0.5f;
		mid += tplane.normal[k] * center;
		half += fabsf(tplane.normal[k]) * extent;
	}
	s32 whole = mid - half > 1 ? 0 : mid + half < -1 ? 1 : -1;
	// The test below puts a fragment that lies in the plane on the back
	// side, solid where the world brush there has skip faces. It can be
	// seen from where it faces, so it goes behind the node's own faces
	// and in front of everything on the far side, like one of them. Not
	// in the alpha pass: that tests depth against the opaque world faces
	// instead of keys, and would draw it over one in the same plane.
	if(!r_alphapass && whole < 0 && R_BEdgesInPlane(pedges, &tplane)){
		r_currentbkey = pnode->key;
		R_RenderBmodelFace(pedges, psurf);
		return;
	}
	for(; whole >= 0 && pedges; pedges = pnextedge){
		pnextedge = pedges->pnext;
		pedges->pnext = psideedges[whole];
		psideedges[whole] = pedges;
	}
	static mvertex_t *pfrontenter, *pfrontexit;
	static s32 ncross; // more than two only through float error
	ncross = 0;
	for(; pedges; pedges = pnextedge){ // clip edges to BSP plane
		pnextedge = pedges->pnext;
		// set the status for the last point as the previous point
		mvertex_t *plastvert = pedges->v[0];
		f32 lastdist = DotProduct(plastvert->position, tplane.normal)
			- tplane.dist;
// CyanBun96: comparison of lastside and side against a low threshold instead of
// 0 prevents corruption in some edge cases(no pun intended)
		s32 lastside = (lastdist > 0.001) ? 0 : 1;
		mvertex_t *pvert = pedges->v[1];
		f32 dist = DotProduct(pvert->position, tplane.normal)
			- tplane.dist;
		s32 side = (dist > 0.001) ? 0 : 1;
		if(side != lastside){
			if(numbverts >= MAX_BMODEL_VERTS) // clipped
				return;
			// generate the clipped vertex. With both in front, the
			// threshold put one behind, and it's at that one
			f32 frac = lastdist / (lastdist - dist);
			if(lastdist > 0 && dist > 0)
				frac = side;
			mvertex_t *ptvert = &pbverts[numbverts++];
			for(s32 k = 0; k < 3; k++)
				ptvert->position[k] = plastvert->position[k] +
					frac * (pvert->position[k] -
					plastvert->position[k]);
			// split into two edges, one on each side, and remember
			// entering and exiting points
			if(numbedges >= (MAX_BMODEL_EDGES - 1)){
				Con_Printf("Out of edges for bmodel\n");
				return;
			}
			R_NewBEdge(&psideedges[lastside], plastvert, ptvert);
			R_NewBEdge(&psideedges[side], ptvert, pvert);
			if(side == 0) // entering for front, exiting for back
				pfrontenter = ptvert;
			else
				pfrontexit = ptvert;
			ncross++;
			split = 1;
		}else{ // add the edge to the appropriate side
			pedges->pnext = psideedges[side];
			psideedges[side] = pedges;
		}
	}
	// if anything was clipped, reconstitute and add the edges along
	// the clip plane to both sides(but in opposite directions)
	if(split){
		if(numbedges >= MAX_BMODEL_EDGES - ncross){
			Con_Printf("Out of edges for bmodel\n");
			return;
		}
		if(ncross > 2){
			if(!R_JoinBCrossings(psideedges, ncross, pnode, psurf))
				return;
		}else{
			R_NewBEdge(&psideedges[0], pfrontexit, pfrontenter);
			R_NewBEdge(&psideedges[1], pfrontenter, pfrontexit);
		}
	}
	for(s32 i = 0; i < 2; i++){ // draw or recurse further
		if(!psideedges[i])
			continue;
		// draw if we reached a non-solid leaf, done if all that's left
		// is a solid leaf and continue down the tree if it's not a leaf
		mnode_t *pn = pnode->children[i];
		// done with this branch if the world traversal didn't get to
		// the node or leaf this pass: it's not in the PVS, off screen
		// or hidden
		if(pn->keypass != r_worldpass)
			continue;
		if(pn->contents < 0){
			if(pn->contents != CONTENTS_SOLID){
				r_currentbkey = ((mleaf_t *) pn)->key;
				R_RenderBmodelFace(psideedges[i], psurf);
			}
		}else if(!split) // the same edges, the same bounds
			R_RecursiveClipBPoly(psideedges[i], pn, psurf, bounds);
		else{
			f32 sidebounds[6];
			R_BEdgeBounds(psideedges[i], sidebounds);
			R_RecursiveClipBPoly(psideedges[i], pn, psurf,
					sidebounds);
		}
	}
}

static bool R_BmodelFaceFront(msurface_t *psurf)
{
	mplane_t *pplane = psurf->plane; // find which node side we are on
	vec_t dot = DotProduct(modelorg, pplane->normal) - pplane->dist;
	return ((psurf->flags&SURF_PLANEBACK) && (dot < -BACKFACE_EPSILON)) ||
		(!(psurf->flags&SURF_PLANEBACK) && (dot > BACKFACE_EPSILON));
}

// The world node R_ClipBmodelFace starts the faces of a brush entity from
static mnode_t *r_bclipnode;

// The deepest node under a world node that a world space box doesn't straddle,
// with the margin of R_RecursiveClipBPoly's whole fragment test. NULL if that
// is solid or the world traversal didn't go into it, then none of what is in
// the box gets drawn.
static mnode_t *R_BoxWorldNode(mnode_t *node, const f32 *mn, const f32 *mx)
{
	while(node->keypass == r_worldpass){
		if(node->contents < 0)
			return node->contents != CONTENTS_SOLID ? node : NULL;
		const mplane_t *p = node->plane;
		f32 lo = -p->dist, hi = -p->dist;
		for(s32 k = 0; k < 3; k++){
			f32 a = p->normal[k] * mn[k], b = p->normal[k] * mx[k];
			lo += a < b ? a : b;
			hi += a < b ? b : a;
		}
		if(lo > 1)
			node = node->children[0];
		else if(hi < -1)
			node = node->children[1];
		else
			return node;
	}
	return NULL;
}

static void R_ClipBmodelFace(model_t *pmodel, msurface_t *psurf,
		s32 clipflags)
{
	if(!R_BmodelFaceFront(psurf))
		return;
	// would be thrown away after clipping it against the BSP anyway
	if(!R_SurfaceInPass(psurf))
		return;
	// hidden behind world faces: no need to clip it against the BSP
	if(r_bmodel_occ_faces && R_CoarseOcclusionTestBmodelFace(psurf)){
		r_bmodel_faces_occluded++;
		return;
	}
	// Cutouts are drawn as whole spans, z-tested against everything
	// opaque, so the pieces that clipping would throw away are hidden all
	// the same, and the pieces don't need the keys of their leaves. The
	// alpha pass sorts translucent faces on those keys though.
	if((psurf->flags & SURF_DRAWCUTOUT) && !r_alphapass){
		R_RenderFace(psurf, clipflags);
		return;
	}
	if(psurf->numedges <= 0)
		Sys_Error("no edges in bmodel");
	numbverts = 0;
	numbedges = psurf->numedges;
	bedge_t *pbedge = pbedges;
	for(s32 j = 0; j < psurf->numedges; j++){
		s32 lindex = pmodel->surfedges[psurf->firstedge + j];
		bool rev = lindex <= 0;
		medge_t *pedge = &pmodel->edges[rev ? -lindex : lindex];
		pbedge[j].v[0] = &r_pcurrentvertbase[pedge->v[rev]];
		pbedge[j].v[1] = &r_pcurrentvertbase[pedge->v[!rev]];
		pbedge[j].pnext = &pbedge[j + 1];
	}
	pbedge[numbedges - 1].pnext = NULL; // mark end of edges
	surf_t *firstsurf = surface_p;
	if(r_bclipnode->contents < 0){ // all of it in a leaf that's drawn
		r_currentbkey = ((mleaf_t *)r_bclipnode)->key;
		R_RenderBmodelFace(pbedge, psurf);
	}else{
		f32 bounds[6];
		R_BEdgeBounds(pbedge, bounds);
		R_RecursiveClipBPoly(pbedge, r_bclipnode, psurf, bounds);
	}
	// faces that come later can be hidden behind it, wherever it got
	// clipped away it's behind world faces anyway
	if(r_bmodel_occ_faces && surface_p > firstsurf)
		R_CoarseOcclusionAddSurface(psurf, surface_p - 1);
}

static void R_UnclippedBmodelFace(SDL_UNUSED model_t *pmodel,
		msurface_t *psurf, s32 clipflags)
{
	if(!R_BmodelFaceFront(psurf))
		return;
	r_currentkey = ((mleaf_t *) currententity->topnode)->key;
	if(r_bmodel_occ_faces && R_CoarseOcclusionTestBmodelFace(psurf)){
		r_bmodel_faces_occluded++;
		return;
	}
	R_RenderFace(psurf, clipflags);
}

typedef void (*bmodelface_f)(model_t *pmodel, msurface_t *psurf,
		s32 clipflags);

// Goes through the faces on the brush model's own BSP tree, skipping the nodes
// whose faces are all off screen or, if the entity isn't rotated, hidden
// behind world faces or only in parts of the world that the world traversal
// didn't go into. The clip planes are in the entity's space, like the
// facebounds. wnode is a world node that has all of them, for faces that get
// clipped against the world BSP if the entity isn't rotated, else NULL.
static void R_BmodelFaceNode(model_t *pmodel, mnode_t *node, s32 clipflags,
		bool occlusion, mnode_t *wnode, bmodelface_f face)
{
	const f32 *b = pmodel->facebounds[node - pmodel->nodes];
	for(s32 i = 0; i < 4; i++){
		if(!(clipflags & (1 << i)))
			continue;
		f32 *n = view_clipplanes[i].normal;
		f32 d = -view_clipplanes[i].dist;
		for(s32 k = 0; k < 3; k++)
			d += n[k] * (n[k] < 0 ? b[k] : b[3 + k]);
		if(d <= 0)
			return; // all on the outside of this plane
	}
	if(occlusion || wnode){
		vec3_t mn, mx;
		for(s32 k = 0; k < 3; k++){
			mn[k] = b[k] + currententity->origin[k];
			mx[k] = b[3 + k] + currententity->origin[k];
		}
		if(occlusion && R_CoarseOcclusionTestBox(mn, mx))
			return;
		// most faces of big brush entities get clipped into leaves
		// that aren't drawn, this finds that before looking at them
		if(wnode && !(wnode = R_BoxWorldNode(wnode, mn, mx)))
			return;
	}
	// the near side first, its faces can hide the ones behind them
	const mplane_t *p = node->plane;
	s32 side = DotProduct(modelorg, p->normal) - p->dist < 0;
	mnode_t *front = node->children[side], *back = node->children[!side];
	if(front->contents >= 0)
		R_BmodelFaceNode(pmodel, front, clipflags, occlusion, wnode,
				face);
	r_bclipnode = wnode ? wnode : currententity->topnode;
	msurface_t *psurf = pmodel->surfaces + node->firstsurface;
	for(u32 i = 0; i < node->numsurfaces; i++, psurf++)
		face(pmodel, psurf, clipflags);
	if(back->contents >= 0)
		R_BmodelFaceNode(pmodel, back, clipflags, occlusion, wnode,
				face);
}

static void R_BmodelFaces(model_t *pmodel, s32 clipflags, bmodelface_f face)
{
	r_bclipnode = currententity->topnode;
	if(pmodel->facetree){
		bool fixed = !currententity->angles[0] &&
			!currententity->angles[1] && !currententity->angles[2];
		mnode_t *wnode = fixed && face == R_ClipBmodelFace ?
			r_bclipnode : NULL;
		R_BmodelFaceNode(pmodel, pmodel->nodes +
				pmodel->hulls[0].firstclipnode, clipflags,
				r_bmodel_occ_faces && fixed, wnode, face);
		return;
	}
	msurface_t *psurf = &pmodel->surfaces[pmodel->firstmodelsurface];
	for(s32 i = 0; i < pmodel->nummodelsurfaces; i++, psurf++)
		face(pmodel, psurf, clipflags);
}

void R_DrawSolidClippedSubmodelPolygons(model_t *pmodel)
{ R_BmodelFaces(pmodel, r_clipflags, R_ClipBmodelFace); }

void R_DrawSubmodelPolygons(model_t *pmodel, s32 clipflags)
{ R_BmodelFaces(pmodel, clipflags, R_UnclippedBmodelFace); }

// The clip planes of clipflags that the node still crosses, -1 if it's off
// screen
Q_SIL s32 R_NodeClipFlags(const mnode_t *node, s32 clipflags)
{
	for(s32 i = 0; i < 4; i++){
		if(!(clipflags & (1 << i)))
			continue; // don't need to clip against it
		// generate accept and reject points
		s32 *pindex = frustum_indexes[i];
		vec3_t acceptpt, rejectpt;
		rejectpt[0] = (f32)node->minmaxs[pindex[0]];
		rejectpt[1] = (f32)node->minmaxs[pindex[1]];
		rejectpt[2] = (f32)node->minmaxs[pindex[2]];
		f64 d = DotProduct(rejectpt, view_clipplanes[i].normal);
		d -= view_clipplanes[i].dist;
		if(d <= 0)
			return -1;
		acceptpt[0] = (f32)node->minmaxs[pindex[3 + 0]];
		acceptpt[1] = (f32)node->minmaxs[pindex[3 + 1]];
		acceptpt[2] = (f32)node->minmaxs[pindex[3 + 2]];
		d = DotProduct(acceptpt, view_clipplanes[i].normal);
		d -= view_clipplanes[i].dist;
		if(d >= 0) // node is entirely on screen
			clipflags &= ~(1 << i);
	}
	return clipflags;
}

// The side of the node's plane the view is on, as a distance
Q_SIL f64 R_NodeDot(const mnode_t *node)
{
	mplane_t *plane = node->plane;
	switch(plane->type){
	case PLANE_X: return modelorg[0] - plane->dist;
	case PLANE_Y: return modelorg[1] - plane->dist;
	case PLANE_Z: return modelorg[2] - plane->dist;
	default: return DotProduct(modelorg, plane->normal) - plane->dist;
	}
}

// The run-ahead job (R_RunAheadNode) skips the nodes the walk found hidden,
// this time if it's fallen behind, or the last time: they likely are again.
// What it doesn't build the walk builds itself.
static SDL_AtomicInt *r_prehidden; // by node, r_worldpass when found hidden
static s32 r_prehiddencap;
static bool r_runahead; // the walk has the job

Q_SIL void R_FoundHidden(mnode_t *node)
{
	if(r_runahead)
		SDL_SetAtomicInt(&r_prehidden[node - cl.worldmodel->nodes],
				r_worldpass);
}

void R_RecursiveWorldNode(mnode_t *node, s32 clipflags)
{
	// Solid leaves other than leaf 0 come from detail and fence brushes
	// (ericw-tools). They can be the only leaves listing the faces they
	// sit on, so they still mark them, as QuakeSpasm does.
	if(node->visframe != r_visframecount)
		return;
	// cull the clipping planes if not trivial accept
	if(clipflags && (clipflags = R_NodeClipFlags(node, clipflags)) < 0)
		return;
	// Hidden behind the world faces drawn so far? Not asked of leaves: the
	// faces they mark are tested one by one anyway, and on big maps the
	// box tests cost more than the marking they save.
	if(node->contents >= 0 && R_CoarseOcclusionTestNode(node)){
		r_nodes_occluded++;
		R_FoundHidden(node);
		return;
	}
	node->keypass = r_worldpass;
	if(node->contents < 0){ // if a leaf node, draw stuff
		mleaf_t *pleaf = (mleaf_t *) node;
		msurface_t **mark = pleaf->firstmarksurface;
		s32 c = pleaf->nummarksurfaces;
		if(c){
			do{
				(*mark)->visframe = r_framecount;
				mark++;
			}while(--c);
		}
		if(node->contents == CONTENTS_SOLID)
			return; // no efrags, and no bmodel is drawn in it
		if(pleaf->efrags) // deal with model fragments in this leaf
			R_StoreEfrags(&pleaf->efrags);
		pleaf->key = r_currentkey;
		r_currentkey++; // all bmodels in a leaf share the same key
	}else{
		// node is just a decision point so go down the apropriate sides
		// find which side of the node we are on
		f64 dot = R_NodeDot(node);
		s32 side = dot >= 0 ? 0 : 1;
		// recurse down the children, front side first
		R_RecursiveWorldNode(node->children[side], clipflags);
		// draw the faces that face us, the cutouts even where no leaf
		// marked them. Seen from behind, a face's edges come out the
		// wrong way round and the edge list draws none of it.
		s32 c = node->numsurfaces;
		if(c){
			msurface_t *surf = cl.worldmodel->surfaces
						+ node->firstsurface;
			if(dot < -BACKFACE_EPSILON || dot > BACKFACE_EPSILON){
				s32 back = dot < 0 ? SURF_PLANEBACK : 0;
				do{
					s32 f = surf->flags;
					if((f & SURF_PLANEBACK) == back &&
						((f & SURF_DRAWCUTOUT) ||
						surf->visframe == r_framecount))
						R_RenderFace(surf, clipflags);
					surf++;
				}while(--c);
			}
		}
		// behind its faces, for brush entity faces lying in its plane
		node->key = r_currentkey++;
		// recurse down the back side
		R_RecursiveWorldNode(node->children[!side], clipflags);
	}
}

// The run ahead: a worker walks the world tree in the same order while the
// world walk goes on, without its occlusion tests, and builds the coarse
// occlusion polygons of the faces it comes to (R_CoarsePreFace). It marks the
// faces of the leaves it goes into itself, they're a few more. Of what it
// reads the walk only changes r_prehidden.
static u32 *r_premark; // by face, the frame a leaf it went into listed it
static s32 r_premarkcap;
static SDL_AtomicInt r_prestop; // the walk is over

static void R_RunAheadNode(mnode_t *node, s32 clipflags)
{
	if(node->visframe != r_visframecount || SDL_GetAtomicInt(&r_prestop))
		return;
	if(clipflags && (clipflags = R_NodeClipFlags(node, clipflags)) < 0)
		return;
	if(node->contents >= 0){ // see R_FoundHidden
		s32 n = node - cl.worldmodel->nodes;
		s32 h = SDL_GetAtomicInt(&r_prehidden[n]);
		if(h && r_worldpass - h <= 1)
			return;
	}
	msurface_t *base = cl.worldmodel->surfaces;
	u32 frame = (u32)r_framecount;
	if(node->contents < 0){
		mleaf_t *pleaf = (mleaf_t *)node;
		msurface_t **mark = pleaf->firstmarksurface;
		for(s32 c = pleaf->nummarksurfaces; c; c--, mark++)
			r_premark[*mark - base] = frame;
		return;
	}
	f64 dot = R_NodeDot(node);
	s32 side = dot >= 0 ? 0 : 1;
	R_RunAheadNode(node->children[side], clipflags);
	s32 c = node->numsurfaces;
	if(c && (dot < -BACKFACE_EPSILON || dot > BACKFACE_EPSILON)){
		s32 back = dot < 0 ? SURF_PLANEBACK : 0;
		msurface_t *surf = base + node->firstsurface;
		for(; c; c--, surf++){
			s32 f = surf->flags;
			if((f & SURF_PLANEBACK) == back && !(f & SURF_DRAWSKY)
				&& ((f & SURF_DRAWCUTOUT) ||
				r_premark[surf - base] == frame))
				R_CoarsePreFace(surf);
		}
	}
	R_RunAheadNode(node->children[!side], clipflags);
}

static void R_RunAheadJob(SDL_UNUSED void *arg)
{ R_RunAheadNode(cl.worldmodel->nodes, 15); }

// p with room for n entries of size each, all 0 if it had to grow. NULL
// without the memory, and the next call tries again.
static void *R_ZeroedArray(void *p, s32 *cap, s32 n, size_t size)
{
	if(*cap >= n)
		return p;
	free(p);
	p = calloc(n, size);
	*cap = p ? n : 0;
	return p;
}

// Starts the run ahead of the world walk if it can, returns its job or -1
static s32 R_RunAheadBegin(const model_t *m)
{
	if(!r_threaded || !R_CoarsePreReady())
		return -1;
	r_premark = R_ZeroedArray(r_premark, &r_premarkcap, m->numsurfaces,
			sizeof(*r_premark));
	r_prehidden = R_ZeroedArray(r_prehidden, &r_prehiddencap,
			m->numnodes, sizeof(*r_prehidden));
	if(!r_premark || !r_prehidden)
		return -1;
	SDL_SetAtomicInt(&r_prestop, 0);
	s32 job = R_JobAdd(R_RunAheadJob, NULL);
	R_JobsPublish();
	R_CoarsePreOn();
	r_runahead = 1;
	return job;
}

static void R_RunAheadEnd(s32 job)
{
	SDL_SetAtomicInt(&r_prestop, 1);
	R_JobCancel(job);
	r_runahead = 0;
}

// Called by R_RenderFace in the main pass for world faces that belong in the
// alpha pass. Each gets a key of its own, to be drawn with later.
void R_RecordAlphaFace(msurface_t *surf, s32 clipflags)
{
	static s32 r_maxalphafaces;
	if(r_numalphafaces == r_maxalphafaces){
		s32 n = r_maxalphafaces ? r_maxalphafaces * 2 : 1024;
		alphaface_t *p = realloc(r_alphafaces, n * sizeof(*p));
		if(!p){
			r_alphafaces_lost = 1;
			return;
		}
		r_alphafaces = p;
		r_maxalphafaces = n;
	}
	alphaface_t *f = &r_alphafaces[r_numalphafaces++];
	f->surf = surf;
	f->key = r_currentkey++;
	f->clipflags = clipflags;
}

void R_RenderWorld()
{
	currententity = &cl_entities[0];
	VectorCopy(r_origin, modelorg);
	model_t *clmodel = currententity->model;
	r_pcurrentvertbase = clmodel->vertexes;
	if(r_alphapass && !r_alphafaces_lost){
		// Walking the tree again would only find the faces the main
		// pass recorded. The node and leaf keys it left stay valid for
		// the brush entities: whatever it skipped is hidden behind
		// opaque world faces.
		s32 endkey = r_currentkey;
		for(s32 i = 0; i < r_numalphafaces; i++){
			r_currentkey = r_alphafaces[i].key;
			R_RenderFace(r_alphafaces[i].surf,
					r_alphafaces[i].clipflags);
		}
		r_currentkey = endkey;
	}else{
		r_numalphafaces = 0;
		r_alphafaces_lost = 0;
		r_worldpass++;
		R_WalkBegin(); // a lone worker takes its emit job first
		s32 job = R_RunAheadBegin(clmodel);
		R_RecursiveWorldNode(clmodel->nodes, 15);
		R_RunAheadEnd(job);
		R_WalkEnd();
	}
}
