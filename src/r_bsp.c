// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

static mvertex_t *pbverts;
static bedge_t *pbedges;
static s32 numbverts, numbedges;
static f32 entity_rotation[3][3];
// The world faces the main pass left for the alpha pass, in BSP order
typedef struct { msurface_t *surf; s32 key, clipflags; } alphaface_t;
static alphaface_t *r_alphafaces;
static s32 r_numalphafaces;
static bool r_alphafaces_lost; // ran out of memory, traverse instead

void R_EntityRotate(vec3_t vec)
{
	vec3_t tvec;
	VectorCopy(vec, tvec);
	vec[0] = DotProduct(entity_rotation[0], tvec);
	vec[1] = DotProduct(entity_rotation[1], tvec);
	vec[2] = DotProduct(entity_rotation[2], tvec);
}

void R_RotateBmodel()
{
	f32 angle = currententity->angles[YAW];
	angle = angle * M_PI * 2 / 360;
	f32 s = sin(angle);
	f32 c = cos(angle);
	f32 temp1[3][3], temp2[3][3], temp3[3][3];
	temp1[0][0] = c; temp1[0][1] = s; temp1[0][2] = 0;
	temp1[1][0] = -s; temp1[1][1] = c; temp1[1][2] = 0;
	temp1[2][0] = 0; temp1[2][1] = 0; temp1[2][2] = 1;
	angle = currententity->angles[PITCH];
	angle = angle * M_PI * 2 / 360;
	s = sin(angle);
	c = cos(angle);
	temp2[0][0] = c; temp2[0][1] = 0; temp2[0][2] = -s;
	temp2[1][0] = 0; temp2[1][1] = 1; temp2[1][2] = 0;
	temp2[2][0] = s; temp2[2][1] = 0; temp2[2][2] = c;
	R_ConcatRotations(temp2, temp1, temp3);
	angle = currententity->angles[ROLL];
	angle = angle * M_PI * 2 / 360;
	s = sin(angle);
	c = cos(angle);
	temp1[0][0] = 1; temp1[0][1] = 0; temp1[0][2] = 0;
	temp1[1][0] = 0; temp1[1][1] = c; temp1[1][2] = s;
	temp1[2][0] = 0; temp1[2][1] = -s; temp1[2][2] = c;
	R_ConcatRotations(temp1, temp3, entity_rotation);
	R_EntityRotate(modelorg); // rotate modelorg and transformation matrix
	R_EntityRotate(vpn);
	R_EntityRotate(vright);
	R_EntityRotate(vup);
	R_TransformFrustum();
}

void R_RecursiveClipBPoly(bedge_t *pedges, mnode_t *pnode, msurface_t *psurf)
{
	bedge_t *psideedges[2], *pnextedge, *ptedge;
	psideedges[0] = psideedges[1] = NULL;
	static bool makeclippededge;
	makeclippededge = 0;
	// transform the BSP plane into model space
	mplane_t *splitplane = pnode->plane;
	mplane_t tplane;
	tplane.dist=splitplane->dist-DotProduct(r_entorigin,splitplane->normal);
	tplane.normal[0] = DotProduct(entity_rotation[0], splitplane->normal);
	tplane.normal[1] = DotProduct(entity_rotation[1], splitplane->normal);
	tplane.normal[2] = DotProduct(entity_rotation[2], splitplane->normal);
	static mvertex_t *pfrontenter, *pfrontexit;
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
			// generate the clipped vertex
			f32 frac = lastdist / (lastdist - dist);
			mvertex_t *ptvert = &pbverts[numbverts++];
			ptvert->position[0] = plastvert->position[0] + frac *
				(pvert->position[0] - plastvert->position[0]);
			ptvert->position[1] = plastvert->position[1] + frac *
				(pvert->position[1] - plastvert->position[1]);
			ptvert->position[2] = plastvert->position[2] + frac *
				(pvert->position[2] - plastvert->position[2]);
			// split into two edges, one on each side, and remember
			// entering and exiting points
			if(numbedges >= (MAX_BMODEL_EDGES - 1)){
				Con_Printf("Out of edges for bmodel\n");
				return;
			}
			ptedge = &pbedges[numbedges];
			ptedge->pnext = psideedges[lastside];
			psideedges[lastside] = ptedge;
			ptedge->v[0] = plastvert;
			ptedge->v[1] = ptvert;
			ptedge = &pbedges[numbedges + 1];
			ptedge->pnext = psideedges[side];
			psideedges[side] = ptedge;
			ptedge->v[0] = ptvert;
			ptedge->v[1] = pvert;
			numbedges += 2;
			if(side == 0){ // entering for front, exiting for back
				pfrontenter = ptvert;
				makeclippededge = 1;
			}else{
				pfrontexit = ptvert;
				makeclippededge = 1;
			}
		}else{ // add the edge to the appropriate side
			pedges->pnext = psideedges[side];
			psideedges[side] = pedges;
		}
	}
	// if anything was clipped, reconstitute and add the edges along
	// the clip plane to both sides(but in opposite directions)
	if(makeclippededge){
		if(numbedges >= (MAX_BMODEL_EDGES - 2)){
			Con_Printf("Out of edges for bmodel\n");
			return;
		}
		ptedge = &pbedges[numbedges];
		ptedge->pnext = psideedges[0];
		psideedges[0] = ptedge;
		ptedge->v[0] = pfrontexit;
		ptedge->v[1] = pfrontenter;
		ptedge = &pbedges[numbedges + 1];
		ptedge->pnext = psideedges[1];
		psideedges[1] = ptedge;
		ptedge->v[0] = pfrontenter;
		ptedge->v[1] = pfrontexit;
		numbedges += 2;
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
		}else 
			R_RecursiveClipBPoly(psideedges[i], pnode->children[i],
					     psurf);
	}
}

static bool R_BmodelFaceFront(msurface_t *psurf)
{
	mplane_t *pplane = psurf->plane; // find which node side we are on
	vec_t dot = DotProduct(modelorg, pplane->normal) - pplane->dist;
	return ((psurf->flags&SURF_PLANEBACK) && (dot < -BACKFACE_EPSILON)) ||
		(!(psurf->flags&SURF_PLANEBACK) && (dot > BACKFACE_EPSILON));
}

static void R_ClipBmodelFace(model_t *pmodel, msurface_t *psurf,
		SDL_UNUSED s32 clipflags)
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
	static mvertex_t r_bmodel_bverts[MAX_BMODEL_VERTS];
	pbverts = r_bmodel_bverts;
	static bedge_t r_bmodel_bedges[MAX_BMODEL_EDGES];
	pbedges = r_bmodel_bedges;
	numbverts = numbedges = 0;
	if(psurf->numedges <= 0)
		Sys_Error("no edges in bmodel");
	bedge_t *pbedge = &pbedges[numbedges];
	numbedges += psurf->numedges;
	medge_t *pedges = pmodel->edges;
	s32 j = 0;
	for(; j < psurf->numedges; j++){
		s32 lindex = pmodel->surfedges[psurf->firstedge + j];
		if(lindex > 0){
			medge_t *pedge = &pedges[lindex];
			pbedge[j].v[0]=&r_pcurrentvertbase[pedge->v[0]];
			pbedge[j].v[1]=&r_pcurrentvertbase[pedge->v[1]];
		}else{
			lindex = -lindex;
			medge_t *pedge = &pedges[lindex];
			pbedge[j].v[0]=&r_pcurrentvertbase[pedge->v[1]];
			pbedge[j].v[1]=&r_pcurrentvertbase[pedge->v[0]];
		}
		pbedge[j].pnext = &pbedge[j + 1];
	}
	pbedge[j - 1].pnext = NULL; // mark end of edges
	R_RecursiveClipBPoly(pbedge, currententity->topnode, psurf);
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
// behind world faces. The clip planes are in the entity's space, like the
// facebounds.
static void R_BmodelFaceNode(model_t *pmodel, mnode_t *node, s32 clipflags,
		bool occlusion, bmodelface_f face)
{
	const f32 *b = node->facebounds;
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
	if(occlusion){
		vec3_t mn, mx;
		for(s32 k = 0; k < 3; k++){
			mn[k] = b[k] + currententity->origin[k];
			mx[k] = b[3 + k] + currententity->origin[k];
		}
		if(R_CoarseOcclusionTestBox(mn, mx))
			return;
	}
	msurface_t *psurf = pmodel->surfaces + node->firstsurface;
	for(u32 i = 0; i < node->numsurfaces; i++, psurf++)
		face(pmodel, psurf, clipflags);
	for(s32 c = 0; c < 2; c++)
		if(node->children[c]->contents >= 0)
			R_BmodelFaceNode(pmodel, node->children[c], clipflags,
					occlusion, face);
}

static void R_BmodelFaces(model_t *pmodel, s32 clipflags, bmodelface_f face)
{
	if(pmodel->facetree){
		bool occlusion = r_bmodel_occ_faces &&
			!currententity->angles[0] &&
			!currententity->angles[1] && !currententity->angles[2];
		R_BmodelFaceNode(pmodel, pmodel->nodes +
				pmodel->hulls[0].firstclipnode, clipflags,
				occlusion, face);
		return;
	}
	msurface_t *psurf = &pmodel->surfaces[pmodel->firstmodelsurface];
	for(s32 i = 0; i < pmodel->nummodelsurfaces; i++, psurf++)
		face(pmodel, psurf, clipflags);
}

void R_DrawSolidClippedSubmodelPolygons(model_t *pmodel)
{
	msurface_t *psurf = &pmodel->surfaces[pmodel->firstmodelsurface];
	texture_t *tex = psurf->texinfo->texture;
	if((psurf->flags&SURF_DRAWSKY) && (tex->width/tex->height != 2))
		return; // avoid drawing broken skies
	R_BmodelFaces(pmodel, r_clipflags, R_ClipBmodelFace);
}

void R_DrawSubmodelPolygons(model_t *pmodel, s32 clipflags)
{ R_BmodelFaces(pmodel, clipflags, R_UnclippedBmodelFace); }

void R_RecursiveWorldNode(mnode_t *node, s32 clipflags)
{
	if(node->contents == CONTENTS_SOLID||node->visframe != r_visframecount)
		return;
	if(clipflags){ // cull the clipping planes if not trivial accept
		for(s32 i = 0; i < 4; i++){
			if(!(clipflags & (1 << i)))
				continue; // don't need to clip against it
			// generate accept and reject points
			s32 *pindex = pfrustum_indexes[i];
			vec3_t acceptpt, rejectpt;
			rejectpt[0] = (f32)node->minmaxs[pindex[0]];
			rejectpt[1] = (f32)node->minmaxs[pindex[1]];
			rejectpt[2] = (f32)node->minmaxs[pindex[2]];
			f64 d=DotProduct(rejectpt,view_clipplanes[i].normal);
			d -= view_clipplanes[i].dist;
			if(d <= 0)
				return;
			acceptpt[0] = (f32)node->minmaxs[pindex[3 + 0]];
			acceptpt[1] = (f32)node->minmaxs[pindex[3 + 1]];
			acceptpt[2] = (f32)node->minmaxs[pindex[3 + 2]];
			d = DotProduct(acceptpt, view_clipplanes[i].normal);
			d -= view_clipplanes[i].dist;
			if(d >= 0) // node is entirely on screen
				clipflags &= ~(1 << i);
		}
	}
	// hidden behind the world faces drawn so far
	if(R_CoarseOcclusionTestNode(node)){
		r_nodes_occluded++;
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
		if(pleaf->efrags) // deal with model fragments in this leaf
			R_StoreEfrags(&pleaf->efrags);
		pleaf->key = r_currentkey;
		r_currentkey++; // all bmodels in a leaf share the same key
	}else{
		// node is just a decision point so go down the apropriate sides
		// find which side of the node we are on
		mplane_t *plane = node->plane;
		f64 dot;
		switch(plane->type){
		case PLANE_X: dot = modelorg[0] - plane->dist; break;
		case PLANE_Y: dot = modelorg[1] - plane->dist; break;
		case PLANE_Z: dot = modelorg[2] - plane->dist; break;
		default:
			dot = DotProduct(modelorg, plane->normal) - plane->dist;
			break;
		}
		s32 side = dot >= 0 ? 0 : 1;
		// recurse down the children, front side first
		R_RecursiveWorldNode(node->children[side], clipflags);
		s32 c = node->numsurfaces; // draw stuff
		if(c){
			msurface_t *surf = cl.worldmodel->surfaces 
						+ node->firstsurface;
			if(dot < -BACKFACE_EPSILON){
				do{
					if((surf->flags & SURF_DRAWCUTOUT) ||
						((surf->flags&SURF_PLANEBACK) &&
						(surf->visframe==r_framecount)))
						R_RenderFace(surf, clipflags);
					surf++;
				}while(--c);
			}else if(dot > BACKFACE_EPSILON){
				do{
					if((surf->flags & SURF_DRAWCUTOUT) ||
						(!(surf->flags&(SURF_PLANEBACK|
						SURF_NODRAWFRONT))&&
						(surf->visframe==r_framecount)))
						R_RenderFace(surf, clipflags);
					surf++;
				}while(--c);
			}
		// all surfaces on the same node share the same sequence number
			r_currentkey++;
		}
		// recurse down the back side
		R_RecursiveWorldNode(node->children[!side], clipflags);
	}
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
		R_RecursiveWorldNode(clmodel->nodes, 15);
	}
	if(r_skyframe == r_framecount)
		R_EmitSkyBox();
}
