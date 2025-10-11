// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

static u32 cacheoffset;
static medge_t tedge;
static medge_t *r_pedge;
static bool r_leftclipped, r_rightclipped;
static bool r_nearzionly;
static mvertex_t r_leftenter, r_leftexit;
static mvertex_t r_rightenter, r_rightexit;
static s32 r_emitted;
static f32 r_nearzi;
static f32 r_u1, r_v1, r_lzi1;
static s32 r_ceilv1;
static bool r_lastvertvalid;
static bool makeleftedge, makerightedge;
extern f32 cur_ent_alpha;

void R_EmitEdge(mvertex_t *pv0, mvertex_t *pv1)
{
	vec3_t local, transformed;
	f32 *world;
	s64 v, v2, ceilv0;
	f32 scale, lzi0, u0, v0;
	if (r_lastvertvalid) {
		u0 = r_u1;
		v0 = r_v1;
		lzi0 = r_lzi1;
		ceilv0 = r_ceilv1;
	} else {
		world = &pv0->position[0]; // transform and project
		VectorSubtract(world, modelorg, local);
		TransformVector(local, transformed);
		if (transformed[2] < NEAR_CLIP)
			transformed[2] = NEAR_CLIP;
		lzi0 = 1.0 / transformed[2];
		// FIXME: build x/yscale into transform?
		scale = xscale * lzi0;
		u0 = (xcenter + scale * transformed[0]);
		if (u0 < r_refdef.fvrectx_adj)
			u0 = r_refdef.fvrectx_adj;
		if (u0 > r_refdef.fvrectright_adj)
			u0 = r_refdef.fvrectright_adj;
		scale = yscale * lzi0;
		v0 = (ycenter - scale * transformed[1]);
		if (v0 < r_refdef.fvrecty_adj)
			v0 = r_refdef.fvrecty_adj;
		if (v0 > r_refdef.fvrectbottom_adj)
			v0 = r_refdef.fvrectbottom_adj;
		ceilv0 = (s32)ceil(v0);
	}
	world = &pv1->position[0];
	VectorSubtract(world, modelorg, local); // transform and project
	TransformVector(local, transformed);
	if (transformed[2] < NEAR_CLIP)
		transformed[2] = NEAR_CLIP;
	r_lzi1 = 1.0 / transformed[2];
	scale = xscale * r_lzi1;
	r_u1 = (xcenter + scale * transformed[0]);
	if (r_u1 < r_refdef.fvrectx_adj)
		r_u1 = r_refdef.fvrectx_adj;
	if (r_u1 > r_refdef.fvrectright_adj)
		r_u1 = r_refdef.fvrectright_adj;
	scale = yscale * r_lzi1;
	r_v1 = (ycenter - scale * transformed[1]);
	if (r_v1 < r_refdef.fvrecty_adj)
		r_v1 = r_refdef.fvrecty_adj;
	if (r_v1 > r_refdef.fvrectbottom_adj)
		r_v1 = r_refdef.fvrectbottom_adj;
	if (r_lzi1 > lzi0)
		lzi0 = r_lzi1;
	if (lzi0 > r_nearzi) // for mipmap finding
		r_nearzi = lzi0;
	// for right edges, all we want is the effect on 1/z
	if (r_nearzionly)
		return;
	r_emitted = 1;
	r_ceilv1 = (s32)ceil(r_v1);
	if (ceilv0 == r_ceilv1) { // create the edge
		// we cache unclipped horizontal edges as fully clipped
		if (cacheoffset != 0x7FFFFFFF) {
			cacheoffset = FULLY_CLIPPED_CACHED |
				(r_framecount & FRAMECOUNT_MASK);
		}
		return; // horizontal edge
	}
	s64 side = ceilv0 > r_ceilv1;
	edge_t *edge = edge_p++;
	edge->owner = r_pedge;
	edge->nearzi = lzi0;
	f64 u, u_step;
	if (side == 0) { // trailing edge (go from p1 to p2)
		v = ceilv0;
		v2 = r_ceilv1 - 1;
		edge->surfs[0] = surface_p - surfaces;
		edge->surfs[1] = 0;
		u_step = ((r_u1 - u0) / (r_v1 - v0));
		u = u0 + ((f32)v - v0) * u_step;
	} else { // leading edge (go from p2 to p1)
		v2 = ceilv0 - 1;
		v = r_ceilv1;
		edge->surfs[0] = 0;
		edge->surfs[1] = surface_p - surfaces;
		u_step = ((u0 - r_u1) / (v0 - r_v1));
		u = r_u1 + ((f32)v - r_v1) * u_step;
	}
	edge->u_step = u_step * 0x100000;
	edge->u = u * 0x100000 + 0xFFFFF;
	// we need to do this to avoid stepping off the edges if a very nearly
	// horizontal edge is less than epsilon above a scan, and numeric error causes
	// it to incorrectly extend to the scan, and the extension of the line goes off
	// the edge of the screen
	// FIXME: is this actually needed?
	if (edge->u < r_refdef.vrect_x_adj_shift20)
		edge->u = r_refdef.vrect_x_adj_shift20;
	if (edge->u > r_refdef.vrectright_adj_shift20)
		edge->u = r_refdef.vrectright_adj_shift20;
	s64 u_check = edge->u; // sort the edge in normally
	if (edge->surfs[0])
		u_check++; // sort trailers after leaders
	// CyanBun96: denser maps like mg1 start begin to chug on the linear
	// edge list traversal as it gets longer. Remembering the last insert
	// point and checking if it fits helps a lot in the best case, and hurts
	// very little in the worst (where performance improvements aren't
	// needed anyway).
	edge_t *pcheck = last_pcheck[v];
	if (!newedges[v] || newedges[v]->u >= u_check) { // Insert at the head
		edge->next = newedges[v];
		newedges[v] = edge;
		last_pcheck[v] = newedges[v]; // cache new head
	} else if (pcheck && pcheck->u < u_check) { // Use cached position
		while (pcheck->next && pcheck->next->u < u_check)
			pcheck = pcheck->next;
		edge->next = pcheck->next;
		pcheck->next = edge;
		last_pcheck[v] = edge; // cache last inserted edge
	} else { // Fallback to normal head traversal
		pcheck = newedges[v];
		while (pcheck->next && pcheck->next->u < u_check)
			pcheck = pcheck->next;
		edge->next = pcheck->next;
		pcheck->next = edge;
		last_pcheck[v] = edge;
	}
	edge->nextremove = removeedges[v2];
	removeedges[v2] = edge;
}

void R_ClipEdge(mvertex_t *pv0, mvertex_t *pv1, clipplane_t *clip)
{
	f32 d0, d1, f;
	mvertex_t clipvert;
	if (clip) {
		do {
			d0 = DotProduct(pv0->position,clip->normal)-clip->dist;
			d1 = DotProduct(pv1->position,clip->normal)-clip->dist;
			if (d0 >= 0) { // point 0 is unclipped
				if (d1 >= 0) {
					continue; // both points are unclipped
				}
				// only point 1 is clipped
				// we don't cache clipped edges
				cacheoffset = 0x7FFFFFFF;
				f = d0 / (d0 - d1);
				clipvert.position[0]=pv0->position[0]+f*
					(pv1->position[0]-pv0->position[0]);
				clipvert.position[1]=pv0->position[1]+f*
					(pv1->position[1]-pv0->position[1]);
				clipvert.position[2]=pv0->position[2]+f*
					(pv1->position[2]-pv0->position[2]);
				if (clip->leftedge) {
					r_leftclipped = 1;
					r_leftexit = clipvert;
				} else if (clip->rightedge) {
					r_rightclipped = 1;
					r_rightexit = clipvert;
				}
				R_ClipEdge(pv0, &clipvert, clip->next);
				return;
			} else { // point 0 is clipped
				if (d1 < 0) { // both points are clipped
					if (!r_leftclipped) // we do cache fully clipped edges
						cacheoffset = FULLY_CLIPPED_CACHED | (r_framecount & FRAMECOUNT_MASK);
					return;
				}
				// only point 0 is clipped
				r_lastvertvalid = 0;
				// we don't cache partially clipped edges
				cacheoffset = 0x7FFFFFFF;
				f = d0 / (d0 - d1);
				clipvert.position[0] = pv0->position[0] + f * (pv1->position[0] - pv0->position[0]);
				clipvert.position[1] = pv0->position[1] + f * (pv1->position[1] - pv0->position[1]);
				clipvert.position[2] = pv0->position[2] + f * (pv1->position[2] - pv0->position[2]);
				if (clip->leftedge) {
					r_leftclipped = 1;
					r_leftenter = clipvert;
				} else if (clip->rightedge) {
					r_rightclipped = 1;
					r_rightenter = clipvert;
				}
				R_ClipEdge(&clipvert, pv1, clip->next);
				return;
			}
		} while ((clip = clip->next) != NULL);
	}
	R_EmitEdge(pv0, pv1); // add the edge
}

void R_EmitCachedEdge()
{
	edge_t *pedge_t=(edge_t*)((uintptr_t)r_edges+r_pedge->cachededgeoffset);
	if (!pedge_t->surfs[0])
		pedge_t->surfs[0] = surface_p - surfaces;
	else
		pedge_t->surfs[1] = surface_p - surfaces;
	if (pedge_t->nearzi > r_nearzi) // for mipmap finding
		r_nearzi = pedge_t->nearzi;
	r_emitted = 1;
}

void R_RenderFace(msurface_t *fa, s32 clipflags)
{
	// Manoel Kasimier - skyboxes - begin
	// Code taken from the ToChriS engine - Author: Vic
	// sky surfaces encountered in the world will cause the
	// environment box surfaces to be emited
	if ((fa->flags & SURF_DRAWSKY) && skybox_name[0]) {
		R_EmitSkyBox();
		return;
	}
	// Manoel Kasimier - skyboxes - end
	if ((s32)r_twopass.value&1 && fa->flags&SURF_WINQUAKE_DRAWTRANSLUCENT) {
		winquake_surface_liquid_alpha=R_LiquidAlphaForFlags(fa->flags);
		if (r_litwater.value) winquake_surface_liquid_alpha = 0.5; // FIXME properly
	} else if (cur_ent_alpha < 1 && r_entalpha.value == 1)
		winquake_surface_liquid_alpha = cur_ent_alpha;
	else winquake_surface_liquid_alpha = 1;
        // Baker: Fully transparent = invisible = don't render
        if (!winquake_surface_liquid_alpha)
                return;
        // Baker: If this is the alpha water pass and we aren't alpha water, get out!
        if (r_wateralphapass && winquake_surface_liquid_alpha == 1)
                return;
        if (!r_wateralphapass && winquake_surface_liquid_alpha < 1) {
                r_foundtranswater = 1;
                return;
        }
	if ((surface_p) >= surf_max) { // skip out if no more surfs
		r_outofsurfaces++;
		return;
	}
	// ditto if not enough edges left
	if ((edge_p + fa->numedges + 4) >= edge_max) {
		r_outofedges += fa->numedges;
		return;
	}
	c_faceclip++;
	clipplane_t *pclip = NULL; // set up clip planes
	s32 i = 3;
	u32 mask = 0x08;
	for (; i >= 0; i--, mask >>= 1) {
		if (clipflags & mask) {
			view_clipplanes[i].next = pclip;
			pclip = &view_clipplanes[i];
		}
	}
	r_emitted = 0; // push the edges through
	r_nearzi = 0;
	r_nearzionly = 0;
	makeleftedge = makerightedge = 0;
	medge_t *pedges, tedge;
	pedges = currententity->model->edges;
	r_lastvertvalid = 0;
	for (i = 0; i < fa->numedges; i++) {
		s32 lindex = currententity->model->surfedges[fa->firstedge + i];
		if (lindex > 0) {
			r_pedge = &pedges[lindex];
			// if the edge is cached, we can just reuse the edge
			if (!insubmodel) {
				if (r_pedge-> cachededgeoffset & FULLY_CLIPPED_CACHED) {
					if ((r_pedge-> cachededgeoffset & FRAMECOUNT_MASK)
							== r_framecount) {
						r_lastvertvalid = 0;
						continue;
					}
				} else {
					if ((((uintptr_t) edge_p - (uintptr_t) r_edges) > r_pedge->cachededgeoffset) &&
							(((edge_t *) ((uintptr_t) r_edges + r_pedge-> cachededgeoffset))-> owner == r_pedge)) {
						R_EmitCachedEdge();
						r_lastvertvalid = 0;
						continue;
					}
				}
			}
			// assume it's cacheable
			cacheoffset = (u32 *) edge_p - (u32 *) r_edges;
			r_leftclipped = r_rightclipped = 0;
			R_ClipEdge(&r_pcurrentvertbase[r_pedge->v[0]],
					&r_pcurrentvertbase[r_pedge->v[1]], pclip);
			r_pedge->cachededgeoffset = cacheoffset;
			if (r_leftclipped)
				makeleftedge = 1;
			if (r_rightclipped)
				makerightedge = 1;
			r_lastvertvalid = 1;
		} else {
			lindex = -lindex;
			r_pedge = &pedges[lindex];
			// if the edge is cached, we can just reuse the edge
			if (!insubmodel) {
				if (r_pedge-> cachededgeoffset & FULLY_CLIPPED_CACHED) {
					if ((r_pedge-> cachededgeoffset & FRAMECOUNT_MASK) == r_framecount) {
						r_lastvertvalid = 0;
						continue;
					}
				} else {
					// it's cached if the cached edge is valid and is owned
					// by this medge_t
					if ((((uintptr_t) edge_p - (uintptr_t) r_edges) > r_pedge->cachededgeoffset) &&
							(((edge_t *) ((uintptr_t) r_edges + r_pedge-> cachededgeoffset))-> owner == r_pedge)) {
						R_EmitCachedEdge();
						r_lastvertvalid = 0;
						continue;
					}
				}
			}
			// assume it's cacheable
			cacheoffset = (u32 *) edge_p - (u32 *) r_edges;
			r_leftclipped = r_rightclipped = 0;
			R_ClipEdge(&r_pcurrentvertbase[r_pedge->v[1]],
					&r_pcurrentvertbase[r_pedge->v[0]], pclip);
			r_pedge->cachededgeoffset = cacheoffset;
			if (r_leftclipped)
				makeleftedge = 1;
			if (r_rightclipped)
				makerightedge = 1;
			r_lastvertvalid = 1;
		}
	}
	// if there was a clip off the left edge, add that edge too
	// FIXME: faster to do in screen space?
	// FIXME: share clipped edges?
	if (makeleftedge) {
		r_pedge = &tedge;
		r_lastvertvalid = 0;
		R_ClipEdge(&r_leftexit, &r_leftenter, pclip->next);
	}
	// if there was a clip off the right edge, get the right r_nearzi
	if (makerightedge) {
		r_pedge = &tedge;
		r_lastvertvalid = 0;
		r_nearzionly = 1;
		R_ClipEdge(&r_rightexit, &r_rightenter,view_clipplanes[1].next);
	}
	// if no edges made it out, return without posting the surface
	if (!r_emitted)
		return;
	r_polycount++;
	surface_p->data = (void *)fa;
	surface_p->nearzi = r_nearzi;
	surface_p->flags = fa->flags;
	surface_p->insubmodel = insubmodel;
	surface_p->spanstate = 0;
	surface_p->entity = currententity;
	surface_p->key = r_currentkey++;
	surface_p->spans = NULL;
	mplane_t *pplane = fa->plane;
	// FIXME: cache this?
	vec3_t p_normal;
	TransformVector(pplane->normal, p_normal);
	// FIXME: cache this?
	f32 distinv = 1.0/(pplane->dist-DotProduct(modelorg, pplane->normal));
	surface_p->d_zistepu = p_normal[0] * xscaleinv * distinv;
	surface_p->d_zistepv = -p_normal[1] * yscaleinv * distinv;
	surface_p->d_ziorigin = p_normal[2] * distinv -
		xcenter * surface_p->d_zistepu - ycenter * surface_p->d_zistepv;
	surface_p++;
}

void R_RenderBmodelFace(bedge_t *pedges, msurface_t *psurf)
{
	if (surface_p >= surf_max) { // skip out if no more surfs
		r_outofsurfaces++;
		return;
	}
	// ditto if not enough edges left
	if ((edge_p + psurf->numedges + 4) >= edge_max) {
		r_outofedges += psurf->numedges;
		return;
	}
	c_faceclip++;
	// this is a dummy to give the caching mechanism someplace to write to
	r_pedge = &tedge;
	clipplane_t *pclip = NULL; // set up clip planes
	s32 i = 3; 
	u32 mask = 0x08;
	for (; i >= 0; i--, mask >>= 1) {
		if (r_clipflags & mask) {
			view_clipplanes[i].next = pclip;
			pclip = &view_clipplanes[i];
		}
	}
	r_emitted = 0; // push the edges through
	r_nearzi = 0;
	r_nearzionly = 0;
	makeleftedge = makerightedge = 0;
	// FIXME: keep clipped bmodel edges in clockwise order so last vertex
	// caching can be used?
	r_lastvertvalid = 0;
	for (; pedges; pedges = pedges->pnext) {
		r_leftclipped = r_rightclipped = 0;
		R_ClipEdge(pedges->v[0], pedges->v[1], pclip);
		if (r_leftclipped)
			makeleftedge = 1;
		if (r_rightclipped)
			makerightedge = 1;
	}
	// if there was a clip off the left edge, add that edge too
	// FIXME: faster to do in screen space?
	// FIXME: share clipped edges?
	if (makeleftedge) {
		r_pedge = &tedge;
		R_ClipEdge(&r_leftexit, &r_leftenter, pclip->next);
	}
	// if there was a clip off the right edge, get the right r_nearzi
	if (makerightedge) {
		r_pedge = &tedge;
		r_nearzionly = 1;
		R_ClipEdge(&r_rightexit, &r_rightenter,
				view_clipplanes[1].next);
	}
	// if no edges made it out, return without posting the surface
	if (!r_emitted)
		return;
	r_polycount++;
	surface_p->data = (void *)psurf;
	surface_p->nearzi = r_nearzi;
	surface_p->flags = psurf->flags;
	surface_p->insubmodel = 1;
	surface_p->spanstate = 0;
	surface_p->entity = currententity;
	surface_p->key = r_currentbkey;
	surface_p->spans = NULL;
	mplane_t *pplane = psurf->plane;
	// FIXME: cache this?
	vec3_t p_normal;
	TransformVector(pplane->normal, p_normal);
	// FIXME: cache this?
	f32 distinv = 1.0/(pplane->dist-DotProduct(modelorg,pplane->normal));
	surface_p->d_zistepu = p_normal[0] * xscaleinv * distinv;
	surface_p->d_zistepv = -p_normal[1] * yscaleinv * distinv;
	surface_p->d_ziorigin = p_normal[2] * distinv -
		xcenter * surface_p->d_zistepu - ycenter * surface_p->d_zistepv;
	surface_p++;
}
