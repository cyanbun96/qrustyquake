// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

// The world walk's emit job uses them on a worker (R_EmitFaceJob)
static Q_TLS u32 cacheoffset;
static Q_TLS medge_t *r_pedge;
static Q_TLS bool r_leftclipped, r_rightclipped;
static Q_TLS bool r_nearzionly;
static Q_TLS mvertex_t r_leftenter, r_leftexit;
static Q_TLS mvertex_t r_rightenter, r_rightexit;
static Q_TLS s32 r_emitted;
static Q_TLS f32 r_nearzi;
static Q_TLS bool r_lastvertvalid;
static Q_TLS bool makeleftedge, makerightedge;
static Q_TLS u8 r_edgerev; // direction of the face currently being walked

// The order WinQuake's compiler used, which matters under chop rounding, and
// the same goes for the order of the multiplications by distinv below
static inline f32 R_PlaneDot(const mplane_t *pplane)
{
	const f32 *n = pplane->normal;
	return (n[2] * modelorg[2] + n[1] * modelorg[1]) + n[0] * modelorg[0];
}

void R_EmitEdge(mvertex_t *pv0, mvertex_t *pv1)
{
	vec3_t local, transformed;
	f32 *world;
	f32 scale, lzi0, u0, v0;
	s64 ceilv0;
	static Q_TLS f32 r_u1, r_v1, r_lzi1;
	static Q_TLS s32 r_ceilv1;
	if(r_lastvertvalid){
		u0 = r_u1;
		v0 = r_v1;
		lzi0 = r_lzi1;
		ceilv0 = r_ceilv1;
	}else{
		world = &pv0->position[0]; // transform and project
		VectorSubtract(world, modelorg, local);
		TransformVector(local, transformed);
		if(transformed[2] < NEAR_CLIP)
			transformed[2] = NEAR_CLIP;
		lzi0 = 1.0 / transformed[2];
		scale = xscale * lzi0;
		u0 = (xcenter + scale * transformed[0]);
		if(u0 < r_refdef.fvrectx_adj)
			u0 = r_refdef.fvrectx_adj;
		if(u0 > r_refdef.fvrectright_adj)
			u0 = r_refdef.fvrectright_adj;
		scale = yscale * lzi0;
		v0 = (ycenter - scale * transformed[1]);
		if(v0 < r_refdef.fvrecty_adj)
			v0 = r_refdef.fvrecty_adj;
		if(v0 > r_refdef.fvrectbottom_adj)
			v0 = r_refdef.fvrectbottom_adj;
		ceilv0 = Q_ceil(v0);
	}
	world = &pv1->position[0];
	VectorSubtract(world, modelorg, local); // transform and project
	TransformVector(local, transformed);
	if(transformed[2] < NEAR_CLIP)
		transformed[2] = NEAR_CLIP;
	r_lzi1 = 1.0 / transformed[2];
	scale = xscale * r_lzi1;
	r_u1 = (xcenter + scale * transformed[0]);
	if(r_u1 < r_refdef.fvrectx_adj)
		r_u1 = r_refdef.fvrectx_adj;
	if(r_u1 > r_refdef.fvrectright_adj)
		r_u1 = r_refdef.fvrectright_adj;
	scale = yscale * r_lzi1;
	r_v1 = (ycenter - scale * transformed[1]);
	if(r_v1 < r_refdef.fvrecty_adj)
		r_v1 = r_refdef.fvrecty_adj;
	if(r_v1 > r_refdef.fvrectbottom_adj)
		r_v1 = r_refdef.fvrectbottom_adj;
	if(r_lzi1 > lzi0)
		lzi0 = r_lzi1;
	if(lzi0 > r_nearzi) // for mipmap finding
		r_nearzi = lzi0;
	if(r_nearzionly) // for right edges, all we want is the effect on 1/z
		return;
	r_emitted = 1;
	r_ceilv1 = Q_ceil(r_v1);
	if(ceilv0 == r_ceilv1){ // create the edge
		// we cache unclipped horizontal edges as fully clipped
		if(cacheoffset != 0x7FFFFFFF){
			cacheoffset = FULLY_CLIPPED_CACHED |
				(r_framecount & FRAMECOUNT_MASK);
		}
		return; // horizontal edge
	}
	s64 side = ceilv0 > r_ceilv1;
	edge_t *edge = edge_p++;
	edge->owner = r_pedge;
	edge->rev = r_edgerev;
	edge->nearzi = lzi0;
	f64 u, u_step;
	s64 v, v2;
	if(side == 0){ // trailing edge(go from p1 to p2)
		v = ceilv0;
		v2 = r_ceilv1 - 1;
		edge->surfs[0] = surface_p - surfaces;
		edge->surfs[1] = 0;
		u_step = ((r_u1 - u0) / (r_v1 - v0));
		u = u0 + ((f32)v - v0) * u_step;
	}else{ // leading edge(go from p2 to p1)
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
	// horizontal edge is less than epsilon above a scan, and numeric error
	// causes it to incorrectly extend to the scan, and the extension of the
	// line goes off the edge of the screen
	if(edge->u < r_refdef.vrect_x_adj_shift20)
		edge->u = r_refdef.vrect_x_adj_shift20;
	if(edge->u > r_refdef.vrectright_adj_shift20)
		edge->u = r_refdef.vrectright_adj_shift20;
	// R_ScanEdges sorts each row once it gets to it
	edge->next = newedges[v];
	newedges[v] = edge;
	edge->vlast = (s32)v2;
	if(r_scanbanded){ // the bands below this one start with it
		edge->vstart = (s32)v;
		if(r_bandof[v2] != r_bandof[v])
			R_CarryEdge(edge, r_bandof[v], r_bandof[v2]);
	}
}

void R_ClipEdge(mvertex_t *pv0, mvertex_t *pv1, clipplane_t *clip)
{
	f32 d0, d1, f;
	mvertex_t clipvert;
	if(clip){
		do{
			d0 = DotProduct(pv0->position,clip->normal)-clip->dist;
			d1 = DotProduct(pv1->position,clip->normal)-clip->dist;
			if(d0 >= 0){ // point 0 is unclipped
				if(d1 >= 0){
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
				if(clip->leftedge){
					r_leftclipped = 1;
					r_leftexit = clipvert;
				}else if(clip->rightedge){
					r_rightclipped = 1;
					r_rightexit = clipvert;
				}
				R_ClipEdge(pv0, &clipvert, clip->next);
				return;
			}else{ // point 0 is clipped
				if(d1 < 0){ // both points are clipped
					if(!r_leftclipped) // fully clipped edge
						cacheoffset=FULLY_CLIPPED_CACHED
						|(r_framecount&FRAMECOUNT_MASK);
					return;
				}
				// only point 0 is clipped
				r_lastvertvalid = 0;
				// we don't cache partially clipped edges
				cacheoffset = 0x7FFFFFFF;
				f = d0 / (d0 - d1);
				clipvert.position[0] = pv0->position[0] + f *
					(pv1->position[0] - pv0->position[0]);
				clipvert.position[1] = pv0->position[1] + f *
					(pv1->position[1] - pv0->position[1]);
				clipvert.position[2] = pv0->position[2] + f *
					(pv1->position[2] - pv0->position[2]);
				if(clip->leftedge){
					r_leftclipped = 1;
					r_leftenter = clipvert;
				}else if(clip->rightedge){
					r_rightclipped = 1;
					r_rightenter = clipvert;
				}
				R_ClipEdge(&clipvert, pv1, clip->next);
				return;
			}
		}while((clip = clip->next) != NULL);
	}
	R_EmitEdge(pv0, pv1); // add the edge
}

// The view clip planes in clipflags, linked in order
static inline clipplane_t *R_ClipPlaneList(s32 clipflags)
{
	clipplane_t *pclip = NULL;
	for(s32 i = 3; i >= 0; i--){
		if(clipflags & (1 << i)){
			view_clipplanes[i].next = pclip;
			pclip = &view_clipplanes[i];
		}
	}
	return pclip;
}

void R_EmitCachedEdge(edge_t *pedge_t)
{
	if(!pedge_t->surfs[0])
		pedge_t->surfs[0] = surface_p - surfaces;
	else
		pedge_t->surfs[1] = surface_p - surfaces;
	if(pedge_t->nearzi > r_nearzi) // for mipmap finding
		r_nearzi = pedge_t->nearzi;
	r_emitted = 1;
}

// Does the face get drawn in this pass? Sets winquake_surface_liquid_alpha.
bool R_SurfaceInPass(msurface_t *fa)
{
	if(fa->flags&SURF_WINQUAKE_DRAWTRANSLUCENT && !currententity->alpha){
		winquake_surface_liquid_alpha=R_LiquidAlphaForFlags(fa->flags);
	}else if(cur_ent_alpha < 1 && r_entalpha.value == 1)
		winquake_surface_liquid_alpha = cur_ent_alpha;
	else winquake_surface_liquid_alpha = 1;
	// Baker: Fully transparent = invisible = don't render
	if(!winquake_surface_liquid_alpha)
		return 0;
// Baker: If this is the alpha water pass and we aren't alpha water, get out!
	if(r_alphapass && winquake_surface_liquid_alpha == 1)
		return 0;
	if(!r_alphapass && winquake_surface_liquid_alpha < 1){
		r_foundtranswater = 1;
		return 0;
	}
	return 1;
}

// The world walk with the workers hands its faces to a job that emits them in
// the same order (R_EmitFaceJob), and goes on with the next ones. The
// occlusion fill of a face has to know that it went in. R_EmitSure, and the
// room left as far as the job has got, can tell for most, else the walk waits.
#define R_EMITQ 4096 // a power of two
#define R_QBATCH 32 // the faces the job gets at a time
typedef struct {
	msurface_t *fa;
	s32 clipflags, key;
	s64 edges; // at most this many edges in the faces before it
} qface_t;
// The pads keep what one side writes off the lines the other one reads
static struct {
	u8 pad0[64];
	qface_t q[R_EMITQ]; // the walk's
	s32 put; // faces queued, head lags behind by up to R_QBATCH
	s32 seen; // what the walk knows of the job: faces done,
	s64 seensurfs, seenedges; // and surface_p, edge_p after them
	s64 edgesum;
	s32 job;
	bool queue; // the faces go to the job
	u8 pad1[64];
	SDL_AtomicInt head, stop; // published by the walk
	u8 pad2[64];
	bool in[R_EMITQ]; // the job's: whether it went in
	SDL_AtomicInt tail, surfs, edges; // faces done, surface_p, edge_p
	u8 pad3[64];
	bool walking; // in the world walk, see R_RenderFace
	clipplane_t planes[16][4]; // as R_ClipPlaneList links them, for
	clipplane_t *clip[16]; // each clipflags, during the walk
	edge_t *edge_p; // the job's edge_p, surface_p and r_polycount at its
	surf_t *surface_p; // start and end, see R_WalkEnd
	s32 polys;
	u8 pad4[64];
} r_q = { .job = -1 };

// The plane of a face as 1/z on the screen. Once, for the edge list and for
// the occlusion fill alike.
static Q_NOINLINE void R_FaceZi(msurface_t *fa, surf_t *s)
{
	mplane_t *pplane = fa->plane;
	vec3_t p_normal;
	TransformVector(pplane->normal, p_normal);
	f32 distinv = 1.0 / (pplane->dist - R_PlaneDot(pplane));
	s->d_zistepu = (p_normal[0] * distinv) * xscaleinv;
	s->d_zistepv = -((p_normal[1] * distinv) * yscaleinv);
	s->d_ziorigin = p_normal[2] * distinv -
		xcenter * s->d_zistepu - ycenter * s->d_zistepv;
}

// Is there room for a surface and the edges of a face? One that doesn't fit
// is counted, the buffers grow by that (V_RenderView).
Q_SIL bool R_FaceFits(s32 numedges)
{
	if(surface_p >= surf_max || edge_p + numedges + 4 >= edge_max){
		r_outofsurfaces++;
		return 0;
	}
	return 1;
}

// The edges and the surface of a face that got past the occlusion test, with
// key, or the next one if that's -1. False if it didn't go in.
Q_SIL bool R_EmitFace(msurface_t *fa, s32 clipflags, s32 key)
{
	if(!R_FaceFits(fa->numedges))
		return 0;
	// the walk's are copies, the job doesn't write the planes it reads
	clipplane_t *planes = r_q.walking ? r_q.planes[clipflags]
		: view_clipplanes;
	clipplane_t *pclip = r_q.walking ? r_q.clip[clipflags]
		: R_ClipPlaneList(clipflags);
	r_emitted = 0; // push the edges through
	r_nearzi = 0;
	r_nearzionly = 0;
	makeleftedge = makerightedge = 0;
	medge_t *pedges, tedge;
	pedges = currententity->model->edges;
	r_lastvertvalid = 0;
	s32 i;
	if(r_showtris.value){ // UNIFIED BSP GEOMETRY HOOK
		mvertex_t *first_vert = NULL;
		mvertex_t *prev_vert = NULL;
		for(i = 0; i < fa->numedges; i++){
			s32 li=currententity->model->surfedges[fa->firstedge+i];
			mvertex_t *vert;
			if(li > 0) vert = &r_pcurrentvertbase[pedges[li].v[0]];
			else vert = &r_pcurrentvertbase[pedges[-li].v[1]];
			if(i == 0){
				first_vert = vert;
			}else{
				// 1. Draw perimeter edge
				R_DrawDebugLine3D(prev_vert->position,
						vert->position);
				// 2. Draw inner triangle diagonal if in mode 1
				if(r_showtris.value == 1 && i > 1 &&
						i < fa->numedges - 1)
					R_DrawDebugLine3D(first_vert->position,
							vert->position);
			}
			prev_vert = vert;
		}
		if(fa->numedges >= 3){ // 3. Close the perimeter loop
			R_DrawDebugLine3D(prev_vert->position,
					first_vert->position);
		}
	}
	const s32 *sedges = currententity->model->surfedges + fa->firstedge;
	mvertex_t *vbase = r_pcurrentvertbase;
	const bool nosub = !insubmodel;
	const u32 fstamp = r_framecount & FRAMECOUNT_MASK;
	for(i = 0; i < fa->numedges; i++){
		s32 lindex = sedges[i];
		s32 rev = lindex <= 0;	// 1 = walks the edge v[1] -> v[0]
		r_pedge = &pedges[rev ? -lindex : lindex];
		// if the edge is cached, we can just reuse the edge
		if(nosub){
			u32 off = r_pedge->cachededgeoffset;
			if(off & FULLY_CLIPPED_CACHED){
				if((off & FRAMECOUNT_MASK) == fstamp){
					r_lastvertvalid = 0;
					continue;
				}
			}else if((uintptr_t)edge_p - (uintptr_t)r_edges > off){
				edge_t *ce = (edge_t *)((uintptr_t)r_edges+off);
				if(ce->owner == r_pedge){
				// Share only if the creator walked the edge the
				// other way and a slot is still free.
					if(ce->rev == (u8)(rev ^ 1) &&
						!(ce->surfs[0]&&ce->surfs[1])){
						R_EmitCachedEdge(ce);
						r_lastvertvalid = 0;
						continue;
					}
				//Same direction or already shared: fall through
				//and build a private edge_t for this face.
				}
			}
		}
		// assume it's cacheable
		cacheoffset = (u32)((uintptr_t)edge_p - (uintptr_t)r_edges);
		r_edgerev = (u8)rev; // R_EmitEdge stores this in the new edge_t
		r_leftclipped = r_rightclipped = 0;
		R_ClipEdge(&vbase[r_pedge->v[rev]], &vbase[r_pedge->v[rev ^ 1]],
				pclip);
		r_pedge->cachededgeoffset = cacheoffset;
		if(r_leftclipped)
			makeleftedge = 1;
		if(r_rightclipped)
			makerightedge = 1;
		r_lastvertvalid = 1;
	}
	// if there was a clip off the left edge, add that edge too
	if(makeleftedge){
		r_pedge = &tedge;
		r_lastvertvalid = 0;
		R_ClipEdge(&r_leftexit, &r_leftenter, pclip->next);
	}
	// if there was a clip off the right edge, get the right r_nearzi
	if(makerightedge){
		r_pedge = &tedge;
		r_lastvertvalid = 0;
		r_nearzionly = 1;
		R_ClipEdge(&r_rightexit, &r_rightenter, planes[1].next);
	}
	// if no edges made it out, return without posting the surface
	if(!r_emitted)
		return 0;
	r_polycount++;
	surface_p->data = (void *)fa;
	surface_p->nearzi = r_nearzi;
	surface_p->flags = fa->flags;
	surface_p->insubmodel = insubmodel;
	surface_p->spanstate = 0;
	surface_p->entity = currententity;
	surface_p->key = key < 0 ? r_currentkey++ : key;
	surface_p->spans = NULL;
	R_FaceZi(fa, surface_p);
	surface_p++;
	return 1;
}

// Will R_EmitFace put the face in, whatever the faces before it left in the
// edge cache? It will once R_EmitEdge gets an edge that no face can have
// cached as fully clipped: one that a clip plane cuts and that keeps a piece
// inside all four (R_ClipEdge), or one inside all four that crosses a row.
// The margins cover rounding other than theirs.
bool R_EmitSure(msurface_t *fa)
{
	const model_t *m = currententity->model;
	const s32 *sedges = m->surfedges + fa->firstedge;
	for(s32 i = 0; i < fa->numedges; i++){
		s32 l = sedges[i];
		const medge_t *e = &m->edges[l < 0 ? -l : l];
		const f32 *p0 = r_pcurrentvertbase[e->v[0]].position;
		const f32 *p1 = r_pcurrentvertbase[e->v[1]].position;
		f32 lo = 0, hi = 1; // the piece inside all four by the margin
		bool cut = 0;
		for(s32 c = 0; c < 4 && lo < hi; c++){
			const clipplane_t *cp = &view_clipplanes[c];
			f32 a0 = DotProduct(p0, cp->normal);
			f32 a1 = DotProduct(p1, cp->normal);
			f32 mg = q_max(fabsf(a0), fabsf(a1)) + fabsf(cp->dist);
			mg = mg * 1e-4f + 1e-3f;
			f32 d0 = a0 - cp->dist, d1 = a1 - cp->dist;
			if(d0 > mg && d1 > mg)
				continue;
			cut = 1;
			if(d0 > mg && d1 < -mg)
				hi = q_min(hi, (d0 - mg) / (d0 - d1));
			else if(d0 < -mg && d1 > mg)
				lo = q_max(lo, (mg - d0) / (d1 - d0));
			else
				hi = 0; // outside, or too close to call
		}
		if(lo >= hi)
			continue;
		if(cut)
			return 1;
		f32 v[2];
		for(s32 k = 0; k < 2; k++){
			vec3_t local, t;
			VectorSubtract(k ? p1 : p0, modelorg, local);
			TransformVector(local, t);
			if(t[2] < NEAR_CLIP)
				t[2] = NEAR_CLIP;
			f32 lzi = 1.0 / t[2];
			f32 scale = yscale * lzi;
			v[k] = ycenter - scale * t[1];
			v[k] = q_max(v[k], r_refdef.fvrecty_adj);
			v[k] = q_min(v[k], r_refdef.fvrectbottom_adj);
		}
		f32 vlo = q_min(v[0], v[1]) + 0.01f;
		f32 vhi = q_max(v[0], v[1]) - 0.01f;
		if(Q_ceil(vlo) < vhi) // a whole row between them
			return 1;
	}
	return 0;
}

static Q_NOINLINE bool R_EmitQueued(const qface_t *f)
{
	return R_EmitFace(f->fa, f->clipflags, f->key);
}

static void R_EmitFaceJob(SDL_UNUSED void *arg)
{
	edge_p = r_q.edge_p;
	surface_p = r_q.surface_p;
	r_polycount = r_q.polys;
	for(s32 t = 0;;){
		if(t == SDL_GetAtomicInt(&r_q.head)){
			if(SDL_GetAtomicInt(&r_q.stop) &&
					t == SDL_GetAtomicInt(&r_q.head))
				break;
			SDL_CPUPauseInstruction();
			continue;
		}
		SDL_MemoryBarrierAcquire();
		r_q.in[t & (R_EMITQ - 1)] =
			R_EmitQueued(&r_q.q[t & (R_EMITQ - 1)]);
		SDL_SetAtomicInt(&r_q.surfs, surface_p - surfaces);
		SDL_SetAtomicInt(&r_q.edges, edge_p - r_edges);
		SDL_MemoryBarrierRelease();
		SDL_SetAtomicInt(&r_q.tail, ++t);
	}
	r_q.edge_p = edge_p;
	r_q.surface_p = surface_p;
	r_q.polys = r_polycount;
}

// The walk can't go on without the job. If no worker has taken it, the walk
// emits the faces itself from now on.
static void R_QueueStall()
{
	if(!R_JobTakeBack(r_q.job))
		return; // it's running
	for(s32 t = SDL_GetAtomicInt(&r_q.tail); t != r_q.put; t++){
		r_q.in[t & (R_EMITQ - 1)] =
			R_EmitQueued(&r_q.q[t & (R_EMITQ - 1)]);
	}
	SDL_SetAtomicInt(&r_q.tail, r_q.put);
	r_q.queue = 0;
	r_q.job = -1;
}

static void R_QueuePublish()
{
	SDL_MemoryBarrierRelease();
	SDL_SetAtomicInt(&r_q.head, r_q.put);
}

// Did queued face i go in? Waits for the job to get there.
static bool R_QueuedIn(s32 i)
{
	R_QueuePublish();
	while(SDL_GetAtomicInt(&r_q.tail) <= i){
		R_QueueStall();
		SDL_CPUPauseInstruction();
	}
	SDL_MemoryBarrierAcquire();
	return r_q.in[i & (R_EMITQ - 1)];
}

// Is there sure to be room (R_FaceFits) for queued face i when the job gets
// to it? Every face before it takes a surface and numedges + 1 edges at most.
// Where the job has got is only looked up when that's not enough.
static bool R_QueueRoomAfter(s32 i, const msurface_t *fa)
{
	s64 surfs = r_q.seensurfs + (i - r_q.seen);
	s64 edges = r_q.seenedges + r_q.q[i & (R_EMITQ - 1)].edges
		- r_q.q[r_q.seen & (R_EMITQ - 1)].edges;
	return surfs < surf_max - surfaces
		&& edges + fa->numedges + 4 < edge_max - r_edges;
}

static bool R_QueueRoom(s32 i, const msurface_t *fa)
{
	if(R_QueueRoomAfter(i, fa))
		return 1;
	s32 t = SDL_GetAtomicInt(&r_q.tail);
	if(t > i)
		return 0; // R_QueuedIn knows
	r_q.seen = t; // the counts may be from after t, that's more
	r_q.seensurfs = SDL_GetAtomicInt(&r_q.surfs);
	r_q.seenedges = SDL_GetAtomicInt(&r_q.edges);
	return R_QueueRoomAfter(i, fa);
}

// Hands the face to the job and fills it in for the occlusion tests if it goes
// in. False if the walk has to emit it itself.
static bool R_QueueFace(msurface_t *fa, s32 clipflags, s32 key)
{
	s32 h = r_q.put;
	if(h - r_q.seen >= R_EMITQ){
		R_QueuePublish();
		while(h - (r_q.seen = SDL_GetAtomicInt(&r_q.tail)) >= R_EMITQ){
			R_QueueStall();
			if(!r_q.queue)
				return 0;
			SDL_CPUPauseInstruction();
		}
		r_q.seensurfs = SDL_GetAtomicInt(&r_q.surfs);
		r_q.seenedges = SDL_GetAtomicInt(&r_q.edges);
	}
	qface_t *f = &r_q.q[h & (R_EMITQ - 1)];
	f->fa = fa;
	f->clipflags = clipflags;
	f->key = key;
	f->edges = r_q.edgesum;
	r_q.edgesum += fa->numedges + 1;
	if(!(++r_q.put & (R_QBATCH - 1)))
		R_QueuePublish();
	if(!r_coarse_active || !R_CoarseOcclusionWillFill(fa))
		return 1;
	if(!(R_CoarseFaceSure(fa) && R_QueueRoom(h, fa)) && !R_QueuedIn(h))
		return 1;
	surf_t zs;
	R_FaceZi(fa, &zs);
	R_CoarseOcclusionAddSurface(fa, &zs);
	return 1;
}

// Around the world walk
void R_WalkBegin()
{
	r_q.walking = 1;
	for(s32 f = 0; f < 16; f++){
		clipplane_t *p = r_q.planes[f], *next = NULL;
		for(s32 i = 3; i >= 0; i--){
			p[i] = view_clipplanes[i];
			p[i].next = next;
			if(f & (1 << i))
				next = &p[i];
		}
		r_q.clip[f] = next;
	}
	if(!r_threaded)
		return;
	r_q.put = r_q.seen = 0;
	r_q.seensurfs = surface_p - surfaces;
	r_q.seenedges = edge_p - r_edges;
	r_q.edgesum = 0;
	r_q.edge_p = edge_p;
	r_q.surface_p = surface_p;
	r_q.polys = r_polycount;
	SDL_SetAtomicInt(&r_q.head, 0);
	SDL_SetAtomicInt(&r_q.tail, 0);
	SDL_SetAtomicInt(&r_q.stop, 0);
	SDL_SetAtomicInt(&r_q.surfs, r_q.seensurfs);
	SDL_SetAtomicInt(&r_q.edges, r_q.seenedges);
	r_q.job = R_JobTry(R_EmitFaceJob, NULL);
	r_q.queue = r_q.job >= 0;
	R_JobsPublish();
}

// The job's edge_p, surface_p and r_polycount are its thread's own: the walk
// takes them over
void R_WalkEnd()
{
	if(r_q.job >= 0){
		R_QueuePublish();
		SDL_SetAtomicInt(&r_q.stop, 1);
		R_JobWait(r_q.job); // runs it here if no worker took it
		r_q.job = -1;
		r_q.queue = 0;
		edge_p = r_q.edge_p;
		surface_p = r_q.surface_p;
		r_polycount = r_q.polys;
	}
	r_q.walking = 0;
}

void R_RenderFace(msurface_t *fa, s32 clipflags)
{
	// Manoel Kasimier - skyboxes - begin
	// Code taken from the ToChriS engine - Author: Vic
	// sky surfaces encountered in the world will cause the
	// environment box surfaces to be emited
	if((fa->flags & SURF_DRAWSKY) && skybox_name[0])
		r_skyframe = r_framecount;
	// Manoel Kasimier - skyboxes - end
	if(!R_SurfaceInPass(fa)){
		if(!r_alphapass && winquake_surface_liquid_alpha > 0 &&
				currententity == &cl_entities[0] &&
				!insubmodel)
			R_RecordAlphaFace(fa, clipflags);
		return;
	}
	c_faceclip++;
	if(r_coarse_active && R_CoarseOcclusionTestSurface(fa)){
		r_coarse_rejected++;
		return;
	}
	// The world walk's faces take their key whether they go in or not, so
	// the walk needn't wait for R_EmitFace
	s32 key = r_q.walking ? r_currentkey++ : -1;
	if(r_q.queue && R_QueueFace(fa, clipflags, key))
		return;
	if(R_EmitFace(fa, clipflags, key) && r_coarse_active)
		R_CoarseOcclusionAddSurface(fa, surface_p - 1);
}

void R_RenderBmodelFace(bedge_t *pedges, msurface_t *psurf)
{
	if(!R_SurfaceInPass(psurf) || !R_FaceFits(psurf->numedges))
		return;
	c_faceclip++;
	static medge_t tedge;
	// this is a dummy to give the caching mechanism someplace to write to
	r_pedge = &tedge;
	clipplane_t *pclip = R_ClipPlaneList(r_clipflags);
	r_emitted = 0; // push the edges through
	r_nearzi = 0;
	r_nearzionly = 0;
	makeleftedge = makerightedge = 0;
	r_lastvertvalid = 0;
	for(; pedges; pedges = pedges->pnext){
		r_leftclipped = r_rightclipped = 0;
		R_ClipEdge(pedges->v[0], pedges->v[1], pclip);
		if(r_leftclipped)
			makeleftedge = 1;
		if(r_rightclipped)
			makerightedge = 1;
	}
	// if there was a clip off the left edge, add that edge too
	if(makeleftedge){
		r_pedge = &tedge;
		R_ClipEdge(&r_leftexit, &r_leftenter, pclip->next);
	}
	// if there was a clip off the right edge, get the right r_nearzi
	if(makerightedge){
		r_pedge = &tedge;
		r_nearzionly = 1;
		R_ClipEdge(&r_rightexit, &r_rightenter,
				view_clipplanes[1].next);
	}
	// if no edges made it out, return without posting the surface
	if(!r_emitted)
		return;
	if((psurf->flags & SURF_DRAWSKY) && skybox_name[0])
		r_skyframe = r_framecount; // see R_RenderFace
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
	vec3_t p_normal;
	TransformVector(pplane->normal, p_normal);
	f32 distinv = 1.0 / (pplane->dist - R_PlaneDot(pplane));
	surface_p->d_zistepu = (xscaleinv * distinv) * p_normal[0];
	surface_p->d_zistepv = -((yscaleinv * distinv) * p_normal[1]);
	surface_p->d_ziorigin = p_normal[2] * distinv -
		xcenter * surface_p->d_zistepu - ycenter * surface_p->d_zistepv;
	surface_p++;
}
