// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.

#include "quakedef.h"

// surfaces are generated in back to front order by the bsp, so if a surf
// pointer is greater than another one, it should be drawn in front
// surfaces[1] is the background, and is used as the active surface stack

// Surfaces that are drawn as whole spans with no occlusion handling: they never
// enter the active surface stack, see R_GenerateSpans()
#define SURF_CUTMASK (SURF_DRAWCUTOUT | SURF_DRAWSKYBOX)
#define CUT_IDLE (-2) // left: not touched on this scanline, not in cut_list
#define CUT_NOLEAD (-1) // left: in cut_list, only trailing edges seen so far

static espan_t *span_p;
static s32 current_iv;
static s64 edge_head_u_shift20, edge_tail_u_shift20;
static edge_t edge_head, edge_tail;
static f32 fv;
static void(*pdrawfunc)(bool step);
static void R_GenerateSpans(bool step);
static void R_GenerateSpansBackward(bool step);
static void R_GenerateSpansAlpha(bool step);

// Cutout/skybox bookkeeping. For every such surface we collect the rightmost
// leading edge and the leftmost trailing edge while walking the active edge
// list once, then emit one span per surface at the end of the scanline. This
// replaces a rescan of the whole edge list for every cutout leading edge.
typedef struct { s32 left, right; } cutspan_t;
static cutspan_t *cut_acc; // indexed by surface number
static u32 *cut_list; // surfaces with any edge on this scanline
static s32 cut_count; // valid entries in cut_list
static s32 cut_surfs; // cutout/skybox surfaces in this frame
static s32 cut_right0; // reset value of cutspan_t.right

// In the alpha pass every surface but the background is drawn like a cutout
static bool R_IsCutSurf(const surf_t *s)
{ return r_alphapass ? s != &surfaces[1] : (s->flags & SURF_CUTMASK) != 0; }

// Decides whether this frame has any cutout/skybox surfaces at all and, if so,
// resets their accumulators. Must run after the surfaces have been emitted.
static void R_CutoutBeginFrame()
{
	cut_surfs = 0;
	cut_count = 0;
	for(surf_t *s = &surfaces[1]; s < surface_p; s++)
		if(R_IsCutSurf(s))
			cut_surfs++;
	if(!cut_surfs)
		return;
	s32 need = r_lsurfs_size + 2;
	static s32 cut_cap; // entries allocated in cut_acc and cut_list
	if(cut_cap < need){
		free(cut_acc);
		free(cut_list);
		cut_acc = malloc(need * sizeof(*cut_acc));
		cut_list = malloc(need * sizeof(*cut_list));
		if(!cut_acc || !cut_list)
			Sys_Error("R_CutoutBeginFrame: out of memory");
		cut_cap = need;
	}
	cut_right0 = (s32)r_refdef.vrectright;
	for(surf_t *s = &surfaces[1]; s < surface_p; s++){
		if(!R_IsCutSurf(s))
			continue;
		cut_acc[s - surfaces].left = CUT_IDLE;
		cut_acc[s - surfaces].right = cut_right0;
	}
}

void R_BeginEdgeFrame()
{
	r_edges = &ledges[0];
	// surface 0 doesn't really exist; it's just a dummy because
	// index 0 is used to indicate no edge attached to surface
	surfaces = &lsurfs[1]; // Point to the first "real" surface
	surf_max = &surfaces[r_lsurfs_size];
	edge_p = r_edges;
	edge_max = &r_edges[r_ledges_size];
	surface_p = &surfaces[2]; // background is surface 1,
	surfaces[1].spans = NULL; // no background spans yet
	surfaces[1].flags = SURF_DRAWBACKGROUND;
	if(r_draworder.value){ //put background behind everything in the world
		pdrawfunc = R_GenerateSpansBackward;
		surfaces[1].key = 0;
		r_currentkey = 1;
	}else{
		pdrawfunc = R_GenerateSpans;
		surfaces[1].key = 0x7FFFFFFF;
		r_currentkey = 0;
	}
	if(r_alphapass)
		pdrawfunc = R_GenerateSpansAlpha;
	memset(newedges + r_refdef.vrect.y, 0, (r_refdef.vrectbottom - 
				r_refdef.vrect.y) * sizeof(*newedges));
}

// R_EmitEdge just pushes the new edges of a row onto newedges[v]. This puts
// them in the order that inserting them one by one into a list sorted on u
// gives (each one before the first with u >= its own, plus one for trailing
// edges): on u, then leading edges before trailing ones, the leading ones
// newest first and the trailing ones oldest first. Every edge gets one key
// that sorts that way, with its index in r_edges in the low bits.
#define NE_IDXBITS 29 // u takes 34 bits at MAXWIDTH 15360
#define NE_IDXMASK (((u64)1 << NE_IDXBITS) - 1)

static void R_SortKeys(u64 *a, u64 *tmp, s64 n)
{
	const s64 run = 16;
	for(s64 r0 = 0; r0 < n; r0 += run){ // insertion sort short runs
		s64 r1 = r0 + run < n ? r0 + run : n;
		for(s64 i = r0 + 1; i < r1; i++){
			u64 k = a[i];
			s64 j = i;
			for(; j > r0 && a[j - 1] > k; j--)
				a[j] = a[j - 1];
			a[j] = k;
		}
	}
	u64 *src = a, *dst = tmp;
	for(s64 w = run; w < n; w *= 2){ // then merge them
		for(s64 l = 0; l < n; l += 2 * w){
			s64 m = l + w < n ? l + w : n;
			s64 r = l + 2 * w < n ? l + 2 * w : n;
			s64 i = l, j = m, o = l;
			while(i < m && j < r)
				dst[o++] = src[i] < src[j] ?
					src[i++] : src[j++];
			while(i < m) dst[o++] = src[i++];
			while(j < r) dst[o++] = src[j++];
		}
		u64 *t = src; src = dst; dst = t;
	}
	if(src != a)
		memcpy(a, src, n * sizeof(*a));
}

static edge_t *R_SortNewEdges(edge_t *list)
{
	if(!list->next)
		return list;
	s64 n = 0;
	static u64 *ne_keys; // ne_cap keys, then ne_cap more for merging
	static s64 ne_cap;
	for(edge_t *e = list; e; e = e->next){
		if(n == ne_cap){
			s64 cap = ne_cap ? ne_cap * 2 : 1024;
			u64 *k = malloc(2 * cap * sizeof(*k));
			if(!k)
				Sys_Error("R_SortNewEdges: out of memory");
			memcpy(k, ne_keys, n * sizeof(*k));
			free(ne_keys);
			ne_keys = k;
			ne_cap = cap;
		}
		u64 idx = (u64)(e - r_edges);
		u64 trail = e->surfs[0] != 0;
		ne_keys[n++] = (u64)e->u << (NE_IDXBITS + 1) |
			trail << NE_IDXBITS |
			(trail ? idx : NE_IDXMASK - idx);
	}
	R_SortKeys(ne_keys, ne_keys + ne_cap, n);
	edge_t *head = NULL, **link = &head;
	for(s64 i = 0; i < n; i++){
		u64 low = ne_keys[i] & NE_IDXMASK;
		bool trail = (ne_keys[i] >> NE_IDXBITS) & 1;
		edge_t *e = r_edges + (trail ? low : NE_IDXMASK - low);
		*link = e;
		link = &e->next;
	}
	*link = NULL;
	return head;
}

// Adds the edges in the linked list edgestoadd, adding them to the edges in the
// linked list edgelist. edgestoadd is assumed to be sorted on u, and non-empty
// (this is actually newedges[v]). edgelist is assumed to be sorted on u, with a
// sentinel at the end(actually, this is the active edge table starting at
// edge_head.next).
static void R_InsertNewEdges(edge_t *edgestoadd, edge_t *edgelist)
{
	edge_t *next_edge;
	do{
		next_edge = edgestoadd->next;
edgesearch:
		if(edgelist->u >= edgestoadd->u)
			goto addedge;
		edgelist = edgelist->next;
		if(edgelist->u >= edgestoadd->u)
			goto addedge;
		edgelist = edgelist->next;
		if(edgelist->u >= edgestoadd->u)
			goto addedge;
		edgelist = edgelist->next;
		if(edgelist->u >= edgestoadd->u)
			goto addedge;
		edgelist = edgelist->next;
		goto edgesearch;
addedge: // insert edgestoadd before edgelist
		edgestoadd->next = edgelist;
		edgestoadd->prev = edgelist->prev;
		edgelist->prev->next = edgestoadd;
		edgelist->prev = edgestoadd;
	}while((edgestoadd = next_edge) != NULL);
}

// Moves one active edge on to the next scanline, or takes it out of the list if
// this was its last scanline, and keeps the list sorted. Called from the span
// generation walk right after the edge's spans have been emitted, so the active
// edge list is only walked once per scanline. prev_u is the (already stepped) u
// of the edge in front of this one; next is e->next, read before any spans were
// emitted. Returns the prev_u to use for next.
static inline s64 R_StepEdge(edge_t *e, edge_t *next, s64 prev_u, s64 iv)
{
	if(e->vlast == iv){ // done with this edge
		next->prev = e->prev;
		e->prev->next = next;
		return prev_u;
	}
	s64 u = e->u + e->u_step;
	e->u = u;
	if(u >= prev_u)
		return u;
	// out of order: pull it out of the edge list and put it back where it
	// belongs. Everything in front of it has already been stepped, and the
	// one right in front is known to be too far right, so start behind it
	next->prev = e->prev;
	e->prev->next = next;
	edge_t *pwedge = e->prev->prev;
	while(pwedge->u > u)
		pwedge = pwedge->prev;
	e->next = pwedge->next;
	e->prev = pwedge;
	e->next->prev = e;
	pwedge->next = e;
	return next->prev->u; // the edge that used to be in front of e
}

static void R_CleanupSpan()
{
	// now that we've reached the right edge of the screen, we're done with
	// any unfinished surfaces, so emit a span for whatever's on top
	surf_t *surf = surfaces[1].next;
	s64 iu = edge_tail_u_shift20;
	if(iu > surf->last_u){
		espan_t *span = span_p++;
		span->u = surf->last_u;
		span->count = iu - span->u;
		span->v = current_iv;
		span->pnext = surf->spans;
		surf->spans = span;
	}
	do{ // reset spanstate for all surfaces in the surface stack
		surf->spanstate = 0;
		surf = surf->next;
	}while(surf != &surfaces[1]);
}

static inline void R_LeadingEdgeBackwards(const edge_t *edge, surf_t *surf)
{
	// it's adding a new surface in, so find the correct place
	// don't start a span if this is an inverted span, with the end edge
	// preceding the start edge(that is, we've already seen the end edge)
	surf_t *surf2;
	s32 iu; // keep here for the OpenBSD compiler
	if(++surf->spanstate == 1){
		surf2 = surfaces[1].next;
		if(surf->key > surf2->key)
			goto newtop;
		// if it's two surfaces on the same plane the one that's already
		// active is in front, so keep going unless it's a bmodel
		if(surf->insubmodel && (surf->key == surf2->key)){
			// must be two bmodels in the same leaf; don't care,
			// because they'll never be farthest anyway
			goto newtop;
		}
continue_search:
		do{
			surf2 = surf2->next;
		}while(surf->key < surf2->key);
		if(surf->key == surf2->key){
		//if it's two surfaces on the same plane, the one that's already
		//active is in front, so keep going unless it's a bmodel
			if(!surf->insubmodel)
				goto continue_search;
		// must be two bmodels in the same leaf; don't care which is
		// really in front, because they'll never be farthest anyway
		}
		goto gotposition;
newtop:
		iu = edge->u >> 20; // emit a span(obscures current top)
		if(iu > surf2->last_u){
			espan_t *span = span_p++;
			span->u = surf2->last_u;
			span->count = iu - span->u;
			span->v = current_iv;
			span->pnext = surf2->spans;
			surf2->spans = span;
		}
		surf->last_u = iu; // set last_u on the new span
gotposition:
		surf->next = surf2; // insert before surf2
		surf->prev = surf2->prev;
		surf2->prev->next = surf;
		surf2->prev = surf;
	}
}

static inline void R_TrailingEdge(surf_t *surf, const edge_t *edge)
{
	// don't generate a span if this is an inverted span, with the end edge
	// preceding the start edge(that is, we haven't seen the start edge yet)
	if(--surf->spanstate == 0){
		if(surf == surfaces[1].next){
			s64 iu = edge->u >> 20; // emit a span
			if(iu > surf->last_u){ // (current top going away)
				espan_t *span = span_p++;
				span->u = surf->last_u;
				span->count = iu - span->u;
				span->v = current_iv;
				span->pnext = surf->spans;
				surf->spans = span;
			} // set last_u on the surface below
			surf->next->last_u = iu;
		}
		surf->prev->next = surf->next;
		surf->next->prev = surf->prev;
	}
}

// Is surf, whose leading edge this is, in front of surf2 on the pixels it
// starts covering? They're two planes, so compare 1/z a pixel to the right of
// the edge: this gives the gap between parallel faces, and the slopes decide
// for faces that meet at the edge. WinQuake took faces within 1% of 1/z as
// meeting and went by the slopes alone, which for parallel faces only depends
// on which side they're seen from, so the far one could win.
static inline bool R_BmodelSurfInFront(const edge_t *edge, const surf_t *surf,
		const surf_t *surf2)
{
	f64 fu = (f32)(edge->u - 0xFFFFF) * (1.0 / 0x100000) + 1;
	f64 newzi = surf->d_ziorigin + fv*surf->d_zistepv + fu*surf->d_zistepu;
	f64 testzi=surf2->d_ziorigin+fv*surf2->d_zistepv+fu*surf2->d_zistepu;
	return newzi >= testzi;
}

// surf is &surfaces[edge->surfs[1]], which the caller has already checked is
// not the "no surface" index 0
static inline void R_LeadingEdge(const edge_t *edge, surf_t *surf)
{
	surf_t *surf2;
	s32 iu; // keep here for the OpenBSD compiler
	// it's adding a new surface in, so find the correct place
	// don't start a span if this is an inverted span, with the end
	// edge preceding the start edge(that is, we've already seen
	// the end edge)
	if(++surf->spanstate == 1){
		surf2 = surfaces[1].next;
		if(surf->key < surf2->key)
			goto newtop;
	//if it's two surfaces on the same plane, the one that's already
	//active is in front, so keep going unless it's a bmodel
		// must be two bmodels in the same leaf; sort on 1/z
		if(surf->insubmodel && (surf->key == surf2->key) &&
				R_BmodelSurfInFront(edge, surf, surf2))
			goto newtop;
		if(!surf->insubmodel){
		// The active surface list is sorted by key, with the background
		// at both ends. Hidden faces tend to have big keys, so when
		// ours is closer to the far end, search from there instead of
		// walking the whole stack. Same insertion point as the forward
		// search below: after every surface whose key is <= ours. (The
		// top surface is a real one here, since our key is not lower
		// than its key, so the far end is a real surface too.)
			surf_t *last = surfaces[1].prev;
			if((s64)surf->key-surf2->key>(s64)last->key-surf->key){
				while(last->key > surf->key)
					last = last->prev;
				surf2 = last->next;
				goto gotposition;
			}
		}
continue_search:
		do{
			surf2 = surf2->next;
		}while(surf->key > surf2->key);
		if(surf->key == surf2->key){
	//if it's two surfaces on the same plane, the one that's already
	//active is in front, so keep going unless it's a bmodel
			if(!surf->insubmodel)
				goto continue_search;
		// must be two bmodels in the same leaf; sort on 1/z
			if(R_BmodelSurfInFront(edge, surf, surf2))
				goto gotposition;
			goto continue_search;
		}
		goto gotposition;
newtop:
		// emit a span(obscures current top)
		iu = edge->u >> 20;
		if(iu > surf2->last_u){
			espan_t *span = span_p++;
			span->u = surf2->last_u;
			span->count = iu - span->u;
			span->v = current_iv;
			span->pnext = surf2->spans;
			surf2->spans = span;
		}
		// set last_u on the new span
		surf->last_u = iu;
gotposition:
		// insert before surf2
		surf->next = surf2;
		surf->prev = surf2->prev;
		surf2->prev->next = surf;
		surf2->prev = surf;
	}
}

static inline void R_CutTrailing(u32 sn, const edge_t *edge)
{
	cutspan_t *acc = &cut_acc[sn];
	s64 iu = edge->u >> 20;
	if(acc->left == CUT_IDLE){
		acc->left = CUT_NOLEAD;
		cut_list[cut_count++] = sn;
	}
	if(iu < acc->right)
		acc->right = iu;
}

static inline void R_CutLeading(u32 sn, const edge_t *edge)
{
	cutspan_t *acc = &cut_acc[sn];
	s64 iu = edge->u >> 20;
	if(iu < 0)
		iu = 0;
	if(acc->left < 0){ // first leading edge
		if(acc->left == CUT_IDLE)
			cut_list[cut_count++] = sn;
		acc->left = iu;
	}else if(iu > acc->left){
		acc->left = iu;
	}
}

// every edge has been seen by now. An inverted span(trailing before
// leading) gives right <= left and is skipped
static void R_EmitCutSpans()
{
	for(s32 i = 0; i < cut_count; i++){
		s32 sn = cut_list[i];
		cutspan_t *acc = &cut_acc[sn];
		// left < 0: only trailing edges seen, nothing to draw
		if(acc->left >= 0 && acc->right > acc->left){
			surf_t *surf = &surfaces[sn];
			espan_t *span = span_p++;
			span->u = acc->left;
			span->count = acc->right - acc->left;
			span->v = current_iv;
			span->pnext = surf->spans;
			surf->spans = span;
		}
		acc->left = CUT_IDLE;
		acc->right = cut_right0;
	}
	cut_count = 0;
}

static void R_GenerateSpans(bool step)
{
	surf_t *const sbase = surfaces; // the global is reloaded after every
				        // store through a pointer otherwise
	// clear active surfaces to just the background surface
	sbase[1].next = sbase[1].prev = &sbase[1];
	sbase[1].last_u = edge_head_u_shift20;
	if(!cut_surfs){ // no cutout or skybox surfaces in this frame
		s64 prev_u = edge_head.u;
		for(edge_t *edge=edge_head.next,*next; edge!=&edge_tail;
								edge=next){
			next = edge->next;
			// A shared(cached) edge can carry a trailing surface in
			// surfs[0] AND a leading one in surfs[1]. Always retire
			// the trailing surface first.
			u32 s0 = edge->surfs[0], s1 = edge->surfs[1];
			if(s0)
				R_TrailingEdge(&sbase[s0], edge);
			if(s1)
				R_LeadingEdge(edge, &sbase[s1]);
			if(step)
				prev_u = R_StepEdge(edge, next, prev_u,
							current_iv);
		}
		R_CleanupSpan();
		return;
	}
	// Cutout and skybox surfaces get a span covering the whole of the
	// surface on this scanline, from its rightmost leading edge to its
	// leftmost trailing edge, with intentional overdraw. They don't take
	// part in the active surface stack, so they only need those two
	// numbers, and we get them as the edges go by.
	s64 prev_u = edge_head.u;
	for(edge_t *edge=edge_head.next,*next; edge!=&edge_tail; edge=next){
		next = edge->next;
		u32 s0 = edge->surfs[0], s1 = edge->surfs[1];
		if(s0){
			surf_t *surf = &sbase[s0];
			if(surf->flags & SURF_CUTMASK)
				R_CutTrailing(s0, edge);
			else
				R_TrailingEdge(surf, edge);
		}
		if(s1){
			surf_t *surf = &sbase[s1];
			if(surf->flags & SURF_CUTMASK)
				R_CutLeading(s1, edge);
			else
				R_LeadingEdge(edge, surf);
		}
		if(step)
			prev_u = R_StepEdge(edge, next, prev_u, current_iv);
	}
	R_EmitCutSpans();
	R_CleanupSpan();
}

// The alpha pass: every surface is translucent, so none of them may hide
// another. They all get whole spans like cutouts, the z-buffer of the opaque
// pass does the occlusion and D_DrawSurfacesAlpha() draws them back to front.
static void R_GenerateSpansAlpha(bool step)
{
	s64 prev_u = edge_head.u;
	for(edge_t *edge=edge_head.next,*next; edge!=&edge_tail; edge=next){
		next = edge->next;
		u32 s0 = edge->surfs[0], s1 = edge->surfs[1];
		if(s0)
			R_CutTrailing(s0, edge);
		if(s1)
			R_CutLeading(s1, edge);
		if(step)
			prev_u = R_StepEdge(edge, next, prev_u, current_iv);
	}
	R_EmitCutSpans();
}

static void R_GenerateSpansBackward(bool step)
{
	surf_t *const sbase = surfaces;
	// clear active surfaces to just the background surface
	sbase[1].next = sbase[1].prev = &sbase[1];
	sbase[1].last_u = edge_head_u_shift20;
	// generate spans
	s64 prev_u = edge_head.u;
	for(edge_t *edge=edge_head.next,*next; edge!=&edge_tail; edge=next){
		next = edge->next;
		u32 s0 = edge->surfs[0], s1 = edge->surfs[1];
		if(s0)
			R_TrailingEdge(&sbase[s0], edge);
		if(s1)
			R_LeadingEdgeBackwards(edge, &sbase[s1]);
		if(step)
			prev_u = R_StepEdge(edge, next, prev_u, current_iv);
	}
	R_CleanupSpan();
}

static void R_FlushSpans(espan_t *basespan_p)
{
	if(r_drawflat.value) D_DrawSurfacesFlat();
	else if(!r_alphapass) D_DrawSurfaces();
	else D_DrawSurfacesAlpha();
	// clear the surface span pointers
	for(surf_t *s = &surfaces[1]; s < surface_p; s++)
		s->spans = NULL;
	span_p = basespan_p;
}

void R_ScanEdges()
{
	static u8 basespans[MAXSPANS * sizeof(espan_t)];
	espan_t *basespan_p = (espan_t *) basespans;
	R_CutoutBeginFrame();
	static espan_t *max_span_p;
	// Leave room for one full scanline of spans: at most one per pixel for
	// the visible surfaces, plus at most one per cutout/skybox surface
	max_span_p = &basespan_p[MAXSPANS - r_refdef.vrect.width - cut_surfs];
	span_p = basespan_p;
	// clear active edges to just background edges around the whole screen
	edge_head.u = (s64)r_refdef.vrect.x << 20;
	edge_head_u_shift20 = edge_head.u >> 20;
	edge_head.u_step = 0;
	edge_head.prev = NULL;
	edge_head.next = &edge_tail;
	edge_head.surfs[0] = 0;
	edge_head.surfs[1] = 1;
	edge_tail.u = ((s64)r_refdef.vrectright << 20) + 0xFFFFF;
	edge_tail_u_shift20 = edge_tail.u >> 20;
	edge_tail.u_step = 0;
	edge_tail.prev = &edge_head;
	static edge_t edge_aftertail;
	edge_tail.next = &edge_aftertail;
	edge_tail.surfs[0] = 1;
	edge_tail.surfs[1] = 0;
	edge_aftertail.u = -1; // force a move
	edge_aftertail.u_step = 0;
	static edge_t edge_sentinel;
	edge_aftertail.next = &edge_sentinel;
	edge_aftertail.prev = &edge_tail;
	edge_sentinel.u = ((u64)2000) << 24;// make sure nothing sorts past this
	edge_sentinel.prev = &edge_aftertail;
	s64 bottom = r_refdef.vrectbottom - 1; // process all scan lines
	s64 iv = r_refdef.vrect.y;
	for(; iv < bottom; iv++){
		current_iv = iv;
		fv = (f32)iv;
		// mark that the head(background start) span is pre-included
		surfaces[1].spanstate = 1;
		if(newedges[iv])
			R_InsertNewEdges(R_SortNewEdges(newedges[iv]),
					edge_head.next);
		(*pdrawfunc)(1); // also steps/removes the edges
		// flush the span list if we can't be sure we have enough spans
		// left for the next scan
		if(span_p >= max_span_p)
			R_FlushSpans(basespan_p);
	}
	// do the last scan(no need to step or sort or remove on the last scan)
	current_iv = iv;
	fv = (f32)iv;
	// mark that the head(background start) span is pre-included
	surfaces[1].spanstate = 1;
	if(newedges[iv])
		R_InsertNewEdges(R_SortNewEdges(newedges[iv]), edge_head.next);
	(*pdrawfunc)(0);
	if(r_drawflat.value) D_DrawSurfacesFlat();
	else if(!r_alphapass) D_DrawSurfaces();
	else D_DrawSurfacesAlpha();
}
