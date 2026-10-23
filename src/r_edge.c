// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.

#include "quakedef.h"

// FIXME the complex cases add new polys on most lines,
// so dont optimize for keeping them the same have multiple free span lists to
// try to get better coherence ? low depth complexity-- 1 to 3 or so this
// breaks spans at every edge, even hidden ones(bad)
// have a sentinal at both ends ?

// surfaces are generated in back to front order by the bsp, so if a surf
// pointer is greater than another one, it should be drawn in front
// surfaces[1] is the background, and is used as the active surface stack

// Surfaces that are drawn as whole spans with no occlusion handling: they never
// enter the active surface stack, see R_GenerateSpans()
#define SURF_CUTMASK (SURF_DRAWCUTOUT | SURF_DRAWSKYBOX)
#define CUT_IDLE (-2) // left: not touched on this scanline, not in cut_list
#define CUT_NOLEAD (-1) // left: in cut_list, only trailing edges seen so far

static espan_t *span_p, *max_span_p;
static edge_t edge_sentinel;
static s32 current_iv;
static s64 edge_head_u_shift20, edge_tail_u_shift20;
static edge_t edge_head, edge_tail;
static edge_t edge_aftertail;
static f32 fv;
static void(*pdrawfunc)();
static void R_GenerateSpans();
static void R_GenerateSpansBackward();

// Cutout/skybox bookkeeping. For every such surface we collect the rightmost
// leading edge and the leftmost trailing edge while walking the active edge
// list once, then emit one span per surface at the end of the scanline. This
// replaces a rescan of the whole edge list for every cutout leading edge.
typedef struct { s32 left, right; } cutspan_t;
static cutspan_t *cut_acc; // indexed by surface number
static u16 *cut_list; // surfaces with any edge on this scanline
static s32 cut_cap; // entries allocated in cut_acc and cut_list
static s32 cut_count; // valid entries in cut_list
static s32 cut_surfs; // cutout/skybox surfaces in this frame
static s32 cut_right0; // reset value of cutspan_t.right

// Decides whether this frame has any cutout/skybox surfaces at all and, if so,
// resets their accumulators. Must run after the surfaces have been emitted.
static void R_CutoutBeginFrame()
{
	cut_surfs = 0;
	cut_count = 0;
	for(surf_t *s = &surfaces[1]; s < surface_p; s++)
		if(s->flags & SURF_CUTMASK)
			cut_surfs++;
	if(!cut_surfs)
		return;
	s32 need = r_lsurfs_size + 2;
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
		if(!(s->flags & SURF_CUTMASK))
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
	memset(newedges + r_refdef.vrect.y, 0, (r_refdef.vrectbottom - 
				r_refdef.vrect.y) * sizeof(*newedges));
	memset(removeedges + r_refdef.vrect.y, 0, (r_refdef.vrectbottom - 
				r_refdef.vrect.y) * sizeof(*removeedges));
	memset(last_pcheck, 0, sizeof(edge_t*)*scr_vrect.height);
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

static void R_RemoveEdges(edge_t *pedge)
{
	do{
		pedge->next->prev = pedge->prev;
		pedge->prev->next = pedge->next;
	}while((pedge = pedge->nextremove) != NULL);
}

// Steps every active edge to the next scanline and keeps the list sorted.
// The u of the previous node is carried in a register instead of being
// reloaded through pedge->prev for every edge.
static void R_StepActiveU(edge_t *pedge)
{
	s64 u, prev_u = pedge->prev->u;
	while(1){
nextedge:
		u = pedge->u + pedge->u_step;
		pedge->u = u;
		if(u < prev_u)
			goto pushback;
		prev_u = u;
		pedge = pedge->next;
		u = pedge->u + pedge->u_step;
		pedge->u = u;
		if(u < prev_u)
			goto pushback;
		prev_u = u;
		pedge = pedge->next;
		u = pedge->u + pedge->u_step;
		pedge->u = u;
		if(u < prev_u)
			goto pushback;
		prev_u = u;
		pedge = pedge->next;
		u = pedge->u + pedge->u_step;
		pedge->u = u;
		if(u < prev_u)
			goto pushback;
		prev_u = u;
		pedge = pedge->next;
		goto nextedge;
pushback:
		if(pedge == &edge_aftertail)
			return;
		// push it back to keep it sorted               
		edge_t *pnext_edge = pedge->next;
		// pull the edge out of the edge list
		pedge->next->prev = pedge->prev;
		pedge->prev->next = pedge->next;
		// find out where the edge goes in the edge list
		edge_t *pwedge = pedge->prev->prev;

		while(pwedge->u > u){
			pwedge = pwedge->prev;
		}
		// put the edge back into the edge list
		pedge->next = pwedge->next;
		pedge->prev = pwedge;
		pedge->next->prev = pedge;
		pwedge->next = pedge;
		pedge = pnext_edge;
		if(pedge == &edge_tail)
			return;
		// pnext_edge's predecessor is now the node that preceded the
		// one we just moved, which has already been stepped
		prev_u = pedge->prev->u;
	}
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
		if(surf->insubmodel && (surf->key == surf2->key)){
		// must be two bmodels in the same leaf; sort on 1/z
			f64 fu = (f32)(edge->u-0xFFFFF) * (1.0/0x100000);
			f64 newzi = surf->d_ziorigin +
				fv*surf->d_zistepv + fu*surf->d_zistepu;
			f64 newzibottom = newzi * 0.99;
			f64 testzi = surf2->d_ziorigin +
				fv*surf2->d_zistepv+fu*surf2->d_zistepu;
			if(newzibottom >= testzi)
				goto newtop;
			f64 newzitop = newzi * 1.01;
			if(newzitop >= testzi &&
				surf->d_zistepu >= surf2->d_zistepu)
					goto newtop;
		}
continue_search:
		do
		{
			surf2 = surf2->next;
		}while(surf->key > surf2->key);
		if(surf->key == surf2->key){
	//if it's two surfaces on the same plane, the one that's already
	//active is in front, so keep going unless it's a bmodel
			if(!surf->insubmodel)
				goto continue_search;
		// must be two bmodels in the same leaf; sort on 1/z
			f64 fu = (f32)(edge->u-0xFFFFF)*(1.0/0x100000);
			f64 newzi = surf->d_ziorigin +
				fv*surf->d_zistepv + fu*surf->d_zistepu;
			f64 newzibottom = newzi * 0.99;
			f64 testzi = surf2->d_ziorigin +
				fv*surf2->d_zistepv+fu*surf2->d_zistepu;
			if(newzibottom >= testzi)
				goto gotposition;
			f64 newzitop = newzi * 1.01;
			if(newzitop >= testzi &&
				surf->d_zistepu >= surf2->d_zistepu)
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

static void R_GenerateSpans()
{
	surf_t *const sbase = surfaces; // the global is reloaded after every
				        // store through a pointer otherwise
	// clear active surfaces to just the background surface
	sbase[1].next = sbase[1].prev = &sbase[1];
	sbase[1].last_u = edge_head_u_shift20;
	if(!cut_surfs){ // no cutout or skybox surfaces in this frame
		for(edge_t *edge=edge_head.next; edge!=&edge_tail;
							edge=edge->next){
			// A shared(cached) edge can carry a trailing surface in
			// surfs[0] AND a leading one in surfs[1]. Always retire
			// the trailing surface first.
			u16 s0 = edge->surfs[0], s1 = edge->surfs[1];
			if(s0)
				R_TrailingEdge(&sbase[s0], edge);
			if(s1)
				R_LeadingEdge(edge, &sbase[s1]);
		}
		R_CleanupSpan();
		return;
	}
	// Cutout and skybox surfaces get a span covering the whole of the
	// surface on this scanline, from its rightmost leading edge to its
	// leftmost trailing edge, with intentional overdraw. They don't take
	// part in the active surface stack, so they only need those two
	// numbers, and we get them as the edges go by.
	for(edge_t *edge=edge_head.next; edge!=&edge_tail; edge=edge->next){
		u16 s0 = edge->surfs[0], s1 = edge->surfs[1];
		if(s0){
			surf_t *surf = &sbase[s0];
			if(surf->flags & SURF_CUTMASK){
				cutspan_t *acc = &cut_acc[s0];
				s64 iu = edge->u >> 20;
				if(acc->left == CUT_IDLE){
					acc->left = CUT_NOLEAD;
					cut_list[cut_count++] = s0;
				}
				if(iu < acc->right)
					acc->right = iu;
			}else{
				R_TrailingEdge(surf, edge);
			}
		}
		if(s1){
			surf_t *surf = &sbase[s1];
			if(surf->flags & SURF_CUTMASK){
				cutspan_t *acc = &cut_acc[s1];
				s64 iu = edge->u >> 20;
				if(iu < 0)
					iu = 0;
				if(acc->left < 0){ // first leading edge
					if(acc->left == CUT_IDLE)
						cut_list[cut_count++] = s1;
					acc->left = iu;
				}else if(iu > acc->left){
					acc->left = iu;
				}
			}else{
				R_LeadingEdge(edge, surf);
			}
		}
	}
	// every edge has been seen by now. An inverted span(trailing before
	// leading) gives right <= left and is skipped
	for(s32 i = 0; i < cut_count; i++){
		s32 sn = cut_list[i];
		cutspan_t *acc = &cut_acc[sn];
		// left < 0: only trailing edges seen, nothing to draw
		if(acc->left >= 0 && acc->right > acc->left){
			surf_t *surf = &sbase[sn];
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
	R_CleanupSpan();
}

static void R_GenerateSpansBackward()
{
	surf_t *const sbase = surfaces;
	// clear active surfaces to just the background surface
	sbase[1].next = sbase[1].prev = &sbase[1];
	sbase[1].last_u = edge_head_u_shift20;
	// generate spans
	for(edge_t *edge=edge_head.next; edge!=&edge_tail; edge=edge->next){
		u16 s0 = edge->surfs[0], s1 = edge->surfs[1];
		if(s0)
			R_TrailingEdge(&sbase[s0], edge);
		if(s1)
			R_LeadingEdgeBackwards(edge, &sbase[s1]);
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
	// Aligned to the cache size. Static because it is ~800KB, too much to
	// keep on the stack(and it lets MAXSPANS be raised safely); this
	// function is not reentrant anyway.
	static _Alignas(CACHE_SIZE) u8 basespans[MAXSPANS * sizeof(espan_t)];
	espan_t *basespan_p = (espan_t *) basespans;
	R_CutoutBeginFrame();
	// Leave room for one full scanline of spans: at most one per pixel for
	// the visible surfaces, plus at most one per cutout/skybox surface
	max_span_p = &basespan_p[MAXSPANS - r_refdef.vrect.width - cut_surfs];
	span_p = basespan_p;
	// clear active edges to just background edges around the whole screen
	// FIXME: most of this only needs to be set up once
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
	edge_tail.next = &edge_aftertail;
	edge_tail.surfs[0] = 1;
	edge_tail.surfs[1] = 0;
	edge_aftertail.u = -1; // force a move
	edge_aftertail.u_step = 0;
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
			R_InsertNewEdges(newedges[iv], edge_head.next);
		(*pdrawfunc) ();
		// flush the span list if we can't be sure we have enough spans
		// left for the next scan
		if(span_p >= max_span_p)
			R_FlushSpans(basespan_p);
		if(removeedges[iv])
			R_RemoveEdges(removeedges[iv]);
		if(edge_head.next != &edge_tail)
			R_StepActiveU(edge_head.next);
	}
	// do the last scan(no need to step or sort or remove on the last scan)
	current_iv = iv;
	fv = (f32)iv;
	// mark that the head(background start) span is pre-included
	surfaces[1].spanstate = 1;
	if(newedges[iv])
		R_InsertNewEdges(newedges[iv], edge_head.next);
	(*pdrawfunc)();
	if(r_drawflat.value) D_DrawSurfacesFlat();
	else if(!r_alphapass) D_DrawSurfaces();
	else D_DrawSurfacesAlpha();
}
