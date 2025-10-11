// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.

#include "quakedef.h"

// FIXME the complex cases add new polys on most lines,
// so dont optimize for keeping them the same have multiple free span lists to
// try to get better coherence ? low depth complexity-- 1 to 3 or so this
// breaks spans at every edge, even hidden ones(bad)
// have a sentinal at both ends ?

// surfaces are generated in back to front order by the bsp, so if a surf
// pointer is greater than another one, it should be drawn in front
// surfaces[1] is the background, and is used as the active surface stack

static espan_t *span_p, *max_span_p;
static edge_t edge_sentinel;
static s32 current_iv;
static s64 edge_head_u_shift20, edge_tail_u_shift20;
static edge_t edge_head, edge_tail;
static edge_t edge_aftertail;
static f32 fv;
static void (*pdrawfunc)();
static edge_t ledges[NUMSTACKEDGES + ((CACHE_SIZE - 1) / sizeof(edge_t)) + 1];
static surf_t lsurfs[NUMSTACKSURFACES + ((CACHE_SIZE - 1) / sizeof(surf_t)) +1];
void R_GenerateSpans();
void R_GenerateSpansBackward();
void R_GenerateSpansCutout();
void R_LeadingEdge(edge_t *edge);
void R_LeadingEdgeBackwards(edge_t *edge);
void R_TrailingEdge(surf_t *surf, edge_t *edge);

void R_BeginEdgeFrame()
{
	r_edges = &ledges[0];
	// surface 0 doesn't really exist; it's just a dummy because
	// index 0 is used to indicate no edge attached to surface
	surfaces = &lsurfs[1]; // Point to the first "real" surface
	surf_max = &surfaces[NUMSTACKSURFACES];
	edge_p = r_edges;
	edge_max = &r_edges[NUMSTACKEDGES];
	surface_p = &surfaces[2]; // background is surface 1,
	surfaces[1].spans = NULL; // no background spans yet
	surfaces[1].flags = SURF_DRAWBACKGROUND;
	if (r_draworder.value) { //put background behind everything in the world
		pdrawfunc = R_GenerateSpansBackward;
		surfaces[1].key = 0;
		r_currentkey = 1;
	} else {
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
// sentinel at the end (actually, this is the active edge table starting at
// edge_head.next).
void R_InsertNewEdges(edge_t *edgestoadd, edge_t *edgelist)
{
	edge_t *next_edge;
	do {
		next_edge = edgestoadd->next;
edgesearch:
		if (edgelist->u >= edgestoadd->u)
			goto addedge;
		edgelist = edgelist->next;
		if (edgelist->u >= edgestoadd->u)
			goto addedge;
		edgelist = edgelist->next;
		if (edgelist->u >= edgestoadd->u)
			goto addedge;
		edgelist = edgelist->next;
		if (edgelist->u >= edgestoadd->u)
			goto addedge;
		edgelist = edgelist->next;
		goto edgesearch;
addedge: // insert edgestoadd before edgelist
		edgestoadd->next = edgelist;
		edgestoadd->prev = edgelist->prev;
		edgelist->prev->next = edgestoadd;
		edgelist->prev = edgestoadd;
	} while ((edgestoadd = next_edge) != NULL);
}

void R_RemoveEdges(edge_t *pedge)
{
	do {
		pedge->next->prev = pedge->prev;
		pedge->prev->next = pedge->next;
	} while ((pedge = pedge->nextremove) != NULL);
}

void R_StepActiveU(edge_t *pedge)
{
	while (1) {
nextedge:
		pedge->u += pedge->u_step;
		if (pedge->u < pedge->prev->u)
			goto pushback;
		pedge = pedge->next;
		pedge->u += pedge->u_step;
		if (pedge->u < pedge->prev->u)
			goto pushback;
		pedge = pedge->next;
		pedge->u += pedge->u_step;
		if (pedge->u < pedge->prev->u)
			goto pushback;
		pedge = pedge->next;
		pedge->u += pedge->u_step;
		if (pedge->u < pedge->prev->u)
			goto pushback;
		pedge = pedge->next;
		goto nextedge;
pushback:
		if (pedge == &edge_aftertail)
			return;
		// push it back to keep it sorted               
		edge_t *pnext_edge = pedge->next;
		// pull the edge out of the edge list
		pedge->next->prev = pedge->prev;
		pedge->prev->next = pedge->next;
		// find out where the edge goes in the edge list
		edge_t *pwedge = pedge->prev->prev;

		while (pwedge->u > pedge->u) {
			pwedge = pwedge->prev;
		}
		// put the edge back into the edge list
		pedge->next = pwedge->next;
		pedge->prev = pwedge;
		pedge->next->prev = pedge;
		pwedge->next = pedge;
		pedge = pnext_edge;
		if (pedge == &edge_tail)
			return;
	}
}

void R_CleanupSpan()
{
	// now that we've reached the right edge of the screen, we're done with
	// any unfinished surfaces, so emit a span for whatever's on top
	surf_t *surf = surfaces[1].next;
	s64 iu = edge_tail_u_shift20;
	if (iu > surf->last_u) {
		espan_t *span = span_p++;
		span->u = surf->last_u;
		span->count = iu - span->u;
		span->v = current_iv;
		span->pnext = surf->spans;
		surf->spans = span;
	}
	do { // reset spanstate for all surfaces in the surface stack
		surf->spanstate = 0;
		surf = surf->next;
	} while (surf != &surfaces[1]);
}

void R_LeadingEdgeBackwards(edge_t *edge)
{
	// it's adding a new surface in, so find the correct place
	surf_t *surf = &surfaces[edge->surfs[1]];
	// don't start a span if this is an inverted span, with the end edge
	// preceding the start edge (that is, we've already seen the end edge)
	surf_t *surf2;
	s32 iu; // keep here for the OpenBSD compiler
	if (++surf->spanstate == 1) {
		surf2 = surfaces[1].next;
		if (surf->key > surf2->key)
			goto newtop;
		// if it's two surfaces on the same plane the one that's already
		// active is in front, so keep going unless it's a bmodel
		if (surf->insubmodel && (surf->key == surf2->key)) {
			// must be two bmodels in the same leaf; don't care,
			// because they'll never be farthest anyway
			goto newtop;
		}
continue_search:
		do {
			surf2 = surf2->next;
		} while (surf->key < surf2->key);
		if (surf->key == surf2->key) {
			// if it's two surfaces on the same plane, the one that's already
			// active is in front, so keep going unless it's a bmodel
			if (!surf->insubmodel)
				goto continue_search;
			// must be two bmodels in the same leaf; don't care which is really
			// in front, because they'll never be farthest anyway
		}
		goto gotposition;
newtop:
		iu = edge->u >> 20; // emit a span (obscures current top)
		if (iu > surf2->last_u) {
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

void R_TrailingEdge(surf_t *surf, edge_t *edge)
{
	// don't generate a span if this is an inverted span, with the end edge
	// preceding the start edge(that is, we haven't seen the start edge yet)
	if (--surf->spanstate == 0) {
		if (surf == surfaces[1].next) {
			s64 iu = edge->u >> 20; // emit a span
			if (iu > surf->last_u) { // (current top going away)
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

void R_LeadingEdge (edge_t *edge)
{
	surf_t *surf, *surf2;
	s32 iu; // keep here for the OpenBSD compiler
	if (edge->surfs[1])
	{
		// it's adding a new surface in, so find the correct place
		surf = &surfaces[edge->surfs[1]];
		// don't start a span if this is an inverted span, with the end
		// edge preceding the start edge (that is, we've already seen the
		// end edge)
		if (++surf->spanstate == 1)
		{
			surf2 = surfaces[1].next;
			if (surf->key < surf2->key)
				goto newtop;
			// if it's two surfaces on the same plane, the one that's already
			// active is in front, so keep going unless it's a bmodel
			if (surf->insubmodel && (surf->key == surf2->key))
			{
				// must be two bmodels in the same leaf; sort on 1/z
				f64 fu = (f32)(edge->u-0xFFFFF) * (1.0/0x100000);
				f64 newzi = surf->d_ziorigin +
					fv*surf->d_zistepv + fu*surf->d_zistepu;
				f64 newzibottom = newzi * 0.99;
				f64 testzi = surf2->d_ziorigin +
					fv*surf2->d_zistepv+fu*surf2->d_zistepu;
				if (newzibottom >= testzi)
					goto newtop;
				f64 newzitop = newzi * 1.01;
				if (newzitop >= testzi &&
					surf->d_zistepu >= surf2->d_zistepu)
						goto newtop;
			}
continue_search:
			do
			{
				surf2 = surf2->next;
			} while (surf->key > surf2->key);
			if (surf->key == surf2->key)
			{
				// if it's two surfaces on the same plane, the one that's already
				// active is in front, so keep going unless it's a bmodel
				if (!surf->insubmodel)
					goto continue_search;
				// must be two bmodels in the same leaf; sort on 1/z
				f64 fu = (f32)(edge->u-0xFFFFF) * (1.0/0x100000);
				f64 newzi = surf->d_ziorigin +
					fv*surf->d_zistepv + fu*surf->d_zistepu;
				f64 newzibottom = newzi * 0.99;
				f64 testzi = surf2->d_ziorigin +
					fv*surf2->d_zistepv+fu*surf2->d_zistepu;
				if (newzibottom >= testzi)
					goto gotposition;
				f64 newzitop = newzi * 1.01;
				if (newzitop >= testzi &&
					surf->d_zistepu >= surf2->d_zistepu)
						goto gotposition;
				goto continue_search;
			}
			goto gotposition;
newtop:
			// emit a span (obscures current top)
			iu = edge->u >> 20;
			if (iu > surf2->last_u)
			{
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
}

void R_GenerateSpans()
{
	// clear active surfaces to just the background surface
	surfaces[1].next = surfaces[1].prev = &surfaces[1];
	surfaces[1].last_u = edge_head_u_shift20;
	// generate spans
	for (edge_t *edge=edge_head.next; edge!=&edge_tail; edge=edge->next) {
		if (surfaces[edge->surfs[1]].flags & SURF_DRAWCUTOUT) {
		// Generates spans for whole cutout surfaces, from the leftmost
		// point of each surface to the rightmost, intentional overdraw
			s32 surfn = edge->surfs[1];
			s64 left = 0;
			s64 right = r_refdef.vrectright;
			for (edge_t *e2=edge_head.next; e2!=&edge_tail; e2=e2->next) {
				if (e2->surfs[1] == surfn && (e2->u>>20) > left)
					left = (e2->u>>20);
				if (e2->surfs[0] == surfn && (e2->u>>20) < right)
					right = (e2->u>>20);
			}
			surf_t *surf = &surfaces[surfn];
			espan_t *span = span_p++;
			span->u = left;
			span->count = right - left;
			span->v = current_iv;
			span->pnext = surf->spans;
			surf->spans = span;
		} else { // Normal span generation that avoids overdraw
			if (edge->surfs[0]) {
				// it has a left surface, so a surface is going away for this span
				surf_t *surf = &surfaces[edge->surfs[0]];
				R_TrailingEdge(surf, edge);
				if (!edge->surfs[1])
					continue;
			}
			R_LeadingEdge(edge);
		}
	}
	R_CleanupSpan();
	return;
}

void R_GenerateSpansBackward()
{
	// clear active surfaces to just the background surface
	surfaces[1].next = surfaces[1].prev = &surfaces[1];
	surfaces[1].last_u = edge_head_u_shift20;
	// generate spans
	for (edge_t *edge=edge_head.next; edge!=&edge_tail; edge=edge->next) {
		if (edge->surfs[0])
			R_TrailingEdge(&surfaces[edge->surfs[0]], edge);
		if (edge->surfs[1])
			R_LeadingEdgeBackwards(edge);
	}
	R_CleanupSpan();
}

// Input: newedges[] array this has links to edges, which have links to surfaces
// Output: Each surface has a linked list of its visible spans
void R_ScanEdges()
{
	// Align the array itself to the cache size
	u8 basespans[MAXSPANS * sizeof(espan_t) + CACHE_SIZE];
	// Pointer to the aligned base of the spans
	espan_t *basespan_p = (espan_t *) basespans;
	// No more pointer adjustment needed because the array is aligned
	max_span_p = &basespan_p[MAXSPANS - r_refdef.vrect.width];
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
	// FIXME: do we need this now that we clamp x in r_draw.c?
	edge_sentinel.u = ((u64)2000) << 24; // make sure nothing sorts past this
	edge_sentinel.prev = &edge_aftertail;
	s64 bottom = r_refdef.vrectbottom - 1; // process all scan lines
	s64 iv = r_refdef.vrect.y;
	for (; iv < bottom; iv++) {
		current_iv = iv;
		fv = (f32)iv;
		// mark that the head (background start) span is pre-included
		surfaces[1].spanstate = 1;
		if (newedges[iv])
			R_InsertNewEdges(newedges[iv], edge_head.next);
		(*pdrawfunc) ();
		// flush the span list if we can't be sure we have enough spans
		// left for the next scan
		if (span_p >= max_span_p) {
			if (r_drawflat.value) D_DrawSurfacesFlat();
			else if((s32)r_twopass.value&1) {
				if (r_pass + r_wateralphapass == 0) D_DrawSurfacesPass1();
				else if (r_pass + r_wateralphapass == 2) D_DrawSurfacesPass3();
			}
			else D_DrawSurfacesPass1();
			// clear the surface span pointers
			for (surf_t *s = &surfaces[1]; s < surface_p; s++)
				s->spans = NULL;
			span_p = basespan_p;
		}
		if (removeedges[iv])
			R_RemoveEdges(removeedges[iv]);
		if (edge_head.next != &edge_tail)
			R_StepActiveU(edge_head.next);
	}
	// do the last scan (no need to step or sort or remove on the last scan)
	current_iv = iv;
	fv = (f32)iv;
	// mark that the head (background start) span is pre-included
	surfaces[1].spanstate = 1;
	if (newedges[iv])
		R_InsertNewEdges(newedges[iv], edge_head.next);
	(*pdrawfunc) ();
	if (r_drawflat.value) D_DrawSurfacesFlat();
	else if((s32)r_twopass.value&1) {
		if (r_pass + r_wateralphapass == 0) D_DrawSurfacesPass1();
		else if (r_pass + r_wateralphapass == 2) D_DrawSurfacesPass3();
	}
	else D_DrawSurfacesPass1();
}
