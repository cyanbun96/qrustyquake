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

// Cutout/skybox bookkeeping. For every such surface we collect the rightmost
// leading edge and the leftmost trailing edge while walking the active edge
// list once, then emit one span per surface at the end of the scanline. This
// replaces a rescan of the whole edge list for every cutout leading edge.
typedef struct { s32 left, right; } cutspan_t;
// The edges that only cutout/skybox surfaces have. Those surfaces only need
// the rightmost leading and the leftmost trailing edge, in whatever order, so
// the edges stay out of the sorted active edge list, where they'd be in the
// way.
typedef struct { s64 u, u_step; u32 surfs[2]; s32 vlast; } cutedge_t;

// What a scan of rows works with: all of them with the surfaces themselves on
// the main thread, or a band of them with its own copies, see R_ScanBand
typedef struct {
	espan_t *span_p;
	s32 iv;
	f32 fv;
	surf_t *sbase; // by number, [1] is the background
	edge_t head, tail, aftertail, sentinel;
	cutspan_t *cut_acc; // indexed by surface number
	u32 *cut_list; // surfaces with any edge on this scanline
	s32 cut_count; // valid entries in cut_list
	s32 cut_cap; // entries allocated in cut_acc and cut_list
	cutedge_t *cut_edges;
	s32 cut_nedges;
	s64 cut_edgecap;
	u64 *ne_keys; // R_SortNewEdges: ne_cap keys, then ne_cap for merging
	s64 ne_cap;
	u32 *touched; // a band's surfaces with spans since its last draw band
	s32 ntouched;
} scan_t;

static scan_t r_scanmain;
// With the workers, a bit per surface that got spans since the last flush, so
// that a flush for a band doesn't have to go through all the surfaces
u64 *r_surfmarks;
u64 *r_flushmarks; // those of all the bands since the surfaces were set up
struct dsurf_s **r_surfdsurf; // how each one's spans get drawn,
SDL_AtomicInt *r_surfscan; // if set up in this scan
u32 r_scan;
static s64 edge_head_u_shift20, edge_tail_u_shift20;
static void(*pdrawfunc)(bool step);
static void R_GenerateSpans(bool step);
static void R_GenerateSpansBackward(bool step);
static void R_GenerateSpansAlpha(bool step);
static void R_GenerateSpansMarked(bool step);
static void R_GenerateSpansBackwardMarked(bool step);
static void R_GenerateSpansAlphaMarked(bool step);
static s32 cut_surfs; // cutout/skybox surfaces in this frame
static s32 cut_right0; // reset value of cutspan_t.right

// In the alpha pass every surface but the background is drawn like a cutout
static bool R_IsCutSurf(const surf_t *s)
{ return r_alphapass ? s != &surfaces[1] : (s->flags & SURF_CUTMASK) != 0; }

// p with room for need of size each, which it keeps
static void *R_Grow(void *p, s64 *cap, s64 need, size_t size)
{
	if(*cap >= need)
		return p;
	s64 n = *cap ? *cap : 256;
	while(n < need)
		n *= 2;
	if(!(p = realloc(p, n * size)))
		Sys_Error("R_Grow: out of memory");
	*cap = n;
	return p;
}

// The surface marks and setups the workers' bands go by, for a scan
static void R_MarksBeginScan()
{
	s32 need = r_lsurfs_size + 2;
	static s32 cap;
	if(cap < need){
		free(r_surfmarks);
		free(r_flushmarks);
		free(r_surfdsurf);
		free(r_surfscan);
		r_surfmarks = calloc((need + 63) / 64, sizeof(*r_surfmarks));
		r_flushmarks = calloc((need + 63) / 64, sizeof(*r_flushmarks));
		r_surfdsurf = malloc(need * sizeof(*r_surfdsurf));
		r_surfscan = calloc(need, sizeof(*r_surfscan));
		if(!r_surfmarks || !r_flushmarks || !r_surfdsurf
				|| !r_surfscan)
			Sys_Error("R_MarksBeginScan: out of memory");
		cap = need;
	}
	if(++r_scan == (u32)-1){ // that's D_BUSY, start over
		memset(r_surfscan, 0, cap * sizeof(*r_surfscan));
		r_scan = 1;
	}
	if(r_threaded){ // cleared again, in case a frame was cut short
		s32 words = (cap + 63) / 64;
		memset(r_surfmarks, 0, words * sizeof(*r_surfmarks));
		memset(r_flushmarks, 0, words * sizeof(*r_flushmarks));
		D_BeginRecords();
	}
}

// Decides whether this frame has any cutout/skybox surfaces at all and, if so,
// resets the accumulators of sc, if any. Must run after the surfaces have been
// emitted.
static void R_CutoutBeginFrame(scan_t *sc)
{
	cut_surfs = 0;
	for(surf_t *s = &surfaces[1]; s < surface_p; s++)
		if(R_IsCutSurf(s))
			cut_surfs++;
	cut_right0 = (s32)r_refdef.vrectright;
	if(!sc)
		return;
	sc->cut_count = 0;
	sc->cut_nedges = 0;
	if(!cut_surfs)
		return;
	s32 need = r_lsurfs_size + 2;
	if(sc->cut_cap < need){
		free(sc->cut_acc);
		free(sc->cut_list);
		sc->cut_acc = malloc(need * sizeof(*sc->cut_acc));
		sc->cut_list = malloc(need * sizeof(*sc->cut_list));
		if(!sc->cut_acc || !sc->cut_list)
			Sys_Error("R_CutoutBeginFrame: out of memory");
		sc->cut_cap = need;
	}
	if(sc->cut_edgecap < r_ledges_size){
		free(sc->cut_edges);
		sc->cut_edges = malloc(r_ledges_size * sizeof(*sc->cut_edges));
		if(!sc->cut_edges)
			Sys_Error("R_CutoutBeginFrame: out of memory");
		sc->cut_edgecap = r_ledges_size;
	}
	for(surf_t *s = &surfaces[1]; s < surface_p; s++){
		if(!R_IsCutSurf(s))
			continue;
		sc->cut_acc[s - surfaces].left = CUT_IDLE;
		sc->cut_acc[s - surfaces].right = cut_right0;
	}
}

// With the workers, the rows are cut into bands, the draw bands, and each one
// is scanned by itself on any thread (R_ScanBand). Then the main thread sets
// up and hands over their spans the way the scan from the top would have.
#define R_MAXSCANBANDS 256
#define R_EDGEBLOCK 4096 // a band's own edges come this many at a time
typedef struct { edge_t *e; s64 u0; } carry_t;
typedef struct { s64 u; edge_t *e; u32 s0, s1; } carrysort_t;
typedef struct {
	s32 y0, y1, job;
	carry_t *carry; // the edges from rows above, with u on their first row
	s64 ncarry, carrycap;
	s32 *rowspans; // the spans each row got
	s64 rowcap;
	scanrec_t *recs; // the surfaces with spans by draw band, by number
	s64 nrecs, reccap;
	s32 *recend; // where each draw band's records end
	s64 nwin, wincap;
	s32 *left; // surface, spanstate: those still in a span after its rows
	s64 nleft, leftcap;
	espan_t **chunks; // the span memory, in use until the frame is drawn
	s64 *chunksize;
	s64 nchunks, chunkcap, chunk;
	espan_t *span_p, *chunk_end;
	s64 frame; // the frame the chunks are in use for
} scanband_t;
// What a thread scans a band with. A band's surfaces get numbers of their own,
// in the order the band meets them, and its edges are copies with those.
typedef struct {
	scan_t sc;
	u64 *map; // by surface: the band's stamp << 32 | its number for it
	s64 mapcap;
	u32 stamp;
	surf_t *ls; // the band's copies of the surfaces, [1] is the background
	u32 *gid; // their numbers in surfaces
	s64 lscap, nls;
	edge_t **eblocks;
	s64 neblocks, eblockcap;
	s32 eb, ei; // the next free edge
	carrysort_t *sort;
	s64 sortcap;
	edge_t **row; // a row's new edges, oldest first
	s64 rowcap;
} scratch_t;
static scanband_t r_sbands[R_MAXSCANBANDS];
static s32 r_nsbands;
static scratch_t r_scratch[R_MAXWORKERS + 1]; // by r_slot
bool r_scanbanded;
s16 r_bandof[MAXHEIGHT];

static void R_ScanBandsBegin()
{
	r_nsbands = 0;
	r_scanbanded = r_threaded && !r_draworder.value;
	if(!r_scanbanded)
		return;
	s32 top = r_refdef.vrect.y, bottom = r_refdef.vrectbottom;
	s32 rows = r_bandrows;
	while((bottom - top) / rows + 2 > R_MAXSCANBANDS)
		rows += r_bandrows;
	for(s32 y = top; y < bottom;){
		s32 y1 = q_min((y / rows + 1) * rows, bottom);
		scanband_t *b = &r_sbands[r_nsbands];
		b->y0 = y;
		b->y1 = y1;
		b->ncarry = 0;
		for(s32 i = y; i < y1; i++)
			r_bandof[i] = r_nsbands;
		r_nsbands++;
		y = y1;
	}
}

// R_EmitEdge: e goes on into the bands after b0 up to b1
void R_CarryEdge(edge_t *e, s32 b0, s32 b1)
{
	for(s32 i = b0 + 1; i <= b1; i++){
		scanband_t *b = &r_sbands[i];
		b->carry = R_Grow(b->carry, &b->carrycap, b->ncarry + 1,
				sizeof(*b->carry));
		b->carry[b->ncarry].e = e;
		b->carry[b->ncarry++].u0 = e->u;
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
	skybox_surf_p = NULL; // until R_EmitSkyBox, if it's drawn
	if(r_draworder.value){ //put background behind everything in the world
		pdrawfunc = r_threaded ? R_GenerateSpansBackwardMarked
			: R_GenerateSpansBackward;
		surfaces[1].key = 0;
		r_currentkey = 1;
	}else{
		pdrawfunc = r_threaded ? R_GenerateSpansMarked
			: R_GenerateSpans;
		surfaces[1].key = 0x7FFFFFFF;
		r_currentkey = 0;
	}
	if(r_alphapass)
		pdrawfunc = r_threaded ? R_GenerateSpansAlphaMarked
			: R_GenerateSpansAlpha;
	memset(newedges + r_refdef.vrect.y, 0, (r_refdef.vrectbottom -
				r_refdef.vrect.y) * sizeof(*newedges));
	R_ScanBandsBegin();
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

static void R_GrowKeys(scan_t *sc, s64 n)
{
	if(sc->ne_cap >= n)
		return;
	s64 cap = sc->ne_cap ? sc->ne_cap : 1024;
	while(cap < n)
		cap *= 2;
	sc->ne_keys = realloc(sc->ne_keys, 2 * cap * sizeof(*sc->ne_keys));
	if(!sc->ne_keys)
		Sys_Error("R_SortNewEdges: out of memory");
	sc->ne_cap = cap;
}

Q_SIL u64 R_NewEdgeKey(const edge_t *e, u64 idx)
{
	u64 trail = e->surfs[0] != 0;
	return (u64)e->u << (NE_IDXBITS + 1) | trail << NE_IDXBITS |
		(trail ? idx : NE_IDXMASK - idx);
}

Q_SIL u64 R_NewEdgeIndex(u64 key)
{
	u64 low = key & NE_IDXMASK;
	return (key >> NE_IDXBITS) & 1 ? low : NE_IDXMASK - low;
}

static edge_t *R_SortNewEdges(scan_t *sc, edge_t *list)
{
	if(!list->next)
		return list;
	s64 n = 0;
	for(edge_t *e = list; e; e = e->next){
		if(n == sc->ne_cap)
			R_GrowKeys(sc, n + 1);
		sc->ne_keys[n++] = R_NewEdgeKey(e, e - r_edges);
	}
	R_SortKeys(sc->ne_keys, sc->ne_keys + sc->ne_cap, n);
	edge_t *head = NULL, **link = &head;
	for(s64 i = 0; i < n; i++){
		edge_t *e = r_edges + R_NewEdgeIndex(sc->ne_keys[i]);
		*link = e;
		link = &e->next;
	}
	*link = NULL;
	return head;
}

// Does the edge only have cutout/skybox surfaces?
Q_SIL bool R_CutOnly(const edge_t *e)
{
	u32 s0 = e->surfs[0], s1 = e->surfs[1];
	return (!s0 || R_IsCutSurf(&surfaces[s0]))
		&& (!s1 || R_IsCutSurf(&surfaces[s1]));
}

// Moves the edges that only cutout/skybox surfaces have from a scanline's new
// edges to cut_edges, returns the others
static edge_t *R_SplitCutEdges(scan_t *sc, edge_t *list)
{
	edge_t *keep = NULL, **link = &keep;
	cutedge_t *ce = sc->cut_edges + sc->cut_nedges;
	for(edge_t *e = list; e; e = e->next){
		if(!R_CutOnly(e)){
			*link = e;
			link = &e->next;
			continue;
		}
		ce->u = e->u;
		ce->u_step = e->u_step;
		ce->surfs[0] = e->surfs[0];
		ce->surfs[1] = e->surfs[1];
		ce->vlast = e->vlast;
		ce++;
	}
	*link = NULL;
	sc->cut_nedges = ce - sc->cut_edges;
	return keep;
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

// Adds the span from u to end (exclusive) on the current scanline. mark 1 sets
// the surface's bit for the workers, 2 lists it for the band's draw band.
static inline void R_EmitSpan(scan_t *sc, surf_t *surf, s64 u, s64 end,
		s32 mark)
{
	espan_t *span = sc->span_p++;
	span->u = u;
	span->count = end - u;
	span->v = sc->iv;
	if(mark && !surf->spans){
		u64 i = surf - sc->sbase;
		if(mark == 1)
			r_surfmarks[i >> 6] |= (u64)1 << (i & 63);
		else
			sc->touched[sc->ntouched++] = i;
	}
	span->pnext = surf->spans;
	surf->spans = span;
}

Q_SIL void R_CleanupSpan(scan_t *sc, s32 mark)
{
	// now that we've reached the right edge of the screen, we're done with
	// any unfinished surfaces, so emit a span for whatever's on top
	surf_t *surf = sc->sbase[1].next;
	if(edge_tail_u_shift20 > surf->last_u)
		R_EmitSpan(sc, surf, surf->last_u, edge_tail_u_shift20, mark);
	do{ // reset spanstate for all surfaces in the surface stack
		surf->spanstate = 0;
		surf = surf->next;
	}while(surf != &sc->sbase[1]);
}

static inline void R_LeadingEdgeBackwards(scan_t *sc, const edge_t *edge,
		surf_t *surf, s32 mark)
{
	// it's adding a new surface in, so find the correct place
	// don't start a span if this is an inverted span, with the end edge
	// preceding the start edge(that is, we've already seen the end edge)
	surf_t *surf2;
	s32 iu; // keep here for the OpenBSD compiler
	if(++surf->spanstate == 1){
		surf2 = sc->sbase[1].next;
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
		if(iu > surf2->last_u)
			R_EmitSpan(sc, surf2, surf2->last_u, iu, mark);
		surf->last_u = iu; // set last_u on the new span
gotposition:
		surf->next = surf2; // insert before surf2
		surf->prev = surf2->prev;
		surf2->prev->next = surf;
		surf2->prev = surf;
	}
}

static inline void R_TrailingEdge(scan_t *sc, surf_t *surf, const edge_t *edge,
		s32 mark)
{
	// don't generate a span if this is an inverted span, with the end edge
	// preceding the start edge(that is, we haven't seen the start edge yet)
	if(--surf->spanstate == 0){
		if(surf == sc->sbase[1].next){
			s64 iu = edge->u >> 20; // emit a span
			if(iu > surf->last_u) // (current top going away)
				R_EmitSpan(sc, surf, surf->last_u, iu, mark);
			// set last_u on the surface below
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
static inline bool R_BmodelSurfInFront(const scan_t *sc, const edge_t *edge,
		const surf_t *surf, const surf_t *surf2)
{
	f32 fv = sc->fv;
	f64 fu = (f32)(edge->u - 0xFFFFF) * (1.0 / 0x100000) + 1;
	f64 newzi = surf->d_ziorigin + fv*surf->d_zistepv + fu*surf->d_zistepu;
	f64 testzi=surf2->d_ziorigin+fv*surf2->d_zistepv+fu*surf2->d_zistepu;
	return newzi >= testzi;
}

// surf is &sbase[edge->surfs[1]], which the caller has already checked is
// not the "no surface" index 0
static inline void R_LeadingEdge(scan_t *sc, const edge_t *edge, surf_t *surf,
		s32 mark)
{
	surf_t *surf2;
	s32 iu; // keep here for the OpenBSD compiler
	// it's adding a new surface in, so find the correct place
	// don't start a span if this is an inverted span, with the end
	// edge preceding the start edge(that is, we've already seen
	// the end edge)
	if(++surf->spanstate == 1){
		surf2 = sc->sbase[1].next;
		if(surf->key < surf2->key)
			goto newtop;
	//if it's two surfaces on the same plane, the one that's already
	//active is in front, so keep going unless it's a bmodel
		// must be two bmodels in the same leaf; sort on 1/z
		if(surf->insubmodel && (surf->key == surf2->key) &&
				R_BmodelSurfInFront(sc, edge, surf, surf2))
			goto newtop;
		if(!surf->insubmodel){
		// The active surface list is sorted by key, with the background
		// at both ends. Hidden faces tend to have big keys, so when
		// ours is closer to the far end, search from there instead of
		// walking the whole stack. Same insertion point as the forward
		// search below: after every surface whose key is <= ours. (The
		// top surface is a real one here, since our key is not lower
		// than its key, so the far end is a real surface too.)
			surf_t *last = sc->sbase[1].prev;
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
			if(R_BmodelSurfInFront(sc, edge, surf, surf2))
				goto gotposition;
			goto continue_search;
		}
		goto gotposition;
newtop:
		// emit a span(obscures current top)
		iu = edge->u >> 20;
		if(iu > surf2->last_u)
			R_EmitSpan(sc, surf2, surf2->last_u, iu, mark);
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

static inline void R_CutTrailing(scan_t *sc, u32 sn, s64 u)
{
	cutspan_t *acc = &sc->cut_acc[sn];
	s64 iu = u >> 20;
	if(acc->left == CUT_IDLE){
		acc->left = CUT_NOLEAD;
		sc->cut_list[sc->cut_count++] = sn;
	}
	if(iu < acc->right)
		acc->right = iu;
}

static inline void R_CutLeading(scan_t *sc, u32 sn, s64 u)
{
	cutspan_t *acc = &sc->cut_acc[sn];
	s64 iu = u >> 20;
	if(iu < 0)
		iu = 0;
	if(acc->left < 0){ // first leading edge
		if(acc->left == CUT_IDLE)
			sc->cut_list[sc->cut_count++] = sn;
		acc->left = iu;
	}else if(iu > acc->left){
		acc->left = iu;
	}
}

// The cutout edges on this scanline, and then on to the next one
Q_SIL void R_CutEdges(scan_t *sc, bool step)
{
	cutedge_t *e = sc->cut_edges;
	s32 n = sc->cut_nedges, iv = sc->iv;
	for(s32 i = 0; i < n;){
		if(e[i].surfs[0])
			R_CutTrailing(sc, e[i].surfs[0], e[i].u);
		if(e[i].surfs[1])
			R_CutLeading(sc, e[i].surfs[1], e[i].u);
		if(step && e[i].vlast == iv){
			e[i] = e[--n]; // done with it, the last one moves in
			continue;
		}
		if(step)
			e[i].u += e[i].u_step;
		i++;
	}
	sc->cut_nedges = n;
}

// every edge has been seen by now. An inverted span(trailing before
// leading) gives right <= left and is skipped
Q_SIL void R_EmitCutSpansImpl(scan_t *sc, s32 mark)
{
	for(s32 i = 0; i < sc->cut_count; i++){
		s32 sn = sc->cut_list[i];
		cutspan_t *acc = &sc->cut_acc[sn];
		// left < 0: only trailing edges seen, nothing to draw
		if(acc->left >= 0 && acc->right > acc->left)
			R_EmitSpan(sc, &sc->sbase[sn], acc->left, acc->right,
					mark);
		acc->left = CUT_IDLE;
		acc->right = cut_right0;
	}
	sc->cut_count = 0;
}

static void R_EmitCutSpans(scan_t *sc){ R_EmitCutSpansImpl(sc, 0); }
static void R_EmitCutSpansMarked(scan_t *sc){ R_EmitCutSpansImpl(sc, 1); }
static void R_EmitCutSpansBand(scan_t *sc){ R_EmitCutSpansImpl(sc, 2); }

Q_SIL void R_EmitCutSpansAs(scan_t *sc, s32 mark)
{
	if(mark == 2) R_EmitCutSpansBand(sc);
	else if(mark) R_EmitCutSpansMarked(sc);
	else R_EmitCutSpans(sc);
}

Q_SIL void R_GenerateSpansImpl(scan_t *sc, bool step, s32 mark)
{
	surf_t *const sbase = sc->sbase; // the global is reloaded after every
				        // store through a pointer otherwise
	const edge_t *const tail = &sc->tail;
	const s32 iv = sc->iv;
	// clear active surfaces to just the background surface
	sbase[1].next = sbase[1].prev = &sbase[1];
	sbase[1].last_u = edge_head_u_shift20;
	if(!cut_surfs){ // no cutout or skybox surfaces in this frame
		s64 prev_u = sc->head.u;
		for(edge_t *edge=sc->head.next,*next; edge!=tail; edge=next){
			next = edge->next;
			// A shared(cached) edge can carry a trailing surface in
			// surfs[0] AND a leading one in surfs[1]. Always retire
			// the trailing surface first.
			u32 s0 = edge->surfs[0], s1 = edge->surfs[1];
			if(s0)
				R_TrailingEdge(sc, &sbase[s0], edge, mark);
			if(s1)
				R_LeadingEdge(sc, edge, &sbase[s1], mark);
			if(step)
				prev_u = R_StepEdge(edge, next, prev_u, iv);
		}
		R_CleanupSpan(sc, mark);
		return;
	}
	// Cutout and skybox surfaces get a span covering the whole of the
	// surface on this scanline, from its rightmost leading edge to its
	// leftmost trailing edge, with intentional overdraw. They don't take
	// part in the active surface stack, so they only need those two
	// numbers, and we get them as the edges go by.
	s64 prev_u = sc->head.u;
	for(edge_t *edge=sc->head.next,*next; edge!=tail; edge=next){
		next = edge->next;
		u32 s0 = edge->surfs[0], s1 = edge->surfs[1];
		if(s0){
			surf_t *surf = &sbase[s0];
			if(surf->flags & SURF_CUTMASK)
				R_CutTrailing(sc, s0, edge->u);
			else
				R_TrailingEdge(sc, surf, edge, mark);
		}
		if(s1){
			surf_t *surf = &sbase[s1];
			if(surf->flags & SURF_CUTMASK)
				R_CutLeading(sc, s1, edge->u);
			else
				R_LeadingEdge(sc, edge, surf, mark);
		}
		if(step)
			prev_u = R_StepEdge(edge, next, prev_u, iv);
	}
	R_CutEdges(sc, step);
	R_EmitCutSpansAs(sc, mark);
	R_CleanupSpan(sc, mark);
}

// The alpha pass: every surface is translucent, so none of them may hide
// another. They all get whole spans like cutouts, the z-buffer of the opaque
// pass does the occlusion and D_DrawSurfacesAlpha() draws them back to front.
Q_SIL void R_GenerateSpansAlphaImpl(scan_t *sc, bool step, s32 mark)
{
	const edge_t *const tail = &sc->tail;
	const s32 iv = sc->iv;
	s64 prev_u = sc->head.u;
	for(edge_t *edge=sc->head.next,*next; edge!=tail; edge=next){
		next = edge->next;
		u32 s0 = edge->surfs[0], s1 = edge->surfs[1];
		if(s0)
			R_CutTrailing(sc, s0, edge->u);
		if(s1)
			R_CutLeading(sc, s1, edge->u);
		if(step)
			prev_u = R_StepEdge(edge, next, prev_u, iv);
	}
	R_CutEdges(sc, step);
	R_EmitCutSpansAs(sc, mark);
}

Q_SIL void R_GenerateSpansBackwardImpl(scan_t *sc, bool step, s32 mark)
{
	surf_t *const sbase = sc->sbase;
	const edge_t *const tail = &sc->tail;
	const s32 iv = sc->iv;
	// clear active surfaces to just the background surface
	sbase[1].next = sbase[1].prev = &sbase[1];
	sbase[1].last_u = edge_head_u_shift20;
	// generate spans
	s64 prev_u = sc->head.u;
	for(edge_t *edge=sc->head.next,*next; edge!=tail; edge=next){
		next = edge->next;
		u32 s0 = edge->surfs[0], s1 = edge->surfs[1];
		if(s0)
			R_TrailingEdge(sc, &sbase[s0], edge, mark);
		if(s1)
			R_LeadingEdgeBackwards(sc, edge, &sbase[s1], mark);
		if(step)
			prev_u = R_StepEdge(edge, next, prev_u, iv);
	}
	R_CleanupSpan(sc, mark);
}

// The span generation, the same marking surfaces for the workers, and for a
// band (R_ScanBand). The main thread's state is known where it's used, so its
// fields are plain globals to the compiler.
static void R_GenerateSpans(bool step)
{ R_GenerateSpansImpl(&r_scanmain, step, 0); }
static void R_GenerateSpansMarked(bool step)
{ R_GenerateSpansImpl(&r_scanmain, step, 1); }
static void R_GenerateSpansBand(scan_t *sc, bool step)
{ R_GenerateSpansImpl(sc, step, 2); }
static void R_GenerateSpansAlpha(bool step)
{ R_GenerateSpansAlphaImpl(&r_scanmain, step, 0); }
static void R_GenerateSpansAlphaMarked(bool step)
{ R_GenerateSpansAlphaImpl(&r_scanmain, step, 1); }
static void R_GenerateSpansAlphaBand(scan_t *sc, bool step)
{ R_GenerateSpansAlphaImpl(sc, step, 2); }
static void R_GenerateSpansBackward(bool step)
{ R_GenerateSpansBackwardImpl(&r_scanmain, step, 0); }
static void R_GenerateSpansBackwardMarked(bool step)
{ R_GenerateSpansBackwardImpl(&r_scanmain, step, 1); }

static void R_DrawSurfaceSpans(s32 act)
{
	if(r_drawflat.value) D_DrawSurfacesFlat(act);
	else if(!r_alphapass) D_DrawSurfaces(act);
	else D_DrawSurfacesAlpha(act);
}

// The workers get the spans of a band once its last scanline is in
static void R_RecordBand(s32 band)
{
	r_recordband = band;
	R_DrawSurfaceSpans(SURFS_RECORD);
	R_BandPublish(band);
}

// The span buffer. Without the workers it's MAXSPANS between flushes, plus
// room for one full scanline of spans: at most one per pixel for the visible
// surfaces, plus at most one per cutout/skybox surface. With them the spans
// stay until the workers are done, the alpha pass goes on after the main pass,
// and the buffer is that many times bigger. They flush at the same points,
// that's when the surfaces get set up and the cache filled (D_DrawSurfaces).
// Scanned in bands, the spans are the bands' own, see R_BandChunk.
static espan_t *basespans;
static s64 span_cap;
static s64 flush_cap; // the buffer without the workers, as it would have grown

// Sets span_p for a scan, returns how many spans make it flush
static s64 R_SpanBuffer(s64 room)
{
	static s64 frame = -1; // of the spans the workers may still use
	if(flush_cap < MAXSPANS + room)
		flush_cap = (MAXSPANS + room) * 3 / 2;
	if(r_scanbanded)
		return flush_cap - room;
	s64 need = r_threaded ? flush_cap * 4 : flush_cap;
	bool cont = r_threaded && frame == r_framecount;
	frame = r_threaded ? r_framecount : -1;
	if(span_cap < need){
		R_ThreadsSync();
		span_cap = need;
		free(basespans);
		basespans = malloc(span_cap * sizeof(*basespans));
		if(!basespans)
			Sys_Error("R_ScanEdges: out of memory");
		cont = 0;
	}
	if(!cont || r_scanmain.span_p + flush_cap + room
			> basespans + span_cap){
		R_ThreadsSync();
		r_scanmain.span_p = basespans;
	}
	return flush_cap - room;
}

// Nothing active yet but the background edges around the whole view
static void R_ScanBeginEdges(scan_t *sc)
{
	sc->head.u = (s64)r_refdef.vrect.x << 20;
	sc->head.u_step = 0;
	sc->head.prev = NULL;
	sc->head.next = &sc->tail;
	sc->head.surfs[0] = 0;
	sc->head.surfs[1] = 1;
	sc->tail.u = ((s64)r_refdef.vrectright << 20) + 0xFFFFF;
	sc->tail.u_step = 0;
	sc->tail.prev = &sc->head;
	sc->tail.next = &sc->aftertail;
	sc->tail.surfs[0] = 1;
	sc->tail.surfs[1] = 0;
	sc->aftertail.u = -1; // force a move
	sc->aftertail.u_step = 0;
	sc->aftertail.next = &sc->sentinel;
	sc->aftertail.prev = &sc->tail;
	sc->sentinel.u = ((u64)2000) << 24; // nothing sorts past this
	sc->sentinel.prev = &sc->aftertail;
}

static void R_ScratchGrow(scratch_t *x)
{ // between scanlines, while no surface is on the stack
	scan_t *sc = &x->sc;
	s64 cap = x->lscap;
	x->ls = R_Grow(x->ls, &cap, x->nls + 1, sizeof(*x->ls));
	cap = x->lscap;
	x->gid = R_Grow(x->gid, &cap, x->nls + 1, sizeof(*x->gid));
	cap = x->lscap;
	sc->cut_acc = R_Grow(sc->cut_acc, &cap, x->nls + 1,
			sizeof(*sc->cut_acc));
	cap = x->lscap;
	sc->cut_list = R_Grow(sc->cut_list, &cap, x->nls + 1,
			sizeof(*sc->cut_list));
	cap = x->lscap;
	sc->touched = R_Grow(sc->touched, &cap, x->nls + 1,
			sizeof(*sc->touched));
	x->lscap = cap;
	sc->sbase = x->ls;
}

// The band's number for surface g, which it copies the first time
static u32 R_BandSurf(scratch_t *x, u32 g)
{
	if(!g)
		return 0;
	u64 m = x->map[g];
	if((u32)(m >> 32) == x->stamp)
		return (u32)m;
	if(x->nls >= x->lscap)
		R_ScratchGrow(x);
	u32 l = x->nls++;
	x->map[g] = (u64)x->stamp << 32 | l;
	x->gid[l] = g;
	const surf_t *f = &surfaces[g];
	surf_t *s = &x->ls[l];
	s->spans = NULL;
	s->spanstate = 0;
	s->key = f->key;
	s->flags = f->flags;
	s->insubmodel = f->insubmodel;
	s->d_ziorigin = f->d_ziorigin;
	s->d_zistepu = f->d_zistepu;
	s->d_zistepv = f->d_zistepv;
	if(r_alphapass ? g != 1 : (f->flags & SURF_CUTMASK) != 0){
		x->sc.cut_acc[l].left = CUT_IDLE;
		x->sc.cut_acc[l].right = cut_right0;
	}
	return l;
}

static void R_ScratchBegin(scratch_t *x)
{
	s64 need = r_lsurfs_size + 2;
	if(x->mapcap < need){
		free(x->map);
		if(!(x->map = calloc(need, sizeof(*x->map))))
			Sys_Error("R_ScanBand: out of memory");
		x->mapcap = need;
		x->stamp = 0;
	}
	if(!++x->stamp){ // a new number for every band scanned
		memset(x->map, 0, x->mapcap * sizeof(*x->map));
		x->stamp = 1;
	}
	x->nls = 1;
	R_BandSurf(x, 1); // the background
	x->eb = -1;
	x->ei = R_EDGEBLOCK;
	x->sc.sbase = x->ls;
	x->sc.cut_count = 0;
	x->sc.cut_nedges = 0;
	x->sc.ntouched = 0;
}

// An edge of the band's own, they stay where they are until the next band
static edge_t *R_BandEdge(scratch_t *x)
{
	if(x->ei == R_EDGEBLOCK){
		if(++x->eb == x->neblocks){
			x->eblocks = R_Grow(x->eblocks, &x->eblockcap,
				x->neblocks + 1, sizeof(*x->eblocks));
			x->eblocks[x->neblocks] =
				malloc(R_EDGEBLOCK * sizeof(edge_t));
			if(!x->eblocks[x->neblocks++])
				Sys_Error("R_ScanBand: out of memory");
		}
		x->ei = 0;
	}
	return &x->eblocks[x->eb][x->ei++];
}

static edge_t *R_BandEdgeOf(scratch_t *x, s64 u, const edge_t *e, u32 s0,
		u32 s1)
{
	edge_t *p = R_BandEdge(x);
	p->u = u;
	p->u_step = e->u_step;
	p->surfs[0] = s0;
	p->surfs[1] = s1;
	p->vlast = e->vlast;
	return p;
}

static void R_BandCutEdge(scratch_t *x, s64 u, const edge_t *e, u32 s0,
		u32 s1)
{
	scan_t *sc = &x->sc;
	sc->cut_edges = R_Grow(sc->cut_edges, &sc->cut_edgecap,
			sc->cut_nedges + 1, sizeof(*sc->cut_edges));
	cutedge_t *ce = &sc->cut_edges[sc->cut_nedges++];
	ce->u = u;
	ce->u_step = e->u_step;
	ce->surfs[0] = s0;
	ce->surfs[1] = s1;
	ce->vlast = e->vlast;
}

// The order of the edges from above on the band's first row. Stepping keeps
// the order of edges that end up on the same u, so the one that came from
// further left is in front, and of those on the same line the newer one, as
// R_InsertNewEdges put it in front. New on the same row, it's R_SortNewEdges'.
static s32 R_CarryCmp(const void *pa, const void *pb)
{
	const carrysort_t *a = pa, *b = pb;
	if(a->u != b->u)
		return a->u < b->u ? -1 : 1;
	const edge_t *x = a->e, *y = b->e;
	if(x->u_step != y->u_step)
		return x->u_step > y->u_step ? -1 : 1;
	if(x->vstart != y->vstart)
		return x->vstart > y->vstart ? -1 : 1;
	bool tx = x->surfs[0] != 0, ty = y->surfs[0] != 0;
	if(tx != ty)
		return tx ? 1 : -1;
	if(!tx)
		return x > y ? -1 : 1;
	return x < y ? -1 : 1;
}

// The band's copies of the edges from above, active as they would be on its
// first row
static void R_BandCarry(scratch_t *x, const scanband_t *b, bool split)
{
	scan_t *sc = &x->sc;
	x->sort = R_Grow(x->sort, &x->sortcap, b->ncarry, sizeof(*x->sort));
	s64 n = 0;
	for(s64 i = 0; i < b->ncarry; i++){
		edge_t *e = b->carry[i].e;
		s64 u = b->carry[i].u0 + (s64)(b->y0 - e->vstart) * e->u_step;
		u32 s0 = R_BandSurf(x, e->surfs[0]);
		u32 s1 = R_BandSurf(x, e->surfs[1]);
		if(split && R_CutOnly(e)){
			R_BandCutEdge(x, u, e, s0, s1);
			continue;
		}
		carrysort_t *c = &x->sort[n++];
		c->u = u;
		c->e = e;
		c->s0 = s0;
		c->s1 = s1;
	}
	qsort(x->sort, n, sizeof(*x->sort), R_CarryCmp);
	edge_t *prev = &sc->head;
	for(s64 i = 0; i < n; i++){
		const carrysort_t *c = &x->sort[i];
		edge_t *p = R_BandEdgeOf(x, c->u, c->e, c->s0, c->s1);
		p->prev = prev;
		prev->next = p;
		prev = p;
	}
	prev->next = &sc->tail;
	sc->tail.prev = prev;
}

// A row's new edges for a band: copies in the order R_SortNewEdges would put
// them, but for the cutout ones, which go to cut_edges. R_EmitEdge puts the
// newest first, so an edge's place counted from the end goes by its index.
static edge_t *R_BandRow(scratch_t *x, edge_t *list, bool split)
{
	scan_t *sc = &x->sc;
	s64 n = 0, k = 0, i = 0;
	for(edge_t *e = list; e; e = e->next)
		n++;
	x->row = R_Grow(x->row, &x->rowcap, n, sizeof(*x->row));
	R_GrowKeys(sc, n);
	for(edge_t *e = list; e; e = e->next, i++){
		u32 s0 = R_BandSurf(x, e->surfs[0]);
		u32 s1 = R_BandSurf(x, e->surfs[1]);
		if(split && R_CutOnly(e)){
			R_BandCutEdge(x, e->u, e, s0, s1);
			continue;
		}
		edge_t *p = R_BandEdgeOf(x, e->u, e, s0, s1);
		x->row[n - 1 - i] = p;
		sc->ne_keys[k++] = R_NewEdgeKey(p, n - 1 - i);
	}
	R_SortKeys(sc->ne_keys, sc->ne_keys + sc->ne_cap, k);
	edge_t *head = NULL, **link = &head;
	for(s64 j = 0; j < k; j++){
		edge_t *p = x->row[R_NewEdgeIndex(sc->ne_keys[j])];
		*link = p;
		link = &p->next;
	}
	*link = NULL;
	return head;
}

static s32 R_RecCmp(const void *pa, const void *pb)
{
	u32 a = ((const scanrec_t *)pa)->surf;
	u32 b = ((const scanrec_t *)pb)->surf;
	return a < b ? -1 : a > b;
}

// The draw band ends: its surfaces and their spans go to the records
static void R_BandWindow(scanband_t *b, scratch_t *x)
{
	scan_t *sc = &x->sc;
	b->recs = R_Grow(b->recs, &b->reccap, b->nrecs + sc->ntouched,
			sizeof(*b->recs));
	scanrec_t *r = b->recs + b->nrecs;
	for(s32 i = 0; i < sc->ntouched; i++){
		surf_t *s = &x->ls[sc->touched[i]];
		r[i].surf = x->gid[sc->touched[i]];
		r[i].spans = s->spans;
		s->spans = NULL;
	}
	qsort(r, sc->ntouched, sizeof(*r), R_RecCmp);
	b->nrecs += sc->ntouched;
	sc->ntouched = 0;
	b->recend = R_Grow(b->recend, &b->wincap, b->nwin + 1,
			sizeof(*b->recend));
	b->recend[b->nwin++] = b->nrecs;
}

// Room for at least room more spans
static void R_BandChunk(scanband_t *b, s64 room)
{
	if(b->frame != r_framecount){ // the workers are done with all of them
		b->frame = r_framecount;
		b->chunk = 0;
	}else
		b->chunk++;
	if(b->chunk == b->nchunks){
		s64 cap = b->chunkcap;
		b->chunks = R_Grow(b->chunks, &cap, b->nchunks + 1,
				sizeof(*b->chunks));
		b->chunksize = R_Grow(b->chunksize, &b->chunkcap,
				b->nchunks + 1, sizeof(*b->chunksize));
		b->chunks[b->nchunks] = NULL;
		b->chunksize[b->nchunks++] = 0;
	}
	s64 i = b->chunk;
	if(b->chunksize[i] < room){
		s64 size = q_max(room * 4, 16384);
		free(b->chunks[i]);
		if(!(b->chunks[i] = malloc(size * sizeof(espan_t))))
			Sys_Error("R_ScanBand: out of memory");
		b->chunksize[i] = size;
	}
	b->span_p = b->chunks[i];
	b->chunk_end = b->chunks[i] + b->chunksize[i];
}

// Scans the rows of a band, starting with the surfaces in left that a band
// above left in the middle of a span (surface, spanstate pairs). The spans of
// each draw band go to the records, and the span counts of each row.
static void R_ScanBand(scanband_t *b, const s32 *left, s64 nleft)
{
	scratch_t *x = &r_scratch[r_slot];
	scan_t *sc = &x->sc;
	bool split = cut_surfs && (r_alphapass || !r_draworder.value);
	R_ScratchBegin(x);
	b->nrecs = b->nwin = b->nleft = 0;
	b->rowspans = R_Grow(b->rowspans, &b->rowcap, b->y1 - b->y0,
			sizeof(*b->rowspans));
	if(b->frame != r_framecount)
		b->chunk_end = NULL;
	R_ScanBeginEdges(sc);
	for(s64 i = 0; i < nleft; i += 2)
		x->ls[R_BandSurf(x, left[i])].spanstate = left[i + 1];
	R_BandCarry(x, b, split);
	for(s32 iv = b->y0; iv < b->y1; iv++){
		sc->iv = iv;
		sc->fv = (f32)iv;
		edge_t *add = newedges[iv];
		if(add)
			add = R_BandRow(x, add, split);
		// one span per pixel at most, and per cutout surface
		s64 room = r_refdef.vrect.width + x->nls;
		if(!b->chunk_end || b->chunk_end - b->span_p < room)
			R_BandChunk(b, room);
		sc->span_p = b->span_p;
		sc->sbase[1].spanstate = 1;
		if(add)
			R_InsertNewEdges(add, sc->head.next);
		bool step = iv + 1 < b->y1;
		if(r_alphapass)
			R_GenerateSpansAlphaBand(sc, step);
		else
			R_GenerateSpansBand(sc, step);
		b->rowspans[iv - b->y0] = sc->span_p - b->span_p;
		b->span_p = sc->span_p;
		if(iv + 1 == b->y1 || (iv + 1) % r_bandrows == 0)
			R_BandWindow(b, x);
	}
	for(s64 l = 2; l < x->nls; l++){ // unbalanced edges, if ever
		if(!x->ls[l].spanstate)
			continue;
		b->left = R_Grow(b->left, &b->leftcap, b->nleft + 2,
				sizeof(*b->left));
		b->left[b->nleft++] = x->gid[l];
		b->left[b->nleft++] = x->ls[l].spanstate;
	}
}

static void R_ScanBandJob(void *arg){ R_ScanBand(arg, NULL, 0); }

// Hands the spans of records r0 to r1 to their surfaces for D_DrawSurfaces:
// those on rows up to f, or all of them if f < 0
static void R_GiveSpans(scanband_t *b, s32 r0, s32 r1, s32 f)
{
	for(s32 i = r0; i < r1; i++){
		scanrec_t *r = &b->recs[i];
		espan_t *s = r->spans, *prev = NULL;
		if(f >= 0) // the newest are first
			for(; s && s->v > f; s = s->pnext)
				prev = s;
		if(!s)
			continue;
		if(prev)
			prev->pnext = NULL;
		else
			r->spans = NULL;
		surfaces[r->surf].spans = s;
		r_surfmarks[r->surf >> 6] |= (u64)1 << (r->surf & 63);
	}
}

// The main pass's record walks between two flushes all go by the same surface
// cache, so they're jobs (D_RecordBandJob), and the main thread only waits for
// them before it sets up at the flush
static s32 *r_recjobs;
static s64 r_nrecjobs, r_recjobcap;

// The spans of records r0 to r1 as R_GiveSpans picks them, to a job of their
// band's record walk
static void R_RecordJob(scanband_t *b, s32 r0, s32 r1, s32 f, s32 band)
{
	drecjob_t *j = R_FrameAlloc(sizeof(*j) + (r1-r0) * sizeof(scanrec_t));
	scanrec_t *out = (scanrec_t *)(j + 1);
	s32 n = 0;
	for(s32 i = r0; i < r1; i++){
		scanrec_t *r = &b->recs[i];
		espan_t *s = r->spans, *prev = NULL;
		if(f >= 0)
			for(; s && s->v > f; s = s->pnext)
				prev = s;
		if(!s)
			continue;
		if(prev)
			prev->pnext = NULL;
		else
			r->spans = NULL;
		out[n].surf = r->surf;
		out[n++].spans = s;
		r_flushmarks[r->surf >> 6] |= (u64)1 << (r->surf & 63);
	}
	if(!n)
		return;
	j->band = band;
	j->n = n;
	j->recs = out;
	r_recjobs = R_Grow(r_recjobs, &r_recjobcap, r_nrecjobs + 1,
			sizeof(*r_recjobs));
	r_recjobs[r_nrecjobs++] = R_JobAdd(D_RecordBandJob, j);
	R_JobsPublish();
}

static void R_WaitRecordJobs()
{
	for(s64 i = 0; i < r_nrecjobs; i++)
		R_JobWait(r_recjobs[i]);
	r_nrecjobs = 0;
	D_StampBlocks(); // before the setup walk evicts any
}

// The bands are scanned by any thread, the main thread takes the first. Then
// it records them and sets up their surfaces in the order of the scan from the
// top, flushing where it would have.
static void R_ScanBands(s64 flush)
{
	bool jobs = !r_alphapass && !r_drawflat.value; // the record walks
	for(s32 i = 1; i < r_nsbands; i++)
		r_sbands[i].job = R_JobAdd(R_ScanBandJob, &r_sbands[i]);
	R_JobsPublish();
	R_ScanBand(&r_sbands[0], NULL, 0);
	s64 since = 0, bottom = r_refdef.vrectbottom - 1;
	for(s32 i = 0; i < r_nsbands; i++){
		scanband_t *b = &r_sbands[i];
		if(i){
			R_JobWait(b->job);
			const scanband_t *a = b - 1;
			if(a->nleft) // it started out in the middle of a span
				R_ScanBand(b, a->left, a->nleft);
		}
		s32 r0 = 0, w = 0;
		for(s32 iv = b->y0; iv < b->y1; iv++){
			since += b->rowspans[iv - b->y0];
			bool last = iv == bottom;
			bool end = iv + 1 == b->y1
				|| (iv + 1) % r_bandrows == 0;
			bool fl = !last && since >= flush;
			s32 band = iv / r_bandrows;
			if((end || fl) && jobs)
				R_RecordJob(b, r0, b->recend[w], end ? -1 : iv,
						band);
			else if(end){
				R_GiveSpans(b, r0, b->recend[w], -1);
				R_RecordBand(band);
			}else if(fl)
				R_GiveSpans(b, r0, b->recend[w], iv);
			if(end)
				r0 = b->recend[w++];
			if(!fl)
				continue;
			if(jobs)
				R_WaitRecordJobs();
			else
				R_RecordBand(band);
			R_DrawSurfaceSpans(SURFS_SETUP);
			since = 0;
		}
	}
	R_WaitRecordJobs();
	R_DrawSurfaceSpans(SURFS_SETUP);
}

void R_ScanEdges()
{
	scan_t *sc = &r_scanmain;
	edge_head_u_shift20 = r_refdef.vrect.x;
	edge_tail_u_shift20 = r_refdef.vrectright;
	R_MarksBeginScan();
	R_CutoutBeginFrame(r_scanbanded ? NULL : sc);
	// r_draworder has no cutouts of its own
	bool split = cut_surfs && (r_alphapass || !r_draworder.value);
	s64 room = r_refdef.vrect.width + cut_surfs;
	s64 flush = R_SpanBuffer(room);
	if(r_scanbanded){
		R_ScanBands(flush);
		return;
	}
	sc->sbase = surfaces;
	espan_t *flush_p = sc->span_p; // the spans since the last flush
	// clear active edges to just background edges around the whole screen
	R_ScanBeginEdges(sc);
	s64 bottom = r_refdef.vrectbottom - 1; // process all scan lines
	for(s64 iv = r_refdef.vrect.y; iv <= bottom; iv++){
		sc->iv = iv;
		sc->fv = (f32)iv;
		// mark that the head(background start) span is pre-included
		surfaces[1].spanstate = 1;
		edge_t *add = newedges[iv];
		if(add && split)
			add = R_SplitCutEdges(sc, add);
		if(add)
			R_InsertNewEdges(R_SortNewEdges(sc, add),
					sc->head.next);
		// also steps/removes the edges, no need on the last scan
		bool last = iv == bottom;
		(*pdrawfunc)(!last);
		if(r_threaded && (last || (iv + 1) % r_bandrows == 0))
			R_RecordBand(iv / r_bandrows);
		// flush the span list if we can't be sure we have enough spans
		// left for the next scan
		if(last || sc->span_p - flush_p < flush)
			continue;
		if(!r_threaded){
			R_DrawSurfaceSpans(SURFS_DRAW);
			flush_p = sc->span_p = basespans;
			continue;
		}
		R_RecordBand(iv / r_bandrows);
		R_DrawSurfaceSpans(SURFS_SETUP);
		if(sc->span_p + flush + 2 * room > basespans + span_cap){
			R_ThreadsSync(); // the workers must be done with them
			sc->span_p = basespans;
		}
		flush_p = sc->span_p;
	}
	R_DrawSurfaceSpans(r_threaded ? SURFS_SETUP : SURFS_DRAW);
}
