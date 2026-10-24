// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

typedef struct { f32 x; f32 y; } coarse_occ_vertex_t;
static f32 winquake_surface_liquid_alpha;
// One bit per tile, one row of 64-bit words per tile row.
static u64 *r_coarse_occlusion_bits;
static s32 r_coarse_occlusion_width;
static s32 r_coarse_occlusion_height;
static s32 r_coarse_occlusion_words; // u64 words per tile row
static f32 *r_coarse_occ_scratch; // 4 arrays of (height + 2) floats
static f32 r_coarse_occ_vw, r_coarse_occ_vh; // viewport size in pixels
static bool r_coarse_active; // set once per frame, see BeginFrame
// Polygon built by the occlusion test, reused when the face is added as an
// occluder right afterwards.
static msurface_t *r_coarse_poly_surf;
static coarse_occ_vertex_t r_coarse_poly[COARSE_OCCLUSION_MAX_VERTS];
static s32 r_coarse_poly_count;
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
static u8 r_edgerev; // direction of the face currently being walked

// r_coarseocclusion enables conservative coarse screen-space occlusion:
// - Static tile coverage bitmap, one bit per tile.
// - Only FULL tiles are recorded.
// - A tile becomes FULL only if all four of its corners are inside the
//   projected convex surface polygon.
// - A surface is rejected only if every tile intersected by its projected
//   polygon is already FULL.
// - Faces with any vertex behind NEAR_CLIP are not considered for coarse
//   occlusion at all. This is conservative around the near plane.
// - Screen-space clipping is performed against r_refdef.vrect.
//
// Implementation notes. Both the test and the fill treat the face as a convex
// polygon and work on whole tile rows ("strips") at once instead of running
// geometry tests per tile:
// - Test: every polygon edge is clipped against each strip it overlaps, which
//   gives the exact horizontal extent of the polygon inside that strip. All
//   tiles in that extent are intersected by the polygon (the polygon is
//   convex, so it is connected inside the strip), and are checked with a
//   couple of word operations.
// - Fill: a tile has all four corners inside a convex polygon exactly when its
//   left/right sides are inside the polygon's horizontal cross-sections at
//   both the top and the bottom line of the tile row. Cross-sections are
//   computed once per tile row boundary, and a whole run of tiles is set at
//   once.
// COARSE_OCCLUSION_MARGIN makes both operations err on the safe side: the
// test looks at slightly more tiles, the fill marks slightly fewer.

void R_CoarseOcclusionBeginFrame()
{
	r_coarse_rejected = 0;
	coccl_enable = (r_coarseocclusion.value == 1) ||
		(r_coarseocclusion.value && cl_entities[0].model &&
		 cl_entities[0].model->bspversion != BSPVERSION);
	r_coarse_active = coccl_enable && !r_showtris.value;
	if(!r_coarse_active)
		return;
	s32 width = (r_refdef.vrect.width + COARSE_OCCLUSION_TILE_SIZE - 1) /
		COARSE_OCCLUSION_TILE_SIZE;
	s32 height = (r_refdef.vrect.height + COARSE_OCCLUSION_TILE_SIZE - 1) /
		COARSE_OCCLUSION_TILE_SIZE;
	s32 words = (width + 63) >> 6;
	if(width != r_coarse_occlusion_width ||
			height != r_coarse_occlusion_height ||
			!r_coarse_occlusion_bits || !r_coarse_occ_scratch){
		free(r_coarse_occlusion_bits);
		free(r_coarse_occ_scratch);
		r_coarse_occlusion_bits = (u64 *)malloc((size_t)words *
				height * sizeof(u64));
		r_coarse_occ_scratch = (f32 *)malloc(4 * (size_t)(height + 2) *
				sizeof(f32));
		if(!r_coarse_occlusion_bits || !r_coarse_occ_scratch){
			free(r_coarse_occlusion_bits);
			free(r_coarse_occ_scratch);
			r_coarse_occlusion_bits = NULL;
			r_coarse_occ_scratch = NULL;
			r_coarse_active = 0;
		}
		r_coarse_occlusion_width = width;
		r_coarse_occlusion_height = height;
		r_coarse_occlusion_words = words;
	}
	r_coarse_occ_vw = (f32)r_refdef.vrect.width;
	r_coarse_occ_vh = (f32)r_refdef.vrect.height;
	r_coarse_poly_surf = NULL;
	if(r_coarse_occlusion_bits)
		memset(r_coarse_occlusion_bits, 0, (size_t)words * height *
				sizeof(u64));
}

void R_CoarseOcclusionClear()
{
	if(r_coarse_occlusion_bits)
		memset(r_coarse_occlusion_bits, 0, (size_t)
			r_coarse_occlusion_words * r_coarse_occlusion_height *
			sizeof(u64));
}


static inline bool R_CoarseOccRunFull(const u64 *row, s32 a, s32 b)
{ // Are all tiles a..b (inclusive) of this row set?
	s32 wa = a >> 6, wb = b >> 6;
	u64 ma = ~(u64)0 << (a & 63);
	u64 mb = ~(u64)0 >> (63 - (b & 63));
	if(wa == wb){
		ma &= mb;
		return (row[wa] & ma) == ma;
	}
	if((row[wa] & ma) != ma)
		return false;
	for(s32 w = wa + 1; w < wb; w++)
		if(row[w] != ~(u64)0)
			return false;
	return (row[wb] & mb) == mb;
}


static inline void R_CoarseOccRunSet(u64 *row, s32 a, s32 b)
{ // Set tiles a..b (inclusive) of this row.
	s32 wa = a >> 6, wb = b >> 6;
	u64 ma = ~(u64)0 << (a & 63);
	u64 mb = ~(u64)0 >> (63 - (b & 63));
	if(wa == wb){
		row[wa] |= ma & mb;
		return;
	}
	row[wa] |= ma;
	for(s32 w = wa + 1; w < wb; w++)
		row[w] = ~(u64)0;
	row[wb] |= mb;
}

// poly is convex, clipped to the viewport and in viewport-relative pixels.
// Returns true when every tile the polygon touches is already FULL.
static bool R_CoarseOccPolyCovered(const coarse_occ_vertex_t *p, s32 n)
{
	const f32 m = COARSE_OCCLUSION_MARGIN;
	const f32 invt = 1.0f / (f32)COARSE_OCCLUSION_TILE_SIZE;
	const f32 tsz = (f32)COARSE_OCCLUSION_TILE_SIZE;
	const s32 tw = r_coarse_occlusion_width;
	const s32 th = r_coarse_occlusion_height;
	const s32 words = r_coarse_occlusion_words;
	const u64 *bits = r_coarse_occlusion_bits;
	f32 *sl = r_coarse_occ_scratch + 2 * (th + 2); // per strip, min x
	f32 *sr = sl + (th + 2); // per strip, max x
	{ // The tile under a vertex is always touched: cheap early out for the
	  // common "visible" case before doing any real work.
		s32 tx = (s32)(p[0].x * invt), ty = (s32)(p[0].y * invt);
		if(tx >= tw) tx = tw - 1;
		if(ty >= th) ty = th - 1;
		if(!((bits[ty * words + (tx >> 6)] >> (tx & 63)) & 1))
			return false;
	}
	f32 ymin = p[0].y, ymax = p[0].y;
	for(s32 i = 1; i < n; i++){
		if(p[i].y < ymin) ymin = p[i].y;
		if(p[i].y > ymax) ymax = p[i].y;
	}
	s32 r0 = (s32)((ymin - m) * invt); // truncation == floor, as it's >= -m
	s32 r1 = (s32)((ymax + m) * invt);
	if(r1 >= th) r1 = th - 1;
	for(s32 r = r0; r <= r1; r++){
		sl[r] = FLT_MAX;
		sr[r] = -FLT_MAX;
	}
	for(s32 i = 0; i < n; i++){
		coarse_occ_vertex_t a = p[i], b = p[i + 1 == n ? 0 : i + 1];
		if(a.y > b.y){
			coarse_occ_vertex_t t = a;
			a = b;
			b = t;
		}
		s32 ra = (s32)((a.y - m) * invt);
		s32 rb = (s32)((b.y + m) * invt);
		if(ra < r0) ra = r0;
		if(rb > r1) rb = r1;
		f32 dy = b.y - a.y;
		if(dy < 1e-5f){ // horizontal: both end points count
			f32 xa = a.x < b.x ? a.x : b.x;
			f32 xb = a.x < b.x ? b.x : a.x;
			for(s32 r = ra; r <= rb; r++){
				if(xa < sl[r]) sl[r] = xa;
				if(xb > sr[r]) sr[r] = xb;
			}
			continue;
		}
		f32 slope = (b.x - a.x) / dy;
		for(s32 r = ra; r <= rb; r++){
			//piece of the edge inside this strip, slightly enlarged
			f32 ylo = (f32)r * tsz - m, yhi = (f32)(r+1) * tsz + m;
			f32 y0 = ylo > a.y ? ylo : a.y;
			f32 y1 = yhi < b.y ? yhi : b.y;
			if(y0 > y1)
				continue;
			f32 x0 = y0 == a.y ? a.x : a.x + (y0 - a.y) * slope;
			f32 x1 = y1 == b.y ? b.x : a.x + (y1 - a.y) * slope;
			if(x0 > x1){
				f32 t = x0;
				x0 = x1;
				x1 = t;
			}
			if(x0 < sl[r]) sl[r] = x0;
			if(x1 > sr[r]) sr[r] = x1;
		}
	}
	bool touched = false;
	for(s32 r = r0; r <= r1; r++){
		if(sl[r] > sr[r])
			continue; // polygon doesn't reach this strip
		s32 a = (s32)((sl[r] - m) * invt);
		s32 b = (s32)((sr[r] + m) * invt);
		if(b >= tw) b = tw - 1;
		if(a > b)
			continue;
		if(!R_CoarseOccRunFull(bits + (size_t)r * words, a, b))
			return false;
		touched = true;
	}
	// Degenerate / off-screen polygons are not treated as occluded.
	return touched;
}

// Marks the tiles that have all four corners inside the convex polygon (same
// preconditions as above), erring on the side of marking fewer.
static void R_CoarseOccPolyFill(const coarse_occ_vertex_t *p, s32 n)
{
	const f32 m = COARSE_OCCLUSION_MARGIN;
	const f32 invt = 1.0f / (f32)COARSE_OCCLUSION_TILE_SIZE;
	const f32 tsz = (f32)COARSE_OCCLUSION_TILE_SIZE;
	const s32 tw = r_coarse_occlusion_width;
	const s32 th = r_coarse_occlusion_height;
	const s32 words = r_coarse_occlusion_words;
	const f32 vw = r_coarse_occ_vw, vh = r_coarse_occ_vh;
	f32 *xl = r_coarse_occ_scratch; // per tile row boundary line
	f32 *xr = xl + (th + 2);
	f32 ymin = p[0].y, ymax = p[0].y;
	for(s32 i = 1; i < n; i++){
		if(p[i].y < ymin) ymin = p[i].y;
		if(p[i].y > ymax) ymax = p[i].y;
	}
	// boundary lines k are at y = min(k * tile, viewport height); those
	// within [ymin, ymax] cross the polygon
	s32 k0 = (s32)ceilf(ymin * invt);
	s32 k1 = (ymax >= vh) ? th : (s32)floorf(ymax * invt);
	if(k1 > th) k1 = th;
	if(k1 - k0 < 1)
		return; // not even one whole tile row
	for(s32 k = k0; k <= k1; k++){
		xl[k] = FLT_MAX;
		xr[k] = -FLT_MAX;
	}
	for(s32 i = 0; i < n; i++){
		coarse_occ_vertex_t a = p[i], b = p[i + 1 == n ? 0 : i + 1];
		// same orientation for both faces sharing an edge
		if(a.y > b.y || (a.y == b.y && a.x > b.x)){
			coarse_occ_vertex_t t = a;
			a = b;
			b = t;
		}
		f32 dy = b.y - a.y;
		s32 ka = (s32)ceilf(a.y * invt);
		if(ka < k0) ka = k0;
		f32 xmin = a.x < b.x ? a.x : b.x;
		f32 xmax = a.x < b.x ? b.x : a.x;
		f32 slope = 0;
		if(dy >= 1e-5f)
			slope = (b.x - a.x) / dy;
		for(s32 k = ka; k <= k1; k++){
			f32 yk = (f32)k * tsz;
			if(yk > vh) yk = vh;
			if(yk > b.y)
				break;
			if(yk < a.y)
				continue;
			f32 x0, x1;
			if(dy == 0){ // lies on the line, covers the whole edge
				x0 = xmin;
				x1 = xmax;
			}else if(dy < 1e-5f){
			    continue; // near-horizontal: skipping only narrows
			}else{
				x0 = yk == a.y ? a.x : yk == b.y ? b.x :
					a.x + (yk - a.y) * slope;
				if(x0 < xmin) x0 = xmin;
				if(x0 > xmax) x0 = xmax;
				x1 = x0;
			}
			if(x0 < xl[k]) xl[k] = x0;
			if(x1 > xr[k]) xr[k] = x1;
		}
	}
	for(s32 r = k0; r < k1; r++){ // tile row r lies between lines r, r+1
		if(xl[r] > xr[r] || xl[r + 1] > xr[r + 1])
			continue;
		f32 l = xl[r] > xl[r + 1] ? xl[r] : xl[r + 1];
		f32 rt = xr[r] < xr[r + 1] ? xr[r] : xr[r + 1];
		// The viewport edges are exact (the polygon was clipped to
		// them), everything else keeps a safety margin.
		s32 ta = l <= 0.0f ? 0 : (s32)ceilf((l + m) * invt);
		s32 tb = rt >= vw ? tw : (s32)((rt - m) * invt);
		if(tb > tw) tb = tw;
		if(ta < tb)
			R_CoarseOccRunSet(r_coarse_occlusion_bits +
				(size_t)r * words, ta, tb - 1);
	}
}

static s32 R_CoarseOccClipVertical(const coarse_occ_vertex_t *in, s32 count,
		coarse_occ_vertex_t *out, f32 clipx, bool keepGreater)
{ // Sutherland-Hodgman polygon clipping against x >= clipx or x <= clipx
	s32 outcount = 0;
	if(count <= 0) return 0;
	for(s32 i = 0; i < count; i++){
		s32 next = (i + 1 == count) ? 0 : i + 1;
		const coarse_occ_vertex_t *a = &in[i];
		const coarse_occ_vertex_t *b = &in[next];
		bool ina = keepGreater ? (a->x >= clipx) : (a->x <= clipx);
		bool inb = keepGreater ? (b->x >= clipx) : (b->x <= clipx);
		if(ina && inb){
			out[outcount++] = *b;
		}else if(ina && !inb){
			f32 dx = b->x - a->x;
			f32 t = (fabsf(dx) > COARSE_OCCLUSION_EPSILON)
				? (clipx - a->x) / dx : 0.0f;
			out[outcount].x = clipx;
			out[outcount].y = a->y + (b->y - a->y) * t;
			outcount++;
		}else if(!ina && inb){
			f32 dx = b->x - a->x;
			f32 t = (fabsf(dx) > COARSE_OCCLUSION_EPSILON)
				? (clipx - a->x) / dx
				: 0.0f;
			out[outcount].x = clipx;
			out[outcount].y = a->y + (b->y - a->y) * t;
			outcount++;
			out[outcount++] = *b;
		}
		if(outcount >= COARSE_OCCLUSION_MAX_VERTS)
			return 0;
	}
	return outcount;
}

static s32 R_CoarseOccClipHorizontal(const coarse_occ_vertex_t *in, s32 count,
                          coarse_occ_vertex_t *out, f32 clipy, bool keepGreater)
{ // Sutherland-Hodgman polygon clipping against y >= clipy or y <= clipy.
	s32 outcount = 0;
	if(count <= 0) return 0;
	for(s32 i = 0; i < count; i++){
		s32 next = (i + 1 == count) ? 0 : i + 1;
		const coarse_occ_vertex_t *a = &in[i];
		const coarse_occ_vertex_t *b = &in[next];
		bool ina = keepGreater ? (a->y >= clipy) : (a->y <= clipy);
		bool inb = keepGreater ? (b->y >= clipy) : (b->y <= clipy);
		if(ina && inb){
			out[outcount++] = *b;
		}else if(ina && !inb){
			f32 dy = b->y - a->y;
			f32 t = (fabsf(dy) > COARSE_OCCLUSION_EPSILON)
				? (clipy - a->y) / dy : 0.0f;
			out[outcount].x = a->x + (b->x - a->x) * t;
			out[outcount].y = clipy;
			outcount++;
		}else if(!ina && inb){
			f32 dy = b->y - a->y;
			f32 t = (fabsf(dy) > COARSE_OCCLUSION_EPSILON)
				? (clipy - a->y) / dy : 0.0f;
			out[outcount].x = a->x + (b->x - a->x) * t;
			out[outcount].y = clipy;
			outcount++;
			out[outcount++] = *b;
		}
		if(outcount >= COARSE_OCCLUSION_MAX_VERTS)
			return 0;
	}
	return outcount;
}

// Build projected, viewport-clipped polygon, in viewport-relative pixels.
// Returns false when we cannot conservatively construct the polygon.
// Important:
// if ANY vertex is behind the near plane, we bail out rather than attempting to
// perform near-plane clipping. This loses optimization opportunities but cannot
// produce a false occlusion result.
static bool R_CoarseOcclusionBuildPolygon(msurface_t *fa,
                              coarse_occ_vertex_t *poly, s32 *outcount)
{
	coarse_occ_vertex_t tmp[COARSE_OCCLUSION_MAX_VERTS];
	if(!fa || fa->numedges < 3) return false;
	if(fa->numedges > COARSE_OCCLUSION_MAX_VERTS) return false;
	model_t *model = currententity->model;
	medge_t *edges = model->edges;
	s32 count = fa->numedges;
	f32 minx = FLT_MAX, maxx = -FLT_MAX, miny = FLT_MAX, maxy = -FLT_MAX;
	for(s32 i = 0; i < count; i++){
		s32 lindex = model->surfedges[fa->firstedge + i];
		mvertex_t *vert;
		vec3_t local;
		vec3_t transformed;
		if(lindex > 0) vert = &r_pcurrentvertbase[edges[lindex].v[0]];
		else vert = &r_pcurrentvertbase[edges[-lindex].v[1]];
		VectorSubtract(vert->position, r_origin, local);
		TransformVector(local, transformed);
		// Do not try to handle near-plane-crossing faces here.
		// R_ClipEdge() remains authoritative for those.
		if(transformed[2] <= NEAR_CLIP) return false;
		f32 zi = 1.0f / transformed[2];
		poly[i].x = xcenter + xscale * zi * transformed[0];
		poly[i].y = ycenter - yscale * zi * transformed[1];
		if(poly[i].x < minx) minx = poly[i].x;
		if(poly[i].x > maxx) maxx = poly[i].x;
		if(poly[i].y < miny) miny = poly[i].y;
		if(poly[i].y > maxy) maxy = poly[i].y;
	}
	// Clip polygon to the actual view rectangle.
	// Use the same viewport represented by r_refdef.vrect.
	f32 left = (f32)r_refdef.vrect.x;
	f32 right = (f32)(r_refdef.vrect.x + r_refdef.vrect.width);
	f32 top = (f32)r_refdef.vrect.y;
	f32 bottom = (f32)(r_refdef.vrect.y + r_refdef.vrect.height);
	// Most faces are small and entirely on screen, don't clip those.
	if(count >= COARSE_OCCLUSION_MAX_VERTS || minx < left || maxx > right
			|| miny < top || maxy > bottom){
		count = R_CoarseOccClipVertical(poly, count, tmp, left, true);
		if(count < 3) return false;
		count = R_CoarseOccClipVertical(tmp, count, poly, right, false);
		if(count < 3) return false;
		count = R_CoarseOccClipHorizontal(poly, count, tmp, top, true);
		if(count < 3) return false;
		count = R_CoarseOccClipHorizontal(tmp, count, poly, bottom,
				false);
		if(count < 3) return false;
	}
	for(s32 i = 0; i < count; i++){
		poly[i].x -= left;
		poly[i].y -= top;
	}
	*outcount = count;
	return true;
}

// Returns true only when the whole projected visible polygon lies inside tiles
// which are already known to be FULL. False is deliberately the default.
bool R_CoarseOcclusionTestSurface(msurface_t *fa)
{
	r_coarse_poly_surf = NULL;
	if(!r_coarse_occlusion_bits) return false;
	// Coarse occlusion is initially restricted to the world.
	// This also avoids entity ordering complications.
	if(currententity != &cl_entities[0]) return false;
	// Never use translucent/sky/cutout geometry as an occluder. The surface
	// test itself could still safely reject against previous opaque
	// coverage, but keeping this restricted makes the initial
	// implementation easier to reason about.
	if(fa->flags & SURF_DRAWSKY) return false;
	if(fa->flags & SURF_DRAWCUTOUT) return false;
	if(fa->flags & SURF_WINQUAKE_DRAWTRANSLUCENT) return false;
	s32 count;
	if(!R_CoarseOcclusionBuildPolygon(fa, r_coarse_poly, &count))
		return false;
	r_coarse_poly_surf = fa; // R_CoarseOcclusionAddSurface can reuse this
	r_coarse_poly_count = count;
	return R_CoarseOccPolyCovered(r_coarse_poly, count);
}

static bool R_CoarseOccluderEligible(msurface_t *fa)
{ // Initial version: only opaque world BSP geometry.
	if(!fa) return 0;
	if(currententity != &cl_entities[0]) return 0;
	if(fa->flags & SURF_DRAWSKY) return 0;
	if(fa->flags & SURF_WINQUAKE_DRAWTRANSLUCENT) return 0;
	if(fa->flags & SURF_DRAWCUTOUT) return 0;
	if(winquake_surface_liquid_alpha < 1.0f) return 0;
	return 1;
}

void R_CoarseOcclusionAddSurface(msurface_t *fa)
{ // Marks only tiles which are COMPLETELY contained by the projected surface
  // polygon. This is intentionally conservative.
	coarse_occ_vertex_t local[COARSE_OCCLUSION_MAX_VERTS];
	coarse_occ_vertex_t *poly;
	s32 count;
	if(!r_coarse_occlusion_bits) return;
	if(!R_CoarseOccluderEligible(fa)) return;
	if(r_coarse_poly_surf == fa){ // built by the test a moment ago
		poly = r_coarse_poly;
		count = r_coarse_poly_count;
	}else{
		poly = local;
		if(!R_CoarseOcclusionBuildPolygon(fa, poly, &count)) return;
	}
	r_coarse_poly_surf = NULL;
	R_CoarseOccPolyFill(poly, count);
}

void R_EmitEdge(mvertex_t *pv0, mvertex_t *pv1)
{
	vec3_t local, transformed;
	f32 *world;
	f32 scale, lzi0, u0, v0;
	s64 ceilv0;
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
		ceilv0 = (s32)ceil(v0);
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
	r_ceilv1 = (s32)ceil(r_v1);
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
	s64 u_check = edge->u; // sort the edge in normally
	if(edge->surfs[0])
		u_check++; // sort trailers after leaders
	// CyanBun96: denser maps like mg1 start begin to chug on the linear
	// edge list traversal as it gets longer. Remembering the last insert
	// point and checking if it fits helps a lot in the best case, and hurts
	// very little in the worst(where performance improvements aren't
	// needed anyway).
	edge_t *pcheck = last_pcheck[v];
	if(!newedges[v] || newedges[v]->u >= u_check){ // Insert at the head
		edge->next = newedges[v];
		newedges[v] = edge;
		last_pcheck[v] = newedges[v]; // cache new head
	}else if(pcheck && pcheck->u < u_check){ // Use cached position
		while(pcheck->next && pcheck->next->u < u_check)
			pcheck = pcheck->next;
		edge->next = pcheck->next;
		pcheck->next = edge;
		last_pcheck[v] = edge; // cache last inserted edge
	}else{ // Fallback to normal head traversal
		pcheck = newedges[v];
		while(pcheck->next && pcheck->next->u < u_check)
			pcheck = pcheck->next;
		edge->next = pcheck->next;
		pcheck->next = edge;
		last_pcheck[v] = edge;
	}
	edge->vlast = (s32)v2;
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

void R_RenderFace(msurface_t *fa, s32 clipflags)
{
	// Manoel Kasimier - skyboxes - begin
	// Code taken from the ToChriS engine - Author: Vic
	// sky surfaces encountered in the world will cause the
	// environment box surfaces to be emited
	if((fa->flags & SURF_DRAWSKY) && skybox_name[0])
		r_skyframe[0] = r_framecount;
	// Manoel Kasimier - skyboxes - end
	if(fa->flags&SURF_WINQUAKE_DRAWTRANSLUCENT && !currententity->alpha){
		winquake_surface_liquid_alpha=R_LiquidAlphaForFlags(fa->flags);
	}else if(cur_ent_alpha < 1 && r_entalpha.value == 1)
		winquake_surface_liquid_alpha = cur_ent_alpha;
	else winquake_surface_liquid_alpha = 1;
	// Baker: Fully transparent = invisible = don't render
	if(!winquake_surface_liquid_alpha)
		return;
// Baker: If this is the alpha water pass and we aren't alpha water, get out!
	if(r_alphapass && winquake_surface_liquid_alpha == 1)
		return;
	if(!r_alphapass && winquake_surface_liquid_alpha < 1){
		r_foundtranswater = 1;
		return;
	}
	if((surface_p) >= surf_max){ // skip out if no more surfs
		r_outofsurfaces++;
		return;
	}
	// ditto if not enough edges left
	if((edge_p + fa->numedges + 4) >= edge_max){
		r_outofedges += fa->numedges;
		return;
	}
	c_faceclip++;
	clipplane_t *pclip = NULL; // set up clip planes
	s32 i = 3;
	u32 mask = 0x08;
	for(; i >= 0; i--, mask >>= 1){
		if(clipflags & mask){
			view_clipplanes[i].next = pclip;
			pclip = &view_clipplanes[i];
		}
	}
	if(r_coarse_active && R_CoarseOcclusionTestSurface(fa)){
		r_coarse_rejected++;
		return;
	}
	r_emitted = 0; // push the edges through
	r_nearzi = 0;
	r_nearzionly = 0;
	makeleftedge = makerightedge = 0;
	medge_t *pedges, tedge;
	pedges = currententity->model->edges;
	r_lastvertvalid = 0;
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
		R_ClipEdge(&r_rightexit, &r_rightenter,view_clipplanes[1].next);
	}
	// if no edges made it out, return without posting the surface
	if(!r_emitted)
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
	vec3_t p_normal;
	TransformVector(pplane->normal, p_normal);
	f32 distinv = 1.0/(pplane->dist-DotProduct(modelorg, pplane->normal));
	surface_p->d_zistepu = p_normal[0] * xscaleinv * distinv;
	surface_p->d_zistepv = -p_normal[1] * yscaleinv * distinv;
	surface_p->d_ziorigin = p_normal[2] * distinv -
		xcenter * surface_p->d_zistepu - ycenter * surface_p->d_zistepv;
	surface_p++;
	if(r_coarse_active)
		R_CoarseOcclusionAddSurface(fa);
}

void R_RenderBmodelFace(bedge_t *pedges, msurface_t *psurf)
{
	if(psurf->flags&SURF_WINQUAKE_DRAWTRANSLUCENT && !currententity->alpha){
	      winquake_surface_liquid_alpha=R_LiquidAlphaForFlags(psurf->flags);
	}else if(cur_ent_alpha < 1 && r_entalpha.value == 1)
		winquake_surface_liquid_alpha = cur_ent_alpha;
	else winquake_surface_liquid_alpha = 1;
	// Baker: Fully transparent = invisible = don't render
	if(!winquake_surface_liquid_alpha)
		return;
// Baker: If this is the alpha water pass and we aren't alpha water, get out!
	if(r_alphapass && winquake_surface_liquid_alpha == 1)
		return;
	if(!r_alphapass && winquake_surface_liquid_alpha < 1){
		r_foundtranswater = 1;
		return;
	}
	if(surface_p >= surf_max){ // skip out if no more surfs
		r_outofsurfaces++;
		return;
	}
	// ditto if not enough edges left
	if((edge_p + psurf->numedges + 4) >= edge_max){
		r_outofedges += psurf->numedges;
		return;
	}
	c_faceclip++;
	// this is a dummy to give the caching mechanism someplace to write to
	r_pedge = &tedge;
	clipplane_t *pclip = NULL; // set up clip planes
	s32 i = 3;
	u32 mask = 0x08;
	for(; i >= 0; i--, mask >>= 1){
		if(r_clipflags & mask){
			view_clipplanes[i].next = pclip;
			pclip = &view_clipplanes[i];
		}
	}
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
	f32 distinv = 1.0/(pplane->dist-DotProduct(modelorg,pplane->normal));
	surface_p->d_zistepu = p_normal[0] * xscaleinv * distinv;
	surface_p->d_zistepv = -p_normal[1] * yscaleinv * distinv;
	surface_p->d_ziorigin = p_normal[2] * distinv -
		xcenter * surface_p->d_zistepu - ycenter * surface_p->d_zistepv;
	surface_p++;
}
