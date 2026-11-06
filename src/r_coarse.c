// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

typedef struct { f32 x; f32 y; } coarse_occ_vertex_t;
// One bit per tile, one row of 64-bit words per tile row.
static u64 *r_coarse_occlusion_bits;
// For every tile that has its bit set: the lowest 1/z that a world face which
// covers the whole tile has inside it, the nearest of those if there are more.
// Only meaningful where the bit is set, so it never needs clearing. Whatever
// has a higher z than this all over the tile is behind a world face there.
static f32 *r_coarse_occlusion_zi;
// The pixels of each tile that opaque world faces cover so far, bit
// (y & 7) * 8 + (x & 7). A tile is FULL once all 64 are covered.
static u64 *r_coarse_occ_pixels;
// The same two for the faces of brush entities, which come after the world
// walk. Only tests with depth look at these: they are kept apart because the
// world faces drawn first would otherwise leave the depth of what's far behind
// the brush entity in the tiles it fills.
static u64 *r_coarse_occ_bpixels;
static f32 *r_coarse_occ_bzi;
static bool r_coarse_occ_bused; // filled in this frame, else to be cleared
static f32 *r_coarse_occ_rows; // 2 arrays of one float per viewport row
static s32 r_coarse_occlusion_width;
static s32 r_coarse_occlusion_height;
static s32 r_coarse_occlusion_words; // u64 words per tile row
static f32 *r_coarse_occ_scratch; // 4 arrays of (height + 2) floats
static f32 r_coarse_occ_vw, r_coarse_occ_vh; // viewport size in pixels
// Polygon built by the occlusion test, reused when the face is added as an
// occluder right afterwards. The clipping can write one vertex past the limit
// before it gives up.
#define COARSE_POLY_SIZE (COARSE_OCCLUSION_MAX_VERTS + 1)
static msurface_t *r_coarse_poly_surf;
static coarse_occ_vertex_t r_coarse_poly[COARSE_POLY_SIZE];
static s32 r_coarse_poly_count;
static f32 r_coarse_poly_zi[2];

// r_coarseocclusion enables conservative coarse screen-space occlusion:
// - Static tile coverage bitmap, one bit per tile.
// - Only FULL tiles are recorded.
// - A tile becomes FULL once every pixel in it is covered by surface
//   polygons, which may take any number of them.
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
// - Fill: the polygon is rasterized the way the edge list does it (pixels are
//   sampled at whole coordinates, left and top edges inclusive), each pixel
//   row is set in the 8 pixel bitmasks of the tiles it crosses.
// COARSE_OCCLUSION_MARGIN makes both operations err on the safe side: the
// test looks at slightly more tiles, the fill marks slightly fewer pixels.

#if COARSE_OCCLUSION_TILE_SIZE != 8
#error "r_coarse_occ_pixels holds 8x8 pixels per tile"
#endif

static void R_CoarseOccClearPixels(u64 *pixels, f32 *zi)
{
	const s32 tw = r_coarse_occlusion_width;
	const s32 th = r_coarse_occlusion_height;
	memset(pixels, 0, (size_t)tw * th * sizeof(u64));
	for(s32 i = 0; i < tw * th; i++)
		zi[i] = FLT_MAX;
	// the pixels of the last tile column and row that are outside the
	// viewport count as covered, so those tiles can still become FULL
	s32 xo = r_refdef.vrect.width & 7, yo = r_refdef.vrect.height & 7;
	if(xo){
		u64 m = 0;
		for(s32 y = 0; y < 8; y++)
			m |= (u64)((0xffu << xo) & 0xffu) << (y * 8);
		for(s32 t = 0; t < th; t++)
			pixels[t * tw + tw - 1] |= m;
	}
	if(yo)
		for(s32 t = 0; t < tw; t++)
			pixels[(th-1) * tw + t] |= ~(u64)0 << (yo*8);
}

// All the arrays in one block, the u64 ones first. NULL bits if it failed.
static void R_CoarseOccAlloc(s32 width, s32 height, s32 words, s32 vh)
{
	static void *block;
	size_t tiles = (size_t)width * height;
	free(block);
	block = malloc(((size_t)words * height + 2 * tiles) * sizeof(u64) +
		(2 * tiles + 4 * (size_t)(height + 2) + 2 * (size_t)vh) *
		sizeof(f32));
	r_coarse_occlusion_bits = block;
	if(!block)
		return;
	r_coarse_occ_pixels = r_coarse_occlusion_bits + (size_t)words * height;
	r_coarse_occ_bpixels = r_coarse_occ_pixels + tiles;
	r_coarse_occlusion_zi = (f32 *)(r_coarse_occ_bpixels + tiles);
	r_coarse_occ_bzi = r_coarse_occlusion_zi + tiles;
	r_coarse_occ_scratch = r_coarse_occ_bzi + tiles;
	r_coarse_occ_rows = r_coarse_occ_scratch + 4 * (size_t)(height + 2);
}

static bool R_OcclusionModeOn(f32 mode)
{ // 1 = on, 2 = on only on maps with a non-vanilla BSP format
	return mode == 1 || (mode && cl_entities[0].model &&
		cl_entities[0].model->bspversion != BSPVERSION);
}

void R_CoarseOcclusionBeginFrame()
{
	r_coarse_rejected = r_nodes_occluded = 0;
	r_bmodels_occluded = r_bmodel_faces_occluded = 0;
	coccl_enable = R_OcclusionModeOn(r_coarseocclusion.value);
	r_coarse_active = coccl_enable && !r_showtris.value;
	r_ent_occ_active = R_OcclusionModeOn(r_entocclusion.value) &&
		!r_showtris.value;
	r_bmodel_occ_active = 0;
	if(!r_coarse_active)
		return;
	s32 width = (r_refdef.vrect.width + COARSE_OCCLUSION_TILE_SIZE - 1) /
		COARSE_OCCLUSION_TILE_SIZE;
	s32 height = (r_refdef.vrect.height + COARSE_OCCLUSION_TILE_SIZE - 1) /
		COARSE_OCCLUSION_TILE_SIZE;
	s32 words = (width + 63) >> 6;
	static s32 r_coarse_occ_rows_n;
	if(width != r_coarse_occlusion_width ||
			height != r_coarse_occlusion_height ||
			r_refdef.vrect.height != r_coarse_occ_rows_n ||
			!r_coarse_occlusion_bits){
		R_CoarseOccAlloc(width, height, words, r_refdef.vrect.height);
		if(!r_coarse_occlusion_bits)
			r_coarse_active = 0;
		r_coarse_occlusion_width = width;
		r_coarse_occlusion_height = height;
		r_coarse_occlusion_words = words;
		r_coarse_occ_rows_n = r_refdef.vrect.height;
	}
	r_coarse_occ_vw = (f32)r_refdef.vrect.width;
	r_coarse_occ_vh = (f32)r_refdef.vrect.height;
	r_coarse_poly_surf = NULL;
	if(r_coarse_occlusion_bits){
		memset(r_coarse_occlusion_bits, 0, (size_t)words * height *
				sizeof(u64));
		R_CoarseOccClearPixels(r_coarse_occ_pixels,
				r_coarse_occlusion_zi);
	}
	r_coarse_occ_bused = 0;
	// r_draworder shows the highest key on top instead of the lowest, which
	// is the opposite of what the brush entity tests rely on
	r_bmodel_occ_active = r_coarse_active && !r_draworder.value &&
		R_OcclusionModeOn(r_bmodelocclusion.value);
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

// Covers pixels x0..x1 (inclusive) of viewport row y. Tiles that get new pixels
// from it take the lowest 1/z that the face (zi plane of surf, in screen
// pixels, and minzi or more everywhere) has in them, if it's the lowest yet.
// brush picks the layer of brush entity faces, which has no FULL bits.
static inline void R_CoarseOccRowFill(s32 y, s32 x0, s32 x1,
		const surf_t *surf, f32 minzi, bool brush)
{
	const s32 tw = r_coarse_occlusion_width;
	const s32 ty = y >> 3;
	u64 *pix = (brush ? r_coarse_occ_bpixels : r_coarse_occ_pixels)
		+ (size_t)ty * tw;
	f32 *zirow = (brush ? r_coarse_occ_bzi : r_coarse_occlusion_zi)
		+ (size_t)ty * tw;
	u64 *bits = r_coarse_occlusion_bits+(size_t)ty*r_coarse_occlusion_words;
	const u32 shift = (u32)(y & 7) * 8;
	// 1/z is linear across the screen, so its minimum over the pixels of a
	// tile is at one of the corner pixels, the u and v parts can be taken
	// separately
	f32 sy = (f32)(r_refdef.vrect.y + (ty << 3));
	f32 zv0 = surf->d_zistepv * sy, zv1 = surf->d_zistepv * (sy + 7.0f);
	f32 zv = surf->d_ziorigin + (zv0 < zv1 ? zv0 : zv1);
	f32 du = surf->d_zistepu * 7.0f;
	s32 t0 = x0 >> 3, t1 = x1 >> 3;
	for(s32 t = t0; t <= t1; t++){
		u32 lo = t == t0 ? (u32)(x0 & 7) : 0;
		u32 hi = t == t1 ? (u32)(x1 & 7) : 7;
		u64 row = (u64)((0xffu >> (7 - hi)) & (0xffu << lo)) << shift;
		u64 fresh = row & ~pix[t];
		if(!fresh)
			continue; // nearer faces got there first
		pix[t] |= fresh;
		f32 zu = surf->d_zistepu * (f32)(r_refdef.vrect.x + (t << 3));
		f32 zi = zv + (du < 0 ? zu + du : zu);
		if(zi < minzi) zi = minzi;
		if(zi < zirow[t])
			zirow[t] = zi;
		if(!brush && pix[t] == ~(u64)0)
			bits[t >> 6] |= (u64)1 << (t & 63);
	}
}

// ceilf and floorf without the libm call
Q_SIL s32 R_CoarseOccCeil(f32 x)  { s32 i = (s32)x; return i + ((f32)i < x); }
Q_SIL s32 R_CoarseOccFloor(f32 x) { s32 i = (s32)x; return i - ((f32)i > x); }

// Are all tiles a..b (inclusive, all known to be set) covered by something
// that is in front of 1/z value zi everywhere in them?
static inline bool R_CoarseOccZiInFront(const f32 *zirow, s32 a, s32 b, f32 zi)
{
	zi *= COARSE_OCCLUSION_ZI_MARGIN;
	for(s32 t = a; t <= b; t++)
		if(zirow[t] <= zi)
			return false;
	return true;
}

// R_CoarseOccRunCovered for tiles the world faces don't cover by themselves
static bool R_CoarseOccRunUnion(s32 r, s32 a, s32 b, f32 zi)
{
	const s32 tw = r_coarse_occlusion_width;
	const f32 *zirow = r_coarse_occlusion_zi + (size_t)r * tw;
	zi *= COARSE_OCCLUSION_ZI_MARGIN;
	const u64 *pix = r_coarse_occ_pixels + (size_t)r * tw;
	const u64 *bpix = r_coarse_occ_bpixels + (size_t)r * tw;
	const f32 *bzi = r_coarse_occ_bzi + (size_t)r * tw;
	for(s32 t = a; t <= b; t++)
		if(((zirow[t] > zi ? pix[t] : 0) | (bzi[t] > zi ? bpix[t] : 0))
				!= ~(u64)0)
			return false;
	return true;
}

// Are tiles a..b (inclusive) of tile row r FULL, and with depth set, also
// covered by faces in front of 1/z value zi? With depth the brush entity faces
// count too, and a tile can be covered by world and brush faces together.
Q_SIL bool R_CoarseOccRunCovered(s32 r, s32 a, s32 b, bool depth, f32 zi)
{
	if(R_CoarseOccRunFull(r_coarse_occlusion_bits +
			(size_t)r * r_coarse_occlusion_words, a, b) &&
			(!depth || R_CoarseOccZiInFront(r_coarse_occlusion_zi +
			(size_t)r * r_coarse_occlusion_width, a, b, zi)))
		return true;
	return depth && r_coarse_occ_bused && R_CoarseOccRunUnion(r, a, b, zi);
}

// poly is convex, clipped to the viewport and in viewport-relative pixels.
// Returns true when every tile the polygon touches is already FULL. With depth
// set, the tiles must also be covered by faces that are in front of 1/z value
// zi, the nearest the tested polygon gets.
static bool R_CoarseOccPolyCovered(const coarse_occ_vertex_t *p, s32 n,
		bool depth, f32 zi)
{
	const f32 m = COARSE_OCCLUSION_MARGIN;
	const f32 invt = 1.0f / (f32)COARSE_OCCLUSION_TILE_SIZE;
	const f32 tsz = (f32)COARSE_OCCLUSION_TILE_SIZE;
	const s32 tw = r_coarse_occlusion_width;
	const s32 th = r_coarse_occlusion_height;
	f32 *sl = r_coarse_occ_scratch + 2 * (th + 2); // per strip, min x
	f32 *sr = sl + (th + 2); // per strip, max x
	{ // The tile under a vertex is always touched: cheap early out for the
	  // common "visible" case before doing any real work.
		s32 tx = (s32)(p[0].x * invt), ty = (s32)(p[0].y * invt);
		if(tx >= tw) tx = tw - 1;
		if(ty >= th) ty = th - 1;
		if(!R_CoarseOccRunCovered(ty, tx, tx, depth, zi))
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
		if(!R_CoarseOccRunCovered(r, a, b, depth, zi))
			return false;
		touched = true;
	}
	// Degenerate / off-screen polygons are not treated as occluded.
	return touched;
}

// Covers the pixels of the convex polygon (same preconditions as above) that
// the edge list would draw it on, erring on the side of covering fewer. surf is
// the face in the edge list and minzi the lowest 1/z of the polygon, for depth.
static void R_CoarseOccPolyFill(const coarse_occ_vertex_t *p, s32 n,
		const surf_t *surf, f32 minzi, bool brush)
{
	const f32 m = COARSE_OCCLUSION_MARGIN;
	const s32 vw = r_refdef.vrect.width, vh = r_refdef.vrect.height;
	f32 *xl = r_coarse_occ_rows; // per pixel row
	f32 *xr = xl + vh;
	f32 ymin = p[0].y, ymax = p[0].y;
	for(s32 i = 1; i < n; i++){
		if(p[i].y < ymin) ymin = p[i].y;
		if(p[i].y > ymax) ymax = p[i].y;
	}
	// rows from ceil(top) to ceil(bottom) - 1, like R_EmitEdge. The
	// viewport edges are exact (the polygon was clipped to them),
	// everything else keeps a safety margin.
	s32 y0 = ymin <= 0.0f ? 0 : R_CoarseOccCeil(ymin + m);
	s32 y1 = R_CoarseOccCeil(ymax - m) - 1;
	if(y1 >= vh) y1 = vh - 1;
	if(y0 > y1)
		return;
	for(s32 y = y0; y <= y1; y++){
		xl[y] = FLT_MAX;
		xr[y] = -FLT_MAX;
	}
	for(s32 i = 0; i < n; i++){
		coarse_occ_vertex_t a = p[i], b = p[i + 1 == n ? 0 : i + 1];
		if(a.y > b.y){
			coarse_occ_vertex_t t = a;
			a = b;
			b = t;
		}
		f32 dy = b.y - a.y;
		if(dy < 1e-5f)
			continue; // the edges next to it reach both its ends
		f32 slope = (b.x - a.x) / dy;
		f32 xmin = a.x < b.x ? a.x : b.x;
		f32 xmax = a.x < b.x ? b.x : a.x;
		s32 ya = R_CoarseOccCeil(a.y), yb = R_CoarseOccFloor(b.y);
		if(ya < y0) ya = y0;
		if(yb > y1) yb = y1;
		for(s32 y = ya; y <= yb; y++){
			f32 x = a.x + ((f32)y - a.y) * slope;
			if(x < xmin) x = xmin;
			if(x > xmax) x = xmax;
			if(x < xl[y]) xl[y] = x;
			if(x > xr[y]) xr[y] = x;
		}
	}
	for(s32 y = y0; y <= y1; y++){
		if(xl[y] > xr[y])
			continue;
		// columns from ceil(left) to ceil(right) - 1
		s32 x0 = xl[y] <= 0.0f ? 0 : R_CoarseOccCeil(xl[y] + m);
		s32 x1 = R_CoarseOccCeil(xr[y] - m) - 1;
		if(xr[y] >= (f32)vw) x1 = vw - 1;
		if(x1 >= vw) x1 = vw - 1;
		if(x0 <= x1)
			R_CoarseOccRowFill(y, x0, x1, surf, minzi, brush);
	}
}

// Sutherland-Hodgman polygon clipping against x (vertical) or y >= clip, or
// <= clip without keepgreater
static s32 R_CoarseOccClip(const coarse_occ_vertex_t *in, s32 count,
		coarse_occ_vertex_t *out, bool vertical, f32 clip,
		bool keepgreater)
{
	s32 outcount = 0;
	for(s32 i = 0; i < count; i++){
		coarse_occ_vertex_t a = in[i];
		coarse_occ_vertex_t b = in[i + 1 == count ? 0 : i + 1];
		f32 ca = vertical ? a.x : a.y, cb = vertical ? b.x : b.y;
		bool ina = keepgreater ? ca >= clip : ca <= clip;
		bool inb = keepgreater ? cb >= clip : cb <= clip;
		if(ina != inb){ // the point where it crosses the line
			f32 d = cb - ca;
			f32 t = fabsf(d) > COARSE_OCCLUSION_EPSILON ?
				(clip - ca) / d : 0.0f;
			coarse_occ_vertex_t *o = &out[outcount++];
			o->x = vertical ? clip : a.x + (b.x - a.x) * t;
			o->y = vertical ? a.y + (b.y - a.y) * t : clip;
		}
		if(inb)
			out[outcount++] = b;
		if(outcount >= COARSE_OCCLUSION_MAX_VERTS)
			return 0;
	}
	return outcount;
}

// Build projected, viewport-clipped polygon, in viewport-relative pixels.
// Returns false when we cannot conservatively construct the polygon. If zi
// isn't NULL it gets the lowest and highest 1/z of the vertices, the farthest
// and the nearest point.
// Important:
// if ANY vertex is behind the near plane, we bail out rather than attempting to
// perform near-plane clipping. This loses optimization opportunities but cannot
// produce a false occlusion result.
static bool R_CoarseOcclusionBuildPolygon(msurface_t *fa,
		coarse_occ_vertex_t *poly, s32 *outcount, f32 *zi)
{
	coarse_occ_vertex_t tmp[COARSE_POLY_SIZE];
	if(!fa || fa->numedges < 3 ||
			fa->numedges > COARSE_OCCLUSION_MAX_VERTS)
		return false;
	model_t *model = currententity->model;
	medge_t *edges = model->edges;
	// R_EmitEdge projects from modelorg, which is the player's position in
	// the space of the entity being drawn (r_origin for the world)
	const f32 *vieworg = currententity == &cl_entities[0]?r_origin:modelorg;
	s32 count = fa->numedges;
	f32 minx = FLT_MAX, maxx = -FLT_MAX, miny = FLT_MAX, maxy = -FLT_MAX;
	f32 nearzi = 0, farzi = FLT_MAX;
	for(s32 i = 0; i < count; i++){
		s32 lindex = model->surfedges[fa->firstedge + i];
		vec3_t local;
		vec3_t transformed;
		// a select, the direction of an edge is a coin toss
		bool back = lindex <= 0;
		mvertex_t *vert = &r_pcurrentvertbase[
			edges[back ? -lindex : lindex].v[back]];
		VectorSubtract(vert->position, vieworg, local);
		TransformVector(local, transformed);
		// Do not try to handle near-plane-crossing faces here.
		// R_ClipEdge() remains authoritative for those.
		if(transformed[2] <= NEAR_CLIP) return false;
		f32 vzi = 1.0f / transformed[2];
		if(vzi > nearzi) nearzi = vzi;
		if(vzi < farzi) farzi = vzi;
		poly[i].x = xcenter + xscale * vzi * transformed[0];
		poly[i].y = ycenter - yscale * vzi * transformed[1];
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
		count = R_CoarseOccClip(poly, count, tmp, 1, left, 1);
		if(count < 3) return false;
		count = R_CoarseOccClip(tmp, count, poly, 1, right, 0);
		if(count < 3) return false;
		count = R_CoarseOccClip(poly, count, tmp, 0, top, 1);
		if(count < 3) return false;
		count = R_CoarseOccClip(tmp, count, poly, 0, bottom, 0);
		if(count < 3) return false;
	}
	for(s32 i = 0; i < count; i++){
		poly[i].x -= left;
		poly[i].y -= top;
	}
	*outcount = count;
	if(zi){
		zi[0] = farzi;
		zi[1] = nearzi;
	}
	return true;
}

// Builds the polygon of a face and tests it, leaving it for
// R_CoarseOcclusionAddSurface
static bool R_CoarseOccTestFace(msurface_t *fa, bool depth)
{
	if(!R_CoarseOcclusionBuildPolygon(fa, r_coarse_poly,
				&r_coarse_poly_count, r_coarse_poly_zi))
		return false;
	r_coarse_poly_surf = fa;
	return R_CoarseOccPolyCovered(r_coarse_poly, r_coarse_poly_count,
			depth, r_coarse_poly_zi[1]);
}

// Returns true only when the whole projected visible polygon lies inside tiles
// which are already known to be FULL. False is deliberately the default.
bool R_CoarseOcclusionTestSurface(msurface_t *fa)
{
	// brush entities have R_CoarseOcclusionTestBmodelFace, and may have
	// left their polygon for R_CoarseOcclusionAddSurface
	if(currententity != &cl_entities[0]) return false;
	r_coarse_poly_surf = NULL;
	if(!r_coarse_occlusion_bits) return false;
	// Translucent/sky/cutout geometry is never used as an occluder (see
	// R_CoarseOccluderEligible), but cutout and translucent faces can
	// still be rejected here. In the main pass everything comes in BSP
	// order, so whatever covered the tiles so far is in front. The alpha
	// pass comes after the whole main pass though, so there the tiles are
	// also covered by what is behind the face, and only the depth of the
	// opaque faces that covered them says whether it's hidden (the
	// z-buffer test that alpha surfaces get).
	if(fa->flags & SURF_DRAWSKY) return false;
	return R_CoarseOccTestFace(fa, r_alphapass);
}

// Opaque world faces, and in the main pass the opaque faces of brush entities.
// Those come after the world walk, when only tests with depth are left, so it
// doesn't matter that they don't come in BSP order.
static bool R_CoarseOccluderEligible(msurface_t *fa)
{
	if(!fa) return 0;
	if(currententity != &cl_entities[0] &&
			(r_alphapass || !r_bmodel_occ_faces))
		return 0;
	if(fa->flags & SURF_DRAWSKY) return 0;
	// the skybox is drawn close by, but sorts behind everything else
	if(fa->flags & SURF_DRAWSKYBOX) return 0;
	if(fa->flags & SURF_WINQUAKE_DRAWTRANSLUCENT) return 0;
	if(fa->flags & SURF_DRAWCUTOUT) return 0;
	if(winquake_surface_liquid_alpha < 1.0f) return 0;
	return 1;
}

void R_CoarseOcclusionAddSurface(msurface_t *fa, const surf_t *surf)
{ // Covers the pixels of the projected surface polygon, a tile is FULL once
  // they're all covered. This is intentionally conservative. surf is the one
  // the face got in the edge list, its 1/z is remembered for the brush entity
  // tests.
	coarse_occ_vertex_t local[COARSE_POLY_SIZE];
	coarse_occ_vertex_t *poly;
	s32 count;
	f32 localzi[2], *zi;
	if(!r_coarse_occlusion_bits) return;
	if(!R_CoarseOccluderEligible(fa)) return;
	if(r_coarse_poly_surf == fa){ // built by the test a moment ago
		poly = r_coarse_poly;
		count = r_coarse_poly_count;
		zi = r_coarse_poly_zi;
	}else{
		poly = local;
		zi = localzi;
		if(!R_CoarseOcclusionBuildPolygon(fa, poly, &count, zi))
			return;
	}
	r_coarse_poly_surf = NULL;
	bool brush = currententity != &cl_entities[0];
	if(brush && !r_coarse_occ_bused){
		R_CoarseOccClearPixels(r_coarse_occ_bpixels, r_coarse_occ_bzi);
		r_coarse_occ_bused = 1;
	}
	R_CoarseOccPolyFill(poly, count, surf, zi[0], brush);
}

// Brush entities go into the edge list with the key of the BSP leaf they are
// in, and the edge scanner shows whichever surface has the lowest key. Along
// any line of sight the BSP order is the depth order, so where a world face is
// nearer than a brush face, its key is lower too and the brush face can't show
// up there. The alpha pass draws its surfaces over the finished frame instead,
// testing every pixel against the z-buffer, so there it's the depth test that
// says a face shows up nowhere.
bool R_CoarseOcclusionTestBmodelFace(msurface_t *fa)
{
	r_coarse_poly_surf = NULL;
	if(!r_bmodel_occ_active || !r_coarse_occlusion_bits) return false;
	if(!r_alphapass && cur_ent_alpha < 1.0f) return false;
	if(currententity == &cl_entities[0]) return false;
	return R_CoarseOccTestFace(fa, 1);
}

// Is the tile under viewport point x, y FULL (and in front of zi, with depth)?
// Yes for points off screen, those don't decide anything.
static inline bool R_CoarseOccPointFull(f32 x, f32 y, bool depth, f32 zi)
{
	if(!(x >= 0 && x < r_coarse_occ_vw && y >= 0 && y < r_coarse_occ_vh))
		return true;
	s32 tx = (s32)x >> 3, ty = (s32)y >> 3;
	return R_CoarseOccRunCovered(ty, tx, tx, depth, zi);
}

// Are all the tiles a world space box touches FULL, as seen from the player's
// point of view? With depth set, they must also be covered by something nearer
// than all of the box.
static bool R_CoarseOccBoxCovered(const f32 *mn, const f32 *mx, bool depth)
{
	const f32 m = COARSE_OCCLUSION_MARGIN;
	const f32 invt = 1.0f / (f32)COARSE_OCCLUSION_TILE_SIZE;
	const s32 tw = r_coarse_occlusion_width, th = r_coarse_occlusion_height;
	// The center is inside the box, so the tile it's on is touched: a
	// cheap early out for the common "visible" case, as in PolyCovered.
	// It only ever says no, so it can round differently from the rest.
	vec3_t c, t;
	for(s32 k = 0; k < 3; k++)
		c[k] = (mn[k] + mx[k]) * 0.5f - r_origin[k];
	TransformVector(c, t);
	if(t[2] <= NEAR_CLIP) return false; // reaches the near plane
	f32 zi = 1.0f / t[2];
	f32 x = xcenter + xscale * zi * t[0] - (f32)r_refdef.vrect.x;
	f32 y = ycenter - yscale * zi * t[1] - (f32)r_refdef.vrect.y;
	// the box's nearest point is at least as near
	if(!R_CoarseOccPointFull(x, y, depth, zi))
		return false;
	// So is the disk across the view through the center, with the radius
	// of the box's smallest half extent: four more points, no division.
	// Most boxes get past the center and then fail the full test below.
	f32 r = mx[0] - mn[0];
	r = mx[1] - mn[1] < r ? mx[1] - mn[1] : r;
	r = mx[2] - mn[2] < r ? mx[2] - mn[2] : r;
	f32 rx = xscale * zi * r * 0.5f, ry = yscale * zi * r * 0.5f;
	if(!R_CoarseOccPointFull(x - rx, y, depth, zi) ||
			!R_CoarseOccPointFull(x + rx, y, depth, zi) ||
			!R_CoarseOccPointFull(x, y - ry, depth, zi) ||
			!R_CoarseOccPointFull(x, y + ry, depth, zi))
		return false;
	// The view space corners are the transformed min corner plus some of
	// the transformed box edges along each axis
	vec3_t base, ax[3];
	VectorSubtract(mn, r_origin, t);
	TransformVector(t, base);
	f32 minz = base[2];
	for(s32 k = 0; k < 3; k++){
		f32 d = mx[k] - mn[k];
		ax[k][0] = d * vright[k];
		ax[k][1] = d * vup[k];
		ax[k][2] = d * vpn[k];
		if(ax[k][2] < 0) minz += ax[k][2];
	}
	if(minz <= NEAR_CLIP) return false; // reaches the near plane
	f32 cx[8], cy[8], cz[8];
	for(s32 i = 0; i < 8; i++){
		cx[i] = base[0], cy[i] = base[1], cz[i] = base[2];
		for(s32 k = 0; k < 3; k++)
			if((i >> k) & 1){
				cx[i] += ax[k][0];
				cy[i] += ax[k][1];
				cz[i] += ax[k][2];
			}
	}
	f32 minx = FLT_MAX, maxx = -FLT_MAX, miny = FLT_MAX, maxy = -FLT_MAX;
	f32 nearzi = 0;
	for(s32 i = 0; i < 8; i++){
		f32 zi = 1.0f / cz[i];
		f32 x = xcenter + xscale * zi * cx[i];
		f32 y = ycenter - yscale * zi * cy[i];
		nearzi = zi > nearzi ? zi : nearzi;
		minx = x < minx ? x : minx;
		maxx = x > maxx ? x : maxx;
		miny = y < miny ? y : miny;
		maxy = y > maxy ? y : maxy;
	}
	minx -= (f32)r_refdef.vrect.x; maxx -= (f32)r_refdef.vrect.x;
	miny -= (f32)r_refdef.vrect.y; maxy -= (f32)r_refdef.vrect.y;
	if(!(maxx >= 0 && minx <= r_coarse_occ_vw && maxy >= 0 &&
			miny <= r_coarse_occ_vh))
		return false; // not on screen at all, nothing to say about it
	if(minx < 0) minx = 0;
	if(miny < 0) miny = 0;
	if(maxx > r_coarse_occ_vw) maxx = r_coarse_occ_vw;
	if(maxy > r_coarse_occ_vh) maxy = r_coarse_occ_vh;
	s32 tx0 = (s32)((minx - m) * invt), tx1 = (s32)((maxx + m) * invt);
	s32 ty0 = (s32)((miny - m) * invt), ty1 = (s32)((maxy + m) * invt);
	if(tx1 >= tw) tx1 = tw - 1;
	if(ty1 >= th) ty1 = th - 1;
	for(s32 r = ty0; r <= ty1; r++)
		if(!R_CoarseOccRunCovered(r, tx0, tx1, depth, nearzi))
			return false;
	return true;
}

// The same as for brush faces, for everything inside a world space box (so
// before the entity's own rotation is applied).
bool R_CoarseOcclusionTestBox(const f32 *mn, const f32 *mx)
{
	if(!r_bmodel_occ_active || !r_coarse_occlusion_bits) return false;
	if(!r_alphapass && cur_ent_alpha < 1.0f) return false;
	return R_CoarseOccBoxCovered(mn, mx, 1);
}

// Can nothing in this world BSP node or leaf show up? Everything in it comes
// after the world faces drawn so far in the BSP order, so like the faces in
// R_CoarseOcclusionTestSurface, it's hidden wherever the tiles are FULL, and
// in the alpha pass only if it's also behind the faces that filled them.
bool R_CoarseOcclusionTestNode(mnode_t *node)
{
	if(!r_coarse_active || !r_coarse_occlusion_bits) return false;
	return R_CoarseOccBoxCovered(node->minmaxs, node->minmaxs + 3,
			r_alphapass);
}
