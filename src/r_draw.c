// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

static u8 *r_coarse_occlusion_bits;
static s32 r_coarse_occlusion_width;
static s32 r_coarse_occlusion_height;
static s32 r_coarse_occlusion_stride;
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

// r_coarseocclusion enables conservative coarse screen-space occlusion:
// - Static tile coverage bitmap.
// - Only FULL tiles are recorded.
// - A tile becomes FULL only if all four of its corners are inside the
//   projected convex surface polygon.
// - A surface is rejected only if every tile intersected by its projected
//   polygon is already FULL.
// - Faces with any vertex behind NEAR_CLIP are not considered for coarse
//   occlusion at all. This is conservative around the near plane.
// - Screen-space clipping is performed against r_refdef.vrect.

void R_CoarseOcclusionBeginFrame()
{ // Call once per rendered frame after r_refdef.vrect is valid.
	s32 width, height, bytes;
	r_coarse_rejected = 0;
	width = (r_refdef.vrect.width + COARSE_OCCLUSION_TILE_SIZE - 1) /
		COARSE_OCCLUSION_TILE_SIZE;
	height = (r_refdef.vrect.height + COARSE_OCCLUSION_TILE_SIZE - 1) /
		COARSE_OCCLUSION_TILE_SIZE;
	bytes = width * height;
	if (width != r_coarse_occlusion_width ||
			height != r_coarse_occlusion_height ||
			!r_coarse_occlusion_bits) {
		free(r_coarse_occlusion_bits);
		r_coarse_occlusion_bits = (u8 *)malloc(bytes);
		r_coarse_occlusion_width = width;
		r_coarse_occlusion_height = height;
		r_coarse_occlusion_stride = width;
	}
	if (r_coarse_occlusion_bits)
		memset(r_coarse_occlusion_bits, 0, bytes);
}

void R_CoarseOcclusionClear()
{
	if(r_coarse_occlusion_bits){
		memset(r_coarse_occlusion_bits, 0,
				r_coarse_occlusion_width *
				r_coarse_occlusion_height);
	}
}

static inline bool R_CoarseOcclusionTileFull(s32 tx, s32 ty)
{
	if ((u32)tx >= (u32)r_coarse_occlusion_width ||
			(u32)ty >= (u32)r_coarse_occlusion_height)
		return false;
	return r_coarse_occlusion_bits[
		ty * r_coarse_occlusion_stride + tx] != 0;
}


static inline void R_CoarseOcclusionSetTile(s32 tx, s32 ty)
{
	if ((u32)tx >= (u32)r_coarse_occlusion_width ||
			(u32)ty >= (u32)r_coarse_occlusion_height)
		return;
	r_coarse_occlusion_bits[ ty * r_coarse_occlusion_stride + tx] = 1;
}

static inline f32 R_CoarseOccCross(const coarse_occ_vertex_t *a,
	const coarse_occ_vertex_t *b, f32 x, f32 y) // 2D cross product: AB x AC
{ return (b->x - a->x) * (y - a->y) - (b->y - a->y) * (x - a->x); }

static bool R_CoarseOccPointInPolygon(const coarse_occ_vertex_t *poly,
                          s32 count, f32 x, f32 y)
{ // Convex polygon, arbitrary winding.
	s32 i;
	f32 sign = 0.0f;
	if (count < 3) return false;
	for (i = 0; i < count; i++) {
		s32 next = (i + 1 == count) ? 0 : i + 1;
		f32 cross = R_CoarseOccCross(&poly[i], &poly[next], x, y);
		if (fabsf(cross) <= COARSE_OCCLUSION_EPSILON)
			continue;
		if (sign == 0.0f)
			sign = cross;
		else if ((sign > 0.0f && cross < 0.0f) ||
				(sign < 0.0f && cross > 0.0f))
			return false;
	}
	return true;
}

static bool R_CoarseOccOnSegment(f32 ax, f32 ay, f32 bx, f32 by, f32 px, f32 py)
{
	if (px < fminf(ax, bx) - COARSE_OCCLUSION_EPSILON ||
		px > fmaxf(ax, bx) + COARSE_OCCLUSION_EPSILON ||
		py < fminf(ay, by) - COARSE_OCCLUSION_EPSILON ||
		py > fmaxf(ay, by) + COARSE_OCCLUSION_EPSILON) return false;
	return true;
}

static bool R_CoarseOccSegmentsIntersect(const coarse_occ_vertex_t *a,
                             const coarse_occ_vertex_t *b,
                             const coarse_occ_vertex_t *c,
                             const coarse_occ_vertex_t *d)
{
	f32 c1, c2, c3, c4;
	c1 = (b->x - a->x) * (c->y - a->y) - (b->y - a->y) * (c->x - a->x);
	c2 = (b->x - a->x) * (d->y - a->y) - (b->y - a->y) * (d->x - a->x);
	c3 = (d->x - c->x) * (a->y - c->y) - (d->y - c->y) * (a->x - c->x);
	c4 = (d->x - c->x) * (b->y - c->y) - (d->y - c->y) * (b->x - c->x);
	if(((c1 > COARSE_OCCLUSION_EPSILON && c2 < -COARSE_OCCLUSION_EPSILON) ||
	   (c1 < -COARSE_OCCLUSION_EPSILON && c2 > COARSE_OCCLUSION_EPSILON)) &&
	   ((c3 > COARSE_OCCLUSION_EPSILON && c4 < -COARSE_OCCLUSION_EPSILON) ||
	   (c3 < -COARSE_OCCLUSION_EPSILON && c4 > COARSE_OCCLUSION_EPSILON)))
		return true;
	if (fabsf(c1) <= COARSE_OCCLUSION_EPSILON &&
		R_CoarseOccOnSegment(a->x, a->y, b->x, b->y, c->x, c->y))
		return true;
	if (fabsf(c2) <= COARSE_OCCLUSION_EPSILON &&
		R_CoarseOccOnSegment(a->x, a->y, b->x, b->y, d->x, d->y))
		return true;
	if (fabsf(c3) <= COARSE_OCCLUSION_EPSILON &&
		R_CoarseOccOnSegment(c->x, c->y, d->x, d->y, a->x, a->y))
		return true;
	if (fabsf(c4) <= COARSE_OCCLUSION_EPSILON &&
		R_CoarseOccOnSegment(c->x, c->y, d->x, d->y, b->x, b->y))
		return true;
	return false;
}

// Exact conservative intersection test between a convex polygon
// and an axis-aligned rectangle.
// False means definitely disjoint.
// True means the polygon and rectangle overlap/touch.
static bool R_CoarseOccPolyIntersectsRect(const coarse_occ_vertex_t *poly,
                              s32 count, f32 x0, f32 y0, f32 x1, f32 y1)
{
	coarse_occ_vertex_t rect[4];
	s32 i;
	rect[0].x = x0; rect[0].y = y0;
	rect[1].x = x1; rect[1].y = y0;
	rect[2].x = x1; rect[2].y = y1;
	rect[3].x = x0; rect[3].y = y1;
	for (i = 0; i < count; i++) { // Polygon vertex inside rectangle.
		if (poly[i].x >= x0 && poly[i].x <= x1 &&
			poly[i].y >= y0 && poly[i].y <= y1) return true;
	}
	for (i = 0; i < 4; i++) // Rectangle corner inside polygon.
		if(R_CoarseOccPointInPolygon(poly, count, rect[i].x, rect[i].y))
			return true;
	for (i = 0; i < count; i++) { // Polygon edge intersects rectangle edge.
		s32 next = (i + 1 == count) ? 0 : i + 1;
		for (s32 j = 0; j < 4; j++) {
			s32 jnext = (j + 1 == 4) ? 0 : j + 1;
			if (R_CoarseOccSegmentsIntersect( &poly[i], &poly[next],
						&rect[j], &rect[jnext]))
				return true;
		}
	}
	return false;
}

// The tile is FULL only if all four of its corners are inside the convex
// polygon. This is intentionally stronger than polygon/tile intersection.
static bool R_CoarseOccTileInsidePolygon(const coarse_occ_vertex_t *poly,
                             s32 count, f32 x0, f32 y0, f32 x1, f32 y1)
{
	if (!R_CoarseOccPointInPolygon(poly, count, x0, y0)) return false;
	if (!R_CoarseOccPointInPolygon(poly, count, x1, y0)) return false;
	if (!R_CoarseOccPointInPolygon(poly, count, x1, y1)) return false;
	if (!R_CoarseOccPointInPolygon(poly, count, x0, y1)) return false;
	return true;
}

static s32 R_CoarseOccClipVertical(const coarse_occ_vertex_t *in,
                        s32 count, coarse_occ_vertex_t *out, f32 clipx,
                        bool keepGreater)
{ // Sutherland-Hodgman polygon clipping against x >= clipx or x <= clipx
	s32 i;
	s32 outcount = 0;
	if (count <= 0) return 0;
	for (i = 0; i < count; i++) {
		s32 next = (i + 1 == count) ? 0 : i + 1;
		const coarse_occ_vertex_t *a = &in[i];
		const coarse_occ_vertex_t *b = &in[next];
		bool ina = keepGreater ? (a->x >= clipx) : (a->x <= clipx);
		bool inb = keepGreater ? (b->x >= clipx) : (b->x <= clipx);
		if (ina && inb) {
			out[outcount++] = *b;
		} else if (ina && !inb) {
			f32 dx = b->x - a->x;
			f32 t = (fabsf(dx) > COARSE_OCCLUSION_EPSILON)
				? (clipx - a->x) / dx : 0.0f;
			out[outcount].x = clipx;
			out[outcount].y = a->y + (b->y - a->y) * t;
			outcount++;
		} else if (!ina && inb) {
			f32 dx = b->x - a->x;
			f32 t = (fabsf(dx) > COARSE_OCCLUSION_EPSILON)
				? (clipx - a->x) / dx
				: 0.0f;
			out[outcount].x = clipx;
			out[outcount].y = a->y + (b->y - a->y) * t;
			outcount++;
			out[outcount++] = *b;
		}
		if (outcount >= COARSE_OCCLUSION_MAX_VERTS)
			return 0;
	}
	return outcount;
}

static s32 R_CoarseOccClipHorizontal(const coarse_occ_vertex_t *in, s32 count,
                          coarse_occ_vertex_t *out, f32 clipy, bool keepGreater)
{ // Sutherland-Hodgman polygon clipping against y >= clipy or y <= clipy.
	s32 i;
	s32 outcount = 0;
	if (count <= 0) return 0;
	for (i = 0; i < count; i++) {
		s32 next = (i + 1 == count) ? 0 : i + 1;
		const coarse_occ_vertex_t *a = &in[i];
		const coarse_occ_vertex_t *b = &in[next];
		bool ina = keepGreater ? (a->y >= clipy) : (a->y <= clipy);
		bool inb = keepGreater ? (b->y >= clipy) : (b->y <= clipy);
		if (ina && inb) {
			out[outcount++] = *b;
		} else if (ina && !inb) {
			f32 dy = b->y - a->y;
			f32 t = (fabsf(dy) > COARSE_OCCLUSION_EPSILON)
				? (clipy - a->y) / dy
				: 0.0f;
			out[outcount].x = a->x + (b->x - a->x) * t;
			out[outcount].y = clipy;
			outcount++;
		} else if (!ina && inb) {
			f32 dy = b->y - a->y;
			f32 t = (fabsf(dy) > COARSE_OCCLUSION_EPSILON)
				? (clipy - a->y) / dy
				: 0.0f;

			out[outcount].x = a->x + (b->x - a->x) * t;
			out[outcount].y = clipy;
			outcount++;

			out[outcount++] = *b;
		}
		if (outcount >= COARSE_OCCLUSION_MAX_VERTS)
			return 0;
	}
	return outcount;
}

// Build projected, viewport-clipped polygon.
// Returns false when we cannot conservatively construct the polygon.
// Important:
// if ANY vertex is behind the near plane, we bail out rather than attempting to
// perform near-plane clipping. This loses optimization opportunities but cannot
// produce a false occlusion result.
static bool R_CoarseOcclusionBuildPolygon(msurface_t *fa,
                              coarse_occ_vertex_t *poly, s32 *outcount)
{
	model_t *model;
	medge_t *edges;
	s32 i;
	s32 count;
	coarse_occ_vertex_t tmp[COARSE_OCCLUSION_MAX_VERTS];
	if (!fa || fa->numedges < 3) return false;
	if (fa->numedges > COARSE_OCCLUSION_MAX_VERTS) return false;
	model = currententity->model;
	edges = model->edges;
	count = fa->numedges;
	for (i = 0; i < count; i++) {
		s32 lindex = model->surfedges[fa->firstedge + i];
		mvertex_t *vert;
		vec3_t local;
		vec3_t transformed;
		f32 zi;
		if (lindex > 0) vert = &r_pcurrentvertbase[edges[lindex].v[0]];
		else vert = &r_pcurrentvertbase[edges[-lindex].v[1]];
		VectorSubtract(vert->position, r_origin, local);
		TransformVector(local, transformed);
		// Do not try to handle near-plane-crossing faces here.
		// R_ClipEdge() remains authoritative for those.
		if (transformed[2] <= NEAR_CLIP) return false;
		zi = 1.0f / transformed[2];
		poly[i].x = xcenter + xscale * zi * transformed[0];
		poly[i].y = ycenter - yscale * zi * transformed[1];
	}
	// Clip polygon to the actual view rectangle.
	// Use the same viewport represented by r_refdef.vrect.
	f32 left = (f32)r_refdef.vrect.x;
	f32 right = (f32)(r_refdef.vrect.x + r_refdef.vrect.width);
	f32 top = (f32)r_refdef.vrect.y;
	f32 bottom = (f32)(r_refdef.vrect.y + r_refdef.vrect.height);
	count = R_CoarseOccClipVertical(poly, count, tmp, left, true);
	if (count < 3) return false;
	count = R_CoarseOccClipVertical(tmp, count, poly, right, false);
	if (count < 3) return false;
	count = R_CoarseOccClipHorizontal(poly, count, tmp, top, true);
	if (count < 3) return false;
	count = R_CoarseOccClipHorizontal(tmp, count, poly, bottom, false);
	if (count < 3) return false;
	*outcount = count;
	return true;
}

// Returns true only when the whole projected visible polygon lies inside tiles
// which are already known to be FULL. False is deliberately the default.
bool R_CoarseOcclusionTestSurface(msurface_t *fa)
{
	coarse_occ_vertex_t poly[COARSE_OCCLUSION_MAX_VERTS];
	s32 count;
	s32 i;
	s32 minx, maxx, miny, maxy;
	s32 tx0, tx1, ty0, ty1;
	bool saw_intersecting_tile = false;
	if (!r_coarse_occlusion_bits) return false;
	// Coarse occlusion is initially restricted to the world.
	// This also avoids entity ordering complications.
	if (currententity != &cl_entities[0]) return false;
	// Never use translucent/sky/cutout geometry as an occluder. The surface
	// test itself could still safely reject against previous opaque
	// coverage, but keeping this restricted makes the initial
	// implementation easier to reason about.
	if (fa->flags & SURF_DRAWSKY) return false;
	if (fa->flags & SURF_DRAWCUTOUT) return false;
	if (fa->flags & SURF_WINQUAKE_DRAWTRANSLUCENT) return false;
	if (!R_CoarseOcclusionBuildPolygon(fa, poly, &count)) return false;
	// Integer pixel-space bounds.
	// Expand slightly so floating-point rounding cannot make us
	// accidentally skip a tile that the polygon actually touches.
	minx = (s32)floorf(poly[0].x) - 1;
	maxx = (s32)ceilf (poly[0].x) + 1;
	miny = (s32)floorf(poly[0].y) - 1;
	maxy = (s32)ceilf (poly[0].y) + 1;
	for (i = 1; i < count; i++) {
		s32 x0 = (s32)floorf(poly[i].x) - 1;
		s32 x1 = (s32)ceilf (poly[i].x) + 1;
		s32 y0 = (s32)floorf(poly[i].y) - 1;
		s32 y1 = (s32)ceilf (poly[i].y) + 1;
		if (x0 < minx) minx = x0;
		if (x1 > maxx) maxx = x1;
		if (y0 < miny) miny = y0;
		if (y1 > maxy) maxy = y1;
	}
	// Convert viewport-space coordinates into tile coordinates.
	tx0 = (minx - r_refdef.vrect.x) / COARSE_OCCLUSION_TILE_SIZE;
	ty0 = (miny - r_refdef.vrect.y) / COARSE_OCCLUSION_TILE_SIZE;
	tx1 = (maxx - r_refdef.vrect.x) / COARSE_OCCLUSION_TILE_SIZE;
	ty1 = (maxy - r_refdef.vrect.y) / COARSE_OCCLUSION_TILE_SIZE;
	if (tx0 < 0) tx0 = 0;
	if (ty0 < 0) ty0 = 0;
	if (tx1 >= r_coarse_occlusion_width) tx1 = r_coarse_occlusion_width - 1;
	if (ty1 >= r_coarse_occlusion_height) ty1 = r_coarse_occlusion_height-1;
	if (tx0 > tx1 || ty0 > ty1) return false;
	for (s32 ty = ty0; ty <= ty1; ty++) {
		for (s32 tx = tx0; tx <= tx1; tx++) {
			f32 x0 = (f32)(r_refdef.vrect.x +
					tx * COARSE_OCCLUSION_TILE_SIZE);
			f32 y0 = (f32)(r_refdef.vrect.y +
					ty * COARSE_OCCLUSION_TILE_SIZE);
			f32 x1 = x0 + COARSE_OCCLUSION_TILE_SIZE;
			f32 y1 = y0 + COARSE_OCCLUSION_TILE_SIZE;
			// Trim edge tiles to the viewport.
			f32 right=(f32)(r_refdef.vrect.x+r_refdef.vrect.width);
			f32 bott =(f32)(r_refdef.vrect.y+r_refdef.vrect.height);
			if (x1 > right)  x1 = right;
			if (y1 > bott) y1 = bott;
			if (!R_CoarseOccPolyIntersectsRect(
				poly, count, x0, y0, x1, y1)) continue;
			saw_intersecting_tile = true;
			// Any intersected tile which is not already completely
			// covered means the face may contribute visible pixels.
			if (!R_CoarseOcclusionTileFull(tx, ty))
				return false;
		}
	}
	// Degenerate / off-screen polygons are not treated as occluded.
	return saw_intersecting_tile;
}


static bool R_CoarseOccluderEligible(msurface_t *fa)
{ // Initial version: only opaque world BSP geometry.
	if (!fa) return 0;
	if (currententity != &cl_entities[0]) return 0;
	if (fa->flags & SURF_DRAWSKY) return 0;
	if (fa->flags & SURF_WINQUAKE_DRAWTRANSLUCENT) return 0;
	if (fa->flags & SURF_DRAWCUTOUT) return 0;
	if (winquake_surface_liquid_alpha < 1.0f) return 0;
	return 1;
}

void R_CoarseOcclusionAddSurface(msurface_t *fa)
{ // Marks only tiles which are COMPLETELY contained by the projected surface
  // polygon. This is intentionally conservative.
	coarse_occ_vertex_t poly[COARSE_OCCLUSION_MAX_VERTS];
	s32 count;
	s32 i;
	s32 minx, maxx, miny, maxy;
	s32 tx0, tx1, ty0, ty1;
	if (!r_coarse_occlusion_bits) return;
	if (!R_CoarseOccluderEligible(fa)) return;
	if (!R_CoarseOcclusionBuildPolygon(fa, poly, &count)) return;
	minx = (s32)floorf(poly[0].x);
	maxx = (s32)ceilf (poly[0].x);
	miny = (s32)floorf(poly[0].y);
	maxy = (s32)ceilf (poly[0].y);
	for (i = 1; i < count; i++) {
		s32 x0 = (s32)floorf(poly[i].x);
		s32 x1 = (s32)ceilf (poly[i].x);
		s32 y0 = (s32)floorf(poly[i].y);
		s32 y1 = (s32)ceilf (poly[i].y);
		if (x0 < minx) minx = x0;
		if (x1 > maxx) maxx = x1;
		if (y0 < miny) miny = y0;
		if (y1 > maxy) maxy = y1;
	}
	tx0 = (minx - r_refdef.vrect.x) / COARSE_OCCLUSION_TILE_SIZE;
	ty0 = (miny - r_refdef.vrect.y) / COARSE_OCCLUSION_TILE_SIZE;
	tx1 = (maxx - r_refdef.vrect.x) / COARSE_OCCLUSION_TILE_SIZE;
	ty1 = (maxy - r_refdef.vrect.y) / COARSE_OCCLUSION_TILE_SIZE;
	if (tx0 < 0) tx0 = 0;
	if (ty0 < 0) ty0 = 0;
	if (tx1 >= r_coarse_occlusion_width)
		tx1 = r_coarse_occlusion_width - 1;
	if (ty1 >= r_coarse_occlusion_height)
		ty1 = r_coarse_occlusion_height - 1;
	if (tx0 > tx1 || ty0 > ty1) return;
	for (s32 ty = ty0; ty <= ty1; ty++) {
		for (s32 tx = tx0; tx <= tx1; tx++) {
			f32 x0 = (f32)(r_refdef.vrect.x +
					tx * COARSE_OCCLUSION_TILE_SIZE);
			f32 y0 = (f32)(r_refdef.vrect.y +
					ty * COARSE_OCCLUSION_TILE_SIZE);
			f32 x1 = x0 + COARSE_OCCLUSION_TILE_SIZE;
			f32 y1 = y0 + COARSE_OCCLUSION_TILE_SIZE;
			f32 right = (f32)(r_refdef.vrect.x +
					r_refdef.vrect.width);
			f32 bottom = (f32)(r_refdef.vrect.y +
					r_refdef.vrect.height);
			if (x1 > right)  x1 = right;
			if (y1 > bottom) y1 = bottom;
			// Only mark a tile if the whole tile is inside the face
			// This is the critical conservative property
			if (R_CoarseOccTileInsidePolygon(
					poly, count, x0, y0, x1, y1))
				R_CoarseOcclusionSetTile(tx, ty);
		}
	}
}

void R_EmitEdge(mvertex_t *pv0, mvertex_t *pv1)
{
	vec3_t local, transformed;
	f32 *world;
	sl64 v, v2, ceilv0;
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
	sl64 side = ceilv0 > r_ceilv1;
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
	// horizontal edge is less than epsilon above a scan, and numeric error
	// causes it to incorrectly extend to the scan, and the extension of the
	// line goes off the edge of the screen
	// FIXME: is this actually needed?
	if (edge->u < r_refdef.vrect_x_adj_shift20)
		edge->u = r_refdef.vrect_x_adj_shift20;
	if (edge->u > r_refdef.vrectright_adj_shift20)
		edge->u = r_refdef.vrectright_adj_shift20;
	sl64 u_check = edge->u; // sort the edge in normally
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
	if ((fa->flags & SURF_DRAWSKY) && skybox_name[0])
		r_skyframe[0] = r_framecount;
	// Manoel Kasimier - skyboxes - end
	if (fa->flags&SURF_WINQUAKE_DRAWTRANSLUCENT && !currententity->alpha) {
		winquake_surface_liquid_alpha=R_LiquidAlphaForFlags(fa->flags);
	} else if (cur_ent_alpha < 1 && r_entalpha.value == 1)
		winquake_surface_liquid_alpha = cur_ent_alpha;
	else winquake_surface_liquid_alpha = 1;
	// Baker: Fully transparent = invisible = don't render
	if (!winquake_surface_liquid_alpha)
		return;
	// Baker: If this is the alpha water pass and we aren't alpha water, get out!
	if (r_alphapass && winquake_surface_liquid_alpha == 1)
		return;
	if (!r_alphapass && winquake_surface_liquid_alpha < 1) {
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
	coccl_enable = (r_coarseocclusion.value == 1)
		|| (r_coarseocclusion.value && 
		cl_entities[0].model->bspversion != BSPVERSION);
	if (!r_showtris.value&&coccl_enable&&R_CoarseOcclusionTestSurface(fa)) {
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
	if (r_showtris.value) { // UNIFIED BSP GEOMETRY HOOK
		mvertex_t *first_vert = NULL;
		mvertex_t *prev_vert = NULL;
		for (i = 0; i < fa->numedges; i++) {
			s32 lindex = currententity->model->surfedges[fa->firstedge + i];
			mvertex_t *vert;
			if (lindex > 0) vert = &r_pcurrentvertbase[pedges[lindex].v[0]];
			else vert = &r_pcurrentvertbase[pedges[-lindex].v[1]];
			if (i == 0) {
				first_vert = vert;
			} else {
				// 1. Draw perimeter edge
				R_DrawDebugLine3D(prev_vert->position, vert->position);
				// 2. Draw inner triangle diagonal if in mode 1
				if (r_showtris.value == 1 && i > 1 && i < fa->numedges - 1) {
					R_DrawDebugLine3D(first_vert->position, vert->position);
				}
			}
			prev_vert = vert;
		}
		if (fa->numedges >= 3) { // 3. Close the perimeter loop
			R_DrawDebugLine3D(prev_vert->position, first_vert->position);
		}
	}
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
	if (R_CoarseOccluderEligible(fa))
	{
		R_CoarseOcclusionAddSurface(fa);
	}
}

void R_RenderBmodelFace(bedge_t *pedges, msurface_t *psurf)
{
	if (psurf->flags&SURF_WINQUAKE_DRAWTRANSLUCENT && !currententity->alpha) {
		winquake_surface_liquid_alpha=R_LiquidAlphaForFlags(psurf->flags);
	} else if (cur_ent_alpha < 1 && r_entalpha.value == 1)
		winquake_surface_liquid_alpha = cur_ent_alpha;
	else winquake_surface_liquid_alpha = 1;
	// Baker: Fully transparent = invisible = don't render
	if (!winquake_surface_liquid_alpha)
		return;
	// Baker: If this is the alpha water pass and we aren't alpha water, get out!
	if (r_alphapass && winquake_surface_liquid_alpha == 1)
		return;
	if (!r_alphapass && winquake_surface_liquid_alpha < 1) {
		r_foundtranswater = 1;
		return;
	}
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
