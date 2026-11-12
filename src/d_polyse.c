// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// d_polyset.c: routines for drawing sets of polygons sharing the same texture
// (used for Alias models)
#include "quakedef.h"

enum { TEXSAMPLE_NEAREST, TEXSAMPLE_DITHERED };
enum { BLEND_OPAQUE, BLEND_MIX, BLEND_STIPPLE };
typedef struct {
	const u8 *colormap;
	s32 blend; // BLEND_*
	s32 mixlevel; // color_mix_lut level, for BLEND_MIX
	s32 level; // D_DitherLevel, for BLEND_STIPPLE
	s32 rgblight; // per-channel coloured lighting instead of the colormap
	const u8 *rgbpix; // D_GetRGBPix of every color, with rgblight
	s32 cutout; // d_cutout
} polysetstyle_t;
typedef struct {
	void *pdest;
	s32 *pz;
	s32 count;
	u8 *ptex;
	s32 sfrac, tfrac, light, zi;
} spanpackage_t;
typedef struct { // the vertices of the edges, indices into r_p
	s32 numleftedges;
	s32 left[3];
	s32 numrightedges;
	s32 right[3];
} edgetable;
// The state of the triangle being drawn, every thread that draws has its own
static Q_TLS s32 r_p[3][6];
static Q_TLS u8 *d_pcolormap;
static Q_TLS s32 d_xdenom;
static Q_TLS const edgetable *pedgetable;
static Q_TLS s32 a_sstepxfrac, a_tstepxfrac, r_lstepx, a_ststepxwhole;
static Q_TLS s32 r_sstepx, r_tstepx, r_lstepy, r_sstepy, r_tstepy;
static Q_TLS s32 r_zistepx, r_zistepy;
static Q_TLS s32 d_aspancount, d_countextrastep;
static Q_TLS spanpackage_t *a_spans;
static Q_TLS spanpackage_t *d_pedgespanpackage;
static Q_TLS u8 *d_pdest, *d_ptex;
static Q_TLS s32 *d_pz;
static Q_TLS s32 d_sfrac, d_tfrac, d_light, d_zi;
static Q_TLS s32 d_ptexextrastep, d_sfracextrastep;
static Q_TLS s32 d_tfracextrastep, d_lightextrastep, d_pdestextrastep;
static Q_TLS s32 d_lightbasestep, d_pdestbasestep, d_ptexbasestep;
static Q_TLS s32 d_sfracbasestep, d_tfracbasestep;
static Q_TLS s32 d_ziextrastep, d_zibasestep;
static Q_TLS s32 d_pzextrastep, d_pzbasestep;
static Q_TLS u8 *skintable[MAX_LBM_HEIGHT];
static Q_TLS s32 ubasestep, errorterm, erroradjustup, erroradjustdown;
// The model being drawn, see D_PolysetBeginModel
static Q_TLS f32 d_entalpha;
static Q_TLS s32 d_rgblight; // coloured light, through d_rgbpix
static Q_TLS const u8 *d_rgbpix;
static Q_TLS s32 d_cutout; // texel to skip: 0xFF if MF_HOLEY, 0x100 (none)

static const edgetable edgetables[12] = {
	{ 1, {0, 2, 0}, 2, {0, 1, 2} },
	{ 2, {1, 0, 2}, 1, {1, 2, 0} },
	{ 1, {0, 2, 0}, 1, {1, 2, 0} },
	{ 1, {1, 0, 0}, 2, {1, 2, 0} },
	{ 2, {0, 2, 1}, 1, {0, 1, 0} },
	{ 1, {2, 1, 0}, 1, {2, 0, 0} },
	{ 1, {2, 1, 0}, 2, {2, 0, 1} },
	{ 2, {2, 1, 0}, 1, {2, 0, 0} },
	{ 1, {1, 0, 0}, 1, {1, 2, 0} },
	{ 1, {2, 1, 0}, 1, {0, 1, 0} },
	{ 1, {1, 0, 0}, 1, {2, 0, 0} },
	{ 1, {0, 2, 0}, 1, {0, 1, 0} },
};

void D_PolysetDrawSpans(spanpackage_t *pspanpackage);
void D_PolysetDrawSpansDithered(spanpackage_t *pspanpackage);
void D_DrawSubdiv();
void D_DrawNonSubdiv();
static void D_DrawNonSubdivBanded();
void D_PolysetRecursiveTriangle(s32 *p1, s32 *p2, s32 *p3);
void D_PolysetSetEdgeTable();
Q_SIL void D_RasterizeAliasPolySmoothImpl(bool banded);
static void D_RasterizeAliasPolySmoothBanded();
void D_PolysetScanLeftEdge(s32 height);

// pix is already lit by the colormap, this only tints it with alias_tint
static inline s32 D_GetRGBPix(s32 pix)
{
	if(pix >= 0xFF) return pix;
	const u8 *c = &CURWORLDPAL[pix * 3];
	s32 rr = (c[0] * alias_tint[0]) >> 8;
	s32 gg = (c[1] * alias_tint[1]) >> 8;
	s32 bb = (c[2] * alias_tint[2]) >> 8;
	if(rr > 255) rr = 255;
	if(gg > 255) gg = 255;
	if(bb > 255) bb = 255;
	return lit_lut[QUANT(rr)][QUANT(gg)][QUANT(bb)];
}

// D_GetRGBPix of all the colors for the current alias_tint, built again when
// that or the palette changes, and once a frame in case the palette was loaded
// again in the same place
static const u8 *D_RGBPixTable()
{
	static u8 tab[256];
	static s32 tint[3];
	static const u8 *pal;
	static s64 frame = -1;
	if(frame != r_framecount || pal != CURWORLDPAL ||
			tint[0] != alias_tint[0] || tint[1] != alias_tint[1] ||
			tint[2] != alias_tint[2]){
		for(s32 i = 0; i < 256; i++)
			tab[i] = D_GetRGBPix(i);
		frame = r_framecount;
		pal = CURWORLDPAL;
		tint[0] = alias_tint[0];
		tint[1] = alias_tint[1];
		tint[2] = alias_tint[2];
	}
	return tab;
}

static void D_PolysetDrawNow()
{
	spanpackage_t spans[(MAXHEIGHT + 1) // +1 for spanpackage marking end
		+ 1 + ((CACHE_SIZE - 1) / sizeof(spanpackage_t)) + 1];
	a_spans = (spanpackage_t *)
		(((uintptr_t) & spans[0] + CACHE_SIZE - 1) & ~(CACHE_SIZE - 1));
	if(r_affinetridesc.drawtype)
		D_DrawSubdiv(); // tests the row of every pixel
	else if(d_bandy0 > 0 || d_bandy1 < (s32)vid.height)
		D_DrawNonSubdivBanded(); // a band of the screen
	else
		D_DrawNonSubdiv();
}

// A pixel of a model drawn by points or by subdivision, its z written
// already. colormap is the row of its light.
Q_SIL void D_PolysetPixel(s32 x, s32 y, const u8 *colormap, s32 texel)
{
	if(texel == d_cutout)
		return;
	s32 pix = colormap[texel];
	if(d_rgblight && texel < 0xE0)
		pix = d_rgbpix[pix];
	u8 *dest = &d_viewbuffer[d_scantable[y] + x];
	if(r_alphastyle.value == 0 && d_entalpha != 1)
		*dest = color_mix_lut[*dest][pix]
			[(s32)(d_entalpha * FOG_LUT_LEVELS)];
	else if(r_alphastyle.value == 1 && d_entalpha != 1){
		if(d_bayer[y & 7][x & 7] < D_DitherLevel(d_entalpha))
			*dest = pix;
	}else
		*dest = pix;
}

static void D_PolysetDrawPoints(finalvert_t *fv, s32 numverts)
{
	for(s32 i = 0; i < numverts; i++, fv++){
		if(fv->v[1] < d_bandy0 || fv->v[1] >= d_bandy1)
			continue;
	// valid triangle coordinates for filling can include the bottom and
	// right clip edges, due to the fill rule; these shouldn't be drawn
		if(fv->v[0]<r_refdef.vrectright&&fv->v[1]<r_refdef.vrectbottom){
			if(fv->v[0] < 0 || fv->v[1] < 0)
				continue; // zbuf segfault
			s32 *zbuf = zspantable[fv->v[1]] + fv->v[0];
			if(zbuf - d_pzbuffer >= (s32)(vid.width * vid.height)
					|| zbuf < d_pzbuffer){
				continue; // zbuf segfault
			}
			s32 z = fv->v[5];
			if(!(z >= *zbuf))
				continue;
			s32 s = fv->v[2] >> 16;
			s32 t = fv->v[3] >> 16;
			if(t < 0 || t >= MAX_LBM_HEIGHT)
				continue;
			u8 *row = skintable[t];
			if(!row)
				continue;
			*zbuf = z;
			D_PolysetPixel(fv->v[0], fv->v[1], (u8 *)acolormap
					+ (fv->v[4] & 0xFF00), row[s]);
		}
	}
}

// The models for the workers (r_thread.c). The command of a model keeps what
// the front end hands D_PolysetDraw and D_PolysetDrawFinalVerts, with copies
// of the vertices and triangles it won't keep, and the bands draw from that.
typedef struct aliasitem_s {
	struct aliasitem_s *next;
	affinetridesc_t desc;
	s32 numverts; // > 0 for the points of D_PolysetDrawFinalVerts
	s32 ymin, ymax;
} aliasitem_t;
typedef struct {
	aliasitem_t *items, **tail;
	void *colormap;
	f32 alpha;
	s32 rgblight;
	u8 rgbpix[256];
	s32 cutout;
	s32 ymin, ymax;
	bool cull; // can skip bands where the z-buffer hides cullrect
	s32 cullrect[5]; // x0, y0, x1, y1 (inclusive) and zbound
} aliascmd_t;
static aliascmd_t *d_alias; // the model being recorded
static finalvert_t *d_modelverts, *d_modelcopy; // its vertices, their copy
static s32 d_modelymin, d_modelymax;

// After R_AliasDrawModel set acolormap, cur_ent_alpha and the coloured light.
// cull: the workers test the rectangle r_entocclusion would, the part of it
// in their bands, against what's in the z-buffer when they get to it.
void D_PolysetBeginModel(bool cull)
{
	d_entalpha = cur_ent_alpha;
	d_rgblight = r_rgblighting.value && colored_aliaslight;
	d_rgbpix = d_rgblight ? D_RGBPixTable() : NULL;
	d_cutout = currententity->model->flags & MF_HOLEY ? 0xFF : 0x100;
	if(r_alphastyle.value == 0 && d_entalpha != 1 && !fog_lut_built)
		R_BuildColorMixLUT(0);
	d_modelverts = NULL;
	d_alias = NULL;
	if(!r_threaded)
		return;
	aliascmd_t *c = d_alias = R_FrameAlloc(sizeof(*c));
	c->items = NULL;
	c->tail = &c->items;
	c->colormap = acolormap;
	c->alpha = d_entalpha;
	c->rgblight = d_rgblight;
	if(d_rgblight)
		memcpy(c->rgbpix, d_rgbpix, sizeof(c->rgbpix));
	c->cutout = d_cutout;
	c->ymin = 0x7fffffff;
	c->ymax = -0x7fffffff;
	s32 *r = c->cullrect;
	c->cull = cull && R_AliasScreenBounds(&r[0], &r[1], &r[2], &r[3],
			&r[4]);
}

// The screen rows of n vertices
static void D_VertRows(const finalvert_t *fv, s32 n, s32 *ymin, s32 *ymax)
{
	*ymin = 0x7fffffff;
	*ymax = -0x7fffffff;
	for(s32 i = 0; i < n; i++){
		*ymin = q_min(*ymin, fv[i].v[1]);
		*ymax = q_max(*ymax, fv[i].v[1]);
	}
}

// The transformed vertices the triangles of the model index. projected: all
// of them, the model is drawn whole.
void D_PolysetModelVerts(finalvert_t *fv, s32 numverts, bool projected)
{
	d_modelverts = fv;
	if(!d_alias)
		return;
	d_modelcopy = R_FrameAlloc(numverts * sizeof(*fv));
	memcpy(d_modelcopy, fv, numverts * sizeof(*fv));
	D_VertRows(fv, projected ? numverts : 0, &d_modelymin, &d_modelymax);
}

static void D_AliasAdd(aliasitem_t *it)
{
	it->next = NULL;
	*d_alias->tail = it;
	d_alias->tail = &it->next;
	d_alias->ymin = q_min(d_alias->ymin, it->ymin);
	d_alias->ymax = q_max(d_alias->ymax, it->ymax);
}

static finalvert_t *D_AliasVerts(finalvert_t *fv, s32 n)
{
	if(fv == d_modelverts)
		return d_modelcopy;
	finalvert_t *copy = R_FrameAlloc(n * sizeof(*fv));
	memcpy(copy, fv, n * sizeof(*fv));
	return copy;
}

static void D_AliasRecordTris()
{
	aliasitem_t *it = R_FrameAlloc(sizeof(*it));
	it->desc = r_affinetridesc;
	it->numverts = 0;
	finalvert_t *fv = r_affinetridesc.pfinalverts;
	mtriangle_t *tri = r_affinetridesc.ptriangles;
	if(r_affinetridesc.numtriangles == 1){ // one of the model's, or
		mtriangle_t *t = R_FrameAlloc(sizeof(*t)); // clipped
		*t = *tri;
		it->desc.ptriangles = t;
		s32 n = 0;
		it->ymin = 0x7fffffff;
		it->ymax = -0x7fffffff;
		for(s32 i = 0; i < 3; i++){
			const finalvert_t *v = &fv[tri->vertindex[i]];
			n = q_max(n, tri->vertindex[i] + 1);
			it->ymin = q_min(it->ymin, v->v[1]);
			it->ymax = q_max(it->ymax, v->v[1]);
		}
		it->desc.pfinalverts = D_AliasVerts(fv, n);
	}else{ // all of the model's
		it->desc.pfinalverts = d_modelcopy;
		it->ymin = d_modelymin;
		it->ymax = d_modelymax;
	}
	D_AliasAdd(it);
}

static void D_AliasCmd(void *arg)
{
	const aliascmd_t *c = arg;
	if(c->cull){
		const s32 *r = c->cullrect;
		s32 y0 = q_max(r[1], d_bandy0), y1 = q_min(r[3], d_bandy1 - 1);
		if(y0 > y1 || R_ZBufferRectBehind(r[0], y0, r[2], y1, r[4]))
			return;
	}
	acolormap = c->colormap;
	d_entalpha = c->alpha;
	d_rgblight = c->rgblight;
	d_rgbpix = c->rgbpix;
	d_cutout = c->cutout;
	for(const aliasitem_t *it = c->items; it; it = it->next){
		if(it->ymax < d_bandy0 || it->ymin >= d_bandy1)
			continue;
		r_affinetridesc = it->desc;
		if(r_affinetridesc.drawtype)
			D_PolysetUpdateTables();
		if(it->numverts)
			D_PolysetDrawPoints(r_affinetridesc.pfinalverts,
					it->numverts);
		else
			D_PolysetDrawNow();
	}
}

void D_PolysetEndModel()
{
	aliascmd_t *c = d_alias;
	d_alias = NULL;
	if(c && c->items)
		R_RowsAdd(c->ymin, c->ymax + 1, D_AliasCmd, c);
}

void D_PolysetDraw()
{
	if(d_alias)
		D_AliasRecordTris();
	else
		D_PolysetDrawNow();
}

void D_PolysetDrawFinalVerts(finalvert_t *fv, s32 numverts)
{
	if(!d_alias){
		D_PolysetDrawPoints(fv, numverts);
		return;
	}
	aliasitem_t *it = R_FrameAlloc(sizeof(*it));
	it->desc = r_affinetridesc;
	it->desc.pfinalverts = D_AliasVerts(fv, numverts);
	it->numverts = numverts;
	D_VertRows(fv, numverts, &it->ymin, &it->ymax);
	D_AliasAdd(it);
}

// Can the triangle have pixels in the rows this thread draws? full: all of it
Q_SIL bool D_TriInBand(const finalvert_t *a, const finalvert_t *b,
		const finalvert_t *c, bool *full)
{
	s32 ymin = q_min(a->v[1], q_min(b->v[1], c->v[1]));
	s32 ymax = q_max(a->v[1], q_max(b->v[1], c->v[1]));
	*full = ymin >= d_bandy0 && ymax <= d_bandy1;
	return ymax >= d_bandy0 && ymin < d_bandy1;
}

void D_DrawSubdiv()
{
	finalvert_t *pfv = r_affinetridesc.pfinalverts;
	mtriangle_t *ptri = r_affinetridesc.ptriangles;
	s32 lnumtriangles = r_affinetridesc.numtriangles;
	for(s32 i = lnumtriangles - 1; i >= 0; i--){
		finalvert_t *i0 = pfv + ptri[i].vertindex[0];
		finalvert_t *i1 = pfv + ptri[i].vertindex[1];
		finalvert_t *i2 = pfv + ptri[i].vertindex[2];
		if(r_showtris.value && r_numdebuglines < MAX_DEBUG_LINES - 3){
			R_DrawDebugLine(i0->v[0],i0->v[1],i1->v[0],i1->v[1],
					r_showtris_color.value);
			R_DrawDebugLine(i1->v[0],i1->v[1],i2->v[0],i2->v[1],
					r_showtris_color.value);
			R_DrawDebugLine(i2->v[0],i2->v[1],i0->v[0],i0->v[1],
					r_showtris_color.value);
		}
		if(((i0->v[1]-i1->v[1])*(i0->v[0]-i2->v[0])-(i0->v[0]-i1->v[0])
					*(i0->v[1]-i2->v[1]))>=0)
			continue;
		d_pcolormap = &((u8 *) acolormap)[i0->v[4] & 0xFF00];
		if(ptri[i].facesfront){
			D_PolysetRecursiveTriangle(i0->v, i1->v, i2->v);
		}else{ // copies, other threads draw from the same vertices
			s32 v0[6], v1[6], v2[6];
			memcpy(v0, i0->v, sizeof(v0));
			memcpy(v1, i1->v, sizeof(v1));
			memcpy(v2, i2->v, sizeof(v2));
			if(i0->flags & ALIAS_ONSEAM)
				v0[2] += r_affinetridesc.seamfixupX16;
			if(i1->flags & ALIAS_ONSEAM)
				v1[2] += r_affinetridesc.seamfixupX16;
			if(i2->flags & ALIAS_ONSEAM)
				v2[2] += r_affinetridesc.seamfixupX16;
			D_PolysetRecursiveTriangle(v0, v1, v2);
		}
	}
}

Q_SIL void D_DrawNonSubdivImpl(bool banded)
{
	finalvert_t *pfv = r_affinetridesc.pfinalverts;
	mtriangle_t *ptri = r_affinetridesc.ptriangles;
	s32 lnumtriangles = r_affinetridesc.numtriangles;
	ptri += lnumtriangles - 1;
	for(s32 i = 0; i < lnumtriangles; i++, ptri--){
		finalvert_t *i0 = pfv + ptri->vertindex[0];
		finalvert_t *i1 = pfv + ptri->vertindex[1];
		finalvert_t *i2 = pfv + ptri->vertindex[2];
		bool full = 1;
		if(banded && !D_TriInBand(i0, i1, i2, &full))
			continue;
		if(r_showtris.value && r_numdebuglines < MAX_DEBUG_LINES - 3){
			R_DrawDebugLine(i0->v[0],i0->v[1],i1->v[0],i1->v[1],
					r_showtris_color.value);
			R_DrawDebugLine(i1->v[0],i1->v[1],i2->v[0],i2->v[1],
					r_showtris_color.value);
			R_DrawDebugLine(i2->v[0],i2->v[1],i0->v[0],i0->v[1],
					r_showtris_color.value);
		}
		d_xdenom=(i0->v[1]-i1->v[1])*(i0->v[0]-i2->v[0])-
			(i0->v[0]-i1->v[0])*(i0->v[1]-i2->v[1]);
		if(d_xdenom >= 0)
			continue;
		memcpy(r_p[0], i0->v, sizeof(r_p[0])); // u, v, s, t, light, iz
		memcpy(r_p[1], i1->v, sizeof(r_p[1]));
		memcpy(r_p[2], i2->v, sizeof(r_p[2]));
		if(!ptri->facesfront){
			if(i0->flags & ALIAS_ONSEAM)
				r_p[0][2] += r_affinetridesc.seamfixupX16;
			if(i1->flags & ALIAS_ONSEAM)
				r_p[1][2] += r_affinetridesc.seamfixupX16;
			if(i2->flags & ALIAS_ONSEAM)
				r_p[2][2] += r_affinetridesc.seamfixupX16;
		}
		D_PolysetSetEdgeTable();
		if(full) // inlined, as without the workers there's only this
			D_RasterizeAliasPolySmoothImpl(0);
		else
			D_RasterizeAliasPolySmoothBanded();
	}
}

void D_DrawNonSubdiv(){ D_DrawNonSubdivImpl(0); }
static void D_DrawNonSubdivBanded(){ D_DrawNonSubdivImpl(1); }

void D_PolysetRecursiveTriangle(s32 *lp1, s32 *lp2, s32 *lp3)
{
	s32 *temp;
	s32 new[6]; // keep here for OpenBSD compiler
	s32 d = lp2[0] - lp1[0];
	if(d < -1 || d > 1)
		goto split;
	d = lp2[1] - lp1[1];
	if(d < -1 || d > 1)
		goto split;
	d = lp3[0] - lp2[0];
	if(d < -1 || d > 1)
		goto split2;
	d = lp3[1] - lp2[1];
	if(d < -1 || d > 1)
		goto split2;
	d = lp1[0] - lp3[0];
	if(d < -1 || d > 1)
		goto split3;
	d = lp1[1] - lp3[1];
	if(d < -1 || d > 1){
split3:
		temp = lp1;
		lp1 = lp3;
		lp3 = lp2;
		lp2 = temp;
		goto split;
	}
	return; // entire tri is filled
split2:
	temp = lp1;
	lp1 = lp2;
	lp2 = lp3;
	lp3 = temp;
split: // split this edge
	new[0] = (lp1[0] + lp2[0]) >> 1;
	new[1] = (lp1[1] + lp2[1]) >> 1;
	new[2] = (lp1[2] + lp2[2]) >> 1;
	new[3] = (lp1[3] + lp2[3]) >> 1;
	new[5] = (lp1[5] + lp2[5]) >> 1;
	if(lp2[1] > lp1[1]) // draw the point if splitting a leading edge
		goto nodraw;
	if((lp2[1] == lp1[1]) && (lp2[0] < lp1[0]))
		goto nodraw;
	if(new[0] < 0 || new[1] < 0 || new[1] >= r_refdef.fvrectbottom
			|| new[0] >= r_refdef.fvrectright)
		goto nodraw; // zbuf segfault
	if(new[1] < d_bandy0 || new[1] >= d_bandy1)
		goto nodraw; // another thread's
	s32 *zbuf = zspantable[new[1]] + new[0];
	if(zbuf - d_pzbuffer >= (s32)(vid.width * vid.height)
			|| zbuf < d_pzbuffer)
		goto nodraw; // zbuf segfault
	s32 z = new[5];
	if(z >= *zbuf){
		*zbuf = z;
		s32 s = new[2] >> 16;
		s32 t = new[3] >> 16;
		if(t < 0 || t >= MAX_LBM_HEIGHT)
			goto nodraw;
		u8 *row = skintable[t];
		if(!row)
			goto nodraw;
		D_PolysetPixel(new[0], new[1], d_pcolormap, row[s]);
	}
nodraw: // recursively continue
	D_PolysetRecursiveTriangle(lp3, lp1, new);
	D_PolysetRecursiveTriangle(lp3, new, lp2);
}

void D_PolysetUpdateTables()
{
	static Q_TLS u8 *skinstart; // what this thread's skintable holds
	static Q_TLS s32 skinwidth;
	if(r_affinetridesc.skinwidth != skinwidth ||
			r_affinetridesc.pskin != skinstart){
		skinwidth = r_affinetridesc.skinwidth;
		skinstart = r_affinetridesc.pskin;
		u8 *s = skinstart; // rows past the skin stay NULL
		s32 rows = r_affinetridesc.skinheight;
		for(s32 i = 0; i < MAX_LBM_HEIGHT; i++, s += skinwidth)
			skintable[i] = i < rows ? s : NULL;
	}
}

Q_SIL void D_PolysetStoreSpanPackage()
{
	d_pedgespanpackage->pdest = d_pdest;
	d_pedgespanpackage->pz = d_pz;
	d_pedgespanpackage->count = d_aspancount;
	d_pedgespanpackage->ptex = d_ptex;
	d_pedgespanpackage->sfrac = d_sfrac;
	d_pedgespanpackage->tfrac = d_tfrac;
	d_pedgespanpackage->light = d_light;
	d_pedgespanpackage->zi = d_zi;
	d_pedgespanpackage++;
}

void D_PolysetScanLeftEdge(s32 height)
{
	do{
		D_PolysetStoreSpanPackage();
		errorterm += erroradjustup;
		if(errorterm >= 0){
			d_pdest += d_pdestextrastep;
			d_pz += d_pzextrastep;
			d_aspancount += d_countextrastep;
			d_ptex += d_ptexextrastep;
			d_sfrac += d_sfracextrastep;
			d_ptex += d_sfrac >> 16;
			d_sfrac &= 0xFFFF;
			d_tfrac += d_tfracextrastep;
			if(d_tfrac & 0x10000){
				d_ptex += r_affinetridesc.skinwidth;
				d_tfrac &= 0xFFFF;
			}
			d_light += d_lightextrastep;
			d_zi += d_ziextrastep;
			errorterm -= erroradjustdown;
		}else{
			d_pdest += d_pdestbasestep;
			d_pz += d_pzbasestep;
			d_aspancount += ubasestep;
			d_ptex += d_ptexbasestep;
			d_sfrac += d_sfracbasestep;
			d_ptex += d_sfrac >> 16;
			d_sfrac &= 0xFFFF;
			d_tfrac += d_tfracbasestep;
			if(d_tfrac & 0x10000){
				d_ptex += r_affinetridesc.skinwidth;
				d_tfrac &= 0xFFFF;
			}
			d_light += d_lightbasestep;
			d_zi += d_zibasestep;
		}
	}while(--height);
}

void D_PolysetSetUpForLineScan(s32 startvertu, s32 startvertv,
		s32 endvertu, s32 endvertv)
{
	errorterm = -1;
	s32 tm = endvertu - startvertu;
	s32 tn = endvertv - startvertv;
	FloorDivMod(tm, tn, &ubasestep, &erroradjustup);
	erroradjustdown = tn;
}

Q_SIL void D_PolysetCalcGradients(s32 skinwidth)
{
	f32 p00_minus_p20 = r_p[0][0] - r_p[2][0];
	f32 p01_minus_p21 = r_p[0][1] - r_p[2][1];
	f32 p10_minus_p20 = r_p[1][0] - r_p[2][0];
	f32 p11_minus_p21 = r_p[1][1] - r_p[2][1];
	f32 xstepdenominv = 1.0 / (f32)d_xdenom;
	f32 ystepdenominv = -xstepdenominv;
// ceil() for light so positive steps are exaggerated, negative steps
// diminished, pushing us away from underflow toward overflow. Underflow is
// very visible, overflow is very unlikely, because of ambient lighting
	f32 t0 = r_p[0][4] - r_p[2][4];
	f32 t1 = r_p[1][4] - r_p[2][4];
	r_lstepx = Q_ceil((t1*p01_minus_p21-t0*p11_minus_p21)*xstepdenominv);
	r_lstepy = Q_ceil((t1*p00_minus_p20-t0*p10_minus_p20)*ystepdenominv);
	t0 = r_p[0][2] - r_p[2][2];
	t1 = r_p[1][2] - r_p[2][2];
	r_sstepx = (s32)((t1*p01_minus_p21-t0*p11_minus_p21)*xstepdenominv);
	r_sstepy = (s32)((t1*p00_minus_p20-t0*p10_minus_p20)*ystepdenominv);
	t0 = r_p[0][3] - r_p[2][3];
	t1 = r_p[1][3] - r_p[2][3];
	r_tstepx = (s32)((t1*p01_minus_p21-t0*p11_minus_p21)*xstepdenominv);
	r_tstepy = (s32)((t1*p00_minus_p20-t0*p10_minus_p20)*ystepdenominv);
	t0 = r_p[0][5] - r_p[2][5];
	t1 = r_p[1][5] - r_p[2][5];
	r_zistepx = (s32)((t1*p01_minus_p21-t0*p11_minus_p21)*xstepdenominv);
	r_zistepy = (s32)((t1*p00_minus_p20-t0*p10_minus_p20)*ystepdenominv);
	a_sstepxfrac = r_sstepx & 0xFFFF;
	a_tstepxfrac = r_tstepx & 0xFFFF;
	a_ststepxwhole = skinwidth * (r_tstepx >> 16) + (r_sstepx >> 16);
}

static Q_NOINLINE polysetstyle_t D_PolysetPixelStyle()
{
	polysetstyle_t st;
	st.colormap = (u8 *)acolormap;
	st.rgblight = d_rgblight;
	st.rgbpix = d_rgbpix;
	st.cutout = d_cutout;
	st.mixlevel = 0;
	st.level = 0;
	if(d_entalpha == 1)
		st.blend = BLEND_OPAQUE;
	else if(r_alphastyle.value == 0){
		st.blend = BLEND_MIX;
		st.mixlevel = (s32)(d_entalpha * FOG_LUT_LEVELS);
	}else if(r_alphastyle.value == 1){
		st.blend = BLEND_STIPPLE;
		st.level = D_DitherLevel(d_entalpha);
	}
	else
		st.blend = BLEND_OPAQUE;
	return st;
}
 

// The skin column and row of ptex (its offset floor-divided by the skin width)
// stepped along with it, and the skin in locals: the pixel stores could
// change anything behind a pointer
typedef struct {
	s32 col, row;
	s32 colstep, rowstep; // whole texels a pixel, 0 <= colstep < skinwidth
	s32 ds, dt; // the dither offsets of this pixel
	s32 dsx, dtx; // xor that swaps them for the next one
	const u8 *skin, *mask;
	s32 skinw, tmax;
} skinwalk_t;

static void D_SkinWalkSetup(skinwalk_t *k)
{ // once a triangle
	s32 skinw = r_affinetridesc.skinwidth;
	s32 sw = r_sstepx >> 16;
	k->colstep = sw % skinw;
	k->rowstep = (r_tstepx >> 16) + sw / skinw;
	if(k->colstep < 0){ k->colstep += skinw; k->rowstep--; }
	k->skin = r_affinetridesc.pskin;
	k->mask = r_affinetridesc.pskinmask;
	k->skinw = skinw;
	k->tmax = r_affinetridesc.skinheight - 1;
}

Q_SIL void D_SkinWalkBegin(skinwalk_t *k, const u8 *ptex, const u8 *pdest)
{
	s32 off = ptex - k->skin;
	k->row = off / k->skinw;
	k->col = off % k->skinw;
	if(k->col < 0){ k->col += k->skinw; k->row--; }
	s32 pixel = pdest - d_viewbuffer;
	s32 y = pixel / screenwidth;
	s32 i = ((pixel - y * screenwidth) & 1) + ((y & 1) << 1);
	k->ds = dither_s[i];
	k->dt = dither_t[i];
	k->dsx = dither_s[i] ^ dither_s[i ^ 1];
	k->dtx = dither_t[i] ^ dither_t[i ^ 1];
}

Q_SIL void D_SkinWalkStep(skinwalk_t *k, s32 scarry, s32 tcarry)
{ // the same carries that step ptex, at most one wrap into the next row
	k->col += k->colstep + scarry;
	k->row += k->rowstep + tcarry;
	if(k->col >= k->skinw){ k->col -= k->skinw; k->row++; }
	k->ds ^= k->dsx;
	k->dt ^= k->dtx;
}

Q_SIL u8 D_PolysetDitheredTexel(const skinwalk_t *k, s32 sfrac, s32 tfrac)
{ // CyanBun96: dithered sampling from Unreal
	s32 col = k->row < 0 ? 0 : k->col; // before the skin reads as its start
	s32 row = k->row < 0 ? 0 : k->row;
	s32 s = ((col << 16) | (sfrac & 0xFFFF)) + k->ds;
	s32 t = ((row << 16) | (tfrac & 0xFFFF)) + k->dt;
	// clamped to the skin, the offsets are <= 0 and col < skinw, so s can
	// only leave it on the left
	s = s < 0 ? 0 : s >> 16;
	t = t < 0 ? 0 : q_min(t >> 16, k->tmax);
	s32 i = s + t * k->skinw;
	// the offset can step off the triangles, onto whatever is around them
	return k->mask[i >> 3] & (1 << (i & 7)) ? k->skin[i]
		: k->skin[k->row * k->skinw + k->col]; // *ptex
}
 
// bay and row are d_bayer's row and the screen row of pdest, for stipple
Q_SIL void D_PolysetPlot(const polysetstyle_t *st, u8 texel, s32 light,
		u8 *pdest, s32 *pz, s32 zi, const u8 *bay, const u8 *row)
{
	if(texel == st->cutout)
		return;
	s32 pix = st->colormap[texel + (light & 0xFF00)];
	if(st->rgblight && texel < 0xE0)
		pix = st->rgbpix[pix];
	switch(st->blend){
	case BLEND_MIX:
		*pdest = color_mix_lut[*pdest][pix][st->mixlevel];
		break;
	case BLEND_STIPPLE:
		if(bay[(pdest - row) & 7] < st->level)
			*pdest = pix;
		break;
	default:
		*pdest = pix;
	}
	*pz = zi;
}
 
Q_SIL void D_PolysetClampSpanLight(s32 *light, s32 *lstep, s32 count)
{ // light is linear across the span, so keeping both ends inside the colormap
  // keeps every pixel inside it. The rounded-up gradients can overshoot by
  // ~1 per step, which adds up on long spans at high resolutions
	if(*light < 0) *light = 0;
	else if(*light > LIGHT_MAX) *light = LIGHT_MAX;
	if(count < 2)
		return;
	s32 end = *light + *lstep * (count - 1);
	if(end > LIGHT_MAX)
		*lstep = (LIGHT_MAX - *light) / (count - 1);
	else if(end < 0)
		*lstep = -*light / (count - 1);
}

Q_SIL void D_PolysetDrawSpan(const spanpackage_t *span, s32 count,
		const polysetstyle_t *st, s32 sample, skinwalk_t *k)
{
	u8 *pdest = span->pdest;
	s32 *pz = span->pz;
	u8 *ptex = span->ptex;
	s32 sfrac = span->sfrac;
	s32 light = span->light;
	s32 lstep = r_lstepx;
	s32 zi = span->zi;
	D_PolysetClampSpanLight(&light, &lstep, count);
	// tfrac:light share a register like in d_polysa.s, so light carries
	// into tfrac there, and tfrac's carry is the step to the next skin row
	u32 tl = ((u32)span->tfrac << 16) | (light & 0xFFFF);
	u32 tlstep = ((u32)a_tstepxfrac << 16) | (lstep & 0xFFFF);
	s32 skinwidth = r_affinetridesc.skinwidth;
	const u8 *bay = NULL, *row = NULL;
	if(st->blend == BLEND_STIPPLE){
		s32 y = (pdest - d_viewbuffer) / screenwidth;
		bay = d_bayer[y & 7];
		row = d_viewbuffer + y * screenwidth;
	}
	if(sample == TEXSAMPLE_DITHERED)
		D_SkinWalkBegin(k, ptex, pdest);
	do{
		if(zi >= *pz){
			u8 texel = sample == TEXSAMPLE_DITHERED ?
			  D_PolysetDitheredTexel(k, sfrac, tl >> 16) : *ptex;
			D_PolysetPlot(st, texel, tl & 0xFFFF, pdest, pz, zi,
					bay, row);
		}
		pdest++;
		pz++;
		zi += r_zistepx;
		ptex += a_ststepxwhole;
		sfrac += a_sstepxfrac;
		s32 scarry = sfrac >> 16;
		ptex += scarry;
		sfrac &= 0xFFFF;
		u32 tlold = tl;
		tl += tlstep;
		s32 tcarry = tl < tlold;
		ptex += skinwidth & -tcarry;
		if(sample == TEXSAMPLE_DITHERED)
			D_SkinWalkStep(k, scarry, tcarry);
	}while(--count);
}
 
Q_SIL s32 D_PolysetNextSpanWidth(const spanpackage_t *span)
{
	s32 width = d_aspancount - span->count;
	errorterm += erroradjustup;
	if(errorterm >= 0){
		d_aspancount += d_countextrastep;
		errorterm -= erroradjustdown;
	}else
		d_aspancount += ubasestep;
	return width;
}
 
Q_SIL s32 D_PolysetSpanInBounds(const spanpackage_t *span, s32 count)
{
	const u8 *dest = span->pdest;
	const s32 *pz = span->pz;
	const u8 *dest_end = d_viewbuffer + vid.width * vid.height;
	const s32 *z_end = d_pzbuffer + vid.width * vid.height;
	return dest >= d_viewbuffer && dest + count <= dest_end &&
		pz >= d_pzbuffer && pz + count <= z_end;
}
 
// The span packages from screen row y on. banded: only the rows of this
// thread's band, for triangles that don't fit into it
Q_SIL void D_PolysetDrawSpansImpl(spanpackage_t *pspanpackage, s32 y,
		s32 sample, bool banded)
{
	const polysetstyle_t st = D_PolysetPixelStyle();
	skinwalk_t k;
	if(sample == TEXSAMPLE_DITHERED)
		D_SkinWalkSetup(&k);
	do{
		if(banded && y >= d_bandy1)
			break;
		s32 count = D_PolysetNextSpanWidth(pspanpackage);
		if(count > 0 && (!banded || y >= d_bandy0)
				&& D_PolysetSpanInBounds(pspanpackage, count))
			D_PolysetDrawSpan(pspanpackage, count, &st, sample,
					&k);
		pspanpackage++;
		y++;
	}while(pspanpackage->count != -999999);
}

void D_PolysetDrawSpans(spanpackage_t *pspanpackage)
{ D_PolysetDrawSpansImpl(pspanpackage, 0, TEXSAMPLE_NEAREST, 0); }

void D_PolysetDrawSpansDithered(spanpackage_t *pspanpackage)
{ D_PolysetDrawSpansImpl(pspanpackage, 0, TEXSAMPLE_DITHERED, 0); }

static Q_NOINLINE void D_PolysetDrawSpansBanded(spanpackage_t *pspanpackage,
		s32 y)
{ D_PolysetDrawSpansImpl(pspanpackage, y, TEXSAMPLE_NEAREST, 1); }

static Q_NOINLINE void D_PolysetDrawSpansDitheredBanded(
		spanpackage_t *pspanpackage, s32 y)
{ D_PolysetDrawSpansImpl(pspanpackage, y, TEXSAMPLE_DITHERED, 1); }

Q_SIL void D_PolysetSpans(spanpackage_t *pspanpackage, s32 y, bool banded)
{
	if(banded && r_dithertex.value)
		D_PolysetDrawSpansDitheredBanded(pspanpackage, y);
	else if(banded) D_PolysetDrawSpansBanded(pspanpackage, y);
	else if(r_dithertex.value) D_PolysetDrawSpansDithered(pspanpackage);
	else D_PolysetDrawSpans(pspanpackage);
}

// sets up the stepping for a left edge of at least two scanlines
Q_SIL void D_PolysetSetUpLeftEdgeSteps()
{
	d_pdestbasestep = screenwidth + ubasestep;
	d_pdestextrastep = d_pdestbasestep + 1;
	d_pzbasestep = d_zwidth + ubasestep;
	d_pzextrastep = d_pzbasestep + 1;
	d_countextrastep = ubasestep + 1;
	// for negative steps in x along left edge, bias toward overflow
	// rather than underflow(sort of turning the floor() we did in
	// the gradient calcs into ceil(), but plus a little bit)
	s32 working_lstepx = ubasestep < 0 ? r_lstepx - 1 : r_lstepx;
	// the extra step is the base step plus one more x step
	s32 sbase = r_sstepy + r_sstepx * ubasestep;
	s32 tbase = r_tstepy + r_tstepx * ubasestep;
	s32 sextra = sbase + r_sstepx;
	s32 textra = tbase + r_tstepx;
	s32 skinwidth = r_affinetridesc.skinwidth;
	d_ptexbasestep = (sbase >> 16) + (tbase >> 16) * skinwidth;
	d_sfracbasestep = sbase & 0xFFFF;
	d_tfracbasestep = tbase & 0xFFFF;
	d_lightbasestep = r_lstepy + working_lstepx * ubasestep;
	d_zibasestep = r_zistepy + r_zistepx * ubasestep;
	d_ptexextrastep = (sextra >> 16) + (textra >> 16) * skinwidth;
	d_sfracextrastep = sextra & 0xFFFF;
	d_tfracextrastep = textra & 0xFFFF;
	d_lightextrastep = d_lightbasestep + working_lstepx;
	d_ziextrastep = d_zibasestep + r_zistepx;
}

// scans out one left edge into span packages; fracmask is 0xFFFF to start
// from the top vertex's s/t fractions, or 0 to drop them like the original
// code does for the bottom edge
Q_SIL void D_PolysetLeftEdge(const s32 *ptop, const s32 *pbottom,
		s32 rightx, s32 height, s32 fracmask)
{
	s32 ystart = ptop[1];
	d_aspancount = ptop[0] - rightx;
	d_ptex = (u8 *) r_affinetridesc.pskin + (ptop[2] >> 16) +
		(ptop[3] >> 16) * r_affinetridesc.skinwidth;
	d_sfrac = ptop[2] & fracmask;
	d_tfrac = ptop[3] & fracmask;
	d_light = ptop[4];
	d_zi = ptop[5];
	d_pdest = (u8 *) d_viewbuffer + ystart * screenwidth + ptop[0];
	d_pz = d_pzbuffer + ystart * d_zwidth + ptop[0];
	if(height == 1){
		D_PolysetStoreSpanPackage();
		return;
	}
	D_PolysetSetUpForLineScan(ptop[0], ptop[1], pbottom[0], pbottom[1]);
	D_PolysetSetUpLeftEdgeSteps();
	D_PolysetScanLeftEdge(height);
}

// banded: the rows of the triangle in this thread's band, it has others too
Q_SIL void D_RasterizeAliasPolySmoothImpl(bool banded)
{
	s32 *plefttop = r_p[pedgetable->left[0]];
	s32 *prighttop = r_p[pedgetable->right[0]];
	s32 *pleftbottom = r_p[pedgetable->left[1]];
	s32 *prightbottom = r_p[pedgetable->right[1]];
	u32 initialleftheight = pleftbottom[1] - plefttop[1];
	u32 initialrightheight = prightbottom[1] - prighttop[1];
	if(initialleftheight>vid.height || initialrightheight>vid.height)
		return; // broken alias model
	s32 y = plefttop[1]; // the row of the first span
	// set the s, t, and light gradients, which are consistent across the
	// triangle because being a triangle, things are affine
	D_PolysetCalcGradients(r_affinetridesc.skinwidth);
	// rasterize the polygon
	// scan out the top(and possibly only) part of the left edge
	d_pedgespanpackage = a_spans;
	D_PolysetLeftEdge(plefttop, pleftbottom, prighttop[0],
			initialleftheight, 0xFFFF);
	// scan out the bottom part of the left edge, if it exists
	if(pedgetable->numleftedges == 2){
		plefttop = pleftbottom;
		pleftbottom = r_p[pedgetable->left[2]];
		D_PolysetLeftEdge(plefttop, pleftbottom, prighttop[0],
				pleftbottom[1] - plefttop[1], 0);
	}
	// scan out the top(and possibly only) part of the right edge, updating
	// the count field
	d_pedgespanpackage = a_spans;
	D_PolysetSetUpForLineScan(prighttop[0], prighttop[1],
			prightbottom[0], prightbottom[1]);
	d_aspancount = 0;
	d_countextrastep = ubasestep + 1;
	s32 originalcount = a_spans[initialrightheight].count;
	a_spans[initialrightheight].count = -999999; // mark end of spanpackages
	D_PolysetSpans(a_spans, y, banded);
	// scan out the bottom part of the right edge, if it exists
	if(pedgetable->numrightedges == 2){
		spanpackage_t *pstart;
		pstart = a_spans + initialrightheight;
		pstart->count = originalcount;
		d_aspancount = prightbottom[0] - prighttop[0];
		prighttop = prightbottom;
		prightbottom = r_p[pedgetable->right[2]];
		s32 height = prightbottom[1] - prighttop[1];
		D_PolysetSetUpForLineScan(prighttop[0], prighttop[1],
				prightbottom[0], prightbottom[1]);
		d_countextrastep = ubasestep + 1;
		//mark end of the spanpackages
		a_spans[initialrightheight + height].count = -999999;
		D_PolysetSpans(pstart, y + initialrightheight, banded);
	}
}

static void D_RasterizeAliasPolySmoothBanded()
{ D_RasterizeAliasPolySmoothImpl(1); }

void D_PolysetSetEdgeTable()
{
	s32 edgetableindex = 0; // assume the vertices are already in top to
		// bottom order determine which edges are right & left, and the
		// order in which to rasterize them
	if(r_p[0][1] >= r_p[1][1]){
		if(r_p[0][1] == r_p[1][1]){
			pedgetable = &edgetables[r_p[0][1]<r_p[2][1] ? 2 : 5];
			return;
		}else
			edgetableindex = 1;
	}
	if(r_p[0][1] == r_p[2][1]){
		pedgetable = edgetableindex ? &edgetables[8] : &edgetables[9];
		return;
	}else if(r_p[1][1] == r_p[2][1]){
		pedgetable = edgetableindex ? &edgetables[10] : &edgetables[11];
		return;
	}
	if(r_p[0][1] > r_p[2][1])
		edgetableindex += 2;
	if(r_p[1][1] > r_p[2][1])
		edgetableindex += 4;
	pedgetable = &edgetables[edgetableindex];
}
