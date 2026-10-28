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
	s32 rgblight; // per-channel coloured lighting instead of the colormap
	const u8 *rgbpix; // D_GetRGBPix of every color, with rgblight
} polysetstyle_t;
typedef struct {
	void *pdest;
	s16 *pz;
	s32 count;
	u8 *ptex;
	s32 sfrac, tfrac, light, zi;
} spanpackage_t;
typedef struct {
	s32 isflattop;
	s32 numleftedges;
	s32 *pleftedgevert0;
	s32 *pleftedgevert1;
	s32 *pleftedgevert2;
	s32 numrightedges;
	s32 *prightedgevert0;
	s32 *prightedgevert1;
	s32 *prightedgevert2;
} edgetable;
static s32 r_p0[6], r_p1[6], r_p2[6];
static u8 *d_pcolormap;
static s32 d_xdenom;
static edgetable *pedgetable;
static s32 a_sstepxfrac, a_tstepxfrac, r_lstepx, a_ststepxwhole;
static s32 r_sstepx, r_tstepx, r_lstepy, r_sstepy, r_tstepy;
static s32 r_zistepx, r_zistepy;
static s32 d_aspancount, d_countextrastep;
static spanpackage_t *a_spans;
static spanpackage_t *d_pedgespanpackage;
static u8 *d_pdest, *d_ptex;
static s16 *d_pz;
static s32 d_sfrac, d_tfrac, d_light, d_zi;
static s32 d_ptexextrastep, d_sfracextrastep;
static s32 d_tfracextrastep, d_lightextrastep, d_pdestextrastep;
static s32 d_lightbasestep, d_pdestbasestep, d_ptexbasestep;
static s32 d_sfracbasestep, d_tfracbasestep;
static s32 d_ziextrastep, d_zibasestep;
static s32 d_pzextrastep, d_pzbasestep;
static u8 *skintable[MAX_LBM_HEIGHT];
static s32 ubasestep, errorterm, erroradjustup, erroradjustdown;

edgetable edgetables[12] = {
	{ 0, 1, r_p0, r_p2, NULL, 2, r_p0, r_p1, r_p2 },
	{ 0, 2, r_p1, r_p0, r_p2, 1, r_p1, r_p2, NULL },
	{ 1, 1, r_p0, r_p2, NULL, 1, r_p1, r_p2, NULL },
	{ 0, 1, r_p1, r_p0, NULL, 2, r_p1, r_p2, r_p0 },
	{ 0, 2, r_p0, r_p2, r_p1, 1, r_p0, r_p1, NULL },
	{ 0, 1, r_p2, r_p1, NULL, 1, r_p2, r_p0, NULL },
	{ 0, 1, r_p2, r_p1, NULL, 2, r_p2, r_p0, r_p1 },
	{ 0, 2, r_p2, r_p1, r_p0, 1, r_p2, r_p0, NULL },
	{ 0, 1, r_p1, r_p0, NULL, 1, r_p1, r_p2, NULL },
	{ 1, 1, r_p2, r_p1, NULL, 1, r_p0, r_p1, NULL },
	{ 1, 1, r_p1, r_p0, NULL, 1, r_p2, r_p0, NULL },
	{ 0, 1, r_p0, r_p2, NULL, 1, r_p0, r_p1, NULL },
};

void D_PolysetDrawSpans(spanpackage_t * pspanpackage);
void D_PolysetCalcGradients(s32 skinwidth);
void D_DrawSubdiv();
void D_DrawNonSubdiv();
void D_PolysetRecursiveTriangle(s32 *p1, s32 *p2, s32 *p3);
void D_PolysetSetEdgeTable();
void D_RasterizeAliasPolySmooth();
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

void D_PolysetDraw()
{
	if(r_alphastyle.value == 0 && cur_ent_alpha != 1 && !fog_lut_built)
		R_BuildColorMixLUT(0);
	spanpackage_t spans[(MAXHEIGHT + 1) // +1 for spanpackage marking end
		+ 1 + ((CACHE_SIZE - 1) / sizeof(spanpackage_t)) + 1];
	a_spans = (spanpackage_t *)
		(((uintptr_t) & spans[0] + CACHE_SIZE - 1) & ~(CACHE_SIZE - 1));
	if(r_affinetridesc.drawtype)
		D_DrawSubdiv();
	else
		D_DrawNonSubdiv();
}

void D_PolysetDrawFinalVerts(finalvert_t *fv, s32 numverts)
{
	for(s32 i = 0; i < numverts; i++, fv++){
	// valid triangle coordinates for filling can include the bottom and
	// right clip edges, due to the fill rule; these shouldn't be drawn
		if(fv->v[0]<r_refdef.vrectright&&fv->v[1]<r_refdef.vrectbottom){
			if(fv->v[0] < 0 || fv->v[1] < 0)
				continue; // zbuf segfault
			s16 *zbuf = zspantable[fv->v[1]] + fv->v[0];
			if(zbuf - d_pzbuffer > (s32)(vid.width * vid.height
					*sizeof(s16)) || zbuf < d_pzbuffer){
				continue; // zbuf segfault
			}
			s32 z = fv->v[5] >> 16;
			if(!(z >= *zbuf))
				continue;
			s32 s = fv->v[2] >> 16;
			s32 t = fv->v[3] >> 16;
			if(t < 0 || t >= MAX_LBM_HEIGHT)
				continue;
			u8 *row = skintable[t];
			if(!row)
				continue;
			s32 pix = row[s];
			*zbuf = z;
			s32 texel = pix;
			pix = ((u8*)acolormap)[texel + (fv->v[4]&0xFF00)];
			if(r_rgblighting.value && colored_aliaslight
					&& texel < 0xE0)
				pix = D_GetRGBPix(pix);
			if(pix == 0xFF) continue;
			if(r_alphastyle.value == 0 && cur_ent_alpha != 1){
				s32 curpix = d_viewbuffer[d_scantable[fv->v[1]]
							+ fv->v[0]];
				d_viewbuffer[d_scantable[fv->v[1]] + fv->v[0]] =
					color_mix_lut[curpix][pix]
					[(s32)((cur_ent_alpha)*FOG_LUT_LEVELS)];
			}
			else if(r_alphastyle.value == 1 && cur_ent_alpha != 1){
				if(D_Dither(&d_viewbuffer[d_scantable[fv->v[1]]
						+ fv->v[0]], cur_ent_alpha))
					d_viewbuffer[d_scantable[fv->v[1]]
						+ fv->v[0]] = pix;
			}
			else d_viewbuffer[d_scantable[fv->v[1]]+fv->v[0]] = pix;
		}
	}
}

void D_DrawSubdiv()
{
	finalvert_t *pfv = r_affinetridesc.pfinalverts;
	mtriangle_t *ptri = r_affinetridesc.ptriangles;
	s32 lnumtriangles = r_affinetridesc.numtriangles;
	for(s32 i = 0; i < lnumtriangles; i++){
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
		}else{
			s32 s0, s1, s2;
			s0 = i0->v[2];
			s1 = i1->v[2];
			s2 = i2->v[2];
			if(i0->flags & ALIAS_ONSEAM)
				i0->v[2] += r_affinetridesc.seamfixupX16;
			if(i1->flags & ALIAS_ONSEAM)
				i1->v[2] += r_affinetridesc.seamfixupX16;
			if(i2->flags & ALIAS_ONSEAM)
				i2->v[2] += r_affinetridesc.seamfixupX16;
			D_PolysetRecursiveTriangle(i0->v, i1->v, i2->v);
			i0->v[2] = s0;
			i1->v[2] = s1;
			i2->v[2] = s2;
		}
	}
}

void D_DrawNonSubdiv()
{
	finalvert_t *pfv = r_affinetridesc.pfinalverts;
	mtriangle_t *ptri = r_affinetridesc.ptriangles;
	s32 lnumtriangles = r_affinetridesc.numtriangles;
	for(s32 i = 0; i < lnumtriangles; i++, ptri++){
		finalvert_t *i0 = pfv + ptri->vertindex[0];
		finalvert_t *i1 = pfv + ptri->vertindex[1];
		finalvert_t *i2 = pfv + ptri->vertindex[2];
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
		r_p0[0] = i0->v[0]; // u
		r_p0[1] = i0->v[1]; // v
		r_p0[2] = i0->v[2]; // s
		r_p0[3] = i0->v[3]; // t
		r_p0[4] = i0->v[4]; // light
		r_p0[5] = i0->v[5]; // iz
		r_p1[0] = i1->v[0];
		r_p1[1] = i1->v[1];
		r_p1[2] = i1->v[2];
		r_p1[3] = i1->v[3];
		r_p1[4] = i1->v[4];
		r_p1[5] = i1->v[5];
		r_p2[0] = i2->v[0];
		r_p2[1] = i2->v[1];
		r_p2[2] = i2->v[2];
		r_p2[3] = i2->v[3];
		r_p2[4] = i2->v[4];
		r_p2[5] = i2->v[5];
		if(!ptri->facesfront){
			if(i0->flags & ALIAS_ONSEAM)
				r_p0[2] += r_affinetridesc.seamfixupX16;
			if(i1->flags & ALIAS_ONSEAM)
				r_p1[2] += r_affinetridesc.seamfixupX16;
			if(i2->flags & ALIAS_ONSEAM)
				r_p2[2] += r_affinetridesc.seamfixupX16;
		}
		D_PolysetSetEdgeTable();
		D_RasterizeAliasPolySmooth();
	}
}

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
	s16 *zbuf = zspantable[new[1]] + new[0];
	if(zbuf - d_pzbuffer > (s32)(vid.width * vid.height * sizeof(s16))
			|| zbuf < d_pzbuffer)
		goto nodraw; // zbuf segfault
	s32 z = new[5] >> 16;
	if(z >= *zbuf){
		*zbuf = z;
		s32 s = new[2] >> 16;
		s32 t = new[3] >> 16;
		if(t < 0 || t >= MAX_LBM_HEIGHT)
			goto nodraw;
		u8 *row = skintable[t];
		if(!row)
			goto nodraw;
		s32 pix = row[s];
		s32 texel = pix;
		pix = d_pcolormap[texel];
		if(r_rgblighting.value && colored_aliaslight && texel < 0xE0)
			pix = D_GetRGBPix(pix);
		if(pix == 0xFF) goto nodraw;
		if(r_alphastyle.value == 0 && cur_ent_alpha != 1){
			s32 curpix = d_viewbuffer[d_scantable[new[1]] + new[0]];
			d_viewbuffer[d_scantable[new[1]] + new[0]] =
				color_mix_lut[curpix][pix]
				[(s32)((cur_ent_alpha)*FOG_LUT_LEVELS)];
		}
		else if(r_alphastyle.value == 1 && cur_ent_alpha != 1){
			if(D_Dither(&d_viewbuffer[d_scantable[new[1]] + new[0]],
						cur_ent_alpha))
				d_viewbuffer[d_scantable[new[1]]+new[0]] = pix;
		}
		else d_viewbuffer[d_scantable[new[1]] + new[0]] = pix;
	}
nodraw: // recursively continue
	D_PolysetRecursiveTriangle(lp3, lp1, new);
	D_PolysetRecursiveTriangle(lp3, new, lp2);
}

void D_PolysetUpdateTables()
{
	static u8 *skinstart;
	static s32 skinwidth;
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
	f64 dm = (f64)tm;
	f64 dn = (f64)tn;
	FloorDivMod(dm, dn, &ubasestep, &erroradjustup);
	erroradjustdown = dn;
}

void D_PolysetCalcGradients(s32 skinwidth)
{
	f32 p00_minus_p20 = r_p0[0] - r_p2[0];
	f32 p01_minus_p21 = r_p0[1] - r_p2[1];
	f32 p10_minus_p20 = r_p1[0] - r_p2[0];
	f32 p11_minus_p21 = r_p1[1] - r_p2[1];
	f32 xstepdenominv = 1.0 / (f32)d_xdenom;
	f32 ystepdenominv = -xstepdenominv;
// ceil() for light so positive steps are exaggerated, negative steps
// diminished, pushing us away from underflow toward overflow. Underflow is
// very visible, overflow is very unlikely, because of ambient lighting
	f32 t0 = r_p0[4] - r_p2[4];
	f32 t1 = r_p1[4] - r_p2[4];
	r_lstepx = (s32)ceil((t1*p01_minus_p21-t0*p11_minus_p21)*xstepdenominv);
	r_lstepy = (s32)ceil((t1*p00_minus_p20-t0*p10_minus_p20)*ystepdenominv);
	t0 = r_p0[2] - r_p2[2];
	t1 = r_p1[2] - r_p2[2];
	r_sstepx = (s32)((t1*p01_minus_p21-t0*p11_minus_p21)*xstepdenominv);
	r_sstepy = (s32)((t1*p00_minus_p20-t0*p10_minus_p20)*ystepdenominv);
	t0 = r_p0[3] - r_p2[3];
	t1 = r_p1[3] - r_p2[3];
	r_tstepx = (s32)((t1*p01_minus_p21-t0*p11_minus_p21)*xstepdenominv);
	r_tstepy = (s32)((t1*p00_minus_p20-t0*p10_minus_p20)*ystepdenominv);
	t0 = r_p0[5] - r_p2[5];
	t1 = r_p1[5] - r_p2[5];
	r_zistepx = (s32)((t1*p01_minus_p21-t0*p11_minus_p21)*xstepdenominv);
	r_zistepy = (s32)((t1*p00_minus_p20-t0*p10_minus_p20)*ystepdenominv);
	a_sstepxfrac = r_sstepx & 0xFFFF;
	a_tstepxfrac = r_tstepx & 0xFFFF;
	a_ststepxwhole = skinwidth * (r_tstepx >> 16) + (r_sstepx >> 16);
}

static polysetstyle_t D_PolysetPixelStyle()
{
	polysetstyle_t st;
	st.colormap = (u8 *)acolormap;
	st.rgblight = r_rgblighting.value && colored_aliaslight;
	st.rgbpix = st.rgblight ? D_RGBPixTable() : NULL;
	st.mixlevel = 0;
	if(cur_ent_alpha == 1)
		st.blend = BLEND_OPAQUE;
	else if(r_alphastyle.value == 0){
		st.blend = BLEND_MIX;
		st.mixlevel = (s32)(cur_ent_alpha * FOG_LUT_LEVELS);
	}else if(r_alphastyle.value == 1)
		st.blend = BLEND_STIPPLE;
	else
		st.blend = BLEND_OPAQUE;
	return st;
}
 

Q_SIL u8 D_PolysetDitheredTexel(const u8 *ptex, s32 sfrac, s32 tfrac,
		const u8 *pdest)
{ // CyanBun96: dithered sampling from Unreal
	const u8 *skin = (u8 *)r_affinetridesc.pskin;
	s32 skinw = r_affinetridesc.skinwidth;
	s32 skinh = r_affinetridesc.skinheight;
	ptrdiff_t off = ptex - skin;
	if(off < 0) off = 0;
	s32 s = ((s32)(off % skinw) << 16) | (sfrac & 0xFFFF);
	s32 t = ((s32)(off / skinw) << 16) | (tfrac & 0xFFFF);
	s32 pixel = (s32)(pdest - (u8 *)screen->pixels);
	s32 y = pixel / vid.width;
	s32 x = pixel - y * vid.width;
	s32 dither = (x & 1) + ((y & 1) << 1);
	s += dither_s[dither];
	t += dither_t[dither];
	s32 smax = (skinw - 1) << 16;
	s32 tmax = (skinh - 1) << 16;
	if(s < 0) s = 0;
	if(t < 0) t = 0;
	if(s >= smax) s = smax;
	if(t >= tmax) t = tmax;
	u8 texel = skin[(s >> 16) + (t >> 16) * skinw];
// The single shotgun model would sometimes get bright blue texels, which should
// be unused as it's the "transparent" color for model skins
// Same with 0xd3 on the torch model
	if(texel == 0xd0 || texel == 0xd3)
		texel = 0;
	return texel;
}
 
Q_SIL void D_PolysetPlot(const polysetstyle_t *st, u8 texel, s32 light,
		u8 *pdest, s16 *pz, s32 zi)
{
	s32 pix = st->colormap[texel + (light & 0xFF00)];
	if(pix == 0xFF)
		return;
	if(st->rgblight && texel < 0xE0)
		pix = st->rgbpix[pix];
	switch(st->blend){
	case BLEND_MIX:
		*pdest = color_mix_lut[*pdest][pix][st->mixlevel];
		break;
	case BLEND_STIPPLE:
		if(D_Dither(pdest, cur_ent_alpha))
			*pdest = pix;
		break;
	default:
		*pdest = pix;
	}
	*pz = zi >> 16;
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
		const polysetstyle_t *st, s32 sample)
{
	u8 *pdest = span->pdest;
	s16 *pz = span->pz;
	u8 *ptex = span->ptex;
	s32 sfrac = span->sfrac;
	s32 tfrac = span->tfrac;
	s32 light = span->light;
	s32 lstep = r_lstepx;
	s32 zi = span->zi;
	D_PolysetClampSpanLight(&light, &lstep, count);
	do{
		if((zi >> 16) >= *pz){
			u8 texel = sample == TEXSAMPLE_DITHERED ?
			  D_PolysetDitheredTexel(ptex,sfrac,tfrac,pdest):*ptex;
			D_PolysetPlot(st, texel, light, pdest, pz, zi);
		}
		pdest++;
		pz++;
		zi += r_zistepx;
		light += lstep;
		ptex += a_ststepxwhole;
		sfrac += a_sstepxfrac;
		ptex += sfrac >> 16;
		sfrac &= 0xFFFF;
		tfrac += a_tstepxfrac;
		if(tfrac & 0x10000){
			ptex += r_affinetridesc.skinwidth;
			tfrac &= 0xFFFF;
		}
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
	const s16 *pz = span->pz;
	const u8 *dest_end = d_viewbuffer + vid.width * vid.height;
	const s16 *z_end = d_pzbuffer + vid.width * vid.height;
	return dest >= d_viewbuffer && dest + count <= dest_end &&
		pz >= d_pzbuffer && pz + count <= z_end;
}
 
Q_SIL void D_PolysetDrawSpansImpl(spanpackage_t *pspanpackage, s32 sample)
{
	const polysetstyle_t st = D_PolysetPixelStyle();
	do{
		s32 count = D_PolysetNextSpanWidth(pspanpackage);
		if(count > 0 && D_PolysetSpanInBounds(pspanpackage, count))
			D_PolysetDrawSpan(pspanpackage, count, &st, sample);
		pspanpackage++;
	}while(pspanpackage->count != -999999);
}
 
void D_PolysetDrawSpans(spanpackage_t *pspanpackage)
{ D_PolysetDrawSpansImpl(pspanpackage, TEXSAMPLE_NEAREST); }
 
void D_PolysetDrawSpansDithered(spanpackage_t *pspanpackage)
{ D_PolysetDrawSpansImpl(pspanpackage, TEXSAMPLE_DITHERED); }

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

void D_RasterizeAliasPolySmooth()
{
	s32 *plefttop = pedgetable->pleftedgevert0;
	s32 *prighttop = pedgetable->prightedgevert0;
	s32 *pleftbottom = pedgetable->pleftedgevert1;
	s32 *prightbottom = pedgetable->prightedgevert1;
	u32 initialleftheight = pleftbottom[1] - plefttop[1];
	u32 initialrightheight = prightbottom[1] - prighttop[1];
	if(initialleftheight>=vid.height || initialrightheight>=vid.height)
		return; // broken alias model
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
		pleftbottom = pedgetable->pleftedgevert2;
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
	if(r_dithertex.value) D_PolysetDrawSpansDithered(a_spans);
	else D_PolysetDrawSpans(a_spans);
	// scan out the bottom part of the right edge, if it exists
	if(pedgetable->numrightedges == 2){
		spanpackage_t *pstart;
		pstart = a_spans + initialrightheight;
		pstart->count = originalcount;
		d_aspancount = prightbottom[0] - prighttop[0];
		prighttop = prightbottom;
		prightbottom = pedgetable->prightedgevert2;
		s32 height = prightbottom[1] - prighttop[1];
		D_PolysetSetUpForLineScan(prighttop[0], prighttop[1],
				prightbottom[0], prightbottom[1]);
		d_countextrastep = ubasestep + 1;
		a_spans[initialrightheight + height].count = -999999;
		if(r_dithertex.value) D_PolysetDrawSpansDithered(pstart);
		else D_PolysetDrawSpans(pstart);//mark end of the spanpackages
	}
}

void D_PolysetSetEdgeTable()
{
	s32 edgetableindex = 0; // assume the vertices are already in top to
		// bottom order determine which edges are right & left, and the
		// order in which to rasterize them
	if(r_p0[1] >= r_p1[1]){
		if(r_p0[1] == r_p1[1]){
		    pedgetable=r_p0[1]<r_p2[1]?&edgetables[2]:&edgetables[5];
			return;
		}else
			edgetableindex = 1;
	}
	if(r_p0[1] == r_p2[1]){
		pedgetable = edgetableindex ? &edgetables[8] : &edgetables[9];
		return;
	}else if(r_p1[1] == r_p2[1]){
		pedgetable = edgetableindex ? &edgetables[10] : &edgetables[11];
		return;
	}
	if(r_p0[1] > r_p2[1])
		edgetableindex += 2;
	if(r_p1[1] > r_p2[1])
		edgetableindex += 4;
	pedgetable = &edgetables[edgetableindex];
}
