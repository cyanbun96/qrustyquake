// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// d_scan.c: Portable C scan-level rasterization code
#include "quakedef.h"

typedef struct {
	u8 *pdest;
	s16 *zp;
	s32 iz;
	s32 s, t;
	s32 sstep, tstep;
	s32 count;
} turbchunk_t;
typedef struct {
	s32 shift;
	s32 count;
	s32 len;
	f32 sdivz, tdivz, zi;
	f32 sdivz_chunk, tdivz_chunk, zi_chunk;
	s32 s, t;
	s32 snext, tnext;
	s32 sstep, tstep;
} spanwalk_t;
typedef struct {
	const u8 *pbase;
	s32 width;
	s32 s, t;
	s32 sstep, tstep;
	s32 s_max, t_max;
	s32 x, y;
	s32 count;
} texcursor_t;

Q_SIL void D_TurbChunk(turbchunk_t *c, s32 trans, s32 filtered, s32 lit,s32 hl);
typedef void(*turbspanfunc_t)(turbchunk_t *);
#define DIST_LUT_SIZE 1024
#define COS_LUT_SIZE 2048
enum { TURBTRANS_OPAQUE, TURBTRANS_MIXED, TURBTRANS_DITHERED };
#define TURBFUNC_FILTERED 3
#define TURBFUNC_LIT 6
#define TURBFUNC_HL 12

#define TURBFUNC_INDEX(trans, filtered, lit, hl) \
((trans)+TURBFUNC_FILTERED*(filtered)+TURBFUNC_LIT*(lit)+TURBFUNC_HL*(hl))

#define TURB_VARIANTS(X) \
	X(OPAQUE, 0, 0, 0) X(MIXED, 0, 0, 0) X(DITHERED, 0, 0, 0) \
	X(OPAQUE, 1, 0, 0) X(MIXED, 1, 0, 0) X(DITHERED, 1, 0, 0) \
	X(OPAQUE, 0, 1, 0) X(MIXED, 0, 1, 0) X(DITHERED, 0, 1, 0) \
	X(OPAQUE, 1, 1, 0) X(MIXED, 1, 1, 0) X(DITHERED, 1, 1, 0) \
	X(OPAQUE, 0, 0, 1) X(MIXED, 0, 0, 1) X(DITHERED, 0, 0, 1) \
	X(OPAQUE, 1, 0, 1) X(MIXED, 1, 0, 1) X(DITHERED, 1, 0, 1) \
	X(OPAQUE, 0, 1, 1) X(MIXED, 0, 1, 1) X(DITHERED, 0, 1, 1) \
	X(OPAQUE, 1, 1, 1) X(MIXED, 1, 1, 1) X(DITHERED, 1, 1, 1)

#define TURB_DEFINE_CHUNKFUNC(tr, filtered, lit, hl) \
	Q_SIL void D_TurbChunk_##tr##_##filtered##_##lit##_##hl(turbchunk_t *c)\
{ D_TurbChunk(c, TURBTRANS_##tr, filtered, lit, hl); }

TURB_VARIANTS(TURB_DEFINE_CHUNKFUNC)

#define TURB_TABLE_ENTRY(tr, filtered, lit, hl) \
		[TURBFUNC_INDEX(TURBTRANS_##tr, filtered, lit, hl)] = \
	D_TurbChunk_##tr##_##filtered##_##lit##_##hl,
	static turbspanfunc_t turbdrawfunc[24] = {
		TURB_VARIANTS(TURB_TABLE_ENTRY)
	};

static u8 *r_turb_pbase;
static s32 *r_turb_turb;
static s16 *pz; // Manoel Kasimier - translucent water
static s32 izi, izistep;
static u8 cutoutbuf[MAXHEIGHT*MAXWIDTH];
static f32 dist_lut[DIST_LUT_SIZE];
static f32 cos_lut[COS_LUT_SIZE];
static f32 dist_lut_max = 8.0f; // max dx*dx + dy*dy
static f32 cos_lut_period = 2.0f * M_PI;
static f32 turb_opacity = 0;

void D_WarpScreen() // this performs a slight compression of the screen at the
{ // same time as the sine warp, to keep the edges from wrapping
	s32 w = r_refdef.vrect.width;
	s32 h = r_refdef.vrect.height;
	f32 wratio = w / (f32)scr_vrect.width;
	f32 hratio = h / (f32)scr_vrect.height;
	u8 *rowptr[MAXHEIGHT + (TURB_AMP2 * 2)];
	s32 column[MAXWIDTH + (TURB_AMP2 * 2)];
	for(s32 v = 0; v < scr_vrect.height + TURB_AMP2 * 2; v++)
		rowptr[v] = d_viewbuffer + (r_refdef.vrect.y * screenwidth) +
		screenwidth * (s32)((f32)v * hratio * h / (h + TURB_AMP2 * 2));
	for(s32 u = 0; u < scr_vrect.width + TURB_AMP2 * 2; u++)
		column[u] = r_refdef.vrect.x + (s32)((f32)u * wratio * w /
				(w + TURB_AMP2 * 2));
	s32 *turb = intsintable + ((s32)(cl.time*TURB_SPEED) & (TURB_CYCLE-1));
	u8 *dest = vid.buffer + scr_vrect.y * vid.width + scr_vrect.x;
	for(s32 v = 0; v < scr_vrect.height; v++, dest += vid.width){
		s32 *col = &column[turb[v]];
		u8 **row = &rowptr[v];
		for(s32 u = 0; u < scr_vrect.width; u += 4){
			dest[u + 0] = row[turb[u + 0]][col[u + 0]];
			dest[u + 1] = row[turb[u + 1]][col[u + 1]];
			dest[u + 2] = row[turb[u + 2]][col[u + 2]];
			dest[u + 3] = row[turb[u + 3]][col[u + 3]];
		}
	}
}

Q_SIL f32 fast_sqrt_dist(f32 x)
{
	f32 f = x * (f32)(DIST_LUT_SIZE - 1) / dist_lut_max;
	s32 i = (s32)f;
	if(i < 0) i = 0;
	if(i >= DIST_LUT_SIZE - 1) i = DIST_LUT_SIZE - 1;
	return dist_lut[i];
}

Q_SIL f32 fast_cos(f32 theta)
{
	f32 n = theta * (f32)COS_LUT_SIZE / cos_lut_period;
	s32 i = (s32)n;
	i &= (COS_LUT_SIZE - 1);
	return cos_lut[i];
}

static void InitHLWaterLUT()
{
	static s32 HLWaterLutInited = 0;
	if(HLWaterLutInited) return;
	for(s32 i = 0; i < DIST_LUT_SIZE; i++){
		f32 x = (f32)i * (dist_lut_max / (DIST_LUT_SIZE - 1));
		dist_lut[i] = sqrtf(x);
	}
	for(s32 i = 0; i < COS_LUT_SIZE; i++){
		f32 t = (f32)i * (cos_lut_period / COS_LUT_SIZE);
		cos_lut[i] = cosf(t);
	}
	HLWaterLutInited = 1;
}

static void HLWarp(f32 *u, f32 *v, f32 time)
{
	f32 x = *u;
	f32 y = *v;
	f32 waveSize = (1.0f / 4.0f) * r_hlwavescale.value;
	f32 waveStrength = 0.02f * r_hlripplescale.value;
	f32 t = time * 3.0f;
	f32 arg1 = y * (1.0f / waveSize) + t;
	f32 arg2 = x * (1.0f / waveSize) + t;
	f32 xu = x + fast_cos(arg1) * waveStrength;
	f32 yv = y + fast_cos(arg2) * waveStrength;
	*u = xu;
	*v = yv;
}

static void HLSingleRipple(f32 *outx, f32 *outy, f32 u, f32 v, f32 px, f32 py,
		f32 size, f32 freq, f32 strength, f32 t)
{
	f32 dx = u - px;
	f32 dy = v - py;
	f32 r2 = dx*dx + dy*dy;
	f32 dist = fast_sqrt_dist(r2);
	f32 f = fast_sqrt_dist(dist) * freq - t;
	f32 wave = fast_cos(f);
	f32 decay = 1.0f - dist / size;
	if(decay < 0.0f)
		decay = 0.0f;
	wave *= decay;
	*outx += dx * wave * strength;
	*outy += dy * wave * strength;
}

static void HLRipple(f32 *u, f32 *v, f32 time)
{
	static const f32 px[3] = {0.2f, 0.1f, 0.7f};
	static const f32 py[3] = {0.8f, 0.2f, 0.5f};
	static const f32 tile[8][2] = {
		{ 1,  0}, { 1,  1}, { 0,  1}, {-1,  1},
		{-1,  0}, {-1, -1}, { 0, -1}, { 1, -1}
	};
	f32 iu = floorf(*u);
	f32 iv = floorf(*v);
	f32 baseu = *u;
	f32 basev = *v;
	f32 ox = 0.0f;
	f32 oy = 0.0f;
	f32 size = 0.8f;
	f32 freq = 200.0f;
	f32 strength = 0.04f * r_hlripplescale.value;
	f32 t = time * 20.0f;
	for(s32 i = 0; i < 3; i++)
		HLSingleRipple(&ox, &oy, baseu, basev, iu + px[i], iv + py[i],
				size, freq, strength, t);
	if(r_hlwaterquality.value)
		for(s32 ti = 0; ti < 8; ti += 2)
		for(s32 i = 0; i < 3; i++)
		HLSingleRipple(&ox, &oy, baseu+tile[ti][0], basev+tile[ti][1],
				iu+px[i], iv+py[i], size, freq, strength, t);

	*u = baseu + ox;
	*v = basev + oy;
}

static inline u8 D_TurbMixLit(u8 pix, u8 lit)
{
	u8 rp = CURWORLDPAL[pix * 3 + 0];
	u8 gp = CURWORLDPAL[pix * 3 + 1];
	u8 bp = CURWORLDPAL[pix * 3 + 2];
	u8 rl = CURWORLDPAL[lit * 3 + 0];
	u8 gl = CURWORLDPAL[lit * 3 + 1];
	u8 bl = CURWORLDPAL[lit * 3 + 2];
	s32 r = rp * rl / 255;
	s32 g = gp * gl / 255;
	s32 b = bp * bl / 255;
	return lit_lut[QUANT(r)][QUANT(g)][QUANT(b)];
}

static inline s32 D_ClampCoord(s32 v, s32 lo, s32 hi)
{ if(v > hi) return hi; if(v < lo) return lo; return v; }

Q_SIL void D_SpanWalkBegin(spanwalk_t *w, const espan_t *pspan, s32 shift)
{
	f32 du = (f32)pspan->u;
	f32 dv = (f32)pspan->v;
	w->shift = shift;
	w->count = pspan->count;
	w->len = 0;
	w->sdivz = d_sdivzorigin + dv*d_sdivzstepv + du*d_sdivzstepu;
	w->tdivz = d_tdivzorigin + dv*d_tdivzstepv + du*d_tdivzstepu;
	w->zi = d_ziorigin + dv * d_zistepv + du * d_zistepu;
	w->sdivz_chunk = d_sdivzstepu * (1 << shift);
	w->tdivz_chunk = d_tdivzstepu * (1 << shift);
	w->zi_chunk = d_zistepu * (1 << shift);
	f32 z = (f32)0x10000 / w->zi;
	w->s = D_ClampCoord((s32)(w->sdivz * z) + sadjust, 0, bbextents);
	w->t = D_ClampCoord((s32)(w->tdivz * z) + tadjust, 0, bbextentt);
	w->snext = w->s;
	w->tnext = w->t;
	w->sstep = 0;
	w->tstep = 0;
}

Q_SIL void D_SpanWalkEvalEnd(spanwalk_t *w)
{
	s32 lo = 1 << w->shift;
	f32 z = (f32)0x10000 / w->zi;
	w->snext = D_ClampCoord((s32)(w->sdivz * z) + sadjust, lo, bbextents);
	w->tnext = D_ClampCoord((s32)(w->tdivz * z) + tadjust, lo, bbextentt);
}

Q_SIL void D_SpanWalkNextChunk(spanwalk_t *w)
{
	w->len = w->count;
	if(w->len >= (1 << w->shift))
		w->len = 1 << w->shift;
	w->count -= w->len;
	w->sstep = 0;
	w->tstep = 0;
	if(w->count){
		w->sdivz += w->sdivz_chunk;
		w->tdivz += w->tdivz_chunk;
		w->zi += w->zi_chunk;
		D_SpanWalkEvalEnd(w);
		w->sstep = (w->snext - w->s) >> w->shift;
		w->tstep = (w->tnext - w->t) >> w->shift;
	}else{
		f32 lastpixel = (f32)(w->len - 1);
		w->sdivz += d_sdivzstepu * lastpixel;
		w->tdivz += d_tdivzstepu * lastpixel;
		w->zi += d_zistepu * lastpixel;
		D_SpanWalkEvalEnd(w);
		if(w->len > 1){
			w->sstep = (w->snext - w->s) / (w->len - 1);
			w->tstep = (w->tnext - w->t) / (w->len - 1);
		}
	}
}

Q_SIL void D_SpanWalkEndChunk(spanwalk_t *w)
{ w->s = w->snext; w->t = w->tnext; }

Q_SIL void D_SpanDepthSetup()
{ izistep = (s32)(d_zistepu * 0x8000 * 0x10000); }

Q_SIL s16 *D_SpanDepthPtr(const espan_t *pspan)
{ return d_pzbuffer + (d_zwidth * pspan->v) + pspan->u; }

Q_SIL s32 D_SpanDepthStart(f32 zi)
{ return (s32)(zi * 0x8000 * 0x10000); }

Q_SIL void D_SpanDepthBegin(const espan_t *pspan, f32 zi)
{ pz = D_SpanDepthPtr(pspan); izi = D_SpanDepthStart(zi); }

Q_SIL s32 D_SpanDepthPass()
{ return *pz <= (izi >> 16); }

Q_SIL void D_SpanDepthStep()
{ izi += izistep; pz++; }

Q_SIL s32 D_DitherIndex(s32 x, s32 y)
{ return(x & 1) + ((y & 1) << 1); }

Q_SIL void D_PixelScreenPos(const u8 *pdest, s32 *x, s32 *y)
{
	s32 pixel_index = (s32)(pdest - (u8 *)screen->pixels);
	*y = pixel_index / scr_vrect.width;
	*x = pixel_index - *y * scr_vrect.width;
}

Q_SIL u8 D_TurbSample(const u8 *pbase, const s32 *turb, f32 time,
		s32 s, s32 t, s32 x, s32 y, s32 filtered, s32 hl)
{
	if(filtered){
		s32 dither_idx = D_DitherIndex(x, y);
		s += dither_s[dither_idx];
		t += dither_t[dither_idx];
	}
	if(hl){
		f32 u = (f32)((s >> 16) & 63) * (1.0f/64.0f);
		f32 v = (f32)((t >> 16) & 63) * (1.0f/64.0f);
		HLWarp(&u, &v, time);
		HLRipple(&u, &v, time);
		s32 su = ((s32)(u * 64.0f)) & 63;
		s32 tv = ((s32)(v * 64.0f)) & 63;
		return *(pbase + (tv << 6) + su);
	}
	s32 sw = ((s + turb[(t>>16)&(TURB_CYCLE-1)]) >> 16) & 63;
	s32 tw = ((t + turb[(s>>16)&(TURB_CYCLE-1)]) >> 16) & 63;
	return *(pbase + (tw << 6) + sw);
}

Q_SIL void D_TurbChunk(turbchunk_t *c, s32 trans, s32 filtered, s32 lit, s32 hl)
{
	u8 *pdest = c->pdest;
	s16 *zp = c->zp;
	s32 iz = c->iz;
	s32 s = c->s;
	s32 t = c->t;
	s32 count = c->count;
	s32 x = 0;
	s32 y = 0;
	const s32 sstep = c->sstep;
	const s32 tstep = c->tstep;
	const s32 izstep = izistep;
	const u8 *pbase = r_turb_pbase;
	const s32 *turb = r_turb_turb;
	const u8 *litbase = (u8 *)litwater_base;
	const f32 time = cl.time;
	const f32 ditherthresh = 1 - turb_opacity;
	const s32 foglut = (trans == TURBTRANS_MIXED) ?
		(s32)(turb_opacity*FOG_LUT_LEVELS) : 0;
	if(filtered)
		D_PixelScreenPos(pdest, &x, &y);
	do{
		s32 draw = 1;
		if(trans != TURBTRANS_OPAQUE){
			draw = *zp <= (iz >> 16);
			if(draw && trans == TURBTRANS_DITHERED)
				draw = D_Dither(pdest, ditherthresh);
		}
		if(draw){
			u8 pix = D_TurbSample(pbase, turb, time, s, t, x, y,
					filtered, hl);
			if(lit)
			    pix = D_TurbMixLit(pix,litbase[pdest-d_viewbuffer]);
			if(trans == TURBTRANS_MIXED)
				pix = color_mix_lut[pix][*pdest][foglut];
			*pdest = pix;
		}
		pdest++;
		if(trans != TURBTRANS_OPAQUE){
			iz += izstep;
			zp++;
		}
		s += sstep;
		t += tstep;
		x++;
	}while(--count > 0);
	c->pdest = pdest;
	c->zp = zp;
	c->iz = iz;
}

static s32 D_TurbFuncIndex()
{
	s32 trans = TURBTRANS_OPAQUE;
	if(r_alphapass)
		trans = r_alphastyle.value ? TURBTRANS_DITHERED:TURBTRANS_MIXED;
	return TURBFUNC_INDEX(trans, r_dithertex.value != 0, lmonly != 0,
			r_hlwater.value != 0);
}

static turbspanfunc_t D_TurbBegin(f32 opacity)
{
	turb_opacity = opacity;
	turbspanfunc_t drawspan = turbdrawfunc[D_TurbFuncIndex()];
	if(r_hlwater.value)
		InitHLWaterLUT();
	if(lmonly)
		R_BuildLitLUT();
	if(r_alphastyle.value == 0 && opacity != 1 && opacity)
		R_BuildColorMixLUT(0);
	r_turb_turb = sintable + ((s32)(cl.time*TURB_SPEED) & (TURB_CYCLE-1));
	r_turb_pbase = (u8 *)cacheblock;
	D_SpanDepthSetup();
	return drawspan;
}

Q_SIL void D_TurbSetChunk(turbchunk_t *c, const spanwalk_t *w)
{
	c->s = w->s & ((TURB_CYCLE << 16) - 1);
	c->t = w->t & ((TURB_CYCLE << 16) - 1);
	c->sstep = w->sstep;
	c->tstep = w->tstep;
	c->count = w->len;
}

static void D_TurbSpan(const espan_t *pspan, turbspanfunc_t drawspan)
{
	spanwalk_t w;
	turbchunk_t c;
	c.pdest = (u8*)((u8*)d_viewbuffer + (screenwidth*pspan->v) + pspan->u);
	D_SpanWalkBegin(&w, pspan, 4); // exact coords every 16px
	c.zp = D_SpanDepthPtr(pspan);
	c.iz = D_SpanDepthStart(w.zi);
	do{
		D_SpanWalkNextChunk(&w);
		D_TurbSetChunk(&c, &w);
		drawspan(&c);
		D_SpanWalkEndChunk(&w);
	}while(w.count > 0);
}

void Turbulent(espan_t *pspan, f32 opacity)
{
	turbspanfunc_t drawspan = D_TurbBegin(opacity);
	do{
		D_TurbSpan(pspan, drawspan);
	}while((pspan = pspan->pnext) != NULL);
}

Q_SIL void D_CursorBegin(texcursor_t *c, const spanwalk_t *w,
		const u8 *pdest, s32 dithered)
{
	c->pbase = (u8 *)cacheblock;
	c->width = cachewidth;
	c->s = w->s;
	c->t = w->t;
	c->sstep = w->sstep;
	c->tstep = w->tstep;
	c->count = w->len;
	c->s_max = 0;
	c->t_max = 0;
	c->x = 0;
	c->y = 0;
	if(dithered){
		c->s_max = (cachewidth - 1) << 16;
		c->t_max = (cacheheight - 1) << 16;
		D_PixelScreenPos(pdest, &c->x, &c->y);
	}
}

Q_SIL u8 D_CursorTexel(const texcursor_t *c, s32 dithered)
{
	s32 s = c->s;
	s32 t = c->t;
	if(dithered){
		// CyanBun96: dithered sampling from Unreal
		s32 idx = D_DitherIndex(c->x, c->y);
		s = D_ClampCoord(s + dither_s[idx], 0, c->s_max);
		t = D_ClampCoord(t + dither_t[idx], 0, c->t_max);
	}
	return c->pbase[(s >> 16) + (t >> 16) * c->width];
}

Q_SIL void D_CursorAdvance(texcursor_t *c)
{
	c->s += c->sstep;
	c->t += c->tstep;
	c->x++;
}

Q_SIL u8 *D_SpanChunkNormal(texcursor_t *c, u8 *pdest, s32 dithered)
{
	do{
		*pdest++ = D_CursorTexel(c, dithered);
		D_CursorAdvance(c);
	}while(--c->count > 0);
	return pdest;
}

Q_SIL u8 *D_SpanChunkCutout(texcursor_t *c, u8 *pdest, s32 dithered)
{
	do{
		if(D_SpanDepthPass()){
			u8 pix = D_CursorTexel(c, dithered);
			cutoutbuf[pdest-d_viewbuffer] = 0;
			if(pix != 0xff){
				*pdest = pix;
				cutoutbuf[pdest-d_viewbuffer] = 1;
			}
		}
		pdest++;
		D_SpanDepthStep();
		D_CursorAdvance(c);
	}while(--c->count > 0);
	return pdest;
}

Q_SIL u8 *D_SpanChunkSkybox(texcursor_t *c, u8 *pdest, s32 dithered)
{
	s32 foglut = r_skyfog.value*FOG_LUT_LEVELS;
	do{
		u8 pix = D_CursorTexel(c, dithered);
		if(fog_density > 0)
			pix = color_mix_lut[pix][fog_pal_index][foglut];
		*pdest++ = pix;
		D_CursorAdvance(c);
	}while(--c->count > 0);
	return pdest;
}

Q_SIL u8 *D_SpanChunkBlend(texcursor_t *c, u8 *pdest, f32 opacity,
		s32 dithered)
{
	s32 foglut = opacity*FOG_LUT_LEVELS;
	if(!fog_lut_built)
		R_BuildColorMixLUT(0);
	do{
		if(D_SpanDepthPass()){
			u8 pix = D_CursorTexel(c, dithered);
			if(pix != 0xff)
				*pdest = color_mix_lut[pix][*pdest][foglut];
		}
		pdest++;
		D_SpanDepthStep();
		D_CursorAdvance(c);
	}while(--c->count > 0);
	return pdest;
}

Q_SIL u8 *D_SpanChunkStipple(texcursor_t *c, u8 *pdest, f32 opacity,
		s32 dithered)
{
	do{
		if(D_SpanDepthPass() && D_Dither(pdest, 1-opacity)){
			u8 pix = D_CursorTexel(c, dithered);
			if(pix != 0xff)
				*pdest = pix;
		}
		pdest++;
		D_SpanDepthStep();
		D_CursorAdvance(c);
	}while(--c->count > 0);
	return pdest;
}

Q_SIL u8 *D_DrawSpanChunk(const spanwalk_t *w, u8 *pdest, s32 type,
		f32 opacity, s32 dithered)
{
	texcursor_t c;
	D_CursorBegin(&c, w, pdest, dithered);
	if(type == SPAN_NORMAL)
		return D_SpanChunkNormal(&c, pdest, dithered);
	if(type == SPAN_CUTOUT)
		return D_SpanChunkCutout(&c, pdest, dithered);
	if(type == SPAN_SKYBOX)
		return D_SpanChunkSkybox(&c, pdest, dithered);
	if(type == SPAN_TRANS){
		if(r_alphastyle.value == 0)
			return D_SpanChunkBlend(&c, pdest, opacity, dithered);
		return D_SpanChunkStipple(&c, pdest, opacity, dithered);
	}
	return pdest;
}

static void D_AllocLitWater()
{
	litwater_base = malloc(vid.width * vid.height);
	if(!litwater_base)
		Sys_Error("Not enough memory for lit water");
}

static u8 *D_SpanDest(const espan_t *pspan)
{
	u8 *base = (u8 *)d_viewbuffer;
	if(lmonly){
		if(!litwater_base)
			D_AllocLitWater();
		base = (u8 *)litwater_base;
	}
	return base + (screenwidth * pspan->v) + pspan->u;
}

Q_SIL void D_DrawSpansImpl(espan_t *pspan, s32 type, f32 opacity, s32 dithered)
{
	s32 depthtested = (type == SPAN_TRANS || type == SPAN_CUTOUT);
	D_SpanDepthSetup();
	do{
		spanwalk_t w;
		u8 *pdest = D_SpanDest(pspan);
		D_SpanWalkBegin(&w, pspan, 3); // exact coords every 8px
		if(depthtested)
			D_SpanDepthBegin(pspan, w.zi);
		do{
			D_SpanWalkNextChunk(&w);
			pdest = D_DrawSpanChunk(&w, pdest, type, opacity,
					dithered);
			D_SpanWalkEndChunk(&w);
		}while(w.count > 0);
	}while((pspan = pspan->pnext) != NULL);
}

void D_DrawSpans(espan_t *pspan, s32 type, f32 opacity)
{ D_DrawSpansImpl(pspan, type, opacity, 0); }

void D_DrawSpansDithered(espan_t *pspan, s32 type, f32 opacity)
{ D_DrawSpansImpl(pspan, type, opacity, 1); }

void D_DrawZSpans(espan_t *pspan)
{
	s32 izistep = (s32)(d_zistepu * 0x8000 * 0x10000);
	do{
		s16 *pdest = d_pzbuffer + (d_zwidth * pspan->v) + pspan->u;
		s32 count = pspan->count;
		f32 du = (f32)pspan->u; // calculate the initial 1/z
		f32 dv = (f32)pspan->v;
		f64 zi = d_ziorigin + dv * d_zistepv + du * d_zistepu;
		s32 izi = (s32)(zi * 0x8000 * 0x10000);
		if((intptr_t) pdest & 0x02){
			*pdest++ = (s16)(izi >> 16);
			izi += izistep;
			count--;
		}
		s32 doublecount = count >> 1;
		if(doublecount > 0){
			do{
				u32 ltemp = izi >> 16;
				izi += izistep;
				ltemp |= izi & 0xFFFF0000;
				izi += izistep;
				*(s32 *)pdest = ltemp;
				pdest += 2;
			}while(--doublecount > 0);
		}
		if(count & 1)
			*pdest = (s16)(izi >> 16);
	}while((pspan = pspan->pnext) != NULL);
}

void D_DrawZSpansTrans(espan_t *pspan)
{
	s32 izistep = (s32)(d_zistepu * 0x8000 * 0x10000);
	do{
		s16 *pdest = d_pzbuffer + (d_zwidth * pspan->v) + pspan->u;
		s32 count = pspan->count;
		f32 du = (f32)pspan->u; // calculate the initial 1/z
		f32 dv = (f32)pspan->v;
		f64 zi = d_ziorigin + dv * d_zistepv + du * d_zistepu;
		s32 izi = (s32)(zi * 0x8000 * 0x10000);
		if((intptr_t) pdest & 0x02){
			if(cutoutbuf[(pdest-d_pzbuffer)] == 1 &&
					*pdest < (s16)(izi >> 16))
				*pdest++ = (s16)(izi >> 16);
			izi += izistep;
			count--;
		}
		s32 doublecount = count >> 1;
		if(doublecount > 0){
			do{
				u32 ltemp = izi >> 16;
				izi += izistep;
				ltemp |= izi & 0xFFFF0000;
				izi += izistep;
				if(cutoutbuf[(pdest-d_pzbuffer)] == 1 &&
						*pdest < (s16)ltemp)
					*(s32 *)pdest = ltemp;
				pdest += 2;
			}while(--doublecount > 0);
		}
		if(count & 1){
			if(cutoutbuf[(pdest-d_pzbuffer)] == 1 &&
					*pdest < (s16)(izi >> 16))
				*pdest = (s16)(izi >> 16);
		}
	}while((pspan = pspan->pnext) != NULL);
}
