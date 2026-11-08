// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// d_scan.c: Portable C scan-level rasterization code
#include "quakedef.h"

typedef struct {
	u8 *pdest;
	s32 *zp;
	s32 iz;
	s32 s, t;
	s32 sstep, tstep;
	s32 ls, lt; // lit water: lightmap cache coords
	s32 lsstep, ltstep;
	s32 count;
	s32 x, y; // screen position, for the dither pattern
} turbchunk_t;
typedef struct {
	s32 shift;
	s32 fine; // extra fraction bits in the steps, see D_DrawSpansImpl
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
	u32 s, t;
	u32 sstep, tstep;
	s32 s_max, t_max;
	s32 x, y;
	s32 count;
} texcursor_t;

Q_SIL void D_TurbChunk(turbchunk_t *c, s32 trans, s32 filtered, s32 lit,s32 hl);
typedef void(*turbspanfunc_t)(turbchunk_t *);
#define DIST_LUT_SIZE 1024
#define COS_LUT_SIZE 2048
enum { TURBTRANS_OPAQUE, TURBTRANS_MIXED, TURBTRANS_DITHERED };
enum { TURBLIT_NONE, TURBLIT_MONO, TURBLIT_RGB };
#define TURBFUNC_FILTERED 3
#define TURBFUNC_LIT 6
#define TURBFUNC_HL 18

#define TURBFUNC_INDEX(trans, filtered, lit, hl) \
((trans)+TURBFUNC_FILTERED*(filtered)+TURBFUNC_LIT*(lit)+TURBFUNC_HL*(hl))

#define TURB_VARIANTS_TR(X, f, l, h) \
	X(OPAQUE, f, l, h) X(MIXED, f, l, h) X(DITHERED, f, l, h)
#define TURB_VARIANTS_F(X, l, h) \
	TURB_VARIANTS_TR(X, 0, l, h) TURB_VARIANTS_TR(X, 1, l, h)
#define TURB_VARIANTS_L(X, h) TURB_VARIANTS_F(X, 0, h) \
	TURB_VARIANTS_F(X, 1, h) TURB_VARIANTS_F(X, 2, h)
#define TURB_VARIANTS(X) TURB_VARIANTS_L(X, 0) TURB_VARIANTS_L(X, 1)

#define TURB_DEFINE_CHUNKFUNC(tr, filtered, lit, hl) \
	Q_SIL void D_TurbChunk_##tr##_##filtered##_##lit##_##hl(turbchunk_t *c)\
{ D_TurbChunk(c, TURBTRANS_##tr, filtered, lit, hl); }

TURB_VARIANTS(TURB_DEFINE_CHUNKFUNC)

#define TURB_TABLE_ENTRY(tr, filtered, lit, hl) \
		[TURBFUNC_INDEX(TURBTRANS_##tr, filtered, lit, hl)] = \
	D_TurbChunk_##tr##_##filtered##_##lit##_##hl,

// The state of the span being drawn, every thread that draws has its own
static Q_TLS u8 *r_turb_pbase;
static Q_TLS turblight_t turblight;
static Q_TLS s32 *r_turb_turb;
static Q_TLS s32 *pz; // Manoel Kasimier - translucent water
static Q_TLS s32 izi, izistep;
static Q_TLS f32 turb_opacity;
static Q_TLS const u8 *d_skyfog; // the skybox texels with the sky fog applied
static u8 *cutoutbuf; // pixels a cutout span gave a texel, one per pixel
static f32 dist_lut[DIST_LUT_SIZE];
static f32 cos_lut[COS_LUT_SIZE];
static f32 dist_lut_max = 8.0f; // max dx*dx + dy*dy
static f32 cos_lut_period = 2.0f * M_PI;

static u8 *warp_rowptr[MAXHEIGHT + (TURB_AMP2 * 2)];
static s32 warp_column[MAXWIDTH + (TURB_AMP2 * 2)];
static s32 *warp_turb;

static void D_WarpRows(s32 v0, s32 v1) // rows of scr_vrect
{
	s32 *turb = warp_turb;
	u8 *dest = vid.buffer + (scr_vrect.y + v0) * vid.width + scr_vrect.x;
	for(s32 v = v0; v < v1; v++, dest += vid.width){
		s32 *col = &warp_column[turb[v]];
		u8 **row = &warp_rowptr[v];
		s32 u = 0;
		for(; u < (scr_vrect.width & ~3); u += 4){
			dest[u + 0] = row[turb[u + 0]][col[u + 0]];
			dest[u + 1] = row[turb[u + 1]][col[u + 1]];
			dest[u + 2] = row[turb[u + 2]][col[u + 2]];
			dest[u + 3] = row[turb[u + 3]][col[u + 3]];
		}
		for(; u < scr_vrect.width; u++)
			dest[u] = row[turb[u]][col[u]];
	}
}

static void D_WarpCmd(SDL_UNUSED void *arg)
{
	s32 h = scr_vrect.height;
	D_WarpRows(CLAMP(0, d_bandy0 - scr_vrect.y, h),
			CLAMP(0, d_bandy1 - scr_vrect.y, h));
}

void D_WarpScreen() // this performs a slight compression of the screen at the
{ // same time as the sine warp, to keep the edges from wrapping
	s32 w = r_refdef.vrect.width;
	s32 h = r_refdef.vrect.height;
	f32 wratio = w / (f32)scr_vrect.width;
	f32 hratio = h / (f32)scr_vrect.height;
	if(r_threaded) // the rows of every band are read
		R_ThreadsHelp();
	for(s32 v = 0; v < scr_vrect.height + TURB_AMP2 * 2; v++)
		warp_rowptr[v] = d_viewbuffer + (r_refdef.vrect.y*screenwidth)
			+ screenwidth * (s32)((f32)v * hratio * h
					/ (h + TURB_AMP2 * 2));
	for(s32 u = 0; u < scr_vrect.width + TURB_AMP2 * 2; u++)
		warp_column[u] = r_refdef.vrect.x + (s32)((f32)u * wratio * w /
				(w + TURB_AMP2 * 2));
	warp_turb = intsintable + ((s32)(cl.time*TURB_SPEED) & (TURB_CYCLE-1));
	if(!r_threaded){
		D_WarpRows(0, scr_vrect.height);
		return;
	}
	for(s32 b = 0; b < r_nbands; b++){
		R_BandAdd(b, D_WarpCmd, NULL);
		R_BandPublish(b);
	}
}

void D_AllocCutoutBuf(s32 size)
{ // calloc leaves big buffers untouched until a cutout is drawn
	free(cutoutbuf);
	if(!(cutoutbuf = calloc(size, 1)))
		Sys_Error("Not enough memory for video mode");
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
#if defined(__x86_64__) || defined(__aarch64__) // one rounding instruction
	f32 iu = floorf(*u);
	f32 iv = floorf(*v);
#else // per pixel, Q_floor beats a slow libm call (x87 on i386)
	f32 iu = Q_floor(*u);
	f32 iv = Q_floor(*v);
#endif
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

static inline s32 D_StepDiv(s32 d, s32 n)
{ // d / n for a span's last run, reciprocal_table_16 in d_draw16.s
	static const u32 r[] = { 0x40000000, 0x2aaaaaaa, 0x20000000,
		0x19999999, 0x15555555, 0x12492492, 0x10000000, 0xe38e38e,
		0xccccccc, 0xba2e8ba, 0xaaaaaaa, 0x9d89d89, 0x9249249,
		0x8888888 };
	return n < 2 ? d : (s32)(((s64)d * r[n - 2]) >> 31);
}

static inline s32 D_ClampCoord(s64 v, s32 lo, s32 hi)
{ if(v > hi) return hi; if(v < lo) return lo; return v; }

Q_SIL void D_SpanWalkBegin(spanwalk_t *w, const espan_t *pspan, s32 shift,
		s32 fine)
{
	f32 du = (f32)pspan->u;
	f32 dv = (f32)pspan->v;
	w->shift = shift;
	w->fine = fine;
	w->count = pspan->count;
	w->len = 0;
	w->sdivz = (du*d_sdivzstepu + dv*d_sdivzstepv) + d_sdivzorigin;
	w->tdivz = (dv*d_tdivzstepv + du*d_tdivzstepu) + d_tdivzorigin;
	w->zi = (dv*d_zistepv + du*d_zistepu) + d_ziorigin;
	w->sdivz_chunk = d_sdivzstepu * (1 << shift);
	w->tdivz_chunk = d_tdivzstepu * (1 << shift);
	w->zi_chunk = d_zistepu * (1 << shift);
	f32 z = (f32)0x10000 / w->zi;
	w->s = D_ClampCoord((s64)(w->sdivz * z) + sadjust, 0, bbextents);
	w->t = D_ClampCoord((s64)(w->tdivz * z) + tadjust, 0, bbextentt);
	w->snext = w->s;
	w->tnext = w->t;
	w->sstep = 0;
	w->tstep = 0;
}

Q_SIL void D_SpanWalkEvalEnd(spanwalk_t *w)
{
	s32 lo = w->fine ? 4096 : 1 << w->shift;
	f32 z = (f32)0x10000 / w->zi;
	w->snext = D_ClampCoord((s64)(w->sdivz * z) + sadjust, lo, bbextents);
	w->tnext = D_ClampCoord((s64)(w->tdivz * z) + tadjust, lo, bbextentt);
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
		w->sstep = (w->snext - w->s) >> (w->shift - w->fine);
		w->tstep = (w->tnext - w->t) >> (w->shift - w->fine);
	}else{
		f32 lastpixel = (f32)(w->len - 1);
		w->sdivz += d_sdivzstepu * lastpixel;
		w->tdivz += d_tdivzstepu * lastpixel;
		w->zi += d_zistepu * lastpixel;
		D_SpanWalkEvalEnd(w);
		if(w->len > 1){
			s32 ds = w->snext - w->s, dt = w->tnext - w->t;
			s32 n = w->len - 1;
			w->sstep = (u32)(w->fine ? D_StepDiv(ds, n) : ds / n)
				<< w->fine;
			w->tstep = (u32)(w->fine ? D_StepDiv(dt, n) : dt / n)
				<< w->fine;
		}
	}
}

Q_SIL void D_SpanWalkEndChunk(spanwalk_t *w)
{ w->s = w->snext; w->t = w->tnext; }

Q_SIL void D_SpanDepthSetup()
{ izistep = (s32)(d_zistepu * 0x8000 * 0x10000); }

Q_SIL s32 *D_SpanDepthPtr(const espan_t *pspan)
{ return d_pzbuffer + (d_zwidth * pspan->v) + pspan->u; }

Q_SIL s32 D_SpanDepthStart(f32 zi)
{ return (s32)(zi * 0x8000 * 0x10000); }

Q_SIL void D_SpanDepthBegin(const espan_t *pspan, f32 zi)
{ pz = D_SpanDepthPtr(pspan); izi = D_SpanDepthStart(zi); }

Q_SIL s32 D_SpanDepthPass()
{ return *pz <= izi; }

Q_SIL void D_SpanDepthStep()
{ izi += izistep; pz++; }

Q_SIL s32 D_DitherIndex(s32 x, s32 y)
{ return(x & 1) + ((y & 1) << 1); }

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
	s32 *zp = c->zp;
	s32 iz = c->iz;
	s32 s = c->s;
	s32 t = c->t;
	s32 count = c->count;
	s32 x = c->x;
	s32 y = c->y;
	s32 ls = c->ls;
	s32 lt = c->lt;
	const s32 sstep = c->sstep;
	const s32 tstep = c->tstep;
	const s32 lsstep = c->lsstep;
	const s32 ltstep = c->ltstep;
	const s32 izstep = izistep;
	const u8 *pbase = r_turb_pbase;
	const s32 *turb = r_turb_turb;
	const u8 *lmbase = turblight.base;
	const s32 lmwidth = turblight.width;
	const u8 *cmap = CURWORLDCMAP;
	const s32 *chr = lit_chan[0], *chg = lit_chan[1], *chb = lit_chan[2];
	const u8 *lut = &lit_lut[0][0][0];
	const f32 time = cl.time;
	const f32 ditherthresh = 1 - turb_opacity;
	const s32 foglut = (trans == TURBTRANS_MIXED) ?
		(s32)(turb_opacity*FOG_LUT_LEVELS) : 0;
	do{
		s32 draw = 1;
		if(trans != TURBTRANS_OPAQUE){
			draw = *zp <= iz;
			if(draw && trans == TURBTRANS_DITHERED)
				draw = D_Dither(pdest, ditherthresh);
		}
		if(draw){
			u8 pix = D_TurbSample(pbase, turb, time, s, t, x, y,
					filtered, hl);
			s32 li = (lt >> 16) * lmwidth + (ls >> 16);
			if(lit == TURBLIT_MONO)
				pix = cmap[(lmbase[li] << 8) + pix];
			else if(lit == TURBLIT_RGB){
				const u8 *l = lmbase + li * 4;
				u8 ir = cmap[(l[0] << 8) + pix];
				u8 ig = cmap[(l[1] << 8) + pix];
				u8 ib = cmap[(l[2] << 8) + pix];
				pix = lut[chr[ir] + chg[ig] + chb[ib]];
			}
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
		if(lit){
			ls += lsstep;
			lt += ltstep;
		}
		x++;
	}while(--count > 0);
	c->pdest = pdest;
	c->zp = zp;
	c->iz = iz;
}

static s32 D_TurbFuncIndex(s32 lit)
{
	s32 trans = TURBTRANS_OPAQUE;
	if(r_alphapass)
		trans = r_alphastyle.value ? TURBTRANS_DITHERED:TURBTRANS_MIXED;
	return TURBFUNC_INDEX(trans, r_dithertex.value != 0, lit,
			r_hlwater.value != 0);
}

static const turbspanfunc_t turbdrawfunc[36] = {
	TURB_VARIANTS(TURB_TABLE_ENTRY)
};

// What Turbulent needs besides the gradients, the tables it uses built.
// ts->light comes first, from D_TurbSetLight.
void D_TurbSetup(turbstate_t *ts, f32 opacity)
{
	if(r_hlwater.value)
		InitHLWaterLUT();
	if(ts->light.mode == TURBLIT_RGB){
		R_BuildLitLUT();
		if(lit_chan_pal != CURWORLDPAL) R_BuildLitChan();
	}
	if(r_alphastyle.value == 0 && opacity != 1 && opacity)
		R_BuildColorMixLUT(0);
	ts->opacity = opacity;
	ts->func = D_TurbFuncIndex(ts->light.mode);
	ts->turb = sintable + ((s32)(cl.time*TURB_SPEED) & (TURB_CYCLE-1));
	ts->pbase = (u8 *)cacheblock;
}

Q_SIL void D_TurbSetChunk(turbchunk_t *c, const spanwalk_t *w)
{
	c->s = w->s & ((TURB_CYCLE << 16) - 1);
	c->t = w->t & ((TURB_CYCLE << 16) - 1);
	c->sstep = w->sstep;
	c->tstep = w->tstep;
	c->count = w->len;
}

void D_TurbSetLight(turblight_t *l, const surfcache_t *cache,
		const msurface_t *pface, s32 mip)
{ // light rows built by R_DrawSurfaceLight, NULL for unlit water
	l->mode = TURBLIT_NONE;
	if(!cache)
		return;
	l->mode = R_LitWaterRGB() ? TURBLIT_RGB : TURBLIT_MONO;
	l->base = cache->data;
	l->width = cache->width;
	l->mip = mip;
	// tiled gradients offset s and t by 0x20000000 plus the texturemins
	// remainder that keeps them aligned to the tile (D_CalcGradients)
	l->soff = 0x20000000
		+ ((pface->texturemins[0] & (TURB_CYCLE-1)) << 16);
	l->toff = 0x20000000
		+ ((pface->texturemins[1] & (TURB_CYCLE-1)) << 16);
	l->smax = ((pface->extents[0] << 16) >> mip) - 1;
	l->tmax = ((pface->extents[1] << 16) >> mip) - 1;
}

Q_SIL void D_TurbSetLightChunk(turbchunk_t *c, const spanwalk_t *w)
{ // lightmap coords follow the unwarped texture coords, clamped like them
	s32 lo = 1 << w->shift;
	s32 mip = turblight.mip;
	s32 smax = turblight.smax;
	s32 tmax = turblight.tmax;
	c->ls = D_ClampCoord((w->s - turblight.soff) >> mip, 0, smax);
	c->lt = D_ClampCoord((w->t - turblight.toff) >> mip, 0, tmax);
	s32 snext = D_ClampCoord((w->snext - turblight.soff) >> mip, lo, smax);
	s32 tnext = D_ClampCoord((w->tnext - turblight.toff) >> mip, lo, tmax);
	c->lsstep = 0;
	c->ltstep = 0;
	if(w->count){
		c->lsstep = (snext - c->ls) >> w->shift;
		c->ltstep = (tnext - c->lt) >> w->shift;
	}else if(w->len > 1){
		c->lsstep = (snext - c->ls) / (w->len - 1);
		c->ltstep = (tnext - c->lt) / (w->len - 1);
	}
}

static void D_TurbSpan(const espan_t *pspan, turbspanfunc_t drawspan)
{
	spanwalk_t w;
	turbchunk_t c;
	c.pdest = (u8*)((u8*)d_viewbuffer + (screenwidth*pspan->v) + pspan->u);
	D_SpanWalkBegin(&w, pspan, 4, 0); // exact coords every 16px
	c.zp = D_SpanDepthPtr(pspan);
	c.iz = D_SpanDepthStart(w.zi);
	c.x = pspan->u;
	c.y = pspan->v;
	do{
		D_SpanWalkNextChunk(&w);
		D_TurbSetChunk(&c, &w);
		if(turblight.mode)
			D_TurbSetLightChunk(&c, &w);
		drawspan(&c);
		c.x += w.len;
		D_SpanWalkEndChunk(&w);
	}while(w.count > 0);
}

void Turbulent(espan_t *pspan, const turbstate_t *ts)
{
	turblight = ts->light;
	turb_opacity = ts->opacity;
	r_turb_turb = ts->turb;
	r_turb_pbase = ts->pbase;
	D_SpanDepthSetup();
	turbspanfunc_t drawspan = turbdrawfunc[ts->func];
	do{
		D_TurbSpan(pspan, drawspan);
	}while((pspan = pspan->pnext) != NULL);
}

Q_SIL void D_CursorBegin(texcursor_t *c, const spanwalk_t *w, s32 dithered,
		s32 x, s32 y)
{
	c->pbase = (u8 *)cacheblock;
	c->width = cachewidth;
	c->s = (u32)w->s << w->fine;
	c->t = (u32)w->t << w->fine;
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
		c->x = x;
		c->y = y;
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

// Dithered texels of a chunk that can't reach past the texture's edges need
// no clamping. A pixel and its neighbour take the two offsets of their row,
// so pairs of pixels step together, each lane with its offset added once.
Q_SIL s32 D_CursorPairs(const texcursor_t *c, s32 *sa, s32 *ta, s32 *sb,
		s32 *tb)
{
	s32 i = D_DitherIndex(c->x, c->y), j = i ^ 1;
	s32 s0 = c->s, s1 = s0 + (c->count - 1) * (s32)c->sstep;
	s32 t0 = c->t, t1 = t0 + (c->count - 1) * (s32)c->tstep;
	if(q_min(s0, s1) + q_min(dither_s[i], dither_s[j]) < 0
	|| q_min(t0, t1) + q_min(dither_t[i], dither_t[j]) < 0
	|| q_max(s0, s1) + q_max(dither_s[i], dither_s[j]) > c->s_max + 0xffff
	|| q_max(t0, t1) + q_max(dither_t[i], dither_t[j]) > c->t_max + 0xffff)
		return 0;
	*sa = s0 + dither_s[i];
	*ta = t0 + dither_t[i];
	*sb = s0 + (s32)c->sstep + dither_s[j];
	*tb = t0 + (s32)c->tstep + dither_t[j];
	return 1;
}

Q_SIL u8 *D_SpanChunkPairs(const texcursor_t *c, u8 *pdest, const u8 *fog,
		s32 sa, s32 ta, s32 sb, s32 tb)
{ // locals only, the byte stores could change anything behind a pointer
	const u8 *pbase = c->pbase;
	s32 width = c->width;
	s32 sstep = 2 * (s32)c->sstep, tstep = 2 * (s32)c->tstep;
	s32 n = c->count;
	for(; n >= 2; n -= 2, pdest += 2){
		u8 a = pbase[(sa >> 16) + (ta >> 16) * width];
		u8 b = pbase[(sb >> 16) + (tb >> 16) * width];
		pdest[0] = fog ? fog[a] : a;
		pdest[1] = fog ? fog[b] : b;
		sa += sstep; ta += tstep;
		sb += sstep; tb += tstep;
	}
	if(n){
		u8 a = pbase[(sa >> 16) + (ta >> 16) * width];
		*pdest++ = fog ? fog[a] : a;
	}
	return pdest;
}

Q_SIL u8 *D_SpanChunkNormal(texcursor_t *c, u8 *pdest, s32 dithered)
{
	s32 sa, ta, sb, tb;
	if(dithered && D_CursorPairs(c, &sa, &ta, &sb, &tb))
		return D_SpanChunkPairs(c, pdest, NULL, sa, ta, sb, tb);
	do{
		*pdest++ = dithered ? D_CursorTexel(c, 1) :
			c->pbase[(c->s >> 20) + (c->t >> 20) * c->width];
		D_CursorAdvance(c);
	}while(--c->count > 0);
	return pdest;
}

Q_SIL u8 *D_SpanChunkCutout(texcursor_t *c, u8 *pdest, s32 dithered)
{ // a local walks cutoutbuf, the byte stores could change the globals
	u8 *cut = cutoutbuf + (pdest - d_viewbuffer);
	do{
		if(D_SpanDepthPass()){
			u8 pix = D_CursorTexel(c, dithered);
			*cut = pix != 0xff;
			if(pix != 0xff) *pdest = pix;
		}
		pdest++;
		cut++;
		D_SpanDepthStep();
		D_CursorAdvance(c);
	}while(--c->count > 0);
	return pdest;
}

// The skybox texels with the sky fog applied, NULL without fog
void D_SkyFogSetup()
{
	static u8 tab[256];
	static s64 frame = -1;
	static s32 pal = -1, lut = -1;
	s32 foglut = CLAMP(0, map_skyfog*FOG_LUT_LEVELS, FOG_LUT_LEVELS-1);
	d_skyfog = NULL;
	if(!(fog_density > 0))
		return;
	if(frame != r_framecount || pal != fog_pal_index || lut != foglut){
		for(s32 i = 0; i < 256; i++)
			tab[i] = color_mix_lut[i][fog_pal_index][foglut];
		frame = r_framecount;
		pal = fog_pal_index;
		lut = foglut;
	}
	d_skyfog = tab;
}

Q_SIL u8 *D_SpanChunkSkybox(texcursor_t *c, u8 *pdest, s32 dithered)
{
	const u8 *fog = d_skyfog;
	s32 sa, ta, sb, tb;
	if(dithered && D_CursorPairs(c, &sa, &ta, &sb, &tb))
		return fog ? D_SpanChunkPairs(c, pdest, fog, sa, ta, sb, tb)
			: D_SpanChunkPairs(c, pdest, NULL, sa, ta, sb, tb);
	if(fog){
		do{
			*pdest++ = fog[D_CursorTexel(c, dithered)];
			D_CursorAdvance(c);
		}while(--c->count > 0);
	}else{
		do{
			*pdest++ = D_CursorTexel(c, dithered);
			D_CursorAdvance(c);
		}while(--c->count > 0);
	}
	return pdest;
}

Q_SIL u8 *D_SpanChunkBlend(texcursor_t *c, u8 *pdest, f32 opacity,
		s32 dithered)
{
	s32 foglut = opacity*FOG_LUT_LEVELS; // D_SetupCached built the table
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
		f32 opacity, s32 dithered, s32 x, s32 y)
{ // x and y are the chunk's screen position, for the dither pattern
	texcursor_t c;
	D_CursorBegin(&c, w, dithered, x, y);
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

static u8 *D_SpanDest(const espan_t *pspan)
{ return (u8 *)d_viewbuffer + (screenwidth * pspan->v) + pspan->u; }

Q_SIL void D_DrawSpansImpl(espan_t *pspan, s32 type, f32 opacity, s32 dithered)
{
	s32 depthtested = (type == SPAN_TRANS || type == SPAN_CUTOUT);
	// d_draw16.s keeps the bits that shifting a 16 pixel run's step drops.
	// 12.20 has room for them, lightmapped surfaces end at 2000 texels
	s32 fine = type == SPAN_NORMAL && !dithered ? 4 : 0;
	D_SpanDepthSetup();
	do{
		spanwalk_t w;
		u8 *pdest = D_SpanDest(pspan);
		s32 x = pspan->u;
		D_SpanWalkBegin(&w, pspan, 4, fine); // exact coords every 16px
		if(depthtested)
			D_SpanDepthBegin(pspan, w.zi);
		do{
			D_SpanWalkNextChunk(&w);
			pdest = D_DrawSpanChunk(&w, pdest, type, opacity,
					dithered, x, pspan->v);
			x += w.len;
			D_SpanWalkEndChunk(&w);
		}while(w.count > 0);
	}while((pspan = pspan->pnext) != NULL);
}

void D_DrawSpans(espan_t *pspan, s32 type, f32 opacity)
{ // a loop of its own for each type
	switch(type){
	case SPAN_NORMAL: D_DrawSpansImpl(pspan, SPAN_NORMAL, 0, 0); break;
	case SPAN_CUTOUT: D_DrawSpansImpl(pspan, SPAN_CUTOUT, 0, 0); break;
	case SPAN_SKYBOX: D_DrawSpansImpl(pspan, SPAN_SKYBOX, 0, 0); break;
	default: D_DrawSpansImpl(pspan, SPAN_TRANS, opacity, 0); break;
	}
}

// The span drawing state set up for a surface, kept for later or elsewhere
void D_SaveGradients(dgrad_t *g)
{
	g->sdivzstepu = d_sdivzstepu; g->tdivzstepu = d_tdivzstepu;
	g->zistepu = d_zistepu;
	g->sdivzstepv = d_sdivzstepv; g->tdivzstepv = d_tdivzstepv;
	g->zistepv = d_zistepv;
	g->sdivzorigin = d_sdivzorigin; g->tdivzorigin = d_tdivzorigin;
	g->ziorigin = d_ziorigin;
	g->sadjust = sadjust; g->tadjust = tadjust;
	g->bbextents = bbextents; g->bbextentt = bbextentt;
	g->cacheblock = cacheblock;
	g->cachewidth = cachewidth; g->cacheheight = cacheheight;
	g->skyfog = d_skyfog;
}

void D_LoadGradients(const dgrad_t *g)
{
	d_sdivzstepu = g->sdivzstepu; d_tdivzstepu = g->tdivzstepu;
	d_zistepu = g->zistepu;
	d_sdivzstepv = g->sdivzstepv; d_tdivzstepv = g->tdivzstepv;
	d_zistepv = g->zistepv;
	d_sdivzorigin = g->sdivzorigin; d_tdivzorigin = g->tdivzorigin;
	d_ziorigin = g->ziorigin;
	sadjust = g->sadjust; tadjust = g->tadjust;
	bbextents = g->bbextents; bbextentt = g->bbextentt;
	cacheblock = g->cacheblock;
	cachewidth = g->cachewidth; cacheheight = g->cacheheight;
	d_skyfog = g->skyfog;
}

// Inlined into D_DrawSurfaces, the dithered loops changed how GCC allocated
// registers for the plain one there (+1.8% instructions)
Q_NOINLINE void D_DrawSpansDithered(espan_t *pspan, s32 type, f32 opacity)
{ D_DrawSpansImpl(pspan, type, opacity, 1); }

void D_DrawZSpans(espan_t *pspan)
{
	s32 izistep = (s32)(d_zistepu * 0x8000 * 0x10000);
	do{
		s32 *pdest = d_pzbuffer + (d_zwidth * pspan->v) + pspan->u;
		s32 count = pspan->count; // local: s32 stores may alias it
		f32 du = (f32)pspan->u; // calculate the initial 1/z
		f32 dv = (f32)pspan->v;
		f64 zi = d_ziorigin + dv * d_zistepv + du * d_zistepu;
		s32 izi = (s32)(zi * 0x8000 * 0x10000);
		for(s32 i = 0; i < count; i++, izi += izistep)
			pdest[i] = izi;
	}while((pspan = pspan->pnext) != NULL);
}

void D_DrawZSpansTrans(espan_t *pspan)
{ // holes stay open, so every pixel is tested on its own
	s32 izistep = (s32)(d_zistepu * 0x8000 * 0x10000);
	do{
		s32 *pdest = d_pzbuffer + (d_zwidth * pspan->v) + pspan->u;
		s32 count = pspan->count;
		f32 du = (f32)pspan->u; // calculate the initial 1/z
		f32 dv = (f32)pspan->v;
		f64 zi = d_ziorigin + dv * d_zistepv + du * d_zistepu;
		s32 izi = (s32)(zi * 0x8000 * 0x10000);
		for(s32 i = 0; i < count; i++, pdest++, izi += izistep)
			if(cutoutbuf[pdest-d_pzbuffer] == 1 && *pdest < izi)
				*pdest = izi;
	}while((pspan = pspan->pnext) != NULL);
}
