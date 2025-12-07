// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// d_scan.c: Portable C scan-level rasterization code
// Water can be: 
// 	unlit/lit
// 	unfiltered/filtered (dithered texture)
// 	opaque/mixed/dithered (transparency)
// 	classic/HL-style
// Any one of those can be toggled, and all combinations have separate drawing
// functions.
#include "quakedef.h"

#define DIST_LUT_SIZE 1024
#define COS_LUT_SIZE 2048

static void D_DrawTurbulentSpan();
static void D_DrawTurbulentSpanMixed();
static void D_DrawTurbulentSpanDithered();
static void D_DrawTurbulentSpanFiltered();
static void D_DrawTurbulentSpanFilteredMixed();
static void D_DrawTurbulentSpanFilteredDithered();
static void D_DrawTurbulentSpanLit();
static void D_DrawTurbulentSpanLitMixed();
static void D_DrawTurbulentSpanLitDithered();
static void D_DrawTurbulentSpanLitFiltered();
static void D_DrawTurbulentSpanLitFilteredMixed();
static void D_DrawTurbulentSpanLitFilteredDithered();
static void D_DrawTurbulentSpanHL();
static void D_DrawTurbulentSpanMixedHL();
static void D_DrawTurbulentSpanDitheredHL();
static void D_DrawTurbulentSpanFilteredHL();
static void D_DrawTurbulentSpanFilteredMixedHL();
static void D_DrawTurbulentSpanFilteredDitheredHL();
static void D_DrawTurbulentSpanLitHL();
static void D_DrawTurbulentSpanLitMixedHL();
static void D_DrawTurbulentSpanLitDitheredHL();
static void D_DrawTurbulentSpanLitFilteredHL();
static void D_DrawTurbulentSpanLitFilteredMixedHL();
static void D_DrawTurbulentSpanLitFilteredDitheredHL();

static void (*turbdrawfunc[24])() = {
D_DrawTurbulentSpan,
D_DrawTurbulentSpanMixed,
D_DrawTurbulentSpanDithered,
D_DrawTurbulentSpanFiltered,
D_DrawTurbulentSpanFilteredMixed,
D_DrawTurbulentSpanFilteredDithered,
D_DrawTurbulentSpanLit,
D_DrawTurbulentSpanLitMixed,
D_DrawTurbulentSpanLitDithered,
D_DrawTurbulentSpanLitFiltered,
D_DrawTurbulentSpanLitFilteredMixed,
D_DrawTurbulentSpanLitFilteredDithered,
D_DrawTurbulentSpanHL,
D_DrawTurbulentSpanMixedHL,
D_DrawTurbulentSpanDitheredHL,
D_DrawTurbulentSpanFilteredHL,
D_DrawTurbulentSpanFilteredMixedHL,
D_DrawTurbulentSpanFilteredDitheredHL,
D_DrawTurbulentSpanLitHL,
D_DrawTurbulentSpanLitMixedHL,
D_DrawTurbulentSpanLitDitheredHL,
D_DrawTurbulentSpanLitFilteredHL,
D_DrawTurbulentSpanLitFilteredMixedHL,
D_DrawTurbulentSpanLitFilteredDitheredHL,
};

static u8 *r_turb_pbase, *r_turb_pdest;
static s32 r_turb_s, r_turb_t, r_turb_sstep, r_turb_tstep;
static s32 *r_turb_turb;
static s32 r_turb_spancount;
static s16 *pz; // Manoel Kasimier - translucent water
static s32 izi, izistep;
static u8 cutoutbuf[MAXHEIGHT*MAXWIDTH];
static f32 dist_lut[DIST_LUT_SIZE];
static f32 cos_lut[COS_LUT_SIZE];
static f32 dist_lut_max = 8.0f; // max dx*dx + dy*dy
static f32 cos_lut_period = 2.0f * M_PI;
static f32 HL_TextureToWaveScale = 1.5f;
static f32 HL_RippleScale = 1.5f;
static s32 HLWaterLutInited = 0;
static f32 turb_opacity = 0;
static s32 flt_y = 0;
static s32 flt_x = 0;
static s32 flt_start_x = 0;
static s32 flt_cur_x = 0;

void D_WarpScreen() // this performs a slight compression of the screen at the
{ // same time as the sine warp, to keep the edges from wrapping
	s32 w = r_refdef.vrect.width;
	s32 h = r_refdef.vrect.height;
	f32 wratio = w / (f32)scr_vrect.width;
	f32 hratio = h / (f32)scr_vrect.height;
	u8 *rowptr[MAXHEIGHT + (AMP2 * 2)];
	s32 column[MAXWIDTH + (AMP2 * 2)];
	for (s32 v = 0; v < scr_vrect.height + AMP2 * 2; v++)
		rowptr[v] = d_viewbuffer + (r_refdef.vrect.y * screenwidth) +
		    screenwidth * (s32)((f32)v * hratio * h / (h + AMP2 * 2));
	for (s32 u = 0; u < scr_vrect.width + AMP2 * 2; u++)
		column[u] = r_refdef.vrect.x + (s32)((f32)u * wratio * w /
						     (w + AMP2 * 2));
	s32 *turb = intsintable + ((s32)(cl.time * SPEED) & (CYCLE - 1));
	u8 *dest = vid.buffer + scr_vrect.y * vid.width + scr_vrect.x;
	for (s32 v = 0; v < scr_vrect.height; v++, dest += vid.width) {
		s32 *col = &column[turb[v]];
		u8 **row = &rowptr[v];
		for (s32 u = 0; u < scr_vrect.width; u += 4) {
			dest[u + 0] = row[turb[u + 0]][col[u + 0]];
			dest[u + 1] = row[turb[u + 1]][col[u + 1]];
			dest[u + 2] = row[turb[u + 2]][col[u + 2]];
			dest[u + 3] = row[turb[u + 3]][col[u + 3]];
		}
	}
}

static inline f32 fast_sqrt_dist(f32 x)
{ // x in [0, 8]
    f32 f = x * (f32)(DIST_LUT_SIZE - 1) / dist_lut_max;
    s32 i = (s32)f;
    if (i < 0) i = 0;
    if (i >= DIST_LUT_SIZE - 1) i = DIST_LUT_SIZE - 1;
    return dist_lut[i];
}

static inline f32 fast_cos(f32 theta)
{ // wrap manually
    f32 n = theta * (f32)COS_LUT_SIZE / cos_lut_period;
    s32 i = (s32)n;
    i &= (COS_LUT_SIZE - 1); // requires power of two table
    return cos_lut[i];
}

static void InitHLWaterLUT()
{
	if(HLWaterLutInited) return;
	for (s32 i = 0; i < DIST_LUT_SIZE; i++) { // sqrt LUT for x in [0, 8]
		f32 x = (f32)i * (dist_lut_max / (DIST_LUT_SIZE - 1));
		dist_lut[i] = sqrtf(x);
	}
	for (s32 i = 0; i < COS_LUT_SIZE; i++) { // cos LUT for θ in [0, 2π]
		f32 t = (f32)i * (cos_lut_period / COS_LUT_SIZE);
		cos_lut[i] = cosf(t);
	}
	HLWaterLutInited = 1;
}

static void HLWarp(f32 *u, f32 *v, f32 time)
{ // classic warp cos(y), cos(x)
	f32 x = *u;
	f32 y = *v;
	f32 waveSize = (1.0f / 8.0f) * HL_TextureToWaveScale;
	f32 waveStrength = 0.01f * HL_RippleScale;
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
{ // Ripple helper
	f32 dx = u - px; // [-2, 2]
	f32 dy = v - py; // [-2, 2]
	f32 r2 = dx*dx + dy*dy; // [0, 8]
	f32 dist = fast_sqrt_dist(r2);
	f32 f = fast_sqrt_dist(dist) * freq - t;
	f32 wave = fast_cos(f);
	f32 decay = 1.0f - dist / size;
	if (decay < 0.0f)
		decay = 0.0f;
	wave *= decay;
	*outx += dx * wave * strength;
	*outy += dy * wave * strength;
}



static void HLRipple(f32 *u, f32 *v, f32 time)
{ // HL three-point ripple, tiled
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
	f32 strength = 0.02f * HL_RippleScale;
	f32 t = time * 20.0f;
	for (s32 i = 0; i < 3; i++) // main tile
		HLSingleRipple(&ox, &oy, baseu, basev, iu + px[i], iv + py[i],
				size, freq, strength, t);
	for (s32 ti = 0; ti < 8; ti += 2) // surrounding tiles
	for (s32 i = 0; i < 3; i++)
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

static void D_DrawTurbulentSpan()
{
	do {
		s32 s=((r_turb_s+r_turb_turb[(r_turb_t>>16)&(CYCLE-1)])>>16)&63;
		s32 t=((r_turb_t+r_turb_turb[(r_turb_s>>16)&(CYCLE-1)])>>16)&63;
		s32 pix = *(r_turb_pbase + (t << 6) + s);
		*r_turb_pdest++ = pix;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanLit()
{
	do {
		s32 s=((r_turb_s+r_turb_turb[(r_turb_t>>16)&(CYCLE-1)])>>16)&63;
		s32 t=((r_turb_t+r_turb_turb[(r_turb_s>>16)&(CYCLE-1)])>>16)&63;
		s32 pix = *(r_turb_pbase + (t << 6) + s);
		s32 lit = *(litwater_base+(r_turb_pdest-d_viewbuffer));
		*r_turb_pdest++ = D_TurbMixLit(pix, lit);
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanFiltered()
{
	do {
		s32 dither_idx = (flt_cur_x & 1) + ((flt_y & 1) << 1);
		s32 s_d = r_turb_s + dither_s[dither_idx];
		s32 t_d = r_turb_t + dither_t[dither_idx];
		s32 s = ((s_d+r_turb_turb[(t_d>>16)&(CYCLE-1)])>>16)&63;
		s32 t = ((t_d+r_turb_turb[(s_d>>16)&(CYCLE-1)])>>16)&63;
		s32 pix = *(r_turb_pbase + (t << 6) + s);
		*r_turb_pdest++ = pix;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
		flt_cur_x++;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanLitFiltered()
{
	do {
		s32 dither_idx = (flt_cur_x & 1) + ((flt_y & 1) << 1);
		s32 s_d = r_turb_s + dither_s[dither_idx];
		s32 t_d = r_turb_t + dither_t[dither_idx];
		s32 s = ((s_d+r_turb_turb[(t_d>>16)&(CYCLE-1)])>>16)&63;
		s32 t = ((t_d+r_turb_turb[(s_d>>16)&(CYCLE-1)])>>16)&63;
		s32 pix = *(r_turb_pbase + (t << 6) + s);
		s32 lit = *(litwater_base+(r_turb_pdest-d_viewbuffer));
		*r_turb_pdest++ = D_TurbMixLit(pix, lit);
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
		flt_cur_x++;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanMixed()
{
	do {
		if (*pz <= (izi >> 16)) {
		s32 s=((r_turb_s+r_turb_turb[(r_turb_t>>16)&(CYCLE-1)])>>16)&63;
		s32 t=((r_turb_t+r_turb_turb[(r_turb_s>>16)&(CYCLE-1)])>>16)&63;
			*r_turb_pdest = color_mix_lut[*(r_turb_pbase+(t<<6)+s)]
			    [*r_turb_pdest][(s32)(turb_opacity*FOG_LUT_LEVELS)];
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanLitMixed()
{
	do {
		if (*pz <= (izi >> 16)) {
		s32 s=((r_turb_s+r_turb_turb[(r_turb_t>>16)&(CYCLE-1)])>>16)&63;
		s32 t=((r_turb_t+r_turb_turb[(r_turb_s>>16)&(CYCLE-1)])>>16)&63;
			s32 pix = *(r_turb_pbase + (t << 6) + s);
			s32 lit = *(litwater_base+(r_turb_pdest-d_viewbuffer));
			pix = D_TurbMixLit(pix, lit);
			*r_turb_pdest = color_mix_lut[pix][*r_turb_pdest]
					[(s32)(turb_opacity*FOG_LUT_LEVELS)];
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanFilteredMixed()
{
	do {
		if (*pz <= (izi >> 16)) {
			s32 dither_idx = (flt_cur_x & 1) + ((flt_y & 1) << 1);
			s32 s_d = r_turb_s + dither_s[dither_idx];
			s32 t_d = r_turb_t + dither_t[dither_idx];
			s32 s = ((s_d+r_turb_turb[(t_d>>16)&(CYCLE-1)])>>16)&63;
			s32 t = ((t_d+r_turb_turb[(s_d>>16)&(CYCLE-1)])>>16)&63;
			*r_turb_pdest = color_mix_lut[*(r_turb_pbase+(t<<6)+s)]
			[*r_turb_pdest][(s32)(turb_opacity*FOG_LUT_LEVELS)];
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
		flt_cur_x++;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanLitFilteredMixed()
{
	do {
		if (*pz <= (izi >> 16)) {
			s32 dither_idx = (flt_cur_x & 1) + ((flt_y & 1) << 1);
			s32 s_d = r_turb_s + dither_s[dither_idx];
			s32 t_d = r_turb_t + dither_t[dither_idx];
			s32 s = ((s_d+r_turb_turb[(t_d>>16)&(CYCLE-1)])>>16)&63;
			s32 t = ((t_d+r_turb_turb[(s_d>>16)&(CYCLE-1)])>>16)&63;
			s32 pix = *(r_turb_pbase + (t << 6) + s);
			s32 lit = *(litwater_base+(r_turb_pdest-d_viewbuffer));
			pix = D_TurbMixLit(pix, lit);
			*r_turb_pdest = color_mix_lut[pix][*r_turb_pdest]
					[(s32)(turb_opacity*FOG_LUT_LEVELS)];
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
		flt_cur_x++;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanDithered()
{
	do {
		if (*pz <= (izi>>16) && D_Dither(r_turb_pdest, 1-turb_opacity)){
		s32 s=((r_turb_s+r_turb_turb[(r_turb_t>>16)&(CYCLE-1)])>>16)&63;
		s32 t=((r_turb_t+r_turb_turb[(r_turb_s>>16)&(CYCLE-1)])>>16)&63;
				*r_turb_pdest = *(r_turb_pbase + (t << 6) + s);
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanLitDithered()
{
	do {
		if (*pz <= (izi>>16) && D_Dither(r_turb_pdest,1-turb_opacity)) {
		s32 s=((r_turb_s+r_turb_turb[(r_turb_t>>16)&(CYCLE-1)])>>16)&63;
		s32 t=((r_turb_t+r_turb_turb[(r_turb_s>>16)&(CYCLE-1)])>>16)&63;
			s32 pix = *(r_turb_pbase + (t << 6) + s);
			s32 lit = *(litwater_base+(r_turb_pdest-d_viewbuffer));
			*r_turb_pdest = D_TurbMixLit(pix, lit);
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanFilteredDithered()
{
	do {
		if (*pz <= (izi>>16) && D_Dither(r_turb_pdest, 1-turb_opacity)){
			s32 dither_idx = (flt_cur_x & 1) + ((flt_y & 1) << 1);
			s32 s_d = r_turb_s + dither_s[dither_idx];
			s32 t_d = r_turb_t + dither_t[dither_idx];
			s32 s = ((s_d+r_turb_turb[(t_d>>16)&(CYCLE-1)])>>16)&63;
			s32 t = ((t_d+r_turb_turb[(s_d>>16)&(CYCLE-1)])>>16)&63;
				*r_turb_pdest = *(r_turb_pbase + (t << 6) + s);
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
		flt_cur_x++;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanLitFilteredDithered()
{
	do {
		if (*pz <= (izi>>16) && D_Dither(r_turb_pdest,1-turb_opacity)) {
			s32 dither_idx = (flt_cur_x & 1) + ((flt_y & 1) << 1);
			s32 s_d = r_turb_s + dither_s[dither_idx];
			s32 t_d = r_turb_t + dither_t[dither_idx];
			s32 s = ((s_d+r_turb_turb[(t_d>>16)&(CYCLE-1)])>>16)&63;
			s32 t = ((t_d+r_turb_turb[(s_d>>16)&(CYCLE-1)])>>16)&63;
			s32 pix = *(r_turb_pbase + (t << 6) + s);
			s32 lit = *(litwater_base+(r_turb_pdest-d_viewbuffer));
			*r_turb_pdest = D_TurbMixLit(pix, lit);
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
		flt_cur_x++;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanHL()
{
	do {
		f32 u = (f32)((r_turb_s >> 16) & 63) * (1.0f/64.0f);
		f32 v = (f32)((r_turb_t >> 16) & 63) * (1.0f/64.0f);
		HLWarp(&u, &v, cl.time);
		HLRipple(&u, &v, cl.time);
		s32 su = ((s32)(u * 64.0f)) & 63;
		s32 tv = ((s32)(v * 64.0f)) & 63;
		*r_turb_pdest++ = r_turb_pbase[(tv << 6) + su];
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanLitHL()
{
	do {
		f32 u = (f32)((r_turb_s >> 16) & 63) * (1.0f/64.0f);
		f32 v = (f32)((r_turb_t >> 16) & 63) * (1.0f/64.0f);
		HLWarp(&u, &v, cl.time);
		HLRipple(&u, &v, cl.time);
		s32 su = ((s32)(u * 64.0f)) & 63;
		s32 tv = ((s32)(v * 64.0f)) & 63;
		s32 pix = r_turb_pbase[(tv << 6) + su];
		s32 lit = *(litwater_base+(r_turb_pdest-d_viewbuffer));
		*r_turb_pdest++ = D_TurbMixLit(pix, lit);
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanFilteredHL()
{
	do {
		s32 dither_idx = (flt_cur_x & 1) + ((flt_y & 1) << 1);
		s32 s_d = r_turb_s + dither_s[dither_idx];
		s32 t_d = r_turb_t + dither_t[dither_idx];
		f32 u = (f32)((s_d >> 16) & 63) * (1.0f/64.0f);
		f32 v = (f32)((t_d >> 16) & 63) * (1.0f/64.0f);
		HLWarp(&u, &v, cl.time);
		HLRipple(&u, &v, cl.time);
		s32 su = ((s32)(u * 64.0f)) & 63;
		s32 tv = ((s32)(v * 64.0f)) & 63;
		*r_turb_pdest++ = r_turb_pbase[(tv << 6) + su];
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
		flt_cur_x++;
	} while (--r_turb_spancount > 0);


}

static void D_DrawTurbulentSpanLitFilteredHL()
{
	do {
		s32 dither_idx = (flt_cur_x & 1) + ((flt_y & 1) << 1);
		s32 s_d = r_turb_s + dither_s[dither_idx];
		s32 t_d = r_turb_t + dither_t[dither_idx];
		f32 u = (f32)((s_d >> 16) & 63) * (1.0f/64.0f);
		f32 v = (f32)((t_d >> 16) & 63) * (1.0f/64.0f);
		HLWarp(&u, &v, cl.time);
		HLRipple(&u, &v, cl.time);
		s32 su = ((s32)(u * 64.0f)) & 63;
		s32 tv = ((s32)(v * 64.0f)) & 63;
		s32 pix = r_turb_pbase[(tv << 6) + su];
		s32 lit = *(litwater_base+(r_turb_pdest-d_viewbuffer));
		*r_turb_pdest++ = D_TurbMixLit(pix, lit);
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
		flt_cur_x++;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanMixedHL()
{
	do {
		if (*pz <= (izi >> 16)) {
			f32 u = (f32)((r_turb_s >> 16) & 63) * (1.0f/64.0f);
			f32 v = (f32)((r_turb_t >> 16) & 63) * (1.0f/64.0f);
			HLWarp(&u, &v, cl.time);
			HLRipple(&u, &v, cl.time);
			s32 su = ((s32)(u * 64.0f)) & 63;
			s32 tv = ((s32)(v * 64.0f)) & 63;
			*r_turb_pdest =color_mix_lut[*(r_turb_pbase+(tv<<6)+su)]
			    [*r_turb_pdest][(s32)(turb_opacity*FOG_LUT_LEVELS)];
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanLitMixedHL()
{
	do {
		if (*pz <= (izi >> 16)) {
			f32 u = (f32)((r_turb_s >> 16) & 63) * (1.0f/64.0f);
			f32 v = (f32)((r_turb_t >> 16) & 63) * (1.0f/64.0f);
			HLWarp(&u, &v, cl.time);
			HLRipple(&u, &v, cl.time);
			s32 su = ((s32)(u * 64.0f)) & 63;
			s32 tv = ((s32)(v * 64.0f)) & 63;
			s32 pix = *(r_turb_pbase + (tv << 6) + su);
			s32 lit = *(litwater_base+(r_turb_pdest-d_viewbuffer));
			pix = D_TurbMixLit(pix, lit);
			*r_turb_pdest = color_mix_lut[pix][*r_turb_pdest]
					[(s32)(turb_opacity*FOG_LUT_LEVELS)];
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanFilteredMixedHL()
{
	do {
		if (*pz <= (izi >> 16)) {
			s32 dither_idx = (flt_cur_x & 1) + ((flt_y & 1) << 1);
			s32 s_d = r_turb_s + dither_s[dither_idx];
			s32 t_d = r_turb_t + dither_t[dither_idx];
			f32 u = (f32)((s_d >> 16) & 63) * (1.0f/64.0f);
			f32 v = (f32)((t_d >> 16) & 63) * (1.0f/64.0f);
			HLWarp(&u, &v, cl.time);
			HLRipple(&u, &v, cl.time);
			s32 su = ((s32)(u * 64.0f)) & 63;
			s32 tv = ((s32)(v * 64.0f)) & 63;
			*r_turb_pdest =color_mix_lut[*(r_turb_pbase+(tv<<6)+su)]
			[*r_turb_pdest][(s32)(turb_opacity*FOG_LUT_LEVELS)];
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
		flt_cur_x++;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanLitFilteredMixedHL()
{
	do {
		if (*pz <= (izi >> 16)) {
			s32 dither_idx = (flt_cur_x & 1) + ((flt_y & 1) << 1);
			s32 s_d = r_turb_s + dither_s[dither_idx];
			s32 t_d = r_turb_t + dither_t[dither_idx];
			f32 u = (f32)((s_d >> 16) & 63) * (1.0f/64.0f);
			f32 v = (f32)((t_d >> 16) & 63) * (1.0f/64.0f);
			HLWarp(&u, &v, cl.time);
			HLRipple(&u, &v, cl.time);
			s32 su = ((s32)(u * 64.0f)) & 63;
			s32 tv = ((s32)(v * 64.0f)) & 63;
			s32 pix = r_turb_pbase[(tv << 6) + su];
			s32 lit = *(litwater_base+(r_turb_pdest-d_viewbuffer));
			pix = D_TurbMixLit(pix, lit);
			*r_turb_pdest = color_mix_lut[pix][*r_turb_pdest]
					[(s32)(turb_opacity*FOG_LUT_LEVELS)];
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
		flt_cur_x++;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanDitheredHL()
{
	do {
		if (*pz <= (izi>>16) && D_Dither(r_turb_pdest, 1-turb_opacity)){
			f32 u = (f32)((r_turb_s >> 16) & 63) * (1.0f/64.0f);
			f32 v = (f32)((r_turb_t >> 16) & 63) * (1.0f/64.0f);
			HLWarp(&u, &v, cl.time);
			HLRipple(&u, &v, cl.time);
			s32 su = ((s32)(u * 64.0f)) & 63;
			s32 tv = ((s32)(v * 64.0f)) & 63;
			*r_turb_pdest = *(r_turb_pbase + (tv << 6) + su);
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanLitDitheredHL()
{
	do {
		if (*pz <= (izi>>16) && D_Dither(r_turb_pdest,1-turb_opacity)) {
			f32 u = (f32)((r_turb_s >> 16) & 63) * (1.0f/64.0f);
			f32 v = (f32)((r_turb_t >> 16) & 63) * (1.0f/64.0f);
			HLWarp(&u, &v, cl.time);
			HLRipple(&u, &v, cl.time);
			s32 su = ((s32)(u * 64.0f)) & 63;
			s32 tv = ((s32)(v * 64.0f)) & 63;
			s32 pix = *(r_turb_pbase + (tv << 6) + su);
			s32 lit = *(litwater_base+(r_turb_pdest-d_viewbuffer));
			*r_turb_pdest = D_TurbMixLit(pix, lit);
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanFilteredDitheredHL()
{
	do {
		if (*pz <= (izi>>16) && D_Dither(r_turb_pdest,1-turb_opacity)) {
			s32 dither_idx = (flt_cur_x & 1) + ((flt_y & 1) << 1);
			s32 s_d = r_turb_s + dither_s[dither_idx];
			s32 t_d = r_turb_t + dither_t[dither_idx];
			f32 u = (f32)((s_d >> 16) & 63) * (1.0f/64.0f);
			f32 v = (f32)((t_d >> 16) & 63) * (1.0f/64.0f);
			HLWarp(&u, &v, cl.time);
			HLRipple(&u, &v, cl.time);
			s32 su = ((s32)(u * 64.0f)) & 63;
			s32 tv = ((s32)(v * 64.0f)) & 63;
			*r_turb_pdest = *(r_turb_pbase + (tv << 6) + su);
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
		flt_cur_x++;
	} while (--r_turb_spancount > 0);
}

static void D_DrawTurbulentSpanLitFilteredDitheredHL()
{
	do {
		if (*pz <= (izi>>16) && D_Dither(r_turb_pdest,1-turb_opacity)) {
			s32 dither_idx = (flt_cur_x & 1) + ((flt_y & 1) << 1);
			s32 s_d = r_turb_s + dither_s[dither_idx];
			s32 t_d = r_turb_t + dither_t[dither_idx];
			f32 u = (f32)((s_d >> 16) & 63) * (1.0f/64.0f);
			f32 v = (f32)((t_d >> 16) & 63) * (1.0f/64.0f);
			HLWarp(&u, &v, cl.time);
			HLRipple(&u, &v, cl.time);
			s32 su = ((s32)(u * 64.0f)) & 63;
			s32 tv = ((s32)(v * 64.0f)) & 63;
			s32 pix = *(r_turb_pbase + (tv << 6) + su);
			s32 lit = *(litwater_base+ (r_turb_pdest-d_viewbuffer));
			*r_turb_pdest = D_TurbMixLit(pix, lit);
		}
		r_turb_pdest++;
		izi += izistep;
		pz++;
		r_turb_s += r_turb_sstep;
		r_turb_t += r_turb_tstep;
		flt_cur_x++;
	} while (--r_turb_spancount > 0);
}

void Turbulent(espan_t *pspan, f32 opacity)
{
	turb_opacity = opacity;
	s32 turb_func_n = r_alphapass ? (r_alphastyle.value ? 2 : 1) : 0;
	turb_func_n += lmonly ? 6 : 0;
	turb_func_n += r_dithertex.value ? 3 : 0;
	turb_func_n += r_hlwater.value ? 12 : 0;
	if(r_hlwater.value) InitHLWaterLUT();
	if(lmonly) R_BuildLitLUT();
	if(r_alphastyle.value == 0) R_BuildColorMixLUT(0);
	void (*pturbdrawfunc)() = turbdrawfunc[turb_func_n];
	r_turb_turb = sintable + ((s32)(cl.time * SPEED) & (CYCLE - 1));
	r_turb_sstep = 0; // keep compiler happy
	r_turb_tstep = 0; // ditto
	r_turb_pbase = (u8 *)cacheblock;
	f32 sdivz16stepu = d_sdivzstepu * 16;
	f32 tdivz16stepu = d_tdivzstepu * 16;
	f32 zi16stepu = d_zistepu * 16;
	izistep = (s32)(d_zistepu * 0x8000 * 0x10000);
	do {
		r_turb_pdest = (u8 *) ((u8 *) d_viewbuffer + (screenwidth * pspan->v) + pspan->u);
		pz = d_pzbuffer + (d_zwidth * pspan->v) + pspan->u; // Manoel Kasimier - translucent water
		s32 count = pspan->count;
		f32 du = (f32)pspan->u; // calculate the initial s/z, t/z,
		f32 dv = (f32)pspan->v; // 1/z, s, and t and clamp
		f32 sdivz = d_sdivzorigin + dv*d_sdivzstepv + du*d_sdivzstepu;
		f32 tdivz = d_tdivzorigin + dv*d_tdivzstepv + du*d_tdivzstepu;
		f32 zi = d_ziorigin + dv * d_zistepv + du * d_zistepu;
		f32 z = (f32)0x10000 / zi; // prescale to 16.16 fixed-point
		izi = (s32)(zi * 0x8000 * 0x10000); // Manoel Kasimier - translucent water
		r_turb_s = (s32)(sdivz * z) + sadjust;
		if (r_turb_s > bbextents)
			r_turb_s = bbextents;
		else if (r_turb_s < 0)
			r_turb_s = 0;
		r_turb_t = (s32)(tdivz * z) + tadjust;
		if (r_turb_t > bbextentt)
			r_turb_t = bbextentt;
		else if (r_turb_t < 0)
			r_turb_t = 0;
		do {
			r_turb_spancount = count; // calculate s and t
			if (count >= 16) // at the far end of the span
				r_turb_spancount = 16;
			count -= r_turb_spancount;
			s32 snext, tnext;
			if (count) {
				// calculate s/z, t/z, zi->fixed s and t at far end of span,
				// calculate s and t steps across span by shifting
				sdivz += sdivz16stepu;
				tdivz += tdivz16stepu;
				zi += zi16stepu;
				z = (f32)0x10000 / zi; // prescale to 16.16 fixed-point
				snext = (s32)(sdivz * z) + sadjust;
				if (snext > bbextents)
					snext = bbextents;
				else if (snext < 16)
					snext = 16; // prevent round-off error on <0 steps from
				// from causing overstepping & running off the
				// edge of the texture
				tnext = (s32)(tdivz * z) + tadjust;
				if (tnext > bbextentt)
					tnext = bbextentt;
				else if (tnext < 16)
					tnext = 16;
				// guard against round-off error on <0 steps
				r_turb_sstep = (snext - r_turb_s) >> 4;
				r_turb_tstep = (tnext - r_turb_t) >> 4;
			} else {
				// calculate s/z, t/z, zi->fixed s and t at last pixel in span (so
				// can't step off polygon), clamp, calculate s and t steps across
				// span by division, biasing steps low so we don't run off the
				// texture
				f32 spancountminus1 = (f32)(r_turb_spancount - 1);
				sdivz += d_sdivzstepu * spancountminus1;
				tdivz += d_tdivzstepu * spancountminus1;
				zi += d_zistepu * spancountminus1;
				z = (f32)0x10000 / zi; // prescale to 16.16 fixed-point
				snext = (s32)(sdivz * z) + sadjust;
				if (snext > bbextents)
					snext = bbextents;
				else if (snext < 16)
					snext = 16; // prevent round-off error on <0 steps from
				// from causing overstepping & running off the
				// edge of the texture
				tnext = (s32)(tdivz * z) + tadjust;
				if (tnext > bbextentt)
					tnext = bbextentt;
				else if (tnext < 16)
					tnext = 16;// guard against round-off error on <0 steps
				if (r_turb_spancount > 1) {
					r_turb_sstep = (snext - r_turb_s) / (r_turb_spancount - 1);
					r_turb_tstep = (tnext - r_turb_t) / (r_turb_spancount - 1);
				}
			}
			r_turb_s = r_turb_s & ((CYCLE << 16) - 1);
			r_turb_t = r_turb_t & ((CYCLE << 16) - 1);
			if(r_dithertex.value) {
				s32 pixel_index = (s32)(r_turb_pdest -
							(u8*)screen->pixels);
				flt_y = pixel_index / scr_vrect.width;
				flt_x = pixel_index - flt_y * scr_vrect.width;
				flt_start_x = flt_x;
				flt_cur_x = flt_start_x;
			}
			pturbdrawfunc();
			r_turb_s = snext;
			r_turb_t = tnext;
		} while (count > 0);
	} while ((pspan = pspan->pnext) != NULL);
}

void D_DrawSpans(espan_t *pspan, s32 type, f32 opacity)
{
	u8 *pbase = (u8 *)cacheblock;
	f32 sdivz8stepu = d_sdivzstepu * 8;
	f32 tdivz8stepu = d_tdivzstepu * 8;
	f32 zi8stepu = d_zistepu * 8;
	izistep = (s32)(d_zistepu * 0x8000 * 0x10000);
	do {
		u8 *pdest = (u8 *)((u8 *) d_viewbuffer +
				      (screenwidth * pspan->v) + pspan->u);
		if (lmonly) {
			if (!litwater_base) {
				litwater_base = malloc(vid.width * vid.height);
				if (!litwater_base)
					Sys_Error("Not enough memory for lit water");
			}
			pdest = (u8 *)((u8 *) litwater_base +
				      (screenwidth * pspan->v) + pspan->u);
		}
		s32 count = pspan->count;
		f32 du = (f32)pspan->u; // calculate the initial s/z, t/z,
		f32 dv = (f32)pspan->v; // 1/z, s, and t and clamp
		f32 sdivz = d_sdivzorigin + dv*d_sdivzstepv + du*d_sdivzstepu;
		f32 tdivz = d_tdivzorigin + dv*d_tdivzstepv + du*d_tdivzstepu;
		f32 zi = d_ziorigin + dv * d_zistepv + du * d_zistepu;
		f32 z = (f32)0x10000 / zi; // prescale to 16.16 fixed-point
		if (type == SPAN_TRANS || type == SPAN_CUTOUT) {
			pz = d_pzbuffer + (d_zwidth * pspan->v) + pspan->u;
			izi = (s32)(zi * 0x8000 * 0x10000);
		}
		s32 s = (s32)(sdivz * z) + sadjust;
		if (s > bbextents)
			s = bbextents;
		else if (s < 0)
			s = 0;
		s32 t = (s32)(tdivz * z) + tadjust;
		if (t > bbextentt)
			t = bbextentt;
		else if (t < 0)
			t = 0;
		do {
			// calculate s and t at the far end of the span
			s32 snext, tnext;
			s32 sstep = 0; // keep compiler happy
			s32 tstep = 0; // ditto
			s32 spancount = count;
			if (count >= 8)
				spancount = 8;
			count -= spancount;
			if (count) {
				// calculate s/z, t/z, zi->fixed s and t at far end of span,
				// calculate s and t steps across span by shifting
				sdivz += sdivz8stepu;
				tdivz += tdivz8stepu;
				zi += zi8stepu;
				z = (f32)0x10000 / zi; // prescale to 16.16 fixed-point
				snext = (s32)(sdivz * z) + sadjust;
				if (snext > bbextents)
					snext = bbextents;
				else if (snext < 8)
					snext = 8; // prevent round-off error on <0 steps from
				// from causing overstepping & running off the
				// edge of the texture
				tnext = (s32)(tdivz * z) + tadjust;
				if (tnext > bbextentt)
					tnext = bbextentt;
				else if (tnext < 8)
					tnext = 8; // guard against round-off error on <0 steps
				sstep = (snext - s) >> 3;
				tstep = (tnext - t) >> 3;
			} else {
				// calculate s/z, t/z, zi->fixed s and t at last pixel in span (so
				// can't step off polygon), clamp, calculate s and t steps across
				// span by division, biasing steps low so we don't run off the
				// texture
				f32 spancountminus1 = (f32)(spancount - 1);
				sdivz += d_sdivzstepu * spancountminus1;
				tdivz += d_tdivzstepu * spancountminus1;
				zi += d_zistepu * spancountminus1;
				z = (f32)0x10000 / zi; // prescale to 16.16 fixed-point
				snext = (s32)(sdivz * z) + sadjust;
				if (snext > bbextents)
					snext = bbextents;
				else if (snext < 8)
					snext = 8; // prevent round-off error on <0 steps from
				// from causing overstepping & running off the
				// edge of the texture
				tnext = (s32)(tdivz * z) + tadjust;
				if (tnext > bbextentt)
					tnext = bbextentt;
				else if (tnext < 8)
					tnext = 8; // guard against round-off error on <0 steps
				if (spancount > 1) {
					sstep = (snext - s) / (spancount - 1);
					tstep = (tnext - t) / (spancount - 1);
				}
			}
			if (type == SPAN_NORMAL) {
				do {
					/*
					// CyanBun96: this is a working implementation of bilinear texture filtering. Don't use. 
					s32 sx = s >> 16;
					s32 sy = t >> 16;
					s32 fracx = (s >> 8) & 0xFF;
					s32 fracy = (t >> 8) & 0xFF;
					// clamp to prevent out-of-bounds
					if (sx < 0) sx = 0;
					if (sy < 0) sy = 0;
					if (sx >= cachewidth - 1) { sx = cachewidth - 2; fracx = 255; }
					if (sy >= cacheheight - 1) { sy = cacheheight - 2; fracy = 255; }
					u8 *row0 = pbase + sy * cachewidth;
					u8 *row1 = row0 + cachewidth;
					u8 i00 = row0[sx];
					u8 i10 = row0[sx + 1];
					u8 i01 = row1[sx];
					u8 i11 = row1[sx + 1];
					s32 lutx = fracx >> 3;   // 0..31 to fit in FOG_LUT_LEVELS
					s32 luty = fracy >> 3;   // 0..31
					if (!fog_lut_built)
						R_BuildColorMixLUT(0);
					// horizontal blends
					u8 top = color_mix_lut[i00][i10][lutx];
					u8 bot = color_mix_lut[i01][i11][lutx];
					// vertical blend
					u8 result = color_mix_lut[top][bot][luty];
					*pdest = result;
					pdest++;
					s += sstep;
					t += tstep;
					*/
					u8 pix = *(pbase + (s >> 16) +
							(t >> 16) * cachewidth);
					*pdest = pix;
					pdest++;
					s += sstep;
					t += tstep;
				} while (--spancount > 0);
			} else if (type == SPAN_CUTOUT) {
				do {
					if (*pz <= (izi >> 16)) {
						u8 pix = *(pbase + (s >> 16) +
							(t >> 16) * cachewidth);
						cutoutbuf[pdest-d_viewbuffer] = 0;
						if (pix != 0xff) {
							*pdest = pix;
							cutoutbuf[pdest-d_viewbuffer] = 1;
						}
					}
					pdest++;
					izi += izistep;
					pz++;
					s += sstep;
					t += tstep;
				} while (--spancount > 0);
			} else if (type == SPAN_SKYBOX) {
				s32 foglut = r_skyfog.value*FOG_LUT_LEVELS;
				do {
					u8 pix = *(pbase + (s >> 16) +
						(t >> 16) * cachewidth);
					if (fog_density > 0)
						pix = color_mix_lut[pix][fog_pal_index][foglut];
					*pdest = pix;
					pdest++;
					s += sstep;
					t += tstep;
				} while (--spancount > 0);
			} else if (type == SPAN_TRANS) {
				s32 foglut = opacity*FOG_LUT_LEVELS;
				if (r_alphastyle.value == 0) {
					if (!fog_lut_built)
						R_BuildColorMixLUT(0);
					do {
						if (*pz <= (izi >> 16)) {
							u8 pix = *(pbase + (s >> 16) + (t >> 16) * cachewidth);
							if (pix != 0xff) {
								pix = color_mix_lut[pix][*pdest][foglut];
								*pdest = pix;
							}
						}
						pdest++;
						izi += izistep;
						pz++;
						s += sstep;
						t += tstep;
					} while (--spancount > 0);
				} else {
					do {
						if (*pz <= (izi >> 16) && D_Dither(pdest, 1-opacity)) {
							u8 pix = *(pbase + (s >> 16) + (t >> 16) * cachewidth);
							if (pix != 0xff) *pdest = pix;
						}
						pdest++;
						izi += izistep;
						pz++;
						s += sstep;
						t += tstep;
					} while (--spancount > 0);
				}
			}
			s = snext;
			t = tnext;
		} while (count > 0);
	} while ((pspan = pspan->pnext) != NULL);
}

void D_DrawSpansDithered(espan_t *pspan, s32 type, f32 opacity)
{ // Separate from the regular function to avoid branching within the loops
	u8 *pbase = (u8 *)cacheblock;
	f32 sdivz8stepu = d_sdivzstepu * 8;
	f32 tdivz8stepu = d_tdivzstepu * 8;
	f32 zi8stepu = d_zistepu * 8;
	izistep = (s32)(d_zistepu * 0x8000 * 0x10000);
	do {
		u8 *pdest = (u8 *)((u8 *) d_viewbuffer +
				      (screenwidth * pspan->v) + pspan->u);
		if (lmonly) {
			if (!litwater_base) {
				litwater_base = malloc(vid.width * vid.height);
				if (!litwater_base)
					Sys_Error("Not enough memory for lit water");
			}
			pdest = (u8 *)((u8 *) litwater_base +
				      (screenwidth * pspan->v) + pspan->u);
		}
		s32 count = pspan->count;
		f32 du = (f32)pspan->u; // calculate the initial s/z, t/z,
		f32 dv = (f32)pspan->v; // 1/z, s, and t and clamp
		f32 sdivz = d_sdivzorigin + dv*d_sdivzstepv + du*d_sdivzstepu;
		f32 tdivz = d_tdivzorigin + dv*d_tdivzstepv + du*d_tdivzstepu;
		f32 zi = d_ziorigin + dv * d_zistepv + du * d_zistepu;
		f32 z = (f32)0x10000 / zi; // prescale to 16.16 fixed-point
		if (type == SPAN_TRANS || type == SPAN_CUTOUT) {
			pz = d_pzbuffer + (d_zwidth * pspan->v) + pspan->u;
			izi = (s32)(zi * 0x8000 * 0x10000);
		}
		s32 s = (s32)(sdivz * z) + sadjust;
		if (s > bbextents)
			s = bbextents;
		else if (s < 0)
			s = 0;
		s32 t = (s32)(tdivz * z) + tadjust;
		if (t > bbextentt)
			t = bbextentt;
		else if (t < 0)
			t = 0;
		do {
			// calculate s and t at the far end of the span
			s32 snext, tnext;
			s32 sstep = 0; // keep compiler happy
			s32 tstep = 0; // ditto
			s32 spancount = count;
			if (count >= 8)
				spancount = 8;
			count -= spancount;
			if (count) {
				// calculate s/z, t/z, zi->fixed s and t at far end of span,
				// calculate s and t steps across span by shifting
				sdivz += sdivz8stepu;
				tdivz += tdivz8stepu;
				zi += zi8stepu;
				z = (f32)0x10000 / zi; // prescale to 16.16 fixed-point
				snext = (s32)(sdivz * z) + sadjust;
				if (snext > bbextents)
					snext = bbextents;
				else if (snext < 8)
					snext = 8; // prevent round-off error on <0 steps from
				// from causing overstepping & running off the
				// edge of the texture
				tnext = (s32)(tdivz * z) + tadjust;
				if (tnext > bbextentt)
					tnext = bbextentt;
				else if (tnext < 8)
					tnext = 8; // guard against round-off error on <0 steps
				sstep = (snext - s) >> 3;
				tstep = (tnext - t) >> 3;
			} else {
				// calculate s/z, t/z, zi->fixed s and t at last pixel in span (so
				// can't step off polygon), clamp, calculate s and t steps across
				// span by division, biasing steps low so we don't run off the
				// texture
				f32 spancountminus1 = (f32)(spancount - 1);
				sdivz += d_sdivzstepu * spancountminus1;
				tdivz += d_tdivzstepu * spancountminus1;
				zi += d_zistepu * spancountminus1;
				z = (f32)0x10000 / zi; // prescale to 16.16 fixed-point
				snext = (s32)(sdivz * z) + sadjust;
				if (snext > bbextents)
					snext = bbextents;
				else if (snext < 8)
					snext = 8; // prevent round-off error on <0 steps from
				// from causing overstepping & running off the
				// edge of the texture
				tnext = (s32)(tdivz * z) + tadjust;
				if (tnext > bbextentt)
					tnext = bbextentt;
				else if (tnext < 8)
					tnext = 8; // guard against round-off error on <0 steps
				if (spancount > 1) {
					sstep = (snext - s) / (spancount - 1);
					tstep = (tnext - t) / (spancount - 1);
				}
			}
			// CyanBun96: dithered sampling from Unreal
			s32 pixel_index = (s32)(pdest - (u8*)screen->pixels);
			s32 y = pixel_index / scr_vrect.width;
			s32 x = pixel_index - y * scr_vrect.width;
			s32 start_x = x;
			s32 cur_x = start_x;
			if (type == SPAN_NORMAL) {
				do {
					s32 dither_idx = (cur_x & 1) + ((y & 1) << 1);
					// Apply dither offset
					s32 s_d = s + dither_s[dither_idx];
					s32 t_d = t + dither_t[dither_idx];
					// Clamp to valid texel range
					if (s_d < 0) s_d = 0;
					if (t_d < 0) t_d = 0;
					s32 s_max = (cachewidth  - 1) << 16;
					s32 t_max = (cacheheight - 1) << 16;
					if (s_d > s_max) s_d = s_max;
					if (t_d > t_max) t_d = t_max;
					u8 pix = *(pbase + (s_d >> 16) + ((t_d >> 16) * cachewidth));
					*pdest++ = pix;
					s += sstep;
					t += tstep;
					cur_x++;
				} while (--spancount > 0);
			} else if (type == SPAN_CUTOUT) {
				do {
					s32 dither_idx = (cur_x & 1) + ((y & 1) << 1);
					s32 s_d = s + dither_s[dither_idx];
					s32 t_d = t + dither_t[dither_idx];
					if (s_d < 0) s_d = 0;
					if (t_d < 0) t_d = 0;
					s32 s_max = (cachewidth  - 1) << 16;
					s32 t_max = (cacheheight - 1) << 16;
					if (s_d > s_max) s_d = s_max;
					if (t_d > t_max) t_d = t_max;
					if (*pz <= (izi >> 16)) {
						u8 pix = *(pbase + (s_d >> 16) + ((t_d >> 16) * cachewidth));
						cutoutbuf[pdest-d_viewbuffer] = 0;
						if (pix != 0xff) {
							*pdest = pix;
							cutoutbuf[pdest-d_viewbuffer] = 1;
						}
					}
					pdest++;
					izi += izistep;
					pz++;
					s += sstep;
					t += tstep;
					cur_x++;
				} while (--spancount > 0);
			} else if (type == SPAN_SKYBOX) {
				s32 foglut = r_skyfog.value*FOG_LUT_LEVELS;
				do {
					s32 dither_idx = (cur_x & 1) + ((y & 1) << 1);
					s32 s_d = s + dither_s[dither_idx];
					s32 t_d = t + dither_t[dither_idx];
					if (s_d < 0) s_d = 0;
					if (t_d < 0) t_d = 0;
					s32 s_max = (cachewidth  - 1) << 16;
					s32 t_max = (cacheheight - 1) << 16;
					if (s_d > s_max) s_d = s_max;
					if (t_d > t_max) t_d = t_max;
					u8 pix = *(pbase + (s_d >> 16) + ((t_d >> 16) * cachewidth));
					if (fog_density > 0)
						pix = color_mix_lut[pix][fog_pal_index][foglut];
					*pdest++ = pix;
					s += sstep;
					t += tstep;
					cur_x++;
				} while (--spancount > 0);
			} else if (type == SPAN_TRANS) {
				s32 foglut = opacity*FOG_LUT_LEVELS;
				if (r_alphastyle.value == 0) {
					if (!fog_lut_built)
						R_BuildColorMixLUT(0);
					do {
						s32 dither_idx = (cur_x & 1) + ((y & 1) << 1);
						s32 s_d = s + dither_s[dither_idx];
						s32 t_d = t + dither_t[dither_idx];
						if (s_d < 0) s_d = 0;
						if (t_d < 0) t_d = 0;
						s32 s_max = (cachewidth  - 1) << 16;
						s32 t_max = (cacheheight - 1) << 16;
						if (s_d > s_max) s_d = s_max;
						if (t_d > t_max) t_d = t_max;
						if (*pz <= (izi >> 16)) {
							u8 pix = *(pbase + (s_d >> 16) + ((t_d >> 16) * cachewidth));
							if (pix != 0xff) {
								pix = color_mix_lut[pix][*pdest][foglut];
								*pdest = pix;
							}
						}
						pdest++;
						izi += izistep;
						pz++;
						s += sstep;
						t += tstep;
						cur_x++;
					} while (--spancount > 0);
				} else {
					do {
						s32 dither_idx = (cur_x & 1) + ((y & 1) << 1);
						s32 s_d = s + dither_s[dither_idx];
						s32 t_d = t + dither_t[dither_idx];
						if (s_d < 0) s_d = 0;
						if (t_d < 0) t_d = 0;
						s32 s_max = (cachewidth  - 1) << 16;
						s32 t_max = (cacheheight - 1) << 16;
						if (s_d > s_max) s_d = s_max;
						if (t_d > t_max) t_d = t_max;
						if (*pz <= (izi >> 16) && D_Dither(pdest, 1-opacity)) {
							u8 pix = *(pbase + (s_d >> 16) + ((t_d >> 16) * cachewidth));
							if (pix != 0xff) *pdest = pix;
						}
						pdest++;
						izi += izistep;
						pz++;
						s += sstep;
						t += tstep;
						cur_x++;
					} while (--spancount > 0);
				}
			}
			s = snext;
			t = tnext;
		} while (count > 0);
	} while ((pspan = pspan->pnext) != NULL);
}

void D_DrawZSpans(espan_t *pspan)
{
	s32 izistep = (s32)(d_zistepu * 0x8000 * 0x10000);
	do {
		s16 *pdest = d_pzbuffer + (d_zwidth * pspan->v) + pspan->u;
		s32 count = pspan->count;
		f32 du = (f32)pspan->u; // calculate the initial 1/z
		f32 dv = (f32)pspan->v;
		f64 zi = d_ziorigin + dv * d_zistepv + du * d_zistepu;
		s32 izi = (s32)(zi * 0x8000 * 0x10000);
		if ((intptr_t) pdest & 0x02) {
			*pdest++ = (s16)(izi >> 16);
			izi += izistep;
			count--;
		}
		s32 doublecount = count >> 1;
		if (doublecount > 0) {
			do {
				u32 ltemp = izi >> 16;
				izi += izistep;
				ltemp |= izi & 0xFFFF0000;
				izi += izistep;
				*(s32 *)pdest = ltemp;
				pdest += 2;
			} while (--doublecount > 0);
		}
		if (count & 1)
			*pdest = (s16)(izi >> 16);
	} while ((pspan = pspan->pnext) != NULL);
}

void D_DrawZSpansTrans(espan_t *pspan)
{
	s32 izistep = (s32)(d_zistepu * 0x8000 * 0x10000);
	do {
		s16 *pdest = d_pzbuffer + (d_zwidth * pspan->v) + pspan->u;
		s32 count = pspan->count;
		f32 du = (f32)pspan->u; // calculate the initial 1/z
		f32 dv = (f32)pspan->v;
		f64 zi = d_ziorigin + dv * d_zistepv + du * d_zistepu;
		s32 izi = (s32)(zi * 0x8000 * 0x10000);
		if ((intptr_t) pdest & 0x02) {
			if(cutoutbuf[(pdest-d_pzbuffer)] == 1 &&
				*pdest < (s16)(izi >> 16))
				*pdest++ = (s16)(izi >> 16);
			izi += izistep;
			count--;
		}
		s32 doublecount = count >> 1;
		if (doublecount > 0) {
			do {
				u32 ltemp = izi >> 16;
				izi += izistep;
				ltemp |= izi & 0xFFFF0000;
				izi += izistep;
				if(cutoutbuf[(pdest-d_pzbuffer)] == 1 &&
					*pdest < (s16)ltemp)
					*(s32 *)pdest = ltemp;
				pdest += 2;
			} while (--doublecount > 0);
		}
		if (count & 1) {
			if(cutoutbuf[(pdest-d_pzbuffer)] == 1 &&
				*pdest < (s16)(izi >> 16))
				*pdest = (s16)(izi >> 16);
		}
	} while ((pspan = pspan->pnext) != NULL);
}
