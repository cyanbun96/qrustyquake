#include "quakedef.h"

static f32 gamma_lut[256];
static s32 color_conv_initialized = 0;
static lab_t lab_palette[256];

f32 lab_f(f32 t) // cube-root if > epsilon, linear otherwise
{ return (t > 0.008856f) ? cbrtf(t) : ((903.3f * t + 16.0f) / 116.0f); }

lab_t rgb_lin_to_lab(f32 r, f32 g, f32 b)
{
	f32 x = r*0.4124f + g*0.3576f + b*0.1805f;
	f32 y = r*0.2126f + g*0.7152f + b*0.0722f;
	f32 z = r*0.0193f + g*0.1192f + b*0.9505f;
	x /= 0.95047f; // normalize by d65 white
	z /= 1.08883f;
	f32 fx = lab_f(x);
	f32 fy = lab_f(y);
	f32 fz = lab_f(z);
	lab_t lab = {116.0f*fy-16.0f, 500.0f*(fx-fy), 200.0f*(fy-fz)};
	return lab;
}

void init_color_conv()
{
	if(color_conv_initialized) return;
	color_conv_initialized = 1;
	for(s32 i = 0; i < 256; ++i){
		f32 v = i / 255.0f;
		gamma_lut[i] = (v > 0.04045f) ?
			powf((v + 0.055f) / 1.055f, 2.4f) : (v / 12.92f);
	}
	for(s32 i = 0, p = 0; i < 256; ++i, p += 3){
		f32 lr = gamma_lut[CURWORLDPAL[p + 0]];
		f32 lg = gamma_lut[CURWORLDPAL[p + 1]];
		f32 lb = gamma_lut[CURWORLDPAL[p + 2]];
		lab_palette[i] = rgb_lin_to_lab(lr, lg, lb);
	}
}

u8 rgbtoi_lab(u8 r, u8 g, u8 b)
{
	init_color_conv();
	lab_t lab0 = rgb_lin_to_lab(gamma_lut[r], gamma_lut[g], gamma_lut[b]);
	f32 bestdist = FLT_MAX;
	u8 besti = 0;
	for(s32 i = 0; i < 256; ++i){
		f32 dl = lab0.l - lab_palette[i].l;
		f32 da = lab0.a - lab_palette[i].a;
		f32 db = lab0.b - lab_palette[i].b;
		f32 dist = dl * dl + da * da + db * db;
		if (dist < bestdist){
			bestdist = dist;
			besti = (u8)i;
			if (!dist) break;
		}
	}
	return besti;
}

u8 rgbtoi(u8 r, u8 g, u8 b)
{
	u8 besti = 0;
	s32 bestdist = 0x7fffffff;
	u8 *p = CURWORLDPAL;
	for(s32 i = 0; i < 256; ++i, p += 3){
		s32 dr = r - p[0]; // squared euclidean distance
		s32 dg = g - p[1];
		s32 db = b - p[2];
		s32 dist = dr*dr + dg*dg + db*db;
		if (dist < bestdist){
			bestdist = dist;
			besti = (u8)i;
			if (!dist) return besti;
		}
	}
	return besti;
}

void R_BuildColorMixLUT(SDL_UNUSED cvar_t *cvar)
{
	if(fog_lut_built) return;
	u8 (*convfunc)(u8,u8,u8) = r_labmixpal.value == 1 ? rgbtoi_lab : rgbtoi;
	for(s32 c1 = 0; c1 < 256; c1++){
	for(s32 c2 = 0; c2 < 256; c2++){
		u8 r1 = CURWORLDPAL[c1*3+0];
		u8 g1 = CURWORLDPAL[c1*3+1];
		u8 b1 = CURWORLDPAL[c1*3+2];
		u8 r2 = CURWORLDPAL[c2*3+0];
		u8 g2 = CURWORLDPAL[c2*3+1];
		u8 b2 = CURWORLDPAL[c2*3+2];
		for(s32 level = 0; level < FOG_LUT_LEVELS; level++){
			f32 factor = (f32)level / (FOG_LUT_LEVELS-1);
			u8 r = (u8)(r1 + factor * (r2 - r1)); // lerp each color
			u8 g = (u8)(g1 + factor * (g2 - g1));
			u8 b = (u8)(b1 + factor * (b2 - b1));
			u8 mixed_index = convfunc(r, g, b);
			color_mix_lut[c1][c2][level] = mixed_index;
	}}}
	fog_lut_built = 1;
}

void R_BuildLitLUT()
{
	if(lit_lut_initialized) return;
	u8 (*convfunc)(u8,u8,u8) = r_labmixpal.value == 1 ? rgbtoi_lab : rgbtoi;
	for(s32 r = 0; r < LIT_LUT_RES; ++r){
	for(s32 g = 0; g < LIT_LUT_RES; ++g){
	for(s32 b = 0; b < LIT_LUT_RES; ++b){
		s32 rr = (r * 255) / (LIT_LUT_RES - 1);
		s32 gg = (g * 255) / (LIT_LUT_RES - 1);
		s32 bb = (b * 255) / (LIT_LUT_RES - 1);
		lit_lut[r][g][b] = convfunc(rr,gg,bb);
	}}}
	lit_lut_initialized = 1;
}

