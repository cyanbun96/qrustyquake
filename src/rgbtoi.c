#include "quakedef.h"

#define NN_CELLS 16 // per axis, in the grids of the nearest color searches
#define NN_NCELLS (NN_CELLS * NN_CELLS * NN_CELLS)
typedef struct { f32 l, a, b; } lab_t;
typedef struct { // the colors that can be nearest somewhere in each cell
	f32 lo[3], step[3]; // where the grid starts, cell size
	u32 first[NN_NCELLS + 1]; // a cell's colors start at list[first[c]]
	u8 *list; // cell by cell, each cell's in palette order
	s32 cap;
} nngrid_t;
static f32 gamma_lut[256];
static lab_t lab_palette[256];
static s32 color_conv_initialized = 0;
static nngrid_t nn_rgb, nn_lab;
static u8 nn_all[256]; // all colors, as a list

f32 lab_f(f32 t) // cube-root if > epsilon, linear otherwise
{ return(t > 0.008856f) ? cbrtf(t) : ((903.3f * t + 16.0f) / 116.0f); }

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
	for(s32 i = 0; i < 256; ++i) nn_all[i] = i;
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

static u8 NN_ScanLab(lab_t lab0, const u8 *list, s32 n)
{ // the first of the listed colors nearest to lab0
	f32 bestdist = FLT_MAX;
	u8 besti = 0;
	for(s32 k = 0; k < n; ++k){
		lab_t *p = &lab_palette[list[k]];
		f32 dl = lab0.l - p->l;
		f32 da = lab0.a - p->a;
		f32 db = lab0.b - p->b;
		f32 dist = dl * dl + da * da + db * db;
		if(dist < bestdist){
			bestdist = dist;
			besti = list[k];
			if(!dist) break;
		}
	}
	return besti;
}

u8 rgbtoi_lab(u8 r, u8 g, u8 b)
{
	init_color_conv();
	lab_t lab0 = rgb_lin_to_lab(gamma_lut[r], gamma_lut[g], gamma_lut[b]);
	return NN_ScanLab(lab0, nn_all, 256);
}

// The first of the listed colors nearest to r, g, b, all of them without a
// list. Inlined, the table builds go through it for every entry.
Q_SIL u8 NN_ScanRGB(u8 r, u8 g, u8 b, const u8 *list, s32 n)
{
	const u8 *pal = CURWORLDPAL;
	u8 besti = 0;
	s32 bestdist = 0x7fffffff;
	for(s32 k = 0; k < n; ++k){
		s32 i = list ? list[k] : k;
		const u8 *p = pal + i * 3;
		s32 dr = r - p[0]; // squared euclidean distance
		s32 dg = g - p[1];
		s32 db = b - p[2];
		s32 dist = dr*dr + dg*dg + db*db;
		if(dist < bestdist){
			bestdist = dist;
			besti = i;
			if(!dist) break;
		}
	}
	return besti;
}

u8 rgbtoi(u8 r, u8 g, u8 b){ return NN_ScanRGB(r, g, b, NULL, 256); }

static void NN_Build(nngrid_t *g, f32 (*p)[3], const f32 *lo, const f32 *hi)
{ // a color can't be the nearest anywhere in a cell if another color is
  // nearer to all of the cell than any of it is to the first one
	u32 n = 0;
	for(s32 a = 0; a < 3; ++a){
		g->lo[a] = lo[a];
		g->step[a] = (hi[a] - lo[a]) / NN_CELLS;
	}
	for(s32 c = 0; c < NN_NCELLS; ++c){
		s32 ci[3] = {c / NN_CELLS / NN_CELLS, c / NN_CELLS % NN_CELLS,
			c % NN_CELLS};
		f32 dnear[256], best = FLT_MAX;
		for(s32 i = 0; i < 256; ++i){
			f32 dn = 0, df = 0;
			for(s32 a = 0; a < 3; ++a){ // a bit bigger cell
				f32 b0 = lo[a] + ci[a] * g->step[a] - 0.01f;
				f32 d0 = b0 - p[i][a];
				f32 d1 = p[i][a] - (b0 + g->step[a] + 0.02f);
				f32 in = d0 > 0 ? d0 : d1 > 0 ? d1 : 0;
				f32 out = q_max(fabsf(d0), fabsf(d1));
				dn += in * in;
				df += out * out;
			}
			dnear[i] = dn;
			if(df < best) best = df;
		}
		best = best * 1.001f + 0.01f; // rounding can't drop a tie
		g->first[c] = n;
		for(s32 i = 0; i < 256; ++i){
			if(dnear[i] > best) continue;
			if(n == (u32)g->cap){
				g->cap = g->cap ? g->cap * 2 : NN_NCELLS * 16;
				if(!(g->list = realloc(g->list, g->cap)))
					Sys_Error("NN_Build: out of memory");
			}
			g->list[n++] = i;
		}
	}
	g->first[NN_NCELLS] = n;
}

static void NN_Prepare()
{ // the grid of the search r_labmixpal picks, for the current palette
	f32 p[256][3];
	if(r_labmixpal.value == 1){
		static const f32 lo[3] = {-1, -100, -120};
		static const f32 hi[3] = {101, 110, 110};
		init_color_conv();
		for(s32 i = 0; i < 256; ++i){
			p[i][0] = lab_palette[i].l;
			p[i][1] = lab_palette[i].a;
			p[i][2] = lab_palette[i].b;
		}
		NN_Build(&nn_lab, p, lo, hi);
		return;
	}
	static const f32 lo[3] = {0, 0, 0}, hi[3] = {256, 256, 256};
	for(s32 i = 0; i < 256 * 3; ++i) p[i / 3][i % 3] = CURWORLDPAL[i];
	NN_Build(&nn_rgb, p, lo, hi);
}

static u8 rgbtoi_nn(u8 r, u8 g, u8 b)
{ // rgbtoi through nn_rgb, which NN_Prepare built
	s32 c = (((r * NN_CELLS) >> 8) * NN_CELLS + ((g * NN_CELLS) >> 8))
		* NN_CELLS + ((b * NN_CELLS) >> 8);
	return NN_ScanRGB(r, g, b, nn_rgb.list + nn_rgb.first[c],
			nn_rgb.first[c + 1] - nn_rgb.first[c]);
}

static u8 rgbtoi_lab_nn(u8 r, u8 g, u8 b)
{ // rgbtoi_lab through nn_lab, which NN_Prepare built
	lab_t lab0 = rgb_lin_to_lab(gamma_lut[r], gamma_lut[g], gamma_lut[b]);
	f32 v[3] = {lab0.l, lab0.a, lab0.b};
	s32 c = 0;
	for(s32 a = 0; a < 3; ++a){
		f32 f = (v[a] - nn_lab.lo[a]) / nn_lab.step[a];
		if(!(f >= 0 && f < NN_CELLS)) // off the grid
			return NN_ScanLab(lab0, nn_all, 256);
		c = c * NN_CELLS + (s32)f;
	}
	return NN_ScanLab(lab0, nn_lab.list + nn_lab.first[c],
			nn_lab.first[c + 1] - nn_lab.first[c]);
}

static void R_MixRGB(u8 (*rgb)[3], const u8 *p1, const u8 *p2)
{ // the color of each level of mixing p2 into p1
	for(s32 level = 0; level < FOG_LUT_LEVELS; level++){
		f32 factor = (f32)level / (FOG_LUT_LEVELS-1);
		for(s32 i = 0; i < 3; ++i) // lerp each color
			rgb[level][i] = (u8)(p1[i] + factor * (p2[i] - p1[i]));
	}
}

void R_BuildColorMixLUT(SDL_UNUSED cvar_t *cvar)
{
	if(fog_lut_built) return;
	NN_Prepare();
	u8(*convfunc)(u8,u8,u8) = r_labmixpal.value==1?rgbtoi_lab_nn:rgbtoi_nn;
	u8 *pal = CURWORLDPAL;
	for(s32 c1 = 0; c1 < 256; c1++){
	for(s32 c2 = 0; c2 < 256; c2++){
		// mixing c2 into c1 is mixing c1 into c2 the other way round,
		// done already when c2 < c1, and its colors are the same ones
		// unless rounding says otherwise (one call site rounds alike)
		u8 rgb[2][FOG_LUT_LEVELS][3];
		for(s32 k = 0; k <= (c2 <= c1); ++k)
			R_MixRGB(rgb[k], pal + 3*(k?c2:c1), pal + 3*(k?c1:c2));
		for(s32 level = 0; level < FOG_LUT_LEVELS; level++){
			s32 m = FOG_LUT_LEVELS - 1 - level;
			u8 *c = rgb[0][level];
			if((c2 < c1 || (c2 == c1 && m < level))
					&& !memcmp(c, rgb[1][m], 3))
				color_mix_lut[c1][c2][level] =
					color_mix_lut[c2][c1][m];
			else color_mix_lut[c1][c2][level] =
				convfunc(c[0], c[1], c[2]);
	}}}
	fog_lut_built = 1;
}

void R_LabMixPalChanged(SDL_UNUSED cvar_t *cvar)
{ // what got mixed the other way gets mixed again when it's needed
	fog_lut_built = lit_lut_initialized = 0;
	fog_initialized = 0; // R_FogMix doesn't check fog_lut_built
	D_FlushCaches(0); // surfaces get lit through lit_lut
}

void R_WorldPalChanged()
{ // what was built from the old colors gets built again when it's needed
	color_conv_initialized = 0;
	refresh_palette++;
	R_LabMixPalChanged(0);
	Draw_PaletteChanged();
	Sky_PaletteChanged();
	Fog_SetPalIndex(0);
}

void R_BuildLitChan()
{
	u8 *pal = CURWORLDPAL;
	for(s32 i = 0; i < 256; ++i){
		lit_chan[0][i] = QUANT(pal[i*3+0]) * LIT_LUT_RES * LIT_LUT_RES;
		lit_chan[1][i] = QUANT(pal[i*3+1]) * LIT_LUT_RES;
		lit_chan[2][i] = QUANT(pal[i*3+2]);
	}
	lit_chan_pal = pal;
}

void R_BuildLitLUT()
{
	if(lit_lut_initialized) return;
	R_BuildLitChan();
	NN_Prepare();
	u8(*convfunc)(u8,u8,u8) = r_labmixpal.value==1?rgbtoi_lab_nn:rgbtoi_nn;
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
