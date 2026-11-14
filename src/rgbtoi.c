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

static void NN_PrepareRGB()
{ // the grid of the RGB search, for the current palette
	f32 p[256][3];
	static const f32 lo[3] = {0, 0, 0}, hi[3] = {256, 256, 256};
	for(s32 i = 0; i < 256 * 3; ++i) p[i / 3][i % 3] = CURWORLDPAL[i];
	NN_Build(&nn_rgb, p, lo, hi);
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
	NN_PrepareRGB();
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

// RGB nearest colours by the top 6 bits of each channel: within 2 of each
// channel, and the palette's own colours map to themselves (the first one,
// if two share a cell)
#define NN_LUTCELL(c) ((c)[0] >> 2 << 12 | (c)[1] >> 2 << 6 | (c)[2] >> 2)
static u8 *nn_lut;
static u8 nn_lut_pal[768]; // what it was built from

static void NN_PrepareLUT()
{
	const u8 *pal = CURWORLDPAL;
	if(nn_lut && !memcmp(nn_lut_pal, pal, 768))
		return;
	if(!nn_lut && !(nn_lut = malloc(64 * 64 * 64)))
		Sys_Error("NN_PrepareLUT: out of memory");
	memcpy(nn_lut_pal, pal, 768);
	NN_PrepareRGB();
	for(s32 i = 0; i < 64 * 64 * 64; i++)
		nn_lut[i] = rgbtoi_nn((i >> 12) * 4 + 2,
				(i >> 6 & 63) * 4 + 2, (i & 63) * 4 + 2);
	for(s32 i = 255; i >= 0; i--)
		nn_lut[NN_LUTCELL(pal + i * 3)] = i;
}

static f32 Dot3(const f32 *a, const f32 *b)
{ return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

static s32 Pal_Dist2(const u8 *a, const u8 *b)
{
	s32 dr = a[0] - b[0], dg = a[1] - b[1], db = a[2] - b[2];
	return dr * dr + dg * dg + db * db;
}

// RGBA pixels to palette indexes in place, Floyd-Steinberg on alternating
// row directions with the error scaled by dither. Where the palette can't
// reach a colour, full diffusion sprinkles far colours over it (blue dots
// on a blue sky): the error a pixel takes is capped, and it only gets
// colours near its own nearest one. The first rows are run twice: without
// error from above, they'd pass a lopsided one down that the next row
// overshoots with a line of far colours. The error that goes past a row's
// ends goes to its edge pixels instead, or they'd make lines at the seams.
#define DITHER_MAXERR 40 // per channel
#define DITHER_SPREAD 64
#define DITHER_WARMUP 2 // rows, even to keep row 0's direction
typedef struct {
	u8 nearby[256][256]; // the colours in each one's spread
	s32 nnear[256];
	u32 outside[256][8]; // bits of the colours out of it
	u32 seenkey[4096]; // the last searches, cell << 8 | n0
	u8 seenc[4096];
} ditherpal_t;

// What each index of a dithered image stands for: the colors that make each
// tile's dither, tiles a 32nd of the image across, average to the tile's
// source, held near the palette's own. Fogging those instead of the palette
// colors undoes some of the hue a dither lost to the palette: the light
// colors of Quake's are pinker than most skies, which fog turns grey.
#define MEANING_TILES 32 // across an image
#define MEANING_HOLD 0.1 // how near the palette's colors the meanings stay
struct dithermeaning_s {
	f64 m[256][256]; // over the tiles, index shares times each other
	f64 v[256][3]; // shares times the tile's mean source color
	f64 u[256]; // the shares
};

dithermeaning_t *R_DitherMeaningBegin()
{
	dithermeaning_t *dm = calloc(1, sizeof(*dm));
	if(!dm)
		Sys_Error("R_DitherMeaningBegin: out of memory");
	return dm;
}

static void Dither_AddTiles(dithermeaning_t *dm, u32 *cnt, f32 *sum,
		s32 ncol)
{ // a row of tiles' index counts and source sums, which it clears
	for(s32 t = 0; t < ncol; t++, cnt += 256, sum += 3){
		u8 used[256];
		f64 a[256];
		s32 n = 0, total = 0;
		for(s32 c = 0; c < 256; c++)
			if(cnt[c]){
				used[n++] = c;
				total += cnt[c];
			}
		for(s32 i = 0; i < n; i++)
			a[i] = (f64)cnt[used[i]] / total;
		for(s32 i = 0; i < n; i++){
			f64 *row = dm->m[used[i]];
			for(s32 j = 0; j < n; j++)
				row[used[j]] += a[i] * a[j];
			for(s32 k = 0; k < 3; k++)
				dm->v[used[i]][k] += a[i] * sum[k] / total;
			dm->u[used[i]] += a[i];
		}
		memset(cnt, 0, 256 * sizeof(*cnt));
		memset(sum, 0, 3 * sizeof(*sum));
	}
}

void R_DitherMeaningEnd(dithermeaning_t *dm, f32 (*q)[3])
{ // q gets the meanings, least squares by a Cholesky solve over the indexes
	const u8 *pal = CURWORLDPAL; // used. Without q it only frees dm.
	s32 idx[256], n = 0;
	f64 b[256][3], *a = NULL;
	if(!q)
		goto done;
	for(s32 c = 0; c < 256; c++){
		for(s32 k = 0; k < 3; k++)
			q[c][k] = pal[c * 3 + k];
		if(dm->u[c] > 0)
			idx[n++] = c;
	}
	if(!n || !(a = malloc(n * n * sizeof(*a))))
		goto done;
	for(s32 i = 0; i < n; i++){
		s32 c = idx[i];
		f64 hold = MEANING_HOLD * dm->u[c];
		for(s32 j = 0; j < n; j++)
			a[i * n + j] = dm->m[c][idx[j]] + (i == j) * hold;
		for(s32 k = 0; k < 3; k++)
			b[i][k] = dm->v[c][k] + hold * pal[c * 3 + k];
	}
	for(s32 j = 0; j < n; j++){ // a = L L^T, L in a's lower half
		f64 d = a[j * n + j];
		for(s32 k = 0; k < j; k++)
			d -= a[j * n + k] * a[j * n + k];
		a[j * n + j] = d = sqrt(q_max(d, 1e-12));
		for(s32 i = j + 1; i < n; i++){
			f64 t = a[i * n + j];
			for(s32 k = 0; k < j; k++)
				t -= a[i * n + k] * a[j * n + k];
			a[i * n + j] = t / d;
		}
	}
	for(s32 i = 0; i < n; i++) // L y = b
		for(s32 k = 0; k < 3; k++){
			for(s32 j = 0; j < i; j++)
				b[i][k] -= a[i * n + j] * b[j][k];
			b[i][k] /= a[i * n + i];
		}
	for(s32 i = n - 1; i >= 0; i--) // L^T x = y
		for(s32 k = 0; k < 3; k++){
			for(s32 j = i + 1; j < n; j++)
				b[i][k] -= a[j * n + i] * b[j][k];
			b[i][k] /= a[i * n + i];
			q[idx[i]][k] = CLAMP(0, b[i][k], 255);
		}
done:
	free(a);
	free(dm);
}

static u8 Dither_Nearby(ditherpal_t *d, u8 n0, s32 cell, const u8 *rgb)
{ // the colour in n0's spread nearest to rgb, which is in nn_lut's cell
	u32 key = cell << 8 | n0, slot = key * 2654435761u >> 20;
	if(d->seenkey[slot] == key)
		return d->seenc[slot];
	const u8 *pal = CURWORLDPAL, *list = d->nearby[n0];
	s32 best = 0x7fffffff;
	u8 c = n0;
	for(s32 j = 0; j < d->nnear[n0]; j++){
		s32 dist = Pal_Dist2(pal + list[j] * 3, rgb);
		if(dist < best){
			best = dist;
			c = list[j];
		}
	}
	d->seenkey[slot] = key;
	return d->seenc[slot] = c;
}

void R_DitherRGBA(u8 *pic, s32 w, s32 h, f32 dither, dithermeaning_t *dm)
{ // errors in 1/256ths, dither out of 256. dm, if any, gets the tiles.
	const u8 *pal = CURWORLDPAL;
	s32 wt = CLAMP(0, dither, 1) * 256;
	s32 tw = q_max(1, w / MEANING_TILES), ncol = (w + tw - 1) / tw;
	ditherpal_t *d = calloc(1, sizeof(ditherpal_t));
	s32 *err = calloc((w + 2) * 6, sizeof(s32)); // this row's, next row's
	u32 *cnt = dm ? calloc(ncol * 256, sizeof(u32)) : NULL; // by tile
	f32 *sum = dm ? calloc(ncol * 3, sizeof(f32)) : NULL;
	if(!d || !err || (dm && (!cnt || !sum)))
		Sys_Error("R_DitherRGBA: out of memory");
	NN_PrepareLUT();
	memset(d->seenkey, 0xff, sizeof(d->seenkey)); // no cell is that high
	for(s32 a = 0; a < 256; a++)
		for(s32 b = 0; b < 256; b++){
			if(Pal_Dist2(pal + a*3, pal + b*3)
					<= DITHER_SPREAD * DITHER_SPREAD)
				d->nearby[a][d->nnear[a]++] = b;
			else d->outside[a][b >> 5] |= 1u << (b & 31);
		}
	for(s32 n = -DITHER_WARMUP; n < h; n++){
		s32 y = n < 0 ? n + DITHER_WARMUP : n;
		s32 *cur = err + (y & 1) * (w + 2) * 3;
		s32 *next = err + (~y & 1) * (w + 2) * 3;
		s32 dir = y & 1 ? -1 : 1; // row 0 goes right, so indexes
		s32 ahead[3] = {0, 0, 0}; // the error for the next pixel along
		for(s32 i = 0; i < w; i++){ // never land on RGBA not read yet
			s32 x = dir > 0 ? i : w - 1 - i;
			u8 src[3], q[3];
			memcpy(src, pic + (y * w + x) * 4, 3);
			u8 n0 = nn_lut[NN_LUTCELL(src)];
			s32 want[3];
			for(s32 k = 0; k < 3; k++){
				want[k] = CLAMP(-DITHER_MAXERR * 256,
					cur[(x + 1) * 3 + k] + ahead[k],
					DITHER_MAXERR * 256) + (src[k] << 8);
				q[k] = CLAMP(0, (want[k] + 128) >> 8, 255);
			}
			s32 cell = NN_LUTCELL(q);
			u8 c = nn_lut[cell];
			if(d->outside[n0][c >> 5] >> (c & 31) & 1)
				c = Dither_Nearby(d, n0, cell, q);
			if(n >= 0)
				pic[y * w + x] = c;
			if(n >= 0 && cnt){
				cnt[x / tw * 256 + c]++;
				for(s32 k = 0; k < 3; k++)
					sum[x / tw * 3 + k] += src[k];
			}
			s32 *below = next + (x + 1) * 3;
			for(s32 k = 0; k < 3; k++){ // the last gets the rest
				s32 r = want[k] - (pal[c*3 + k] << 8);
				r = r * wt >> 8;
				s32 r7 = r * 7 >> 4, r3 = r * 3 >> 4;
				s32 r5 = r * 5 >> 4;
				ahead[k] = r7;
				below[k - dir * 3] += r3;
				below[k] += r5;
				below[k + dir * 3] += r - r7 - r3 - r5;
			}
		}
		s32 *end = next + (dir > 0 ? w : 1) * 3;
		s32 *start = next + (dir > 0 ? 1 : w) * 3;
		for(s32 k = 0; k < 3; k++){ // what went past the ends, to them
			end[k] += ahead[k] + end[k + dir * 3];
			start[k] += start[k - dir * 3];
		}
		memset(cur, 0, (w + 2) * 3 * sizeof(s32));
		if(n >= 0 && cnt && ((y + 1) % tw == 0 || y == h - 1))
			Dither_AddTiles(dm, cnt, sum, ncol);
	}
	free(err);
	free(d);
	free(cnt);
	free(sum);
}

static f32 Lin(f32 v)
{ // 0-255 sRGB to 0-1 linear light
	v /= 255;
	return v > 0.04045f ? powf((v + 0.055f) / 1.055f, 2.4f) : v / 12.92f;
}

static f32 Form(f32 (*w)[3], const f32 *x, const f32 *y)
{ return Dot3(x, w[0]) * y[0] + Dot3(x, w[1]) * y[1] + Dot3(x, w[2]) * y[2]; }

static void Lab_Jacobian(const f32 *lin, f32 (*j)[3])
{ // how L, a and b change with linear r, g and b at lin
	static const f32 m[3][3] = {{0.4124f, 0.3576f, 0.1805f},
		{0.2126f, 0.7152f, 0.0722f}, {0.0193f, 0.1192f, 0.9505f}};
	static const f32 white[3] = {0.95047f, 1, 1.08883f};
	f32 d[3];
	for(s32 i = 0; i < 3; i++){ // lab_f's slope at x, y and z, by r g b
		f32 t = Dot3(m[i], lin) / white[i];
		d[i] = (t > 0.008856f ? 1 / (3 * cbrtf(t * t)) : 903.3f / 116)
			/ white[i];
	}
	for(s32 k = 0; k < 3; k++){
		j[0][k] = 116 * d[1] * m[1][k];
		j[1][k] = 500 * (d[0] * m[0][k] - d[1] * m[1][k]);
		j[2][k] = 200 * (d[1] * m[1][k] - d[2] * m[2][k]);
	}
}

// The palette color, or two and how many 64ths are the second, whose dither
// looks nearest rgb: the mean of their light nearest it in Lab, as Lab goes
// near rgb. Fog chains measure in sRGB, which on a grey can miss its hue by
// several Lab steps. The dither's grain counts too, its steps of lightness
// fully and of color a tenth: grain this fine shows little color, and it's
// colors far apart in hue that mix into a grey of the right tint.
#define NEAR_CANDS 8 // the colors nearest rgb it tries
#define NEAR_GRAIN 0.5f // what grain costs against error
#define NEAR_GRAINCOLOR 0.1f // of that for color, not lightness
void R_NearestDither(const f32 *rgb, u8 *lo, u8 *hi, s32 *t)
{
	const u8 *pal = CURWORLDPAL;
	s32 cand[NEAR_CANDS], n = 0;
	f32 dist[NEAR_CANDS], lin[3], j[3][3], w[3][3], e[NEAR_CANDS][3];
	init_color_conv();
	for(s32 i = 0; i < 256; i++){ // the nearest, nearest first
		f32 v[3] = {pal[i * 3] - rgb[0], pal[i * 3 + 1] - rgb[1],
			pal[i * 3 + 2] - rgb[2]};
		f32 dd = Dot3(v, v);
		if(n == NEAR_CANDS && dd >= dist[NEAR_CANDS - 1])
			continue;
		s32 k = n < NEAR_CANDS ? n++ : NEAR_CANDS - 1;
		for(; k > 0 && dist[k - 1] > dd; k--){
			dist[k] = dist[k - 1];
			cand[k] = cand[k - 1];
		}
		dist[k] = dd;
		cand[k] = i;
	}
	for(s32 k = 0; k < 3; k++)
		lin[k] = Lin(rgb[k]);
	Lab_Jacobian(lin, j);
	for(s32 a = 0; a < 3; a++) // Lab distance squared, as a form on light
		for(s32 b = 0; b < 3; b++)
			w[a][b] = j[0][a] * j[0][b] + j[1][a] * j[1][b]
				+ j[2][a] * j[2][b];
	f32 best = FLT_MAX;
	for(s32 i = 0; i < n; i++){
		for(s32 k = 0; k < 3; k++)
			e[i][k] = gamma_lut[pal[cand[i] * 3 + k]] - lin[k];
		f32 cost = Form(w, e[i], e[i]);
		if(cost < best){
			best = cost;
			*lo = *hi = cand[i];
			*t = 0;
		}
	}
	for(s32 a = 0; a < n; a++)
	for(s32 b = a + 1; b < n; b++){
		if(Pal_Dist2(pal + cand[a] * 3, pal + cand[b] * 3)
				> DITHER_SPREAD * DITHER_SPREAD)
			continue;
		f32 ab[3], err[3];
		for(s32 k = 0; k < 3; k++)
			ab[k] = e[b][k] - e[a][k];
		f32 len = Form(w, ab, ab);
		if(!(len > 0))
			continue;
		s32 tt = CLAMP(0, -Form(w, e[a], ab) / len * 64 + 0.5f, 64);
		if(tt <= 0 || tt >= 64)
			continue;
		f32 r = tt / 64.0f;
		for(s32 k = 0; k < 3; k++)
			err[k] = e[a][k] + r * ab[k];
		const lab_t *pa = &lab_palette[cand[a]];
		const lab_t *pb = &lab_palette[cand[b]];
		f32 dl = pb->l - pa->l, da = pb->a - pa->a, db = pb->b - pa->b;
		f32 cost = Form(w, err, err) + NEAR_GRAIN * r * (1 - r)
			* (dl * dl + NEAR_GRAINCOLOR * (da * da + db * db));
		if(cost < best){
			best = cost;
			*lo = cand[a];
			*hi = cand[b];
			*t = tt;
		}
	}
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
	Fog_PaletteChanged(); // r_fogdither 0 mixes through lit_lut
	D_FlushCaches(0); // surfaces get lit through lit_lut
}

void R_WorldPalChanged()
{ // what was built from the old colors gets built again when it's needed
	color_conv_initialized = 0;
	refresh_palette++;
	R_LabMixPalChanged(0);
	Draw_PaletteChanged();
	Sky_PaletteChanged();
	Fog_PaletteChanged();
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
