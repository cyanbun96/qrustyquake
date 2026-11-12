#include "quakedef.h"

// Fog as Ironwail draws it: 1 - exp(-(density * 0.85 / 64 * distance)^2) of
// the fog color over every pixel, by its distance from the eye. With
// r_fogdepthcorrection 0 it's QuakeSpasm's, density / 64 along the view axis.
// In the palette each color's way to the fog color is a chain of palette
// colors, and a pixel between two links gets one of them by an 8x8 Bayer
// threshold, so that the mean of the dither follows the way.

// fog_level is by the top 16 bits of a float, 2^31 / distance: the z-buffer
// holds 2^31 / z, and for distance from the eye that's times fog_cos, the
// cosine of the angle between the pixel's direction and the view axis.
#define FOG_GRAY 0.3f // the fog color of a map that only sets a density
#define FOG_LEVELS 64 // fog amounts in fog_pair, none to all
#define FOG_DEPTHS (161 << 7) // floats up to 2^32, where (u32)-1 rounds to
#define FOG_SPANS 8 // stretches of a way, each gives the colors nearest it
#define FOG_SPANLINKS 2 // as links
#define FOG_LINKS (1 + FOG_SPANS * FOG_SPANLINKS + 2) // and c and the end pair
#define FOG_ENDLEAD (1 - 4.0f / (FOG_LEVELS - 1)) // where the end pair starts
#define FOG_NOISE 0.2f // what the dither's noise costs in a chain
#define FOG_PACE 0.004f // what a link going faster or slower than it costs

typedef struct { f32 d, r, g, b; } fog_t;
// A link: palette color c placed at s on a color's way to the fog, e off it.
// The end of the way has c -1 and the mean of the end pair as its rgb.
typedef struct { f32 s, e[3], rgb[3]; s32 c; } foglink_t;
static fog_t fog_to = {0, FOG_GRAY, FOG_GRAY, FOG_GRAY}; // where a fade ends
static fog_t fog_from = {0, FOG_GRAY, FOG_GRAY, FOG_GRAY}; // and began
static f64 fade_time, fade_done;
static u8 fog_rgb[3]; // the fog color fog_pair was built for
static bool fog_built; // fog_pair is for fog_rgb and the palette
static f32 fog_level_d = -1; // the density fog_level was built for
static bool fog_radial; // fog_level is by distance from the eye
static s32 fog_skylevel; // the fog amount of sky pixels
static u8 fog_level[FOG_DEPTHS]; // fog amount by depth
static u32 fog_pair[256][FOG_LEVELS]; // low | high << 8 | threshold << 16
static f32 fog_palf[3][256]; // the palette's channels, for Fog_Links
static f32 *fog_cos; // by pixel of the view, for distance from the eye
// Sky gets map_skyfog instead of fog by distance, but it keeps its real z to
// hide what's behind it. So this remembers the z it wrote: a pixel shows sky,
// or a translucent surface over it, where the z-buffer still holds that value.
// It's kept xor fog_skykey, a key for the frame, so that what an earlier frame
// left there is a 2^-32 chance.
static u32 *fog_skyz = 0;
static s32 fog_skyz_size = 0;
static u32 fog_skykey;
static u32 fog_skyrow[MAXHEIGHT]; // frame this row last got sky on
static u32 lfsr = 0x1337; // lfsr_random's state, non-zero

u32 lfsr_random()
{
	lfsr ^= lfsr >> 7;
	lfsr ^= lfsr << 9;
	lfsr ^= lfsr >> 13;
	return lfsr;
}

// Where the fog is now, fades included
static fog_t Fog_Faded()
{
	if(fade_done <= cl.time)
		return fog_to;
	f32 f = (fade_done - cl.time) / fade_time;
	fog_t m = {f * fog_from.d + (1 - f) * fog_to.d,
		f * fog_from.r + (1 - f) * fog_to.r,
		f * fog_from.g + (1 - f) * fog_to.g,
		f * fog_from.b + (1 - f) * fog_to.b};
	return m;
}

void Fog_Update(f32 d, f32 r, f32 g, f32 b, f32 time) // from Ironwail
{ // a fade starts from where the fog is, even halfway through another one
	fog_t to = {d, r, g, b};
	if(time > 0)
		fog_from = Fog_Faded();
	fog_to = to;
	fade_time = time;
	fade_done = cl.time + time;
}

void Fog_FogCommand_f() // from Ironwail
{
	f32 d = fog_to.d, r = fog_to.r, g = fog_to.g, b = fog_to.b, t = 0;
	switch(Cmd_Argc()){
	default:
	case 1:
		Con_Printf("usage:\n");
		Con_Printf("   fog <density>\n");
		Con_Printf("   fog <red> <green> <blue>\n");
		Con_Printf("   fog <density> <red> <green> <blue>\n");
		Con_Printf("current values:\n");
		Con_Printf("   \"density\" is \"%f\"\n", fog_to.d);
		Con_Printf("   \"red\" is \"%f\"\n", fog_to.r);
		Con_Printf("   \"green\" is \"%f\"\n", fog_to.g);
		Con_Printf("   \"blue\" is \"%f\"\n", fog_to.b);
		return;
	case 3: // and a fade time
		t = atof(Cmd_Argv(2));
		// fallthrough
	case 2:
		d = atof(Cmd_Argv(1));
		break;
	case 4:
		r = atof(Cmd_Argv(1));
		g = atof(Cmd_Argv(2));
		b = atof(Cmd_Argv(3));
		break;
	case 6: // mg1's trigger_fog sends its delay as the fade time
		t = atof(Cmd_Argv(5));
		// fallthrough
	case 5:
		d = atof(Cmd_Argv(1));
		r = atof(Cmd_Argv(2));
		g = atof(Cmd_Argv(3));
		b = atof(Cmd_Argv(4));
		break;
	}
	Fog_Update(q_max(d, 0), CLAMP(0, r, 1), CLAMP(0, g, 1), CLAMP(0, b, 1),
			t);
}

void Fog_ParseWorldspawn() // from Quakespasm
{ // called at map load
	static const fog_t none = {0, FOG_GRAY, FOG_GRAY, FOG_GRAY};
	c8 key[128], value[4096];
	fog_to = fog_from = none;
	fade_time = fade_done = 0;
	const c8 *data = COM_Parse(cl.worldmodel->entities);
	if(!data || com_token[0] != '{') return; // error
	while(COM_ParseEntityPair(&data, key, value))
		if(!strcmp("fog", key))
			sscanf(value, "%f %f %f %f", &fog_to.d, &fog_to.r,
					&fog_to.g, &fog_to.b);
}

void Fog_PaletteChanged()
{ fog_built = 0; }

static f32 Fog_Dot(const f32 *a, const f32 *b)
{ return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

static u32 Fog_Pair(s32 lo, s32 hi, s32 t)
{
	if(t >= 64)
		lo = hi;
	if(t <= 0 || t >= 64)
		return lo | lo << 8;
	return lo | hi << 8 | t << 16;
}

// The two palette colors and threshold whose dither comes nearest the fog
// color, or one color if that's nearer. Every color's last level, so that all
// fog looks the same once it's thick enough to hide what's in it.
static u32 Fog_EndPair(const f32 *fog)
{
	const u8 *pal = CURWORLDPAL;
	s32 cand[8], n = 0;
	f32 dist[8];
	for(s32 i = 0; i < 256; i++){ // the 8 nearest, nearest first
		f32 v[3] = {pal[i * 3] - fog[0], pal[i * 3 + 1] - fog[1],
			pal[i * 3 + 2] - fog[2]};
		f32 d = Fog_Dot(v, v);
		if(n == 8 && d >= dist[7])
			continue;
		s32 j = n < 8 ? n++ : 7;
		for(; j > 0 && dist[j - 1] > d; j--){
			dist[j] = dist[j - 1];
			cand[j] = cand[j - 1];
		}
		dist[j] = d;
		cand[j] = i;
	}
	f32 best = dist[0];
	s32 lo = cand[0], hi = lo, t = 0;
	for(s32 i = 0; i < n; i++)
	for(s32 j = i + 1; j < n; j++){
		f32 pf[3], pq[3];
		for(s32 k = 0; k < 3; k++){
			pf[k] = fog[k] - pal[cand[i] * 3 + k];
			pq[k] = pal[cand[j] * 3 + k] - pal[cand[i] * 3 + k];
		}
		f32 len = Fog_Dot(pq, pq);
		if(!len)
			continue;
		s32 tt = CLAMP(0, Fog_Dot(pf, pq) / len * 64 + 0.5f, 64);
		f32 r = tt / 64.0f, e[3];
		for(s32 k = 0; k < 3; k++)
			e[k] = r * pq[k] - pf[k];
		f32 cost = Fog_Dot(e, e) + FOG_NOISE * r * (1 - r) * len;
		if(cost < best && tt > 0 && tt < 64){
			best = cost;
			lo = cand[i];
			hi = cand[j];
			t = tt;
		}
	}
	return Fog_Pair(lo, hi, t);
}

// Palette color i as a link on the way from color a to the fog, ab long,
// placed where it's nearest the way but not past most
static foglink_t Fog_Link(s32 i, const f32 *a, const f32 *ab, f32 len,
		f32 most)
{
	const u8 *pal = CURWORLDPAL;
	f32 v[3];
	for(s32 k = 0; k < 3; k++)
		v[k] = pal[i * 3 + k] - a[k];
	f32 s = q_min(Fog_Dot(v, ab) / len, most);
	foglink_t l = {s, {0, 0, 0}, {0, 0, 0}, i};
	for(s32 k = 0; k < 3; k++){
		l.e[k] = v[k] - l.s * ab[k];
		l.rgb[k] = pal[i * 3 + k];
	}
	return l;
}

// The links a chain from color c to the fog can use, a is c and ab the way:
// in each of FOG_SPANS stretches of the way the palette colors nearest it, and
// the colors of the end pair, a few levels from the end so that they can lead
// into it. c first, the rest in the order they lie along the way. Returns how
// many.
static s32 Fog_Links(foglink_t *l, s32 c, const f32 *a, const f32 *ab,
		u32 end)
{
	f32 len = Fog_Dot(ab, ab), s[256], d[256];
	f32 dist[FOG_SPANS][FOG_SPANLINKS];
	s32 cand[FOG_SPANS][FOG_SPANLINKS], ncand[FOG_SPANS] = {0}, n = 1;
	foglink_t start = {0, {0, 0, 0}, {a[0], a[1], a[2]}, c};
	l[0] = start;
	if(!(len > 0))
		return 1;
	for(s32 i = 0; i < 256; i++){ // where each color lies, a loop that
		f32 v[3]; // vectorizes, Fog_Link places the few that are kept
		for(s32 k = 0; k < 3; k++)
			v[k] = fog_palf[k][i] - a[k];
		s[i] = Fog_Dot(v, ab) / len;
		f32 t = q_min(s[i], 1);
		for(s32 k = 0; k < 3; k++)
			v[k] -= t * ab[k];
		d[i] = Fog_Dot(v, v);
	}
	for(s32 i = 0; i < 256; i++){
		if(i == c || s[i] <= 0) // c starts the way
			continue;
		if(i == (s32)(end & 255) || i == (s32)(end >> 8 & 255)){
			l[n++] = Fog_Link(i, a, ab, len, FOG_ENDLEAD);
			continue;
		}
		s32 b = q_min(s[i] * FOG_SPANS, FOG_SPANS - 1), *m = &ncand[b];
		if(*m == FOG_SPANLINKS && d[i] >= dist[b][*m - 1])
			continue;
		s32 j = *m < FOG_SPANLINKS ? (*m)++ : *m - 1;
		for(; j > 0 && dist[b][j - 1] > d[i]; j--){
			dist[b][j] = dist[b][j - 1];
			cand[b][j] = cand[b][j - 1];
		}
		dist[b][j] = d[i];
		cand[b][j] = i;
	}
	for(s32 b = 0; b < FOG_SPANS; b++)
		for(s32 j = 0; j < ncand[b]; j++)
			l[n++] = Fog_Link(cand[b][j], a, ab, len, 1);
	for(s32 i = 2; i < n; i++){
		foglink_t t = l[i];
		s32 j = i;
		for(; j > 1 && l[j - 1].s > t.s; j--)
			l[j] = l[j - 1];
		l[j] = t;
	}
	return n;
}

// What the dither from link p to link q costs: how far its mean strays from
// the way, its noise and how far its pace is from the way's
static f32 Fog_LinkCost(const foglink_t *p, const foglink_t *q, const f32 *ab)
{
	f32 ds = q->s - p->s, pq[3], g[3];
	for(s32 k = 0; k < 3; k++){
		pq[k] = q->rgb[k] - p->rgb[k];
		g[k] = pq[k] - ab[k] * ds;
	}
	f32 e = Fog_Dot(p->e, p->e) + Fog_Dot(p->e, q->e)
		+ Fog_Dot(q->e, q->e);
	return ds * (e / 3 + FOG_NOISE * Fog_Dot(pq, pq) / 6)
		+ FOG_PACE * Fog_Dot(g, g) / ds;
}

// The same for keeping link p's color from its place to the last level,
// where the end pair takes over at once
static f32 Fog_HoldCost(const foglink_t *p, const foglink_t *end,
		const f32 *ab)
{
	f32 h = 1 - p->s, ab2 = Fog_Dot(ab, ab), jump[3];
	for(s32 k = 0; k < 3; k++)
		jump[k] = end->rgb[k] - p->rgb[k];
	return h * Fog_Dot(p->e, p->e) - h * h * Fog_Dot(p->e, ab)
		+ h * h * h / 3 * ab2 + FOG_PACE * h * ab2
		+ FOG_PACE * (FOG_LEVELS - 1) * Fog_Dot(jump, jump);
}

// Color c's fog levels: the chain of links that costs the least, then the
// two links around each level and how far it is between them. A chain can
// reach the end pair from one of its colors, the threshold going from that
// color's to the pair's, or hold its last color until the end pair.
static void Fog_BuildChain(u32 *pair, s32 c, const f32 *fog, u32 end)
{
	const u8 *pal = CURWORLDPAL;
	s32 elo = end & 255, ehi = end >> 8 & 255, et = end >> 16;
	foglink_t l[FOG_LINKS], e = {1, {0, 0, 0}, {0, 0, 0}, -1};
	f32 a[3], ab[3], cost[FOG_LINKS];
	s32 from[FOG_LINKS], chain[FOG_LINKS];
	for(s32 k = 0; k < 3; k++){
		a[k] = pal[c * 3 + k];
		ab[k] = fog[k] - a[k];
		e.rgb[k] = pal[elo * 3 + k]
			+ et / 64.0f * (pal[ehi * 3 + k] - pal[elo * 3 + k]);
		e.e[k] = e.rgb[k] - fog[k];
	}
	s32 n = Fog_Links(l, c, a, ab, end), last = 0;
	bool toend = 0;
	f32 best = FLT_MAX;
	cost[0] = 0;
	for(s32 i = 0; i < n; i++){
		if(i){
			cost[i] = FLT_MAX;
			for(s32 j = 0; j < i; j++){
				if(l[j].s >= l[i].s)
					continue;
				f32 t = cost[j]
					+ Fog_LinkCost(&l[j], &l[i], ab);
				if(t < cost[i]){
					cost[i] = t;
					from[i] = j;
				}
			}
		}
		f32 t = cost[i] + Fog_HoldCost(&l[i], &e, ab);
		if(t < best){
			best = t;
			last = i;
			toend = 0;
		}
		if((l[i].c != elo && l[i].c != ehi) || l[i].s >= 1)
			continue;
		t = cost[i] + Fog_LinkCost(&l[i], &e, ab);
		if(t < best){
			best = t;
			last = i;
			toend = 1;
		}
	}
	s32 len = 0;
	for(s32 i = last; i; i = from[i])
		chain[len++] = i;
	chain[len] = 0; // c, chain runs backwards
	for(s32 lev = 0; lev < FOG_LEVELS - 1; lev++){
		f32 s = (f32)lev / (FOG_LEVELS - 1);
		while(len && l[chain[len - 1]].s <= s)
			len--;
		const foglink_t *p = &l[chain[len]];
		if(len){
			const foglink_t *q = &l[chain[len - 1]];
			s32 t = (s - p->s) / (q->s - p->s) * 64 + 0.5f;
			pair[lev] = Fog_Pair(p->c, q->c, t);
		}else if(toend){
			f32 u = (s - p->s) / (1 - p->s); // p's color to end
			s32 t = p->c == elo ? u * et + 0.5f
				: 64 - u * (64 - et) + 0.5f;
			pair[lev] = Fog_Pair(elo, ehi, t);
		}else
			pair[lev] = p->c;
	}
	pair[FOG_LEVELS - 1] = end;
}

static void Fog_BuildPairs()
{
	f32 fog[3] = {fog_rgb[0], fog_rgb[1], fog_rgb[2]};
	for(s32 i = 0; i < 256 * 3; i++)
		fog_palf[i % 3][i / 3] = CURWORLDPAL[i];
	u32 end = Fog_EndPair(fog);
	for(s32 c = 0; c < 256; c++)
		Fog_BuildChain(fog_pair[c], c, fog, end);
	fog_built = 1;
}

static void Fog_BuildLevels(f32 density, bool radial)
{
	f64 k = density * (radial ? 0.85 : 1) / 64;
	for(s32 i = 0; i < FOG_DEPTHS; i++){
		u32 b = (u32)i << 16 | 0x8000; // the middle of the floats here
		f32 f;
		memcpy(&f, &b, sizeof(f));
		f64 kr = k * 2147483648.0 / f; // f is 2^31 / distance
		// past these ends it's levels 0 and 63 without exp
		f64 a = kr < 0.06 ? 0 : kr > 3 ? 1 : 1 - exp(-kr * kr);
		fog_level[i] = a * (FOG_LEVELS - 1) + 0.5;
	}
	fog_level_d = density;
	fog_radial = radial;
}

// worldspawn colors can be past 1, and r_fogbrightness takes them further
static u8 Fog_Channel(f32 c)
{ return CLAMP(0, Q_rint(c * 255 * r_fogbrightness.value), 255); }

void R_FogSetupFrame()
{ // after R_ThreadsBeginFrame and before anything is drawn
	fog_t f = Fog_Faded();
	if(r_lockfog.value){
		f.d = r_lockfogd.value;
		f.r = r_lockfogr.value;
		f.g = r_lockfogg.value;
		f.b = r_lockfogb.value;
	}
	fog_density = r_nofog.value || !(f.d > 0) ? 0 : f.d;
	if(!fog_density)
		return;
	u8 rgb[3] = {Fog_Channel(f.r), Fog_Channel(f.g), Fog_Channel(f.b)};
	bool radial = r_fogdepthcorrection.value != 0;
	if(!fog_built || memcmp(rgb, fog_rgb, 3)){
		memcpy(fog_rgb, rgb, 3);
		Fog_BuildPairs();
	}
	if(fog_level_d != fog_density || fog_radial != radial)
		Fog_BuildLevels(fog_density, radial);
	fog_skylevel = CLAMP(0, map_skyfog, 1) * (FOG_LEVELS - 1) + 0.5f;
	fog_skykey = (u32)r_framecount * 0x9e3779b9; // the same once in 2^32
}

// w pixels of a view row from screen column x. bay is the row's thresholds,
// cos the pixels' fog_cos. Inlined so that radial gets its own loop.
Q_SIL void R_FogPixels(u8 *pdest, const s32 *zb, s32 x, s32 w,
		const u8 *bay, const f32 *cos, s32 radial)
{ // in locals, the stores to pdest could change any global
	const u8 *level = fog_level;
	const u32 *pair = &fog_pair[0][0];
	for(s32 i = 0; i < w; i++){
		f32 f = (u32)zb[i]; // the void's negative z comes out near
		if(radial)
			f *= cos[i];
		u32 b;
		memcpy(&b, &f, sizeof(b));
		u32 e = pair[pdest[i] * FOG_LEVELS + level[b >> 16]];
		pdest[i] = bay[(x + i) & 7] < e >> 16 ? e >> 8 : e;
	}
}

static void R_FogRow(u8 *pdest, const s32 *zb, s32 x, s32 w, const u8 *bay,
		const f32 *cos)
{
	if(fog_radial) R_FogPixels(pdest, zb, x, w, bay, cos, 1);
	else R_FogPixels(pdest, zb, x, w, bay, NULL, 0);
}

// w sky pixels from screen column x
static void R_FogSky(u8 *pdest, s32 x, s32 w, const u8 *bay)
{
	const u32 *pair = &fog_pair[0][fog_skylevel];
	for(s32 i = 0; i < w; i++){
		u32 e = pair[pdest[i] * FOG_LEVELS];
		pdest[i] = bay[(x + i) & 7] < e >> 16 ? e >> 8 : e;
	}
}

void R_FogSkyAlloc()
{ // before R_FogMarkSky, on the main thread
	if(!fog_density) return;
	s32 size = vid.width * vid.height;
	if(fog_skyz_size != size){
		R_ThreadsSync();
		free(fog_skyz);
		fog_skyz = calloc(size, sizeof(*fog_skyz));
		if(!fog_skyz) Sys_Error("R_FogMarkSky: out of memory");
		fog_skyz_size = size;
	}
}

// Fogs the sky as it's drawn, so that translucent surfaces, which don't write
// z, go over fogged sky and R_DrawFog leaves them be
void R_FogMarkSky(const espan_t *span)
{ // has to come after the span's z and colors were written
	if(!fog_density) return;
	u32 key = fog_skykey;
	for(; span; span = span->pnext){ // locals, a u32 store could change
		s32 ofs = d_zwidth * span->v + span->u; // span->count
		s32 n = span->count;
		const s32 *zb = d_pzbuffer + ofs;
		u32 *sky = fog_skyz + ofs;
		for(s32 i = 0; i < n; i++)
			sky[i] = (u32)zb[i] ^ key;
		fog_skyrow[span->v] = (u32)r_framecount;
		if(fog_skylevel)
			R_FogSky(d_viewbuffer + screenwidth * span->v
				+ span->u, span->u, n, d_bayer[span->v & 7]);
	}
}

// R_FogRow for a row with sky on it, the runs of sky are fogged already
static void R_FogRowAroundSky(u8 *pdest, const s32 *zb, s32 x, s32 w,
		const u8 *bay, const f32 *cos, const u32 *sky)
{
	u32 key = fog_skykey;
	for(s32 i = 0; i < w;){
		bool issky = sky[i] == ((u32)zb[i] ^ key);
		s32 end = i + 1;
		while(end < w && (sky[end] == ((u32)zb[end] ^ key)) == issky)
			end++;
		if(!issky)
			R_FogRow(pdest + i, zb + i, x + i, end - i, bay,
					cos ? cos + i : NULL);
		i = end;
	}
}

// Rows y0 to y1 - 1 of the view, before any warp
static void R_FogRows(s32 y0, s32 y1)
{
	const vrect_t *r = &r_refdef.vrect;
	y0 = q_max(y0, r->y);
	y1 = q_min(y1, r->y + r->height);
	for(s32 v = y0; v < y1; v++){
		s32 ofs = v * vid.width + r->x; // the view and z buffers are
		u8 *pdest = d_viewbuffer + ofs; // vid.width wide, the view
		const s32 *zb = d_pzbuffer + ofs; // can be narrower
		const u8 *bay = d_bayer[v & 7];
		const f32 *cos = fog_radial ? fog_cos + (v - r->y) * r->width
			: NULL;
		if(fog_skyrow[v] == (u32)r_framecount)
			R_FogRowAroundSky(pdest, zb, r->x, r->width, bay, cos,
					fog_skyz + ofs);
		else
			R_FogRow(pdest, zb, r->x, r->width, bay, cos);
	}
}

static void R_FogCmd(SDL_UNUSED void *arg)
{ R_FogRows(d_bandy0, d_bandy1); }

// fog_cos for the view, when it has changed
static void R_FogCosines()
{
	static f32 xc, yc, xs, ys;
	static s32 w, h;
	const vrect_t *r = &r_refdef.vrect;
	if(fog_cos && w == r->width && h == r->height && xc == xcenter
			&& yc == ycenter && xs == xscale && ys == yscale)
		return;
	w = r->width; h = r->height;
	xc = xcenter; yc = ycenter; xs = xscale; ys = yscale;
	free(fog_cos);
	if(!(fog_cos = malloc(sizeof(*fog_cos) * w * h)))
		Sys_Error("R_DrawFog: out of memory");
	for(s32 v = 0; v < h; v++){
		f32 ty = (r->y + v - ycenter) / yscale;
		for(s32 u = 0; u < w; u++){
			f32 tx = (r->x + u - xcenter) / xscale;
			fog_cos[v * w + u] = 1 / sqrtf(1 + tx * tx + ty * ty);
		}
	}
}

void R_DrawFog()
{
	if(!fog_density) return;
	const vrect_t *r = &r_refdef.vrect;
	if(fog_radial)
		R_FogCosines();
	sb_updates = 0; // draw sbar over fog
	if(r_threaded)
		R_RowsAdd(r->y, r->y + r->height, R_FogCmd, NULL);
	else
		R_FogRows(r->y, r->y + r->height);
}
