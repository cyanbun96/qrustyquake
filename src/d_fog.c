#include "quakedef.h"

#define FOG_FACTOR_LUT_DOWNSAMPLE 4
#define FOG_FACTOR_LUT_SIZE (65536 / FOG_FACTOR_LUT_DOWNSAMPLE)
#define FOG_MAX_DENSITY 0.9999
#define FOG_MAX_MIX 0.99f // stays a valid level when zadj's 1.005 gets on top

static f32 fog_red; // CyanBun96: we store the actual RGB values in these,
static f32 fog_green; // but they get quantized to a single index in the color
static f32 fog_blue; // palette before use, stored in fog_pal_index
static s32 *randarr = 0; // the bias is only ever used as a whole number
static f32 *fog_factor_lut = 0;
static f32 *fog_mix_lut = 0; // fog_factor_lut * r_fogfactor, for mix styles
static f32 old_r_fogfactor = -1234;
static f32 old_r_fogscale = -1234;
static f32 old_fog_density = -1234;
static f32 old_noisebias = -1234;
static f32 *zadjarr = 0;
// Sky gets map_skyfog instead of fog by distance, but it keeps its real z to
// hide what's behind it. So this remembers the z it wrote: a pixel still shows
// sky where the z-buffer holds that value. Low 32 bits z, high 32 the frame.
static u64 *fog_skyz = 0;
static s32 fog_skyz_size = 0;
static u32 fog_skyrow[MAXHEIGHT]; // frame this row last got sky on

static u8 Fog_Channel(f32 c) // worldspawn colours can be past 1 and so can
{ return CLAMP(0, c * 255.0f * r_fogbrightness.value, 255); } // the brightness

void Fog_SetPalIndex(SDL_UNUSED cvar_t *cvar)
{
	fog_pal_index = rgbtoi_lab(Fog_Channel(fog_red),
			Fog_Channel(fog_green), Fog_Channel(fog_blue));
}

void Fog_Update(f32 d, f32 r, f32 g, f32 b){
	if(r_lockfog.value){
		d = r_lockfogd.value;
		r = r_lockfogr.value;
		g = r_lockfogg.value;
		b = r_lockfogb.value;
	}
	fog_density = CLAMP(0,d,FOG_MAX_DENSITY);
	fog_red = CLAMP(0,r,1);
	fog_green = CLAMP(0,g,1);
	fog_blue = CLAMP(0,b,1);
	Fog_SetPalIndex(0);
	vid.recalc_refdef = 1;
	fog_initialized = 0;
}

void Fog_FogCommand_f() // yanked from Quakespasm, mostly
{
	f32 d, r, g, b;
	switch(Cmd_Argc()){
	default:
	case 1:
		Con_Printf("usage:\n");
		Con_Printf("   fog <density>\n");
		Con_Printf("   fog <red> <green> <blue>\n");
		Con_Printf("   fog <density> <red> <green> <blue>\n");
		Con_Printf("current values:\n");
		Con_Printf("   \"density\" is \"%f\"\n", fog_density);
		Con_Printf("   \"red\" is \"%f\"\n", fog_red);
		Con_Printf("   \"green\" is \"%f\"\n", fog_green);
		Con_Printf("   \"blue\" is \"%f\"\n", fog_blue);
		return;
	case 2:
	case 3: // QuakeSpasm's fade time is ignored
		d = atof(Cmd_Argv(1));
		r = fog_red;
		g = fog_green;
		b = fog_blue;
		break;
	case 4:
		d = fog_density;
		r = atof(Cmd_Argv(1));
		g = atof(Cmd_Argv(2));
		b = atof(Cmd_Argv(3));
		break;
	case 5:
		d = atof(Cmd_Argv(1));
		r = atof(Cmd_Argv(2));
		g = atof(Cmd_Argv(3));
		b = atof(Cmd_Argv(4));
		break;
	}
	Fog_Update(d, r, g, b);
}

void Fog_ParseWorldspawn() // from Quakespasm
{ // called at map load
	c8 key[128], value[4096];
	if(r_lockfog.value){
		fog_density = CLAMP(0, r_lockfogd.value, FOG_MAX_DENSITY);
		fog_red = r_lockfogr.value;
		fog_green = r_lockfogg.value;
		fog_blue = r_lockfogb.value;
		Fog_SetPalIndex(0);
		return;
	}else fog_density = fog_red = fog_green = fog_blue = 0;
	const c8 *data = COM_Parse(cl.worldmodel->entities);
	if(!data || com_token[0] != '{') return; // error
	while(COM_ParseEntityPair(&data, key, value))
		if(!strcmp("fog", key))
			sscanf(value, "%f %f %f %f", &fog_density, &fog_red,
					&fog_green, &fog_blue);
	fog_density = CLAMP(0, fog_density, FOG_MAX_DENSITY);
	Fog_SetPalIndex(0);
}

static void R_InitDepthCorrection()
{ // the array is normalized to [0.8,1] with 0.8 in center and 1 in corners
	static s32 zadjarr_w = -1;
	static s32 zadjarr_h = -1;
	if(zadjarr_w==scr_vrect.width&&zadjarr_h==scr_vrect.height) return;
	zadjarr_w = scr_vrect.width; zadjarr_h = scr_vrect.height;
	zadjarr = realloc(zadjarr, sizeof(f32)*zadjarr_w*zadjarr_h);
	s32 xmid = scr_vrect.width/2;
	s32 ymid = scr_vrect.height/2;
	f32 max = xmid*xmid + ymid*ymid;
	for(s32 i = 0; i < scr_vrect.height; ++i){
	for(s32 j = 0; j < scr_vrect.width; ++j){
		f32 dist = (j-xmid)*(j-xmid) + (i-ymid)*(i-ymid);
		f32 noise = ((f32)((s32)lfsr_random() % 64))/(64.0*100.0)-0.005;
		zadjarr[i*scr_vrect.width+j]=(dist/max)*(dist/max)/5+0.8+noise;
	}}
}

u32 lfsr_random()
{
	static u32 lfsr = 0x1337; // non-zero seed
	lfsr^=lfsr>>7;
	lfsr^=lfsr<<9;
	lfsr^=lfsr>>13;
	return lfsr;
}

void R_InitFog(f32 noisebias)
{
	if(!randarr)randarr = malloc(RANDARR_SIZE*sizeof(s32));
	for(s32 i = 0; i < RANDARR_SIZE; ++i){ // fog bias array
		f32 r = (lfsr_random() & 0xFFFF) / 65535.0f; // == [0,1]
		randarr[i] = r * noisebias;
	}
	old_noisebias = noisebias;
	if(!fog_lut_built) R_BuildColorMixLUT(0);
	Fog_SetPalIndex(0);
	fog_initialized = 1;
}

void R_InitFogLUT()
{
	if(!fog_factor_lut)
		fog_factor_lut = malloc(FOG_FACTOR_LUT_SIZE*sizeof(f32));
	if(!fog_mix_lut)
		fog_mix_lut = malloc(FOG_FACTOR_LUT_SIZE*sizeof(f32));
	f32 fog_scale = 32.0f * r_fogscale.value;
	f32 fog_scale_inv = 1.0f / fog_scale;
	f32 fd_sq = (1.0f - fog_density) * (1.0f - fog_density);
	for(s32 i = 0; i < FOG_FACTOR_LUT_SIZE; ++i){
		s32 z = (i * FOG_FACTOR_LUT_DOWNSAMPLE) - 32768;
		s32 z_clamped = (z < fog_scale) ? fog_scale : z;
		f32 z_norm = z_clamped * fog_scale_inv;
		fog_factor_lut[i] = expf(-fd_sq * (z_norm * z_norm));
		// the tail of expf goes through denormals, which are very slow
		// to multiply, and is far below anything visible anyway
		if(fog_factor_lut[i] < 1e-6f)
			fog_factor_lut[i] = 0;
		fog_mix_lut[i] = fog_factor_lut[i] * r_fogfactor.value;
		if(fog_mix_lut[i] < 1e-6f)
			fog_mix_lut[i] = 0;
		if(fog_mix_lut[i] > FOG_MAX_MIX)
			fog_mix_lut[i] = FOG_MAX_MIX;
	}
	old_r_fogfactor = r_fogfactor.value;
	old_r_fogscale = r_fogscale.value;
	old_fog_density = fog_density;
}

static inline f32 compute_fog_lut(const f32 *lut, s32 z)
{
	z = z + 32768; // convert signed to unsigned
	if(z < 0) z = 0;
	if(z >= 65536) z = 65535;
	s32 idx = z / FOG_FACTOR_LUT_DOWNSAMPLE;
	if(idx >= FOG_FACTOR_LUT_SIZE) idx = FOG_FACTOR_LUT_SIZE - 1;
	return lut[idx];
}

static void R_FogMix(u8 *pdest, const s32 *zb, s32 w, s32 n, const f32 *zadj)
{
	const s32 *bias = randarr;
	const f32 *lut = fog_mix_lut;
	const u8 *mix = &color_mix_lut[0][fog_pal_index][0];
	const s32 stride = sizeof(color_mix_lut[0]);
	for(s32 i = 0; i < w; ++i){
		s32 z = (zb[i] >> 16) + bias[(n - i) & (RANDARR_SIZE - 1)]
			+ 32768;
		if(z < 0) z = 0;
		if(z > 65535) z = 65535;
		f32 ffactor = lut[z / FOG_FACTOR_LUT_DOWNSAMPLE];
		if(zadj)
			ffactor *= zadj[i];
		s32 lut_idx = (s32)(ffactor * FOG_LUT_LEVELS);
		pdest[i] = mix[pdest[i] * stride + lut_idx];
	}
}

// The fog styles but mix: 0 noisy, 1 dither, 2 dither + noise. Inlined so
// that each style, with and without depth correction, gets its own loop.
Q_SIL void R_FogPixels(u8 *pdest, const s32 *zb, s32 w, s32 n,
		const f32 *zadj, s32 style)
{
	// in locals, the stores to pdest could change any global
	const s32 *rnd = randarr;
	const f32 *lut = fog_factor_lut, fogfactor = r_fogfactor.value;
	const u8 fog = fog_pal_index;
	for(s32 i = 0; i < w; ++i){
		s32 bias = rnd[(n-i)&(RANDARR_SIZE-1)];
		f32 ffactor = compute_fog_lut(lut, (zb[i]>>16) + bias)
			* fogfactor;
		if(zadj)
			ffactor *= zadj[i];
		if(style == 2)
			ffactor -= bias*0.05;
		if(style == 1 ? D_Dither(&pdest[i], ffactor) :
				(lfsr_random() & 0xFFFF) / 65535.0f < ffactor)
			pdest[i] = fog;
	}
}

// w pixels of one view row, n is how many of the view's pixels are left
static void R_FogRow(u8 *pdest, const s32 *zb, s32 w, s32 n, const f32 *zadj)
{
	switch((s32)r_fogstyle.value + (zadj ? 10 : 0)){
	case 0: R_FogPixels(pdest, zb, w, n, NULL, 0); break;
	case 1: R_FogPixels(pdest, zb, w, n, NULL, 1); break;
	case 2: R_FogPixels(pdest, zb, w, n, NULL, 2); break;
	case 10: R_FogPixels(pdest, zb, w, n, zadj, 0); break;
	case 11: R_FogPixels(pdest, zb, w, n, zadj, 1); break;
	case 12: R_FogPixels(pdest, zb, w, n, zadj, 2); break;
	default: R_FogMix(pdest, zb, w, n, zadj); break; // 3: mix
	}
}

void R_FogMarkSky(const espan_t *span)
{ // has to come after the span's z was written
	if(!fog_density || r_nofog.value) return;
	s32 size = vid.width * vid.height;
	if(fog_skyz_size != size){
		free(fog_skyz);
		fog_skyz = calloc(size, sizeof(*fog_skyz));
		if(!fog_skyz) Sys_Error("R_FogMarkSky: out of memory");
		fog_skyz_size = size;
	}
	u64 stamp = (u64)r_framecount << 32;
	for(; span; span = span->pnext){
		s32 ofs = d_zwidth * span->v + span->u;
		for(s32 i = 0; i < span->count; i++)
			fog_skyz[ofs + i] = stamp | (u32)d_pzbuffer[ofs + i];
		fog_skyrow[span->v] = (u32)r_framecount;
	}
}

// R_FogRow for a row with sky on it, on the runs of pixels between the sky
static void R_FogRowAroundSky(u8 *pdest, const s32 *zb, s32 w, s32 n,
		const f32 *zadj, const u64 *sky)
{
	u64 stamp = (u64)r_framecount << 32;
	for(s32 x = 0; x < w;){
		if(sky[x] == (stamp | (u32)zb[x])){
			x++;
			continue;
		}
		s32 end = x + 1;
		while(end < w && sky[end] != (stamp | (u32)zb[end]))
			end++;
		R_FogRow(pdest + x, zb + x, end - x, n - x,
				zadj ? zadj + x : NULL);
		x = end;
	}
}

void R_DrawFog()
{
	if(!fog_density || r_nofog.value) return;
	f32 noisebias = 10 * r_fognoise.value;
	if(!fog_initialized || old_noisebias!=noisebias)R_InitFog(noisebias);
	if(old_r_fogscale!=r_fogscale.value || old_fog_density!=fog_density
			|| old_r_fogfactor != r_fogfactor.value)
		R_InitFogLUT();
	if(r_fogdepthcorrection.value) R_InitDepthCorrection();
	sb_updates = 0; // draw sbar over fog
	s32 w = scr_vrect.width;
	s32 ofs = scr_vrect.y * vid.width + scr_vrect.x;
	u8 *pdest = vid.buffer + ofs; // the screen and z buffer are vid.width
	const s32 *zb = d_pzbuffer + ofs; // wide, the view can be narrower
	const f32 *zadj = r_fogdepthcorrection.value ? zadjarr : NULL;
	s32 v = scr_vrect.y;
	for(s32 n = w * scr_vrect.height; n > 0; n -= w, v++){
		if(fog_skyrow[v] == (u32)r_framecount)
			R_FogRowAroundSky(pdest, zb, w, n, zadj,
					fog_skyz + (zb - d_pzbuffer));
		else
			R_FogRow(pdest, zb, w, n, zadj);
		pdest += vid.width;
		zb += vid.width;
		if(zadj) zadj += w;
	}
}
