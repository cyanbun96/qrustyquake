#include "quakedef.h"

#define FOG_FACTOR_LUT_DOWNSAMPLE 4
#define FOG_FACTOR_LUT_SIZE (65536 / FOG_FACTOR_LUT_DOWNSAMPLE)

static f32 fog_red; // CyanBun96: we store the actual RGB values in these,
static f32 fog_green; // but they get quantized to a single index in the color
static f32 fog_blue; // palette before use, stored in fog_pal_index
static f32 *randarr = 0;
static u32 lfsr = 0x1337; // non-zero seed
static f32 *fog_factor_lut = 0;
static f32 old_r_fogscale = -1234;
static f32 old_fog_density = -1234;
static f32 old_noisebias = -1234;
static f32 *zadjarr = 0;
static s32 zadjarr_w = -1;
static s32 zadjarr_h = -1;

void Fog_SetPalIndex(SDL_UNUSED cvar_t *cvar)
{
	fog_pal_index = rgbtoi_lab(fog_red * 255.0f*r_fogbrightness.value,
				   fog_green*255.0f*r_fogbrightness.value,
				   fog_blue *255.0f*r_fogbrightness.value);
}

void Fog_Update(f32 d, f32 r, f32 g, f32 b) {
	if(r_lockfog.value){
		d = r_lockfogd.value;
		r = r_lockfogr.value;
		g = r_lockfogg.value;
		b = r_lockfogb.value;
	}
	if     (d < 0.0f) d = 0.0f;
	if     (r < 0.0f) r = 0.0f;
	else if(r > 1.0f) r = 1.0f;
	if     (g < 0.0f) g = 0.0f;
	else if(g > 1.0f) g = 1.0f;
	if     (b < 0.0f) b = 0.0f;
	else if(b > 1.0f) b = 1.0f;
	fog_density = d;
	fog_red = r;
	fog_green = g;
	fog_blue = b;
	Fog_SetPalIndex(0);
	vid.recalc_refdef = 1;
	fog_initialized = 0;
}

void Fog_FogCommand_f () // yanked from Quakespasm, mostly
{ // TODO time transitions
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
		d = Q_atof(Cmd_Argv(1));
		r = fog_red;
		g = fog_green;
		b = fog_blue;
		break;
	case 3: //TEST
		d = Q_atof(Cmd_Argv(1));
		r = fog_red;
		g = fog_green;
		b = fog_blue;
		break;
	case 4:
		d = fog_density;
		r = Q_atof(Cmd_Argv(1));
		g = Q_atof(Cmd_Argv(2));
		b = Q_atof(Cmd_Argv(3));
		break;
	case 5:
		d = Q_atof(Cmd_Argv(1));
		r = Q_atof(Cmd_Argv(2));
		g = Q_atof(Cmd_Argv(3));
		b = Q_atof(Cmd_Argv(4));
		break;
	}
	Fog_Update(d, r, g, b);
}

void Fog_ParseWorldspawn () // from Quakespasm
{ // called at map load
	s8 key[128], value[4096];
	if(r_lockfog.value){
		fog_density = r_lockfogd.value;
		fog_red = r_lockfogr.value;
		fog_green = r_lockfogg.value;
		fog_blue = r_lockfogb.value;
		Fog_SetPalIndex(0);
		return;
	} else fog_density = fog_red = fog_green = fog_blue = 0;
	const s8 *data = COM_Parse(cl.worldmodel->entities);
	if(!data || com_token[0] != '{') return; // error
	while(1){
		if(!data) return; // error
		if(com_token[0] == '}') break; // end of worldspawn
		if(com_token[0] == '_')q_strlcpy(key, com_token+1, sizeof(key));
		else q_strlcpy(key, com_token, sizeof(key));
		while(key[0] && key[strlen(key)-1] == ' ') // no trailing spaces
			key[strlen(key)-1] = 0;
		data = COM_ParseEx(data, CPE_ALLOWTRUNC);
		if(!data) return; // error
		q_strlcpy(value, com_token, sizeof(value));
		if(!strcmp("fog", key))
			sscanf(value, "%f %f %f %f", &fog_density, &fog_red,
					&fog_green, &fog_blue);
	}
	Fog_SetPalIndex(0);
}

static void R_InitDepthCorrection()
{ // the array is normalized to [0.8,1] with 0.8 in center and 1 in corners
	if(zadjarr_w==scr_vrect.width&&zadjarr_h==scr_vrect.height) return;
	zadjarr_w = scr_vrect.width; zadjarr_h = scr_vrect.height;
	zadjarr = realloc(zadjarr, sizeof(f32)*zadjarr_w*zadjarr_h);
	s32 xmid = scr_vrect.width/2;
	s32 ymid = scr_vrect.height/2;
	f32 max = xmid*xmid + ymid*ymid;
	for(s32 i = 0; i < scr_vrect.height; ++i){
	for(s32 j = 0; j < scr_vrect.width; ++j){
		f32 dist = (j-xmid)*(j-xmid) + (i-ymid)*(i-ymid);
		zadjarr[i*scr_vrect.width + j] = (dist/max)*(dist/max)/5+0.8;
	}}
}

u32 lfsr_random() {lfsr^=lfsr>>7; lfsr^=lfsr<<9; lfsr^=lfsr>>13; return lfsr;}

void R_InitFog(f32 noisebias)
{
	if(!randarr)randarr = malloc(RANDARR_SIZE*sizeof(f32));
	for(s32 i = 0; i < RANDARR_SIZE; ++i) { // fog bias array
		randarr[i] = (lfsr_random() & 0xFFFF) / 65535.0f; // == [0,1]
		randarr[i] *= noisebias;
	}
	old_noisebias = noisebias;
	if(!fog_lut_built) R_BuildColorMixLUT(0);
	Fog_SetPalIndex(0);
	fog_initialized = 1;
}

void R_InitFogLUT() {
	if(!fog_factor_lut)
		fog_factor_lut = malloc(FOG_FACTOR_LUT_SIZE*sizeof(f32));
	f32 fog_scale = 32.0f * r_fogscale.value;
	f32 fog_scale_inv = 1.0f / fog_scale;
	f32 fd_sq = (1.0f - fog_density) * (1.0f - fog_density);
	for (s32 i = 0; i < FOG_FACTOR_LUT_SIZE; ++i) {
		s32 z = (i * FOG_FACTOR_LUT_DOWNSAMPLE) - 32768;
		s32 z_clamped = (z < fog_scale) ? fog_scale : z;
		f32 z_norm = z_clamped * fog_scale_inv;
		fog_factor_lut[i] = expf(-fd_sq * (z_norm * z_norm));
	}
	old_r_fogscale = r_fogscale.value;
	old_fog_density = fog_density;
}

static inline f32 compute_fog_lut(s32 z) {
	z = z + 32768; // convert signed to unsigned
	if (z < 0) z = 0;
	if (z >= 65536) z = 65535;
	s32 idx = z / FOG_FACTOR_LUT_DOWNSAMPLE;
	if (idx >= FOG_FACTOR_LUT_SIZE) idx = FOG_FACTOR_LUT_SIZE - 1;
	return fog_factor_lut[idx];
}

void R_DrawFog(){
	if(!fog_density || r_nofog.value) return;
	f32 noisebias = 10 * r_fognoise.value;
	if(!fog_initialized || old_noisebias!=noisebias)R_InitFog(noisebias);
	if(old_r_fogscale!=r_fogscale.value||old_fog_density!=fog_density)
		R_InitFogLUT();
	if(r_fogdepthcorrection.value) R_InitDepthCorrection();
	sb_updates = 0; // draw sbar over fog
	s32 j = 0;
	u8 *pdest = screen->pixels;
	s32 area = scr_vrect.width * scr_vrect.height;
	switch((s32)r_fogstyle.value + (r_fogdepthcorrection.value?10:0)){
	case 0: // noisy
		for(s32 i = 0; i < area; ++i){
			s32 bias = randarr[(area-j++)%RANDARR_SIZE];
			f32 ffactor = compute_fog_lut(d_pzbuffer[i]+bias)*r_fogfactor.value;
			if((lfsr_random() & 0xFFFF) / 65535.0f < ffactor)
				pdest[i] = fog_pal_index;
		} break;
	case 2: // dither + noise
		for(s32 i = 0; i < area; ++i){
			s32 bias = randarr[(area-j++)%RANDARR_SIZE];
			f32 ffactor = compute_fog_lut(d_pzbuffer[i]+bias)*r_fogfactor.value;
			ffactor -= bias*0.05;
			if((lfsr_random() & 0xFFFF) / 65535.0f < ffactor)
				pdest[i] = fog_pal_index;
		} break;
	case 1: // dither
		for(s32 i = 0; i < area; ++i){
			s32 bias = randarr[(area-j++)%RANDARR_SIZE];
			f32 ffactor = compute_fog_lut(d_pzbuffer[i]+bias)*r_fogfactor.value;
			if(D_Dither(&pdest[i], ffactor)) pdest[i]=fog_pal_index;
		} break;
	default: case 3: // mix
		for(s32 i = 0; i < area; ++i){
			s32 bias = randarr[(area-j++)%RANDARR_SIZE];
			f32 ffactor = compute_fog_lut(d_pzbuffer[i]+bias)*r_fogfactor.value;
			u8 pix = pdest[i];
			s32 lut_idx = (s32)(ffactor * FOG_LUT_LEVELS);
			pdest[i] = color_mix_lut[pix][fog_pal_index][lut_idx];
		} break;
	case 10: // depth-corrected, noisy
		for(s32 i = 0; i < area; ++i){
			s32 bias = randarr[(area-j++)%RANDARR_SIZE];
			f32 ffactor = compute_fog_lut(d_pzbuffer[i]+bias)*r_fogfactor.value*zadjarr[i];
			if((lfsr_random() & 0xFFFF) / 65535.0f < ffactor)
				pdest[i] = fog_pal_index;
		} break;
	case 12: // depth-corrected, dither + noise
		for(s32 i = 0; i < area; ++i){
			s32 bias = randarr[(area-j++)%RANDARR_SIZE];
			f32 ffactor = compute_fog_lut(d_pzbuffer[i]+bias)*r_fogfactor.value*zadjarr[i];
			ffactor -= bias*0.05;
			if((lfsr_random() & 0xFFFF) / 65535.0f < ffactor)
				pdest[i] = fog_pal_index;
		} break;
	case 11: // depth-corrected, dither
		for(s32 i = 0; i < area; ++i){
			s32 bias = randarr[(area-j++)%RANDARR_SIZE];
			f32 ffactor = compute_fog_lut(d_pzbuffer[i]+bias)*r_fogfactor.value*zadjarr[i];
			if(D_Dither(&pdest[i], ffactor)) pdest[i]=fog_pal_index;
		} break;
	case 13: // depth-corrected, mix
		for(s32 i = 0; i < area; ++i){
			s32 bias = randarr[(area-j++)%RANDARR_SIZE];
			f32 ffactor = compute_fog_lut(d_pzbuffer[i]+bias)*r_fogfactor.value*zadjarr[i];
			u8 pix = pdest[i];
			s32 lut_idx = (s32)(ffactor * FOG_LUT_LEVELS);
			pdest[i] = color_mix_lut[pix][fog_pal_index][lut_idx];
		} break;
	}
}
