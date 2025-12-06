// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// r_surf.c: surface-related refresh code
#include "quakedef.h"

static s32 lightleft, blocksize, sourcetstep;
static s32 lightleft_g;
static s32 lightleft_b;
static s32 lightright, lightleftstep, lightrightstep, blockdivshift;
static s32 lightright_g, lightleftstep_g, lightrightstep_g;
static s32 lightright_b, lightleftstep_b, lightrightstep_b;
static unsigned blockdivmask;
static void *prowdestbase;
static u8 *pbasesource;
static s32 surfrowbytes; // used by ASM files
static unsigned *r_lightptr;
static unsigned *r_lightptr_g;
static unsigned *r_lightptr_b;
static s32 r_stepback;
static s32 r_lightwidth;
static s32 r_numhblocks, r_numvblocks;
static u8 *r_source, *r_sourcemax;
static u32 blocklights[18 * 18];
static u32 blocklights_g[18 * 18];
static u32 blocklights_b[18 * 18];

void R_DrawSurfaceBlock(s32 miplvl);
void R_DrawSurfaceBlockRGBMono(s32 miplvl);
void R_DrawSurfaceBlockRGB(s32 miplvl);

void R_AddDynamicLights()
{
	//johnfitz -- lit support via lordhavoc
	f32 cred, cgreen, cblue, brightness;
	//johnfitz
	msurface_t *surf = r_drawsurf.surf;
	s32 smax = (surf->extents[0] >> 4) + 1;
	s32 tmax = (surf->extents[1] >> 4) + 1;
	mtexinfo_t *tex = surf->texinfo;
	for (s32 lnum = 0; lnum < MAX_DLIGHTS; lnum++) {
		if (!(surf->dlightbits & (1 << lnum)))
			continue; // not lit by this light
		f32 rad = cl_dlights[lnum].radius;
		f32 dist = DotProduct(cl_dlights[lnum].origin,
				surf->plane->normal) - surf->plane->dist;
		rad -= fabs(dist);
		f32 minlight = cl_dlights[lnum].minlight;
		if (rad < minlight)
			continue;
		minlight = rad - minlight;
		vec3_t impact, local;
		for (s32 i = 0; i < 3; i++)
			impact[i] = cl_dlights[lnum].origin[i] -
				surf->plane->normal[i] * dist;
		local[0] = DotProduct(impact, tex->vecs[0]) + tex->vecs[0][3];
		local[1] = DotProduct(impact, tex->vecs[1]) + tex->vecs[1][3];
		local[0] -= surf->texturemins[0];
		local[1] -= surf->texturemins[1];
		//johnfitz -- lit support via lordhavoc
		u32 *bl = blocklights;
		u32 *bl_g = blocklights_g;
		u32 *bl_b = blocklights_b;
                cred = cl_dlights[lnum].color[0] * 256.0f;
                cgreen = cl_dlights[lnum].color[1] * 256.0f;
                cblue = cl_dlights[lnum].color[2] * 256.0f;
                //johnfitz
		for (s32 t = 0; t < tmax; t++) {
			s32 td = local[1] - t * 16;
			if (td < 0)
				td = -td;
			for (s32 s = 0; s < smax; s++) {
				s32 sd = local[0] - s * 16;
				if (sd < 0)
					sd = -sd;
				dist = sd > td ? sd + td / 2 : td + sd / 2;
				if (dist < minlight) {
					//johnfitz -- lit support via lordhavoc
					brightness = rad - dist;
					bl[0] += (s32) (brightness * cred);
					bl_g[0] += (s32) (brightness * cgreen);
					bl_b[0] += (s32) (brightness * cblue);
				}
				bl += 1;
				bl_g += 1;
				bl_b += 1;
				//johnfitz
			}
		}
	}
}


void R_BuildLightMap()
{ // Combine and scale multiple lightmaps into the 8.8 format in blocklights
	msurface_t *surf = r_drawsurf.surf;
	s32 smax = surf->extents[0] / 16 + 1;
	s32 tmax = surf->extents[1] / 16 + 1;
	s32 size = smax * tmax;
	u8 *lightmap = surf->samples;
	if (r_fullbright.value || !cl.worldmodel->lightdata) {
		memset(blocklights,   0, size * sizeof(*blocklights));
		memset(blocklights_g, 0, size * sizeof(*blocklights_g));
		memset(blocklights_b, 0, size * sizeof(*blocklights_b));
		return;
	}
	for (s32 i = 0; i < size; i++) { // clear to ambient
		blocklights[i] = r_refdef.ambientlight << 8;
		blocklights_g[i] = r_refdef.ambientlight << 8;
		blocklights_b[i] = r_refdef.ambientlight << 8;
	}
	if (lightmap) // add all the lightmaps
		for (s32 maps = 0; maps < MAXLIGHTMAPS && surf->styles[maps] != 255; maps++) {
			u32 scale = r_drawsurf.lightadj[maps];	// 8.8 fraction         
			u32 *bl = blocklights;
			u32 *bl_g = blocklights_g;
			u32 *bl_b = blocklights_b;
			for (s32 i = 0; i < size; i++) {
				*bl++ += *lightmap++ * scale;
				*bl_g++ += *lightmap++ * scale;
				*bl_b++ += *lightmap++ * scale;
			}
		}
	if ((u32)surf->dlightframe == r_framecount) //add all the dynamic lights
		R_AddDynamicLights();
	for (s32 i = 0; i < size; i++) { // bound, invert, and shift
		s32 t = (255 * 256 - (s32)blocklights[i]) >> (8 - VID_CBITS);
		if (t < 64)
			t = 64;
		blocklights[i] = t;
	}
	for (s32 i = 0; i < size; i++) { // Green
		s32 t = (255 * 256 - (s32)blocklights_g[i]) >> (8 - VID_CBITS);
		if (t < 64)
			t = 64;
		blocklights_g[i] = t;
	}
	for (s32 i = 0; i < size; i++) { // Blue
		s32 t = (255 * 256 - (s32)blocklights_b[i]) >> (8 - VID_CBITS);
		if (t < 64)
			t = 64;
		blocklights_b[i] = t;
	}
}

texture_t *R_TextureAnimation(texture_t *base)
{ // Returns the proper texture for a given time and base texture
	if (currententity->frame && base->alternate_anims)
		base = base->alternate_anims;
	if (!base->anim_total)
		return base;
	s32 reletive = (s32)(cl.time * 10) % base->anim_total;
	s32 count = 0;
	while (base->anim_min > reletive || base->anim_max <= reletive) {
		base = base->anim_next;
		if (!base)
			Sys_Error("R_TextureAnimation: broken cycle");
		if (++count > 100)
			Sys_Error("R_TextureAnimation: infinite cycle");
	}
	return base;
}

s32 R_DrawSurface()
{
	R_BuildLightMap(); // calculate the lightings
	surfrowbytes = r_drawsurf.rowbytes;
	texture_t *mt = r_drawsurf.texture;
	r_source = (u8 *) mt + mt->offsets[r_drawsurf.surfmip];
	// the fractional light values should range from 0 to (VID_GRADES - 1)
	// << 16 from a source range of 0 - 255
	s32 texwidth = mt->width >> r_drawsurf.surfmip;
	blocksize = 16 >> r_drawsurf.surfmip;
	blockdivshift = 4 - r_drawsurf.surfmip;
	blockdivmask = (1 << blockdivshift) - 1;
	r_lightwidth = (r_drawsurf.surf->extents[0] >> 4) + 1;
	r_numhblocks = r_drawsurf.surfwidth >> blockdivshift;
	r_numvblocks = r_drawsurf.surfheight >> blockdivshift;
	// TODO: only needs to be set when there is a display settings change
	s32 horzblockstep = blocksize;
	s32 smax = mt->width >> r_drawsurf.surfmip;
	s32 twidth = texwidth;
	s32 tmax = mt->height >> r_drawsurf.surfmip;
	if (!smax || !tmax) {
		Con_DPrintf("Division by 0 in R_DrawSurface\n");
		return 1;
	}
	sourcetstep = texwidth;
	r_stepback = tmax * twidth;
	r_sourcemax = r_source + (tmax * smax);
	s32 soffset = r_drawsurf.surf->texturemins[0];
	s32 basetoffset = r_drawsurf.surf->texturemins[1];
	// << 16 components are to guarantee positive values for %
	soffset = ((soffset >> r_drawsurf.surfmip) + (smax << 16)) % smax;
	u8 *basetptr = &r_source[((((basetoffset>>r_drawsurf.surfmip)
				+ (tmax << 16)) % tmax) * twidth)];
	u8 *pcolumndest = r_drawsurf.surfdat;
	for (s32 u = 0; u < r_numhblocks; u++) {
		r_lightptr = blocklights + u;
		r_lightptr_g = blocklights_g + u;
		r_lightptr_b = blocklights_b + u;
		prowdestbase = pcolumndest;
		pbasesource = basetptr + soffset;
		if (!lit_loaded) R_DrawSurfaceBlock(r_drawsurf.surfmip);
		else if (!r_rgblighting.value) R_DrawSurfaceBlockRGBMono(r_drawsurf.surfmip);
		else R_DrawSurfaceBlockRGB(r_drawsurf.surfmip);
		soffset = soffset + blocksize;
		if (soffset >= smax)
			soffset = 0;
		pcolumndest += horzblockstep;
	}
	return 0;
}

void R_DrawSurfaceBlock(s32 miplvl)
{
	u8 *psource = pbasesource;
	u8 *prowdest = prowdestbase;
	for (s32 v = 0; v < r_numvblocks; v++) {
		// FIXME: make these locals?
		// FIXME: use delta rather than both right and left, like ASM?
		lightleft = r_lightptr[0];
		lightright = r_lightptr[1];
		r_lightptr += r_lightwidth;
		lightleftstep = (r_lightptr[0] - lightleft) >> (4-miplvl);
		lightrightstep = (r_lightptr[1] - lightright) >> (4-miplvl);
		for (s32 i = 0; i < 1 << (4-miplvl); i++) {
			s32 lighttemp = lightleft - lightright;
			s32 lightstep = lighttemp >> (4-miplvl);
			s32 light = lightright;
			for (s32 b = (1 << (4-miplvl)) - 1; b >= 0; b--) {
				u8 pix = psource[b];
				prowdest[b] = CURWORLDCMAP[(light&0xFF00)+pix];
				light += lightstep;
			}
			psource += sourcetstep;
			lightright += lightrightstep;
			lightleft += lightleftstep;
			prowdest += surfrowbytes;
		}
		if (psource >= r_sourcemax)
			psource -= r_stepback;
	}
}

void R_DrawSurfaceBlockRGBMono(s32 miplvl)
{
	u8 *psource = pbasesource;
	u8 *prowdest = prowdestbase;
	for (s32 v = 0; v < r_numvblocks; v++) {
		// FIXME: make these locals?
		// FIXME: use delta rather than both right and left, like ASM?
		lightleft = r_lightptr[0]; // Red
		lightright = r_lightptr[1];
		r_lightptr += r_lightwidth;
		lightleftstep = (r_lightptr[0] - lightleft) >> (4-miplvl);
		lightrightstep = (r_lightptr[1] - lightright) >> (4-miplvl);
		lightleft_g = r_lightptr_g[0]; // Green
		lightright_g = r_lightptr_g[1];
		r_lightptr_g += r_lightwidth;
		lightleftstep_g = (r_lightptr_g[0] - lightleft_g) >> (4-miplvl);
		lightrightstep_g = (r_lightptr_g[1] - lightright_g) >> (4-miplvl);
		lightleft_b = r_lightptr_b[0]; // Blue
		lightright_b = r_lightptr_b[1];
		r_lightptr_b += r_lightwidth;
		lightleftstep_b = (r_lightptr_b[0] - lightleft_b) >> (4-miplvl);
		lightrightstep_b = (r_lightptr_b[1] - lightright_b) >> (4-miplvl);
		for (s32 i = 0; i < 1 << (4-miplvl); i++) {
			s32 lighttemp = lightleft - lightright; // Red
			s32 lightstep = lighttemp >> (4-miplvl);
			s32 light = lightright;
			s32 light_g = lightright_g;
			s32 light_b = lightright_b;
			for (s32 b = (1 << (4-miplvl)) - 1; b >= 0; b--) {
				u8 pix = psource[b];
				s32 lightavg = ((light&0xFF00) + (light_g&0xFF00) + (light_b&0xFF00))/3;
				prowdest[b] = CURWORLDCMAP[(lightavg & 0xFF00) + pix];
				light += lightstep;
			}
			psource += sourcetstep;
			lightright += lightrightstep;
			lightleft += lightleftstep;
			lightright_g += lightrightstep_g;
			lightleft_g += lightleftstep_g;
			lightright_b += lightrightstep_b;
			lightleft_b += lightleftstep_b;
			prowdest += surfrowbytes;
		}
		if (psource >= r_sourcemax)
			psource -= r_stepback;
	}
}

void R_DrawSurfaceBlockRGB(s32 miplvl)
{
	if (!lit_lut_initialized) R_BuildLitLUT();
	u8 *psource = pbasesource;
	u8 *prowdest = prowdestbase;
	for (s32 v = 0; v < r_numvblocks; v++) {
		// FIXME: make these locals?
		// FIXME: use delta rather than both right and left, like ASM?
		lightleft = r_lightptr[0]; // Red
		lightright = r_lightptr[1];
		r_lightptr += r_lightwidth;
		lightleftstep = (r_lightptr[0] - lightleft) >> (4-miplvl);
		lightrightstep = (r_lightptr[1] - lightright) >> (4-miplvl);
		lightleft_g = r_lightptr_g[0]; // Green
		lightright_g = r_lightptr_g[1];
		r_lightptr_g += r_lightwidth;
		lightleftstep_g = (r_lightptr_g[0] - lightleft_g) >> (4-miplvl);
		lightrightstep_g = (r_lightptr_g[1] - lightright_g) >> (4-miplvl);
		lightleft_b = r_lightptr_b[0]; // Blue
		lightright_b = r_lightptr_b[1];
		r_lightptr_b += r_lightwidth;
		lightleftstep_b = (r_lightptr_b[0] - lightleft_b) >> (4-miplvl);
		lightrightstep_b = (r_lightptr_b[1] - lightright_b) >> (4-miplvl);
		for (s32 i = 0; i < 1 << (4-miplvl); i++) {
			s32 lighttemp = lightleft - lightright; // Red
			s32 lightstep = lighttemp >> (4-miplvl);
			s32 light = lightright;
			s32 lighttemp_g = lightleft_g - lightright_g; // Green
			s32 lightstep_g = lighttemp_g >> (4-miplvl);
			s32 light_g = lightright_g;
			s32 lighttemp_b = lightleft_b - lightright_b; // Blue
			s32 lightstep_b = lighttemp_b >> (4-miplvl);
			s32 light_b = lightright_b;
			for (s32 b = (1 << (4-miplvl)) - 1; !lmonly && b >= 0; b--) {
				u8 tex = psource[b];
				u8 ir = CURWORLDCMAP[(light   & 0xFF00) + tex];
				u8 ig = CURWORLDCMAP[(light_g & 0xFF00) + tex];
				u8 ib = CURWORLDCMAP[(light_b & 0xFF00) + tex];
				u8 r_ = CURWORLDPAL[ir * 3 + 0];
				u8 g_ = CURWORLDPAL[ig * 3 + 1];
				u8 b_ = CURWORLDPAL[ib * 3 + 2];
				prowdest[b] = lit_lut[QUANT(r_)][QUANT(g_)][QUANT(b_)];
				light += lightstep;
				light_g += lightstep_g;
				light_b += lightstep_b;
			}
			for (s32 b = (1 << (4-miplvl)) - 1; lmonly && b >= 0; b--) {
				u8 ir = CURWORLDCMAP[(light   & 0xFF00)+15];
				u8 ig = CURWORLDCMAP[(light_g & 0xFF00)+15];
				u8 ib = CURWORLDCMAP[(light_b & 0xFF00)+15];
				u8 r_ = CURWORLDPAL[ir * 3 + 0];
				u8 g_ = CURWORLDPAL[ig * 3 + 1];
				u8 b_ = CURWORLDPAL[ib * 3 + 2];
				prowdest[b] = lit_lut[QUANT(r_)][QUANT(g_)][QUANT(b_)];
				light += lightstep;
				light_g += lightstep_g;
				light_b += lightstep_b;
			}
			psource += sourcetstep;
			lightright += lightrightstep;
			lightleft += lightleftstep;
			lightright_g += lightrightstep_g;
			lightleft_g += lightleftstep_g;
			lightright_b += lightrightstep_b;
			lightleft_b += lightleftstep_b;
			prowdest += surfrowbytes;
		}
		if (psource >= r_sourcemax)
			psource -= r_stepback;
	}
}

void R_GenTurbTile(u8 *pbasetex, void *pdest)
{
	s32 *turb = sintable + ((s32)(cl.time * SPEED) & (CYCLE - 1));
	u8 *pd = (u8 *) pdest;
	for (s32 i = 0; i < TILE_SIZE; i++) {
		for (s32 j = 0; j < TILE_SIZE; j++) {
			s32 s = (((j << 16) + turb[i & (CYCLE-1)]) >> 16) & 63;
			s32 t = (((i << 16) + turb[j & (CYCLE-1)]) >> 16) & 63;
			*pd++ = *(pbasetex + (t << 6) + s);
		}
	}
}
