// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

static vec3_t transformed_modelorg;
static vec3_t world_transformed_modelorg;

static s32 D_MipLevelForScale(f32 scale)
{
	s32 lmiplevel;
	scale *= r_mipscale.value;
	if(scale >= d_scalemip[0]) lmiplevel = 0;
	else if(scale >= d_scalemip[1]) lmiplevel = 1;
	else if(scale >= d_scalemip[2]) lmiplevel = 2;
	else lmiplevel = 3;
	if(lmiplevel < d_minmip) lmiplevel = d_minmip;
	return lmiplevel;
}

static void D_DrawSolidSurface(surf_t *surf, s32 color)
{
	r_drawnpolycount++;
	for(espan_t *span = surf->spans; span; span = span->pnext)
		memset(d_viewbuffer + screenwidth * span->v + span->u, color,
				span->count);
}

// 16.16 texture offset as id computed it, in s64 so that huge maps (|s| > 32k)
// don't overflow before it cancels out. x87 rounds the sum to 24 bits but not
// the integer going into it
static inline s64 D_Adjust(const vec3_t org, const vec3_t axis, f32 ofs,
		f32 t, s64 mins)
{
	f32 d = org[0] * axis[0] + (org[1] * axis[1] + org[2] * axis[2]);
	d = d * 0x10000 + 0.5f;
	return (s64)(f32)((f64)((s64)d - mins) + ofs * t);
}

static void D_CalcGradients(msurface_t *pface, s32 tiledextents)
{
	f32 mipscale = 1.0 / (f32)(1 << miplevel);
	vec3_t p_saxis, p_taxis, p_temp;
	TransformVector(pface->texinfo->vecs[0], p_saxis);
	TransformVector(pface->texinfo->vecs[1], p_taxis);
	f32 t = xscaleinv * mipscale;
	d_sdivzstepu = p_saxis[0] * t;
	d_tdivzstepu = p_taxis[0] * t;
	t = yscaleinv * mipscale;
	d_sdivzstepv = -p_saxis[1] * t;
	d_tdivzstepv = -p_taxis[1] * t;
	d_sdivzorigin = p_saxis[2] * mipscale - xcenter * d_sdivzstepu -
						ycenter * d_sdivzstepv;
	d_tdivzorigin = p_taxis[2] * mipscale - xcenter * d_tdivzstepu -
						ycenter * d_tdivzstepv;
	VectorScale(transformed_modelorg, mipscale, p_temp);
	t = 0x10000 * mipscale;
	mtexinfo_t *tex = pface->texinfo;
	if(!tiledextents){
		sadjust = D_Adjust(p_temp, p_saxis, tex->vecs[0][3], t,
			((s64)pface->texturemins[0] * 0x10000) >> miplevel);
		tadjust = D_Adjust(p_temp, p_taxis, tex->vecs[1][3], t,
			((s64)pface->texturemins[1] * 0x10000) >> miplevel);
		// -1(-epsilon) so we never wander off the edge of the texture
		bbextents = ((pface->extents[0] << 16) >> miplevel) - 1;
		bbextentt = ((pface->extents[1] << 16) >> miplevel) - 1;
	}
	else{ // Tiled textures(water, sky) always have miplevel of 0
		// WinQuake's -8192 mins, then whole tiles to stay near them
		s64 mins = (s64)-8192 * 0x10000;
		sadjust = D_Adjust(p_temp, p_saxis, tex->vecs[0][3], t, mins)
			- (s64)(pface->texturemins[0] & -TURB_CYCLE) * 0x10000;
		tadjust = D_Adjust(p_temp, p_taxis, tex->vecs[1][3], t, mins)
			- (s64)(pface->texturemins[1] & -TURB_CYCLE) * 0x10000;
		bbextents = 0x3fffffff;
		bbextentt = 0x3fffffff;
	}
}

void D_DrawSurfacesFlat()
{
	r_drawnpolycount++;
	for(surf_t *s = &surfaces[1]; s < surface_p; s++){
		if(!s->spans) continue;
		d_zistepu = s->d_zistepu;
		d_zistepv = s->d_zistepv;
		d_ziorigin = s->d_ziorigin;
		D_DrawSolidSurface(s, (uintptr_t) s->data & 0xFF);
		D_DrawZSpans(s->spans);
	}
}

static void D_SetZGradients(const surf_t *s)
{
	d_zistepu = s->d_zistepu;
	d_zistepv = s->d_zistepv;
	d_ziorigin = s->d_ziorigin;
}

static void D_DrawZSpansFar(espan_t *spans)
{
	d_zistepu = 0; // set up gradient for background surf that places it
	d_zistepv = 0; // effectively at infinity distance from the viewpoint
	d_ziorigin = -0.9;
	D_DrawZSpans(spans);
}

static void D_DrawSky(surf_t *s, msurface_t *pface)
{
	r_drawnpolycount++;
	s32 skyi = R_SkyIndexForTexture(pface->texinfo->texture);
	if(!r_skies[skyi].made && !r_skysmooth.value)
		R_MakeSky(pface->texinfo->texture);
	if(fog_density>0){
		if(map_skyfog>=1) D_DrawSkyScansOnlyFog(s->spans);
		else D_DrawSkyScansFog(s->spans, pface);
	}
	else D_DrawSkyScans(s->spans, pface);
	D_DrawZSpans(s->spans);
}

static void D_DrawSkybox(surf_t *s)
{ // Manoel Kasimier
	r_drawnpolycount++;
	msurface_t *pface = s->data;
	texture_t *tex = pface->texinfo->texture;
	cacheblock = (u8*)(r_skypixels[tex->offsets[0]]);
	cachewidth = tex->width;
	cacheheight = tex->height;
	miplevel = 0;
	D_SetZGradients(s);
	D_CalcGradients(pface, 0);
	if(fog_density > 0 && !fog_lut_built) R_BuildColorMixLUT(0);
	if(r_dithertex.value)
		D_DrawSpansDithered(s->spans, SPAN_SKYBOX, 0);
	else
		D_DrawSpans(s->spans, SPAN_SKYBOX, 0);
	D_DrawZSpansFar(s->spans);
	R_FogMarkSky(s->spans);
}

// Faces too big for a surface cache block at the requested mip level
// (vanilla limits them to 256, GL ports to 2000) use a coarser mip instead
static void D_FitMipLevel(msurface_t *pface)
{
	while(miplevel < MIPLEVELS-1 && ((pface->extents[0] >> miplevel) > 256
			|| (pface->extents[1] >> miplevel) > 256))
		miplevel++;
}

static s32 D_SurfMipLevel(const surf_t *s, msurface_t *pface)
{ // for the nearest point of the surface
	f32 scale = s->nearzi * scale_for_mip * pface->texinfo->mipadjust;
	if(s->flags & SURF_DRAWCUTOUT)
		scale *= r_cutoutmipscale.value;
	return D_MipLevelForScale(scale);
}

// Opaque, cutout and translucent faces, from the surface cache. Inlined so
// that every caller gets the span drawing for its own type.
Q_SIL void D_DrawCachedSurf(surf_t *s, msurface_t *pface, s32 type,
		f32 opacity)
{
	miplevel = D_SurfMipLevel(s, pface);
	D_FitMipLevel(pface);
	surfcache_t *pcurrentcache = D_CacheSurface(pface, miplevel);
	if(pcurrentcache == NULL){ s->spans = 0; return; }
	r_drawnpolycount++;
	cacheblock = (u8*) pcurrentcache->data;
	cachewidth = pcurrentcache->width;
	cacheheight = pcurrentcache->height;
	D_CalcGradients(pface, 0);
	if(r_dithertex.value && !miplevel)
		D_DrawSpansDithered(s->spans, type, opacity);
	else
		D_DrawSpans(s->spans, type, opacity);
	if(type == SPAN_NORMAL)
		D_DrawZSpans(s->spans);
	else if(type == SPAN_CUTOUT)
		D_DrawZSpansTrans(s->spans);
}

static void D_DrawWater(surf_t *s, msurface_t *pface, f32 opacity)
{ // Manoel Kasimier
	r_drawnpolycount++;
	cacheblock = (u8*)((u8*)pface->texinfo->texture
			+ pface->texinfo->texture->offsets[0]);
	cachewidth = 64;
	cacheheight = 64;
	D_CalcGradients(pface, 1);
	Turbulent(s->spans, opacity);
	if(!r_alphapass) D_DrawZSpans(s->spans);
}

static void D_FitLitWaterMip(msurface_t *pface)
{ // RGB light rows take 4 bytes per texel, D_SCAlloc caps a block at 64k
	s32 bpp = R_LitWaterRGB() ? 4 : 1;
	D_FitMipLevel(pface);
	while(miplevel < MIPLEVELS-1 && (pface->extents[0] >> miplevel)
			* (pface->extents[1] >> miplevel) * bpp > 0x10000)
		miplevel++;
}

static void D_DrawTurb(surf_t *s, msurface_t *pface, f32 opacity)
{
	miplevel = 0;
	if(!s->entity->model->haslitwater || !r_litwater.value
			|| (pface->flags & SURF_DRAWTILED)){ // unlit
		D_DrawWater(s, pface, opacity);
		return;
	}
	// the cache holds only light rows, Turbulent applies them per pixel
	miplevel = D_SurfMipLevel(s, pface);
	D_FitLitWaterMip(pface);
	r_lightonly = 1; // tells D_CacheSurface it's caching lit water
	surfcache_t *pcurrentcache = D_CacheSurface(pface, miplevel);
	r_lightonly = 0;
	s32 lightmip = miplevel;
	miplevel = 0; // the water texture itself is tiled at full res
	D_TurbSetLight(pcurrentcache, pface, lightmip); // unlit if NULL
	D_DrawWater(s, pface, opacity);
	D_TurbSetLight(NULL, NULL, 0);
}

static void D_SwitchSubModelOn(surf_t *s)
{
	currententity = s->entity;
	vec3_t local_modelorg;
	VectorSubtract(r_origin, currententity->origin, local_modelorg);
	TransformVector(local_modelorg, transformed_modelorg);
	R_RotateBmodel();
}

static void D_SwitchSubModelOff()
{
	currententity = &cl_entities[0];
	VectorCopy(world_transformed_modelorg, transformed_modelorg);
	VectorCopy(base_vpn, vpn);
	VectorCopy(base_vup, vup);
	VectorCopy(base_vright, vright);
	VectorCopy(base_modelorg, modelorg);
	R_TransformFrustum();
}

// Brush entity faces are drawn in the space of their entity. Those of one
// entity come one after another, so it only switches when the entity changes.
// NULL goes back to the world.
static void D_SwitchEntity(surf_t *s)
{
	entity_t *e = s && s->insubmodel ? s->entity : &cl_entities[0];
	if(e == currententity)
		return;
	if(currententity != &cl_entities[0])
		D_SwitchSubModelOff();
	if(e != &cl_entities[0])
		D_SwitchSubModelOn(s);
}

static void D_BeginSurfaces()
{
	currententity = &cl_entities[0];
	TransformVector(modelorg, transformed_modelorg);
	VectorCopy(transformed_modelorg, world_transformed_modelorg);
}

static bool D_BrokenSurf(surf_t *s, msurface_t *pface)
{
	if(pface && pface->extents[0] && pface->extents[1])
		return 0;
	if(pface)Con_DPrintf("Broken surface extents %hd %hd\n",
			pface->extents[0], pface->extents[1]);
	else Con_DPrintf("Broken surface\n");
	s->spans = NULL;
	return 1;
}

// Cutouts are drawn over all the other surfaces, they're z-tested. Faces too
// big to draw are drawn black with the rest.
static bool D_CutoutPass(const surf_t *s)
{ return (s->flags & (SURF_DRAWCUTOUT | SURF_NOTEXTURE)) == SURF_DRAWCUTOUT; }

void D_DrawSurfaces()
{
	D_BeginSurfaces();
	if(skybox_name[0] && r_skyframe == r_framecount){
		for(surf_t *s = skybox_surf_p; s < surface_p; s++){
			if(s->spans && (s->flags & SURF_DRAWSKYBOX)){
				D_DrawSkybox(s);
				s->spans = NULL;
			}
		}
	}
	bool foundcutouts = 0;
	for(surf_t *s = &surfaces[1]; s < surface_p; s++){
		if(!s->spans) continue;
		msurface_t *pface = s->data;
		if(s->flags & SURF_DRAWBACKGROUND){
			if(skybox_name[0]) continue;
			D_DrawSolidSurface(s, (s32)r_clearcolor.value & 0xFF);
			D_DrawZSpansFar(s->spans);
			s->spans = NULL;
			continue;
		}
		if(D_CutoutPass(s)){
			foundcutouts = 1;
			continue;
		}
		if(D_BrokenSurf(s, pface))
			continue;
		D_SetZGradients(s);
		D_SwitchEntity(s);
		miplevel = 0;
		if(s->flags & SURF_NOTEXTURE){
			Con_DPrintf("D_DrawSurfaces: drawing "
					"SURF_NOTEXTURE\n");
			D_DrawSolidSurface(s, 0);
			D_DrawZSpans(s->spans);
		}else if(s->flags & SURF_DRAWSKY){
			if(!skybox_name[0]) D_DrawSky(s, pface);
			else D_DrawZSpans(s->spans);
			R_FogMarkSky(s->spans);
		}else if(s->flags & SURF_DRAWTURB){
			// R_SurfaceInPass kept the translucent ones for later,
			// everything that made it here is opaque
			D_DrawTurb(s, pface, 0);
		}else
			D_DrawCachedSurf(s, pface, SPAN_NORMAL, 0);
		s->spans = NULL;
	}
	for(surf_t *s = &surfaces[1]; foundcutouts && s < surface_p; s++){
		if(!s->spans || !D_CutoutPass(s)) continue;
		msurface_t *pface = s->data;
		if(D_BrokenSurf(s, pface))
			continue;
		D_SetZGradients(s);
		D_SwitchEntity(s);
		D_DrawCachedSurf(s, pface, SPAN_CUTOUT, 0);
	}
	D_SwitchEntity(NULL);
}

// Translucent surfaces get all of their pixels, so they have to be blended
// back to front. A lower key is nearer, faces with keys of their own sort
// exactly. Brush entity faces share the key of their leaf, between those the
// one with the farther nearest point goes first.
static s32 D_AlphaSurfCmp(const void *pa, const void *pb)
{
	const surf_t *a = *(surf_t *const *)pa, *b = *(surf_t *const *)pb;
	if(a->key != b->key)
		return (a->key < b->key) == !r_draworder.value ? 1 : -1;
	if(a->nearzi != b->nearzi)
		return a->nearzi > b->nearzi ? 1 : -1;
	return a < b ? -1 : a > b; // keep the emission order
}

void D_DrawSurfacesAlpha()
{
	D_BeginSurfaces();
	static surf_t **order;
	static s32 order_cap;
	if(order_cap < r_lsurfs_size){
		free(order);
		order = malloc(r_lsurfs_size * sizeof(*order));
		if(!order)
			Sys_Error("D_DrawSurfacesAlpha: out of memory");
		order_cap = r_lsurfs_size;
	}
	s32 count = 0;
	for(surf_t *s = &surfaces[1]; s < surface_p; s++)
		if(s->spans)
			order[count++] = s;
	qsort(order, count, sizeof(*order), D_AlphaSurfCmp);
	for(s32 i = 0; i < count; i++){
		surf_t *s = order[i];
		msurface_t *pface = s->data;
		D_SetZGradients(s);
		D_SwitchEntity(s);
		if(s->flags & SURF_DRAWTURB){
			f32 opacity = 1;
			if(s->entity->alpha && r_entalpha.value)
				opacity -= (f32)s->entity->alpha / 255;
			else if(s->flags & SURF_WINQUAKE_DRAWTRANSLUCENT)
				opacity -= R_LiquidAlphaForFlags(s->flags);
			D_DrawTurb(s, pface, opacity);
		}else if(s->entity->alpha && r_entalpha.value == 1)
			D_DrawCachedSurf(s, pface, SPAN_TRANS,
					1 - (f32)s->entity->alpha / 255);
	}
	D_SwitchEntity(NULL);
}
