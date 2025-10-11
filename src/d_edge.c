// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

static s32 miplevel;
static vec3_t transformed_modelorg;
static vec3_t world_transformed_modelorg;

s32 D_MipLevelForScale(f32 scale)
{
	s32 lmiplevel;
	scale *= r_mipscale.value;
	if (scale >= d_scalemip[0]) lmiplevel = 0;
	else if (scale >= d_scalemip[1]) lmiplevel = 1;
	else if (scale >= d_scalemip[2]) lmiplevel = 2;
	else lmiplevel = 3;
	if (lmiplevel < d_minmip) lmiplevel = d_minmip;
	return lmiplevel;
}

void D_DrawSolidSurface(surf_t *surf, s32 color)
{ // FIXME: clean this up
	s32 pix = (color << 24) | (color << 16) | (color << 8) | color;
	for (espan_t *span = surf->spans; span; span = span->pnext) {
		u8 *pdest = (u8 *) d_viewbuffer + screenwidth * span->v;
		s32 u = span->u;
		s32 u2 = span->u + span->count - 1;
		((u8 *) pdest)[u] = pix;
		if (u2 - u < 8)
			for (u++; u <= u2; u++)
				((u8 *) pdest)[u] = pix;
		else {
			for (u++; u & 3; u++)
				((u8 *) pdest)[u] = pix;
			u2 -= 4;
			for (; u <= u2; u += 4)
				*(s32 *)((u8 *) pdest + u) = pix;
			u2 += 4;
			for (; u <= u2; u++)
				((u8 *) pdest)[u] = pix;
		}
	}
}

void D_CalcGradients(msurface_t *pface)
{
	f32 mipscale = 1.0 / (f32)(1 << miplevel);
	vec3_t p_saxis, p_taxis;
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
	vec3_t p_temp1;
	VectorScale(transformed_modelorg, mipscale, p_temp1);
	t = 0x10000 * mipscale;
	sadjust = ((s32) (DotProduct(p_temp1, p_saxis) * 0x10000 + 0.5)) -
				((pface->texturemins[0] << 16) >> miplevel)
				+ pface->texinfo->vecs[0][3] * t;
	tadjust = ((s32) (DotProduct(p_temp1, p_taxis) * 0x10000 + 0.5)) -
				((pface->texturemins[1] << 16) >> miplevel)
				+ pface->texinfo->vecs[1][3] * t;
	// -1 (-epsilon) so we never wander off the edge of the texture
	bbextents = ((pface->extents[0] << 16) >> miplevel) - 1;
	bbextentt = ((pface->extents[1] << 16) >> miplevel) - 1;
}

void D_DrawSurfacesFlat()
{
	for (surf_t *s = &surfaces[1]; s < surface_p; s++) {
		if (!s->spans) continue;
		d_zistepu = s->d_zistepu;
		d_zistepv = s->d_zistepv;
		d_ziorigin = s->d_ziorigin;
		D_DrawSolidSurface(s, (uintptr_t) s->data & 0xFF);
		D_DrawZSpans(s->spans);
	}
}

static void D_DrawSky(surf_t *s)
{
	if (!r_skymade) R_MakeSky();
	D_DrawSkyScans8(s->spans);
	D_DrawZSpans(s->spans);
}

static void D_DrawSkybox(surf_t *s, msurface_t *pface)
{ // Manoel Kasimier
	extern u8 r_skypixels[6][SKYBOX_MAX_SIZE*SKYBOX_MAX_SIZE];
	cacheblock = (u8 *)(r_skypixels[pface->texinfo->texture->offsets[0]]);
	cachewidth = pface->texinfo->texture->width;
	d_zistepu = s->d_zistepu;
	d_zistepv = s->d_zistepv;
	d_ziorigin = s->d_ziorigin;
	D_CalcGradients (pface);
	if(fog_density > 0 && !fog_lut_built) build_color_mix_lut(0);
	D_DrawSpans(s->spans, SPAN_SKYBOX, 0);
	d_zistepu = 0; // set up gradient for background surf that places it
	d_zistepv = 0; // effectively at infinity distance from the viewpoint
	d_ziorigin = -0.9;
	D_DrawZSpans (s->spans);
}

static void D_DrawBackground(surf_t *s)
{
	d_zistepu = 0; // set up gradient for background surf that places it
	d_zistepv = 0; // effectively at infinity distance from the viewpoint
	d_ziorigin = -0.9;
	D_DrawSolidSurface(s, (s32)r_clearcolor.value & 0xFF);
	D_DrawZSpans(s->spans);
}

static void D_DrawUnlitWater(surf_t *s, msurface_t *pface)
{ // Manoel Kasimier
	cacheblock = (u8 *) ((u8 *) pface->texinfo->texture + pface->texinfo->texture->offsets[0]);
	cachewidth = 64;
	D_CalcGradients(pface);
	f32 opacity = 1;
	if (s->entity && s->entity->alpha && r_entalpha.value == 1)
		opacity -= (f32)s->entity->alpha/255;
	else if (s->flags & SURF_WINQUAKE_DRAWTRANSLUCENT)
		opacity -= R_LiquidAlphaForFlags(s->flags);
	Turbulent8(s->spans, opacity);
	if (!r_alphapass) D_DrawZSpans(s->spans);
}

static void D_DrawTransSurf(surf_t *s, msurface_t *pface)
{
	surfcache_t *pcurrentcache = D_CacheSurface(pface, miplevel);
	cacheblock = (u8 *) pcurrentcache->data;
	cachewidth = pcurrentcache->width;
	D_CalcGradients(pface);
	f32 opacity = 1 - (f32)s->entity->alpha / 255;
	D_DrawSpans(s->spans, SPAN_TRANS, opacity);
}

static void D_DrawLitWater(surf_t *s, msurface_t *pface)
{ // FIXME this is horrible.
	miplevel = D_MipLevelForScale(s->nearzi * scale_for_mip
			* pface->texinfo->mipadjust);
	lmonly = 1; // this is how we know it's lit water that we're drawing
	surfcache_t *pcurrentcache = D_CacheSurface(pface, miplevel);
	cacheblock = (u8 *) pcurrentcache->data;
	cachewidth = pcurrentcache->width;
	D_CalcGradients(pface);
	D_DrawSpans(s->spans, SPAN_NORMAL, 0); // draw the lightmap to a separate buffer
	miplevel = 0;
	cacheblock = (u8 *) pface->texinfo->texture + pface->texinfo->texture->offsets[0];
	cachewidth = 64;
	D_CalcGradients(pface);
	f32 opacity = 1;
	if (s->entity && s->entity->alpha && r_entalpha.value == 1)
		opacity -= (f32)s->entity->alpha/255;
	else if (s->flags & SURF_WINQUAKE_DRAWTRANSLUCENT)
		opacity -= R_LiquidAlphaForFlags(s->flags);
	Turbulent8(s->spans, opacity);
	if (!r_alphapass) D_DrawZSpans(s->spans);
	lmonly = 0;
}

static void D_DrawCutoutSurf(surf_t *s, msurface_t *pface)
{
	surfcache_t *pcurrentcache = D_CacheSurface(pface, miplevel);
	cacheblock = (u8 *) pcurrentcache->data;
	cachewidth = pcurrentcache->width;
	D_CalcGradients(pface);
	D_DrawSpans(s->spans, SPAN_CUTOUT, 0);
	D_DrawZSpansTrans(s->spans);
}

static void D_DrawNormalSurf(surf_t *s, msurface_t *pface)
{
	surfcache_t *pcurrentcache = D_CacheSurface(pface, miplevel);
	cacheblock = (u8 *) pcurrentcache->data;
	cachewidth = pcurrentcache->width;
	D_CalcGradients(pface);
	D_DrawSpans(s->spans, SPAN_NORMAL, 0);
	D_DrawZSpans(s->spans);
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

void D_DrawSurfaces()
{
	currententity = &cl_entities[0];
	TransformVector(modelorg, transformed_modelorg);
	VectorCopy(transformed_modelorg, world_transformed_modelorg);
	s32 foundcutouts = 0;
	for (surf_t *s = &surfaces[1]; s < surface_p; s++) {
		if (!s->spans) continue;
		msurface_t *pface = s->data;
		if (s->flags & SURF_DRAWBACKGROUND) {
			if (skybox_name[0]) continue;
			D_DrawBackground(s);
			continue;
		}
		if (!pface || !pface->extents[0] || !pface->extents[1]) {
			if(pface)Con_DPrintf("Broken surface extents %hd %hd\n",
					pface->extents[0], pface->extents[1]);
			else Con_DPrintf("Broken surface\n");
			s->spans = NULL;
			continue;
		}
		d_zistepu = s->d_zistepu;
		d_zistepv = s->d_zistepv;
		d_ziorigin = s->d_ziorigin;
		if(s->flags&(SURF_DRAWTURB|SURF_DRAWBACKGROUND|SURF_DRAWSKYBOX))
			miplevel = 0;
		else miplevel = D_MipLevelForScale(s->nearzi * scale_for_mip
				* pface->texinfo->mipadjust);
		if (s->insubmodel) D_SwitchSubModelOn(s);
		if (s->flags & SURF_DRAWSKY) {
			D_DrawSky(s);
		} else if (s->flags & SURF_DRAWSKYBOX) {
			D_DrawSkybox(s, pface);
		} else if (s->flags & SURF_DRAWTURB &&
			(!s->entity->model->haslitwater || !r_litwater.value)) {
			D_DrawUnlitWater(s, pface);
		} else if (s->flags & SURF_DRAWCUTOUT) {
			foundcutouts = 1;
		} else D_DrawNormalSurf(s, pface);
		if (s->insubmodel) D_SwitchSubModelOff();
		if (!(s->flags & SURF_DRAWCUTOUT))
			s->spans = NULL;
	}
	if (!foundcutouts) return;
	for (surf_t *s = &surfaces[1]; s < surface_p; s++) {
		if (!s->spans) continue;
		msurface_t *pface = s->data;
		if (!pface || !pface->extents[0] || !pface->extents[1]) {
			if(pface)Con_DPrintf("Broken surface extents %hd %hd\n",
					pface->extents[0], pface->extents[1]);
			else Con_DPrintf("Broken surface\n");
			s->spans = NULL;
			continue;
		}
		d_zistepu = s->d_zistepu;
		d_zistepv = s->d_zistepv;
		d_ziorigin = s->d_ziorigin;
		miplevel = D_MipLevelForScale(s->nearzi * scale_for_mip
				* pface->texinfo->mipadjust);
		if (s->insubmodel) D_SwitchSubModelOn(s);
		D_DrawCutoutSurf(s, pface);
		if (s->insubmodel) D_SwitchSubModelOff();
	}
}

void D_DrawSurfacesAlpha()
{
	currententity = &cl_entities[0];
	TransformVector(modelorg, transformed_modelorg);
	VectorCopy(transformed_modelorg, world_transformed_modelorg);
	for (surf_t *s = &surfaces[1]; s < surface_p; s++) {
		if (!s->spans) continue;
		msurface_t *pface = s->data;
		u64 is_ent = (ul64)s->entity & 0xffff000ull;
		if(!(s->flags&SURF_DRAWTURB) && !is_ent) continue;
		if (pface == 0) continue;
		d_zistepu = s->d_zistepu;
		d_zistepv = s->d_zistepv;
		d_ziorigin = s->d_ziorigin;
		miplevel = 0;
		if (s->insubmodel) D_SwitchSubModelOn(s);
		if (s->flags & SURF_DRAWTURB && (!s->entity->model->haslitwater || !r_litwater.value)) {
			D_DrawUnlitWater(s, pface);
		} else if (is_ent && s->entity->alpha && r_entalpha.value == 1) {
			D_DrawTransSurf(s, pface);
		} else if (s->flags & SURF_DRAWTURB && s->entity->model->haslitwater && r_litwater.value) {
			D_DrawLitWater(s, pface);
		}
		if (s->insubmodel) D_SwitchSubModelOff();
	}
}
