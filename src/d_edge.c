// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

// How the spans of a surface get drawn. D_SetupSurf works it out and
// D_DrawSurf needs nothing else, so with r_threads the workers draw the spans
// of each band while the main thread goes on (r_thread.c). The setup of a
// band's surfaces can be a job too (D_RecordBandJob), so what it works with is
// the thread's own.
// DS_SKY and DS_SKYZ mark the sky, R_FogMarkSky fogs the colors there. DS_SKYZ
// only writes z, over a skybox's colors, or over nothing under full sky fog,
// which hides the box so it isn't drawn. The skybox pass writes no z and isn't
// marked: every pixel of it gets a surface or the background drawn over it in
// the main pass. What sky faces don't cover is the void, the background
// (DS_FAR) draws it in r_clearcolor, and R_DrawFog leaves its far z alone.
// DS_LATER has to have its block built first.
enum { DS_NONE, DS_SOLID, DS_FAR, DS_SKY, DS_SKYZ, DS_SKYBOX, DS_TURB,
	DS_CACHED, DS_LATER };
typedef struct dsurf_s {
	s32 kind;
	s32 type; // SPAN_* of DS_CACHED
	s32 color; // DS_SOLID and DS_FAR
	bool dithered;
	bool zbuffer; // DS_TURB, opaque water writes its z
	bool early; // can be drawn before the setup walk comes to it, D_Surf
	u32 walk; // the setup walk that last came to it
	f32 opacity;
	dgrad_t g;
	msurface_t *pface;
	surfcache_t *cache; // the block this was set up with, if any
	surfcache_t **cachespot; // where the surface keeps it
	texture_t *texture; // what was in it
	vec3_t sky[3]; // the view of DS_SKY
	turbstate_t turb;
} dsurf_t;
typedef struct { const dsurf_t *d; espan_t *spans; } dsurfcmd_t;

static Q_TLS vec3_t transformed_modelorg;
static Q_TLS vec3_t world_transformed_modelorg;
static Q_TLS vec3_t d_vpn, d_vright, d_vup; // the view in d_entity's space
static Q_TLS entity_t *d_entity; // the one being set up
static Q_TLS bool d_nobuild; // set up only what needs no cache block built
static Q_TLS s32 d_npolys; // set up, for r_drawnpolycount
static SDL_AtomicInt d_jobpolys; // of the jobs
static SDL_Mutex *d_lazy; // for what's built the first time it's needed
static u32 d_walk; // counts the setup walks
static Q_TLS bool d_injob; // D_RecordBandJob
// The cache blocks the jobs found, which the main thread marks used once
// they're done (D_StampBlocks). Two jobs can look up the same block: a brush
// face cut by the world's leaves is a surface for each piece.
typedef struct { surfcache_t **c; s32 n, cap; } dfound_t;
static dfound_t d_found[R_MAXWORKERS + 1];
#define D_BUSY (-1) // r_surfscan: the surface is being set up

// Tables and the console aren't the jobs' to change at the same time
static void D_Lazy(bool lock)
{
	if(!d_lazy)
		return;
	if(lock) SDL_LockMutex(d_lazy);
	else SDL_UnlockMutex(d_lazy);
}

// TransformVector in d_entity's space
static void D_Transform(const vec3_t in, vec3_t out)
{
	out[0] = DotProduct(in, d_vright);
	out[1] = DotProduct(in, d_vup);
	out[2] = DotProduct(in, d_vpn);
}

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

static void D_DrawSolidSurface(espan_t *span, s32 color)
{
	for(; span; span = span->pnext)
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
	D_Transform(pface->texinfo->vecs[0], p_saxis);
	D_Transform(pface->texinfo->vecs[1], p_taxis);
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
		// WinQuake's -8192 mins, then whole tiles to stay near them.
		// The coords wrap around, no extents (D_WalkCoord).
		s64 mins = (s64)-8192 * 0x10000;
		sadjust = D_Adjust(p_temp, p_saxis, tex->vecs[0][3], t, mins)
			- (s64)(pface->texturemins[0] & -TURB_CYCLE) * 0x10000;
		tadjust = D_Adjust(p_temp, p_taxis, tex->vecs[1][3], t, mins)
			- (s64)(pface->texturemins[1] & -TURB_CYCLE) * 0x10000;
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
	d_ziorigin = ZI_BACKGROUND;
	D_DrawZSpans(spans);
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

static void D_FitLitWaterMip(msurface_t *pface)
{ // RGB light rows take 4 bytes per texel, D_SCAlloc caps a block at 64k
	s32 bpp = R_LitWaterRGB() ? 4 : 1;
	D_FitMipLevel(pface);
	while(miplevel < MIPLEVELS-1 && (pface->extents[0] >> miplevel)
			* (pface->extents[1] >> miplevel) * bpp > 0x10000)
		miplevel++;
}

static void D_SwitchSubModelOn(surf_t *s)
{ // as R_RotateBmodel turns the view
	d_entity = s->entity;
	vec3_t local_modelorg;
	VectorSubtract(r_origin, d_entity->origin, local_modelorg);
	D_Transform(local_modelorg, transformed_modelorg);
	f32 m[3][3];
	R_EntityRotation(d_entity, m);
	R_RotateBy(m, d_vpn);
	R_RotateBy(m, d_vright);
	R_RotateBy(m, d_vup);
}

static void D_SwitchSubModelOff()
{
	d_entity = &cl_entities[0];
	VectorCopy(world_transformed_modelorg, transformed_modelorg);
	VectorCopy(base_vpn, d_vpn);
	VectorCopy(base_vup, d_vup);
	VectorCopy(base_vright, d_vright);
}

static Q_NOINLINE void D_SwitchTo(entity_t *e, surf_t *s)
{
	if(d_entity != &cl_entities[0])
		D_SwitchSubModelOff();
	if(e != &cl_entities[0])
		D_SwitchSubModelOn(s);
}

// Brush entity faces are drawn in the space of their entity. Those of one
// entity come one after another, so it only switches when the entity changes.
// NULL goes back to the world.
Q_SIL void D_SwitchEntity(surf_t *s)
{
	entity_t *e = s && s->insubmodel ? s->entity : &cl_entities[0];
	if(e != d_entity)
		D_SwitchTo(e, s);
}

static void D_BeginSurfaces()
{
	d_entity = &cl_entities[0];
	VectorCopy(base_vpn, d_vpn);
	VectorCopy(base_vup, d_vup);
	VectorCopy(base_vright, d_vright);
	D_Transform(base_modelorg, transformed_modelorg);
	VectorCopy(transformed_modelorg, world_transformed_modelorg);
}

static bool D_BrokenSurf(surf_t *s, msurface_t *pface)
{
	if(pface && pface->extents[0] && pface->extents[1])
		return 0;
	D_Lazy(1);
	if(pface)Con_DPrintf("Broken surface extents %hd %hd\n",
			pface->extents[0], pface->extents[1]);
	else Con_DPrintf("Broken surface\n");
	D_Lazy(0);
	s->spans = NULL;
	return 1;
}

// Cutouts are drawn over all the other surfaces, they're z-tested. Faces too
// big to draw are drawn black with the rest.
static bool D_CutoutPass(const surf_t *s)
{ return (s->flags & (SURF_DRAWCUTOUT | SURF_NOTEXTURE)) == SURF_DRAWCUTOUT; }

// r_lerplightstyles keeps a block while its light is close enough. A block
// that's only close gets built again when it's lost, so it can't be drawn from
// before the setup walk comes to its surface.
static bool D_CacheExact(const surfcache_t *c)
{
	const s32 *a = c->lightadj, *b = r_drawsurf.lightadj;
	return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
}

// Full sky fog hides all of the sky, its colors aren't drawn
bool D_SkyAllFog()
{ return fog_density && map_skyfog >= 1; }

static void D_SetupSky(dsurf_t *d, msurface_t *pface)
{
	d_npolys++;
	sky_t *sky = &r_skies[R_SkyIndexForTexture(pface->texinfo->texture)];
	d->kind = D_SkyAllFog() ? DS_SKYZ : DS_SKY;
	D_Lazy(1);
	if(!sky->made && !r_skysmooth.value)
		R_MakeSky(pface->texinfo->texture);
	if(d->kind == DS_SKY && r_skysmooth.value && sky->blend)
		R_BuildColorMixLUT(0);
	D_Lazy(0);
	d->pface = pface;
	VectorCopy(d_vpn, d->sky[0]);
	VectorCopy(d_vright, d->sky[1]);
	VectorCopy(d_vup, d->sky[2]);
}

static void D_SetupSkybox(dsurf_t *d, surf_t *s)
{ // Manoel Kasimier
	d_npolys++;
	msurface_t *pface = s->data;
	texture_t *tex = pface->texinfo->texture;
	cacheblock = (u8*)(r_skypixels[tex->offsets[0]]);
	cachewidth = tex->width;
	cacheheight = tex->height;
	miplevel = 0;
	D_CalcGradients(pface, 0);
	d->kind = DS_SKYBOX;
	d->dithered = r_dithertex.value != 0;
}

// The surface's block at miplevel, built if it has to be. While a band is
// recorded it's only looked up, and without one the surface waits for the
// setup walk.
static surfcache_t *D_SurfBlock(dsurf_t *d, msurface_t *pface)
{
	if(!d_nobuild)
		return D_CacheSurface(pface, miplevel, d_entity->frame);
	surfcache_t *c = D_CacheLookup(pface, miplevel, d_entity->frame);
	if(!c){
		d->kind = DS_LATER;
		d->early = 0;
	}else if(d_injob){
		dfound_t *f = &d_found[r_slot];
		if(f->n == f->cap){
			f->cap = f->cap ? f->cap * 2 : 256;
			f->c = realloc(f->c, f->cap * sizeof(*f->c));
			if(!f->c)
				Sys_Error("D_SurfBlock: out of memory");
		}
		f->c[f->n++] = c;
	}else if(r_threaded)
		c->frame = r_framecount & 0x7fffffff; // see D_BlockUsed
	return c;
}

// The block the surface is drawn from, see D_SurfStillValid
static void D_KeepBlock(dsurf_t *d, surfcache_t *c, msurface_t *pface, s32 mip)
{
	d->cache = c;
	d->cachespot = &pface->cachespots[mip];
	d->texture = c->texture;
	d->early = D_CacheExact(c);
}

// Opaque, cutout and translucent faces, from the surface cache
static void D_SetupCached(dsurf_t *d, surf_t *s, msurface_t *pface, s32 type,
		f32 opacity)
{
	miplevel = D_SurfMipLevel(s, pface);
	D_FitMipLevel(pface);
	surfcache_t *pcurrentcache = D_SurfBlock(d, pface);
	if(pcurrentcache == NULL)
		return;
	d_npolys++;
	cacheblock = (u8*) pcurrentcache->data;
	cachewidth = pcurrentcache->width;
	cacheheight = pcurrentcache->height;
	D_CalcGradients(pface, 0);
	if(type == SPAN_TRANS && r_alphastyle.value == 0){
		D_Lazy(1);
		R_BuildColorMixLUT(0);
		D_Lazy(0);
	}
	d->kind = DS_CACHED;
	d->type = type;
	d->opacity = opacity;
	d->dithered = r_dithertex.value && !miplevel;
	D_KeepBlock(d, pcurrentcache, pface, miplevel);
}

static void D_SetupWater(dsurf_t *d, msurface_t *pface, f32 opacity)
{ // Manoel Kasimier
	d_npolys++;
	cacheblock = (u8*)((u8*)pface->texinfo->texture
			+ pface->texinfo->texture->offsets[0]);
	cachewidth = 64;
	cacheheight = 64;
	D_CalcGradients(pface, 1);
	D_Lazy(1);
	D_TurbSetup(&d->turb, opacity);
	D_Lazy(0);
	d->kind = DS_TURB;
	d->zbuffer = !r_alphapass;
}

static void D_SetupTurb(dsurf_t *d, surf_t *s, msurface_t *pface, f32 opacity)
{
	miplevel = 0;
	if(!s->entity->model->haslitwater || !r_litwater.value
			|| (pface->flags & SURF_DRAWTILED)){ // unlit
		D_TurbSetLight(&d->turb.light, NULL, NULL, 0);
		D_SetupWater(d, pface, opacity);
		return;
	}
	// the cache holds only light rows, Turbulent applies them per pixel
	miplevel = D_SurfMipLevel(s, pface);
	D_FitLitWaterMip(pface);
	r_lightonly = 1; // tells D_CacheSurface it's caching lit water
	surfcache_t *pcurrentcache = D_SurfBlock(d, pface);
	r_lightonly = 0;
	if(d->kind == DS_LATER)
		return;
	s32 lightmip = miplevel;
	miplevel = 0; // the water texture itself is tiled at full res
	D_TurbSetLight(&d->turb.light, pcurrentcache, pface, lightmip);
	D_SetupWater(d, pface, opacity); // unlit without a block
	if(pcurrentcache)
		D_KeepBlock(d, pcurrentcache, pface, lightmip);
}

// Is what s was set up with still there? A block in the surface cache can be
// taken for another surface, or rebuilt for another frame of its texture.
static bool D_SurfStillValid(const dsurf_t *d)
{
	return !d->cache || (*d->cachespot == d->cache
			&& d->cache->texture == d->texture);
}

enum { PASS_MAIN, PASS_CUTOUT, PASS_SKYBOX, PASS_ALPHA, PASS_FLAT };

// Works out how the spans of s get drawn. For the workers the result stays
// with the surface for the rest of its pass, so the flushes after the first
// only check it, the jobs of other bands too. Without them it's drawn right
// away from the draw state.
static dsurf_t *D_SetupSurf(surf_t *s, s32 pass)
{
	static dsurf_t now;
	dsurf_t *d = &now;
	dsurf_t **slot = NULL;
	s32 i = s - surfaces;
	if(r_threaded){
		slot = &r_surfdsurf[i];
		for(;;){ // a job of another band may be at it
			s32 st = SDL_GetAtomicInt(&r_surfscan[i]);
			if(st == D_BUSY){
				SDL_CPUPauseInstruction();
				continue;
			}
			if((u32)st == r_scan){
				SDL_MemoryBarrierAcquire();
				if(*slot && D_SurfStillValid(*slot))
					return *slot;
			}
			if(SDL_CompareAndSwapAtomicInt(&r_surfscan[i], st,
					D_BUSY))
				break;
		}
		d = R_FrameAlloc(sizeof(*d));
	}
	d->kind = DS_NONE;
	d->cache = NULL;
	d->early = 1;
	d->walk = d_walk;
	msurface_t *pface = s->data;
	if(pass == PASS_FLAT){
		d_npolys++;
		D_SetZGradients(s);
		d->kind = DS_SOLID;
		d->color = (uintptr_t)s->data & 0xFF;
	}else if(pass == PASS_SKYBOX){
		D_SetZGradients(s);
		D_SetupSkybox(d, s);
	}else if(s->flags & SURF_DRAWBACKGROUND){
		d_npolys++;
		d->kind = DS_FAR;
		d->color = (s32)r_clearcolor.value & 0xFF;
	}else if(pass == PASS_ALPHA){
		D_SetZGradients(s);
		D_SwitchEntity(s);
		if(s->flags & SURF_DRAWTURB){
			f32 opacity = 1;
			if(s->entity->alpha && r_entalpha.value)
				opacity -= (f32)s->entity->alpha / 255;
			else if(s->flags & SURF_WINQUAKE_DRAWTRANSLUCENT)
				opacity -= R_LiquidAlphaForFlags(s->flags);
			D_SetupTurb(d, s, pface, opacity);
		}else if(s->entity->alpha && r_entalpha.value == 1)
			D_SetupCached(d, s, pface, SPAN_TRANS,
					1 - (f32)s->entity->alpha / 255);
	}else if(!D_BrokenSurf(s, pface)){
		D_SetZGradients(s);
		D_SwitchEntity(s);
		miplevel = 0;
		if(pass == PASS_CUTOUT)
			D_SetupCached(d, s, pface, SPAN_CUTOUT, 0);
		else if(s->flags & SURF_NOTEXTURE){
			D_Lazy(1);
			Con_DPrintf("D_DrawSurfaces: drawing "
					"SURF_NOTEXTURE\n");
			D_Lazy(0);
			d_npolys++;
			d->kind = DS_SOLID;
			d->color = 0;
		}else if(s->flags & SURF_DRAWSKY){
			if(!skybox_name[0]) D_SetupSky(d, pface);
			else d->kind = DS_SKYZ;
		}else if(s->flags & SURF_DRAWTURB){
			// R_SurfaceInPass kept the translucent ones for later,
			// everything that made it here is opaque
			D_SetupTurb(d, s, pface, 0);
		}else
			D_SetupCached(d, s, pface, SPAN_NORMAL, 0);
	}
	// the sky and translucent surfaces mark what they fog, the marks are
	// there already for the jobs (D_BeginRecords)
	if(d->kind == DS_SKY || d->kind == DS_SKYZ || pass == PASS_ALPHA)
		R_FogMarkAlloc();
	// Kept for the rest of the scan. Without a block it's tried again, as
	// the cache does, and DS_LATER gets set up by the walk.
	if(slot){
		s32 st = 0;
		if(d->kind != DS_LATER){
			*slot = d->kind == DS_NONE ? NULL : d;
			D_SaveGradients(&d->g);
			st = r_scan;
		}
		SDL_MemoryBarrierRelease();
		SDL_SetAtomicInt(&r_surfscan[i], st);
	}
	return d;
}

// A call for each type, so that each gets its own inlined loop
Q_SIL void D_DrawSurfSpans(const dsurf_t *d, espan_t *spans, s32 type)
{
	f32 opacity = type == SPAN_TRANS ? d->opacity : 0;
	if(d->dithered)
		D_DrawSpansDithered(spans, type, opacity);
	else
		D_DrawSpans(spans, type, opacity);
}

// The spans of one surface, on any thread
static void D_DrawSurf(const dsurf_t *d, espan_t *spans, bool load)
{
	if(load)
		D_LoadGradients(&d->g);
	switch(d->kind){
	case DS_SOLID:
		D_DrawSolidSurface(spans, d->color);
		D_DrawZSpans(spans);
		break;
	case DS_FAR:
		D_DrawSolidSurface(spans, d->color);
		D_DrawZSpansFar(spans);
		break;
	case DS_SKY:
		D_DrawSkyScans(spans, d->pface, d->sky);
		D_DrawZSpans(spans);
		R_FogMarkSky(spans);
		break;
	case DS_SKYZ:
		D_DrawZSpans(spans);
		R_FogMarkSky(spans);
		break;
	case DS_TURB:
		Turbulent(spans, &d->turb);
		if(d->zbuffer)
			D_DrawZSpans(spans);
		break;
	case DS_SKYBOX:
		D_DrawSurfSpans(d, spans, SPAN_SKYBOX);
		break;
	case DS_CACHED:
		if(d->type == SPAN_NORMAL){
			D_DrawSurfSpans(d, spans, SPAN_NORMAL);
			D_DrawZSpans(spans);
		}else if(d->type == SPAN_CUTOUT){
			D_DrawSurfSpans(d, spans, SPAN_CUTOUT);
			D_DrawZSpansTrans(spans);
		}else
			D_DrawSurfSpans(d, spans, SPAN_TRANS);
		break;
	}
}

static void D_SurfCmd(void *arg)
{
	const dsurfcmd_t *c = arg;
	D_DrawSurf(c->d, c->spans, 1);
}

static void D_SurfCmdAdd(s32 band, const dsurf_t *d, espan_t *spans)
{
	dsurfcmd_t *c = R_FrameAlloc(sizeof(*c));
	c->d = d;
	c->spans = spans;
	R_BandAdd(band, D_SurfCmd, c);
}

// The surfaces with spans get walked in D_DrawSurfaces' order in one of three
// ways. Without the workers they're set up and drawn. With them the scanline
// flush for a band hands the workers what it can (D_Surf) and keeps the spans
// of the rest. The setup walk at each flush of the span buffer sets those up,
// in the order they'd have without the workers: they fill the surface cache,
// and what it keeps depends on that order (r_lerplightstyles keeps old light,
// the RAM icon, the cache growing). Each band keeps its own records, so its
// walk can be a job (D_RecordBandJob).
typedef struct { u32 surf; espan_t *spans; } srecord_t;
typedef struct {
	srecord_t *r; // the spans the band got since the last setup
	s32 n, cap, next; // the ones up to next are handed over
} brecords_t;
static brecords_t d_records[R_MAXBANDS];

void D_BeginRecords()
{ // in case a frame was cut short
	for(s32 b = 0; b < R_MAXBANDS; b++)
		d_records[b].n = d_records[b].next = 0;
	for(s32 i = 0; i <= R_MAXWORKERS; i++)
		d_found[i].n = 0;
	if(!d_lazy)
		d_lazy = SDL_CreateMutex();
	R_FogMarkAlloc(); // before any job marks with it
}

static void D_RecordSurf(s32 band, surf_t *s, espan_t *spans)
{
	brecords_t *b = &d_records[band];
	if(b->n == b->cap){
		b->cap = b->cap ? b->cap * 2 : 256;
		b->r = realloc(b->r, b->cap * sizeof(*b->r));
		if(!b->r)
			Sys_Error("D_RecordSurf: out of memory");
	}
	b->r[b->n].surf = s - surfaces;
	b->r[b->n++].spans = spans;
}

// Can n be drawn before the setup walk is over? Set up in this scan, with
// what it'd have when the walk comes to it.
static bool D_SetUpEarly(u32 n)
{
	const dsurf_t *d = r_surfdsurf[n];
	return (u32)SDL_GetAtomicInt(&r_surfscan[n]) == r_scan && (!d ||
		(D_SurfStillValid(d) && (d->early || d->walk == d_walk)));
}

// Hands the bands what they recorded, as far as the surfaces are set up: in
// the order of each band, so a band stops at its first surface that isn't
static void D_PublishRecords(bool all)
{
	for(s32 band = 0; band < r_nbands; band++){
		brecords_t *b = &d_records[band];
		for(; b->next < b->n; b->next++){
			srecord_t *r = &b->r[b->next];
			if(!all && !D_SetUpEarly(r->surf))
				break;
			const dsurf_t *d = r_surfdsurf[r->surf];
			if(d && d->kind != DS_NONE)
				D_SurfCmdAdd(band, d, r->spans);
		}
		if(all)
			b->n = b->next = 0;
		R_BandPublish(band);
	}
}

// A block of the surface cache drawn from in this frame is about to change:
// what's recorded from it is drawn first
void D_SurfCacheHazard()
{
	if(!r_threaded)
		return;
	D_PublishRecords(0);
	R_ThreadsSync();
}

// A band's surfaces that need no cache block built go to the workers right
// away. The others wait for the setup walk at the flush, and so does whatever
// has to come after them: the cutouts come after all of the main pass, and the
// alpha pass in its order. The main pass faces don't overlap.
typedef struct { bool waitmain, waitrest; } dwait_t;
static dwait_t d_wait; // of the main thread's walk

static void D_Record(s32 band, surf_t *s, espan_t *spans, s32 pass,
		dwait_t *w)
{
	bool can = pass == PASS_CUTOUT ? !w->waitmain && !w->waitrest
		: pass == PASS_ALPHA ? !w->waitrest : 1;
	const dsurf_t *d = NULL;
	if(can){
		d_nobuild = 1;
		d = D_SetupSurf(s, pass);
		d_nobuild = 0;
	}
	if(d && d->early){
		if(d->kind != DS_NONE && spans)
			D_SurfCmdAdd(band, d, spans);
		return;
	}
	if(pass == PASS_MAIN)
		w->waitmain = 1;
	else
		w->waitrest = 1;
	D_RecordSurf(band, s, spans);
}

Q_SIL void D_Surf(surf_t *s, s32 pass, s32 act)
{
	if(act == SURFS_RECORD){
		D_Record(r_recordband, s, s->spans, pass, &d_wait);
		return;
	}
	dsurf_t *d = D_SetupSurf(s, pass);
	if(act == SURFS_SETUP)
		d->walk = d_walk; // its records can go now
	else if(d->kind != DS_NONE && s->spans)
		D_DrawSurf(d, s->spans, 0);
}

static const u8 debruijn[64] = {
	0, 1, 48, 2, 57, 49, 28, 3, 61, 58, 50, 42, 38, 29, 17, 4,
	62, 55, 59, 36, 53, 51, 43, 22, 45, 39, 33, 30, 24, 18, 12, 5,
	63, 47, 56, 27, 60, 41, 37, 16, 54, 35, 52, 21, 44, 32, 23, 11,
	46, 26, 40, 15, 34, 20, 31, 10, 25, 14, 19, 9, 13, 8, 7, 6 };

// The next surface from i on that the walk is for, or end: the ones with spans
// without the workers, else those a band or the bands since the last setup got
static s32 D_NextSpanSurf(s32 i, s32 end, s32 act)
{
	if(act == SURFS_DRAW){
		const surf_t *s = surfaces;
		while(i < end && !s[i].spans)
			i++;
		return i;
	}
	const u64 *marks = act == SURFS_SETUP ? r_flushmarks : r_surfmarks;
	if(i >= end)
		return end;
	s32 w = i >> 6;
	u64 bits = marks[w] & (~(u64)0 << (i & 63));
	while(!bits){
		if(++w > (end - 1) >> 6)
			return end;
		bits = marks[w];
	}
	i = w * 64 + debruijn[((bits & -bits) * 0x03f79d71b4cb0a89ULL) >> 58];
	return i < end ? i : end;
}

#define FOR_SPAN_SURFS(i, from) \
	for(s32 i = D_NextSpanSurf(from, nsurfs, act); i < nsurfs; \
			i = D_NextSpanSurf(i + 1, nsurfs, act))

static void D_EndSurfaces(s32 nsurfs, s32 act)
{ // start the next flush clean
	s32 words = (nsurfs + 63) >> 6;
	r_drawnpolycount += d_npolys + SDL_SetAtomicInt(&d_jobpolys, 0);
	d_npolys = 0;
	if(act == SURFS_SETUP){
		memset(r_flushmarks, 0, words * sizeof(*r_flushmarks));
		D_PublishRecords(1);
		return;
	}
	FOR_SPAN_SURFS(i, 0)
		surfaces[i].spans = NULL;
	if(act == SURFS_DRAW)
		return;
	for(s32 w = 0; w < words; w++){
		r_flushmarks[w] |= r_surfmarks[w];
		r_surfmarks[w] = 0;
	}
}

// Without spans to look at when setting up, the cutout and skybox surfaces
// are told apart by their flags alone
Q_SIL bool D_Pending(const surf_t *s, s32 act)
{ return act == SURFS_SETUP || s->spans; }

// The per-band state of the surface walks with the workers
static void D_BeginWalk(s32 act)
{
	d_wait.waitmain = d_wait.waitrest = 0;
	d_walk += act == SURFS_SETUP;
}

void D_DrawSurfacesFlat(s32 act)
{
	D_BeginWalk(act);
	if(act != SURFS_RECORD)
		r_drawnpolycount++;
	s32 nsurfs = surface_p - surfaces;
	FOR_SPAN_SURFS(i, 1)
		D_Surf(&surfaces[i], PASS_FLAT, act);
	D_EndSurfaces(nsurfs, act);
}

void D_DrawSurfaces(s32 act)
{
	D_BeginWalk(act);
	D_BeginSurfaces();
	s32 nsurfs = surface_p - surfaces;
	bool skybox = skybox_surf_p != NULL;
	if(skybox){
		FOR_SPAN_SURFS(i, skybox_surf_p - surfaces){
			surf_t *s = &surfaces[i];
			if(D_Pending(s, act) && (s->flags & SURF_DRAWSKYBOX)){
				D_Surf(s, PASS_SKYBOX, act);
				s->spans = NULL;
			}
		}
	}
	bool foundcutouts = 0;
	FOR_SPAN_SURFS(i, 1){
		surf_t *s = &surfaces[i];
		if(!D_Pending(s, act) || (act == SURFS_SETUP && skybox
				&& (s->flags & SURF_DRAWSKYBOX)))
			continue;
		if(D_CutoutPass(s)){
			foundcutouts = 1;
			continue;
		}
		D_Surf(s, PASS_MAIN, act);
		s->spans = NULL;
	}
	if(foundcutouts){
		FOR_SPAN_SURFS(i, 1){
			surf_t *s = &surfaces[i];
			if(D_Pending(s, act) && D_CutoutPass(s))
				D_Surf(s, PASS_CUTOUT, act);
		}
	}
	D_SwitchEntity(NULL);
	D_EndSurfaces(nsurfs, act);
}

// D_DrawSurfaces' record walk of the main pass for one band, on any thread.
// The record walks between two setup walks all go by the same surface cache,
// so they can go on at the same time, each with what its band got: n records
// sorted by surface.
void D_RecordBandJob(void *arg)
{
	const drecjob_t *j = arg;
	const scanrec_t *r = j->recs;
	dwait_t w = {0, 0};
	d_injob = 1;
	D_BeginSurfaces();
	u32 sky0 = skybox_surf_p ? (u32)(skybox_surf_p - surfaces) : ~0u;
	bool cutouts = 0;
	for(s32 i = 0; i < j->n; i++){
		surf_t *s = &surfaces[r[i].surf];
		if(r[i].surf >= sky0 && (s->flags & SURF_DRAWSKYBOX))
			D_Record(j->band, s, r[i].spans, PASS_SKYBOX, &w);
	}
	for(s32 i = 0; i < j->n; i++){
		surf_t *s = &surfaces[r[i].surf];
		if(r[i].surf >= sky0 && (s->flags & SURF_DRAWSKYBOX))
			continue;
		if(D_CutoutPass(s))
			cutouts = 1;
		else
			D_Record(j->band, s, r[i].spans, PASS_MAIN, &w);
	}
	for(s32 i = 0; cutouts && i < j->n; i++){
		surf_t *s = &surfaces[r[i].surf];
		if(D_CutoutPass(s))
			D_Record(j->band, s, r[i].spans, PASS_CUTOUT, &w);
	}
	R_BandPublish(j->band);
	SDL_AddAtomicInt(&d_jobpolys, d_npolys);
	d_npolys = 0;
	d_injob = 0;
}

// The record walks are done: the blocks they found are used in this frame
void D_StampBlocks()
{
	for(s32 i = 0; i <= R_MAXWORKERS; i++){
		for(s32 k = 0; k < d_found[i].n; k++)
			d_found[i].c[k]->frame = r_framecount & 0x7fffffff;
		d_found[i].n = 0;
	}
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

void D_DrawSurfacesAlpha(s32 act)
{
	D_BeginWalk(act);
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
	s32 count = 0, nsurfs = surface_p - surfaces;
	FOR_SPAN_SURFS(i, 1)
		order[count++] = &surfaces[i];
	qsort(order, count, sizeof(*order), D_AlphaSurfCmp);
	for(s32 i = 0; i < count; i++)
		D_Surf(order[i], PASS_ALPHA, act);
	D_SwitchEntity(NULL);
	D_EndSurfaces(nsurfs, act);
}
