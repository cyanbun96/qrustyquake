// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// d_surf.c: rasterization driver surface heap manager
#include "quakedef.h"

static f32 surfscale;                                           
static u64 sc_size;
static surfcache_t *sc_base;
s32 D_SurfaceCacheForRes(s32 width, s32 height)
{
	s32 size;
	if (COM_CheckParm("-surfcachesize")) {
		size = Q_atoi(com_argv[COM_CheckParm("-surfcachesize")+1])*1024;
		return size;
	}
	size = SURFCACHE_SIZE_AT_320X200;
	s32 pix = width * height;
	if (pix > 64000)
		size += (pix - 64000) * 3;
	return size;
}

void D_ClearCacheGuard()
{
	u8 *s = (u8 *) sc_base + sc_size;
	for (s32 i = 0; i < GUARDSIZE; i++)
		s[i] = (u8) i;
}

void D_InitCaches(void *buffer, s32 size)
{
	Con_Printf("%ik surface cache\n", size / 1024);
	sc_size = size - GUARDSIZE;
	sc_base = (surfcache_t *) buffer;
	sc_rover = sc_base;
	sc_base->next = NULL;
	sc_base->owner = NULL;
	sc_base->size = sc_size;
	D_ClearCacheGuard();
	d_pzbuffer_size = size;
}

void D_FlushCaches(SDL_UNUSED cvar_t *cvar)
{
	surfcache_t *c;
	if (!sc_base)
		return;
	for (c = sc_base; c; c = c->next)
		if (c->owner)
			*c->owner = NULL;
	sc_rover = sc_base;
	sc_base->next = NULL;
	sc_base->owner = NULL;
	sc_base->size = sc_size;
}

void D_AllocCaches() // allocate z buffer and surface cache
{ // CyanBun96: reallocs the cache if it's too small
	s32 chunk = vid.width * vid.height * sizeof(*d_pzbuffer);
	s32 cachesize = d_pzbuffer_size + 
		D_SurfaceCacheForRes(vid.width, vid.height);
	chunk += cachesize;
	D_FlushCaches(0);
	d_pzbuffer = realloc(d_pzbuffer, chunk);
	if(d_pzbuffer == NULL) Sys_Error("Not enough memory for surf cache\n");
	u8 *cache = (u8*)d_pzbuffer+vid.width*vid.height*sizeof(*d_pzbuffer);
	D_InitCaches(cache, cachesize);
}

surfcache_t *D_SCAlloc(s32 width, uintptr_t size)
{
	if ((width < 0) || (width > 256)) {
		Con_DPrintf("D_SCAlloc: bad cache width %d\n", width);
		return NULL;
	}
	if ((size <= 0) || (size > 0x10000)) {
		Con_DPrintf("D_SCAlloc: bad cache size %d\n", size);
		return NULL;
	}
	size = (uintptr_t) & ((surfcache_t *) 0)->data[size];
	size = (size + 3) & ~3;
	if (size > sc_size) {
		Con_DPrintf("D_SCAlloc: %i > cache size", size);
		return NULL;
	}
	// if there is not size bytes after the rover, reset to the start
	bool wrapped_this_time = 0;
	if (!sc_rover || (u64)((u8 *) sc_rover - (u8 *) sc_base)
			> sc_size - size) {
		if (sc_rover) {
			wrapped_this_time = 1;
		}
		sc_rover = sc_base;
	}
	// colect and free surfcache_t blocks until rover block is large enough
	surfcache_t *new = sc_rover;
	if (sc_rover->owner)
		*sc_rover->owner = NULL;
	while (new->size < (s32)size) {
		sc_rover = sc_rover->next; // free another
		if (!sc_rover)
			Sys_Error("D_SCAlloc: hit the end of memory");
		if (sc_rover->owner)
			*sc_rover->owner = NULL;

		new->size += sc_rover->size;
		new->next = sc_rover->next;
	}
	if (new->size - size > 256) { // create a fragment out of any leftovers
		sc_rover = (surfcache_t *) ((u8 *) new + size);
		sc_rover->size = new->size - size;
		sc_rover->next = new->next;
		sc_rover->width = 0;
		sc_rover->owner = NULL;
		new->next = sc_rover;
		new->size = size;
	} else
		sc_rover = new->next;
	new->width = width;
	if (width > 0) // DEBUG
		new->height = (size - sizeof(*new) + sizeof(new->data)) / width;
	new->owner = NULL; // should be set properly after return
	if (d_roverwrapped && (wrapped_this_time || sc_rover >=d_initial_rover))
			r_cache_thrash = 1;
	else if (wrapped_this_time)
		d_roverwrapped = 1;
	return new;
}

surfcache_t *D_CacheSurface(msurface_t *surface, s32 miplevel)
{
	// if the surface is animating or flashing, flush the cache
	r_drawsurf.texture = R_TextureAnimation(surface->texinfo->texture);
	r_drawsurf.lightadj[0] = d_lightstylevalue[surface->styles[0]];
	r_drawsurf.lightadj[1] = d_lightstylevalue[surface->styles[1]];
	r_drawsurf.lightadj[2] = d_lightstylevalue[surface->styles[2]];
	r_drawsurf.lightadj[3] = d_lightstylevalue[surface->styles[3]];
	// see if the cache holds apropriate data
	surfcache_t *cache = surface->cachespots[miplevel];
	if (cache && !cache->dlight && surface->dlightframe != r_framecount
	    && cache->texture == r_drawsurf.texture
	    && cache->lightadj[0] == r_drawsurf.lightadj[0]
	    && cache->lightadj[1] == r_drawsurf.lightadj[1]
	    && cache->lightadj[2] == r_drawsurf.lightadj[2]
	    && cache->lightadj[3] == r_drawsurf.lightadj[3])
		return cache;
	surfscale = 1.0 / (1 << miplevel); // determine shape of surface
	r_drawsurf.surfmip = miplevel;
	r_drawsurf.surfwidth = surface->extents[0] >> miplevel;
	r_drawsurf.rowbytes = r_drawsurf.surfwidth;
	r_drawsurf.surfheight = surface->extents[1] >> miplevel;
	// allocate memory if needed
	if (!cache) // if a texture just animated, don't reallocate it
	{
		cache = D_SCAlloc(r_drawsurf.surfwidth,
				  r_drawsurf.surfwidth * r_drawsurf.surfheight);
		if (cache == NULL) return NULL;
		surface->cachespots[miplevel] = cache;
		cache->owner = &surface->cachespots[miplevel];
		cache->mipscale = surfscale;
	}
	cache->dlight = (surface->dlightframe == r_framecount);
	r_drawsurf.surfdat = (u8 *) cache->data;
	cache->texture = r_drawsurf.texture;
	cache->lightadj[0] = r_drawsurf.lightadj[0];
	cache->lightadj[1] = r_drawsurf.lightadj[1];
	cache->lightadj[2] = r_drawsurf.lightadj[2];
	cache->lightadj[3] = r_drawsurf.lightadj[3];
	r_drawsurf.surf = surface; // draw and light the surface texture
	c_surf++;
	if (R_DrawSurface()) return NULL;
	return surface->cachespots[miplevel];
}
