// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// d_surf.c: rasterization driver surface heap manager
#include "quakedef.h"
#include <stddef.h> // Required for offsetof

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
{
    // CyanBun96: reallocs the cache if it's too small
    s32 chunk = vid.width * vid.height * sizeof(*d_pzbuffer);
    s32 cachesize = d_pzbuffer_size + D_SurfaceCacheForRes(vid.width, vid.height);
    chunk += cachesize;

    D_FlushCaches(0);

    // FIX: Standard safe realloc pattern
    void *tmp_ptr = realloc(d_pzbuffer, chunk);
    if (tmp_ptr == NULL) {
        Sys_Error("Not enough memory for surf cache\n");
    }
    d_pzbuffer = tmp_ptr;

    u8 *cache = (u8*)d_pzbuffer + vid.width * vid.height * sizeof(*d_pzbuffer);
    D_InitCaches(cache, cachesize);
    vid.recalc_refdef = 1;
}

surfcache_t *D_SCAlloc(s32 width, uintptr_t data_size)
{
    // 1. Sanity Checks
    if ((width < 0) || (width > 256)) {
        Con_DPrintf("D_SCAlloc: bad cache width %d\n", width);
        return NULL;
    }
    if ((data_size <= 0) || (data_size > 0x10000)) {
        Con_DPrintf("D_SCAlloc: bad cache size %d\n", data_size);
        return NULL;
    }

    // 2. Calculate total block size (Header + Data + Alignment)
    // FIX: Use offsetof for clean calculation of struct overhead
    uintptr_t total_size = offsetof(surfcache_t, data) + data_size;
    total_size = (total_size + 3) & ~3; // Align to 4 bytes

    if (total_size > sc_size) {
        Con_DPrintf("D_SCAlloc: %i > cache size", total_size);
        return NULL;
    }

    // 3. Rover Management (Ring Buffer Logic)
    // If there is not enough space after the rover, wrap to the start.
    bool wrapped_this_time = 0;
    if (!sc_rover || (u64)((u8 *) sc_rover - (u8 *) sc_base) > sc_size - total_size) {
        if (sc_rover) {
            wrapped_this_time = 1;
        }
        sc_rover = sc_base;
    }

    // 4. Allocation Loop
    // Collect and free surfcache_t blocks until the rover block is large enough
    // FIX: Renamed 'new' to 'block' for C++ compatibility
    surfcache_t *block = sc_rover;
    
    if (sc_rover->owner)
        *sc_rover->owner = NULL;

    while (block->size < (s32)total_size) {
        sc_rover = sc_rover->next; // free the next block in line
        if (!sc_rover)
            Sys_Error("D_SCAlloc: hit the end of memory");
        
        if (sc_rover->owner)
            *sc_rover->owner = NULL;

        block->size += sc_rover->size;
        block->next = sc_rover->next;
    }

    // 5. Fragmentation
    // If the found block is significantly larger than needed, split it.
    if (block->size - total_size > 256) { 
        sc_rover = (surfcache_t *) ((u8 *) block + total_size);
        sc_rover->size = block->size - total_size;
        sc_rover->next = block->next;
        sc_rover->width = 0; // 0 width marks this as an empty fragment
        sc_rover->owner = NULL;
        
        block->next = sc_rover;
        block->size = total_size;
    } else {
        sc_rover = block->next;
    }

    // 6. Finalize Block
    block->width = width;
    block->owner = NULL; // should be set properly by caller

    // FIX: Cleaned up the "DEBUG" calculation.
    // We only calculate height if width > 0 (to avoid Div/0 on fragments).
    // Using offsetof ensures we accurately calculate available data payload.
    if (width > 0) {
        size_t available_payload = block->size - offsetof(surfcache_t, data);
        block->height = available_payload / width;
    } else {
        block->height = 0;
    }

    // 7. Thrashing Detection
    if (d_roverwrapped && (wrapped_this_time || sc_rover >= d_initial_rover))
        r_cache_thrash = 1;
    else if (wrapped_this_time)
        d_roverwrapped = 1;

    return block;
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
