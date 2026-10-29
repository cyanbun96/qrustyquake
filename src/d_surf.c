// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// d_surf.c: rasterization driver surface heap manager
#include "quakedef.h"

// A cache that gets recycled faster than this many frames rebuilds most
// surfaces before they get a chance to be reused, so it grows
#define SC_MINTURNOVER 30
#define SC_MAXSIZE (256 * 1024 * 1024)

static u64 sc_size;
static surfcache_t *sc_base;
static s64 sc_wrapframe; // when the rover last went back to the start
s32 D_SurfaceCacheForRes(s32 width, s32 height)
{
	s32 size;
	if(COM_CheckParm("-surfcachesize")){
		size = atoi(com_argv[COM_CheckParm("-surfcachesize")+1])*1024;
		return size;
	}
	size = 600*1024; // surfcache size at 320x200
	s32 pix = width * height;
	if(pix > 64000)
		size += (pix - 64000) * 3;
	return size;
}

void D_ClearCacheGuard()
{
	u8 *s = (u8*)sc_base + sc_size;
	for(s32 i = 0; i < GUARDSIZE; i++)
		s[i] = (u8)i;
}

void D_InitCaches(void *buffer, s32 size)
{
	Con_DPrintf("%ik surface cache\n", size / 1024);
	sc_size = size - GUARDSIZE;
	sc_base = (surfcache_t *) buffer;
	sc_rover = sc_base;
	sc_base->next = NULL;
	sc_base->owner = NULL;
	sc_base->size = sc_size;
	sc_wrapframe = r_framecount;
	D_ClearCacheGuard();
	d_pzbuffer_size = size;
}

void D_FlushCaches(SDL_UNUSED cvar_t *cvar)
{
	if(!sc_base)
		return;
	for(surfcache_t *c = sc_base; c; c = c->next)
		if(c->owner){
			uintptr_t owner_addr = (uintptr_t)c->owner;
			if(owner_addr < (uintptr_t)host_parms.membase ||
					owner_addr >= (uintptr_t)((uintptr_t)
					host_parms.membase+host_parms.memsize)){
			  Con_DPrintf("D_FlushCaches: corrupt owner pointer\n");
			}
			*c->owner = NULL;
		}
	sc_rover = sc_base;
	sc_base->next = NULL;
	sc_base->owner = NULL;
	sc_base->size = sc_size;
	sc_wrapframe = r_framecount;
}

void D_AllocCaches() // allocate z buffer and surface cache
{ // CyanBun96: reallocs the cache if it's too small
	s32 chunk = vid.width * vid.height * sizeof(*d_pzbuffer);
	s32 cachesize = D_SurfaceCacheForRes(vid.width, vid.height);
	if(d_pzbuffer_size){ // thrashing, double it
		if(d_pzbuffer_size >= SC_MAXSIZE)
			return;
		cachesize = d_pzbuffer_size * 2;
		if(cachesize > SC_MAXSIZE)
			cachesize = SC_MAXSIZE;
	}
	chunk += cachesize;
	D_FlushCaches(0);
	d_pzbuffer = realloc(d_pzbuffer, chunk);
	memset(d_pzbuffer, 0, chunk);
	if(d_pzbuffer == NULL) Sys_Error("Not enough memory for surf cache\n");
	u8 *cache = (u8*)d_pzbuffer+vid.width*vid.height*sizeof(*d_pzbuffer);
	D_InitCaches(cache, cachesize);
	vid.recalc_refdef = 1;
}

surfcache_t *D_SCAlloc(s32 width, uintptr_t size)
{
	if((width < 0) || (width > 256)){
		Con_DPrintf("D_SCAlloc: bad cache width %d\n", width);
		return NULL;
	}
	if((size <= 0) || (size > 0x10000)){
		Con_DPrintf("D_SCAlloc: bad cache size %d\n", size);
		return NULL;
	}
	size = (uintptr_t) & ((surfcache_t *) 0)->data[size];
	size = (size + 3) & ~3;
	if(size > sc_size){
		Con_DPrintf("D_SCAlloc: %i > cache size", size);
		return NULL;
	}
	// if there is not size bytes after the rover, reset to the start
	bool wrapped_this_time = 0;
	if(!sc_rover || (u64)((u8 *) sc_rover - (u8 *) sc_base)
			> sc_size - size){
		if(sc_rover){
			wrapped_this_time = 1;
			if(r_framecount - sc_wrapframe < SC_MINTURNOVER)
				r_cache_thrash = 1;
			sc_wrapframe = r_framecount;
		}
		sc_rover = sc_base;
	}
	// colect and free surfcache_t blocks until rover block is large enough
	surfcache_t *new = sc_rover;
	if(sc_rover->owner)
		*sc_rover->owner = NULL;
	uintptr_t cache_lo = (uintptr_t)sc_base;
	uintptr_t cache_hi = cache_lo + sc_size;
	while(new->size < (s32)size){
		sc_rover = sc_rover->next; // free another
		if(!sc_rover){
			Con_DPrintf("D_SCAlloc: hit the end of memory\n");
			return NULL;
		}
		uintptr_t rover_addr = (uintptr_t)sc_rover;
		if(rover_addr < cache_lo || rover_addr >= cache_hi){
			Con_DPrintf("D_SCAlloc: corrupt sc_rover\n");
			return NULL;
		}
		if(sc_rover->size <= 0 || (u64)sc_rover->size > sc_size){
			Con_DPrintf("D_SCAlloc: corrupt block size\n");
			return NULL;
		}
		if(sc_rover->next){ // validate next
			uintptr_t next_addr = (uintptr_t)sc_rover->next;
			if(next_addr < cache_lo || next_addr >= cache_hi){
			    Con_DPrintf("D_SCAlloc: corrupt next pointer\n");
				return NULL;
			}
		}
		if(sc_rover->owner){//validate owner(against memory, not cache)
			uintptr_t owner_addr = (uintptr_t)sc_rover->owner;
			if(owner_addr < (uintptr_t)host_parms.membase ||
					owner_addr >= (uintptr_t)((uintptr_t)
					host_parms.membase+host_parms.memsize)){
			    Con_DPrintf("D_SCAlloc: corrupt owner pointer\n");
				return NULL;
			}
			*sc_rover->owner = NULL;
		}
		new->size += sc_rover->size;
		new->next = sc_rover->next;
	}
	if(new->size - size > 256){ // create a fragment out of any leftovers
		sc_rover = (surfcache_t *) ((u8 *) new + size);
		sc_rover->size = new->size - size;
		sc_rover->next = new->next;
		sc_rover->width = 0;
		sc_rover->owner = NULL;
		new->next = sc_rover;
		new->size = size;
	}else
		sc_rover = new->next;
	new->width = width;
	if(width > 0) // DEBUG
		new->height = (size - sizeof(*new) + sizeof(new->data)) / width;
	new->owner = NULL; // should be set properly after return
	if(d_roverwrapped && (wrapped_this_time || sc_rover >=d_initial_rover))
			r_cache_thrash = 1;
	else if(wrapped_this_time)
		d_roverwrapped = 1;
	return new;
}

static s32 D_LightAdjMatches(msurface_t *surf, surfcache_t *cache)
{ // lerped lightstyles change every frame; keep the cached surface while
  // the summed light error stays under one colormap shade at its peak
	s32 *a = cache->lightadj, *b = r_drawsurf.lightadj;
	if(!r_lerplightstyles.value)
		return a[0]==b[0] && a[1]==b[1] && a[2]==b[2] && a[3]==b[3];
	s32 err = 0;
	for(s32 m = 0; m < MAXLIGHTMAPS; m++)
		err += abs(a[m] - b[m]) * surf->lightpeak[m];
	return err < 1 << (16 - VID_CBITS); // one shade of blocklight
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
	if(cache && !cache->dlight && surface->dlightframe != r_framecount
	    && cache->texture == r_drawsurf.texture
	    && D_LightAdjMatches(surface, cache))
		return cache;
	static f32 surfscale;
	surfscale = 1.0 / (1 << miplevel); // determine shape of surface
	r_drawsurf.surfmip = miplevel;
	r_drawsurf.surfwidth = surface->extents[0] >> miplevel;
	r_drawsurf.rowbytes = r_drawsurf.surfwidth;
	if(r_lightonly && R_LitWaterRGB())
		r_drawsurf.rowbytes *= 4;
	r_drawsurf.surfheight = surface->extents[1] >> miplevel;
	// allocate memory if needed
	if(!cache){ // if a texture just animated, don't reallocate it
		cache = D_SCAlloc(r_drawsurf.surfwidth,
				  r_drawsurf.rowbytes * r_drawsurf.surfheight);
		if(cache == NULL) return NULL;
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
	if(R_DrawSurface()) return NULL;
	return surface->cachespots[miplevel];
}
