// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
#include "quakedef.h"

void Z_Free(void *ptr) { free(ptr); }

void *Z_Malloc(s32 size)
{
	void *ptr = malloc(size);
	memset(ptr, 0, size);
	return ptr;
}

void *Z_Realloc(void *ptr, s32 size)
{ //FIXME this does not zero out the new memory
	return realloc(ptr, size);
	if(!ptr) return Z_Malloc(size);
}

s8 *Z_Strdup(const s8 *s)
{
	size_t sz = strlen(s) + 1;
	s8 *ptr = (s8 *) Z_Malloc(sz);
	memcpy(ptr, s, sz);
	return ptr;
}

void Hunk_Print_f() { return; }

void *Hunk_AllocName(s32 size, const s8 *name)
{
	size = sizeof(hunk_t) + ((size+15)&~15);
	void *p = malloc(size);
	memset(p, 0, size);
	return p;
}

void Hunk_Check() { return; /*check passed ok!*/ }

void *Hunk_Alloc(s32 size)
{ return Hunk_AllocName(size, "unknown"); }

s32 Hunk_LowMark() { return 0; }

void Hunk_FreeToLowMark(s32 mark) { }
s32 Hunk_HighMark() { return 0; }
void Hunk_FreeToHighMark(SDL_UNUSED s32 mark) { }
void *Hunk_HighAllocName(SDL_UNUSED s32 size, SDL_UNUSED const s8 *name) { }
void Cache_Report() { }

void *Hunk_TempAlloc(s32 size)
{
	size = (size+15)&~15;
	return malloc(size);
}

s8 *Hunk_Strdup(const s8 *s, const s8 *name)
{
	size_t sz = strlen(s) + 1;
	s8 *ptr = (s8 *) Hunk_AllocName(sz, name);
	memcpy(ptr, s, sz);
	return ptr;
}

cache_system_t *Cache_TryAlloc(s32 size, SDL_UNUSED bool nobottom)
{
	cache_system_t *p = malloc(size);
	if(!p) Sys_Error("Cache_TryAlloc: failed to allocate %d\n", size);
	return p;
}

void Cache_Flush()
{ // Throw everything out, so new data will be demand cached
	return; //FIXME this is particularly bad
}

void Cache_Init()
{
	Cmd_AddCommand("flush", Cache_Flush);
}

void Cache_Free(cache_user_t *c, SDL_UNUSED bool freetextures)
{ // Frees the memory
	if(!c->data) Sys_Error("Cache_Free: not allocated");
	free(c->data);
	c->data = NULL;
}

void *Cache_Check(cache_user_t *c) { return c->data; }

void *Cache_Alloc(cache_user_t *c, s32 size, SDL_UNUSED const s8 *name)
{
	if(c->data) Sys_Error("Cache_Alloc: already allocated");
	c->data = malloc(size);
	return c->data;
}
