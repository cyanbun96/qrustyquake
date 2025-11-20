// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
#include "quakedef.h"

static mem_journal_t *journal_head = 0;
static mem_journal_t *journal_tail = 0;

void *Q_Malloc(u64 size, cache_user_t *cache_user, s8 *name)
{
	if(sv_cheats.value == 8) Con_DPrintf("malloc %d\n", size);
	mem_journal_t *entry = malloc(sizeof(mem_journal_t));
	if(journal_tail) {
		journal_tail->next = entry;
	}
	journal_tail = entry;
	if(!journal_head) journal_head = entry;
	void *ptr = malloc(size);
	entry->addr = (u64)ptr;
	entry->size = size;
	entry->cache_user = cache_user;
	Q_strncpy(entry->name, name?name:"Unknown", 15);
	entry->next = 0;
	return ptr;
}

void Q_Free(void *ptr)
{
	if(sv_cheats.value == 8) Con_DPrintf("free\n");
	mem_journal_t *entry = 0;
	mem_journal_t *last = 0;
	for(mem_journal_t *search = journal_head; search; search =search->next){
		if(search->addr == (u64)ptr) {
			entry = search;
			break;
		}
		last = search;
	}
	if(!entry) Sys_Error("Q_Free: called on an unknown allocation");
	free((void*)(entry->addr));
	if(entry == journal_tail) {
		journal_tail = last;
		last->next = 0;
	} else last->next = entry->next;
	free(entry);
}

void *Q_Realloc(void *ptr, u64 size, cache_user_t *cache_user, s8 *name)
{
	//if(sv_cheats.value == 9) __asm__("int3");
	//if(sv_cheats.value == 8) Con_DPrintf("realloc\n");
	if(!ptr) return Q_Malloc(size, cache_user, name);
	mem_journal_t *entry = 0;
	for(mem_journal_t *search = journal_head; search; search =search->next){
		if(search->addr == (u64)ptr) {
			entry = search;
			break;
		}
	}
	if(!entry) Sys_Error("Q_Realloc: called on an unknown allocation");
	if(entry->size == size) return ptr;
	ptr = realloc(ptr, size);
	entry->addr = (u64)ptr;
	entry->size = size;
	entry->cache_user = cache_user;
	if(name)Q_strncpy(entry->name, name, 15);
	return ptr;
}

void Mem_Journal_Show() {
	u64 tot = 0, mallocs = 0;
	for(mem_journal_t *entry = journal_head; entry; entry = entry->next) {
		Con_Printf("%x %u %x %s\n", entry->addr, entry->size,
				entry->cache_user, entry->name);
		tot += entry->size;
		mallocs++;
	}
	Con_Printf("Total: %d.%d.%d\n", tot/1000000, (tot/1000)%1000, tot%1000);
	Con_Printf("Mallocs: %d\n", mallocs);
}

void Z_Free(void *ptr) { Q_Free(ptr); }

void *Z_Malloc(s32 size)
{
	void *ptr = Q_Malloc(size, 0, 0);
	memset(ptr, 0, size);
	return ptr;
}

void *Z_Realloc(void *ptr, s32 size)
{ //FIXME this does not zero out the new memory
	return Q_Realloc(ptr, size, 0, 0);
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

void *Hunk_AllocName(s32 size, s8 *name)
{
	size = sizeof(hunk_t) + ((size+15)&~15);
	void *p = Q_Malloc(size, 0, name);
	memset(p, 0, size);
	return p;
}

void Hunk_Check() { return; /*check passed ok!*/ }

void *Hunk_Alloc(s32 size)
{ return Hunk_AllocName(size, "unknown"); }

void Hunk_FreeToLowMark(s32 mark)
{
	Con_DPrintf("NOT FREEING (LOWMARK)\n");
	//if(sv_cheats.value == 9) __asm__("int3");
}

void Hunk_FreeToHighMark(SDL_UNUSED s32 mark)
{
	Con_DPrintf("NOT FREEING (HIGHMARK)\n");
}

s32 Hunk_LowMark() { Con_DPrintf("Hunk_LowMark DEPRECATED\n"); return 0; }

s32 Hunk_HighMark() { Con_DPrintf("Hunk_HighMark DEPRECATED\n"); return 0; }

void *Hunk_HighAllocName(SDL_UNUSED s32 size, SDL_UNUSED s8 *name) { }
void Cache_Report() { }

void *Hunk_TempAlloc(s32 size)
{
	size = (size+15)&~15;
	return Q_Malloc(size, 0, 0);
}

s8 *Hunk_Strdup(const s8 *s, s8 *name)
{
	size_t sz = strlen(s) + 1;
	s8 *ptr = (s8 *) Hunk_AllocName(sz, name);
	memcpy(ptr, s, sz);
	return ptr;
}

void Cache_Free(cache_user_t *c, SDL_UNUSED bool freetextures)
{
	if(!c->data) Sys_Error("Cache_Free: not allocated");
	Q_Free(c->data);
	c->data = NULL;
}

void Cache_Flush()
{ // Throw everything out, so new data will be demand cached
	mem_journal_t *last = 0;
	for(mem_journal_t *entry = journal_head; entry;) {
		mem_journal_t *next = entry->next;
		if(entry->cache_user) {
			Cache_Free(entry->cache_user, 0);
		}
		entry = next;
	}
}

void *Cache_Check(cache_user_t *c) { return c->data; }

void *Cache_Alloc(cache_user_t *c, s32 size, s8 *name)
{
	if(c->data) Sys_Error("Cache_Alloc: already allocated");
	c->data = Q_Malloc(size, c, name);
	return c->data;
}
