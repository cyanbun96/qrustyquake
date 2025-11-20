// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
#include "quakedef.h"

void *Q_Malloc(u64 size, cache_user_t *cache_user, s32 type, s8 *name)
{
	if(sv_cheats.value == 8) Con_DPrintf("malloc %d\n", size);
	mem_journal_t *entry = malloc(sizeof(mem_journal_t));
	if(journal_tail) {
		journal_tail->next = entry;
	}
	journal_tail = entry;
	if(!journal_head) journal_head = entry;
	void *ptr = malloc(size);
	memset(ptr, 0, size);
	entry->addr = (u64)ptr;
	entry->size = size;
	entry->cache_user = cache_user;
	entry->type = type;
	Q_strncpy(entry->name, name?name:"Unknown", 15);
	entry->name[15] = 0;
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

void *Q_Realloc(void *ptr, u64 size, cache_user_t *cache_user, s32 type,s8*name)
{
	//if(sv_cheats.value == 9) __asm__("int3");
	//if(sv_cheats.value == 8) Con_DPrintf("realloc\n");
	if(!ptr) return Q_Malloc(size, cache_user, type, name);
	mem_journal_t *entry = 0;
	for(mem_journal_t *search = journal_head; search; search =search->next){
		if(search->addr == (u64)ptr) {
			entry = search;
			break;
		}
	}
	if(!entry) return Q_Malloc(size, cache_user, type, name);
	if(entry->size == size) return ptr;
	ptr = realloc(ptr, size);
	if(entry->size < size) memset(ptr + entry->size, 0, size - entry->size);
	entry->addr = (u64)ptr;
	entry->size = size;
	entry->cache_user = cache_user;
	if(name)Q_strncpy(entry->name, name, 15);
	return ptr;
}

void Mem_Journal_Show() {
	u64 tot = 0, mallocs = 0;
	for(mem_journal_t *entry = journal_head; entry; entry = entry->next) {
		Con_Printf("%x %u %x %d %s\n", entry->addr, entry->size,
				entry->cache_user, entry->type, entry->name);
		tot += entry->size;
		mallocs++;
	}
	Con_Printf("Total: %d.%d.%d\n", tot/1000000, (tot/1000)%1000, tot%1000);
	Con_Printf("Mallocs: %d\n", mallocs);
}

void *Z_Realloc(void *ptr, s32 size)
{
	return Q_Realloc(ptr, size, 0, 0, 0);
	if(!ptr) return Q_Malloc(size, 0, 0, "z_realloc");
}

s8 *Z_Strdup(const s8 *s)
{
	size_t sz = strlen(s) + 1;
	s8 *ptr = (s8 *) Q_Malloc(sz, 0, 0, "z_strdup");
	memcpy(ptr, s, sz);
	return ptr;
}

void Hunk_Print_f() { return; }

void *Hunk_AllocName(s32 size, s8 *name)
{
	size = sizeof(hunk_t) + ((size+15)&~15);
	void *p = Q_Malloc(size, 0, 1, name);
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

s32 Hunk_LowMark() { Con_DPrintf("Hunk_LowMark DEPRECATED\n"); return 0; }

void Cache_Report() { }

void *Hunk_TempAlloc(s32 size)
{
	size = (size+15)&~15;
	return Q_Malloc(size, 0, 1, 0);
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
	c->data = Q_Malloc(size, c, 1, name);
	return c->data;
}
