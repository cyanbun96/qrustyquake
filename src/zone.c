// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
#include "quakedef.h"

// ZONE MEMORY ALLOCATION
// There is never any space between memblocks, and there will never be two
// contiguous free memblocks.
// The rover can be left pointing at a non-empty block
// The zone calls are pretty much only used for small strings and structures,
// all big things are allocated on the hunk.

static memzone_t *mainzone;
static u8 *hunk_base;
static s32 hunk_size;
static s32 hunk_low_used;
static s32 hunk_high_used;
static bool hunk_tempactive;
static s32 hunk_tempmark;
static cache_system_t cache_head;

void Cache_FreeLow(s32 new_low_hunk);
void Cache_FreeHigh(s32 new_high_hunk);
cache_system_t *Cache_TryAlloc(s32 size, bool nobottom);
void Cache_Free(cache_user_t *c, bool freetextures);

void Z_Free(void *ptr)
{
	if(!ptr) Sys_Error("Z_Free: NULL pointer");
	memblock_t *block = (memblock_t *) ( (u8 *)ptr - sizeof(memblock_t));
	if(block->id != ZONEID)
		Sys_Error("Z_Free: freed a pointer without ZONEID");
	if(block->tag == 0)
		Sys_Error("Z_Free: freed a freed pointer");
	block->tag = 0; // mark as free
	memblock_t *other = block->prev;
	if(!other->tag) { // merge with previous free block
		other->size += block->size;
		other->next = block->next;
		other->next->prev = other;
		if(block == mainzone->rover)
			mainzone->rover = other;
		block = other;
	}
	other = block->next;
	if(!other->tag) { // merge the next free block onto the end
		block->size += other->size;
		block->next = other->next;
		block->next->prev = block;
		if(other == mainzone->rover)
			mainzone->rover = block;
	}
}


static void *Z_TagMalloc(s32 size, s32 tag)
{
	memblock_t *start, *rover, *newblock, *base;
	if(!tag) Sys_Error("Z_TagMalloc: tried to use a 0 tag");
	// scan through the block list looking for the first free block
	// of sufficient size
	size += sizeof(memblock_t); // account for size of block header
	size += 4; // space for memory trash tester
	size = (size + 7) & ~7; // align to 8-u8 boundary
	base = rover = mainzone->rover;
	start = base->prev;
	do {
		if(rover == start)return 0; //scaned all the way around the list
		if(rover->tag) base = rover = rover->next;
		else rover = rover->next;
	} while(base->tag || base->size < size);
	// found a block big enough
	s32 extra = base->size - size;
	if(extra > MINFRAGMENT)
	{ // there will be a free fragment after the allocated block
		newblock = (memblock_t *) ((u8 *)base + size );
		newblock->size = extra;
		newblock->tag = 0; // free block
		newblock->prev = base;
		newblock->id = ZONEID;
		newblock->next = base->next;
		newblock->next->prev = newblock;
		base->next = newblock;
		base->size = size;
	}
	base->tag = tag; // no longer a free block
	mainzone->rover = base->next; // next allocation will start looking here
	base->id = ZONEID;
	// marker for memory trash testing
	*(s32 *)((u8 *)base + base->size - 4) = ZONEID;
	return(void *) ((u8 *)base + sizeof(memblock_t));
}

static void Z_CheckHeap()
{
	for(memblock_t *block = mainzone->blocklist.next;; block = block->next){
		if(block->next == &mainzone->blocklist)
			break; // all blocks have been hit
		if((u8 *)block + block->size != (u8 *)block->next)
	Sys_Error("Z_CheckHeap: block size does not touch the next block");
		if(block->next->prev != block)
	Sys_Error("Z_CheckHeap: next block doesn't have proper back link");
		if(!block->tag && !block->next->tag)
			Sys_Error("Z_CheckHeap: two consecutive free blocks");
	}
}

void *Z_Malloc(s32 size)
{
	Z_CheckHeap(); // DEBUG
	void *buf = Z_TagMalloc(size, 1);
	if(!buf) Sys_Error("Z_Malloc: failed on allocation of %i bytes", size);
	Q_memset(buf, 0, size);
	return buf;
}

void *Z_Realloc(void *ptr, s32 size)
{
	if(!ptr) return Z_Malloc(size);
	memblock_t *block = (memblock_t *) ((u8 *) ptr - sizeof(memblock_t));
	if(block->id != ZONEID)
		Sys_Error("Z_Realloc: realloced a pointer without ZONEID");
	if(block->tag == 0)
		Sys_Error("Z_Realloc: realloced a freed pointer");
	s32 old_size = block->size;
	old_size -= (4 + (s32)sizeof(memblock_t)); // see Z_TagMalloc()
	void *old_ptr = ptr;
	Z_Free(ptr);
	ptr = Z_TagMalloc(size, 1);
	if(!ptr) Sys_Error("Z_Realloc: failed on allocation of %i bytes", size);
	if(ptr != old_ptr) memmove(ptr, old_ptr, q_min(old_size, size));
	if(old_size < size) memset((u8 *)ptr + old_size, 0, size - old_size);
	return ptr;
}

s8 *Z_Strdup(const s8 *s)
{
	size_t sz = strlen(s) + 1;
	s8 *ptr = (s8 *) Z_Malloc(sz);
	memcpy(ptr, s, sz);
	return ptr;
}

void Z_Print(memzone_t *zone)
{
	Con_Printf("zone size: %i location: %p\n",mainzone->size,mainzone);
	for(memblock_t *block = zone->blocklist.next;; block = block->next) {
		Con_Printf("block:%p size:%7i tag:%3i\n",
				block, block->size, block->tag);
		if(block->next == &zone->blocklist)
			break; // all blocks have been hit
		if( (u8 *)block + block->size != (u8 *)block->next)
		Con_Printf("ERROR: block size does not touch the next block\n");
		if( block->next->prev != block)
		Con_Printf("ERROR: next block doesn't have proper back link\n");
		if(!block->tag && !block->next->tag)
			Con_Printf("ERROR: two consecutive free blocks\n");
	}
}

void Hunk_Check()
{ // Run consistancy and sentinel trahing checks
	for(hunk_t *h=(hunk_t *)hunk_base; (u8 *)h != hunk_base+hunk_low_used;){
		if(h->sentinel != HUNK_SENTINEL)
			Sys_Error("Hunk_Check: trashed sentinel");
		if(h->size < (s32) sizeof(hunk_t) ||
				h->size + (u8 *)h - hunk_base > hunk_size)
			Sys_Error("Hunk_Check: bad size");
		h = (hunk_t *)((u8 *)h+h->size);
	}
}

// If "all" is specified, every single allocation is printed.
// Otherwise, allocations with the same name will be totaled up before printing.
void Hunk_Print(bool all)
{
	s8 name[HUNKNAME_LEN];
	s32 count = 0;
	s32 sum = 0;
	s32 totalblocks = 0;
	hunk_t *h = (hunk_t *)hunk_base;
	hunk_t *endlow = (hunk_t *)(hunk_base + hunk_low_used);
	hunk_t *starthigh = (hunk_t *)(hunk_base + hunk_size - hunk_high_used);
	hunk_t *endhigh = (hunk_t *)(hunk_base + hunk_size);
	Con_Printf(" :%8i total hunk size\n", hunk_size);
	Con_Printf("-------------------------\n");
	while(1) {
		if(h == endlow){ // skip to the high hunk if done with low hunk
			Con_Printf("-------------------------\n");
	Con_Printf(" :%8i REMAINING\n",hunk_size-hunk_low_used-hunk_high_used);
			Con_Printf("-------------------------\n");
			h = starthigh;
		}
		if(h == endhigh) break; // if totally done, break
		if(h->sentinel != HUNK_SENTINEL) // run consistancy checks
			Sys_Error("Hunk_Check: trashed sentinel");
		if(h->size < (s32) sizeof(hunk_t) || h->size +
				(u8 *)h - hunk_base > hunk_size)
			Sys_Error("Hunk_Check: bad size");
		hunk_t *next = (hunk_t *)((u8 *)h+h->size);
		count++;
		totalblocks++;
		sum += h->size;
		memcpy(name, h->name, HUNKNAME_LEN); // print the single block
		if(all) Con_Printf("%8p :%8i %8s\n",h, h->size, name);
		if(next == endlow || next == endhigh ||
				strncmp(h->name, next->name, HUNKNAME_LEN - 1)){
			if(!all) Con_Printf(" :%8i %8s(TOTAL)\n",sum, name);
			count = 0;
			sum = 0;
		}
		h = next;
	}
	Con_Printf("-------------------------\n");
	Con_Printf("%8i total blocks\n", totalblocks);

}

void Hunk_Print_f() { Hunk_Print(0); }

void *Hunk_AllocName(s32 size, const s8 *name)
{
	if(size < 0) Sys_Error("Hunk_Alloc: bad size: %i", size);
	size = sizeof(hunk_t) + ((size+15)&~15);
	if(hunk_size - hunk_low_used - hunk_high_used < size)
		Sys_Error("Hunk_Alloc: failed on %i bytes",size);
	hunk_t *h = (hunk_t *)(hunk_base + hunk_low_used);
	hunk_low_used += size;
	Cache_FreeLow(hunk_low_used);
	memset(h, 0, size);
	h->size = size;
	h->sentinel = HUNK_SENTINEL;
	q_strlcpy(h->name, name, HUNKNAME_LEN);
	return(void *)(h+1);
}

void *Hunk_Alloc(s32 size)
{ return Hunk_AllocName(size, "unknown"); }

s32 Hunk_LowMark()
{ return hunk_low_used; }

void Hunk_FreeToLowMark(s32 mark)
{
	if(mark < 0 || mark > hunk_low_used)
		Sys_Error("Hunk_FreeToLowMark: bad mark %i", mark);
	memset(hunk_base + mark, 0, hunk_low_used - mark);
	hunk_low_used = mark;
}

s32 Hunk_HighMark()
{
	if(hunk_tempactive) {
		hunk_tempactive = 0;
		Hunk_FreeToHighMark(hunk_tempmark);
	}
	return hunk_high_used;
}

void Hunk_FreeToHighMark(s32 mark)
{
	if(hunk_tempactive) {
		hunk_tempactive = 0;
		Hunk_FreeToHighMark(hunk_tempmark);
	}
	if(mark < 0 || mark > hunk_high_used)
		Sys_Error("Hunk_FreeToHighMark: bad mark %i", mark);
	memset(hunk_base + hunk_size - hunk_high_used, 0, hunk_high_used-mark);
	hunk_high_used = mark;
}

void *Hunk_HighAllocName(s32 size, const s8 *name)
{
	if(size < 0) Sys_Error("Hunk_HighAllocName: bad size: %i", size);
	if(hunk_tempactive) {
		Hunk_FreeToHighMark(hunk_tempmark);
		hunk_tempactive = 0;
	}
	size = sizeof(hunk_t) + ((size+15)&~15);
	if(hunk_size - hunk_low_used - hunk_high_used < size) {
		Con_Printf("Hunk_HighAlloc: failed on %i bytes\n",size);
		return NULL;
	}
	hunk_high_used += size;
	Cache_FreeHigh(hunk_high_used);
	hunk_t *h = (hunk_t *)(hunk_base + hunk_size - hunk_high_used);
	memset(h, 0, size);
	h->size = size;
	h->sentinel = HUNK_SENTINEL;
	q_strlcpy(h->name, name, HUNKNAME_LEN);
	return(void *)(h+1);
}


void *Hunk_TempAlloc(s32 size) // Return space from the top of the hunk
{
	size = (size+15)&~15;
	if(hunk_tempactive) {
		Hunk_FreeToHighMark(hunk_tempmark);
		hunk_tempactive = 0;
	}
	hunk_tempmark = Hunk_HighMark();
	void *buf = Hunk_HighAllocName(size, "temp");
	hunk_tempactive = 1;
	return buf;
}

s8 *Hunk_Strdup(const s8 *s, const s8 *name)
{
	size_t sz = strlen(s) + 1;
	s8 *ptr = (s8 *) Hunk_AllocName(sz, name);
	memcpy(ptr, s, sz);
	return ptr;
}

void Cache_Move( cache_system_t *c)
{
	// we are clearing up space at the bottom, so only allocate it late
	cache_system_t *new_cs = Cache_TryAlloc(c->size, 1);
	if(new_cs) {
		Q_memcpy( new_cs+1, c+1, c->size - sizeof(cache_system_t) );
		new_cs->user = c->user;
		Q_memcpy(new_cs->name, c->name, sizeof(new_cs->name));
		Cache_Free(c->user, 0); //johnfitz -- added second argument
		new_cs->user->data = (void *)(new_cs+1);
	} else Cache_Free(c->user, 1); // tough luck...
}


void Cache_FreeLow(s32 new_low_hunk)
{ // Throw things out until the hunk can be expanded to the given point
	while(1) {
		cache_system_t *c = cache_head.next;
		if(c == &cache_head) return; // nothing in cache at all
		if((u8 *)c >= hunk_base + new_low_hunk)
			return; // there is space to grow the hunk
		Cache_Move( c ); // reclaim the space
	}
}


void Cache_FreeHigh(s32 new_high_hunk)
{ // Throw things out until the hunk can be expanded to the given point
	cache_system_t *prev = NULL;
	while(1) {
		cache_system_t *c = cache_head.prev;
		if(c == &cache_head)
			return; // nothing in cache at all
		if((u8 *)c + c->size <= hunk_base + hunk_size - new_high_hunk)
			return; // there is space to grow the hunk
		if(c == prev)
			Cache_Free(c->user, 1); // didn't move out of the way
		else {
			Cache_Move(c); // try to move it
			prev = c;
		}
	}
}

void Cache_UnlinkLRU(cache_system_t *cs)
{
	if(!cs->lru_next || !cs->lru_prev)
		Sys_Error("Cache_UnlinkLRU: NULL link");
	cs->lru_next->lru_prev = cs->lru_prev;
	cs->lru_prev->lru_next = cs->lru_next;
	cs->lru_prev = cs->lru_next = NULL;
}

void Cache_MakeLRU(cache_system_t *cs)
{
	if(cs->lru_next || cs->lru_prev)Sys_Error("Cache_MakeLRU: active link");
	cache_head.lru_next->lru_prev = cs;
	cs->lru_next = cache_head.lru_next;
	cs->lru_prev = &cache_head;
	cache_head.lru_next = cs;
}

// Looks for a free block of memory between the high and low hunk marks
// Size should already include the header and padding
cache_system_t *Cache_TryAlloc(s32 size, bool nobottom)
{
	cache_system_t *cs, *new_cs;
	// is the cache completely empty?
	if(!nobottom && cache_head.prev == &cache_head) {
		if(hunk_size - hunk_high_used - hunk_low_used < size)
		Sys_Error("Cache_TryAlloc: %i is greater then free hunk", size);
		new_cs = (cache_system_t *) (hunk_base + hunk_low_used);
		memset(new_cs, 0, sizeof(*new_cs));
		new_cs->size = size;
		cache_head.prev = cache_head.next = new_cs;
		new_cs->prev = new_cs->next = &cache_head;
		Cache_MakeLRU(new_cs);
		return new_cs;
	}
	// search from the bottom up for space
	new_cs = (cache_system_t *) (hunk_base + hunk_low_used);
	cs = cache_head.next;
	do {
		if(!nobottom || cs != cache_head.next) {
			if( (u8 *)cs - (u8 *)new_cs >= size) { // found space
				memset(new_cs, 0, sizeof(*new_cs));
				new_cs->size = size;
				new_cs->next = cs;
				new_cs->prev = cs->prev;
				cs->prev->next = new_cs;
				cs->prev = new_cs;
				Cache_MakeLRU(new_cs);
				return new_cs;
			}
		}
		// continue looking
		new_cs = (cache_system_t *)((u8 *)cs + cs->size);
		cs = cs->next;
	} while(cs != &cache_head);
	// try to allocate one at the very end
	if(hunk_base + hunk_size - hunk_high_used - (u8 *)new_cs >= size) {
		memset(new_cs, 0, sizeof(*new_cs));
		new_cs->size = size;
		new_cs->next = &cache_head;
		new_cs->prev = cache_head.prev;
		cache_head.prev->next = new_cs;
		cache_head.prev = new_cs;
		Cache_MakeLRU(new_cs);
		return new_cs;
	}
	return NULL; // couldn't allocate
}

void Cache_Flush()
{ // Throw everything out, so new data will be demand cached
	while(cache_head.next != &cache_head)
		Cache_Free( cache_head.next->user, 1); // reclaim the space
}

void Cache_Print()
{
	for(cache_system_t *cd=cache_head.next; cd!=&cache_head; cd=cd->next)
		Con_Printf("%8i : %s\n", cd->size, cd->name);
}

void Cache_Report()
{
	Con_DPrintf("%4.1f megabyte data cache\n",
		(hunk_size-hunk_high_used-hunk_low_used) / (f32)(1024*1024));
}

void Cache_Init()
{
	cache_head.next = cache_head.prev = &cache_head;
	cache_head.lru_next = cache_head.lru_prev = &cache_head;
	Cmd_AddCommand("flush", Cache_Flush);
}


void Cache_Free(cache_user_t *c, SDL_UNUSED bool freetextures)
{ // Frees the memory and removes it from the LRU list
	if(!c->data) Sys_Error("Cache_Free: not allocated");
	cache_system_t *cs = ((cache_system_t *)c->data) - 1;
	cs->prev->next = cs->next;
	cs->next->prev = cs->prev;
	cs->next = cs->prev = NULL;
	c->data = NULL;
	Cache_UnlinkLRU(cs);
}

void *Cache_Check(cache_user_t *c)
{
	if(!c->data) return NULL;
	cache_system_t *cs = ((cache_system_t *)c->data) - 1;
	Cache_UnlinkLRU(cs); // move to head of LRU
	Cache_MakeLRU(cs);
	return c->data;
}

void *Cache_Alloc(cache_user_t *c, s32 size, const s8 *name)
{
	if(c->data) Sys_Error("Cache_Alloc: already allocated");
	if(size <= 0) Sys_Error("Cache_Alloc: size %i", size);
	size = (size + sizeof(cache_system_t) + 15) & ~15;
	while(1) { // find memory for it
		cache_system_t *cs = Cache_TryAlloc(size, 0);
		if(cs) {
			q_strlcpy(cs->name, name, CACHENAME_LEN);
			c->data = (void *)(cs+1);
			cs->user = c;
			break;
		}
		// free the least recently used cahedat
		if(cache_head.lru_prev == &cache_head)
			Sys_Error("Cache_Alloc: out of memory");
		Cache_Free(cache_head.lru_prev->user, 1);
	}
	return Cache_Check(c);
}

static void Memory_InitZone(memzone_t *zone, s32 size)
{
	memblock_t *block;
	// set the entire zone to one free block
	zone->blocklist.next = zone->blocklist.prev = block =
		(memblock_t *)((u8 *)zone + sizeof(memzone_t));
	zone->blocklist.tag = 1; // in use block
	zone->blocklist.id = 0;
	zone->blocklist.size = 0;
	zone->rover = block;
	block->prev = block->next = &zone->blocklist;
	block->tag = 0; // free block
	block->id = ZONEID;
	block->size = size - sizeof(memzone_t);
}

void Memory_Init(void *buf, s32 size)
{
	s32 zonesize = DYNAMIC_SIZE;
	hunk_base = (u8 *) buf;
	hunk_size = size;
	hunk_low_used = 0;
	hunk_high_used = 0;
	Cache_Init();
	s32 p = COM_CheckParm("-zone");
	if(p) {
		if(p < com_argc-1) zonesize = Q_atoi(com_argv[p+1]) * 1024;
		else
	    Sys_Error("Memory_Init: you must specify a size in KB after -zone");
	}
	mainzone = (memzone_t *) Hunk_AllocName(zonesize, "zone" );
	Memory_InitZone(mainzone, zonesize);
	Cmd_AddCommand("hunk_print", Hunk_Print_f); //johnfitz
}
