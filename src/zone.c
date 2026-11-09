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

typedef struct cache_system_s {
	s32 size; // including this header
	cache_user_t *user;
	c8 name[CACHENAME_LEN];
	struct cache_system_s *prev, *next;
	struct cache_system_s *lru_prev, *lru_next; // for LRU flushing
} cache_system_t;
typedef struct memblock_s {
	s32 size; // including the header and possibly tiny fragments
	s32 tag; // a tag of 0 is a free block
	s32 id; // should be ZONEID
	s32 pad; // pad to 64 bit boundary
	struct memblock_s *next, *prev;
} memblock_t;
typedef struct {
	s32 size; // total bytes malloced, including header
	memblock_t blocklist; // start / end cap for linked list
	memblock_t *rover;
} memzone_t;
typedef struct { // the hunk grows by these, each followed by its memory
	s32 base; // hunk offset of the first byte
	s32 size;
	s32 used; // by low hunk blocks
	s32 pad; // keep the memory 16 byte aligned
} hunkseg_t;
#define MAX_SEGMENTS 8
#define SEG_MEM(seg) ((u8 *)((seg) + 1))
#define SEG_END(seg) (SEG_MEM(seg) + (seg)->size)
#define LASTSEG (hunk_segs[hunk_numsegs-1])
static memzone_t *mainzone;
static hunkseg_t *hunk_segs[MAX_SEGMENTS];
static s32 hunk_numsegs;
static s32 hunk_low_used; // offsets run through the segments in order
static cache_system_t cache_head;

void Cache_FreeLow(s32 new_low_hunk);
cache_system_t *Cache_TryAlloc(s32 size, bool nobottom);
void Cache_Free(cache_user_t *c);

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


static void *Z_TagMalloc(s32 size)
{
	memblock_t *start, *rover, *newblock, *base;
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
	if(extra > 64/*minfragment*/)
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
	base->tag = 1; // no longer a free block
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
	Z_CheckHeap();
	void *buf = Z_TagMalloc(size);
	if(!buf) Sys_Error("Z_Malloc: failed on allocation of %i bytes", size);
	memset(buf, 0, size);
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
	ptr = Z_TagMalloc(size);
	if(!ptr) Sys_Error("Z_Realloc: failed on allocation of %i bytes", size);
	if(ptr != old_ptr) memmove(ptr, old_ptr, q_min(old_size, size));
	if(old_size < size) memset((u8 *)ptr + old_size, 0, size - old_size);
	return ptr;
}

c8 *Z_Strdup(const c8 *s)
{
	size_t sz = strlen(s) + 1;
	c8 *ptr = (c8 *) Z_Malloc(sz);
	memcpy(ptr, s, sz);
	return ptr;
}

static s32 Hunk_Size()
{ return LASTSEG->base + LASTSEG->size; }

static s32 Hunk_SegForOfs(s32 ofs)
{ // an offset at the end of a segment is the start of the next one
	s32 i = 0;
	while(i < hunk_numsegs - 1 && ofs >= hunk_segs[i+1]->base) i++;
	return i;
}

static hunk_t *Hunk_Block(hunkseg_t *seg, s32 ofs)
{ // Run consistancy and sentinel trahing checks
	hunk_t *h = (hunk_t *)(SEG_MEM(seg) + ofs);
	if(h->sentinel != HUNK_SENTINEL)
		Sys_Error("Hunk_Check: trashed sentinel");
	if(h->size < (s32) sizeof(hunk_t) || h->size > seg->used - ofs)
		Sys_Error("Hunk_Check: bad size");
	return h;
}

void Hunk_Check()
{
	for(s32 i = 0; i < hunk_numsegs; i++)
		for(s32 ofs = 0; ofs < hunk_segs[i]->used;)
			ofs += Hunk_Block(hunk_segs[i], ofs)->size;
}

bool Hunk_Contains(const void *p)
{
	for(s32 i = 0; i < hunk_numsegs; i++)
		if((uintptr_t)p - (uintptr_t)SEG_MEM(hunk_segs[i])
				< (uintptr_t)hunk_segs[i]->size)
			return 1;
	return 0;
}

// Allocations with the same name are totaled up before printing.
void Hunk_Print_f()
{
	const c8 *name = NULL;
	s32 sum = 0;
	s32 totalblocks = 0;
	Con_Printf(" :%8i total hunk size\n", Hunk_Size());
	for(s32 i = 0; i < hunk_numsegs && hunk_numsegs > 1; i++)
		Con_Printf(" :%8i segment %i\n", hunk_segs[i]->size, i);
	Con_Printf("-------------------------\n");
	for(s32 i = 0; i < hunk_numsegs; i++)
		for(s32 ofs = 0; ofs < hunk_segs[i]->used;){
			hunk_t *h = Hunk_Block(hunk_segs[i], ofs);
			if(name && strncmp(name, h->name, HUNKNAME_LEN - 1)){
				Con_Printf(" :%8i %8s(TOTAL)\n", sum, name);
				sum = 0;
			}
			name = h->name;
			sum += h->size;
			totalblocks++;
			ofs += h->size;
		}
	if(name) Con_Printf(" :%8i %8s(TOTAL)\n", sum, name);
	Con_Printf("-------------------------\n");
	Con_Printf(" :%8i REMAINING\n", Hunk_Size() - hunk_low_used);
	Con_Printf("-------------------------\n");
	Con_Printf("%8i total blocks\n", totalblocks);
}

static void Hunk_AddSegment(s32 size)
{ // Grow instead of failing, like Ironwail. A segment is twice the last one,
  // or less if that much memory isn't there.
	s32 base = Hunk_Size();
	s64 left = (s64)INT32_MAX - base; // offsets are s32
	if(hunk_numsegs == MAX_SEGMENTS || size > left)
		Sys_Error("Hunk_Alloc: failed on %i bytes", size);
	s32 newsize = q_max(q_min(2 * (s64)LASTSEG->size, left), size);
	hunkseg_t *seg;
	while(!(seg = malloc(sizeof(hunkseg_t) + newsize)) && newsize > size)
		newsize = q_max(newsize / 4 * 3, size);
	if(!seg) Sys_Error("Hunk_Alloc: failed on %i bytes", size);
	Cache_Flush(); // the cache only lives in the last segment
	seg->base = base;
	seg->size = newsize;
	seg->used = 0;
	hunk_segs[hunk_numsegs++] = seg;
	Con_DPrintf("Hunk grown to %i MB\n", Hunk_Size() >> 20);
}

void *Hunk_AllocInternal(s32 size, const c8 *name, bool clear)
{
	if(size < 0) Sys_Error("Hunk_Alloc: bad size: %i", size);
	size = sizeof(hunk_t) + ((size+15)&~15);
	s32 i = Hunk_SegForOfs(hunk_low_used);
	for(; i < hunk_numsegs && size > hunk_segs[i]->base
			+ hunk_segs[i]->size - hunk_low_used; i++)
		hunk_low_used = hunk_segs[i]->base + hunk_segs[i]->size;
	if(i == hunk_numsegs) Hunk_AddSegment(size);
	hunkseg_t *seg = hunk_segs[i];
	hunk_t *h = (hunk_t *)(SEG_MEM(seg) + hunk_low_used - seg->base);
	hunk_low_used += size;
	seg->used = hunk_low_used - seg->base;
	Cache_FreeLow(hunk_low_used);
	if(clear)memset(h, 0, size);
	h->size = size;
	h->sentinel = HUNK_SENTINEL;
	SDL_strlcpy(h->name, name, HUNKNAME_LEN);
	return(void *)(h+1);
}

void *Hunk_AllocNoFill(s32 size)
{ return Hunk_AllocInternal(size, "nofill", 0); }

void *Hunk_AllocName(s32 size, const c8 *name)
{ return Hunk_AllocInternal(size, name, 1); }

void *Hunk_Alloc(s32 size)
{ return Hunk_AllocName(size, "unknown"); }

s32 Hunk_LowMark()
{ return hunk_low_used; }

void Hunk_FreeToLowMark(s32 mark)
{
	if(mark < 0 || mark > hunk_low_used)
		Sys_Error("Hunk_FreeToLowMark: bad mark %i", mark);
	for(s32 i = Hunk_SegForOfs(mark); i < hunk_numsegs; i++){
		hunkseg_t *seg = hunk_segs[i];
		s32 used = q_max(mark - seg->base, 0);
		if(seg->used > used)
			memset(SEG_MEM(seg) + used, 0, seg->used - used);
		seg->used = used;
	}
	hunk_low_used = mark;
}

bool Hunk_Restart(s32 mark)
{ // If the blocks since mark went into two segments, frees them and moves to
  // the start of the next segment, so they can be built again in one piece.
	s32 i = Hunk_SegForOfs(mark);
	if(hunk_low_used <= hunk_segs[i]->base + hunk_segs[i]->size)
		return 0;
	Hunk_FreeToLowMark(mark);
	hunk_low_used = hunk_segs[i+1]->base;
	return 1;
}

void Cache_Move( cache_system_t *c)
{
	R_ThreadsSync(); // the renderer's threads may be drawing from it
	// we are clearing up space at the bottom, so only allocate it late
	cache_system_t *new_cs = Cache_TryAlloc(c->size, 1);
	if(new_cs) {
		memcpy( new_cs+1, c+1, c->size - sizeof(cache_system_t) );
		new_cs->user = c->user;
		memcpy(new_cs->name, c->name, sizeof(new_cs->name));
		Cache_Free(c->user);
		new_cs->user->data = (void *)(new_cs+1);
	} else Cache_Free(c->user); // tough luck...
}


static u8 *Cache_Bottom(s32 low_hunk)
{ // the cache lives in the last segment, above the low hunk
	return SEG_MEM(LASTSEG) + q_max(low_hunk - LASTSEG->base, 0);
}

void Cache_FreeLow(s32 new_low_hunk)
{ // Throw things out until the hunk can be expanded to the given point
	u8 *bottom = Cache_Bottom(new_low_hunk);
	while(1) {
		cache_system_t *c = cache_head.next;
		if(c == &cache_head) return; // nothing in cache at all
		if((u8 *)c >= bottom)
			return; // there is space to grow the hunk
		Cache_Move( c ); // reclaim the space
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

// Looks for a free block of memory between the low hunk and the end
// Size should already include the header and padding
cache_system_t *Cache_TryAlloc(s32 size, bool nobottom)
{
	cache_system_t *cs, *new_cs;
	// is the cache completely empty?
	if(!nobottom && cache_head.prev == &cache_head) {
		if(SEG_END(LASTSEG) - Cache_Bottom(hunk_low_used) < size)
			Hunk_AddSegment(size);
		new_cs = (cache_system_t *) Cache_Bottom(hunk_low_used);
		memset(new_cs, 0, sizeof(*new_cs));
		new_cs->size = size;
		cache_head.prev = cache_head.next = new_cs;
		new_cs->prev = new_cs->next = &cache_head;
		Cache_MakeLRU(new_cs);
		return new_cs;
	}
	// search from the bottom up for space
	new_cs = (cache_system_t *) Cache_Bottom(hunk_low_used);
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
	if(SEG_END(LASTSEG) - (u8 *)new_cs >= size) {
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
		Cache_Free(cache_head.next->user); // reclaim the space
}

void Cache_Report()
{
	Con_DPrintf("%4.1f megabyte data cache\n", (SEG_END(LASTSEG)
		- Cache_Bottom(hunk_low_used)) / (f32)(1024*1024));
}

void Cache_Init()
{
	cache_head.next = cache_head.prev = &cache_head;
	cache_head.lru_next = cache_head.lru_prev = &cache_head;
	Cmd_AddCommand("flush", Cache_Flush);
}


void Cache_Free(cache_user_t *c)
{ // Frees the memory and removes it from the LRU list
	if(!c->data) Sys_Error("Cache_Free: not allocated");
	R_ThreadsSync(); // the renderer's threads may be drawing from it
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

void *Cache_Alloc(cache_user_t *c, s32 size, const c8 *name)
{
	if(c->data) Sys_Error("Cache_Alloc: already allocated");
	if(size <= 0) Sys_Error("Cache_Alloc: size %i", size);
	size = (size + sizeof(cache_system_t) + 15) & ~15;
	while(1) { // find memory for it
		cache_system_t *cs = Cache_TryAlloc(size, 0);
		if(cs) {
			SDL_strlcpy(cs->name, name, CACHENAME_LEN);
			c->data = (void *)(cs+1);
			cs->user = c;
			break;
		}
		// free the least recently used cahedat
		if(cache_head.lru_prev == &cache_head)
			Sys_Error("Cache_Alloc: out of memory");
		Cache_Free(cache_head.lru_prev->user);
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
	s32 zonesize = 4*1024*1024;
	hunk_segs[0] = buf;
	hunk_segs[0]->base = 0;
	hunk_segs[0]->size = size - sizeof(hunkseg_t);
	hunk_segs[0]->used = 0;
	hunk_numsegs = 1;
	hunk_low_used = 0;
	Cache_Init();
	s32 p = COM_CheckParm("-zone");
	if(p) {
		if(p < com_argc-1) zonesize = atoi(com_argv[p+1]) * 1024;
		else
	    Sys_Error("Memory_Init: you must specify a size in KB after -zone");
	}
	// Z_Malloc clears each block, so the zone's pages stay untouched until
	// it's used
	mainzone = (memzone_t *) Hunk_AllocInternal(zonesize, "zone", 0);
	Memory_InitZone(mainzone, zonesize);
	Cmd_AddCommand("hunk_print", Hunk_Print_f); //johnfitz
}
