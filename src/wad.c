// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
#include "quakedef.h"

typedef struct {
	c8 identification[4]; // should be WAD2 or 2DAW
	s32 numlumps;
	s32 infotableofs;
} wadinfo_t;
typedef struct {
	s32 filepos;
	s32 disksize;
	s32 size; // uncompressed
	c8 type;
	c8 compression;
	c8 pad1, pad2;
	c8 name[16]; // must be null terminated
} lumpinfo_t;
static s32 wad_numlumps;
static lumpinfo_t *wad_lumps;
static u8 *wad_base = NULL;

// Lowercases name and pads with spaces and a terminating 0 to the length of
// lumpinfo_t->name.
// Used so lumpname lookups can proceed rapidly by comparing 4 chars at a time
// Space padding is so names can be printed nicely in tables.
// Can safely be performed in place.
void W_CleanupName(const c8 *in, c8 *out)
{
	s32 i = 0;
	for(; i < 16; i++){
		s32 c = in[i];
		if(!c) break;
		if(c >= 'A' && c <= 'Z')
			c += ('a' - 'A');
		out[i] = c;
	}
	for(; i < 16; i++) out[i] = 0;
}

void W_LoadWadFile()
{ //johnfitz -- filename is now hard-coded for honesty
	const c8 *filename = "gfx.wad";
	if(wad_base) free(wad_base);
	wad_base = COM_LoadMallocFile(filename, NULL);
	if(!wad_base)
		Sys_Error("W_LoadWadFile: couldn't load %s\n\n"
"Basedir is: %s\n\n"
"Check that this has an " GAMENAME " subdirectory containing pak0.pak and pak1.pak, "
"or use the -basedir command-line option to specify another directory.",
				filename, com_basedir);
	wadinfo_t *header = (wadinfo_t *)wad_base;
	if(header->identification[0]!='W'||header->identification[1]!='A'||
		header->identification[2]!='D'||header->identification[3]!='2')
		Sys_Error("Wad file %s doesn't have WAD2 id\n",filename);
	wad_numlumps = (s32)(header->numlumps);
	s32 infotableofs = (s32)(header->infotableofs);
	wad_lumps = (lumpinfo_t *)(wad_base + infotableofs);
	lumpinfo_t *lump_p = wad_lumps;
	for(s32 i = 0; i < wad_numlumps; i++, lump_p++){
		lump_p->filepos = (s32)(lump_p->filepos);
		lump_p->size = (s32)(lump_p->size);
		W_CleanupName(lump_p->name, lump_p->name);
	}
}

static lumpinfo_t *W_GetLumpinfo(lumpinfo_t *lumps, s32 numlumps,const c8 *name)
{
	c8 clean[16];
	W_CleanupName(name, clean);
	lumpinfo_t *lump_p = lumps;
	for(s32 i = 0; i < numlumps; i++, lump_p++){
		if(!strcmp(clean, lump_p->name))
			return lump_p;
	}
	Con_SafePrintf("W_GetLumpinfo: %s not found\n", name);
	return NULL;
}

void *W_GetLumpName(const c8 *name)
{
	lumpinfo_t *lump = W_GetLumpinfo(wad_lumps, wad_numlumps, name);
	if(!lump) return NULL; //johnfitz
	return(void *)(wad_base + lump->filepos);
}
