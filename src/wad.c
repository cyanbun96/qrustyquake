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
#define TYP_QPIC 66 // a lump that is a qpic_t
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

static bool W_LumpOK(const lumpinfo_t *l, s64 size)
{ // inside the file, and so are the pixels its users read
	if(l->filepos < 0 || l->disksize < 0
			|| (s64)l->filepos + l->disksize > size)
		return 0;
	if(!strcmp(l->name, "conchars")) // raw, Draw_Character reads 128x128
		return l->disksize >= 128 * 128;
	if(l->type != TYP_QPIC)
		return 1;
	const qpic_t *pic = (qpic_t *)(wad_base + l->filepos);
	return l->disksize >= 8 && pic->width >= 0 && pic->height >= 0
		&& (s64)pic->width * pic->height + 8 <= l->disksize;
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
	if(com_filesize < (s32)sizeof(wadinfo_t)
			|| memcmp(header->identification, "WAD2", 4))
		Sys_Error("Wad file %s doesn't have WAD2 id\n",filename);
	wad_numlumps = header->numlumps;
	s32 infotableofs = header->infotableofs;
	s64 size = com_filesize;
	if(infotableofs < (s32)sizeof(wadinfo_t) || wad_numlumps < 0
			|| infotableofs + wad_numlumps
			* (s64)sizeof(lumpinfo_t) > size)
		Sys_Error("Wad file %s has a bad lump table\n", filename);
	wad_lumps = (lumpinfo_t *)(wad_base + infotableofs);
	lumpinfo_t *lump_p = wad_lumps;
	for(s32 i = 0; i < wad_numlumps; i++, lump_p++){
		lump_p->name[15] = 0;
		W_CleanupName(lump_p->name, lump_p->name);
		if(!W_LumpOK(lump_p, size))
			Sys_Error("Wad file %s has a bad lump %s\n",
					filename, lump_p->name);
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
