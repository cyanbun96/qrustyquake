// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.

#include "quakedef.h"

int wad_numlumps;
lumpinfo_t *wad_lumps;
byte *wad_base;

void SwapPic(qpic_t *pic)
{
	pic->width = LittleLong(pic->width);
	pic->height = LittleLong(pic->height);
}

// Lowercases name and pads with spaces and a terminating 0 to the length of
// lumpinfo_t->name.
// Used so lumpname lookups can proceed rapidly by comparing 4 chars at a time
// Space padding is so names can be printed nicely in tables.
// Can safely be performed in place.
void W_CleanupName(char *in, char *out)
{
	int i = 0;
	for (; i < 16; i++) {
		char c = in[i];
		if (!c)
			break;
		if (c >= 'A' && c <= 'Z')
			c += ('a' - 'A');
		out[i] = c;
	}
	for (; i < 16; i++)
		out[i] = 0;
}

void W_LoadWadFile()
{
	wad_base = COM_LoadHunkFile("gfx.wad");
	if (!wad_base)
		Sys_Error("W_LoadWadFile: couldn't load %s", "gfx.wad");
	wadinfo_t *header = (wadinfo_t *) wad_base;
	if (header->identification[0] != 'W'
	 || header->identification[1] != 'A'
	 || header->identification[2] != 'D'
	 || header->identification[3] != '2')
		Sys_Printf("Wad file %s doesn't have WAD2 id\n", "gfx.wad");
	wad_numlumps = LittleLong(header->numlumps);
	int infotableofs = LittleLong(header->infotableofs);
	wad_lumps = (lumpinfo_t *) (wad_base + infotableofs);
	lumpinfo_t *lump_p = wad_lumps;
	for (int i = 0; i < wad_numlumps; i++, lump_p++) {
		lump_p->filepos = LittleLong(lump_p->filepos);
		lump_p->size = LittleLong(lump_p->size);
		W_CleanupName(lump_p->name, lump_p->name);
		if (lump_p->type == TYP_QPIC)
			SwapPic((qpic_t *) (wad_base + lump_p->filepos));
	}
}

lumpinfo_t *W_GetLumpinfo(char *name)
{
	char clean[16];
	W_CleanupName(name, clean);
	lumpinfo_t *lump_p = wad_lumps;
	for (int i = 0; i < wad_numlumps; i++, lump_p++)
		if (!strcmp(clean, lump_p->name))
			return lump_p;
	Sys_Error("W_GetLumpinfo: %s not found", name);
	return NULL;
}

void *W_GetLumpName(char *name)
{
	lumpinfo_t *lump = W_GetLumpinfo(name);
	return (void *)(wad_base + lump->filepos);
}
