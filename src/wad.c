// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
#include "quakedef.h"

static s32 wad_numlumps;
static lumpinfo_t *wad_lumps;
static u8 *wad_base = NULL;

void SwapPic(qpic_t *pic);

// Lowercases name and pads with spaces and a terminating 0 to the length of
// lumpinfo_t->name.
// Used so lumpname lookups can proceed rapidly by comparing 4 chars at a time
// Space padding is so names can be printed nicely in tables.
// Can safely be performed in place.
void W_CleanupName(const s8 *in, s8 *out)
{
	s32 i = 0;
	for(; i < 16; i++) {
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
	const s8 *filename = WADFILENAME;
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
	wad_numlumps = LittleLong(header->numlumps);
	s32 infotableofs = LittleLong(header->infotableofs);
	wad_lumps = (lumpinfo_t *)(wad_base + infotableofs);
	lumpinfo_t *lump_p = wad_lumps;
	for(s32 i = 0; i < wad_numlumps; i++, lump_p++) {
		lump_p->filepos = LittleLong(lump_p->filepos);
		lump_p->size = LittleLong(lump_p->size);
		W_CleanupName(lump_p->name, lump_p->name);
		if(lump_p->type == TYP_QPIC)
			SwapPic((qpic_t *)(wad_base + lump_p->filepos));
	}
}

static lumpinfo_t *W_GetLumpinfo(lumpinfo_t *lumps, s32 numlumps,const s8 *name)
{
	s8 clean[16];
	W_CleanupName(name, clean);
	lumpinfo_t *lump_p = lumps;
	for(s32 i = 0; i < numlumps; i++, lump_p++) {
		if(!strcmp(clean, lump_p->name))
			return lump_p;
	}
	Con_SafePrintf("W_GetLumpinfo: %s not found\n", name);
	return NULL;
}

void *W_GetLumpName(const s8 *name)
{
	lumpinfo_t *lump = W_GetLumpinfo(wad_lumps, wad_numlumps, name);
	if(!lump) return NULL; //johnfitz
	return(void *)(wad_base + lump->filepos);
}

void *W_GetLumpNum(s32 num)
{
	if(num < 0 || num > wad_numlumps)
		Sys_Error("W_GetLumpNum: bad number: %i", num);
	lumpinfo_t *lump = wad_lumps + num;
	return(void *)(wad_base + lump->filepos);
}

static bool W_OpenWadFile(const s8 *filename, fshandle_t *fh)
{
	FILE *f;
	s64 length = (s64)COM_FOpenFile(filename, &f, NULL);
	if(length == -1) return 0;
	fh->file = f;
	fh->start = ftell(f);
	fh->pos = 0;
	fh->length = length;
	fh->pak = file_from_pak;
	return 1;
}

static wad_t *W_AddWadFile(const s8 *name, fshandle_t *fh)
{
	wadinfo_t header;
	FS_fread((void *)&header, 1, sizeof(header), fh);
	s32 id=LittleLong(header.identification[0]|(header.identification[1]<<8)
		|(header.identification[2]<<16)|(header.identification[3]<<24));
	if(id != WADID && id != WADID_VALVE) {
		printf("%s is not a valid WAD\n", name);
		return NULL;
	}
	s32 numlumps = LittleLong(header.numlumps);
	s32 infotableofs = LittleLong(header.infotableofs);
	if(numlumps < 0 || infotableofs < 0) {
		printf("%s is not a valid WAD(%i lumps, %i info table offset)\n"
				, name, numlumps, infotableofs);
		return NULL;
	}
	if(!numlumps) {
		printf("WAD file %s has no lumps, ignored\n", name);
		return NULL;
	}
	lumpinfo_t *lumps = (lumpinfo_t *)malloc(numlumps * sizeof(lumpinfo_t));
	FS_fseek(fh, infotableofs, SEEK_SET);
	FS_fread(lumps, 1, numlumps * sizeof(lumpinfo_t), fh);
	// parse the directory
	lumpinfo_t *info = lumps;
	for(s32 i = 0; i < numlumps; i++, info++) {
		W_CleanupName(info->name, info->name);
		info->filepos = LittleLong(info->filepos);
		info->size = LittleLong(info->size);
		s32 disksize = LittleLong(info->disksize);

		if(info->filepos + info->size > fh->length &&
				!(info->filepos + disksize > fh->length))
			info->size = disksize;
		// ensure lump sanity
		if(info->filepos < 0 || info->size < 0 || info->filepos
				+ info->size > fh->length) {
			if(info->filepos > fh->length || info->size < 0) {

printf("WAD file %s lump \"%.16s\" begins %li bytes beyond end of WAD\n"
		, name, info->name, info->filepos - fh->length);

				info->filepos = 0;
				info->size = q_max(0,info->size-info->filepos);
			} else { 

printf( "WAD file %s lump \"%.16s\" extends %li bytes beyond end of WAD(lump size is %i)\n"
, name, info->name, (info->filepos + info->size) - fh->length, info->size);

				info->size = q_max(0,info->size-info->filepos);
			}
		}
	}
	wad_t *wad = (wad_t *)malloc(sizeof(wad_t));
	q_strlcpy(wad->name, name, sizeof(wad->name));
	wad->id = id;
	wad->fh = *fh;
	wad->numlumps = numlumps;
	wad->lumps = lumps;
	Con_DPrintf("%s\n", name);
	return wad;
}

wad_t *W_LoadWadList(const s8 *names)
{
	s8 *newnames = q_strdup(names);
	wad_t *wad, *wads = NULL;
	s8 filename[MAX_QPATH];
	for(s8 *name = newnames; name && *name;) {
		s8 *e = strchr(name, ';');
		if(e) *e++ = 0;
		// remove all of the leading garbage left by the map editor
		COM_FileBase(name, filename, sizeof(filename));
		COM_AddExtension(filename, ".wad", sizeof(filename));
		fshandle_t fh;
		if(!W_OpenWadFile(filename, &fh)) { // try the "gfx" directory
			memmove(filename + 4, filename, sizeof(filename) - 4);
			memcpy(filename, "gfx/", 4);
			filename[sizeof(filename) - 1] = 0;
			if(!W_OpenWadFile(filename, &fh)) {
				name = e;
				continue;
			}
		}
		wad = W_AddWadFile(filename, &fh);
		if(wad) {
			wad->next = wads;
			wads = wad;
		} else FS_fclose(&fh);
		name = e;
	}
	free(newnames);
	return wads;
}

void W_FreeWadList(wad_t *wads)
{
	while(wads) {
		FS_fclose(&wads->fh);
		free(wads->lumps);
		wad_t *next = wads->next;
		free(wads);
		wads = next;
	}
}

lumpinfo_t *W_GetLumpinfoList(wad_t *wads, const s8 *name, wad_t **out_wad)
{
	while(wads) {
		lumpinfo_t *info=W_GetLumpinfo(wads->lumps,wads->numlumps,name);
		if(info) {
			*out_wad = wads;
			return info;
		}
		wads = wads->next;
	}
	return NULL;
}

void SwapPic(qpic_t *pic)
{
	pic->width = LittleLong(pic->width);
	pic->height = LittleLong(pic->height);
}
