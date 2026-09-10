// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2010-2011 O. Sezer <sezero@users.sourceforge.net>
// GPLv3 See LICENSE for details.
// snd_mem.c: sound caching
#include "quakedef.h"

static u8 *data_p;
static u8 *iff_end;
static u8 *last_chunk;
static u8 *iff_data;
static s32 iff_chunk_len;
wavinfo_t GetWavinfo(const s8 *name, u8 *wav, s32 wavlength);

static void ResampleSfx(sfx_t *sfx, s32 inrate, s32 inwidth, u8 *data)
{
	sfxcache_t *sc = (sfxcache_t *) Cache_Check(&sfx->cache);
	if(!sc) return;
	f32 stepscale = (f32)inrate / shm->speed; //this is usually 0.5, 1, or 2
	s32 outcount = sc->length / stepscale;
	sc->length = outcount;
	if(sc->loopstart != -1) sc->loopstart = sc->loopstart / stepscale;
	sc->speed = shm->speed;
	if(loadas8bit.value) sc->width = 1;
	else sc->width = inwidth;
	sc->stereo = 0;
	// resample / decimate to the current source rate
	if(stepscale == 1 && inwidth == 1 && sc->width == 1) {
		// fast special case
		for(s32 i = 0; i < outcount; i++)
			((s8 *)sc->data)[i] = (s32)((u8)(data[i]) - 128);
	} else {
		// general case
		// samplefrac can overflow 2**31 with very big sounds, see below
		int64_t samplefrac = 0;
		s32 fracstep = (s32)(stepscale * 256);
		for(s32 i = 0; i < outcount; i++) {
			s32 srcsample = (s32)(samplefrac >> 8);
			samplefrac += fracstep; // int64_t to prevent overflow
			s32 sample = inwidth == 2 ? 
				LittleShort(((s16 *)data)[srcsample]):
				(s32)((u8)(data[srcsample]) - 128) << 8;
			if(sc->width == 2) ((s16 *)sc->data)[i] = sample;
			else ((s8 *)sc->data)[i] = sample >> 8;
		}
	}
}

sfxcache_t *S_LoadSound(sfx_t *s)
{
	s8 namebuffer[256];
	u8 stackbuf[1*1024]; // avoid dirtying the cache heap
	sfxcache_t *sc = (sfxcache_t *) Cache_Check(&s->cache);
	if(sc) return sc; // see if still in memory
	q_strlcpy(namebuffer, "sound/", sizeof(namebuffer)); // load it in
	q_strlcat(namebuffer, s->name, sizeof(namebuffer));
	u8 *data = COM_LoadStackFile(namebuffer, stackbuf, sizeof(stackbuf), 0);
	if(!data) {
		Con_Printf("Couldn't load %s\n", namebuffer);
		return NULL;
	}
	wavinfo_t info = GetWavinfo(s->name, data, com_filesize);
	if(info.channels != 1) {
		Con_Printf("%s is a stereo sample\n",s->name);
		return NULL;
	}
	if(info.width != 1 && info.width != 2) {
		Con_Printf("%s is not 8 or 16 bit\n", s->name);
		return NULL;
	}
	f32 stepscale = (f32)info.rate / shm->speed;
	s32 len = info.samples / stepscale;
	len = len * info.width * info.channels;
	if(info.samples == 0 || len == 0) {
		Con_Printf("%s has zero samples\n", s->name);
		return NULL;
	}
	sc = (sfxcache_t*)Cache_Alloc(&s->cache,len+sizeof(sfxcache_t),s->name);
	if(!sc) return NULL;
	sc->length = info.samples;
	sc->loopstart = info.loopstart;
	sc->speed = info.rate;
	sc->width = info.width;
	sc->stereo = info.channels;
	ResampleSfx(s, sc->speed, sc->width, data + info.dataofs);
	return sc;
}

static s16 GetLittleShort()
{
	s16 val = *data_p;
	val = val + (*(data_p+1)<<8);
	data_p += 2;
	return val;
}

static s32 GetLittleLong()
{
	s32 val = *data_p;
	val = val + (*(data_p+1)<<8);
	val = val + (*(data_p+2)<<16);
	val = val + (*(data_p+3)<<24);
	data_p += 4;
	return val;
}

static void FindNextChunk(const s8 *name)
{
	while(1) {
		if(last_chunk + 8 >= iff_end) {
			data_p = NULL; // Need at least 8 bytes for a chunk
			return;
		}
		data_p = last_chunk + 4;
		iff_chunk_len = GetLittleLong();
		if(iff_chunk_len < 0 || iff_chunk_len > iff_end - data_p) {
			data_p = NULL;
			return;
		}
		last_chunk = data_p + ((iff_chunk_len + 1) & ~1);
		data_p -= 8;
		if(!strncmp((s8 *)data_p, name, 4)) return;
	}
}

static void FindChunk(const s8 *name)
{
	last_chunk = iff_data;
	FindNextChunk(name);
}

wavinfo_t GetWavinfo(const s8 *name, u8 *wav, s32 wavlength)
{
	wavinfo_t info;
	memset(&info, 0, sizeof(info));
	if(!wav) return info;
	iff_data = wav;
	iff_end = wav + wavlength;
	FindChunk("RIFF");
	if(!(data_p && !strncmp((s8 *)data_p + 8, "WAVE", 4))) {
		Con_Printf("%s missing RIFF/WAVE chunks\n", name);
		return info;
	}
	iff_data = data_p + 12;
	FindChunk("fmt ");
	if(!data_p) {
		Con_Printf("%s is missing fmt chunk\n", name);
		return info;
	}
	data_p += 8;
	s32 format = GetLittleShort();
	if(format != WAV_FORMAT_PCM) {
		Con_Printf("%s is not Microsoft PCM format\n", name);
		return info;
	}
	info.channels = GetLittleShort();
	info.rate = GetLittleLong();
	data_p += 4 + 2;
	s32 i = GetLittleShort();
	if(i != 8 && i != 16) return info;
	info.width = i / 8;
	FindChunk("cue ");
	if(data_p) {
		data_p += 32;
		info.loopstart = GetLittleLong();
		//if the next chunk is a LIST chunk look for a cue length marker
		FindNextChunk("LIST");
		if(data_p) {
			if(!strncmp((s8 *)data_p + 28, "mark", 4)) {
				data_p += 24;
				i = GetLittleLong(); // samples in loop
				info.samples = info.loopstart + i;
			}
		}
	} else info.loopstart = -1;
	FindChunk("data");
	if(!data_p) {
		Con_Printf("%s is missing data chunk\n", name);
		return info;
	}
	data_p += 4;
	s32 samples = GetLittleLong() / info.width;
	if(info.samples) {
		if(samples < info.samples)
			Sys_Error("%s has a bad loop length", name);
	} else info.samples = samples;
	if(info.loopstart >= info.samples) {
		printf("%s has loop start >= end\n", name);
		info.loopstart = -1;
		info.samples = samples;
	}
	info.dataofs = data_p - wav;
	return info;
}
