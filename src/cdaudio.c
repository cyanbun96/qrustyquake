// Copyright (C) 1996-2001 Id Software, Inc.
// Copyright (C) 2002-2009 John Fitzgibbons and others
// Copyright (C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
#include "quakedef.h"

#ifndef AVAIL_SDL3MIXER
static void CDAudio_NotImplemented_f(void)
{ Con_Printf("Not implemented in this build, sorry!\n"); }
void CDAudio_Play(SDL_UNUSED u8 track, SDL_UNUSED bool looping) { }
void CDAudio_Stop() { }
void CDAudio_Pause() { }
void CDAudio_Resume() { }
void CDAudio_Update() { }
void CDAudio_Shutdown() { }
bool CDAudio_Init()
{
	Con_Printf("Music unavailable\n");
	Cmd_AddCommand("music", CDAudio_NotImplemented_f);
	Cmd_AddCommand("music_stop", CDAudio_NotImplemented_f);
	Cmd_AddCommand("music_pause", CDAudio_NotImplemented_f);
	Cmd_AddCommand("music_resume", CDAudio_NotImplemented_f);
	return false;
}
#else

// erysdren: note that i've placed these in a particular order, so that when
// they're iterated in CDAudio_Play() the most common formats for Quake mod
// music will come up first
static struct music_format {
	const char *name;
	size_t num_extensions;
	const char **extensions;
	bool required;
} music_formats[] = {
	{"OGG", 1, (const char *[]){".ogg"}, true},
	{"OPUS", 2, (const char *[]){".ogg", ".opus"}, true},
	{"MP3", 1, (const char *[]){".mp3"}, true},
	{"FLAC", 1, (const char *[]){".flac"}, false},
	{"MID", 2, (const char *[]){".mid", ".midi"}, false},
	{"MOD", 1, (const char *[]){".mod"}, false},
	{"WAVPACK", 2, (const char *[]){".wav", ".wv"}, false}
};

static MIX_Mixer *mixer = NULL;
static MIX_Audio *current_music = NULL;
static s8 current_name[MAX_OSPATH];
static u8 *loaded_file = NULL;
static float last_volume = -1;

void BGM_Play(s8 *musicname, bool looping)
{
	s8 filename[MAX_OSPATH];
	u8 *file = NULL;
	SDL_IOStream *io = NULL;
	MIX_Audio *music = NULL;
	if(!musicname || !*musicname) {
		Con_DPrintf("null music file name\n");
		return;
	}
	CDAudio_Stop();
	for(s32 i = 0; i < Q_COUNTOF(music_formats); i++) {
		for(s32 j = 0; j < music_formats[i].num_extensions; j++) {
			q_snprintf(filename, sizeof(filename), "music/%s%s",
				musicname, music_formats[i].extensions[j]);
			file = COM_LoadMallocFile(filename, NULL);
			if(file) goto found;
		}
	}
	if(!file) {
		Con_Printf("file for %s not found\n", filename);
		return;
	}
found:
	loaded_file = file;
	io = SDL_IOFromConstMem(file, com_filesize);
	if(!io) {
		Con_Printf("failed to create IOStream for %s: %s\n", filename, SDL_GetError());
		return;
	}
	music = MIX_LoadAudio_IO(mixer, io, true, true);
	if (!music) {
		Con_Printf("failed to load %s: %s\n", filename, SDL_GetError());
		return;
	}
	CDAudio_Stop();
	current_music = music;
	if (!MIX_PlayAudio(mixer, current_music)) {
		Con_Printf("failed to play %s: %s\n", filename, SDL_GetError());
		return;
	}
	q_strlcpy(current_name, musicname, MAX_OSPATH);
}

static void BGM_Play_f()
{
	if(Cmd_Argc() == 2)
		BGM_Play(Cmd_Argv(1), 1);
	else {
		if(current_music) {
			Con_Printf("Playing %s, use 'music <musicfile>' to change\n", current_name);
		} else Con_Printf("music <musicfile>\n");
	}
}

void CDAudio_Play(u8 track, bool looping)
{
	s8 name[16];
	q_snprintf(name, sizeof(name), "track%02d", (s32)track);
	BGM_Play(name, looping);
}

void CDAudio_Stop()
{
	MIX_StopAllTracks(mixer, 0);
	if (loaded_file) free(loaded_file);
	loaded_file = NULL;
	memset(current_name, 0, sizeof(current_name));
	current_music = NULL;
}
void CDAudio_Update()
{
	if (bgmvolume.value < 0) Cvar_SetQuick(&bgmvolume, "0");
	if (bgmvolume.value > 1) Cvar_SetQuick(&bgmvolume, "1");
	if (last_volume != bgmvolume.value) {
		last_volume = bgmvolume.value;
		MIX_SetMasterGain(mixer, bgmvolume.value);
	}
}

bool CDAudio_Init()
{
	if (safemode || COM_CheckParm("-nosound")) return;
	MIX_Init();
	mixer = MIX_CreateMixerDevice(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, 0);
	if (mixer) Con_Printf("SDL Mixer initialized\n");
	else {
		Con_Printf("SDL Mixer initialization failed: %s\n", SDL_GetError());
		return false;
	}
	Cmd_AddCommand("music", BGM_Play_f);
	Cmd_AddCommand("music_stop", CDAudio_Stop);
	Cmd_AddCommand("music_pause", CDAudio_Pause);
	Cmd_AddCommand("music_resume", CDAudio_Resume);
	return true;
}

void CDAudio_Pause()
{ MIX_PauseAllTracks(mixer); }

void CDAudio_Resume()
{ MIX_ResumeAllTracks(mixer); }

void CDAudio_Shutdown()
{ MIX_Quit(); }

#endif // AVAIL_SDL3MIXER
