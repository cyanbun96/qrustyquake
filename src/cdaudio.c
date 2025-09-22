// Copyright (C) 1996-2001 Id Software, Inc.
// Copyright (C) 2002-2009 John Fitzgibbons and others
// Copyright (C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
#include "quakedef.h"

#ifndef AVAIL_SDL3MIXER
static void CDAudio_NotImplemented_f(void)
{
	Con_Printf("not implemented in this build, sorry!\n");
}

void CDAudio_Play(SDL_UNUSED u8 track, SDL_UNUSED bool looping)
{

}

void CDAudio_Stop()
{

}

void CDAudio_Pause()
{

}

void CDAudio_Resume()
{

}

void CDAudio_Update()
{

}

bool CDAudio_Init()
{
	Con_Printf("Music unavailable\n");
	Cmd_AddCommand("music", CDAudio_NotImplemented_f);
	Cmd_AddCommand("music_stop", CDAudio_NotImplemented_f);
	Cmd_AddCommand("music_pause", CDAudio_NotImplemented_f);
	Cmd_AddCommand("music_resume", CDAudio_NotImplemented_f);
	return false;
}

void CDAudio_Shutdown()
{

}
#else
#ifndef ASIZE
#define ASIZE(a) (sizeof(a)/sizeof(a[0]))
#endif

// erysdren: note that i've placed these in a particular order, so that when
// they're iterated in CDAudio_Play() the most common formats for Quake mod
// music will come up first
static struct music_format {
	u32 flag;
	const char *name;
	size_t num_extensions;
	const char **extensions;
	bool required;
} music_formats[] = {
	{0, "OGG", 1, (const char *[]){".ogg"}, true},
	{0, "OPUS", 2, (const char *[]){".ogg", ".opus"}, true},
	{0, "MP3", 1, (const char *[]){".mp3"}, true},
	{0, "FLAC", 1, (const char *[]){".flac"}, false},
	{0, "MID", 2, (const char *[]){".mid", ".midi"}, false},
	{0, "MOD", 1, (const char *[]){".mod"}, false},
	{0, "WAVPACK", 2, (const char *[]){".wav", ".wv"}, false}
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
	size_t i, j;

	if (!musicname || !*musicname)
	{
		Con_DPrintf("null music file name\n");
		return;
	}

	CDAudio_Stop();
	for (i = 0; i < ASIZE(music_formats); i++)
	{
		for (j = 0; j < music_formats[i].num_extensions; j++)
		{
			q_snprintf(filename, sizeof(filename), "music/%s%s", musicname, music_formats[i].extensions[j]);
			file = COM_LoadMallocFile(filename, NULL);
			if (file)
				goto found;
		}
	}

	if (!file)
	{
		Con_Printf("file for %s not found\n", filename);
		return;
	}

found:

	loaded_file = file;
	io = SDL_IOFromConstMem(file, com_filesize);
	if (!io)
	{
		Con_Printf("failed to create IOStream for %s: %s\n", filename, SDL_GetError());
		return;
	}

	music = MIX_LoadAudio_IO(mixer, io, true, true);
	if (!music)
	{
		Con_Printf("failed to load %s: %s\n", filename, SDL_GetError());
		return;
	}

	CDAudio_Stop();

	current_music = music;

	if (!MIX_PlayAudio(mixer, current_music))
	{
		Con_Printf("failed to play %s: %s\n", filename, SDL_GetError());
		return;
	}

	q_strlcpy(current_name, musicname, MAX_OSPATH);
}

static void BGM_Play_f()
{
	if (Cmd_Argc() == 2)
		BGM_Play (Cmd_Argv(1), 1);
	else {
		if (current_music) {
			Con_Printf ("Playing %s, use 'music <musicfile>' to change\n", current_name);
		} else Con_Printf ("music <musicfile>\n");
	}
}

void CDAudio_Play(u8 track, bool looping)
{
	char filename[MAX_QPATH];
	u8 *file = NULL;
	SDL_IOStream *io = NULL;
	MIX_Audio *music = NULL;
	size_t i, j;

	CDAudio_Stop();
	for (i = 0; i < ASIZE(music_formats); i++)
	{
		for (j = 0; j < music_formats[i].num_extensions; j++)
		{
			q_snprintf(filename, sizeof(filename), "music/track%02d%s", (int)track, music_formats[i].extensions[j]);
			file = COM_LoadMallocFile(filename, NULL);
			if (file)
				goto found;
		}
	}

	if (!file)
	{
		Con_Printf("file for track %02d not found\n", (int)track);
		return;
	}

found:

	loaded_file = file;
	io = SDL_IOFromConstMem(file, com_filesize);
	if (!io)
	{
		Con_Printf("failed to create IOStream for track %02d: %s\n", (int)track, SDL_GetError());
		return;
	}

	music = MIX_LoadAudio_IO(mixer, io, true, true);
	if (!music)
	{
		Con_Printf("failed to load track %02d: %s\n", (int)track, SDL_GetError());
		return;
	}

	current_music = music;

	if (!MIX_PlayAudio(mixer, current_music))
	{
		Con_Printf("failed to play track %02d: %s\n", (int)track, SDL_GetError());
		return;
	}

	q_snprintf(current_name, sizeof(current_name), "track%02d", (int)track);
}

void CDAudio_Stop()
{
	MIX_StopAllTracks(mixer, 0);
	if (loaded_file) free(loaded_file);
	loaded_file = NULL;
	memset(current_name, 0, sizeof(current_name));
	current_music = NULL;
}

void CDAudio_Pause()
{
	MIX_PauseAllTracks(mixer);
}

void CDAudio_Resume()
{
	MIX_ResumeAllTracks(mixer);
}

void CDAudio_Update()
{
	if (bgmvolume.value < 0)
		Cvar_SetQuick(&bgmvolume, "0");
	if (bgmvolume.value > 1)
		Cvar_SetQuick(&bgmvolume, "1");
	if (last_volume != bgmvolume.value) {
		last_volume = bgmvolume.value;
		MIX_SetMasterGain(mixer, bgmvolume.value);
	}
}

bool CDAudio_Init()
{
	size_t i;
	MIX_Init();
	mixer = MIX_CreateMixerDevice(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, 0);
	/*for (i = 0; i < ASIZE(music_formats); i++)
	{
		MIX_InitFlags flags = MIX_Init(music_formats[i].flag);
		if (!(flags & music_formats[i].flag))
		{
			Con_Printf("Failed to initialize SDL Mixer format \"%s\": %s\n", music_formats[i].name, SDL_GetError());
			if (music_formats[i].required)
				return false;
		}
		Con_Printf("Initialized SDL Mixer format \"%s\"\n", music_formats[i].name);
	}

	if (!MIX_PlayAudio(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, NULL))
	{
		Con_Printf("Failed to SDL Mixer audio output device: %s\n", SDL_GetError());
		return false;
	}*/

	if (mixer) Con_Printf("SDL Mixer initialized\n");
	else Con_Printf("SDL Mixer initialization failed: %s\n", SDL_GetError());

	Cmd_AddCommand("music", BGM_Play_f);
	Cmd_AddCommand("music_stop", CDAudio_Stop);
	Cmd_AddCommand("music_pause", CDAudio_Pause);
	Cmd_AddCommand("music_resume", CDAudio_Resume);

	return true;
}

void CDAudio_Shutdown()
{
	MIX_Quit();
}
#endif // AVAIL_SDL3MIXER
