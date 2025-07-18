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
	{MIX_INIT_OGG, "OGG", 1, (const char *[]){".ogg"}, true},
	{MIX_INIT_OPUS, "OPUS", 2, (const char *[]){".ogg", ".opus"}, true},
	{MIX_INIT_MP3, "MP3", 1, (const char *[]){".mp3"}, true},
	{MIX_INIT_FLAC, "FLAC", 1, (const char *[]){".flac"}, false},
	{MIX_INIT_MID, "MID", 2, (const char *[]){".mid", ".midi"}, false},
	{MIX_INIT_MOD, "MOD", 1, (const char *[]){".mod"}, false},
	{MIX_INIT_WAVPACK, "WAVPACK", 2, (const char *[]){".wav", ".wv"}, false}
};

static Mix_Music *current_music = NULL;
static s8 current_name[MAX_OSPATH];
static u8 *loaded_file = NULL;

void BGM_Play(s8 *musicname, bool looping)
{
	s8 filename[MAX_OSPATH];
	u8 *file = NULL;
	SDL_IOStream *io = NULL;
	Mix_Music *music = NULL;
	size_t i, j;

	if (!musicname || !*musicname)
	{
		Con_DPrintf("null music file name\n");
		return;
	}

	if (loaded_file) free(loaded_file);
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

	music = Mix_LoadMUS_IO(io, true);
	if (!music)
	{
		Con_Printf("failed to load %s: %s\n", filename, SDL_GetError());
		return;
	}

	CDAudio_Stop();

	current_music = music;

	if (!Mix_PlayMusic(current_music, looping ? -1 : 0))
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
	Mix_Music *music = NULL;
	size_t i, j;

	if (loaded_file) free(loaded_file);
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

	music = Mix_LoadMUS_IO(io, true);
	if (!music)
	{
		Con_Printf("failed to load track %02d: %s\n", (int)track, SDL_GetError());
		return;
	}

	CDAudio_Stop();

	current_music = music;

	if (!Mix_PlayMusic(current_music, looping ? -1 : 0))
	{
		Con_Printf("failed to play track %02d: %s\n", (int)track, SDL_GetError());
		return;
	}

	q_snprintf(current_name, sizeof(current_name), "track%02d", (int)track);
}

void CDAudio_Stop()
{
	if (current_music)
		Mix_FreeMusic(current_music);
	current_music = NULL;
	current_name[0] = 0;
}

void CDAudio_Pause()
{
	Mix_PauseMusic();
}

void CDAudio_Resume()
{
	Mix_ResumeMusic();
}

void CDAudio_Update()
{
	static float last_volume = 0;

	if (bgmvolume.value < 0)
		Cvar_SetQuick(&bgmvolume, "0");
	if (bgmvolume.value > 1)
		Cvar_SetQuick(&bgmvolume, "1");

	if (last_volume != bgmvolume.value)
		Mix_VolumeMusic(bgmvolume.value * MIX_MAX_VOLUME);
}

bool CDAudio_Init()
{
	size_t i;
	for (i = 0; i < ASIZE(music_formats); i++)
	{
		MIX_InitFlags flags = Mix_Init(music_formats[i].flag);
		if (!(flags & music_formats[i].flag))
		{
			Con_Printf("Failed to initialize SDL Mixer format \"%s\": %s\n", music_formats[i].name, SDL_GetError());
			if (music_formats[i].required)
				return false;
		}
		Con_Printf("Initialized SDL Mixer format \"%s\"\n", music_formats[i].name);
	}

	if (!Mix_OpenAudio(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, NULL))
	{
		Con_Printf("Failed to SDL Mixer audio output device: %s\n", SDL_GetError());
		return false;
	}

	Con_Printf("SDL Mixer initialized\n");

	Cmd_AddCommand("music", BGM_Play_f);
	Cmd_AddCommand("music_stop", CDAudio_Stop);
	Cmd_AddCommand("music_pause", CDAudio_Pause);
	Cmd_AddCommand("music_resume", CDAudio_Resume);

	return true;
}

void CDAudio_Shutdown()
{
	CDAudio_Stop();
	Mix_Quit();
}
#endif // AVAIL_SDL3MIXER
