// Copyright(C) 1999-2005 Id Software, Inc.
// Copyright(C) 2005-2012 O.Sezer <sezero@users.sourceforge.net>
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
// snd_sdl.c - SDL audio driver for Hexen II: Hammer of Thyrion(uHexen2)
// based on implementations found in the quakeforge and ioquake3 projects.
#include "quakedef.h"

static s32 buffersize;

static void SDLCALL paint_audio(SDL_UNUSED void *unused, Uint8 *stream, s32 len)
{
	if(!shm) { // shouldn't happen, but just in case
		memset(stream, 0, len);
		return;
	}
	s32 pos = (shm->samplepos * (shm->samplebits / 8));
	if(pos >= buffersize)
		shm->samplepos = pos = 0;
	s32 tobufend = buffersize - pos; // bytes to buffer's end
	s32 len1 = len;
	s32 len2 = 0;
	if(len1 > tobufend) {
		len1 = tobufend;
		len2 = len - len1;
	}
	memcpy(stream, shm->buffer + pos, len1);
	if(len2 <= 0) {
		shm->samplepos += (len1 / (shm->samplebits / 8));
	} else { // wraparound?
		if(!stream || !len1 || !len2 || !shm->buffer)
		{memset(stream,0,len);return;} // fix for WINE
		memcpy(stream + len1, shm->buffer, len2);
		shm->samplepos = (len2 / (shm->samplebits / 8));
	}
	if(shm->samplepos >= buffersize)
		shm->samplepos = 0;
}

void SDLCALL paint_audio_new(void *userdata, SDL_AudioStream *stream, s32 additional_amount, SDL_UNUSED s32 total_amount)
{
	if (additional_amount > 0) {
		Uint8 *data = SDL_stack_alloc(Uint8, additional_amount);
		if (data) {
			paint_audio(userdata, data, additional_amount);
			SDL_PutAudioStreamData(stream, data, additional_amount);
			SDL_stack_free(data);
		}
	}
}

bool SNDDMA_Init(dma_t *dma)
{
	s8 drivername[128];
	if(!SDL_InitSubSystem(SDL_INIT_AUDIO)) {
		Con_Printf("Couldn't init SDL audio: %s\n", SDL_GetError());
		return 0;
	}
	const SDL_AudioSpec desired = { SDL_AUDIO_U8, 2, 11025 }; // TODO configurable sample rates
	s32 samples = 256;
	SDL_AudioStream *stream = SDL_OpenAudioDeviceStream(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, &desired, paint_audio_new, 0);
	SDL_ResumeAudioDevice(SDL_GetAudioStreamDevice(stream));
	memset((void *) dma, 0, sizeof(dma_t));
	shm = dma;
	// Fill the audio DMA information block
	// Since we passed NULL as the 'obtained' spec to SDL_OpenAudio(),
	// SDL will convert to hardware format for us if needed, hence we
	// directly use the desired values here.
	shm->samplebits = (desired.format & 0xFF); // first u8 of format is bits
	shm->signed8 = (desired.format == SDL_AUDIO_S8);
	shm->speed = desired.freq;
	shm->channels = desired.channels;
	s32 tmp = (samples * desired.channels) * 10;
	if(tmp & (tmp - 1)) { // make it a power of two
		s32 val = 1;
		while(val < tmp) val <<= 1;
		tmp = val;
	}
	shm->samples = tmp;
	shm->samplepos = 0;
	shm->submission_chunk = 1;
	Con_Printf("SDL audio spec: %d Hz, %d samples, %d channels\n",
			desired.freq, samples, desired.channels);
	const s8 *driver = SDL_GetCurrentAudioDriver();
	const s8 *device = SDL_GetAudioDeviceName(0);
	q_snprintf(drivername, sizeof(drivername), "%s - %s",
			driver != NULL ? driver : (s8*)"(UNKNOWN)",
			device != NULL ? device : (s8*)"(UNKNOWN)");
	buffersize = shm->samples * (shm->samplebits / 8);
	Con_Printf("SDL audio driver: %s, %d bytes buffer\n",
			drivername, buffersize);
	shm->buffer = (u8 *) calloc(1, buffersize);
	if(!shm->buffer) {
		//SDL_CloseAudio();
		SDL_QuitSubSystem(SDL_INIT_AUDIO);
		shm = NULL;
		Con_Printf("Failed allocating memory for SDL audio\n");
		return 0;
	}
	//SDL_PauseAudio(0);
	return 1;
}

void SNDDMA_Shutdown()
{
	if(shm) {
		Con_Printf("Shutting down SDL sound\n");
		//SDL_CloseAudio();
		SDL_QuitSubSystem(SDL_INIT_AUDIO);
		if(shm->buffer) free(shm->buffer);
		shm->buffer = NULL;
		shm = NULL;
	}
}

s32 SNDDMA_GetDMAPos() { return shm->samplepos; }
void SNDDMA_LockBuffer() { /*SDL_LockAudio();*/ }
void SNDDMA_Submit() { /*SDL_UnlockAudio();*/ }
void SNDDMA_BlockSound() { /*SDL_PauseAudio(1);*/ }
void SNDDMA_UnblockSound(){ /*SDL_PauseAudio(0);*/ }

