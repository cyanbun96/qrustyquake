// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.

// sound.h -- client sound i/o functions

#define DEFAULT_SOUND_PACKET_VOLUME 255
#define DEFAULT_SOUND_PACKET_ATTENUATION 1.0
#define MAX_CHANNELS 128
#define MAX_DYNAMIC_CHANNELS 8

typedef struct
{
	int left;
	int right;
} portable_samplepair_t;

typedef struct sfx_s
{
	char name[MAX_QPATH];
	cache_user_t cache;
} sfx_t;

typedef struct
{
	int length;
	int loopstart;
	int speed;
	int width;
	int stereo;
	byte data[1]; // variable sized
} sfxcache_t;

typedef struct
{
	qboolean gamealive;
	qboolean soundalive;
	qboolean splitbuffer;
	int channels;
	int samples; // mono samples in buffer
	int submission_chunk; // don't mix less than this #
	int samplepos; // in mono samples
	int samplebits;
	int speed;
	unsigned char *buffer;
} dma_t;

typedef struct
{
	sfx_t *sfx; // sfx number
	int leftvol; // 0-255 volume
	int rightvol; // 0-255 volume
	int end; // end time in global paintsamples
	int pos; // sample position in sfx
	int looping; // where to loop, -1 = no looping
	int entnum; // to allow overriding a specific sound
	int entchannel; //
	vec3_t origin; // origin of sound effect
	vec_t dist_mult; // distance multiplier (attenuation/clipK)
	int master_vol; // 0-255 master volume
} channel_t;

typedef struct
{
	int rate;
	int width;
	int channels;
	int loopstart;
	int samples;
	int dataofs; // chunk starts this many bytes from file start
} wavinfo_t;

// User-setable variables
extern channel_t channels[MAX_CHANNELS];
// 0 to MAX_DYNAMIC_CHANNELS-1 = normal entity sounds
// MAX_DYNAMIC_CHANNELS to MAX_DYNAMIC_CHANNELS + NUM_AMBIENTS -1 = water, etc
// MAX_DYNAMIC_CHANNELS + NUM_AMBIENTS to total_channels = static sounds
extern int total_channels;
// Fake dma is a synchronous faking of the DMA progress used for
// isolating performance in the renderer. The fakedma_updates is
// number of times S_Update() is called per second.
extern qboolean fakedma;
extern int fakedma_updates;
extern int paintedtime;
extern vec3_t listener_origin;
extern vec3_t listener_forward;
extern vec3_t listener_right;
extern vec3_t listener_up;
extern volatile dma_t *shm;
extern volatile dma_t sn;
extern vec_t sound_nominal_clip_dist;
extern cvar_t loadas8bit;
extern cvar_t bgmvolume;
extern cvar_t volume;
extern qboolean snd_initialized;
extern int snd_blocked;

void S_LocalSound(char *s);
sfxcache_t *S_LoadSound(sfx_t *s);
wavinfo_t GetWavinfo(char *name, byte *wav, int wavlength);
void SND_InitScaletable();
void SNDDMA_Submit();
void S_Init();
void S_Startup();
void S_Shutdown();
void S_StartSound(int entnum, int entchannel, sfx_t *sfx, vec3_t origin,
		float fvol, float attenuation);
void S_StaticSound(sfx_t *sfx, vec3_t origin, float vol, float attenuation);
void S_StopSound(int entnum, int entchannel);
void S_StopAllSounds(qboolean clear);
void S_ClearBuffer();
void S_Update(vec3_t origin, vec3_t v_forward, vec3_t v_right, vec3_t v_up);
void S_ExtraUpdate();
sfx_t *S_PrecacheSound(char *sample);
void S_TouchSound(char *sample);
void S_PaintChannels(int endtime);
void S_InitPaintChannels();
channel_t *SND_PickChannel(int entnum, int entchannel); // picks a channel based
// on priorities, empty slots, number of channels
void SND_Spatialize(channel_t *ch); // spatializes a channel
qboolean SNDDMA_Init(); // initializes cycling through a DMA buffer and returns
// information on it
int SNDDMA_GetDMAPos(); // gets the current DMA position
void SNDDMA_Shutdown(); // shutdown the DMA xfer.
