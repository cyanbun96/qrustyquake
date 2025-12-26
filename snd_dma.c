// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2007-2008 Kristian Duske
// Copyright(C) 2010-2011 O. Sezer <sezero@users.sourceforge.net>
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
// snd_dma.c -- main control for any streaming sound output device
#include "quakedef.h"

static void S_Play();
static void S_PlayVol();
static void S_SoundList();
static void S_Update_();
static void S_StopAllSoundsC();
static s32 snd_blocked = 0;
static bool snd_initialized = 0;
static dma_t sn;
static vec3_t listener_origin;
static vec3_t listener_forward;
static vec3_t listener_right;
static vec3_t listener_up;
static sfx_t *known_sfx = NULL; // hunk allocated [MAX_SFX]
static s32 num_sfx;
static sfx_t *ambient_sfx[NUM_AMBIENTS];
static bool sound_started = 0;
static s32 soundtime; // sample PAIRS

static void S_SoundInfo_f()
{
	if(!sound_started || !shm) {
		Con_Printf("sound system not started\n");
		return;
	}
	Con_Printf("%d bit, %s, %d Hz\n", shm->samplebits,
			(shm->channels == 2) ? "stereo" : "mono", shm->speed);
	Con_Printf("%5d samples\n", shm->samples);
	Con_Printf("%5d samplepos\n", shm->samplepos);
	Con_Printf("%5d submission_chunk\n", shm->submission_chunk);
	Con_Printf("%5d total_channels\n", total_channels);
	Con_Printf("%p dma buffer\n", shm->buffer);
}

static void SND_Callback_sfxvolume(SDL_UNUSED cvar_t *cvar) { SND_InitScaletable(); }

static void SND_Callback_snd_filterquality(SDL_UNUSED cvar_t *cvar)
{
	if(snd_filterquality.value < 1 || snd_filterquality.value > 5)
	{
		Con_Printf("snd_filterquality must be between 1 and 5\n");
		Cvar_SetQuick(&snd_filterquality, SND_FILTERQUALITY_DEFAULT);
	}
}

void S_Startup()
{
	if(!snd_initialized) return;
	sound_started = SNDDMA_Init(&sn);
	if(!sound_started) Con_Printf("Failed initializing sound\n");
	else Con_Printf("Audio: %d bit, %s, %d Hz\n", shm->samplebits,
			(shm->channels == 2) ? "stereo" : "mono", shm->speed);
}

void S_Init()
{
	s32 i;
	if(snd_initialized) {
		Con_Printf("Sound is already initialized\n");
		return;
	}
	Cvar_RegisterVariable(&nosound);
	Cvar_RegisterVariable(&sfxvolume);
	Cvar_RegisterVariable(&precache);
	Cvar_RegisterVariable(&loadas8bit);
	Cvar_RegisterVariable(&bgmvolume);
	Cvar_RegisterVariable(&ambient_level);
	Cvar_RegisterVariable(&ambient_fade);
	Cvar_RegisterVariable(&snd_noextraupdate);
	Cvar_RegisterVariable(&snd_show);
	Cvar_RegisterVariable(&_snd_mixahead);
	Cvar_RegisterVariable(&sndspeed);
	Cvar_RegisterVariable(&snd_mixspeed);
	Cvar_RegisterVariable(&snd_filterquality);
	if(safemode || COM_CheckParm("-nosound")) return;
	Con_Printf("\nSound Initialization\n");
	Cmd_AddCommand("play", S_Play);
	Cmd_AddCommand("playvol", S_PlayVol);
	Cmd_AddCommand("stopsound", S_StopAllSoundsC);
	Cmd_AddCommand("soundlist", S_SoundList);
	Cmd_AddCommand("soundinfo", S_SoundInfo_f);
	i = COM_CheckParm("-sndspeed");
	if(i && i < com_argc-1) Cvar_SetQuick(&sndspeed, com_argv[i + 1]);
	i = COM_CheckParm("-mixspeed");
	if(i && i < com_argc-1) Cvar_SetQuick(&snd_mixspeed, com_argv[i + 1]);
	if(host_parms.memsize < 0x800000) {
		Cvar_SetQuick(&loadas8bit, "1");
		Con_Printf("loading all sounds as 8bit\n");
	}
	Cvar_SetCallback(&sfxvolume, SND_Callback_sfxvolume);
	Cvar_SetCallback(&snd_filterquality, &SND_Callback_snd_filterquality);
	SND_InitScaletable();
	known_sfx = (sfx_t *) Hunk_AllocName(MAX_SFX*sizeof(sfx_t), "sfx_t");
	num_sfx = 0;
	snd_initialized = 1;
	S_Startup();
	if(sound_started == 0) return;
	ambient_sfx[AMBIENT_WATER] = S_PrecacheSound("ambience/water1.wav");
	ambient_sfx[AMBIENT_SKY] = S_PrecacheSound("ambience/wind2.wav");
	S_StopAllSounds(1);
}

void S_Shutdown()
{
	if(!sound_started) return;
	sound_started = 0;
	snd_blocked = 0;
	SNDDMA_Shutdown();
	shm = NULL;
}

static sfx_t *S_FindName(const s8 *name)
{
	if(!name) Sys_Error("S_FindName: NULL");
	if(strlen(name) >= MAX_QPATH) Sys_Error("Sound name too s64: %s", name);
	s32 i = 0;
	for(; i < num_sfx; i++) // see if already loaded
		if(!strcmp(known_sfx[i].name, name))
			return &known_sfx[i];
	if(num_sfx == MAX_SFX) Sys_Error("S_FindName: out of sfx_t");
	sfx_t *sfx = &known_sfx[i];
	q_strlcpy(sfx->name, name, sizeof(sfx->name));
	num_sfx++;
	return sfx;
}

void S_TouchSound(const s8 *name)
{
	if(!sound_started) return;
	sfx_t *sfx = S_FindName(name);
	Cache_Check(&sfx->cache);
}

sfx_t *S_PrecacheSound(const s8 *name)
{
	if(!sound_started || nosound.value) return NULL;
	sfx_t *sfx = S_FindName(name);
	if(precache.value) S_LoadSound(sfx); // cache it in
	return sfx;
}

channel_t *SND_PickChannel(s32 entnum, s32 entchannel)
{ // picks a channel based on priorities, empty slots, number of channels
	s32 first_to_die = -1; // Check for replacement sound
	s32 life_left = 0x7fffffff; // or find the best one to replace
	for(s32 ch_idx=NUM_AMBIENTS; ch_idx < NUM_AMBIENTS+MAX_DYNAMIC_CHANNELS;
			ch_idx++){
		if(entchannel != 0 // channel 0 never overrides
			&& snd_channels[ch_idx].entnum == entnum
			&& (snd_channels[ch_idx].entchannel == entchannel
				|| entchannel == -1) ) {
			first_to_die = ch_idx;
			break; // always override sound from same entity
		}
		// don't let monster sounds override player sounds
		if(snd_channels[ch_idx].entnum == cl.viewentity
			&& entnum != cl.viewentity && snd_channels[ch_idx].sfx)
			continue;
		if(snd_channels[ch_idx].end - paintedtime < life_left) {
			life_left = snd_channels[ch_idx].end - paintedtime;
			first_to_die = ch_idx;
		}
	}
	if(first_to_die == -1) return NULL;
	if(snd_channels[first_to_die].sfx) snd_channels[first_to_die].sfx=NULL;
	return &snd_channels[first_to_die];
}

void SND_Spatialize(channel_t *ch)
{ // spatializes a channel
	vec_t dot, dist, lscale, rscale, scale;
	vec3_t source_vec;
	// anything coming from the view entity will always be full volume
	if(ch->entnum == cl.viewentity) {
		ch->leftvol = ch->master_vol;
		ch->rightvol = ch->master_vol;
		return;
	}
	// calculate stereo seperation and distance attenuation
	VectorSubtract(ch->origin, listener_origin, source_vec);
	dist = VectorNormalize(source_vec) * ch->dist_mult;
	dot = DotProduct(listener_right, source_vec);
	if(shm->channels == 1) { rscale = 1.0; lscale = 1.0;
	} else { rscale = 1.0 + dot; lscale = 1.0 - dot; }
	// add in distance effect
	scale = (1.0 - dist) * rscale;
	ch->rightvol = (s32) (ch->master_vol * scale);
	if(ch->rightvol < 0) ch->rightvol = 0;
	scale = (1.0 - dist) * lscale;
	ch->leftvol = (s32) (ch->master_vol * scale);
	if(ch->leftvol < 0) ch->leftvol = 0;
}

void S_StartSound(s32 en, s32 ech, sfx_t *sfx, vec3_t origin, f32 fvol, f32 atn)
{ // Start a sound effect
	if(!sound_started) return;
	if(!sfx) return;
	if(nosound.value) return;
	channel_t *target_chan = SND_PickChannel(en, ech); // pick a channel
	if(!target_chan) return;
	memset(target_chan, 0, sizeof(*target_chan)); // spatialize
	VectorCopy(origin, target_chan->origin);
	target_chan->dist_mult = atn / sound_nominal_clip_dist;
	target_chan->master_vol = (s32) (fvol * 255);
	target_chan->entnum = en;
	target_chan->entchannel = ech;
	SND_Spatialize(target_chan);
	if(!target_chan->leftvol && !target_chan->rightvol) return;//not audible
	sfxcache_t *sc = S_LoadSound(sfx); // new channel
	if(!sc) {
		target_chan->sfx = NULL;
		return; // couldn't load the sound's data
	}
	target_chan->sfx = sfx;
	target_chan->pos = 0.0;
	target_chan->end = paintedtime + sc->length;
// if an identical sound has also been started this frame, offset the pos
// a bit to keep it from just making the first one louder
	channel_t *check = &snd_channels[NUM_AMBIENTS];
	for(s32 ch_idx=NUM_AMBIENTS; ch_idx < NUM_AMBIENTS+MAX_DYNAMIC_CHANNELS;
			ch_idx++, check++) {
		if(check == target_chan) continue;
		if(check->sfx == sfx && !check->pos) {
			s32 skip = 0.1 * shm->speed; // 0.1 * sc->speed
			if(skip > sc->length) skip = sc->length;
			if(skip > 0) skip = rand() % skip;
			target_chan->pos += skip;
			target_chan->end -= skip;
			break;
		}
	}
}

void S_StopSound(s32 entnum, s32 entchannel)
{
	for(s32 i = 0; i < MAX_DYNAMIC_CHANNELS; i++)
		if(snd_channels[i].entnum == entnum
			&& snd_channels[i].entchannel == entchannel) {
			snd_channels[i].end = 0;
			snd_channels[i].sfx = NULL;
			return;
		}
}

void S_StopAllSounds(bool clear)
{
	if(!sound_started) return;
	total_channels = MAX_DYNAMIC_CHANNELS + NUM_AMBIENTS; // no statics
	for(s32 i = 0; i < MAX_CHANNELS; i++)
		if(snd_channels[i].sfx) snd_channels[i].sfx = NULL;
	memset(snd_channels, 0, MAX_CHANNELS * sizeof(channel_t));
	if(clear) S_ClearBuffer();
}

static void S_StopAllSoundsC() { S_StopAllSounds(1); }

void S_ClearBuffer()
{
	if(!sound_started || !shm) return;
	SNDDMA_LockBuffer();
	if(! shm->buffer) return;
	s_rawend = 0;
	s32 clear = 0;
	if(shm->samplebits == 8 && !shm->signed8) clear = 0x80;
	memset(shm->buffer, clear, shm->samples * shm->samplebits / 8);
	SNDDMA_Submit();
}

void S_StaticSound(sfx_t *sfx, vec3_t origin, f32 vol, f32 attenuation)
{
	if(!sfx) return;
	if(total_channels == MAX_CHANNELS) {
		Con_Printf("total_channels == MAX_CHANNELS\n");
		return;
	}
	channel_t *ss = &snd_channels[total_channels];
	total_channels++;
	sfxcache_t *sc = S_LoadSound(sfx);
	if(!sc) return;
	if(sc->loopstart == -1) {
		Con_Printf("Sound %s not looped\n", sfx->name);
		return;
	}
	ss->sfx = sfx;
	VectorCopy(origin, ss->origin);
	ss->master_vol = (s32)vol;
	ss->dist_mult = (attenuation / 64) / sound_nominal_clip_dist;
	ss->end = paintedtime + sc->length;
	SND_Spatialize(ss);
}

static void S_UpdateAmbientSounds()
{
	if(cls.state != ca_connected) return; // no ambients when disconnected
	if(!cl.worldmodel) return; // calc ambient sound levels
	mleaf_t *l = Mod_PointInLeaf(listener_origin, cl.worldmodel);
	s32 vol = 0, amb_ch = 0;
	if(!l || !ambient_level.value) {
		for(amb_ch = 0; amb_ch < NUM_AMBIENTS; amb_ch++)
			snd_channels[amb_ch].sfx = NULL;
		return;
	}
	for(amb_ch = 0; amb_ch < NUM_AMBIENTS; amb_ch++) {
		channel_t *chan = &snd_channels[amb_ch];
		chan->sfx = ambient_sfx[amb_ch];
		vol = (s32)(ambient_level.value*l->ambient_sound_level[amb_ch]);
		if(vol < 8) vol = 0;
		if(chan->master_vol < vol) { // don't adjust volume too fast
		    chan->master_vol+=(s32)(host_frametime*ambient_fade.value);
		    if(chan->master_vol > vol) chan->master_vol = vol;
		} else if(chan->master_vol > vol) {
		    chan->master_vol-=(s32)(host_frametime*ambient_fade.value);
		    if(chan->master_vol < vol) chan->master_vol = vol;
		}
		chan->leftvol = chan->rightvol = chan->master_vol;
	}
}

void S_Update(vec3_t origin, vec3_t forward, vec3_t right, vec3_t up)
{ // Called once each time through the main loop
	s32 i, j;
	if(!sound_started || (snd_blocked > 0)) return;
	VectorCopy(origin, listener_origin);
	VectorCopy(forward, listener_forward);
	VectorCopy(right, listener_right);
	VectorCopy(up, listener_up);
	S_UpdateAmbientSounds(); // update general area ambient sound sources
	channel_t *combine = NULL;
	// update spatialization for static and dynamic sounds
	channel_t *ch = snd_channels + NUM_AMBIENTS;
	for(i = NUM_AMBIENTS; i < total_channels; i++, ch++) {
		if(!ch->sfx) continue;
		SND_Spatialize(ch); // respatialize channel
		if(!ch->leftvol && !ch->rightvol) continue;
		// try to combine static sounds with a previous channel of the
		// same sound effect so we don't mix five torches every frame
		if(i < MAX_DYNAMIC_CHANNELS + NUM_AMBIENTS) continue;
		// see if it can just use the last one
		if(combine && combine->sfx == ch->sfx) {
			combine->leftvol += ch->leftvol;
			combine->rightvol += ch->rightvol;
			ch->leftvol = ch->rightvol = 0;
			continue;
		}
		// search for one
		combine = snd_channels + MAX_DYNAMIC_CHANNELS + NUM_AMBIENTS;
		for(j = MAX_DYNAMIC_CHANNELS+NUM_AMBIENTS; j < i; j++,combine++)
			if(combine->sfx == ch->sfx) break;
		if(j == total_channels) combine = NULL;
		if(combine != ch) {
			combine->leftvol += ch->leftvol;
			combine->rightvol += ch->rightvol;
			ch->leftvol = ch->rightvol = 0;
		}
	}
	if(snd_show.value) { // debugging output
		s32 total = 0;
		ch = snd_channels;
		for(i = 0; i < total_channels; i++, ch++) {
			if(ch->sfx && (ch->leftvol || ch->rightvol) ) {
				Con_Printf("%3i %3i %s\n", ch->leftvol,
						ch->rightvol, ch->sfx->name);
				total++;
			}
		}
		Con_Printf("----(%i)----\n", total);
	}
	S_Update_();
}

static void GetSoundtime()
{
	static s32 buffers;
	static s32 oldsamplepos;
	s32 fullsamples = shm->samples / shm->channels;
	// it is possible to miscount buffers if it has wrapped twice between
	// calls to S_Update. Oh well.
	s32 samplepos = SNDDMA_GetDMAPos();
	if(samplepos < oldsamplepos) {
		buffers++; // buffer wrapped
		if(paintedtime > 0x40000000){//chop things to avoid 32bit limits
			buffers = 0;
			paintedtime = fullsamples;
			S_StopAllSounds(1);
		}
	}
	oldsamplepos = samplepos;
	soundtime = buffers*fullsamples + samplepos/shm->channels;
}

void S_ExtraUpdate()
{
	if(snd_noextraupdate.value) return; // don't pollute timings
	S_Update_();
}

static void S_Update_()
{
	u32 endtime;
	s32 samps;
	if(!sound_started || (snd_blocked > 0))
		return;
	SNDDMA_LockBuffer();
	if(! shm->buffer)
		return;
	// Updates DMA time
	GetSoundtime();
	// check to make sure that we haven't overshot
	if(paintedtime < soundtime)
	{
		// Con_Printf("S_Update_ : overflow\n");
		paintedtime = soundtime;
	}
	// mix ahead of current position
	endtime = soundtime + (u32)(_snd_mixahead.value * shm->speed);
	samps = shm->samples >> (shm->channels - 1);
	endtime = q_min(endtime, (u32)(soundtime + samps));
	S_PaintChannels(endtime);
	SNDDMA_Submit();
}

void S_BlockSound()
{
	if(sound_started && snd_blocked == 0){
		snd_blocked = 1;
		S_ClearBuffer();
		if(shm) SNDDMA_BlockSound();
	}
}

void S_UnblockSound()
{
	if(!sound_started || !snd_blocked) return;
	if(snd_blocked == 1){ // --snd_blocked == 0 
		snd_blocked = 0;
		SNDDMA_UnblockSound();
		S_ClearBuffer();
	}
}

static void S_Play()
{
	static s32 hash = 345;
	s8 name[256];
	s32 i = 1;
	while(i < Cmd_Argc()) {
		q_strlcpy(name, Cmd_Argv(i), sizeof(name));
		if(!strrchr(Cmd_Argv(i), '.'))
			q_strlcat(name, ".wav", sizeof(name));
		sfx_t *sfx = S_PrecacheSound(name);
		S_StartSound(hash++, 0, sfx, listener_origin, 1.0, 1.0);
		i++;
	}
}

static void S_PlayVol()
{
	static s32 hash = 543;
	s8 name[256];
	s32 i = 1;
	while(i < Cmd_Argc()) {
		q_strlcpy(name, Cmd_Argv(i), sizeof(name));
		if(!strrchr(Cmd_Argv(i), '.'))
			q_strlcat(name, ".wav", sizeof(name));
		sfx_t *sfx = S_PrecacheSound(name);
		f32 vol = atof(Cmd_Argv(i + 1));
		S_StartSound(hash++, 0, sfx, listener_origin, vol, 1.0);
		i += 2;
	}
}

static void S_SoundList()
{
	s32 i = 0;
	s32 total = 0;
	for(sfx_t *sfx = known_sfx; i < num_sfx; i++, sfx++) {
		sfxcache_t *sc = (sfxcache_t *) Cache_Check(&sfx->cache);
		if(!sc) continue;
		s32 size = sc->length*sc->width*(sc->stereo + 1);
		total += size;
		if(sc->loopstart >= 0) Con_SafePrintf("L");
		else Con_SafePrintf(" ");
		Con_SafePrintf("(%2db) %6i : %s\n",sc->width*8,size,sfx->name);
	}
	Con_Printf("%i sounds, %i bytes\n", num_sfx, total);
}

void S_LocalSound(const s8 *name)
{
	if(nosound.value) return;
	if(!sound_started) return;
	sfx_t *sfx = S_PrecacheSound(name);
	if(!sfx) { Con_Printf("S_LocalSound: can't cache %s\n", name); return; }
	S_StartSound(cl.viewentity, -1, sfx, vec3_origin, 1, 1);
}
