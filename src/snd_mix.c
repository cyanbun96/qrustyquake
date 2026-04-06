// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2010-2011 O. Sezer <sezero@users.sourceforge.net>
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
// snd_mix.c -- portable code to mix sounds for snd_dma.c
#include "quakedef.h"

static portable_samplepair_t paintbuffer[PAINTBUFFER_SIZE];
static s32 snd_scaletable[32][256];
static s32 *snd_p, snd_linear_count;
static s16 *snd_out;
static s32 snd_vol;

static void Snd_WriteLinearBlastStereo16()
{
	for(s32 i = 0; i < snd_linear_count; i += 2) {
		s32 val = snd_p[i] / 256;
		if(val > 0x7fff) snd_out[i] = 0x7fff;
		else if(val < (s16)0x8000) snd_out[i] = (s16)0x8000;
		else snd_out[i] = val;
		val = snd_p[i+1] / 256;
		if(val > 0x7fff) snd_out[i+1] = 0x7fff;
		else if(val < (s16)0x8000) snd_out[i+1] = (s16)0x8000;
		else snd_out[i+1] = val;
	}
}

static void S_TransferStereo16(s32 endtime)
{
	snd_p = (s32 *) paintbuffer;
	s32 lpaintedtime = paintedtime;
	while(lpaintedtime < endtime) {
		// handle recirculating buffer issues
		s32 lpos = lpaintedtime & ((shm->samples >> 1) - 1);
		snd_out = (s16 *)shm->buffer + (lpos << 1);
		snd_linear_count = (shm->samples >> 1) - lpos;
		if(lpaintedtime + snd_linear_count > endtime)
			snd_linear_count = endtime - lpaintedtime;
		snd_linear_count <<= 1;
		// write a linear blast of samples
		Snd_WriteLinearBlastStereo16();
		snd_p += snd_linear_count;
		lpaintedtime += (snd_linear_count >> 1);
	}
}

static void S_TransferPaintBuffer(s32 endtime)
{
	if(shm->samplebits == 16 && shm->channels == 2) {
		S_TransferStereo16(endtime);
		return;
	}
	s32 *p = (s32 *)paintbuffer;
	s32 count = (endtime - paintedtime) * shm->channels;
	s32 out_mask = shm->samples - 1;
	s32 out_idx = paintedtime * shm->channels & out_mask;
	s32 step = 3 - shm->channels;
	if(shm->samplebits == 16) {
		s16 *out = (s16 *)shm->buffer;
		while(count--) {
			s32 val = *p / 256;
			p+= step;
			if(val > 0x7fff) val = 0x7fff;
			else if(val < (s16)0x8000) val = (s16)0x8000;
			out[out_idx] = val;
			out_idx = (out_idx + 1) & out_mask;
		}
	} else if(shm->samplebits == 8 && !shm->signed8) {
		u8 *out = shm->buffer;
		while(count--) {
			s32 val = *p / 256;
			p+= step;
			if(val > 0x7fff) val = 0x7fff;
			else if(val < (s16)0x8000) val = (s16)0x8000;
			out[out_idx] = (val / 256) + 128;
			out_idx = (out_idx + 1) & out_mask;
		}
	} else if(shm->samplebits == 8) {
		s8 *out = (s8 *) shm->buffer;
		while(count--) {
			s32 val = *p / 256;
			p+= step;
			if(val > 0x7fff) val = 0x7fff;
			else if(val < (s16)0x8000) val = (s16)0x8000;
			out[out_idx] = (val / 256);
			out_idx = (out_idx + 1) & out_mask;
		}
	}
}

// Makes a lowpass filter kernel, from equation 16-4 in
// "The Scientist and Engineer's Guide to Digital Signal Processing"
// M is the kernel size(not counting the center point), must be even
// kernel has room for M+1 floats
// f_c is the filter cutoff frequency, as a fraction of the samplerate
static void S_MakeBlackmanWindowKernel(f32 *kernel, s32 M, f32 f_c)
{
	for(s32 i = 0; i <= M; i++) {
		if(i == M/2)
			kernel[i] = 2 * M_PI * f_c;
		else {
			kernel[i]=(sin(2*M_PI*f_c*(i-M/2.0))/(i-(M/2.0)))
				*(0.42-0.5*cos(2*M_PI*i/(f64)M)
				+0.08*cos(4*M_PI*i/(f64)M));
		}
	}
	// normalize the kernel so all of the values sum to 1
	f32 sum = 0;
	for(s32 i = 0; i <= M; i++) sum += kernel[i];
	for(s32 i = 0; i <= M; i++) kernel[i] /= sum;
}

static void S_UpdateFilter(filter_t *filter, s32 M, f32 f_c)
{
	if(filter->f_c != f_c || filter->M != M) {
		if(filter->memory != NULL) free(filter->memory);
		if(filter->kernel != NULL) free(filter->kernel);
		filter->M = M;
		filter->f_c = f_c;
		filter->parity = 0;
		// M + 1 rounded up to the next multiple of 16
		filter->kernelsize = (M + 1) + 16 - ((M + 1) % 16);
		filter->memory = (f32 *)calloc(filter->kernelsize, sizeof(f32));
		filter->kernel = (f32 *)calloc(filter->kernelsize, sizeof(f32));
		S_MakeBlackmanWindowKernel(filter->kernel, M, f_c);
	}
}

// Lowpass-filter the given buffer containing 44100Hz audio.
// As an optimization, it decimates the audio to 11025Hz(setting every sample
// position that's not a multiple of 4 to 0), then convoluting with the filter
// kernel is 4x faster, because we can skip 3/4 of the input samples that are
// known to be 0 and skip 3/4 of the filter kernel.
static void S_ApplyFilter(filter_t *filter, s32 *data, s32 stride, s32 count)
{
	const s32 kernelsize = filter->kernelsize;
	const f32 *kernel = filter->kernel;
	f32 *input = (f32 *) malloc(sizeof(f32) * (filter->kernelsize + count));
	// set up the input buffer
	// memory holds the previous filter->kernelsize samples of input.
	memcpy(input, filter->memory, filter->kernelsize * sizeof(f32));
	for(s32 i = 0; i < count; i++)
		input[filter->kernelsize+i] = data[i*stride]/(32768.0*256.0);
	//copy out the last filter->kernelsize samples to 'memory' for next time
	memcpy(filter->memory, input + count, filter->kernelsize * sizeof(f32));
	// apply the filter
	s32 parity = filter->parity;
	for(s32 i = 0; i < count; i++) {
		const f32 *input_plus_i = input + i;
		f32 val[4] = {0, 0, 0, 0};
		for(s32 j = (4 - parity) % 4; j < kernelsize; j+=16) {
			val[0] += kernel[j] * input_plus_i[j];
			val[1] += kernel[j+4] * input_plus_i[j+4];
			val[2] += kernel[j+8] * input_plus_i[j+8];
			val[3] += kernel[j+12] * input_plus_i[j+12];
		}
		// 4.0 factor is to increase volume by 12 dB; this is to make up
		// the volume drop caused by the zero-filling this filter does.
		data[i*stride] = (val[0]+val[1]+val[2]+val[3])*(32768.0*1024.0);
		parity = (parity + 1) % 4;
	}
	filter->parity = parity;
	free(input);
}

// lowpass filters 24-bit integer samples in 'data' (stored in 32-bit ints).
// assumes 44100Hz sample rate, and lowpasses at around 5kHz
// memory should be a zero-filled filter_t struct
static void S_LowpassFilter(s32 *data, s32 stride, s32 count, filter_t *memory)
{
	s32 M;
	f32 bw;
	switch((s32)snd_filterquality.value) {
	case 1: M = 126; bw = 0.900; break;
	case 2: M = 150; bw = 0.915; break;
	case 3: M = 174; bw = 0.930; break;
	case 4: M = 198; bw = 0.945; break;
	case 5: default: M = 222; bw = 0.960; break;
	}
	f32 f_c = (bw * 11025 / 2.0) / 44100.0;
	S_UpdateFilter(memory, M, f_c);
	S_ApplyFilter(memory, data, stride, count);
}

static void SND_PaintChannelFrom8
	(channel_t *ch, sfxcache_t *sc, s32 endtime, s32 paintbufferstart);
static void SND_PaintChannelFrom16
	(channel_t *ch, sfxcache_t *sc, s32 endtime, s32 paintbufferstart);

void S_PaintChannels(s32 endtime)
{
	snd_vol = sfxvolume.value * 256;
	while(paintedtime < endtime) {
		s32 end = endtime; // if paintbuffer is smaller than DMA buffer
		if(endtime - paintedtime > PAINTBUFFER_SIZE)
			end = paintedtime + PAINTBUFFER_SIZE;
		// clear the paint buffer
		memset(paintbuffer, 0, (end - paintedtime)
				* sizeof(portable_samplepair_t));
		// paint in the channels.
		channel_t *ch = snd_channels;
		for(s32 i = 0; i < total_channels; i++, ch++) {
			if(!ch->sfx) continue;
			if(!ch->leftvol && !ch->rightvol) continue;
			sfxcache_t *sc = S_LoadSound(ch->sfx);
			if(!sc) continue;
			s32 ltime = paintedtime;
			while(ltime < end) { // paint up to end
				s32 count;
				if(ch->end < end) count = ch->end - ltime;
				else count = end - ltime;
				if(count > 0) {
			  // the last param to SND_PaintChannelFrom is the index
			  // to start painting to in the paintbuffer, usually 0.
					if(sc->width ==1) SND_PaintChannelFrom8(
					    ch, sc, count, ltime - paintedtime);
					else SND_PaintChannelFrom16(
					    ch, sc, count, ltime - paintedtime);
					ltime += count;
				}
				// if at end of loop, restart
				if(ltime >= ch->end) {
					if(sc->loopstart >= 0) {
					 ch->pos = sc->loopstart;
					 ch->end = ltime + sc->length - ch->pos;
					} else { // channel just stopped
						ch->sfx = NULL;
						break;
					}
				}
			}
		}
		// clip each sample to 0dB, then reduce by 6dB(to leave some
		// headroom for the lowpass filter and the music). the lowpass
		// will smooth out the clipping
		for(s32 i = 0; i < end-paintedtime; i++) {
			paintbuffer[i].left = CLAMP(-32768 * 256,
					paintbuffer[i].left, 32767 * 256) / 2;
			paintbuffer[i].right = CLAMP(-32768 * 256,
					paintbuffer[i].right, 32767 * 256) / 2;
		}
		// apply a lowpass filter
		if(sndspeed.value == 11025 && shm->speed == 44100) {
			static filter_t memory_l, memory_r;
			S_LowpassFilter((s32 *)paintbuffer, 2,
					end - paintedtime, &memory_l);
			S_LowpassFilter(((s32 *)paintbuffer) + 1, 2,
					end - paintedtime, &memory_r);
		}
		// paint in the music
		if(s_rawend >= paintedtime) { //copy from streaming sound source
			s32 stop = (end < s_rawend) ? end : s_rawend;
			for(s32 i = paintedtime; i < stop; i++) {
				s32 s = i & (MAX_RAW_SAMPLES - 1);
				// lower music by 6db to match sfx
				paintbuffer[i - paintedtime].left +=
					s_rawsamples[s].left / 2;
				paintbuffer[i - paintedtime].right +=
					s_rawsamples[s].right / 2;
			}
		}
		S_TransferPaintBuffer(end); // transfer according to DMA format
		paintedtime = end;
	}
}

void SND_InitScaletable()
{
	for(s32 i = 0; i < 32; i++) {
		s32 scale = i * 8 * 256 * sfxvolume.value;
		for(s32 j = 0; j < 256; j++)
			snd_scaletable[i][j] = ((j<128) ? j : j - 256) * scale;
	}
}

static void SND_PaintChannelFrom8
	(channel_t *ch, sfxcache_t *sc, s32 count, s32 paintbufferstart)
{
	if(ch->leftvol > 255) ch->leftvol = 255;
	if(ch->rightvol > 255) ch->rightvol = 255;
	s32 *lscale = snd_scaletable[ch->leftvol >> 3];
	s32 *rscale = snd_scaletable[ch->rightvol >> 3];
	u8 *sfx = (u8 *)sc->data + ch->pos;
	for(s32 i = 0; i < count; i++) {
		s32 data = sfx[i];
		paintbuffer[paintbufferstart + i].left += lscale[data];
		paintbuffer[paintbufferstart + i].right += rscale[data];
	}
	ch->pos += count;
}

static void SND_PaintChannelFrom16
	(channel_t *ch, sfxcache_t *sc, s32 count, s32 paintbufferstart)
{
	s32 leftvol = ch->leftvol * snd_vol;
	s32 rightvol = ch->rightvol * snd_vol;
	leftvol /= 256;
	rightvol /= 256;
	s16 *sfx = (s16 *)sc->data + ch->pos;
	for(s32 i = 0; i < count; i++) {
		s32 data = sfx[i];
		// this was causing integer overflow as observed in quakespasm
		// with the warpspasm mod moved >>8 to left/right volume above.
		s32 left = data * leftvol;
		s32 right = data * rightvol;
		paintbuffer[paintbufferstart + i].left += left;
		paintbuffer[paintbufferstart + i].right += right;
	}
	ch->pos += count;
}
