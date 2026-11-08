// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

static Q_TLS vec3_t sky_view[3]; // vpn, vright, vup of the surface

static void D_Sky_uv_To_st(s32 u, s32 v, s32 *s, s32 *t)
{
	const f32 *fwd = sky_view[0], *rt = sky_view[1], *up = sky_view[2];
	vec3_t end;
	if(r_refdef.fov_x == 90 && !r_fovmode.value){ // WinQuake's dome
		f32 size = q_max(r_refdef.vrect.width, r_refdef.vrect.height);
		f32 wu = 8192.0f * (f32)(u - ((s32)vid.width >> 1)) / size;
		f32 wv = 8192.0f * (f32)(((s32)vid.height >> 1) - v) / size;
		end[0] = 4096 * fwd[0] + wu * rt[0] + wv * up[0];
		end[1] = 4096 * fwd[1] + wu * rt[1] + wv * up[1];
		end[2] = 4096 * fwd[2] + wu * rt[2] + wv * up[2];
	}else{
		f32 du = (u - xcenter) / xscale;
		f32 dv = (ycenter - v) / yscale;
		end[0] = fwd[0] + du * rt[0] + dv * up[0];
		end[1] = fwd[1] + du * rt[1] + dv * up[1];
		end[2] = fwd[2] + du * rt[2] + dv * up[2];
	}
	end[2] *= 3;
	VectorNormalize(end);
	f32 temp = skytime * skyspeed;
	*s = (s32)((temp + 6 * (SKYSIZE / 2 - 1) * end[0]) * 65536.0f);
	*t = (s32)((temp + 6 * (SKYSIZE / 2 - 1) * end[1]) * 65536.0f);
}

Q_SIL void D_SkyScans(espan_t *pspan, msurface_t *pface, s32 fog,
		s32 smooth, s32 blend)
{ // smooth samples the layers apart, the front one at a fractional offset
  // blend mixes the front 50% over the back instead of masking (Q64)
	sky_t *sky = &r_skies[R_SkyIndexForTexture(pface->texinfo->texture)];
	u8 *skysrc = smooth ? sky->source + 128 : sky->source;
	u8 *front = sky->bottomsky;
	s32 shift = smooth ? (s32)(skytime * skyspeed * 65536.0f) : 0;
	s32 mix = 0; // D_SetupSky built the mix table
	if(fog)
		mix = CLAMP(0, (FOG_LUT_LEVELS-1)*map_skyfog,
			FOG_LUT_LEVELS-1);
	do{
		u8 *pdest = (u8*)d_viewbuffer+(screenwidth*pspan->v)+pspan->u;
		s32 count = pspan->count;
		s32 u = pspan->u; // calculate the initial s & t
		s32 v = pspan->v;
		s32 s, t, snext = 0, tnext = 0;
		D_Sky_uv_To_st(u, v, &s, &t);
		do{
			s32 spancount = count>=SKY_SPAN_MAX?SKY_SPAN_MAX:count;
			count -= spancount;
			s32 sstep = 0;
			s32 tstep = 0;
			if(count){
				u += spancount;
				// calculate s and t at far end of span,
				// s and t steps across span by shifting
				D_Sky_uv_To_st(u, v, &snext, &tnext);
				sstep = (snext - s) >> SKY_SPAN_SHIFT;
				tstep = (tnext - t) >> SKY_SPAN_SHIFT;
			}else{
				// calculate s and t at last pixel in span,
				// s and t steps across span by division
				s32 spancountminus1 = (f32)(spancount - 1);
				if(spancountminus1 > 0){
					u += spancountminus1;
					D_Sky_uv_To_st(u, v, &snext, &tnext);
					sstep = (snext - s) / spancountminus1;
					tstep = (tnext - t) / spancountminus1;
				}
			}
			do{
				u8 c = skysrc[((t & R_SKY_TMASK)
						>> 8)+((s&R_SKY_SMASK) >> 16)];
				if(smooth){
					u8 f = front[(((t+shift)&R_SKY_TMASK)
						>> 9) + (((s+shift)
						&R_SKY_SMASK) >> 16)];
					c = blend ? color_mix_lut[f][c]
						[FOG_LUT_LEVELS/2] : f ? f : c;
				}
				if(fog) c = color_mix_lut[c]
					[fog_pal_index][mix];
				*pdest++ = c;
				s += sstep;
				t += tstep;
			}while(--spancount > 0);
			s = snext;
			t = tnext;
		}while(count > 0);
	}while((pspan = pspan->pnext) != NULL);
}

static bool D_SkyBlend(msurface_t *pface)
{ return r_skies[R_SkyIndexForTexture(pface->texinfo->texture)].blend; }

static void D_DrawSkyScansOnlyFog(espan_t *pspan)
{
	do{
		u8 *pdest = (u8 *)((u8 *) d_viewbuffer +
					  (screenwidth * pspan->v) + pspan->u);
		s32 count = pspan->count;
		do{
			s32 spancount = count>=SKY_SPAN_MAX?SKY_SPAN_MAX:count;
			count -= spancount;
			do{
				*pdest++ = fog_pal_index;
			}while(--spancount > 0);
		}while(count > 0);
	}while((pspan = pspan->pnext) != NULL);
}

// fog: 0 none, 1 the sky fogged by map_skyfog, 2 all fog. view is the vpn,
// vright and vup of the surface, brush entities have their own
void D_DrawSkyScans(espan_t *pspan, msurface_t *pface, const vec3_t *view,
		s32 fog)
{
	for(s32 i = 0; i < 3; i++)
		VectorCopy(view[i], sky_view[i]);
	if(fog == 2) D_DrawSkyScansOnlyFog(pspan);
	else if(fog){
		if(!r_skysmooth.value) D_SkyScans(pspan, pface, 1, 0, 0);
		else if(D_SkyBlend(pface)) D_SkyScans(pspan, pface, 1, 1, 1);
		else D_SkyScans(pspan, pface, 1, 1, 0);
	}else{
		if(!r_skysmooth.value) D_SkyScans(pspan, pface, 0, 0, 0);
		else if(D_SkyBlend(pface)) D_SkyScans(pspan, pface, 0, 1, 1);
		else D_SkyScans(pspan, pface, 0, 1, 0);
	}
}
