// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

void D_Sky_uv_To_st(s32 u, s32 v, s32 *s, s32 *t)
{
	f32 du = (u - xcenter) / xscale;
	f32 dv = (ycenter - v) / yscale;
	vec3_t end;
	end[0] = vpn[0] + du * vright[0] + dv * vup[0];
	end[1] = vpn[1] + du * vright[1] + dv * vup[1];
	end[2] = vpn[2] + du * vright[2] + dv * vup[2];
	end[2] *= 3;
	VectorNormalize(end);
	f32 temp = skytime * skyspeed;
	*s = (s32)((temp + 6 * (SKYSIZE / 2 - 1) * end[0]) * 65536.0f);
	*t = (s32)((temp + 6 * (SKYSIZE / 2 - 1) * end[1]) * 65536.0f);
}

void D_DrawSkyScans(espan_t *pspan, msurface_t *pface)
{
	s32 skyi = R_SkyIndexForTexture(pface->texinfo->texture);
	do {
		u8 *pdest = (u8 *)((u8 *) d_viewbuffer +
					  (screenwidth * pspan->v) + pspan->u);
		s32 count = pspan->count;
		s32 u = pspan->u; // calculate the initial s & t
		s32 v = pspan->v;
		s32 s, t, snext = 0, tnext = 0;
		D_Sky_uv_To_st(u, v, &s, &t);
		do {
			s32 spancount = count >= SKY_SPAN_MAX ?
				SKY_SPAN_MAX : count;
			count -= spancount;
			s32 sstep = 0;
			s32 tstep = 0;
			if (count) {
				u += spancount;
				// calculate s and t at far end of span,
				// calculate s and t steps across span by shifting
				D_Sky_uv_To_st(u, v, &snext, &tnext);
				sstep = (snext - s) >> SKY_SPAN_SHIFT;
				tstep = (tnext - t) >> SKY_SPAN_SHIFT;
			} else {
				// calculate s and t at last pixel in span,
				// calculate s and t steps across span by division
				s32 spancountminus1 = (f32)(spancount - 1);
				if (spancountminus1 > 0) {
					u += spancountminus1;
					D_Sky_uv_To_st(u, v, &snext, &tnext);
					sstep = (snext - s) / spancountminus1;
					tstep = (tnext - t) / spancountminus1;
				}
			}
			do {
				*pdest++ = r_skysource[skyi][((t & R_SKY_TMASK) >> 8)+
						((s & R_SKY_SMASK) >> 16)];
				s += sstep;
				t += tstep;
			} while (--spancount > 0);
			s = snext;
			t = tnext;
		} while (count > 0);
	} while ((pspan = pspan->pnext) != NULL);
}

void D_DrawSkyScansFog(espan_t *pspan, msurface_t *pface)
{
	R_BuildColorMixLUT(0);
	s32 skyi = R_SkyIndexForTexture(pface->texinfo->texture);
	s32 mix = CLAMP(0, (FOG_LUT_LEVELS - 1)*r_skyfog.value, FOG_LUT_LEVELS);
	do {
		u8 *pdest = (u8 *)((u8 *) d_viewbuffer +
					  (screenwidth * pspan->v) + pspan->u);
		s32 count = pspan->count;
		s32 u = pspan->u; // calculate the initial s & t
		s32 v = pspan->v;
		s32 s, t, snext = 0, tnext = 0;
		D_Sky_uv_To_st(u, v, &s, &t);
		do {
			s32 spancount = count >= SKY_SPAN_MAX ?
				SKY_SPAN_MAX : count;
			count -= spancount;
			s32 sstep = 0;
			s32 tstep = 0;
			if (count) {
				u += spancount;
				// calculate s and t at far end of span,
				// calculate s and t steps across span by shifting
				D_Sky_uv_To_st(u, v, &snext, &tnext);
				sstep = (snext - s) >> SKY_SPAN_SHIFT;
				tstep = (tnext - t) >> SKY_SPAN_SHIFT;
			} else {
				// calculate s and t at last pixel in span,
				// calculate s and t steps across span by division
				s32 spancountminus1 = (f32)(spancount - 1);
				if (spancountminus1 > 0) {
					u += spancountminus1;
					D_Sky_uv_To_st(u, v, &snext, &tnext);
					sstep = (snext - s) / spancountminus1;
					tstep = (tnext - t) / spancountminus1;
				}
			}
			do {
				u8 c = r_skysource[skyi][((t & R_SKY_TMASK) >> 8)+ 
						((s & R_SKY_SMASK) >> 16)];
				*pdest++ = color_mix_lut[c][fog_pal_index][mix];
				s += sstep;
				t += tstep;
			} while (--spancount > 0);
			s = snext;
			t = tnext;
		} while (count > 0);
	} while ((pspan = pspan->pnext) != NULL);
}

void D_DrawSkyScansOnlyFog(espan_t *pspan)
{
	do {
		u8 *pdest = (u8 *)((u8 *) d_viewbuffer +
					  (screenwidth * pspan->v) + pspan->u);
		s32 count = pspan->count;
		do {
			s32 spancount = count>=SKY_SPAN_MAX?SKY_SPAN_MAX:count;
			count -= spancount;
			do {
				*pdest++ = fog_pal_index;
			} while (--spancount > 0);
		} while (count > 0);
	} while ((pspan = pspan->pnext) != NULL);
}
