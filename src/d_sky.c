// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

void D_Sky_uv_To_st(s32 u, s32 v, s32 *s, s32 *t)
{
	f32 temp = r_refdef.vrect.width >= r_refdef.vrect.height ?
		(f32)r_refdef.vrect.width : (f32)r_refdef.vrect.height;
	f32 wu = 8192.0 * (f32)(u - ((s32)vid.width >> 1)) / temp;
	f32 wv = 8192.0 * (f32)(((s32)vid.height >> 1) - v) / temp;
	vec3_t end;
	end[0] = 4096 * vpn[0] + wu * vright[0] + wv * vup[0];
	end[1] = 4096 * vpn[1] + wu * vright[1] + wv * vup[1];
	end[2] = 4096 * vpn[2] + wu * vright[2] + wv * vup[2];
	end[2] *= 3;
	VectorNormalize(end);
	temp = skytime * skyspeed; // TODO: add D_SetupFrame & set this there
	*s = (s32)((temp + 6 * (SKYSIZE / 2 - 1) * end[0]) * 0x10000);
	*t = (s32)((temp + 6 * (SKYSIZE / 2 - 1) * end[1]) * 0x10000);
}

void D_DrawSkyScans(espan_t *pspan)
{
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
				*pdest++ = r_skysource[((t & R_SKY_TMASK) >> 8)+
						((s & R_SKY_SMASK) >> 16)];
				s += sstep;
				t += tstep;
			} while (--spancount > 0);
			s = snext;
			t = tnext;
		} while (count > 0);
	} while ((pspan = pspan->pnext) != NULL);
}
