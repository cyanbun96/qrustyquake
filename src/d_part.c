// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// d_part.c: software driver module for drawing particles
#include "quakedef.h"

void D_DrawParticle(particle_t *pparticle)
{
	vec3_t local, transformed;
	VectorSubtract(pparticle->org, r_origin, local); // transform point
	transformed[0] = DotProduct(local, r_pright);
	transformed[1] = DotProduct(local, r_pup);
	transformed[2] = DotProduct(local, r_ppn);
	if(transformed[2] < PARTICLE_Z_CLIP)
		return;
	f32 zi = 1.0 / transformed[2]; // project the point
	s32 u = (s32)(xcenter + zi * transformed[0] + 0.5);
	s32 v = (s32)(ycenter - zi * transformed[1] + 0.5);
	s32 izi = (s32)(zi * 0x8000);
	s32 pix;
	if(r_particlesize.value){
		pix = r_particlesize.value;
	} else {
		pix = izi >> d_pix_shift;
		if(pix < d_pix_min) pix = d_pix_min;
		else if(pix > d_pix_max) pix = d_pix_max;
	}
	if(r_particlescale.value > 0)
		pix *= r_particlescale.value;
	s32 half = pix >> 1;
	s32 x0 = u - half;
	s32 y0 = v - half;
	s32 x1 = x0 + pix;
	s32 y1 = y0 + pix;
	if(x1 <= d_vrectx || x0 >= r_refdef.vrectright ||
	    y1 <= d_vrecty || y0 >= r_refdef.vrectbottom)
		return;
	s32 cx0 = x0 < d_vrectx ? d_vrectx : x0;
	s32 cy0 = y0 < d_vrecty ? d_vrecty : y0;
	s32 cx1 = x1 > r_refdef.vrectright ? r_refdef.vrectright : x1;
	s32 cy1 = y1 > r_refdef.vrectbottom ? r_refdef.vrectbottom : y1;
	s32 draw_w = cx1 - cx0;
	s32 draw_h = cy1 - cy0;
	s16 *pz = d_pzbuffer +(d_zwidth * cy0) + cx0;
	u8 *pdest = d_viewbuffer + d_scantable[cy0] + cx0;
	if(r_particlestyle.value == 1){ // circle
		s32 cy = cy0;
		s32 dy = cy - v;
		s32 err = dy*dy;
		s32 rr = half*half;
		for(s32 count = draw_h; count;
		    count--, cy++, pz += d_zwidth, pdest += screenwidth){
			dy = cy - v;
			err = dy*dy;
			if(err > rr) continue;
			s32 dx = (s32)(sqrtf((f32)(rr - err)));
			s32 sx = u - dx;
			s32 ex = u + dx + 1;
			if(sx < cx0) sx = cx0;
			if(ex > cx1) ex = cx1;
			for(s32 x = sx - cx0; x < ex - cx0; x++){
				if(pz[x] <= izi){
					pz[x] = izi;
					pdest[x] = pparticle->color;
				}
			}
		}
	} else { // square
		for(s32 count = draw_h;
		    count; count--, pz += d_zwidth, pdest += screenwidth){
			for(s32 i = 0; i < draw_w; i++){
				if(pz[i] <= izi){
					pz[i] = izi;
					pdest[i] = pparticle->color;
				}
			}
		}
	}
}
