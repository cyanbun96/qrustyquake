// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// d_part.c: software driver module for drawing particles
#include "quakedef.h"

typedef struct { // a particle on the screen, see D_ProjectParticle
	s32 cx0, cy0, cx1, cy1; // the square, clipped to the view
	s32 u, v, half, izi;
	s32 color;
} partdraw_t;
static partdraw_t *d_parts; // the frame's particles, for the workers
static s32 d_nparts, d_partcap;

#define PART_CIRCLE(TEST, COLOR) { \
	s32 cy = cy0; \
	s32 rr = half*half; \
	for(s32 count = draw_h; count; \
	    count--, cy++, pz += d_zwidth, pdest += screenwidth){ \
		s32 dy = cy - v; \
		s32 err = dy*dy; \
		if(err > rr) continue; \
		s32 dx = (s32)(sqrtf((f32)(rr - err))); \
		s32 sx = u - dx; \
		s32 ex = u + dx + 1; \
		if(sx < cx0) sx = cx0; \
		if(ex > cx1) ex = cx1; \
		for(s32 x = sx - cx0; x < ex - cx0; x++){ \
			if(TEST){ \
				pz[x] = izi; \
				pdest[x] = COLOR; \
			} \
		} \
	} \
}
#define PART_SQUARE(TEST, COLOR) { \
	for(s32 count = draw_h; \
	    count; count--, pz += d_zwidth, pdest += screenwidth){ \
		for(s32 x = 0; x < draw_w; x++){ \
			if(TEST){ \
				pz[x] = izi; \
				pdest[x] = COLOR; \
			} \
		} \
	} \
}

static s32 D_ParticleAlpha()
{
	s32 part_alpha = r_particlealpha.value * FOG_LUT_LEVELS;
	return part_alpha > FOG_LUT_LEVELS ? FOG_LUT_LEVELS : part_alpha;
}

// Where the particle lands, false if it's not on the screen
static bool D_ProjectParticle(const particle_t *pparticle, partdraw_t *pd)
{
	s32 part_alpha = D_ParticleAlpha();
	if(part_alpha <= 0)
		return 0;
	vec3_t local, transformed;
	VectorSubtract(pparticle->org, r_origin, local); // transform point
	transformed[0] = DotProduct(local, r_pright);
	transformed[1] = DotProduct(local, r_pup);
	transformed[2] = DotProduct(local, r_ppn);
	if(transformed[2] < 8.0) // particle_z_clip
		return 0;
	f32 zi = 1.0 / transformed[2]; // project the point
	s32 u = (s32)(xcenter + zi * transformed[0] + 0.5);
	s32 v = (s32)(ycenter - zi * transformed[1] + 0.5);
	s32 izi = (s32)(zi * 0x8000 * 0x10000);
	s32 pix;
	if(r_particlesize.value){
		pix = r_particlesize.value;
	}else{
		pix = (izi >> 16) >> d_pix_shift;
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
		return 0;
	pd->cx0 = x0 < d_vrectx ? d_vrectx : x0;
	pd->cy0 = y0 < d_vrecty ? d_vrecty : y0;
	pd->cx1 = x1 > r_refdef.vrectright ? r_refdef.vrectright : x1;
	pd->cy1 = y1 > r_refdef.vrectbottom ? r_refdef.vrectbottom : y1;
	pd->u = u;
	pd->v = v;
	pd->half = half;
	pd->izi = izi;
	pd->color = pparticle->color;
	if(part_alpha < FOG_LUT_LEVELS && !r_alphastyle.value
			&& !fog_lut_built)
		R_BuildColorMixLUT(0);
	return 1;
}

// The rows of the particle that this thread draws
static void D_FillParticle(const partdraw_t *pd)
{
	s32 part_alpha = D_ParticleAlpha();
	s32 level = D_DitherLevel((f32)part_alpha / FOG_LUT_LEVELS);
	s32 cx0 = pd->cx0, cx1 = pd->cx1;
	s32 cy0 = q_max(pd->cy0, d_bandy0);
	s32 cy1 = q_min(pd->cy1, d_bandy1);
	if(cy0 >= cy1)
		return;
	s32 u = pd->u, v = pd->v, half = pd->half, izi = pd->izi;
	s32 draw_w = cx1 - cx0;
	s32 draw_h = cy1 - cy0;
	s32 *pz = d_pzbuffer +(d_zwidth * cy0) + cx0;
	u8 *pdest = d_viewbuffer + d_scantable[cy0] + cx0;
#define MIX_PIXEL color_mix_lut[pdest[x]][pd->color][part_alpha]
// the row is cy1 - count in both shapes
#define DITHER d_bayer[(cy1 - count) & 7][(cx0 + x) & 7] < level
	if(r_particlestyle.value == 1 && part_alpha == FOG_LUT_LEVELS){// circle
		PART_CIRCLE(pz[x] <= izi, pd->color);
	}else if(part_alpha == FOG_LUT_LEVELS){ // square
		PART_SQUARE(pz[x] <= izi, pd->color);
	}else if(r_particlestyle.value == 1 && r_alphastyle.value){
		// circle(dither alpha)
		PART_CIRCLE(pz[x] <= izi && DITHER, pd->color);
	}else if(r_alphastyle.value){ // square(dither alpha)
		PART_SQUARE(pz[x] <= izi && DITHER, pd->color);
	}else if(r_particlestyle.value == 1){ // circle(mix alpha)
		PART_CIRCLE(pz[x] <= izi, MIX_PIXEL);
	}else{ // square(mix alpha)
		PART_SQUARE(pz[x] <= izi, MIX_PIXEL);
	}
#undef MIX_PIXEL
#undef DITHER
}

void D_DrawParticle(particle_t *pparticle)
{
	partdraw_t pd;
	if(!D_ProjectParticle(pparticle, &pd))
		return;
	if(!r_threaded){
		D_FillParticle(&pd);
		return;
	}
	if(d_nparts == d_partcap){ // kept for D_EndParticles
		d_partcap = d_partcap ? d_partcap * 2 : 1024;
		d_parts = realloc(d_parts, d_partcap * sizeof(*d_parts));
		if(!d_parts)
			Sys_Error("D_DrawParticle: out of memory");
	}
	d_parts[d_nparts++] = pd;
}

typedef struct { s32 count; partdraw_t *parts; } partband_t;

static void D_ParticleCmd(void *arg)
{
	const partband_t *pb = arg;
	for(s32 i = 0; i < pb->count; i++)
		D_FillParticle(&pb->parts[i]);
}

// The particles of each band in the order they came, for the workers
void D_EndParticles()
{
	if(!d_nparts)
		return;
	partband_t *bands = R_FrameAlloc(r_nbands * sizeof(*bands));
	memset(bands, 0, r_nbands * sizeof(*bands));
	for(s32 i = 0; i < d_nparts; i++)
		for(s32 b = d_parts[i].cy0 / r_bandrows;
		    b <= (d_parts[i].cy1 - 1) / r_bandrows; b++)
			bands[b].count++;
	for(s32 b = 0; b < r_nbands; b++){
		if(!bands[b].count)
			continue;
		bands[b].parts = R_FrameAlloc(bands[b].count
				* sizeof(*bands[b].parts));
		bands[b].count = 0;
	}
	for(s32 i = 0; i < d_nparts; i++)
		for(s32 b = d_parts[i].cy0 / r_bandrows;
		    b <= (d_parts[i].cy1 - 1) / r_bandrows; b++)
			bands[b].parts[bands[b].count++] = d_parts[i];
	for(s32 b = 0; b < r_nbands; b++){
		if(!bands[b].count)
			continue;
		R_BandAdd(b, D_ParticleCmd, &bands[b]);
		R_BandPublish(b);
	}
	d_nparts = 0;
}
