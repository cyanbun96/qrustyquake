// d_part.c: software driver module for drawing particles

#include "quakedef.h"

void D_DrawParticle(particle_t *pparticle)
{
	vec3_t local, transformed;
	VectorSubtract(pparticle->org, r_origin, local); // transform point
	transformed[0] = DotProduct(local, r_pright);
	transformed[1] = DotProduct(local, r_pup);
	transformed[2] = DotProduct(local, r_ppn);
	if (transformed[2] < PARTICLE_Z_CLIP)
		return;
	f32 zi = 1.0 / transformed[2]; // project the point
	s32 u = (s32)(xcenter + zi * transformed[0] + 0.5); // FIXME: preadjust xcenter and ycenter
	s32 v = (s32)(ycenter - zi * transformed[1] + 0.5);
	
	// Clipping Check 1: Top/Left bounds
	if ((v > d_vrectbottom_particle) || (u > d_vrectright_particle)
		|| (v < d_vrecty) || (u < d_vrectx)) 
		return;
	
	s16 *pz = d_pzbuffer + (d_zwidth * v) + u;
	u8 *pdest = d_viewbuffer + d_scantable[v] + u;
	s32 izi = (s32)(zi * 0x8000);
	s32 pix;
	
	if (r_particlescale.value) {
		pix = r_particlescale.value;
	} else {
		pix = izi >> d_pix_shift;
		if (pix < d_pix_min) pix = d_pix_min;
		else if (pix > d_pix_max) pix = d_pix_max;
	}

	// FIX: Use u/v directly. Do not recalculate from pointers.
	// We must clip against the absolute bottom/right edges of the viewport.
	s32 pix_y = v;
	s32 pix_x = u;
	
	// Calculate absolute bounds
	s32 bound_bottom = scr_vrect.y + scr_vrect.height;
	s32 bound_right  = scr_vrect.x + scr_vrect.width;

	// Bottom Clipping
	if (pix_y + pix > bound_bottom) {
		pix = bound_bottom - pix_y;
		if (pix <= 0) return; // Entirely clipped
	}

	// Right Clipping
	s32 pix_w = pix;
	if (pix_x + pix_w > bound_right) {
		pix_w = bound_right - pix_x;
		if (pix_w <= 0) return; // Entirely clipped
	}

	// Safety clamp for buffer width (handle odd edge cases)
	s32 max_w = screenwidth - pix_x;
	if (pix_w > max_w) pix_w = max_w;

	// Draw loop
	for (s32 count=pix; count; count--, pz+=d_zwidth, pdest+=screenwidth) {
		for (s32 i = 0; i < pix_w; i++) {
			if (pz[i] <= izi) {
				pz[i] = izi;
				pdest[i] = pparticle->color;
			}
		}
	}
}
