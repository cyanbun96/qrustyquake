// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

static s32 clip_current;
static vec5_t clip_verts[2][MAXWORKINGVERTS];
static s32 sprite_width, sprite_height;

void R_RotateSprite(f32 beamlength)
{
	vec3_t vec;

	if (beamlength == 0.0)
		return;

	VectorScale(r_spritedesc.vpn, -beamlength, vec);
	VectorAdd(r_entorigin, vec, r_entorigin);
	VectorSubtract(modelorg, vec, modelorg);
}


s32 R_ClipSpriteFace(s32 nump, clipplane_t *pclipplane) // Clips the winding at
{ // clip_verts[clip_current] and changes clip_current, throws out the back side
	f32 clipdist = pclipplane->dist;
	f32 *pclipnormal = pclipplane->normal;
	f32 *in, *outstep;
	if (clip_current) { // calc dists
		in = clip_verts[1][0];
		outstep = clip_verts[0][0];
		clip_current = 0;
	} else {
		in = clip_verts[0][0];
		outstep = clip_verts[1][0];
		clip_current = 1;
	}
	f32 *instep = in;
	f32 dists[MAXWORKINGVERTS + 1];
	for (s32 i = 0; i < nump; i++, instep += sizeof(vec5_t) / sizeof(f32))
		dists[i] = DotProduct(instep, pclipnormal) - clipdist;
	dists[nump] = dists[0]; // handle wraparound case
	Q_memcpy(instep, in, sizeof(vec5_t));
	instep = in; // clip the winding
	s32 outcount = 0;
	for (s32 i = 0; i < nump; i++, instep += sizeof(vec5_t)/sizeof(f32)) {
		if (dists[i] >= 0) {
			Q_memcpy(outstep, instep, sizeof(vec5_t));
			outstep += sizeof(vec5_t) / sizeof(f32);
			outcount++;
		}
		if (dists[i] == 0 || dists[i + 1] == 0
			|| (dists[i] > 0) == (dists[i + 1] > 0))
			continue;
		// split it into a new vertex
		f32 frac = dists[i] / (dists[i] - dists[i + 1]);
		f32 *vert2 = instep + sizeof(vec5_t) / sizeof(f32);
		outstep[0] = instep[0] + frac * (vert2[0] - instep[0]);
		outstep[1] = instep[1] + frac * (vert2[1] - instep[1]);
		outstep[2] = instep[2] + frac * (vert2[2] - instep[2]);
		outstep[3] = instep[3] + frac * (vert2[3] - instep[3]);
		outstep[4] = instep[4] + frac * (vert2[4] - instep[4]);
		outstep += sizeof(vec5_t) / sizeof(f32);
		outcount++;
	}
	return outcount;
}

void R_SetupAndDrawSprite()
{
	if (DotProduct(r_spritedesc.vpn, modelorg) >= 0) // backface cull
		return;
	// build the sprite poster in worldspace
	vec3_t left, up, right, down, transformed, local;
	VectorScale(r_spritedesc.vright,r_spritedesc.pspriteframe->right,right);
	VectorScale(r_spritedesc.vup, r_spritedesc.pspriteframe->up, up);
	VectorScale(r_spritedesc.vright, r_spritedesc.pspriteframe->left, left);
	VectorScale(r_spritedesc.vup, r_spritedesc.pspriteframe->down, down);
	vec5_t *pverts = clip_verts[0];
	pverts[0][0] = r_entorigin[0] + up[0] + left[0];
	pverts[0][1] = r_entorigin[1] + up[1] + left[1];
	pverts[0][2] = r_entorigin[2] + up[2] + left[2];
	pverts[0][3] = 0;
	pverts[0][4] = 0;
	pverts[1][0] = r_entorigin[0] + up[0] + right[0];
	pverts[1][1] = r_entorigin[1] + up[1] + right[1];
	pverts[1][2] = r_entorigin[2] + up[2] + right[2];
	pverts[1][3] = sprite_width;
	pverts[1][4] = 0;
	pverts[2][0] = r_entorigin[0] + down[0] + right[0];
	pverts[2][1] = r_entorigin[1] + down[1] + right[1];
	pverts[2][2] = r_entorigin[2] + down[2] + right[2];
	pverts[2][3] = sprite_width;
	pverts[2][4] = sprite_height;
	pverts[3][0] = r_entorigin[0] + down[0] + left[0];
	pverts[3][1] = r_entorigin[1] + down[1] + left[1];
	pverts[3][2] = r_entorigin[2] + down[2] + left[2];
	pverts[3][3] = 0;
	pverts[3][4] = sprite_height;
	s32 nump = 4; // clip to the frustum in worldspace
	clip_current = 0;
	for (s32 i = 0; i < 4; i++) {
		nump = R_ClipSpriteFace(nump, &view_clipplanes[i]);
		if (nump < 3)
			return;
		if (nump >= MAXWORKINGVERTS)
			Sys_Error("R_SetupAndDrawSprite: too many points");
	}
	// transform vertices into viewspace and project
	f32 *pv = &clip_verts[clip_current][0][0];
	r_spritedesc.nearzi = -999999;
	emitpoint_t outverts[MAXWORKINGVERTS + 1];
	for (s32 i = 0; i < nump; i++) {
		VectorSubtract(pv, r_origin, local);
		TransformVector(local, transformed);
		if (transformed[2] < NEAR_CLIP)
			transformed[2] = NEAR_CLIP;
		emitpoint_t *pout = &outverts[i];
		pout->zi = 1.0 / transformed[2];
		if (pout->zi > r_spritedesc.nearzi)
			r_spritedesc.nearzi = pout->zi;
		pout->s = pv[3];
		pout->t = pv[4];
		f32 scale = xscale * pout->zi;
		pout->u = (xcenter + scale * transformed[0]);
		scale = yscale * pout->zi;
		pout->v = (ycenter - scale * transformed[1]);
		pv += sizeof(vec5_t) / sizeof(*pv);
	}
	r_spritedesc.nump = nump;
	r_spritedesc.pverts = outverts;
	D_DrawSprite();
}

mspriteframe_t *R_GetSpriteframe(msprite_t *psprite)
{
	s32 frame = currententity->frame;
	if ((frame >= psprite->numframes) || (frame < 0)) {
		Con_DPrintf("R_DrawSprite: no such frame %d\n", frame);
		frame = 0;
	}
	mspritegroup_t *pspritegroup;
	mspriteframe_t *pspriteframe;
	if (psprite->frames[frame].type == SPR_SINGLE)
		pspriteframe = psprite->frames[frame].frameptr;
	else {
		pspritegroup = (mspritegroup_t*)psprite->frames[frame].frameptr;
		f32 *pintervals = pspritegroup->intervals;
		s32 numframes = pspritegroup->numframes;
		f32 fullinterval = pintervals[numframes - 1];
		f32 time = cl.time + currententity->syncbase;
		// when loading in Mod_LoadSpriteGroup, we guaranteed all interval values
		// are positive, so we don't have to worry about division by 0
		f32 targettime =
			time - ((s32)(time / fullinterval)) * fullinterval;
		s32 i = 0;
		for (; i < (numframes - 1); i++)
			if (pintervals[i] > targettime)
				break;
		pspriteframe = pspritegroup->frames[i];
	}
	return pspriteframe;
}

void R_DrawSprite()
{
	msprite_t *psprite = currententity->model->cache.data;
	r_spritedesc.pspriteframe = R_GetSpriteframe(psprite);
	sprite_width = r_spritedesc.pspriteframe->width;
	sprite_height = r_spritedesc.pspriteframe->height;
	// TODO: make this caller-selectable
	if (psprite->type == SPR_FACING_UPRIGHT) {
		// generate the sprite's axes, with vup straight up in worldspace, and
		// r_spritedesc.vright perpendicular to modelorg.
		// This will not work if the view direction is very close to straight up or
		// down, because the cross product will be between two nearly parallel
		// vectors and starts to approach an undefined state, so we don't draw if
		// the two vectors are less than 1 degree apart
		vec3_t tvec;
		tvec[0] = -modelorg[0];
		tvec[1] = -modelorg[1];
		tvec[2] = -modelorg[2];
		VectorNormalize(tvec);
		f32 dot = tvec[2]; // same as DotProduct (tvec, r_spritedesc.vup) because
		// r_spritedesc.vup is 0, 0, 1
		if ((dot > 0.999848) || (dot < -0.999848)) // cos(1 degree) = 0.999848
			return;
		r_spritedesc.vup[0] = 0;
		r_spritedesc.vup[1] = 0;
		r_spritedesc.vup[2] = 1;
		r_spritedesc.vright[0] = tvec[1];
		// CrossProduct(r_spritedesc.vup, -modelorg,
		r_spritedesc.vright[1] = -tvec[0];
		//              r_spritedesc.vright)
		r_spritedesc.vright[2] = 0;
		VectorNormalize(r_spritedesc.vright);
		r_spritedesc.vpn[0] = -r_spritedesc.vright[1];
		r_spritedesc.vpn[1] = r_spritedesc.vright[0];
		r_spritedesc.vpn[2] = 0;
		// CrossProduct (r_spritedesc.vright, r_spritedesc.vup,
		//  r_spritedesc.vpn)
	} else if (psprite->type == SPR_VP_PARALLEL) {
		// generate the sprite's axes, completely parallel to the viewplane. There
		// are no problem situations, because the sprite is always in the same
		// position relative to the viewer
		for (s32 i = 0; i < 3; i++) {
			r_spritedesc.vup[i] = vup[i];
			r_spritedesc.vright[i] = vright[i];
			r_spritedesc.vpn[i] = vpn[i];
		}
	} else if (psprite->type == SPR_VP_PARALLEL_UPRIGHT) {
		f32 dot = vpn[2];
		if ((dot > 0.999848) || (dot < -0.999848))
			return;
		r_spritedesc.vup[0] = 0;
		r_spritedesc.vup[1] = 0;
		r_spritedesc.vup[2] = 1;
		r_spritedesc.vright[0] = vpn[1];
		r_spritedesc.vright[1] = -vpn[0];
		r_spritedesc.vright[2] = 0;
		VectorNormalize(r_spritedesc.vright);
		r_spritedesc.vpn[0] = -r_spritedesc.vright[1];
		r_spritedesc.vpn[1] = r_spritedesc.vright[0];
		r_spritedesc.vpn[2] = 0;
	} else if (psprite->type == SPR_ORIENTED) {
		// generate the sprite's axes, according to the sprite's world orientation
		AngleVectors(currententity->angles, r_spritedesc.vpn,
			     r_spritedesc.vright, r_spritedesc.vup);
	} else if (psprite->type == SPR_VP_PARALLEL_ORIENTED) {
		// generate the sprite's axes, parallel to the viewplane, but rotated in
		// that plane around the center according to the sprite entity's roll
		// angle. So vpn stays the same, but vright and vup rotate
		f32 angle = currententity->angles[ROLL] * (M_PI * 2 / 360);
		f32 sr = sin(angle);
		f32 cr = cos(angle);
		for (s32 i = 0; i < 3; i++) {
			r_spritedesc.vpn[i] = vpn[i];
			r_spritedesc.vright[i] = vright[i] * cr + vup[i] * sr;
			r_spritedesc.vup[i] = vright[i] * -sr + vup[i] * cr;
		}
	} else
		Sys_Error("R_DrawSprite: Bad sprite type %d", psprite->type);
	R_RotateSprite(psprite->beamlength);
	R_SetupAndDrawSprite();
}
