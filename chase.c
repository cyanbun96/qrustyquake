// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// chase.c -- chase camera code
#include "quakedef.h"

static vec3_t chase_dest;

void Chase_Init()
{
	Cvar_RegisterVariable((struct cvar_s *)&chase_back);
	Cvar_RegisterVariable((struct cvar_s *)&chase_up);
	Cvar_RegisterVariable((struct cvar_s *)&chase_right);
	Cvar_RegisterVariable((struct cvar_s *)&chase_active);
}

void TraceLine(vec3_t start, vec3_t end, vec3_t impact)
{
	trace_t tr;
	memset(&tr, 0, sizeof(tr));
	SV_RecursiveHullCheck(cl.worldmodel->hulls, 0, 0, 1, start, end, &tr);
	VectorCopy(tr.endpos, impact);
}

void Chase_Update()
{
	vec3_t forward, up, right, dest, stop;
	// if can't see player, reset
	AngleVectors(cl.viewangles, forward, right, up);
	for(s32 i = 0; i < 3; i++) // calc exact destination
		chase_dest[i] = r_refdef.vieworg[i]
		    - forward[i] * chase_back.value
		    - right[i] * chase_right.value;
	chase_dest[2] = r_refdef.vieworg[2] + chase_up.value;
	// find the spot the player is looking at
	VectorMA(r_refdef.vieworg, 4096, forward, dest);
	TraceLine(r_refdef.vieworg, dest, stop);
	// calculate pitch to look at the same spot from camera
	VectorSubtract(stop, r_refdef.vieworg, stop);
	f32 dist = DotProduct(stop, forward);
	if(dist < 1) dist = 1;
	r_refdef.viewangles[PITCH] = -atan(stop[2] / dist) / M_PI * 180;
	VectorCopy(chase_dest, r_refdef.vieworg); // move towards destination
}
