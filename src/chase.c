// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// Copyright (C) 2002-2009 John Fitzgibbons and others
// Copyright (C) 2010-2014 QuakeSpasm developers
// chase.c -- chase camera code
#include "quakedef.h"

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
{ // TODO: stay at least 8 units away from all walls in this leaf
	vec3_t forward, up, right;
	vec3_t ideal, crosshair, temp;
	AngleVectors(cl.viewangles, forward, right, up);
	// calc ideal camera location before checking for walls
	for(s32 i = 0; i < 3; i++)
		ideal[i] = cl.viewent.origin[i]
			- forward[i]*chase_back.value
			+ right[i]*chase_right.value;
	ideal[2] = cl.viewent.origin[2] + chase_up.value;
	// make sure camera is not in or behind a wall
	TraceLine(r_refdef.vieworg, ideal, temp);
	if(VectorLength(temp)) VectorCopy(temp, ideal);
	// place camera
	VectorCopy(ideal, r_refdef.vieworg);
	// find the spot the player is looking at
	VectorMA(cl.viewent.origin, 1<<20, forward, temp);
	TraceLine(cl.viewent.origin, temp, crosshair);
	// calculate camera angles to look at the same spot
	VectorSubtract(crosshair, r_refdef.vieworg, temp);
	VectorAngles(temp, r_refdef.viewangles);
	if(r_refdef.viewangles[PITCH]==90 || r_refdef.viewangles[PITCH]==-90)
		r_refdef.viewangles[YAW] = cl.viewangles[YAW];
}
