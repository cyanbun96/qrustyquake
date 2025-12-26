// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
// sv_move.c -- monster movement
#include "quakedef.h"

static s32 c_yes, c_no;

bool SV_CheckBottom(edict_t *ent) // Returns 0 if any part of the bottom of the
{ // entity is off an edge that is not a staircase.
	vec3_t mins, maxs, start, stop;
	f32 mid, bottom;
	VectorAdd(ent->v.origin, ent->v.mins, mins);
	VectorAdd(ent->v.origin, ent->v.maxs, maxs);
	// if all of the points under the corners are solid world, don't bother
	// with tougher checks. the corners must be within 16 of the midpoint
	start[2] = mins[2] - 1;
	for(s32 x = 0; x <= 1; x++)
	for(s32 y = 0; y <= 1; y++) {
		start[0] = x ? maxs[0] : mins[0];
		start[1] = y ? maxs[1] : mins[1];
		if(SV_PointContents(start) != CONTENTS_SOLID) goto realcheck;
	}
	c_yes++;
	return 1; // we got out easy
realcheck:
	c_no++;
	start[2] = mins[2]; // check it for real...
	start[0] = stop[0] = (mins[0] + maxs[0])*0.5; // the midpoint must be
	start[1] = stop[1] = (mins[1] + maxs[1])*0.5; // within 16 of the bottom
	stop[2] = start[2] - 2*STEPSIZE;
	trace_t trace = SV_Move(start, vec3_origin, vec3_origin, stop, 1, ent);
	if(trace.fraction == 1.0) return 0;
	mid = bottom = trace.endpos[2];
	for(s32 x = 0; x <= 1; x++) // corners must be within 16 of the midpoint
	for(s32 y = 0; y <= 1; y++) {
		start[0] = stop[0] = x ? maxs[0] : mins[0];
		start[1] = stop[1] = y ? maxs[1] : mins[1];
		trace = SV_Move(start, vec3_origin, vec3_origin, stop, 1, ent);
		if(trace.fraction != 1.0 && trace.endpos[2] > bottom)
			bottom = trace.endpos[2];
		if(trace.fraction == 1.0 || mid - trace.endpos[2] > STEPSIZE)
			return 0;
	}
	c_yes++;
	return 1;
}

// Called by monster program code.
// The move will be adjusted for slopes and stairs, but if the move isn't
// possible, no move is done, 0 is returned, and
// pr_global_struct->trace_normal is set to the normal of the blocking wall
bool SV_movestep(edict_t *ent, vec3_t move, bool relink)
{
	vec3_t oldorg, neworg, end;
	VectorCopy(ent->v.origin, oldorg);
	VectorAdd(ent->v.origin, move, neworg);
	if((s32)ent->v.flags & (FL_SWIM|FL_FLY)){//flying monsters don't step up
		// try one move with vertical motion, then one without
		for(s32 i = 0; i < 2; i++) {
			VectorAdd(ent->v.origin, move, neworg);
			edict_t *enemy = PROG_TO_EDICT(ent->v.enemy);
			if(i == 0 && enemy != sv.edicts) {
				f32 dz = ent->v.origin[2] - 
				    PROG_TO_EDICT(ent->v.enemy)->v.origin[2];
				if(dz > 40) neworg[2] -= 8;
				if(dz < 30) neworg[2] += 8;
			}
			trace_t trace = SV_Move(ent->v.origin, ent->v.mins,
					ent->v.maxs, neworg, 0, ent);
			if(trace.fraction == 1) {
				if(((s32)ent->v.flags & FL_SWIM) &&
				 SV_PointContents(trace.endpos)==CONTENTS_EMPTY)
					return 0; // swim monster left water
				VectorCopy(trace.endpos, ent->v.origin);
				if(relink) SV_LinkEdict(ent, 1);
				return 1;
			}
			if(enemy == sv.edicts) break;
		}
		return 0;
	}
	// push down from a step height above the wished position
	neworg[2] += STEPSIZE;
	VectorCopy(neworg, end);
	end[2] -= STEPSIZE*2;
	trace_t trace = SV_Move(neworg, ent->v.mins, ent->v.maxs, end, 0, ent);
	if(trace.allsolid) return 0;
	if(trace.startsolid) {
		neworg[2] -= STEPSIZE;
		trace = SV_Move(neworg, ent->v.mins, ent->v.maxs, end, 0, ent);
		if(trace.allsolid || trace.startsolid) return 0;
	}
	if(trace.fraction == 1) {
		// if monster had the ground pulled out, go ahead and fall
		if( (s32)ent->v.flags & FL_PARTIALGROUND ) {
			VectorAdd(ent->v.origin, move, ent->v.origin);
			if(relink) SV_LinkEdict(ent, 1);
			ent->v.flags = (s32)ent->v.flags & ~FL_ONGROUND;
			return 1;
		}
		return 0; // walked off an edge
	}
	// check point traces down for dangling corners
	VectorCopy(trace.endpos, ent->v.origin);
	if(!SV_CheckBottom(ent)) {
		if((s32)ent->v.flags & FL_PARTIALGROUND){
			// entity had floor mostly pulled out from underneath
			// it and is trying to correct
			if(relink) SV_LinkEdict(ent, 1);
			return 1;
		}
		VectorCopy(oldorg, ent->v.origin);
		return 0;
	}
	if((s32)ent->v.flags & FL_PARTIALGROUND)
		ent->v.flags = (s32)ent->v.flags & ~FL_PARTIALGROUND;
	ent->v.groundentity = EDICT_TO_PROG(trace.ent);
	if(relink) SV_LinkEdict(ent, 1); // the move is ok
	return 1;
}

bool SV_StepDirection(edict_t *ent, f32 yaw, f32 dist)
{ // Turns to the movement direction and walks the current distance if facing it
	vec3_t move, oldorigin;
	ent->v.ideal_yaw = yaw;
	PF_changeyaw();
	yaw = yaw*M_PI*2 / 360;
	move[0] = cos(yaw)*dist;
	move[1] = sin(yaw)*dist;
	move[2] = 0;
	VectorCopy(ent->v.origin, oldorigin);
	if(SV_movestep(ent, move, 0)) {
		f32 delta = ent->v.angles[YAW] - ent->v.ideal_yaw;
		if(delta > 45 && delta < 315) {
			VectorCopy(oldorigin, ent->v.origin);
		} // not turned far enough, so don't take the step
		SV_LinkEdict(ent, 1);
		return 1;
	}
	SV_LinkEdict(ent, 1);
	return 0;
}

void SV_FixCheckBottom(edict_t *ent)
{ ent->v.flags = (s32)ent->v.flags | FL_PARTIALGROUND; }

void SV_NewChaseDir(edict_t *actor, edict_t *enemy, f32 dist)
{
	f32 d[3];
	f32 tdir;
	f32 olddir = anglemod( (s32)(actor->v.ideal_yaw/45)*45 );
	f32 turnaround = anglemod(olddir - 180);
	f32 deltax = enemy->v.origin[0] - actor->v.origin[0];
	f32 deltay = enemy->v.origin[1] - actor->v.origin[1];
	if(deltax>10)       d[1]= 0;
	else if(deltax<-10) d[1]= 180;
	else                d[1]= DI_NODIR;
	if(deltay<-10)      d[2]= 270;
	else if(deltay>10)  d[2]= 90;
	else                d[2]= DI_NODIR;
	if(d[1] != DI_NODIR && d[2] != DI_NODIR) { // try direct route
		if(d[1] == 0) tdir = d[2] == 90 ? 45 : 315;
		else          tdir = d[2] == 90 ? 135 : 215;
		if(tdir!=turnaround&&SV_StepDirection(actor,tdir,dist)) return;
	}
	if(((rand()&3)&1)||abs((s32)deltay)>abs((s32)deltax)){
		tdir=d[1]; // try other directions
		d[1]=d[2];
		d[2]=tdir;
	}
	if(d[1]!=DI_NODIR&&d[1]!=turnaround&&SV_StepDirection(actor,d[1],dist))
		return;
	if(d[2]!=DI_NODIR&&d[2]!=turnaround&&SV_StepDirection(actor,d[2],dist))
		return;
	// there is no direct path to the player, so pick another direction
	if(olddir!=DI_NODIR && SV_StepDirection(actor, olddir, dist)) return;
	if(rand()&1){ // randomly determine direction of search
		for(tdir = 0; tdir<=315; tdir += 45)
			if(tdir!=turnaround&&SV_StepDirection(actor,tdir,dist))
				return;
	} else {
		for(tdir = 315; tdir >= 0; tdir -= 45)
			if(tdir!=turnaround&&SV_StepDirection(actor,tdir,dist))
				return;
	}
	if(turnaround!=DI_NODIR&&SV_StepDirection(actor,turnaround,dist))
		return;
	actor->v.ideal_yaw = olddir; // can't move
		// if a bridge was pulled out from underneath a monster, it may
		// not have a valid standing position at all
	if(!SV_CheckBottom(actor)) SV_FixCheckBottom(actor);
}

bool SV_CloseEnough(edict_t *ent, edict_t *goal, f32 dist)
{
	for(s32 i = 0; i < 3; i++) {
		if(goal->v.absmin[i] > ent->v.absmax[i] + dist) return 0;
		if(goal->v.absmax[i] < ent->v.absmin[i] - dist) return 0;
	}
	return 1;
}

void SV_MoveToGoal()
{
	edict_t *ent = PROG_TO_EDICT(pr_global_struct->self);
	edict_t *goal = PROG_TO_EDICT(ent->v.goalentity);
	f32 dst = G_FLOAT(OFS_PARM0);
	if(!((s32)ent->v.flags & (FL_ONGROUND|FL_FLY|FL_SWIM))) {
		G_FLOAT(OFS_RETURN) = 0;
		return;
	}
	// if the next step hits the enemy, return immediately
	if(PROG_TO_EDICT(ent->v.enemy)!=sv.edicts&&SV_CloseEnough(ent,goal,dst))
		return;
	// bump around...
	if((rand()&3)==1 || !SV_StepDirection(ent, ent->v.ideal_yaw, dst))
		SV_NewChaseDir(ent, goal, dst);
}
