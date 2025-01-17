// Copyright (C) 1996-2001 Id Software, Inc.
// Copyright (C) 2002-2009 John Fitzgibbons and others
// GPLv3 See LICENSE for details.

#include "quakedef.h"

#define MOVE_EPSILON 0.01
#define STOP_EPSILON 0.1
#define MAX_CLIP_PLANES 5
#define STEPSIZE 18

/*
pushmove objects do not obey gravity, and do not interact with each other or
trigger fields, but block normal movement and push normal objects when they move
onground is set for toss objects when they come to a complete rest. it is set
for steping or walking objects

doors, plats, etc are SOLID_BSP, and MOVETYPE_PUSH
bonus items are SOLID_TRIGGER touch, and MOVETYPE_TOSS
corpses are SOLID_NOT and MOVETYPE_TOSS
crates are SOLID_BBOX and MOVETYPE_TOSS
walking monsters are SOLID_SLIDEBOX and MOVETYPE_STEP
flying/floating monsters are SOLID_SLIDEBOX and MOVETYPE_FLY

solid_edge items only clip against bsp models.
*/

cvar_t sv_friction = { "sv_friction", "4", false, true, 0, NULL };
cvar_t sv_stopspeed = { "sv_stopspeed", "100", false, false, 0, NULL };
cvar_t sv_gravity = { "sv_gravity", "800", false, true, 0, NULL };
cvar_t sv_maxvelocity = { "sv_maxvelocity", "2000", false, false, 0, NULL };
cvar_t sv_nostep = { "sv_nostep", "0", false, false, 0, NULL };

void SV_Physics_Toss(edict_t * ent);

void SV_CheckVelocity(edict_t *ent)
{
	for (int i = 0; i < 3; i++) { // bound velocity
		if (IS_NAN(ent->v.velocity[i])) {
			Con_Printf("Got a NaN velocity on %s\n",
					PR_GetString(ent->v.classname));
			ent->v.velocity[i] = 0;
		}
		if (IS_NAN(ent->v.origin[i])) {
			Con_Printf("Got a NaN origin on %s\n",
					PR_GetString(ent->v.classname));
			ent->v.origin[i] = 0;
		}
		if (ent->v.velocity[i] > sv_maxvelocity.value)
			ent->v.velocity[i] = sv_maxvelocity.value;
		else if (ent->v.velocity[i] < -sv_maxvelocity.value)
			ent->v.velocity[i] = -sv_maxvelocity.value;
	}
}

// Runs thinking code if time. There is some play in the exact time the think
// function will be called, because it is called before any movement is done
// in a frame. Not used for pushmove objects, because they must be exact.
// Returns false if the entity removed itself.
qboolean SV_RunThink(edict_t *ent)
{
	float thinktime = ent->v.nextthink;
	if (thinktime <= 0 || thinktime > sv.time + host_frametime)
		return true;
	if (thinktime < sv.time)
		thinktime = sv.time; // don't let things stay in the past.
				     // it is possible to start that way
				     // by a trigger with a local time.
	float oldframe = ent->v.frame; //johnfitz
	ent->v.nextthink = 0;
	pr_global_struct->time = thinktime;
	pr_global_struct->self = EDICT_TO_PROG(ent);
	pr_global_struct->other = EDICT_TO_PROG(sv.edicts);
	PR_ExecuteProgram(ent->v.think);
	//johnfitz -- PROTOCOL_FITZQUAKE
	//capture interval to nextthink here and send it to client for better
	//lerp timing, but only if interval is not 0.1 (which client assumes)
	ent->sendinterval = false;
	if (!ent->free && ent->v.nextthink && (ent->v.movetype == MOVETYPE_STEP
				|| ent->v.frame != oldframe)) {
		int i = Q_rint((ent->v.nextthink - thinktime) * 255); //johnfitz
		if (i >= 0 && i < 256 && i != 25 && i != 26) //25 and 26 are close enough to 0.1 to not send
			ent->sendinterval = true;
	} //johnfitz
	return !ent->free;
}

void SV_Impact(edict_t *e1, edict_t *e2)
{ // Two entities have touched, so run their touch functions
	int old_self, old_other;
	old_self = pr_global_struct->self;
	old_other = pr_global_struct->other;
	pr_global_struct->time = sv.time;
	if (e1->v.touch && e1->v.solid != SOLID_NOT) {
		pr_global_struct->self = EDICT_TO_PROG(e1);
		pr_global_struct->other = EDICT_TO_PROG(e2);
		PR_ExecuteProgram(e1->v.touch);
	}
	if (e2->v.touch && e2->v.solid != SOLID_NOT) {
		pr_global_struct->self = EDICT_TO_PROG(e2);
		pr_global_struct->other = EDICT_TO_PROG(e1);
		PR_ExecuteProgram(e2->v.touch);
	}
	pr_global_struct->self = old_self;
	pr_global_struct->other = old_other;
}

// Slide off of the impacting object
// returns the blocked flags (1 = floor, 2 = step / wall)
int ClipVelocity(vec3_t in, vec3_t normal, vec3_t out, float overbounce)
{
	int blocked = 0;
	if (normal[2] > 0)
		blocked |= 1; // floor
	if (!normal[2])
		blocked |= 2; // step
	float backoff = DotProduct(in, normal) * overbounce;
	for (int i = 0; i < 3; i++) {
		float change = normal[i] * backoff;
		out[i] = in[i] - change;
		if (out[i] > -STOP_EPSILON && out[i] < STOP_EPSILON)
			out[i] = 0;
	}
	return blocked;
}

// The basic solid body movement clip that slides along multiple planes
// Returns the clipflags if the velocity was modified (hit something solid)
// 1 = floor
// 2 = wall / step
// 4 = dead stop
// If steptrace is not NULL, the trace of any vertical wall hit will be stored
int SV_FlyMove(edict_t *ent, float time, trace_t *steptrace)
{
	vec3_t planes[MAX_CLIP_PLANES];
	int numbumps = 4;
	int blocked = 0;
	vec3_t primal_velocity, original_velocity, new_velocity;
	VectorCopy(ent->v.velocity, original_velocity);
	VectorCopy(ent->v.velocity, primal_velocity);
	int numplanes = 0;
	float time_left = time;
	for (int bumpcount = 0; bumpcount < numbumps; bumpcount++) {
		if (!ent->v.velocity[0] && !ent->v.velocity[1]
				&& !ent->v.velocity[2])
			break;
		vec3_t end;
		for (int i = 0; i < 3; i++)
			end[i] = ent->v.origin[i]+time_left*ent->v.velocity[i];
		trace_t trace = SV_Move(ent->v.origin, ent->v.mins, ent->v.maxs,
				end, false, ent);
		if (trace.allsolid) { // entity is trapped in another solid
			VectorCopy(vec3_origin, ent->v.velocity);
			return 3;
		}
		if (trace.fraction > 0) { // actually covered some distance
			VectorCopy(trace.endpos, ent->v.origin);
			VectorCopy(ent->v.velocity, original_velocity);
			numplanes = 0;
		}
		if (trace.fraction == 1)
			break; // moved the entire distance
		if (!trace.ent)
			Sys_Error("SV_FlyMove: !trace.ent");
		if (trace.plane.normal[2] > 0.7) {
			blocked |= 1; // floor
			if (trace.ent->v.solid == SOLID_BSP) {
				ent->v.flags = (int)ent->v.flags | FL_ONGROUND;
				ent->v.groundentity = EDICT_TO_PROG(trace.ent);
			}
		}
		if (!trace.plane.normal[2]) {
			blocked |= 2; // step
			if (steptrace)
				*steptrace = trace; // save for player extrafriction
		}
		// run the impact function
		SV_Impact(ent, trace.ent);
		if (ent->free)
			break; // removed by the impact function
		time_left -= time_left * trace.fraction;
		// cliped to another plane
		if (numplanes >= MAX_CLIP_PLANES) { // this shouldn't really happen
			VectorCopy(vec3_origin, ent->v.velocity);
			return 3;
		}
		VectorCopy(trace.plane.normal, planes[numplanes]);
		numplanes++;
		// modify original_velocity so it parallels all of the clip planes
		int i = 0;
		for (; i < numplanes; i++) {
			ClipVelocity(original_velocity, planes[i], new_velocity,
					1);
			int j = 0;
			for (; j < numplanes; j++)
				if (j != i) {
					if (DotProduct(new_velocity, planes[j])
							< 0)
						break; // not ok
				}
			if (j == numplanes)
				break;
		}
		if (i != numplanes) { // go along this plane
			VectorCopy(new_velocity, ent->v.velocity);
		} else { // go along the crease
			if (numplanes != 2) {
				VectorCopy(vec3_origin, ent->v.velocity);
				return 7;
			}
			vec3_t dir;
			CrossProduct(planes[0], planes[1], dir);
			float d = DotProduct(dir, ent->v.velocity);
			VectorScale(dir, d, ent->v.velocity);
		}
		// if original velocity is against the original velocity, stop
		// dead to avoid tiny occilations in sloping corners
		if (DotProduct(ent->v.velocity, primal_velocity) <= 0) {
			VectorCopy(vec3_origin, ent->v.velocity);
			return blocked;
		}
	}
	return blocked;
}

void SV_AddGravity(edict_t *ent)
{
	float ent_gravity;
	eval_t *val;
	val = GetEdictFieldValue(ent, "gravity");
	if (val && val->_float)
		ent_gravity = val->_float;
	else
		ent_gravity = 1.0;
	ent->v.velocity[2] -= ent_gravity * sv_gravity.value * host_frametime;
}

trace_t SV_PushEntity(edict_t *ent, vec3_t push)
{ // Does not change the entities velocity at all
	vec3_t end;
	VectorAdd(ent->v.origin, push, end);
	trace_t trace;
	if (ent->v.movetype == MOVETYPE_FLYMISSILE)
		trace = SV_Move(ent->v.origin, ent->v.mins, ent->v.maxs, end,
					MOVE_MISSILE, ent);
	else if (ent->v.solid == SOLID_TRIGGER || ent->v.solid == SOLID_NOT)
		// only clip against bmodels
		trace = SV_Move(ent->v.origin, ent->v.mins, ent->v.maxs, end,
					MOVE_NOMONSTERS, ent);
	else
		trace = SV_Move(ent->v.origin, ent->v.mins, ent->v.maxs, end,
					MOVE_NORMAL, ent);
	VectorCopy(trace.endpos, ent->v.origin);
	SV_LinkEdict(ent, true);
	if (trace.ent)
		SV_Impact(ent, trace.ent);
	return trace;
}

void SV_PushMove(edict_t *pusher, float movetime)
{
	int i, e;
	edict_t *check, *block;
	vec3_t mins, maxs, move;
	vec3_t entorig, pushorig;
	int num_moved;
	edict_t **moved_edict; //johnfitz -- dynamically allocate
	vec3_t *moved_from; //johnfitz -- dynamically allocate
	int mark; //johnfitz
	if (!pusher->v.velocity[0] && !pusher->v.velocity[1]
			&& !pusher->v.velocity[2]) {
		pusher->v.ltime += movetime;
		return;
	}
	for (i = 0; i < 3; i++) {
		move[i] = pusher->v.velocity[i] * movetime;
		mins[i] = pusher->v.absmin[i] + move[i];
		maxs[i] = pusher->v.absmax[i] + move[i];
	}
	VectorCopy(pusher->v.origin, pushorig);
	// move the pusher to it's final position
	VectorAdd(pusher->v.origin, move, pusher->v.origin);
	pusher->v.ltime += movetime;
	SV_LinkEdict(pusher, false);
	//johnfitz -- dynamically allocate
	mark = Hunk_LowMark();
	moved_edict = Hunk_Alloc(sv.num_edicts * sizeof(edict_t *));
	moved_from = Hunk_Alloc(sv.num_edicts * sizeof(vec3_t));
	//johnfitz
	// see if any solid entities are inside the final position
	num_moved = 0;
	check = NEXT_EDICT(sv.edicts);
	for (e = 1; e < sv.num_edicts; e++, check = NEXT_EDICT(check)) {
		if (check->free)
			continue;
		if (check->v.movetype == MOVETYPE_PUSH
				|| check->v.movetype == MOVETYPE_NONE
				|| check->v.movetype == MOVETYPE_NOCLIP)
			continue;
		// if the entity is standing on the pusher, it will definately be moved
		if (!(((int)check->v.flags & FL_ONGROUND)
					&& PROG_TO_EDICT(check->v.groundentity) == pusher)) {
			if (check->v.absmin[0] >= maxs[0]
					|| check->v.absmin[1] >= maxs[1]
					|| check->v.absmin[2] >= maxs[2]
					|| check->v.absmax[0] <= mins[0]
					|| check->v.absmax[1] <= mins[1]
					|| check->v.absmax[2] <= mins[2])
				continue;
			// see if the ent's bbox is inside the pusher's final position
			if (!SV_TestEntityPosition(check))
				continue;
		}
		// remove the onground flag for non-players
		if (check->v.movetype != MOVETYPE_WALK)
			check->v.flags = (int)check->v.flags & ~FL_ONGROUND;
		VectorCopy(check->v.origin, entorig);
		VectorCopy(check->v.origin, moved_from[num_moved]);
		moved_edict[num_moved] = check;
		num_moved++;
		// try moving the contacted entity
		pusher->v.solid = SOLID_NOT;
		SV_PushEntity(check, move);
		pusher->v.solid = SOLID_BSP;
		// if it is still inside the pusher, block
		block = SV_TestEntityPosition(check);
		if (block) { // fail the move
			if (check->v.mins[0] == check->v.maxs[0])
				continue;
			if (check->v.solid == SOLID_NOT || check->v.solid == SOLID_TRIGGER) { // corpse
				check->v.mins[0] = check->v.mins[1] = 0;
				VectorCopy(check->v.mins, check->v.maxs);
				continue;
			}
			VectorCopy(entorig, check->v.origin);
			SV_LinkEdict(check, true);
			VectorCopy(pushorig, pusher->v.origin);
			SV_LinkEdict(pusher, false);
			pusher->v.ltime -= movetime;
			// if the pusher has a "blocked" function, call it
			// otherwise, just stay in place until the obstacle is gone
			if (pusher->v.blocked) {
				pr_global_struct->self = EDICT_TO_PROG(pusher);
				pr_global_struct->other = EDICT_TO_PROG(check);
				PR_ExecuteProgram(pusher->v.blocked);
			}
			// move back any entities we already moved
			for (i = 0; i < num_moved; i++) {
				VectorCopy(moved_from[i],
						moved_edict[i]->v.origin);
				SV_LinkEdict(moved_edict[i], false);
			}
			Hunk_FreeToLowMark(mark); //johnfitz
			return;
		}
	}
	Hunk_FreeToLowMark(mark); //johnfitz
}

void SV_Physics_Pusher(edict_t *ent)
{
	float oldltime = ent->v.ltime;
	float thinktime = ent->v.nextthink;
	float movetime;
	if (thinktime < ent->v.ltime + host_frametime) {
		movetime = thinktime - ent->v.ltime;
		if (movetime < 0)
			movetime = 0;
	} else
		movetime = host_frametime;
	if (movetime) {
		SV_PushMove(ent, movetime); // advances ent->v.ltime if not blocked
	}
	if (thinktime > oldltime && thinktime <= ent->v.ltime) {
		ent->v.nextthink = 0;
		pr_global_struct->time = sv.time;
		pr_global_struct->self = EDICT_TO_PROG(ent);
		pr_global_struct->other = EDICT_TO_PROG(sv.edicts);
		PR_ExecuteProgram(ent->v.think);
		if (ent->free)
			return;
	}
}

// This is a big hack to try and fix the rare case of getting stuck in the
// world clipping hull.
void SV_CheckStuck(edict_t *ent)
{
	if (!SV_TestEntityPosition(ent)) {
		VectorCopy(ent->v.origin, ent->v.oldorigin);
		return;
	}
	vec3_t org;
	VectorCopy(ent->v.origin, org);
	VectorCopy(ent->v.oldorigin, ent->v.origin);
	if (!SV_TestEntityPosition(ent)) {
		Con_DPrintf("Unstuck.\n");
		SV_LinkEdict(ent, true);
		return;
	}
	for (int z = 0; z < 18; z++)
		for (int i = -1; i <= 1; i++)
			for (int j = -1; j <= 1; j++) {
				ent->v.origin[0] = org[0] + i;
				ent->v.origin[1] = org[1] + j;
				ent->v.origin[2] = org[2] + z;
				if (!SV_TestEntityPosition(ent)) {
					Con_DPrintf("Unstuck.\n");
					SV_LinkEdict(ent, true);
					return;
				}
			}
	VectorCopy(org, ent->v.origin);
	Con_DPrintf("player is stuck.\n");
}

qboolean SV_CheckWater(edict_t *ent)
{
	vec3_t point;
	point[0] = ent->v.origin[0];
	point[1] = ent->v.origin[1];
	point[2] = ent->v.origin[2] + ent->v.mins[2] + 1;
	ent->v.waterlevel = 0;
	ent->v.watertype = CONTENTS_EMPTY;
	int cont = SV_PointContents(point);
	if (cont <= CONTENTS_WATER) {
		ent->v.watertype = cont;
		ent->v.waterlevel = 1;
		point[2] =
			ent->v.origin[2] + (ent->v.mins[2] + ent->v.maxs[2]) * 0.5;
		cont = SV_PointContents(point);
		if (cont <= CONTENTS_WATER) {
			ent->v.waterlevel = 2;
			point[2] = ent->v.origin[2] + ent->v.view_ofs[2];
			cont = SV_PointContents(point);
			if (cont <= CONTENTS_WATER)
				ent->v.waterlevel = 3;
		}
	}
	return ent->v.waterlevel > 1;
}

void SV_WallFriction(edict_t *ent, trace_t *trace)
{
	vec3_t forward, right, up;
	AngleVectors(ent->v.v_angle, forward, right, up);
	float d = DotProduct(trace->plane.normal, forward);
	d += 0.5;
	if (d >= 0)
		return;
	// cut the tangential velocity
	float i = DotProduct(trace->plane.normal, ent->v.velocity);
	vec3_t into, side;
	VectorScale(trace->plane.normal, i, into);
	VectorSubtract(ent->v.velocity, into, side);
	ent->v.velocity[0] = side[0] * (1 + d);
	ent->v.velocity[1] = side[1] * (1 + d);
}

// Player has come to a dead stop, possibly due to the problem with limited
// float precision at some angle joins in the BSP hull.
// Try fixing by pushing one pixel in each direction.
// This is a hack, but in the interest of good gameplay...
int SV_TryUnstick(edict_t *ent, vec3_t oldvel)
{
	vec3_t oldorg;
	vec3_t dir;
	trace_t steptrace;
	VectorCopy(ent->v.origin, oldorg);
	VectorCopy(vec3_origin, dir);
	for (int i = 0; i < 8; i++) {
		switch (i) { // try pushing a little in an axial direction
			case 0: dir[0] =  2; dir[1] =  0; break;
			case 1: dir[0] =  0; dir[1] =  2; break;
			case 2: dir[0] = -2; dir[1] =  0; break;
			case 3: dir[0] =  0; dir[1] = -2; break;
			case 4: dir[0] =  2; dir[1] =  2; break;
			case 5: dir[0] = -2; dir[1] =  2; break;
			case 6: dir[0] =  2; dir[1] = -2; break;
			case 7: dir[0] = -2; dir[1] = -2; break;
		}
		SV_PushEntity(ent, dir);
		ent->v.velocity[0] = oldvel[0]; // retry the original move
		ent->v.velocity[1] = oldvel[1];
		ent->v.velocity[2] = 0;
		int clip = SV_FlyMove(ent, 0.1, &steptrace);
		if (fabs(oldorg[1] - ent->v.origin[1]) > 4
				|| fabs(oldorg[0] - ent->v.origin[0]) > 4) {
			return clip;
		}
		// go back to the original pos and try again
		VectorCopy(oldorg, ent->v.origin);
	}
	VectorCopy(vec3_origin, ent->v.velocity);
	return 7; // still not moving
}

void SV_WalkMove(edict_t *ent)
{ // Only used by players
	vec3_t upmove, downmove;
	vec3_t oldorg, oldvel;
	vec3_t nosteporg, nostepvel;
	trace_t steptrace, downtrace;
	// do a regular slide move unless it looks like you ran into a step
	int oldonground = (int)ent->v.flags & FL_ONGROUND;
	ent->v.flags = (int)ent->v.flags & ~FL_ONGROUND;
	VectorCopy(ent->v.origin, oldorg);
	VectorCopy(ent->v.velocity, oldvel);
	int clip = SV_FlyMove(ent, host_frametime, &steptrace);
	if (!(clip & 2))
		return; // move didn't block on a step
	if (!oldonground && ent->v.waterlevel == 0)
		return; // don't stair up while jumping
	if (ent->v.movetype != MOVETYPE_WALK)
		return; // gibbed by a trigger
	if (sv_nostep.value)
		return;
	if ((int)sv_player->v.flags & FL_WATERJUMP)
		return;
	VectorCopy(ent->v.origin, nosteporg);
	VectorCopy(ent->v.velocity, nostepvel);
	// try moving up and forward to go up a step
	VectorCopy(oldorg, ent->v.origin); // back to start pos
	VectorCopy(vec3_origin, upmove);
	VectorCopy(vec3_origin, downmove);
	upmove[2] = STEPSIZE;
	downmove[2] = -STEPSIZE + oldvel[2] * host_frametime;
	// move up
	SV_PushEntity(ent, upmove); // FIXME: don't link?
				    // move forward
	ent->v.velocity[0] = oldvel[0];
	ent->v.velocity[1] = oldvel[1];
	ent->v.velocity[2] = 0;
	clip = SV_FlyMove(ent, host_frametime, &steptrace);
	// check for stuckness, possibly due to the limited precision of floats
	// in the clipping hulls
	if (clip) {
		if (fabs(oldorg[1]-ent->v.origin[1]) < 0.03125 && fabs(oldorg[0]
			- ent->v.origin[0]) < 0.03125) { // stepping up didn't make any progress
			clip = SV_TryUnstick(ent, oldvel);
		}
	}
	if (clip & 2) // extra friction based on view angle
		SV_WallFriction(ent, &steptrace);
	// move down
	downtrace = SV_PushEntity(ent, downmove); // FIXME: don't link?
	if (downtrace.plane.normal[2] > 0.7) {
		if (ent->v.solid == SOLID_BSP) {
			ent->v.flags = (int)ent->v.flags | FL_ONGROUND;
			ent->v.groundentity = EDICT_TO_PROG(downtrace.ent);
		}
	} else {
		// if the push down didn't end up on good ground, use the move without
		// the step up. This happens near wall / slope combinations, and can
		// cause the player to hop up higher on a slope too steep to climb
		VectorCopy(nosteporg, ent->v.origin);
		VectorCopy(nostepvel, ent->v.velocity);
	}
}

void SV_Physics_Client(edict_t *ent, int num)
{ // Player character actions
	if (!svs.clients[num - 1].active)
		return; // unconnected slot

	pr_global_struct->time = sv.time; // call standard client pre-think
	pr_global_struct->self = EDICT_TO_PROG(ent);
	PR_ExecuteProgram(pr_global_struct->PlayerPreThink);
	SV_CheckVelocity(ent); // do a move

	switch ((int)ent->v.movetype) { // decide which move function to call
		case MOVETYPE_NONE:
			if (!SV_RunThink(ent))
				return;
			break;
		case MOVETYPE_WALK:
			if (!SV_RunThink(ent))
				return;
			if (!SV_CheckWater(ent) &&
					!((int)ent->v.flags & FL_WATERJUMP))
				SV_AddGravity(ent);
			SV_CheckStuck(ent);
			SV_WalkMove(ent);
			break;
		case MOVETYPE_TOSS:
		case MOVETYPE_BOUNCE:
			SV_Physics_Toss(ent);
			break;
		case MOVETYPE_FLY:
			if (!SV_RunThink(ent))
				return;
			SV_FlyMove(ent, host_frametime, NULL);
			break;
		case MOVETYPE_NOCLIP:
			if (!SV_RunThink(ent))
				return;
			VectorMA(ent->v.origin, host_frametime, ent->v.velocity,
					ent->v.origin);
			break;
		default:
			Sys_Error("SV_Physics_client: bad movetype %i",
					(int)ent->v.movetype);
	}
	SV_LinkEdict(ent, true); // call standard player post-think
	pr_global_struct->time = sv.time;
	pr_global_struct->self = EDICT_TO_PROG(ent);
	PR_ExecuteProgram(pr_global_struct->PlayerPostThink);
}

void SV_Physics_None(edict_t *ent) // Non moving objects can only think
{
	SV_RunThink(ent); // regular thinking
}


void SV_Physics_Noclip(edict_t *ent)
{ // A moving object that doesn't obey physics
	if (!SV_RunThink(ent)) // regular thinking
		return;
	VectorMA(ent->v.angles, host_frametime, ent->v.avelocity,ent->v.angles);
	VectorMA(ent->v.origin, host_frametime, ent->v.velocity, ent->v.origin);
	SV_LinkEdict(ent, false);
}

void SV_CheckWaterTransition(edict_t *ent)
{
	int cont = SV_PointContents(ent->v.origin);
	if (!ent->v.watertype) { // just spawned here
		ent->v.watertype = cont;
		ent->v.waterlevel = 1;
		return;
	}
	if (cont <= CONTENTS_WATER) {
		if (ent->v.watertype == CONTENTS_EMPTY) { // just crossed into water
			SV_StartSound(ent, 0, "misc/h2ohit1.wav", 255, 1);
		}
		ent->v.watertype = cont;
		ent->v.waterlevel = 1;
	} else {
		if (ent->v.watertype != CONTENTS_EMPTY) { // just crossed into water
			SV_StartSound(ent, 0, "misc/h2ohit1.wav", 255, 1);
		}
		ent->v.watertype = CONTENTS_EMPTY;
		ent->v.waterlevel = cont;
	}
}


void SV_Physics_Toss(edict_t *ent)
{ // Toss, bounce, and fly movement. When onground, do nothing.
	trace_t trace;
	vec3_t move;
	float backoff;
	if (!SV_RunThink(ent)) // regular thinking
		return;
	if (((int)ent->v.flags & FL_ONGROUND)) // if onground ret without moving
		return;
	SV_CheckVelocity(ent);
	if (ent->v.movetype != MOVETYPE_FLY // add gravity
			&& ent->v.movetype != MOVETYPE_FLYMISSILE)
		SV_AddGravity(ent);
	VectorMA(ent->v.angles, host_frametime, ent->v.avelocity, // move angles
			ent->v.angles);
	VectorScale(ent->v.velocity, host_frametime, move); // move origin
	trace = SV_PushEntity(ent, move);
	if (trace.fraction == 1)
		return;
	if (ent->free)
		return;
	if (ent->v.movetype == MOVETYPE_BOUNCE)
		backoff = 1.5;
	else
		backoff = 1;
	ClipVelocity(ent->v.velocity, trace.plane.normal, ent->v.velocity,
			backoff);

	if (trace.plane.normal[2] > 0.7) { // stop if on ground
		if (ent->v.velocity[2] < 60
				|| ent->v.movetype != MOVETYPE_BOUNCE) {
			ent->v.flags = (int)ent->v.flags | FL_ONGROUND;
			ent->v.groundentity = EDICT_TO_PROG(trace.ent);
			VectorCopy(vec3_origin, ent->v.velocity);
			VectorCopy(vec3_origin, ent->v.avelocity);
		}
	}
	SV_CheckWaterTransition(ent); // check for in water
}

// Monsters freefall when they don't have a ground entity, otherwise
// all movement is done with discrete steps.
// This is also used for objects that have become still on the ground, but
// will fall if the floor is pulled out from under them.
void SV_Physics_Step(edict_t *ent)
{
	if (!((int)ent->v.flags & (FL_ONGROUND | FL_FLY | FL_SWIM))) { // freefall if not onground
		qboolean hitsound = (ent->v.velocity[2]<sv_gravity.value*-0.1);
		SV_AddGravity(ent);
		SV_CheckVelocity(ent);
		SV_FlyMove(ent, host_frametime, NULL);
		SV_LinkEdict(ent, true);
		if ((int)ent->v.flags & FL_ONGROUND && hitsound) // just hit ground
			SV_StartSound(ent, 0, "demon/dland2.wav",255,1);
	}
	SV_RunThink(ent); // regular thinking
	SV_CheckWaterTransition(ent);
}

void SV_Physics()
{
	// let the progs know that a new frame has started
	pr_global_struct->self = EDICT_TO_PROG(sv.edicts);
	pr_global_struct->other = EDICT_TO_PROG(sv.edicts);
	pr_global_struct->time = sv.time;
	PR_ExecuteProgram(pr_global_struct->StartFrame);
	edict_t *ent = sv.edicts; // treat each object in turn
	for (int i = 0; i < sv.num_edicts; i++, ent = NEXT_EDICT(ent)) {
		if (ent->free)
			continue;
		if (pr_global_struct->force_retouch) {
			SV_LinkEdict(ent, true); // force retouch even for stationary
		}
		if (i > 0 && i <= svs.maxclients)
			SV_Physics_Client(ent, i);
		else if (ent->v.movetype == MOVETYPE_PUSH)
			SV_Physics_Pusher(ent);
		else if (ent->v.movetype == MOVETYPE_NONE)
			SV_Physics_None(ent);
		else if (ent->v.movetype == MOVETYPE_NOCLIP)
			SV_Physics_Noclip(ent);
		else if (ent->v.movetype == MOVETYPE_STEP)
			SV_Physics_Step(ent);
		else if (ent->v.movetype == MOVETYPE_TOSS
				|| ent->v.movetype == MOVETYPE_BOUNCE
				|| ent->v.movetype == MOVETYPE_FLY
				|| ent->v.movetype == MOVETYPE_FLYMISSILE)
			SV_Physics_Toss(ent);
		else
			Sys_Error("SV_Physics: bad movetype %i",
					(int)ent->v.movetype);
	}
	if (pr_global_struct->force_retouch)
		pr_global_struct->force_retouch--;
	sv.time += host_frametime;
}
