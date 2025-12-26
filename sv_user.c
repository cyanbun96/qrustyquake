// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
// sv_user.c -- server code for moving users
#include "quakedef.h"

static vec3_t sv_fwd, sv_rgt, sv_up;

void SV_SetIdealPitch()
{
	f32 z[MAX_FORWARD];
	if(!((s32)sv_player->v.flags & FL_ONGROUND)) return;
	f32 angleval = sv_player->v.angles[YAW] * M_PI*2 / 360;
	f32 sinval = sin(angleval);
	f32 cosval = cos(angleval);
	s32 i = 0;
	for(; i < MAX_FORWARD; i++) {
		vec3_t top, bottom;
		top[0] = sv_player->v.origin[0] + cosval*(i+3)*12;
		top[1] = sv_player->v.origin[1] + sinval*(i+3)*12;
		top[2] = sv_player->v.origin[2] + sv_player->v.view_ofs[2];
		bottom[0] = top[0];
		bottom[1] = top[1];
		bottom[2] = top[2] - 160;
	     trace_t tr=SV_Move(top,vec3_origin,vec3_origin,bottom,1,sv_player);
		if(tr.allsolid) return; // looking at a wall
		if(tr.fraction == 1) return; // near a dropoff
		z[i] = top[2] + tr.fraction*(bottom[2]-top[2]);
	}
	s32 dir = 0;
	s32 steps = 0;
	for(s32 j = 1; j < i; j++) {
		s32 step = z[j] - z[j-1];
		if(step > -ON_EPSILON && step < ON_EPSILON) continue;
		if(dir && ( step-dir > ON_EPSILON || step-dir < -ON_EPSILON ) )
			return; // mixed changes
		steps++;
		dir = step;
	}
	if(!dir) {
		sv_player->v.idealpitch = 0;
		return;
	}
	if(steps < 2) return;
	sv_player->v.idealpitch = -dir * sv_idealpitchscale.value;
}

void SV_UserFriction()
{
	f32 *vel = velocity;
	f32 speed = sqrt(vel[0]*vel[0] + vel[1]*vel[1]);
	if(!speed) return;
	// if the leading edge is over a dropoff, increase friction
	vec3_t start, stop;
	start[0] = stop[0] = origin[0] + vel[0]/speed*16;
	start[1] = stop[1] = origin[1] + vel[1]/speed*16;
	start[2] = origin[2] + sv_player->v.mins[2];
	stop[2] = start[2] - 34;
	trace_t trace=SV_Move(start,vec3_origin,vec3_origin,stop,1,sv_player);
	f32 friction = sv_friction.value;
	if(trace.fraction == 1.0) friction *= sv_edgefriction.value;
	// apply friction
	f32 control = speed < sv_stopspeed.value ? sv_stopspeed.value : speed;
	f32 newspeed = speed - host_frametime*control*friction;
	if(newspeed < 0) newspeed = 0;
	newspeed /= speed;
	vel[0] = vel[0] * newspeed;
	vel[1] = vel[1] * newspeed;
	vel[2] = vel[2] * newspeed;
}

void SV_Accelerate(f32 wishspeed, const vec3_t wishdir)
{
	f32 currentspeed = DotProduct(velocity, wishdir);
	f32 addspeed = wishspeed - currentspeed;
	if(addspeed <= 0) return;
	f32 accelspeed = sv_accelerate.value*host_frametime*wishspeed;
	if(accelspeed > addspeed) accelspeed = addspeed;
	for(s32 i = 0; i < 3; i++) velocity[i] += accelspeed*wishdir[i];
}

void SV_AirAccelerate(f32 wishspeed, vec3_t wishveloc)
{
	f32 wishspd = VectorNormalize(wishveloc);
	if(wishspd > 30) wishspd = 30;
	f32 currentspeed = DotProduct(velocity, wishveloc);
	f32 addspeed = wishspd - currentspeed;
	if(addspeed <= 0) return;
	f32 accelspeed = sv_accelerate.value*wishspeed * host_frametime;
	if(accelspeed > addspeed) accelspeed = addspeed;
	for(s32 i = 0; i < 3; i++) velocity[i] += accelspeed*wishveloc[i];
}

void DropPunchAngle()
{
	f32 len = VectorNormalize(sv_player->v.punchangle);
	len -= 10*host_frametime;
	if(len < 0) len = 0;
	VectorScale(sv_player->v.punchangle, len, sv_player->v.punchangle);
}

void SV_WaterMove()
{
	AngleVectors(sv_player->v.v_angle, sv_fwd, sv_rgt, sv_up); //user intent
	vec3_t wishvel;
	for(s32 i = 0; i < 3; i++)
		wishvel[i] = sv_fwd[i]*cmd.forwardmove + sv_rgt[i]*cmd.sidemove;
	if(!cmd.forwardmove && !cmd.sidemove && !cmd.upmove)
		wishvel[2] -= 60; // drift towards bottom
	else wishvel[2] += cmd.upmove;
	f32 wishspeed = VectorLength(wishvel);
	if(wishspeed > sv_maxspeed.value) {
		VectorScale(wishvel, sv_maxspeed.value/wishspeed, wishvel);
		wishspeed = sv_maxspeed.value;
	}
	wishspeed *= 0.7;
	f32 speed = VectorLength(velocity); // water friction
	f32 newspeed = 0;
	if(speed) {
		newspeed = speed - host_frametime * speed * sv_friction.value;
		if(newspeed < 0) newspeed = 0;
		VectorScale(velocity, newspeed/speed, velocity);
	}
	// water acceleration
	if(!wishspeed) return;
	f32 addspeed = wishspeed - newspeed;
	if(addspeed <= 0) return;
	VectorNormalize(wishvel);
	f32 accelspeed = sv_accelerate.value * wishspeed * host_frametime;
	if(accelspeed > addspeed)
		accelspeed = addspeed;
	for(s32 i = 0; i < 3; i++)
		velocity[i] += accelspeed * wishvel[i];
}

void SV_WaterJump()
{
	if(sv.time > sv_player->v.teleport_time || !sv_player->v.waterlevel) {
		sv_player->v.flags = (s32)sv_player->v.flags & ~FL_WATERJUMP;
		sv_player->v.teleport_time = 0;
	}
	sv_player->v.velocity[0] = sv_player->v.movedir[0];
	sv_player->v.velocity[1] = sv_player->v.movedir[1];
}

void SV_NoclipMove()
{ // johnfitz - new, alternate noclip. old noclip is still handled in SV_AirMove
	AngleVectors(sv_player->v.v_angle, sv_fwd, sv_rgt, sv_up);
	velocity[0] = sv_fwd[0]*cmd.forwardmove + sv_rgt[0]*cmd.sidemove;
	velocity[1] = sv_fwd[1]*cmd.forwardmove + sv_rgt[1]*cmd.sidemove;
	velocity[2] = sv_fwd[2]*cmd.forwardmove + sv_rgt[2]*cmd.sidemove;
	velocity[2] += cmd.upmove*2; //doubled to match running speed
	if(VectorLength(velocity) > sv_maxspeed.value) {
		VectorNormalize(velocity);
		VectorScale(velocity, sv_maxspeed.value, velocity);
	}
}

void SV_AirMove()
{
	AngleVectors(sv_player->v.angles, sv_fwd, sv_rgt, sv_up);
	f32 fmove = cmd.forwardmove;
	f32 smove = cmd.sidemove;
	// hack to not let you back into teleporter
	if(sv.time < sv_player->v.teleport_time && fmove < 0) fmove = 0;
	vec3_t wishvel, wishdir;
	for(s32 i = 0; i < 3; i++)
		wishvel[i] = sv_fwd[i]*fmove + sv_rgt[i]*smove;
	if((s32)sv_player->v.movetype != MOVETYPE_WALK) wishvel[2] = cmd.upmove;
	else wishvel[2] = 0;
	VectorCopy(wishvel, wishdir);
	f32 wishspeed = VectorNormalize(wishdir);
	if(wishspeed > sv_maxspeed.value) {
		VectorScale(wishvel, sv_maxspeed.value/wishspeed, wishvel);
		wishspeed = sv_maxspeed.value;
	}
	if(sv_player->v.movetype == MOVETYPE_NOCLIP) { // noclip
		VectorCopy(wishvel, velocity);
	} else if(onground) {
		SV_UserFriction();
		SV_Accelerate(wishspeed, wishdir);
	} else // not on ground, so little effect on velocity
		SV_AirAccelerate(wishspeed, wishvel);
}

void SV_ClientThink() // the move fields specify an intended velocity in pix/sec
{ // the angle fields specify an exact angular motion in degrees
	if(sv_player->v.movetype == MOVETYPE_NONE) return;
	onground = (s32)sv_player->v.flags & FL_ONGROUND;
	origin = sv_player->v.origin;
	velocity = sv_player->v.velocity;
	DropPunchAngle();
	// if dead, behave differently
	if(sv_player->v.health <= 0) return;
	// angles
	// show 1/3 the pitch angle and all the roll angle
	cmd = host_client->cmd;
	angles = sv_player->v.angles;
	vec3_t v_angle;
	VectorAdd(sv_player->v.v_angle, sv_player->v.punchangle, v_angle);
	angles[ROLL] = V_CalcRoll(sv_player->v.angles, sv_player->v.velocity)*4;
	if(!sv_player->v.fixangle) {
		angles[PITCH] = -v_angle[PITCH]/3;
		angles[YAW] = v_angle[YAW];
	}
	if((s32)sv_player->v.flags & FL_WATERJUMP) {
		SV_WaterJump();
		return;
	}
	if(sv_player->v.movetype == MOVETYPE_NOCLIP && sv_altnoclip.value)//walk
		SV_NoclipMove();
     else if(sv_player->v.waterlevel>=2&&sv_player->v.movetype!=MOVETYPE_NOCLIP)
		SV_WaterMove();
	else SV_AirMove();
}

void SV_ReadClientMove(usercmd_t *move)
{
	host_client->ping_times[host_client->num_pings%NUM_PING_TIMES]
		= sv.time - MSG_ReadFloat(); // read ping time
	host_client->num_pings++;
	vec3_t angle;
	for(s32 i = 0; i < 3; i++) // read current angles
		if(sv.protocol == PROTOCOL_NETQUAKE)
			angle[i] = MSG_ReadAngle(sv.protocolflags);
		else
			angle[i] = MSG_ReadAngle16(sv.protocolflags);
	VectorCopy(angle, host_client->edict->v.v_angle);
	move->forwardmove = MSG_ReadShort(); // read movement
	move->sidemove = MSG_ReadShort();
	move->upmove = MSG_ReadShort();
	s32 bits = MSG_ReadByte(); // read buttons
	host_client->edict->v.button0 = bits & 1;
	host_client->edict->v.button2 = (bits & 2)>>1;
	s32 i = MSG_ReadByte();
	if(i) host_client->edict->v.impulse = i;
}

bool SV_ReadClientMessage()
{ // Returns 0 if the client should be killed
	const s8 *s;
	s32 ret;
	do {
nextmsg:
		ret = NET_GetMessage(host_client->netconnection);
		if(ret == -1) {
		    Sys_Printf("SV_ReadClientMessage: NET_GetMessage failed\n");
			return 0;
		}
		if(!ret) return 1;
		MSG_BeginReading();
		while(1) {
			if(!host_client->active)
				return 0; // a command caused an error
			if(msg_badread) {
				Sys_Printf("SV_ReadClientMessage: badread\n");
				return 0;
			}
			s32 ccmd = MSG_ReadChar();
			switch(ccmd) {
			case -1:
				goto nextmsg; // end of message
			default:
		       Sys_Printf("SV_ReadClientMessage: unknown command s8\n");
				return 0;
			case clc_nop:
				break;
			case clc_stringcmd:
				s = MSG_ReadString();
				ret = 0;
				if(q_strncasecmp(s, "status", 6) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "god", 3) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "notarget", 8) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "fly", 3) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "name", 4) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "noclip", 6) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "setpos", 6) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "say", 3) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "say_team", 8) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "tell", 4) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "color", 5) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "kill", 4) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "pause", 5) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "spawn", 5) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "begin", 5) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "prespawn", 8) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "kick", 4) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "ping", 4) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "give", 4) == 0)
					ret = 1;
				else if(q_strncasecmp(s, "ban", 3) == 0)
					ret = 1;
				if(ret == 1)
					Cmd_ExecuteString(s, src_client);
				else
					Con_DPrintf("%s tried to %s\n",
							host_client->name, s);
				break;
			case clc_disconnect:
				return 0;
			case clc_move:
				SV_ReadClientMove(&host_client->cmd);
				break;
			}
		}
	} while(ret == 1);
	return 1;
}

void SV_RunClients()
{
	s32 i = 0;
	for(host_client = svs.clients ; i<svs.maxclients ; i++, host_client++) {
		if(!host_client->active) continue;
		sv_player = host_client->edict;
		if(!SV_ReadClientMessage()) {
			SV_DropClient(0); // client misbehaved...
			continue;
		}
		if(!host_client->spawned) {
			// clear client movement until a new packet is received
			memset(&host_client->cmd, 0, sizeof(host_client->cmd));
			continue;
		}
		// always pause in single player if in console or menus
		if(!sv.paused && (svs.maxclients > 1 || key_dest == key_game) )
			SV_ClientThink();
	}
}
