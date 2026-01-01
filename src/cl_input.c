// Copyright (C) 1996-2001 Id Software, Inc.
// Copyright (C) 2002-2009 John Fitzgibbons and others
// Copyright (C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
#include "quakedef.h"

static s32 in_impulse;
static kbutton_t in_left, in_right, in_forward, in_back, in_klook;
static kbutton_t in_lookup, in_lookdown, in_moveleft, in_moveright;
static kbutton_t in_speed, in_use, in_jump, in_attack;
static kbutton_t in_up, in_down;

//Continuous button event tracking is complicated by the fact that two different
//input sources (say, mouse button 1 and the control key) can both press the
//same button, but the button should only be released when both of the
//pressing key have been released.
//When a key event issues a button command (+forward, +attack, etc), it appends
//its key number as a parameter to the command so it can be matched up with
//the release.
//state bit 0 is the current state of the key
//state bit 1 is edge triggered on the up to down transition
//state bit 2 is edge triggered on the down to up transition

void KeyDown (kbutton_t *b)
{
	const s8 *c = Cmd_Argv(1); // typed manually in console
	s32 k = c[0] ? atoi(c) : -1; // for continuous down
	if(k == b->down[0] || k == b->down[1]) return; // repeating key
	if(!b->down[0]) b->down[0] = k;
	else if(!b->down[1]) b->down[1] = k;
	else {
		Con_Printf ("Three keys down for a button!\n");
		return;
	}
	if(b->state & 1) return; // still down
	b->state |= 1 + 2; // down + impulse down
}

void KeyUp (kbutton_t *b)
{
	const s8 *c = Cmd_Argv(1);
	if(!c[0]){ //typed manually in console, assume for unsticking, clear all
		b->down[0] = b->down[1] = 0;
		b->state = 4; // impulse up
		return;
	}
	s32 k = atoi(c);
	if(b->down[0] == k) b->down[0] = 0;
	else if(b->down[1] == k) b->down[1] = 0;
	else return; // key up without coresponding down (menu pass through)
	if(b->down[0]||b->down[1]) return; // some other key is holding it down
	if(!(b->state & 1)) return; // still up (this should not happen)
	b->state &= ~1; // now up
	b->state |= 4; // impulse up
}

void IN_KLookDown (){KeyDown(&in_klook);}
void IN_KLookUp (){KeyUp(&in_klook);}
void IN_MLookDown (){KeyDown(&in_mlook);}
void IN_MLookUp (){ KeyUp(&in_mlook);
	if( !(in_mlook.state&1) && lookspring.value)
		V_StartPitchDrift(); }
void IN_UpDown(){KeyDown(&in_up);}
void IN_UpUp(){KeyUp(&in_up);}
void IN_DownDown(){KeyDown(&in_down);}
void IN_DownUp(){KeyUp(&in_down);}
void IN_LeftDown(){KeyDown(&in_left);}
void IN_LeftUp(){KeyUp(&in_left);}
void IN_RightDown(){KeyDown(&in_right);}
void IN_RightUp(){KeyUp(&in_right);}
void IN_ForwardDown(){KeyDown(&in_forward);}
void IN_ForwardUp(){KeyUp(&in_forward);}
void IN_BackDown(){KeyDown(&in_back);}
void IN_BackUp(){KeyUp(&in_back);}
void IN_LookupDown(){KeyDown(&in_lookup);}
void IN_LookupUp(){KeyUp(&in_lookup);}
void IN_LookdownDown(){KeyDown(&in_lookdown);}
void IN_LookdownUp(){KeyUp(&in_lookdown);}
void IN_MoveleftDown(){KeyDown(&in_moveleft);}
void IN_MoveleftUp(){KeyUp(&in_moveleft);}
void IN_MoverightDown(){KeyDown(&in_moveright);}
void IN_MoverightUp(){KeyUp(&in_moveright);}
void IN_SpeedDown(){KeyDown(&in_speed);}
void IN_SpeedUp(){KeyUp(&in_speed);}
void IN_StrafeDown(){KeyDown(&in_strafe);}
void IN_StrafeUp(){KeyUp(&in_strafe);}
void IN_AttackDown(){KeyDown(&in_attack);}
void IN_AttackUp(){KeyUp(&in_attack);}
void IN_UseDown (){KeyDown(&in_use);}
void IN_UseUp (){KeyUp(&in_use);}
void IN_JumpDown (){KeyDown(&in_jump);}
void IN_JumpUp (){KeyUp(&in_jump);}
void IN_Impulse (){in_impulse=Q_atoi(Cmd_Argv(1));}

f32 CL_KeyState (kbutton_t *key)
{
	bool impulsedown = key->state & 2;
	bool impulseup = key->state & 4;
	bool down = key->state & 1;
	f32 val = 0;
	if(impulsedown && !impulseup){
		if(down) val = 0.5; // pressed and held this frame
		else val = 0; // I_Error ();
	}
	if(impulseup && !impulsedown){
		if(down) val = 0; // I_Error ();
		else val = 0; // released this frame
	}
	if(!impulsedown && !impulseup){
		if(down) val = 1.0; // held the entire frame
		else val = 0; // up the entire frame
	}
	if(impulsedown && impulseup){
		if(down) val = 0.75; // released and re-pressed this frame
		else val = 0.25; // pressed and released this frame
	}
	key->state &= 1; // clear impulses
	return val;
}

void CL_AdjustAngles ()
{ // Moves the local angle positions
	f32 speed = host_frametime;
	if((in_speed.state & 1) ^ (cl_alwaysrun.value != 0.0))
		speed *= cl_anglespeedkey.value;
	if(!(in_strafe.state & 1)){
		cl.viewangles[YAW] -= speed*cl_yawspeed.value
			* CL_KeyState (&in_right);
		cl.viewangles[YAW] += speed*cl_yawspeed.value
			* CL_KeyState (&in_left);
		cl.viewangles[YAW] = anglemod(cl.viewangles[YAW]);
	}
	if(in_klook.state & 1){
		V_StopPitchDrift ();
		cl.viewangles[PITCH] -= speed*cl_pitchspeed.value
			* CL_KeyState (&in_forward);
		cl.viewangles[PITCH] += speed*cl_pitchspeed.value
			* CL_KeyState (&in_back);
	}
	f32 up = CL_KeyState (&in_lookup);
	f32 down = CL_KeyState(&in_lookdown);
	cl.viewangles[PITCH] -= speed*cl_pitchspeed.value * up;
	cl.viewangles[PITCH] += speed*cl_pitchspeed.value * down;
	if(up || down) V_StopPitchDrift ();
	if(cl.viewangles[PITCH] > cl_maxpitch.value) //johnfitz -- variable
		cl.viewangles[PITCH] = cl_maxpitch.value; // pitch clamping
	if(cl.viewangles[PITCH] < cl_minpitch.value)
		cl.viewangles[PITCH] = cl_minpitch.value;
	if(cl.viewangles[ROLL] > 50) cl.viewangles[ROLL] = 50;
	if(cl.viewangles[ROLL] < -50) cl.viewangles[ROLL] = -50;
}


void CL_BaseMove (usercmd_t *cmd)
{ // Send the intended movement message to the server
	if(cls.signon != SIGNONS) return;
	Q_memset (cmd, 0, sizeof(*cmd));
	if(in_strafe.state & 1){
		cmd->sidemove += cl_sidespeed.value * CL_KeyState (&in_right);
		cmd->sidemove -= cl_sidespeed.value * CL_KeyState (&in_left);
	}
	cmd->sidemove += cl_sidespeed.value * CL_KeyState (&in_moveright);
	cmd->sidemove -= cl_sidespeed.value * CL_KeyState (&in_moveleft);
	cmd->upmove += cl_upspeed.value * CL_KeyState (&in_up);
	cmd->upmove -= cl_upspeed.value * CL_KeyState (&in_down);
	if(!(in_klook.state & 1)){
		cmd->forwardmove += cl_forwardspeed.value
			* CL_KeyState (&in_forward);
		cmd->forwardmove -= cl_backspeed.value * CL_KeyState (&in_back);
	}
	if((in_speed.state & 1) ^ (cl_alwaysrun.value != 0.0)){
		cmd->forwardmove *= cl_movespeedkey.value;
		cmd->sidemove *= cl_movespeedkey.value;
		cmd->upmove *= cl_movespeedkey.value;
	}
}

void CL_SendMove (const usercmd_t *cmd)
{
	u8 data[128];
	sizebuf_t buf = {0, 0, data, 128, 0};
	if(cmd){
		cl.cmd = *cmd; // send the movement message
		MSG_WriteByte (&buf, clc_move);
		MSG_WriteFloat (&buf,cl.mtime[0]);//so server can get ping times
		for(s32 i = 0; i < 3; i++) //johnfitz -- 16-bit angles
			if(cl.protocol == PROTOCOL_NETQUAKE)
				MSG_WriteAngle (&buf, cl.viewangles[i],
						cl.protocolflags);
			else MSG_WriteAngle16 (&buf, cl.viewangles[i],
						cl.protocolflags);
		MSG_WriteShort (&buf, cmd->forwardmove);
		MSG_WriteShort (&buf, cmd->sidemove);
		MSG_WriteShort (&buf, cmd->upmove);
		s32 bits = 0; // send button bits
		if(in_attack.state & 3) bits |= 1;
		in_attack.state &= ~2;
		if(in_jump.state & 3) bits |= 2;
		in_jump.state &= ~2;
		MSG_WriteByte (&buf, bits);
		MSG_WriteByte (&buf, in_impulse);
		in_impulse = 0;
	}
	if(cls.demoplayback) return; // deliver the message
	// allways dump the first two message, because it may contain leftover
	// inputs from the last level
	if(++cl.movemessages <= 2) return;
	if(NET_SendUnreliableMessage (cls.netcon, &buf) == -1) {
		Con_Printf ("CL_SendMove: lost server connection\n");
		CL_Disconnect ();
	}
}

void CL_InitInput ()
{
	Cmd_AddCommand ("+moveup",IN_UpDown);
	Cmd_AddCommand ("-moveup",IN_UpUp);
	Cmd_AddCommand ("+movedown",IN_DownDown);
	Cmd_AddCommand ("-movedown",IN_DownUp);
	Cmd_AddCommand ("+left",IN_LeftDown);
	Cmd_AddCommand ("-left",IN_LeftUp);
	Cmd_AddCommand ("+right",IN_RightDown);
	Cmd_AddCommand ("-right",IN_RightUp);
	Cmd_AddCommand ("+forward",IN_ForwardDown);
	Cmd_AddCommand ("-forward",IN_ForwardUp);
	Cmd_AddCommand ("+back",IN_BackDown);
	Cmd_AddCommand ("-back",IN_BackUp);
	Cmd_AddCommand ("+lookup", IN_LookupDown);
	Cmd_AddCommand ("-lookup", IN_LookupUp);
	Cmd_AddCommand ("+lookdown", IN_LookdownDown);
	Cmd_AddCommand ("-lookdown", IN_LookdownUp);
	Cmd_AddCommand ("+strafe", IN_StrafeDown);
	Cmd_AddCommand ("-strafe", IN_StrafeUp);
	Cmd_AddCommand ("+moveleft", IN_MoveleftDown);
	Cmd_AddCommand ("-moveleft", IN_MoveleftUp);
	Cmd_AddCommand ("+moveright", IN_MoverightDown);
	Cmd_AddCommand ("-moveright", IN_MoverightUp);
	Cmd_AddCommand ("+speed", IN_SpeedDown);
	Cmd_AddCommand ("-speed", IN_SpeedUp);
	Cmd_AddCommand ("+attack", IN_AttackDown);
	Cmd_AddCommand ("-attack", IN_AttackUp);
	Cmd_AddCommand ("+use", IN_UseDown);
	Cmd_AddCommand ("-use", IN_UseUp);
	Cmd_AddCommand ("+jump", IN_JumpDown);
	Cmd_AddCommand ("-jump", IN_JumpUp);
	Cmd_AddCommand ("impulse", IN_Impulse);
	Cmd_AddCommand ("+klook", IN_KLookDown);
	Cmd_AddCommand ("-klook", IN_KLookUp);
	Cmd_AddCommand ("+mlook", IN_MLookDown);
	Cmd_AddCommand ("-mlook", IN_MLookUp);
}

