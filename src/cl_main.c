// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// cl_main.c  -- client main loop
#include "quakedef.h"

static efrag_t cl_efrags[MAX_EFRAGS];

void CL_ClearState()
{
	if (!sv.active) Host_ClearMemory();
	memset(&cl, 0, sizeof(cl)); // wipe the entire cl structure
	SZ_Clear(&cls.message);
	memset(cl_efrags, 0, sizeof(cl_efrags)); // clear other arrays
	memset(cl_entities, 0, sizeof(cl_entities));
	memset(cl_dlights, 0, sizeof(cl_dlights));
	memset(cl_lightstyle, 0, sizeof(cl_lightstyle));
	memset(cl_temp_entities, 0, sizeof(cl_temp_entities));
	memset(cl_beams, 0, sizeof(cl_beams));
	// allocate the efrags and chain together into a free list
	cl.free_efrags = cl_efrags;
	s32 i = 0;
	for (; i < MAX_EFRAGS - 1; i++)
		cl.free_efrags[i].entnext = &cl.free_efrags[i + 1];
	cl.free_efrags[i].entnext = NULL;
}

void CL_Disconnect() // Sends a disconnect message to the server
{ // This is also called on Host_Error, so it shouldn't cause any errors
	S_StopAllSounds(1); // stop sounds (especially looping!)
	if (cls.demoplayback) // if running a local server, shut it down
		CL_StopPlayback();
	else if (cls.state == ca_connected) {
		if (cls.demorecording) CL_Stop_f();
		Con_DPrintf("Sending clc_disconnect\n");
		SZ_Clear(&cls.message);
		MSG_WriteByte(&cls.message, clc_disconnect);
		NET_SendUnreliableMessage(cls.netcon, &cls.message);
		SZ_Clear(&cls.message);
		NET_Close(cls.netcon);
		cls.state = ca_disconnected;
		if (sv.active) Host_ShutdownServer(0);
	}
	cls.demoplayback = cls.timedemo = 0;
	cls.signon = 0;
}

void CL_Disconnect_f()
{
	CL_Disconnect();
	if (sv.active) Host_ShutdownServer(0);
}

void CL_EstablishConnection(s8 *host)
{ // Host should be either "local" or a net address to be passed on
	if (cls.state == ca_dedicated || cls.demoplayback) return;
	CL_Disconnect();
	cls.netcon = NET_Connect(host);
	if (!cls.netcon) Host_Error("CL_Connect: connect failed\n");
	Con_DPrintf("CL_EstablishConnection: connected to %s\n", host);
	cls.demonum = -1; // not in the demo loop now
	cls.state = ca_connected;
	cls.signon = 0; // need all the signon messages before playing
}

void CL_SignonReply()
{ // An svc_signonnum has been received, perform a client side setup
	s8 str[8192];
	Con_DPrintf("CL_SignonReply: %i\n", cls.signon);
	switch (cls.signon) {
	case 1:
		MSG_WriteByte(&cls.message, clc_stringcmd);
		MSG_WriteString(&cls.message, "prespawn");
		break;
	case 2:
		MSG_WriteByte(&cls.message, clc_stringcmd);
		MSG_WriteString(&cls.message,
				va("name \"%s\"\n", cl_name.string));
		MSG_WriteByte(&cls.message, clc_stringcmd);
		MSG_WriteString(&cls.message,
				va("color %i %i\n", ((s32)cl_color.value) >> 4,
				   ((s32)cl_color.value) & 15));
		MSG_WriteByte(&cls.message, clc_stringcmd);
		sprintf(str, "spawn %s", cls.spawnparms);
		MSG_WriteString(&cls.message, str);
		break;
	case 3:
		MSG_WriteByte(&cls.message, clc_stringcmd);
		MSG_WriteString(&cls.message, "begin");
		Cache_Report(); // print remaining memory
		break;
	case 4:
		SCR_EndLoadingPlaque(); // allow normal screen updates
		break;
	}
}

void CL_NextDemo()
{ // Called to play the next demo in the demo loop
	s8 str[1024];
	if (cls.demonum == -1) return; // don't play demos
	SCR_BeginLoadingPlaque();
	if (!cls.demos[cls.demonum][0] || cls.demonum == MAX_DEMOS) {
		cls.demonum = 0;
		if (!cls.demos[cls.demonum][0]) {
			Con_Printf("No demos listed with startdemos\n");
			cls.demonum = -1;
			return;
		}
	}
	sprintf(str, "playdemo %s\n", cls.demos[cls.demonum]);
	Cbuf_InsertText(str);
	cls.demonum++;
}

void CL_PrintEntities_f()
{
	entity_t *ent = cl_entities;
	for (s32 i = 0; i < cl.num_entities; i++, ent++) {
		Con_Printf("%3i:", i);
		if (!ent->model) { Con_Printf("EMPTY\n"); continue; }
		Con_Printf("%s:%2i  (%5.1f,%5.1f,%5.1f) [%5.1f %5.1f %5.1f]\n",
			   ent->model->name, ent->frame, ent->origin[0],
			   ent->origin[1], ent->origin[2], ent->angles[0],
			   ent->angles[1], ent->angles[2]);
	}
}

dlight_t *CL_AllocDlight (s32 key)
{
        s32 i;
        dlight_t *dl;
        if (key) { // first look for an exact key match
                dl = cl_dlights;
                for (i=0 ; i<MAX_DLIGHTS ; i++, dl++) {
                        if (dl->key == key) {
                                memset (dl, 0, sizeof(*dl));
                                dl->key = key;
                                dl->color[0] = dl->color[1] = dl->color[2] = 1; //johnfitz -- lit support via lordhavoc
                                return dl;
                        }
                }
        }
        dl = cl_dlights; // then look for anything else
        for (i=0 ; i<MAX_DLIGHTS ; i++, dl++) {
                if (dl->die < cl.time) {
                        memset (dl, 0, sizeof(*dl));
                        dl->key = key;
                        dl->color[0] = dl->color[1] = dl->color[2] = 1; //johnfitz -- lit support via lordhavoc
                        return dl;
                }
        }
        dl = &cl_dlights[0];
        memset (dl, 0, sizeof(*dl));
        dl->key = key;
        dl->color[0] = dl->color[1] = dl->color[2] = 1; //johnfitz -- lit support via lordhavoc
        return dl;
}

void CL_DecayLights()
{
	f32 time = cl.time - cl.oldtime;
	dlight_t *dl = cl_dlights;
	for (s32 i = 0; i < MAX_DLIGHTS; i++, dl++) {
		if (dl->die < cl.time || !dl->radius)
			continue;
		dl->radius -= time * dl->decay;
		if (dl->radius < 0)
			dl->radius = 0;
	}
}

f32 CL_LerpPoint() // Determines the fraction between the last two messages
{ // that the objects should be put at.
	f32 f = cl.mtime[0] - cl.mtime[1];
	if (!f || cls.timedemo || (sv.active && !host_netinterval)) {
		cl.time = cl.mtime[0];
		return 1;
	}
	if (f > 0.1) { // dropped packet, or start of demo
		cl.mtime[1] = cl.mtime[0] - 0.1;
		f = 0.1;
	}
	f32 frac = (cl.time - cl.mtime[1]) / f;
	if (frac < 0) {
		if (frac < -0.01) cl.time = cl.mtime[1];
		frac = 0;
	} else if (frac > 1) {
		if (frac > 1.01) cl.time = cl.mtime[0];
		frac = 1;
	}
	if (cl_nolerp.value) return 1; //johnfitz -- better nolerp behavior
	return frac;
}

void CL_RelinkEntities()
{
	f32 frac = CL_LerpPoint(); // determine partial update time
	cl_numvisedicts = 0;
	for (s32 i = 0; i < 3; i++) // interpolate player info
		cl.velocity[i] = cl.mvelocity[1][i] +
		    frac * (cl.mvelocity[0][i] - cl.mvelocity[1][i]);
	if (cls.demoplayback) {
		for (s32 j = 0; j < 3; j++) { // interpolate the angles
			f32 d = cl.mviewangles[0][j] - cl.mviewangles[1][j];
			if (d > 180) d -= 360;
			else if (d < -180) d += 360;
			cl.viewangles[j] = cl.mviewangles[1][j] + frac * d;
		}
	}
	f32 bobjrotate = anglemod(100 * cl.time);
	// start on the entity after the world
	entity_t *ent = cl_entities + 1;
	for (s32 i = 1; i < cl.num_entities; i++, ent++) {
		if (!ent->model) { // empty slot
			if (ent->forcelink)
				R_RemoveEfrags(ent); // just became empty
			continue;
		} // if the object wasn't included in the last packet, remove it
		if (ent->msgtime != cl.mtime[0]) {
			ent->model = NULL;
			continue;
		}
		vec3_t oldorg;
		VectorCopy(ent->origin, oldorg);
		if (ent->forcelink) { // the entity was not updated in the last 
					// message so move to the final spot
			VectorCopy(ent->msg_origins[0], ent->origin);
			VectorCopy(ent->msg_angles[0], ent->angles);
		} else { // if the delta is large assume a teleport, don't lerp
			f32 f = frac;
			vec3_t delta;
			for (s32 j = 0; j < 3; j++) {
				delta[j] =
				    ent->msg_origins[0][j] -
				    ent->msg_origins[1][j];
				if (delta[j] > 100 || delta[j] < -100)
					f = 1;//assume teleportation, not motion
			}
			// interpolate the origin and angles
			for (s32 j = 0; j < 3; j++) {
				ent->origin[j] =
				    ent->msg_origins[1][j] + f * delta[j];
				f32 d = ent->msg_angles[0][j] -
				    ent->msg_angles[1][j];
				if (d > 180)
					d -= 360;
				else if (d < -180)
					d += 360;
				ent->angles[j] = ent->msg_angles[1][j] + f * d;
			}
		}
		if (ent->model->flags & EF_ROTATE) // rotate bin objects locally
			ent->angles[1] = bobjrotate;
		if (ent->effects & EF_BRIGHTFIELD)
			R_EntityParticles(ent);
		dlight_t *dl;
		if (ent->effects & EF_MUZZLEFLASH) {
			vec3_t fv, rv, uv;
			dl = CL_AllocDlight(i);
			VectorCopy(ent->origin, dl->origin);
			dl->origin[2] += 16;
			AngleVectors(ent->angles, fv, rv, uv);
			VectorMA(dl->origin, 18, fv, dl->origin);
			dl->radius = 200 + (rand() & 31);
			dl->minlight = 32;
			dl->die = cl.time + 0.1;
		}
		if (ent->effects & EF_BRIGHTLIGHT) {
			dl = CL_AllocDlight(i);
			VectorCopy(ent->origin, dl->origin);
			dl->origin[2] += 16;
			dl->radius = 400 + (rand() & 31);
			dl->die = cl.time + 0.001;
		}
		if (ent->effects & EF_DIMLIGHT) {
			dl = CL_AllocDlight(i);
			VectorCopy(ent->origin, dl->origin);
			dl->radius = 200 + (rand() & 31);
			dl->die = cl.time + 0.001;
		}
		if (ent->model->flags & EF_GIB)
			R_RocketTrail(oldorg, ent->origin, 2);
		else if (ent->model->flags & EF_ZOMGIB)
			R_RocketTrail(oldorg, ent->origin, 4);
		else if (ent->model->flags & EF_TRACER)
			R_RocketTrail(oldorg, ent->origin, 3);
		else if (ent->model->flags & EF_TRACER2)
			R_RocketTrail(oldorg, ent->origin, 5);
		else if (ent->model->flags & EF_ROCKET) {
			R_RocketTrail(oldorg, ent->origin, 0);
			dl = CL_AllocDlight(i);
			VectorCopy(ent->origin, dl->origin);
			dl->radius = 200;
			dl->die = cl.time + 0.01;
		} else if (ent->model->flags & EF_GRENADE)
			R_RocketTrail(oldorg, ent->origin, 1);
		else if (ent->model->flags & EF_TRACER3)
			R_RocketTrail(oldorg, ent->origin, 6);
		ent->forcelink = 0;
		if (i == cl.viewentity && !chase_active.value)
			continue;
		if (cl_numvisedicts < MAX_VISEDICTS) {
			cl_visedicts[cl_numvisedicts] = ent;
			cl_numvisedicts++;
		}
	}
}

s32 CL_ReadFromServer()
{ // Read all incoming data from the server
	cl.oldtime = cl.time;
	cl.time += host_frametime;
	s32 ret;
	do {
		ret = CL_GetMessage();
		if (ret == -1)
			Host_Error("CL_ReadFromServer: lost server connection");
		if (!ret)
			break;
		cl.last_received_message = realtime;
		CL_ParseServerMessage();
	} while (ret && cls.state == ca_connected);
	if (cl_shownet.value)
		Con_Printf("\n");
	CL_RelinkEntities();
	CL_UpdateTEnts(); // bring the links up to date
	return 0;
}

void CL_AccumulateCmd () // Spike: split from CL_SendCmd, to do clientside
{ // viewangle changes separately from outgoing packets.
    if (cls.signon == SIGNONS) {
        CL_AdjustAngles (); //basic keyboard looking
        IN_Move (&cl.pendingcmd); // accumulate movement from other devices
    }
}
void CL_SendCmd ()
{
	usercmd_t cmd;
	if (cls.state != ca_connected)
		return;
	if (cls.signon == SIGNONS) {
		CL_BaseMove (&cmd); // get basic movement from keyboard
		// allow mice or other external controllers to add to the move
		cmd.forwardmove += cl.pendingcmd.forwardmove;
		cmd.sidemove    += cl.pendingcmd.sidemove;
		cmd.upmove      += cl.pendingcmd.upmove;
		CL_SendMove (&cmd); // send the unreliable message
	}
	else
		CL_SendMove (NULL);
	memset(&cl.pendingcmd, 0, sizeof(cl.pendingcmd));
	if (cls.demoplayback) {
		SZ_Clear (&cls.message);
		return;
	}
	if (!cls.message.cursize) // send the reliable message
		return; // no message at all
	if (!NET_CanSendMessage (cls.netcon)) {
		Con_DPrintf ("CL_SendCmd: can't send\n");
		return;
	}
	if (NET_SendMessage (cls.netcon, &cls.message) == -1)
		Host_Error ("CL_SendCmd: lost server connection");
	SZ_Clear (&cls.message);
}

void CL_Init()
{
	SZ_Alloc(&cls.message, 1024);
	CL_InitInput();
	CL_InitTEnts();
	Cvar_RegisterVariable(&cl_name);
	Cvar_RegisterVariable(&cl_color);
	Cvar_RegisterVariable(&cl_upspeed);
	Cvar_RegisterVariable(&cl_forwardspeed);
	Cvar_RegisterVariable(&cl_backspeed);
	Cvar_RegisterVariable(&cl_sidespeed);
	Cvar_RegisterVariable(&cl_movespeedkey);
	Cvar_RegisterVariable(&cl_yawspeed);
	Cvar_RegisterVariable(&cl_pitchspeed);
	Cvar_RegisterVariable(&cl_anglespeedkey);
	Cvar_RegisterVariable(&cl_shownet);
	Cvar_RegisterVariable(&cl_nolerp);
	Cvar_RegisterVariable(&lookspring);
	Cvar_RegisterVariable(&lookstrafe);
	Cvar_RegisterVariable(&sensitivity);
	Cvar_RegisterVariable(&jlooksens);
	Cvar_RegisterVariable(&jmovesens);
	Cvar_RegisterVariable(&jdeadzone);
	Cvar_RegisterVariable(&jtriggerthresh);
	Cvar_RegisterVariable(&jmoveaxisx);
	Cvar_RegisterVariable(&jmoveaxisy);
	Cvar_RegisterVariable(&jlookaxisx);
	Cvar_RegisterVariable(&jlookaxisy);
	Cvar_RegisterVariable(&jtrigaxis1);
	Cvar_RegisterVariable(&jtrigaxis2);
	Cvar_RegisterVariable(&jenterbutton);
	Cvar_RegisterVariable(&jescapebutton);
	Cvar_RegisterVariable(&joysticknum);
	Cvar_SetCallback(&joysticknum, IN_InitJoystick);
	Cvar_RegisterVariable(&m_pitch);
	Cvar_RegisterVariable(&m_yaw);
	Cvar_RegisterVariable(&m_forward);
	Cvar_RegisterVariable(&m_side);
	Cvar_RegisterVariable (&cl_maxpitch); //johnfitz -- variable pitch clamp
	Cvar_RegisterVariable (&cl_minpitch); //johnfitz -- variable pitch clamp
	Cmd_AddCommand("entities", CL_PrintEntities_f);
	Cmd_AddCommand("disconnect", CL_Disconnect_f);
	Cmd_AddCommand("record", CL_Record_f);
	Cmd_AddCommand("stop", CL_Stop_f);
	Cmd_AddCommand("playdemo", CL_PlayDemo_f);
	Cmd_AddCommand("timedemo", CL_TimeDemo_f);
}
