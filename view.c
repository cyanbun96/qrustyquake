// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// view.c -- player eye positioning
#include "quakedef.h"
// The view is allowed to move slightly from it's 1 position for bobbing,
// but if it exceeds 8 pixels linear distance(spherical, not box), the list of
// entities sent from the server may not include everything in the pvs,
// especially when crossing a water boudnary.

static u8 gammatable[256];
static f32 v_dmg_time, v_dmg_roll, v_dmg_pitch;
static vec3_t forward, right, up;
static cshift_t cshift_empty = { { 130, 80, 50 }, 0 }; // Palette flashes 
static cshift_t cshift_water = { { 130, 80, 50 }, 128 };
static cshift_t cshift_slime = { { 0, 25, 5 }, 150 };
static cshift_t cshift_lava = { { 255, 80, 0 }, 150 };

f32 V_CalcRoll(vec3_t angles, vec3_t velocity)
{ // Used by view and sv_user
	AngleVectors(angles, forward, right, up);
	f32 side = DotProduct(velocity, right);
	f32 sign = side < 0 ? -1 : 1;
	side = fabs(side);
	f32 value = cl_rollangle.value;
	side=side<cl_rollspeed.value?side*value/cl_rollspeed.value:value;
	return side * sign;
}

f32 V_CalcBob()
{
	f32 cycle=cl.time-(s32)(cl.time/cl_bobcycle.value)*cl_bobcycle.value;
	cycle /= cl_bobcycle.value;
	if(cycle < cl_bobup.value) cycle = M_PI * cycle / cl_bobup.value;
	else cycle=M_PI+M_PI*(cycle-cl_bobup.value)/(1.0-cl_bobup.value);
	// bob is proportional to velocity in the xy plane
	// (don't count Z, or jumping messes it up)
	f32 bob = sqrt(cl.velocity[0] * cl.velocity[0] + cl.velocity[1]
			* cl.velocity[1]) * cl_bob.value;
	bob = bob * 0.3 + bob * 0.7 * sin(cycle);
	if(bob > 4) bob = 4;
	else if(bob < -7) bob = -7;
	return bob;
}

void V_StartPitchDrift()
{
	if(cl.laststop == cl.time)
		return; // something else is keeping it from drifting
	if(cl.nodrift || !cl.pitchvel){
		cl.pitchvel = v_centerspeed.value;
		cl.nodrift = 0;
		cl.driftmove = 0;
	}
}

void V_StopPitchDrift()
{
	cl.laststop = cl.time;
	cl.nodrift = 1;
	cl.pitchvel = 0;
}

// Moves the client pitch angle towards cl.idealpitch sent by the server.
// If the user is adjusting pitch manually, either with lookup/lookdown,
// mlook and mouse, or klook and keyboard, pitch drifting is constantly stopped.
// Drifting is enabled when the center view key is hit, mlook is released and
// lookspring is non 0, or when 
void V_DriftPitch()
{
	if(noclip_anglehack || !cl.onground || cls.demoplayback){
		cl.driftmove = 0;
		cl.pitchvel = 0;
		return;
	}
	if(cl.nodrift){ // don't count small mouse motion
		if(fabs(cl.cmd.forwardmove) < cl_forwardspeed.value)
			cl.driftmove = 0;
		else cl.driftmove += host_frametime;
		if(cl.driftmove > v_centermove.value)
			V_StartPitchDrift();
		return;
	}
	f32 delta = cl.idealpitch - cl.viewangles[PITCH];
	if(!delta){
		cl.pitchvel = 0;
		return;
	}
	f32 move = host_frametime * cl.pitchvel;
	cl.pitchvel += host_frametime * v_centerspeed.value;
	if(delta > 0){
		if(move > delta){
			cl.pitchvel = 0;
			move = delta;
		}
		cl.viewangles[PITCH] += move;
	} else if(delta < 0){
		if(move > -delta){
			cl.pitchvel = 0;
			move = -delta;
		}
		cl.viewangles[PITCH] -= move;
	}
}

void BuildGammaTable(f32 g)
{
	if(g == 1.0){
		for(s32 i = 0; i < 256; i++)
			gammatable[i] = i;
		return;
	}
	for(s32 i = 0; i < 256; i++){
		s32 inf = 255 * pow((i + 0.5) / 255.5, g) + 0.5;
		if(inf < 0) inf = 0;
		if(inf > 255) inf = 255;
		gammatable[i] = inf;
	}
}

bool V_CheckGamma()
{
	static f32 oldgammavalue;
	if(v_gamma.value == oldgammavalue) return 0;
	oldgammavalue = v_gamma.value;
	BuildGammaTable(v_gamma.value);
	vid.recalc_refdef = 1; // force a surface cache flush
	return 1;
}

void V_ParseDamage()
{
	vec3_t from;
	vec3_t forward, right, up;
	s32 armor = MSG_ReadByte();
	s32 blood = MSG_ReadByte();
	for(s32 i = 0; i < 3; i++)
		from[i] = MSG_ReadCoord(cl.protocolflags);
	f32 count = blood*0.5 + armor*0.5;
	if(count < 10) count = 10;
	cl.faceanimtime = cl.time + 0.2; // but sbar face into pain frame
	cl.cshifts[CSHIFT_DAMAGE].percent += 3*count;
	if(cl.cshifts[CSHIFT_DAMAGE].percent < 0)
		cl.cshifts[CSHIFT_DAMAGE].percent = 0;
	if(cl.cshifts[CSHIFT_DAMAGE].percent > 150)
		cl.cshifts[CSHIFT_DAMAGE].percent = 150;
	if(armor > blood) {
		cl.cshifts[CSHIFT_DAMAGE].destcolor[0] = 200;
		cl.cshifts[CSHIFT_DAMAGE].destcolor[1] = 100;
		cl.cshifts[CSHIFT_DAMAGE].destcolor[2] = 100;
	} else if(armor) {
		cl.cshifts[CSHIFT_DAMAGE].destcolor[0] = 220;
		cl.cshifts[CSHIFT_DAMAGE].destcolor[1] = 50;
		cl.cshifts[CSHIFT_DAMAGE].destcolor[2] = 50;
	} else {
		cl.cshifts[CSHIFT_DAMAGE].destcolor[0] = 255;
		cl.cshifts[CSHIFT_DAMAGE].destcolor[1] = 0;
		cl.cshifts[CSHIFT_DAMAGE].destcolor[2] = 0;
	}
	entity_t *ent = &cl_entities[cl.viewentity];//calculate view angle kicks
	VectorSubtract(from, ent->origin, from);
	VectorNormalize(from);
	AngleVectors(ent->angles, forward, right, up);
	f32 side = DotProduct(from, right);
	v_dmg_roll = count*side*v_kickroll.value;
	side = DotProduct(from, forward);
	v_dmg_pitch = count*side*v_kickpitch.value;
	v_dmg_time = v_kicktime.value;
}

void V_cshift_f()
{
	cshift_empty.destcolor[0] = atoi(Cmd_Argv(1));
	cshift_empty.destcolor[1] = atoi(Cmd_Argv(2));
	cshift_empty.destcolor[2] = atoi(Cmd_Argv(3));
	cshift_empty.percent = atoi(Cmd_Argv(4));
}

void V_BonusFlash_f()
{ // When you run over an item, the server sends this command
	cl.cshifts[CSHIFT_BONUS].destcolor[0] = 215;
	cl.cshifts[CSHIFT_BONUS].destcolor[1] = 186;
	cl.cshifts[CSHIFT_BONUS].destcolor[2] = 69;
	cl.cshifts[CSHIFT_BONUS].percent = 50;
}


void V_SetContentsColor(s32 contents)
{ // Underwater, lava, etc each has a color shift
	switch(contents){
	case CONTENTS_EMPTY:
	case CONTENTS_SOLID: cl.cshifts[CSHIFT_CONTENTS] = cshift_empty; break;
	case CONTENTS_LAVA:  cl.cshifts[CSHIFT_CONTENTS] = cshift_lava; break;
	case CONTENTS_SLIME: cl.cshifts[CSHIFT_CONTENTS] = cshift_slime; break;
	default:             cl.cshifts[CSHIFT_CONTENTS] = cshift_water; }
}

void V_CalcPowerupCshift()
{
	if(cl.items & IT_QUAD){
		cl.cshifts[CSHIFT_POWERUP].destcolor[0] = 0;
		cl.cshifts[CSHIFT_POWERUP].destcolor[1] = 0;
		cl.cshifts[CSHIFT_POWERUP].destcolor[2] = 255;
		cl.cshifts[CSHIFT_POWERUP].percent = 30;
	} else if(cl.items & IT_SUIT){
		cl.cshifts[CSHIFT_POWERUP].destcolor[0] = 0;
		cl.cshifts[CSHIFT_POWERUP].destcolor[1] = 255;
		cl.cshifts[CSHIFT_POWERUP].destcolor[2] = 0;
		cl.cshifts[CSHIFT_POWERUP].percent = 20;
	} else if(cl.items & IT_INVISIBILITY){
		cl.cshifts[CSHIFT_POWERUP].destcolor[0] = 100;
		cl.cshifts[CSHIFT_POWERUP].destcolor[1] = 100;
		cl.cshifts[CSHIFT_POWERUP].destcolor[2] = 100;
		cl.cshifts[CSHIFT_POWERUP].percent = 100;
	} else if(cl.items & IT_INVULNERABILITY){
		cl.cshifts[CSHIFT_POWERUP].destcolor[0] = 255;
		cl.cshifts[CSHIFT_POWERUP].destcolor[1] = 255;
		cl.cshifts[CSHIFT_POWERUP].destcolor[2] = 0;
		cl.cshifts[CSHIFT_POWERUP].percent = 30;
	} else cl.cshifts[CSHIFT_POWERUP].percent = 0;
}

void V_UpdatePalette()
{
	V_CalcPowerupCshift();
	f32 frametime = fabs(cl.time - cl.oldtime);
	bool new = 0;
	for(s32 i = 0; i < NUM_CSHIFTS; i++){
		if(cl.cshifts[i].percent != cl.prev_cshifts[i].percent){
			new = 1;
			cl.prev_cshifts[i].percent = cl.cshifts[i].percent;
		}
		for(s32 j = 0; j < 3; j++)
			if(cl.cshifts[i].destcolor[j] !=
					cl.prev_cshifts[i].destcolor[j]){
				new = 1;
				cl.prev_cshifts[i].destcolor[j] =
					cl.cshifts[i].destcolor[j];
			}
	}
	// drop the damage value
	cl.cshifts[CSHIFT_DAMAGE].percent -= frametime * 150;
	if(cl.cshifts[CSHIFT_DAMAGE].percent <= 0)
		cl.cshifts[CSHIFT_DAMAGE].percent = 0;
	// drop the bonus value
	cl.cshifts[CSHIFT_BONUS].percent -= frametime * 100;
	if(cl.cshifts[CSHIFT_BONUS].percent <= 0)
		cl.cshifts[CSHIFT_BONUS].percent = 0;
	bool force = V_CheckGamma();
	if(!new && !force && !refresh_palette)
		return;
	refresh_palette--;
	u8 pal[768];
	u8 *basepal = CURWORLDPAL;
	u8 *newpal = pal;
	for(s32 i = 0; i < 256; i++){
		s32 r = basepal[0];
		s32 g = basepal[1];
		s32 b = basepal[2];
		basepal += 3;
		for(s32 j = 0; j < NUM_CSHIFTS; j++){
			r += (cl.cshifts[j].percent *
					(cl.cshifts[j].destcolor[0] - r)) >> 8;
			g += (cl.cshifts[j].percent *
					(cl.cshifts[j].destcolor[1] - g)) >> 8;
			b += (cl.cshifts[j].percent *
					(cl.cshifts[j].destcolor[2] - b)) >> 8;
		}
		newpal[0] = gammatable[r];
		newpal[1] = gammatable[g];
		newpal[2] = gammatable[b];
		newpal += 3;
	}
	VID_SetPalette(pal, screen);
}

f32 angledelta(f32 a) // View rendering 
{
	a = anglemod(a);
	if(a > 180) a -= 360;
	return a;
}

void CalcGunAngle()
{
	static f32 oldyaw = 0;
	static f32 oldpitch = 0;
	f32 yaw = r_refdef.viewangles[YAW];
	f32 pitch = -r_refdef.viewangles[PITCH];
	yaw = angledelta(yaw - r_refdef.viewangles[YAW]) * 0.4;
	if(yaw > 10) yaw = 10;
	if(yaw < -10) yaw = -10;
	pitch = angledelta(-pitch - r_refdef.viewangles[PITCH]) * 0.4;
	if(pitch > 10) pitch = 10;
	if(pitch < -10) pitch = -10;
	f32 move = host_frametime * 20;
	if(yaw > oldyaw){
		if(oldyaw + move < yaw) yaw = oldyaw + move;
	} else if(oldyaw - move > yaw) yaw = oldyaw - move;
	if(pitch > oldpitch){
		if(oldpitch + move < pitch) pitch = oldpitch + move;
	} else if(oldpitch - move > pitch) pitch = oldpitch - move;
	oldyaw = yaw;
	oldpitch = pitch;
	cl.viewent.angles[YAW] = r_refdef.viewangles[YAW] + yaw;
	cl.viewent.angles[PITCH] = -(r_refdef.viewangles[PITCH] + pitch);
	cl.viewent.angles[ROLL] -= v_idlescale.value * sin(cl.time *
			v_iroll_cycle.value) * v_iroll_level.value;
	cl.viewent.angles[PITCH] -= v_idlescale.value * sin(cl.time *
			v_ipitch_cycle.value) * v_ipitch_level.value;
	cl.viewent.angles[YAW] -= v_idlescale.value * sin(cl.time *
			v_iyaw_cycle.value) * v_iyaw_level.value;
}

void V_BoundOffsets()
{
	entity_t *ent = &cl_entities[cl.viewentity];
	// absolutely bound refresh reletive to entity clipping hull
	// so the view can never be inside a solid wall
	if(r_refdef.vieworg[0] < ent->origin[0] - 14)
		r_refdef.vieworg[0] = ent->origin[0] - 14;
	else if(r_refdef.vieworg[0] > ent->origin[0] + 14)
		r_refdef.vieworg[0] = ent->origin[0] + 14;
	if(r_refdef.vieworg[1] < ent->origin[1] - 14)
		r_refdef.vieworg[1] = ent->origin[1] - 14;
	else if(r_refdef.vieworg[1] > ent->origin[1] + 14)
		r_refdef.vieworg[1] = ent->origin[1] + 14;
	if(r_refdef.vieworg[2] < ent->origin[2] - 22)
		r_refdef.vieworg[2] = ent->origin[2] - 22;
	else if(r_refdef.vieworg[2] > ent->origin[2] + 30)
		r_refdef.vieworg[2] = ent->origin[2] + 30;
}

void V_AddIdle()
{ // Idle swaying
	r_refdef.viewangles[ROLL] += v_idlescale.value
		* sin(cl.time * v_iroll_cycle.value) * v_iroll_level.value;
	r_refdef.viewangles[PITCH] += v_idlescale.value
		* sin(cl.time * v_ipitch_cycle.value) * v_ipitch_level.value;
	r_refdef.viewangles[YAW] += v_idlescale.value
		* sin(cl.time * v_iyaw_cycle.value) * v_iyaw_level.value;
}


void V_CalcViewRoll()
{ // Roll is induced by movement and damage
	f32 side = V_CalcRoll(cl_entities[cl.viewentity].angles, cl.velocity);
	r_refdef.viewangles[ROLL] += side;
	if(v_dmg_time > 0){
	     r_refdef.viewangles[ROLL]+=v_dmg_time/v_kicktime.value*v_dmg_roll;
	    r_refdef.viewangles[PITCH]+=v_dmg_time/v_kicktime.value*v_dmg_pitch;
		v_dmg_time -= host_frametime;
	}
	if(cl.stats[STAT_HEALTH] <= 0)
		r_refdef.viewangles[ROLL] = 80; // dead view angle
}

void V_CalcIntermissionRefdef()
{
	entity_t *ent = &cl_entities[cl.viewentity]; // player model
	entity_t *view = &cl.viewent; // weapon model
	VectorCopy(ent->origin, r_refdef.vieworg);
	VectorCopy(ent->angles, r_refdef.viewangles);
	view->model = NULL;
	f32 old = v_idlescale.value; // allways idle in intermission
	v_idlescale.value = 1;
	V_AddIdle();
	v_idlescale.value = old;
}

void V_CalcRefdef()
{
	static f32 oldz = 0;
	V_DriftPitch();
	entity_t *ent = &cl_entities[cl.viewentity];
	entity_t *view = &cl.viewent;
	// transform the view offset by the model's matrix to get the offset
	// from model origin for the view
	vec3_t angles;
	ent->angles[YAW]=cl.viewangles[YAW];//the model should face the view dir
	ent->angles[PITCH]=-cl.viewangles[PITCH];
	f32 bob = V_CalcBob();
	VectorCopy(ent->origin, r_refdef.vieworg); // refresh position
	r_refdef.vieworg[2] += cl.viewheight + bob;
	// never let it sit exactly on a node line, because a water plane can
	// dissapear when viewed with the eye exactly on it. the server protocol
	// only specifies to 1/16 pixel, so add 1/32 in each axis
	r_refdef.vieworg[0] += 1.0 / 32;
	r_refdef.vieworg[1] += 1.0 / 32;
	r_refdef.vieworg[2] += 1.0 / 32;
	VectorCopy(cl.viewangles, r_refdef.viewangles);
	V_CalcViewRoll();
	V_AddIdle();
	// offsets
	angles[PITCH] = -ent->angles[PITCH]; // because entity pitches are
	angles[YAW] = ent->angles[YAW]; // actually backward
	angles[ROLL] = ent->angles[ROLL];
	vec3_t forward, right, up;
	AngleVectors(angles, forward, right, up);
	for(s32 i = 0; i < 3; i++)
		r_refdef.vieworg[i] += scr_ofsx.value * forward[i]
			+ scr_ofsy.value * right[i] + scr_ofsz.value * up[i];
	V_BoundOffsets();
	VectorCopy(cl.viewangles, view->angles); // set up gun position
	CalcGunAngle();
	VectorCopy(ent->origin, view->origin);
	view->origin[2] += cl.viewheight;
	for(s32 i = 0; i < 3; i++)
		view->origin[i] += forward[i] * bob * 0.4;
	view->origin[2] += bob;
	// fudge position around to keep amount of weapon visible
	// roughly equal with different FOV
	if(scr_viewsize.value == 110) view->origin[2] += 1;
	else if(scr_viewsize.value == 100) view->origin[2] += 2;
	else if(scr_viewsize.value == 90) view->origin[2] += 1;
	else if(scr_viewsize.value == 80) view->origin[2] += 0.5;
	view->model = cl.model_precache[cl.stats[STAT_WEAPON]];
	view->frame = cl.stats[STAT_WEAPONFRAME];
	view->colormap = CURWORLDCMAP;
	// set up the refresh position
	VectorAdd(r_refdef.viewangles, cl.punchangle, r_refdef.viewangles);
	// smooth out stair step ups
	if(cl.onground && ent->origin[2] - oldz > 0){
		f32 steptime;
		steptime = cl.time - cl.oldtime;
		if(steptime < 0) steptime = 0;
		oldz += steptime * 80;
		if(oldz > ent->origin[2]) oldz = ent->origin[2];
		if(ent->origin[2] - oldz > 12) oldz = ent->origin[2] - 12;
		r_refdef.vieworg[2] += oldz - ent->origin[2];
		view->origin[2] += oldz - ent->origin[2];
	} else oldz = ent->origin[2];
	if(chase_active.value) Chase_Update();
}

void V_AllocLedges()
{
	r_ledges_size *= 1.2;
	r_ledges_size += NUMSTACKEDGES;
	Con_DPrintf("Reallocing ledges %d\n", r_ledges_size);
	s32 alloc_sz = (r_ledges_size+((CACHE_SIZE-1)/sizeof(edge_t))+1)
				* sizeof(edge_t);
	ledges = realloc(ledges, alloc_sz);

}

void V_AllocLsurfs()
{
	r_lsurfs_size *= 1.2;
	r_lsurfs_size += NUMSTACKSURFACES;
	Con_DPrintf("Reallocing lsurfs %d\n", r_lsurfs_size);
	s32 alloc_sz = (r_lsurfs_size+((CACHE_SIZE-1)/sizeof(surf_t))+1)
				* sizeof(surf_t);
	lsurfs = realloc(lsurfs, alloc_sz);
}

void V_RenderView() // The player's clipping box goes from(-16 -16 -24) to(16
{// 16 32) from the entity origin so any view position inside that will be valid
	if(con_forcedup) return;
	if(cl.maxclients > 1){ // don't allow cheats in multiplayer
		Cvar_Set("scr_ofsx", "0");
		Cvar_Set("scr_ofsy", "0");
		Cvar_Set("scr_ofsz", "0");
	}
	if(cl.intermission) V_CalcIntermissionRefdef();
	else if(!cl. paused) V_CalcRefdef();
	if(!cl_entities[0].model || !cl.worldmodel)
		Sys_Error("R_RenderView: NULL worldmodel");
	R_PushDlights();
	R_RenderView();
	if(r_dspeeds.value) R_PrintDSpeeds();
	if(r_reportsurfout.value && r_outofsurfaces)
		Con_Printf("short %d surfaces\n", r_outofsurfaces);
	if(r_reportedgeout.value && r_outofedges)
		Con_Printf("short roughly %d edges\n", r_outofedges * 2 / 3);
	if(r_outofedges) V_AllocLedges();
	if(r_outofsurfaces) V_AllocLsurfs();
	if(!crosshair.value) return;
	drawlayer = lyr_crosshair.value;
	Draw_CharacterScaled(scr_vrect.x + scr_vrect.width / 2 - uiscale*4
			+ cl_crossx.value, scr_vrect.y + scr_vrect.height / 2
			+ cl_crossy.value, cl_crosschar.value, uiscale);
	drawlayer = lyr_main.value;
}

void V_Init()
{
	Cmd_AddCommand("v_cshift", V_cshift_f);
	Cmd_AddCommand("bf", V_BonusFlash_f);
	Cmd_AddCommand("centerview", V_StartPitchDrift);
	Cvar_RegisterVariable(&v_centermove);
	Cvar_RegisterVariable(&v_centerspeed);
	Cvar_RegisterVariable(&v_iyaw_cycle);
	Cvar_RegisterVariable(&v_iroll_cycle);
	Cvar_RegisterVariable(&v_ipitch_cycle);
	Cvar_RegisterVariable(&v_iyaw_level);
	Cvar_RegisterVariable(&v_iroll_level);
	Cvar_RegisterVariable(&v_ipitch_level);
	Cvar_RegisterVariable(&v_idlescale);
	Cvar_RegisterVariable(&crosshair);
	Cvar_RegisterVariable(&cl_crossx);
	Cvar_RegisterVariable(&cl_crossy);
	Cvar_RegisterVariable(&cl_crosschar);
	Cvar_RegisterVariable(&gl_cshiftpercent);
	Cvar_RegisterVariable(&scr_ofsx);
	Cvar_RegisterVariable(&scr_ofsy);
	Cvar_RegisterVariable(&scr_ofsz);
	Cvar_RegisterVariable(&cl_rollspeed);
	Cvar_RegisterVariable(&cl_rollangle);
	Cvar_RegisterVariable(&cl_bob);
	Cvar_RegisterVariable(&cl_bobcycle);
	Cvar_RegisterVariable(&cl_bobup);
	Cvar_RegisterVariable(&v_kicktime);
	Cvar_RegisterVariable(&v_kickroll);
	Cvar_RegisterVariable(&v_kickpitch);
	BuildGammaTable(1.0); // no gamma yet
	Cvar_RegisterVariable(&v_gamma);
	V_AllocLedges();
	V_AllocLsurfs();
}
