// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
// cl_tent.c -- client side temporary entities

#include "quakedef.h"

static s32 num_temp_entities;
static sfx_t *cl_sfx_wizhit;
static sfx_t *cl_sfx_knighthit;
static sfx_t *cl_sfx_tink1;
static sfx_t *cl_sfx_ric1;
static sfx_t *cl_sfx_ric2;
static sfx_t *cl_sfx_ric3;
static sfx_t *cl_sfx_r_exp3;

void CL_InitTEnts()
{
	cl_sfx_wizhit = S_PrecacheSound("wizard/hit.wav");
	cl_sfx_knighthit = S_PrecacheSound("hknight/hit.wav");
	cl_sfx_tink1 = S_PrecacheSound("weapons/tink1.wav");
	cl_sfx_ric1 = S_PrecacheSound("weapons/ric1.wav");
	cl_sfx_ric2 = S_PrecacheSound("weapons/ric2.wav");
	cl_sfx_ric3 = S_PrecacheSound("weapons/ric3.wav");
	cl_sfx_r_exp3 = S_PrecacheSound("weapons/r_exp3.wav");
}

void CL_ParseBeam(model_t *m)
{
	s32 ent = MSG_ReadShort();
	vec3_t start, end;
	start[0] = MSG_ReadCoord(cl.protocolflags);
	start[1] = MSG_ReadCoord(cl.protocolflags);
	start[2] = MSG_ReadCoord(cl.protocolflags);
	end[0] = MSG_ReadCoord(cl.protocolflags);
	end[1] = MSG_ReadCoord(cl.protocolflags);
	end[2] = MSG_ReadCoord(cl.protocolflags);
	s32 i;
	beam_t *b; // override any beam with the same entity
	for(i = 0, b = cl_beams; i < MAX_BEAMS; i++, b++)
		if(b->entity == ent) {
			b->entity = ent;
			b->model = m;
			b->endtime = cl.time + 0.2;
			VectorCopy(start, b->start);
			VectorCopy(end, b->end);
			return;
		}
	for(i = 0, b = cl_beams; i < MAX_BEAMS; i++, b++) // find a free beam
		if(!b->model || b->endtime < cl.time) {
			b->entity = ent;
			b->model = m;
			b->endtime = cl.time + 0.2;
			VectorCopy(start, b->start);
			VectorCopy(end, b->end);
			return;
		}
}

void CL_ParseTEnt()
{
	dlight_t *dl;
	vec3_t pos;
	s32 type = MSG_ReadByte();
	switch(type) {
	case TE_WIZSPIKE: // spike hitting wall
		pos[0] = MSG_ReadCoord(cl.protocolflags);
		pos[1] = MSG_ReadCoord(cl.protocolflags);
		pos[2] = MSG_ReadCoord(cl.protocolflags);
		R_RunParticleEffect(pos, vec3_origin, 20, 30);
		S_StartSound(-1, 0, cl_sfx_wizhit, pos, 1, 1);
		break;
	case TE_KNIGHTSPIKE: // spike hitting wall
		pos[0] = MSG_ReadCoord(cl.protocolflags);
		pos[1] = MSG_ReadCoord(cl.protocolflags);
		pos[2] = MSG_ReadCoord(cl.protocolflags);
		R_RunParticleEffect(pos, vec3_origin, 226, 20);
		S_StartSound(-1, 0, cl_sfx_knighthit, pos, 1, 1);
		break;
	case TE_SPIKE: // spike hitting wall
		pos[0] = MSG_ReadCoord(cl.protocolflags);
		pos[1] = MSG_ReadCoord(cl.protocolflags);
		pos[2] = MSG_ReadCoord(cl.protocolflags);
		R_RunParticleEffect(pos, vec3_origin, 0, 10);
		if(rand() % 5)
			S_StartSound(-1, 0, cl_sfx_tink1, pos, 1, 1);
		else{
			s32 rnd = rand() & 3;
			if(rnd==1)S_StartSound(-1,0,cl_sfx_ric1,pos,1,1);
			else if(rnd==2)S_StartSound(-1,0,cl_sfx_ric2,pos,1,1);
			else S_StartSound(-1,0,cl_sfx_ric3,pos,1,1);
		}
		break;
	case TE_SUPERSPIKE: // super spike hitting wall
		pos[0] = MSG_ReadCoord(cl.protocolflags);
		pos[1] = MSG_ReadCoord(cl.protocolflags);
		pos[2] = MSG_ReadCoord(cl.protocolflags);
		R_RunParticleEffect(pos, vec3_origin, 0, 20);
		if(rand() % 5)
			S_StartSound(-1, 0, cl_sfx_tink1, pos, 1, 1);
		else{
			s32 rnd = rand() & 3;
			if(rnd==1)S_StartSound(-1,0,cl_sfx_ric1,pos,1,1);
			else if(rnd==2)S_StartSound(-1,0,cl_sfx_ric2,pos,1,1);
			else S_StartSound(-1,0,cl_sfx_ric3,pos,1,1);
		}
		break;
	case TE_GUNSHOT: // bullet hitting wall
		pos[0] = MSG_ReadCoord(cl.protocolflags);
		pos[1] = MSG_ReadCoord(cl.protocolflags);
		pos[2] = MSG_ReadCoord(cl.protocolflags);
		R_RunParticleEffect(pos, vec3_origin, 0, 20);
		break;
	case TE_EXPLOSION: // rocket explosion
		pos[0] = MSG_ReadCoord(cl.protocolflags);
		pos[1] = MSG_ReadCoord(cl.protocolflags);
		pos[2] = MSG_ReadCoord(cl.protocolflags);
		R_ParticleExplosion(pos);
		dl = CL_AllocDlight(0);
		VectorCopy(pos, dl->origin);
		dl->radius = 350;
		dl->die = cl.time + 0.5;
		dl->decay = 300;
		S_StartSound(-1, 0, cl_sfx_r_exp3, pos, 1, 1);
		break;
	case TE_TAREXPLOSION: // tarbaby explosion
		pos[0] = MSG_ReadCoord(cl.protocolflags);
		pos[1] = MSG_ReadCoord(cl.protocolflags);
		pos[2] = MSG_ReadCoord(cl.protocolflags);
		R_BlobExplosion(pos);
		S_StartSound(-1, 0, cl_sfx_r_exp3, pos, 1, 1);
		break;
	case TE_LIGHTNING1: // lightning bolts
		CL_ParseBeam(Mod_ForName("progs/bolt.mdl", 1));
		break;
	case TE_LIGHTNING2: // lightning bolts
		CL_ParseBeam(Mod_ForName("progs/bolt2.mdl", 1));
		break;
	case TE_LIGHTNING3: // lightning bolts
		CL_ParseBeam(Mod_ForName("progs/bolt3.mdl", 1));
		break;
	case TE_BEAM: // grappling hook beam
		CL_ParseBeam(Mod_ForName("progs/beam.mdl", 1));
		break;
	case TE_LAVASPLASH:
		pos[0] = MSG_ReadCoord(cl.protocolflags);
		pos[1] = MSG_ReadCoord(cl.protocolflags);
		pos[2] = MSG_ReadCoord(cl.protocolflags);
		R_LavaSplash(pos);
		break;
	case TE_TELEPORT:
		pos[0] = MSG_ReadCoord(cl.protocolflags);
		pos[1] = MSG_ReadCoord(cl.protocolflags);
		pos[2] = MSG_ReadCoord(cl.protocolflags);
		R_TeleportSplash(pos);
		break;
	case TE_EXPLOSION2: // color mapped explosion
		pos[0] = MSG_ReadCoord(cl.protocolflags);
		pos[1] = MSG_ReadCoord(cl.protocolflags);
		pos[2] = MSG_ReadCoord(cl.protocolflags);
		s32 colorStart = MSG_ReadByte();
		s32 colorLength = MSG_ReadByte();
		R_ParticleExplosion2(pos, colorStart, colorLength);
		dl = CL_AllocDlight(0);
		VectorCopy(pos, dl->origin);
		dl->radius = 350;
		dl->die = cl.time + 0.5;
		dl->decay = 300;
		S_StartSound(-1, 0, cl_sfx_r_exp3, pos, 1, 1);
		break;
	default:
		Sys_Error("CL_ParseTEnt: bad type");
	}
}

entity_t *CL_NewTempEntity()
{
	if(cl_numvisedicts == MAX_VISEDICTS) return NULL;
	if(num_temp_entities == MAX_TEMP_ENTITIES) return NULL;
	entity_t *ent = &cl_temp_entities[num_temp_entities];
	memset(ent, 0, sizeof(*ent));
	num_temp_entities++;
	cl_visedicts[cl_numvisedicts] = ent;
	cl_numvisedicts++;
	ent->scale = ENTSCALE_DEFAULT;
	ent->colormap = CURWORLDCMAP;
	return ent;
}

void CL_UpdateTEnts()
{
	num_temp_entities = 0;
	srand((s32)(cl.time * 1000)); //johnfitz -- freeze beams when paused
	s32 i;
	beam_t *b;
	for(i = 0, b = cl_beams; i < MAX_BEAMS; i++, b++) { // update lightning
		if(!b->model || b->endtime < cl.time) continue;
		// if coming from the player, update the start position
		if(b->entity == cl.viewentity)
			VectorCopy(cl_entities[cl.viewentity].origin, b->start);
		// calculate pitch and yaw
		vec3_t dist, org;
		VectorSubtract(b->end, b->start, dist);
		f32 yaw, pitch;
		if(dist[1] == 0 && dist[0] == 0) {
			yaw = 0;
			pitch = dist[2] > 0 ? 90 : 270;
		} else {
			yaw = (s32) (atan2(dist[1], dist[0]) * 180 / M_PI);
			if(yaw < 0) yaw += 360;
			f32 forward = sqrt(dist[0]*dist[0] + dist[1]*dist[1]);
			pitch = (s32) (atan2(dist[2], forward) * 180 / M_PI);
			if(pitch < 0) pitch += 360;
		}
		VectorCopy(b->start, org); // add new entities for the lightning
		f32 d = VectorNormalize(dist);
		while(d > 0) {
			entity_t *ent = CL_NewTempEntity();
			if(!ent) return;
			VectorCopy(org, ent->origin);
			ent->model = b->model;
			ent->angles[0] = pitch;
			ent->angles[1] = yaw;
			ent->angles[2] = rand()%360;
			for(s32 j = 0; j < 3; j++)
				org[j] += dist[j]*30;
			d -= 30;
		}
	}
}
