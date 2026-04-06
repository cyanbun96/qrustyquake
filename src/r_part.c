// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

static s32 ramp1[8] = { 0x6f, 0x6d, 0x6b, 0x69, 0x67, 0x65, 0x63, 0x61 };
static s32 ramp2[8] = { 0x6f, 0x6e, 0x6d, 0x6c, 0x6b, 0x6a, 0x68, 0x66 };
static s32 ramp3[8] = { 0x6d, 0x6b, 6, 5, 4, 3 };
static f32 beamlength = 16;
static particle_t *active_particles, *free_particles;
static particle_t *particles;
static s32 r_numparticles;
static vec3_t avelocities[NUMVERTEXNORMALS];

void R_InitParticles()
{
	s32 i = COM_CheckParm("-particles");
	if (i) {
		r_numparticles = (s32)(Q_atoi(com_argv[i + 1]));
		if (r_numparticles < ABSOLUTE_MIN_PARTICLES)
			r_numparticles = ABSOLUTE_MIN_PARTICLES;
	} else
		r_numparticles = MAX_PARTICLES;
	particles = (particle_t *)
		Hunk_AllocName(r_numparticles*sizeof(particle_t), "particles");
}

void R_EntityParticles(entity_t *ent)
{
	f32 dist = 64;
	if (!avelocities[0][0]) {
		for (s32 i = 0; i < NUMVERTEXNORMALS; i++) {
			avelocities[i][0] = (rand() & 255) * 0.01;
			avelocities[i][1] = (rand() & 255) * 0.01;
			avelocities[i][2] = (rand() & 255) * 0.01;
		}
	}
	for (s32 i = 0; i < NUMVERTEXNORMALS; i++) {
		f32 angle = cl.time * avelocities[i][0];
		f32 sy = sin(angle);
		f32 cy = cos(angle);
		angle = cl.time * avelocities[i][1];
		f32 sp = sin(angle);
		f32 cp = cos(angle);
		angle = cl.time * avelocities[i][2];
		vec3_t forward;
		forward[0] = cp * cy;
		forward[1] = cp * sy;
		forward[2] = -sp;
		if (!free_particles)
			return;
		particle_t *p = free_particles;
		free_particles = p->next;
		p->next = active_particles;
		active_particles = p;
		p->die = cl.time + 0.01;
		p->color = 0x6f;
		p->type = pt_explode;
		p->org[0] = ent->origin[0] + r_avertexnormals[i][0] * dist
			+ forward[0] * beamlength;
		p->org[1] = ent->origin[1] + r_avertexnormals[i][1] * dist
			+ forward[1] * beamlength;
		p->org[2] = ent->origin[2] + r_avertexnormals[i][2] * dist
			+ forward[2] * beamlength;
	}
}

void R_ClearParticles()
{
	free_particles = &particles[0];
	active_particles = NULL;
	for (s32 i = 0; i < r_numparticles; i++)
		particles[i].next = &particles[i + 1];
	particles[r_numparticles - 1].next = NULL;
}

void R_ReadPointFile_f()
{
	s8 name[MAX_OSPATH];
	sprintf(name, "maps/%s.pts", sv.name);
	FILE *f;
	COM_FOpenFile(name, &f, NULL);
	if (!f) {
		Con_Printf("couldn't open %s\n", name);
		return;
	}
	Con_Printf("Reading %s...\n", name);
	s32 c = 0;
	for (;;) {
		vec3_t org;
		s32 r = fscanf(f, "%f %f %f\n", &org[0], &org[1], &org[2]);
		if (r != 3)
			break;
		c++;
		if (!free_particles) {
			Con_Printf("Not enough free particles\n");
			break;
		}
		particle_t *p = free_particles;
		free_particles = p->next;
		p->next = active_particles;
		active_particles = p;
		p->die = 99999;
		p->color = (-c) & 15;
		p->type = pt_static;
		VectorCopy(vec3_origin, p->vel);
		VectorCopy(org, p->org);
	}
	fclose(f);
	Con_Printf("%i points read\n", c);
}

void R_ParseParticleEffect ()
{
        vec3_t org, dir;
        for (s32 i=0 ; i<3 ; i++)
                org[i] = MSG_ReadCoord (cl.protocolflags);
        for (s32 i=0 ; i<3 ; i++)
                dir[i] = MSG_ReadChar () * (1.0/16);
        s32 msgcount = MSG_ReadByte ();
        s32 color = MSG_ReadByte ();
        s32 count = msgcount == 255 ? 1024 : msgcount;
        R_RunParticleEffect (org, dir, color, count);
}

// CyanBun96: TODO a lot of shared code between these. consolidate?
void R_ParticleExplosion(vec3_t org)
{
	for (s32 i = 0; i < 1024; i++) {
		if (!free_particles)
			return;
		particle_t *p = free_particles;
		free_particles = p->next;
		p->next = active_particles;
		active_particles = p;
		p->die = cl.time + 5;
		p->color = ramp1[0];
		p->ramp = rand() & 3;
		if (i & 1) {
			p->type = pt_explode;
			for (s32 j = 0; j < 3; j++) {
				p->org[j] = org[j] + ((rand() % 32) - 16);
				p->vel[j] = (rand() % 512) - 256;
			}
		} else {
			p->type = pt_explode2;
			for (s32 j = 0; j < 3; j++) {
				p->org[j] = org[j] + ((rand() % 32) - 16);
				p->vel[j] = (rand() % 512) - 256;
			}
		}
	}
}

void R_ParticleExplosion2(vec3_t org, s32 colorStart, s32 colorLength)
{
	s32 colorMod = 0;
	for (s32 i = 0; i < 512; i++) {
		if (!free_particles)
			return;
		particle_t *p = free_particles;
		free_particles = p->next;
		p->next = active_particles;
		active_particles = p;
		p->die = cl.time + 0.3;
		p->color = colorStart + (colorMod % colorLength);
		colorMod++;
		p->type = pt_blob;
		for (s32 j = 0; j < 3; j++) {
			p->org[j] = org[j] + ((rand() % 32) - 16);
			p->vel[j] = (rand() % 512) - 256;
		}
	}
}

void R_BlobExplosion(vec3_t org)
{
	for (s32 i = 0; i < 1024; i++) {
		if (!free_particles)
			return;
		particle_t *p = free_particles;
		free_particles = p->next;
		p->next = active_particles;
		active_particles = p;
		p->die = cl.time + 1 + (rand() & 8) * 0.05;
		if (i & 1) {
			p->type = pt_blob;
			p->color = 66 + rand() % 6;
			for (s32 j = 0; j < 3; j++) {
				p->org[j] = org[j] + ((rand() % 32) - 16);
				p->vel[j] = (rand() % 512) - 256;
			}
		} else {
			p->type = pt_blob2;
			p->color = 150 + rand() % 6;
			for (s32 j = 0; j < 3; j++) {
				p->org[j] = org[j] + ((rand() % 32) - 16);
				p->vel[j] = (rand() % 512) - 256;
			}
		}
	}
}

void R_RunParticleEffect(vec3_t org, vec3_t dir, s32 color, s32 count)
{
	for (s32 i = 0; i < count; i++) {
		if (!free_particles)
			return;
		particle_t *p = free_particles;
		free_particles = p->next;
		p->next = active_particles;
		active_particles = p;
		if (count == 1024) { // rocket explosion
			p->die = cl.time + 5;
			p->color = ramp1[0];
			p->ramp = rand() & 3;
			if (i & 1) {
				p->type = pt_explode;
				for (s32 j = 0; j < 3; j++) {
					p->org[j] =
					    org[j] + ((rand() % 32) - 16);
					p->vel[j] = (rand() % 512) - 256;
				}
			} else {
				p->type = pt_explode2;
				for (s32 j = 0; j < 3; j++) {
					p->org[j] =
					    org[j] + ((rand() % 32) - 16);
					p->vel[j] = (rand() % 512) - 256;
				}
			}
		} else {
			p->die = cl.time + 0.1 * (rand() % 5);
			p->color = (color & ~7) + (rand() & 7);
			p->type = pt_slowgrav;
			for (s32 j = 0; j < 3; j++) {
				p->org[j] = org[j] + ((rand() & 15) - 8);
				p->vel[j] = dir[j] * 15;
			}
		}
	}
}

void R_LavaSplash(vec3_t org)
{
	for (s32 i = -16; i < 16; i++)
		for (s32 j = -16; j < 16; j++)
			for (s32 k = 0; k < 1; k++) {
				if (!free_particles)
					return;
				particle_t *p = free_particles;
				free_particles = p->next;
				p->next = active_particles;
				active_particles = p;
				p->die = cl.time + 2 + (rand() & 31) * 0.02;
				p->color = 224 + (rand() & 7);
				p->type = pt_slowgrav;
				vec3_t dir;
				dir[0] = j * 8 + (rand() & 7);
				dir[1] = i * 8 + (rand() & 7);
				dir[2] = 256;
				p->org[0] = org[0] + dir[0];
				p->org[1] = org[1] + dir[1];
				p->org[2] = org[2] + (rand() & 63);
				VectorNormalize(dir);
				f32 vel = 50 + (rand() & 63);
				VectorScale(dir, vel, p->vel);
			}
}

void R_TeleportSplash(vec3_t org)
{
	for (s32 i = -16; i < 16; i += 4)
		for (s32 j = -16; j < 16; j += 4)
			for (s32 k = -24; k < 32; k += 4) {
				if (!free_particles)
					return;
				particle_t *p = free_particles;
				free_particles = p->next;
				p->next = active_particles;
				active_particles = p;
				p->die = cl.time + 0.2 + (rand() & 7) * 0.02;
				p->color = 7 + (rand() & 7);
				p->type = pt_slowgrav;
				vec3_t dir;
				dir[0] = j * 8;
				dir[1] = i * 8;
				dir[2] = k * 8;
				p->org[0] = org[0] + i + (rand() & 3);
				p->org[1] = org[1] + j + (rand() & 3);
				p->org[2] = org[2] + k + (rand() & 3);
				VectorNormalize(dir);
				f32 vel = 50 + (rand() & 63);
				VectorScale(dir, vel, p->vel);
			}
}

void R_RocketTrail(vec3_t start, vec3_t end, s32 type)
{
	static s32 tracercount;
	vec3_t vec;
	VectorSubtract(end, start, vec);
	f32 len = VectorNormalize(vec);
	s32 dec;
	if (type < 128)
		dec = 3;
	else {
		dec = 1;
		type -= 128;
	}
	while (len > 0) {
		len -= dec;
		if (!free_particles)
			return;
		particle_t *p = free_particles;
		free_particles = p->next;
		p->next = active_particles;
		active_particles = p;
		VectorCopy(vec3_origin, p->vel);
		p->die = cl.time + 2;
		switch (type) {
		case 0: // rocket trail
			p->ramp = (rand() & 3);
			p->color = ramp3[(s32)p->ramp];
			p->type = pt_fire;
			for (s32 j = 0; j < 3; j++)
				p->org[j] = start[j] + ((rand() % 6) - 3);
			break;

		case 1: // smoke smoke
			p->ramp = (rand() & 3) + 2;
			p->color = ramp3[(s32)p->ramp];
			p->type = pt_fire;
			for (s32 j = 0; j < 3; j++)
				p->org[j] = start[j] + ((rand() % 6) - 3);
			break;
		case 2: // blood
			p->type = pt_grav;
			p->color = 67 + (rand() & 3);
			for (s32 j = 0; j < 3; j++)
				p->org[j] = start[j] + ((rand() % 6) - 3);
			break;
		case 3:
		case 5: // tracer
			p->die = cl.time + 0.5;
			p->type = pt_static;
			if (type == 3)
				p->color = 52 + ((tracercount & 4) << 1);
			else
				p->color = 230 + ((tracercount & 4) << 1);
			tracercount++;
			VectorCopy(start, p->org);
			if (tracercount & 1) {
				p->vel[0] = 30 * vec[1];
				p->vel[1] = 30 * -vec[0];
			} else {
				p->vel[0] = 30 * -vec[1];
				p->vel[1] = 30 * vec[0];
			}
			break;
		case 4: // slight blood
			p->type = pt_grav;
			p->color = 67 + (rand() & 3);
			for (s32 j = 0; j < 3; j++)
				p->org[j] = start[j] + ((rand() % 6) - 3);
			len -= 3;
			break;
		case 6: // voor trail
			p->color = 9 * 16 + 8 + (rand() & 3);
			p->type = pt_static;
			p->die = cl.time + 0.3;
			for (s32 j = 0; j < 3; j++)
				p->org[j] = start[j] + ((rand() & 15) - 8);
			break;
		}
		VectorAdd(start, vec, start);
	}
}

void R_DrawParticles()
{
	VectorScale(vright, xscaleshrink, r_pright);
	VectorScale(vup, yscaleshrink, r_pup);
	VectorCopy(vpn, r_ppn);
	f32 frametime = cl.time - cl.oldtime;
	f32 time3 = frametime * 15;
	f32 time2 = frametime * 10;
	f32 time1 = frametime * 5;
	f32 grav = frametime * sv_gravity.value * 0.05;
	f32 dvel = 4 * frametime;
	particle_t *kill;
	for (;;) {
		kill = active_particles;
		if (kill && kill->die < cl.time) {
			active_particles = kill->next;
			kill->next = free_particles;
			free_particles = kill;
			continue;
		}
		break;
	}
	for (particle_t *p = active_particles; p; p = p->next) {
		for (;;) {
			kill = p->next;
			if (kill && kill->die < cl.time) {
				p->next = kill->next;
				kill->next = free_particles;
				free_particles = kill;
				continue;
			}
			break;
		}
		D_DrawParticle(p);
		p->org[0] += p->vel[0] * frametime;
		p->org[1] += p->vel[1] * frametime;
		p->org[2] += p->vel[2] * frametime;
		switch (p->type) {
		case pt_static:
			break;
		case pt_fire:
			p->ramp += time1;
			if (p->ramp >= 6)
				p->die = -1;
			else
				p->color = ramp3[(s32)p->ramp];
			p->vel[2] += grav;
			break;

		case pt_explode:
			p->ramp += time2;
			if (p->ramp >= 8)
				p->die = -1;
			else
				p->color = ramp1[(s32)p->ramp];
			for (s32 i = 0; i < 3; i++)
				p->vel[i] += p->vel[i] * dvel;
			p->vel[2] -= grav;
			break;
		case pt_explode2:
			p->ramp += time3;
			if (p->ramp >= 8)
				p->die = -1;
			else
				p->color = ramp2[(s32)p->ramp];
			for (s32 i = 0; i < 3; i++)
				p->vel[i] -= p->vel[i] * frametime;
			p->vel[2] -= grav;
			break;
		case pt_blob:
			for (s32 i = 0; i < 3; i++)
				p->vel[i] += p->vel[i] * dvel;
			p->vel[2] -= grav;
			break;
		case pt_blob2:
			for (s32 i = 0; i < 2; i++)
				p->vel[i] -= p->vel[i] * dvel;
			p->vel[2] -= grav;
			break;
		case pt_grav:
		case pt_slowgrav:
			p->vel[2] -= grav;
			break;
		}
	}
}
