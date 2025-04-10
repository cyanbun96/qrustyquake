// Copyright (C) 1996-2001 Id Software, Inc.
// Copyright (C) 2002-2005 John Fitzgibbons and others
// Copyright (C) 2007-2008 Kristian Duske
// GPLv3 See LICENSE for details.

#define MOVE_NORMAL 0
#define MOVE_NOMONSTERS 1
#define MOVE_MISSILE 2

typedef struct
{
	vec3_t normal;
	float dist;
} plane_t;

typedef struct
{
	qboolean allsolid; // if true, plane is not valid
	qboolean startsolid; // if true, the initial point was in a solid area
	qboolean inopen, inwater;
	float fraction; // time completed, 1.0 = didn't hit anything
	vec3_t endpos; // final position
	plane_t plane; // surface normal at impact
	edict_t *ent; // entity the surface is on
} trace_t;

void SV_ClearWorld();
// called after the world model has been loaded, before linking any entities

void SV_UnlinkEdict(edict_t *ent);
// call before removing an entity, and before trying to move one,
// so it doesn't clip against itself
// flags ent->v.modified

void SV_LinkEdict(edict_t *ent, qboolean touch_triggers);
// Needs to be called any time an entity changes origin, mins, maxs, or solid
// flags ent->v.modified
// sets ent->v.absmin and ent->v.absmax
// if touchtriggers, calls prog functions for the intersected triggers

int SV_PointContents(vec3_t p);
int SV_TruePointContents(vec3_t p);
// returns the CONTENTS_* value from the world at the given point.
// does not check any entities at all
// the non-true version remaps the water current contents to content_water

edict_t *SV_TestEntityPosition(edict_t *ent);
trace_t SV_Move(vec3_t start, vec3_t mins, vec3_t maxs,
		vec3_t end, int type, edict_t *passedict);
// mins and maxs are relative
// if the entire move stays in a solid volume, trace.allsolid will be set
// if the starting point is in a solid, it will be allowed to move out
// to an open area
// nomonsters is used for line of sight or edge testing, where monsters
// shouldn't be considered solid objects
// passedict is explicitly excluded from clipping checks(normally NULL)
