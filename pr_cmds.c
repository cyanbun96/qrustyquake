// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
#include "quakedef.h"

static s8 pr_string_temp[STRINGTEMP_BUFFERS][STRINGTEMP_LENGTH];
static u8 pr_string_tempindex = 0;
static u8 *checkpvs; //ericw -- changed to malloc
static s32 checkpvs_capacity;
static s32 c_invis, c_notvis;

static s8 *PR_GetTempString()
{ return pr_string_temp[(STRINGTEMP_BUFFERS-1) & ++pr_string_tempindex]; }

static const s8* PF_GetStringArg(s32 idx, void* userdata)
{
	if(userdata) idx += *(s32*)userdata;
	if(idx < 0 || idx >= pr_argc) return "";
	return LOC_GetString(G_STRING(OFS_PARM0 + idx * 3));
}

static s8 *PF_VarString(s32 first)
{
	static s8 out[1024];
	out[0] = 0;
	if(first >= pr_argc) return out;
	const s8 *format = LOC_GetString(G_STRING((OFS_PARM0 + first * 3)));
	size_t s = 0;
	if(LOC_HasPlaceholders(format)) {
		s32 offset = first + 1;
		s = LOC_Format(format,PF_GetStringArg,&offset,out,sizeof(out));
	} else for(s32 i = first; i < pr_argc; i++) {
	    s=q_strlcat(out,LOC_GetString(G_STRING(OFS_PARM0+i*3)),sizeof(out));
		if(s >= sizeof(out)) {
			printf("PF_VarString: overflow(string truncated)\n");
			return out;
		}
	}
	return out;
}

static void PF_error()
{ // This is a TERMINAL error, which will kill off the entire server.
	s8 *s = PF_VarString(0);
	Con_Printf("======SERVER ERROR in %s:\n%s\n",
			PR_GetString(pr_xfunction->s_name), s);
	edict_t *ed = PROG_TO_EDICT(pr_global_struct->self);
	ED_Print(ed);
	Host_Error("Program error");
}


static void PF_objerror() // Dumps out self, then an error message. The program
{ // is aborted and self is removed, but the level can continue.
	s8 *s = PF_VarString(0);
	Con_Printf("======OBJECT ERROR in %s:\n%s\n",
			PR_GetString(pr_xfunction->s_name), s);
	edict_t *ed = PROG_TO_EDICT(pr_global_struct->self);
	ED_Print(ed);
	ED_Free(ed);
}

static void PF_makevectors()
{ // Writes new values for v_forward, v_up, and v_right based on angles
	AngleVectors(G_VECTOR(OFS_PARM0), pr_global_struct->v_forward,
			pr_global_struct->v_right, pr_global_struct->v_up);
}

// This is the only valid way to move an object without using the physics of the
// world(setting velocity and waiting). Directly changing origin will not set
// internal links correctly, so clipping would be messed up. This should be
// called when an object is spawned, and then only if it is teleported.
static void PF_setorigin()
{
	edict_t *e = G_EDICT(OFS_PARM0);
	f32 *org = G_VECTOR(OFS_PARM1);
	VectorCopy(org, e->v.origin);
	SV_LinkEdict(e, 0);
}

static void SetMinMaxSize(edict_t *e, f32 *minvec, f32 *maxvec, SDL_UNUSED bool rotate)
{
	vec3_t rmin, rmax;
	for(s32 i = 0; i < 3; i++)
		if(minvec[i] > maxvec[i]) PR_RunError("backwards mins/maxs");
	VectorCopy(minvec, rmin);
	VectorCopy(maxvec, rmax);
	VectorCopy(rmin, e->v.mins); // set derived values
	VectorCopy(rmax, e->v.maxs);
	VectorSubtract(maxvec, minvec, e->v.size);
	SV_LinkEdict(e, 0);
}

static void PF_setsize() // the size box is rotated by the current angle
{
	edict_t *e = G_EDICT(OFS_PARM0);
	f32 *minvec = G_VECTOR(OFS_PARM1);
	f32 *maxvec = G_VECTOR(OFS_PARM2);
	SetMinMaxSize(e, minvec, maxvec, 0);
}

static void PF_setmodel()
{
	edict_t *e = G_EDICT(OFS_PARM0);
	const s8 *m = G_STRING(OFS_PARM1);
	// check to see if model was properly precached
	s32 i = 0;
	const s8 **check = sv.model_precache;
	for(; *check; i++, check++)
		if(!strcmp(*check, m)) break;
	if(!*check) PR_RunError("no precache: %s", m);
	e->v.model = PR_SetEngineString(*check);
	e->v.modelindex = i;
	model_t *mod = sv.models[ (s32)e->v.modelindex];
	if(mod){ //johnfitz -- correct physics cullboxes for bmodels
		if(mod->type == mod_brush)
			SetMinMaxSize(e, mod->clipmins, mod->clipmaxs, 1);
		else SetMinMaxSize(e, mod->mins, mod->maxs, 1);
	} else SetMinMaxSize(e, vec3_origin, vec3_origin, 1);
}

static void PF_bprint() // broadcast print to everyone on server
{ s8 *s = PF_VarString(0); SV_BroadcastPrintf("%s", s); }

static void PF_sprint() // single print to a specific client
{
	s32 entnum = G_EDICTNUM(OFS_PARM0);
	s8 *s = PF_VarString(1);
	if(entnum < 1 || entnum > svs.maxclients) {
		Con_Printf("tried to sprint to a non-client\n");
		return;
	}
	client_t *client = &svs.clients[entnum-1];
	MSG_WriteChar(&client->message,svc_print);
	MSG_WriteString(&client->message, s );
}

static void PF_centerprint() // single centerprint to a specific client
{
	s32 entnum = G_EDICTNUM(OFS_PARM0);
	s8 *s = PF_VarString(1);
	if(entnum < 1 || entnum > svs.maxclients) {
		Con_Printf("tried to sprint to a non-client\n");
		return;
	}
	client_t *client = &svs.clients[entnum-1];
	MSG_WriteChar(&client->message,svc_centerprint);
	MSG_WriteString(&client->message, s);
}

static void PF_normalize()
{
	vec3_t newvalue;
	f32 *v1 = G_VECTOR(OFS_PARM0);
	f64 new_temp = (f64)v1[0]*v1[0] + (f64)v1[1]*v1[1] + (f64)v1[2]*v1[2];
	new_temp = sqrt(new_temp);
	if(new_temp == 0) newvalue[0] = newvalue[1] = newvalue[2] = 0;
	else {
		new_temp = 1 / new_temp;
		newvalue[0] = v1[0] * new_temp;
		newvalue[1] = v1[1] * new_temp;
		newvalue[2] = v1[2] * new_temp;
	}
	VectorCopy(newvalue, G_VECTOR(OFS_RETURN));
}

static void PF_vlen()
{
	f32 *v1 = G_VECTOR(OFS_PARM0);
	f64 new_temp = (f64)v1[0]*v1[0] + (f64)v1[1]*v1[1] + (f64)v1[2]*v1[2];
	new_temp = sqrt(new_temp);
	G_FLOAT(OFS_RETURN) = new_temp;
}

static void PF_vectoyaw()
{
	f32 *value1 = G_VECTOR(OFS_PARM0);
	f32 yaw;
	if(value1[1] == 0 && value1[0] == 0) yaw = 0;
	else {
		yaw = (s32) (atan2(value1[1], value1[0]) * 180 / M_PI);
		if(yaw < 0) yaw += 360;
	}
	G_FLOAT(OFS_RETURN) = yaw;
}

static void PF_vectoangles()
{
	f32 yaw, pitch;
	f32 *value1 = G_VECTOR(OFS_PARM0);
	if(value1[1] == 0 && value1[0] == 0) {
		yaw = 0;
		pitch = value1[2] > 0 ? 90 : 270;
	} else {
		yaw = (s32) (atan2(value1[1], value1[0]) * 180 / M_PI);
		if(yaw < 0) yaw += 360;
		f32 forward = sqrt(value1[0]*value1[0] + value1[1]*value1[1]);
		pitch = (s32) (atan2(value1[2], forward) * 180 / M_PI);
		if(pitch < 0) pitch += 360;
	}
	G_FLOAT(OFS_RETURN+0) = pitch;
	G_FLOAT(OFS_RETURN+1) = yaw;
	G_FLOAT(OFS_RETURN+2) = 0;
}

static void PF_random() // Returns a number from 0 <= num < 1
{
	f32 num = (rand() & 0x7fff) / ((f32)0x7fff);
	G_FLOAT(OFS_RETURN) = num;
}

static void PF_particle()
{
	f32 *org = G_VECTOR(OFS_PARM0);
	f32 *dir = G_VECTOR(OFS_PARM1);
	f32 color = G_FLOAT(OFS_PARM2);
	f32 count = G_FLOAT(OFS_PARM3);
	SV_StartParticle(org, dir, color, count);
}

static void PF_ambientsound()
{
	s32 large = 0; //johnfitz -- PROTOCOL_FITZQUAKE
	f32 *pos = G_VECTOR(OFS_PARM0);
	const s8 *samp = G_STRING(OFS_PARM1);
	f32 vol = G_FLOAT(OFS_PARM2);
	f32 attenuation = G_FLOAT(OFS_PARM3);
	// check to see if samp was properly precached
	s32 soundnum = 0;
	const s8 **check = sv.sound_precache;
	for(; *check; check++, soundnum++)
		if(!strcmp(*check, samp)) break;
	if(!*check) { Con_Printf("no precache: %s\n", samp); return; }
	if(soundnum > 255) {
		if(sv.protocol == PROTOCOL_NETQUAKE) return;
		else large = 1;
	}
	SV_ReserveSignonSpace(17);
	// add an svc_spawnambient command to the level signon packet
	if(large) MSG_WriteByte(sv.signon,svc_spawnstaticsound2);
	else MSG_WriteByte(sv.signon,svc_spawnstaticsound);
	for(s32 i = 0; i < 3; i++)
		MSG_WriteCoord(sv.signon, pos[i], sv.protocolflags);
	if(large) MSG_WriteShort(sv.signon, soundnum);
	else MSG_WriteByte(sv.signon, soundnum);
	MSG_WriteByte(sv.signon, vol*255);
	MSG_WriteByte(sv.signon, attenuation*64);
}

// Each entity can have eight independant sound sources, like voice,
// weapon, feet, etc.
// Channel 0 is an auto-allocate channel, the others override anything
// already running on that entity/channel pair.
// An attenuation of 0 will play full volume everywhere in the level.
// Larger attenuations will drop off.
static void PF_sound()
{
	edict_t *entity = G_EDICT(OFS_PARM0);
	s32 channel = G_FLOAT(OFS_PARM1);
	const s8 *sample = G_STRING(OFS_PARM2);
	s32 volume = G_FLOAT(OFS_PARM3) * 255;
	f32 attenuation = G_FLOAT(OFS_PARM4);
	SV_StartSound(entity, channel, sample, volume, attenuation);
}

static void PF_break()
{ Con_Printf("break statement\n"); *(s32 *)-4 = 0; } // dump to debugger

// Used for use tracing and shot targeting
// Traces are blocked by bbox and exact bsp entityes, and also slide box
// entities if the tryents flag is set.
static void PF_traceline()
{
	f32 *v1 = G_VECTOR(OFS_PARM0);
	f32 *v2 = G_VECTOR(OFS_PARM1);
	s32 nomonsters = G_FLOAT(OFS_PARM2);
	edict_t *ent = G_EDICT(OFS_PARM3);
	// FIXME: Why do we hit this with certain progs.dat?
	if(developer.value){
		if(isnan(v1[0]) || isnan(v1[1]) || isnan(v1[2]) ||
			isnan(v2[0]) || isnan(v2[1]) || isnan(v2[2])){
	printf("NAN in traceline:\nv1(%f %f %f) v2(%f %f %f)\nentity %d\n",
		v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], NUM_FOR_EDICT(ent));
		}
	}
	if(isnan(v1[0]) || isnan(v1[1]) || isnan(v1[2]))
		v1[0] = v1[1] = v1[2] = 0;
	if(isnan(v2[0]) || isnan(v2[1]) || isnan(v2[2]))
		v2[0] = v2[1] = v2[2] = 0;
	trace_t trace = SV_Move(v1,vec3_origin,vec3_origin,v2,nomonsters,ent);
	pr_global_struct->trace_allsolid = trace.allsolid;
	pr_global_struct->trace_startsolid = trace.startsolid;
	pr_global_struct->trace_fraction = trace.fraction;
	pr_global_struct->trace_inwater = trace.inwater;
	pr_global_struct->trace_inopen = trace.inopen;
	VectorCopy(trace.endpos, pr_global_struct->trace_endpos);
	VectorCopy(trace.plane.normal, pr_global_struct->trace_plane_normal);
	pr_global_struct->trace_plane_dist = trace.plane.dist;
	if(trace.ent) pr_global_struct->trace_ent = EDICT_TO_PROG(trace.ent);
	else pr_global_struct->trace_ent = EDICT_TO_PROG(sv.edicts);
}

static s32 PF_newcheckclient(s32 check)
{
	s32 i;
	edict_t *ent;
	if(check < 1) check = 1; // cycle to the next one
	if(check > svs.maxclients) check = svs.maxclients;
	if(check == svs.maxclients) i = 1;
	else i = check + 1;
	for(;; i++) {
		if(i == svs.maxclients+1) i = 1;
		ent = EDICT_NUM(i);
		if(i == check) break; // didn't find anything else
		if(ent->free) continue;
		if(ent->v.health <= 0) continue;
		if((s32)ent->v.flags & FL_NOTARGET) continue;
		break; // anything that is a client, or has a client as an enemy
	}
	vec3_t org; // get the PVS for the entity
	VectorAdd(ent->v.origin, ent->v.view_ofs, org);
	mleaf_t *leaf = Mod_PointInLeaf(org, sv.worldmodel);
	u8 *pvs = Mod_LeafPVS(leaf, sv.worldmodel);
	s32 pvsbytes = (sv.worldmodel->numleafs+7)>>3;
	if(checkpvs == NULL || pvsbytes > checkpvs_capacity) {
		checkpvs_capacity = pvsbytes;
		checkpvs = (u8 *) realloc(checkpvs, checkpvs_capacity);
		if(!checkpvs)
Sys_Error("PF_newcheckclient: realloc() failed on %d bytes", checkpvs_capacity);
	}
	memcpy(checkpvs, pvs, pvsbytes);
	return i;
}

// Returns a client(or object that has a client enemy) that would be a
// valid target.
// If there are more than one valid options, they are cycled each frame
// If(self.origin + self.viewofs) is not in the PVS of the current target,
// it is not returned at all.
static void PF_checkclient()
{
	// find a new check if on a new frame
	if(sv.time - sv.lastchecktime >= 0.1) {
		sv.lastcheck = PF_newcheckclient(sv.lastcheck);
		sv.lastchecktime = sv.time;
	}
	edict_t *ent = EDICT_NUM(sv.lastcheck);
	if(ent->free || ent->v.health<=0){ //return check if it might be visible
		RETURN_EDICT(sv.edicts);
		return;
	}
	// if current entity can't possibly see the check entity, return 0
	edict_t *self = PROG_TO_EDICT(pr_global_struct->self);
	vec3_t view;
	VectorAdd(self->v.origin, self->v.view_ofs, view);
	mleaf_t *leaf = Mod_PointInLeaf(view, sv.worldmodel);
	s32 l = (leaf - sv.worldmodel->leafs) - 1;
	if((l < 0) || !(checkpvs[l>>3] & (1 << (l & 7)))) {
		c_notvis++;
		RETURN_EDICT(sv.edicts);
		return;
	}
	c_invis++; // might be able to see it
	RETURN_EDICT(ent);
}

static void PF_stuffcmd()
{ // Sends text over to the client's execution buffer
	s32 entnum = G_EDICTNUM(OFS_PARM0);
	if(entnum < 1 || entnum > svs.maxclients)
		PR_RunError("Parm 0 not a client");
	const s8 *str = G_STRING(OFS_PARM1);
	client_t *old = host_client;
	host_client = &svs.clients[entnum-1];
	Host_ClientCommands("%s", str);
	host_client = old;
}

static void PF_localcmd()
{ // Sends text over to the client's execution buffer
	const s8 *str = G_STRING(OFS_PARM0);
	Cbuf_AddText(str);
}

static void PF_cvar()
{
	const s8 *str = G_STRING(OFS_PARM0);
	G_FLOAT(OFS_RETURN) = Cvar_VariableValue(str);
}

static void PF_cvar_set()
{
	const s8 *var = G_STRING(OFS_PARM0);
	const s8 *val = G_STRING(OFS_PARM1);
	Cvar_Set(var, val);
}


static void PF_findradius()
{ // Returns a chain of entities that have origins within a spherical area
	edict_t *chain = (edict_t *)sv.edicts;
	f32 *org = G_VECTOR(OFS_PARM0);
	f32 rad = G_FLOAT(OFS_PARM1);
	rad *= rad;
	edict_t *ent = NEXT_EDICT(sv.edicts);
	for(s32 i = 1; i < sv.num_edicts; i++, ent = NEXT_EDICT(ent)) {
		f32 d, lensq;
		if(ent->free) continue;
		if(ent->v.solid == SOLID_NOT) continue;
		d=org[0]-(ent->v.origin[0]+(ent->v.mins[0]+ent->v.maxs[0])*0.5);
		lensq = d * d;
		if(lensq > rad) continue;
		d=org[1]-(ent->v.origin[1]+(ent->v.mins[1]+ent->v.maxs[1])*0.5);
		lensq += d * d;
		if(lensq > rad) continue;
		d=org[2]-(ent->v.origin[2]+(ent->v.mins[2]+ent->v.maxs[2])*0.5);
		lensq += d * d;
		if(lensq > rad) continue;
		ent->v.chain = EDICT_TO_PROG(chain);
		chain = ent;
	}
	RETURN_EDICT(chain);
}

static void PF_dprint() { Con_DPrintf("%s",PF_VarString(0)); }

static void PF_ftos()
{
	f32 v = G_FLOAT(OFS_PARM0);
	s8 *s = PR_GetTempString();
	if(v == (s32)v) sprintf(s, "%d",(s32)v);
	else sprintf(s, "%5.1f",v);
	G_INT(OFS_RETURN) = PR_SetEngineString(s);
}

static void PF_fabs()
{
	f32 v;
	v = G_FLOAT(OFS_PARM0);
	G_FLOAT(OFS_RETURN) = fabs(v);
}

static void PF_vtos()
{
	s8 *s = PR_GetTempString();
	sprintf(s, "'%5.1f %5.1f %5.1f'", G_VECTOR(OFS_PARM0)[0],
			G_VECTOR(OFS_PARM0)[1], G_VECTOR(OFS_PARM0)[2]);
	G_INT(OFS_RETURN) = PR_SetEngineString(s);
}

static void PF_Spawn()
{
	edict_t *ed = ED_Alloc();
	RETURN_EDICT(ed);
}

static void PF_Remove()
{
	edict_t *ed;
	ed = G_EDICT(OFS_PARM0);
	ED_Free(ed);
}

static void PF_Find()
{
	s32 e = G_EDICTNUM(OFS_PARM0);
	s32 f = G_INT(OFS_PARM1);
	const s8 *s = G_STRING(OFS_PARM2);
	if(!s) PR_RunError("PF_Find: bad search string");
	for(e++; e < sv.num_edicts; e++) {
		edict_t *ed = EDICT_NUM(e);
		if(ed->free) continue;
		const s8 *t = E_STRING(ed,f);
		if(!t) continue;
		if(!strcmp(t,s)) {
			RETURN_EDICT(ed);
			return;
		}
	}
	RETURN_EDICT(sv.edicts);
}

static void PR_CheckEmptyString(const s8 *s)
{ if(s[0] <= ' ') PR_RunError("Bad string"); }

static void PF_precache_file() // only used to copy files with qcc, does nothing
{ G_INT(OFS_RETURN) = G_INT(OFS_PARM0); }

static void PF_precache_sound()
{
	if(sv.state != ss_loading)
PR_RunError("PF_Precache_*: Precache can only be done in spawn functions");
	const s8 *s = G_STRING(OFS_PARM0);
	G_INT(OFS_RETURN) = G_INT(OFS_PARM0);
	PR_CheckEmptyString(s);
	for(s32 i = 0; i < MAX_SOUNDS; i++) {
		if(!sv.sound_precache[i]) { sv.sound_precache[i] = s; return; }
		if(!strcmp(sv.sound_precache[i], s)) return;
	}
	PR_RunError("PF_precache_sound: overflow");
}

static void PF_precache_model()
{
	if(sv.state != ss_loading)
PR_RunError("PF_Precache_*: Precache can only be done in spawn functions");
	const s8 *s = G_STRING(OFS_PARM0);
	G_INT(OFS_RETURN) = G_INT(OFS_PARM0);
	PR_CheckEmptyString(s);
	for(s32 i = 0; i < MAX_MODELS; i++) {
		if(!sv.model_precache[i]) {
			sv.model_precache[i] = s;
			sv.models[i] = Mod_ForName(s, 1);
			return;
		}
		if(!strcmp(sv.model_precache[i], s)) return;
	}
	PR_RunError("PF_precache_model: overflow");
}

static void PF_coredump() { ED_PrintEdicts(); }
static void PF_traceon() { pr_trace = 1; }
static void PF_traceoff() { pr_trace = 0; }
static void PF_eprint() { ED_PrintNum(G_EDICTNUM(OFS_PARM0)); }

static void PF_walkmove()
{
	edict_t *ent = PROG_TO_EDICT(pr_global_struct->self);
	f32 yaw = G_FLOAT(OFS_PARM0);
	f32 dist = G_FLOAT(OFS_PARM1);
	if(!((s32)ent->v.flags & (FL_ONGROUND|FL_FLY|FL_SWIM))) {
		G_FLOAT(OFS_RETURN) = 0;
		return;
	}
	yaw = yaw * M_PI * 2 / 360;
	vec3_t move;
	move[0] = cos(yaw) * dist;
	move[1] = sin(yaw) * dist;
	move[2] = 0;
	// save program state, because SV_movestep may call other progs
	dfunction_t *oldf = pr_xfunction;
	s32 oldself = pr_global_struct->self;
	G_FLOAT(OFS_RETURN) = SV_movestep(ent, move, 1);
	pr_xfunction = oldf; // restore program state
	pr_global_struct->self = oldself;
}

static void PF_droptofloor()
{
	edict_t *ent = PROG_TO_EDICT(pr_global_struct->self);
	vec3_t end;
	VectorCopy(ent->v.origin, end);
	end[2] -= 256;
	trace_t trace=SV_Move(ent->v.origin,ent->v.mins,ent->v.maxs,end,0,ent);
	if(trace.fraction == 1 || trace.allsolid) G_FLOAT(OFS_RETURN) = 0;
	else {
		VectorCopy(trace.endpos, ent->v.origin);
		SV_LinkEdict(ent, 0);
		ent->v.flags = (s32)ent->v.flags | FL_ONGROUND;
		ent->v.groundentity = EDICT_TO_PROG(trace.ent);
		G_FLOAT(OFS_RETURN) = 1;
	}
}

static void PF_lightstyle()
{
	s32 style = G_FLOAT(OFS_PARM0);
	const s8 *val = G_STRING(OFS_PARM1);
	// bounds check to avoid clobbering sv struct
	if(style < 0 || style >= MAX_LIGHTSTYLES) {
		printf("PF_lightstyle: invalid style %d\n", style);
		return;
	}
	sv.lightstyles[style] = val; // change the string in sv
	// send message to all clients on this server
	if(sv.state != ss_active) return;
	client_t *client = svs.clients;
	s32 j = 0;
	for(; j < svs.maxclients; j++, client++) {
		if(client->active || client->spawned) {
			MSG_WriteChar(&client->message, svc_lightstyle);
			MSG_WriteChar(&client->message, style);
			MSG_WriteString(&client->message, val);
		}
	}
}

static void PF_rint()
{
	f32 f = G_FLOAT(OFS_PARM0);
	if(f > 0) G_FLOAT(OFS_RETURN) = (s32)(f + 0.5);
	else G_FLOAT(OFS_RETURN) = (s32)(f - 0.5);
}

static void PF_floor() { G_FLOAT(OFS_RETURN) = floor(G_FLOAT(OFS_PARM0)); }
static void PF_ceil() { G_FLOAT(OFS_RETURN) = ceil(G_FLOAT(OFS_PARM0)); }

static void PF_checkbottom()
{
	edict_t *ent = G_EDICT(OFS_PARM0);
	G_FLOAT(OFS_RETURN) = SV_CheckBottom(ent);
}

static void PF_pointcontents()
{
	f32 *v = G_VECTOR(OFS_PARM0);
	G_FLOAT(OFS_RETURN) = SV_PointContents(v);
}

static void PF_nextent()
{
	s32 i = G_EDICTNUM(OFS_PARM0);
	while(1) {
		i++;
		if(i == sv.num_edicts) {
			RETURN_EDICT(sv.edicts);
			return;
		}
		edict_t *ent = EDICT_NUM(i);
		if(!ent->free) {
			RETURN_EDICT(ent);
			return;
		}
	}
}


static void PF_aim() // Pick a vector for the player to shoot along
{
	vec3_t start, dir, end, bestdir;
	edict_t *ent = G_EDICT(OFS_PARM0);
	/*speed*/(void)G_FLOAT(OFS_PARM1);
	VectorCopy(ent->v.origin, start);
	start[2] += 20;
	// try sending a trace straight
	VectorCopy(pr_global_struct->v_forward, dir);
	VectorMA(start, 2048, dir, end);
	trace_t tr = SV_Move(start, vec3_origin, vec3_origin, end, 0, ent);
	if(tr.ent && tr.ent->v.takedamage == DAMAGE_AIM && (!teamplay.value
			|| ent->v.team <= 0 || ent->v.team != tr.ent->v.team)) {
		VectorCopy(pr_global_struct->v_forward, G_VECTOR(OFS_RETURN));
		return;
	}
	VectorCopy(dir, bestdir); // try all possible entities
	f32 bestdist = sv_aim.value;
	edict_t *bestent = NULL;
	edict_t *check = NEXT_EDICT(sv.edicts);
	for(s32 i = 1; i < sv.num_edicts; i++, check = NEXT_EDICT(check)) {
		if(check->v.takedamage != DAMAGE_AIM) continue;
		if(check == ent) continue;
		if(teamplay.value &&ent->v.team>0 &&ent->v.team==check->v.team)
			continue; // don't aim at teammate
		for(s32 j = 0; j < 3; j++)
			end[j] = check->v.origin[j] +
				0.5*(check->v.mins[j]+check->v.maxs[j]);
		VectorSubtract(end, start, dir);
		VectorNormalize(dir);
		f32 dist = DotProduct(dir, pr_global_struct->v_forward);
		if(dist < bestdist) continue; // to far to turn
		tr = SV_Move(start, vec3_origin, vec3_origin, end, 0, ent);
		if(tr.ent == check) { // can shoot at this one
			bestdist = dist;
			bestent = check;
		}
	}
	if(bestent) {
		VectorSubtract(bestent->v.origin, ent->v.origin, dir);
		f32 dist = DotProduct(dir, pr_global_struct->v_forward);
		VectorScale(pr_global_struct->v_forward, dist, end);
		end[2] = dir[2];
		VectorNormalize(end);
		VectorCopy(end, G_VECTOR(OFS_RETURN));
	}
	else VectorCopy(bestdir, G_VECTOR(OFS_RETURN));
}

void PF_changeyaw()
{ // This was a major timewaster in progs, so it was converted to C
	edict_t *ent;
	f32 ideal, current, move, speed;
	ent = PROG_TO_EDICT(pr_global_struct->self);
	current = anglemod( ent->v.angles[1] );
	ideal = ent->v.ideal_yaw;
	speed = ent->v.yaw_speed;
	if(current == ideal)
		return;
	move = ideal - current;
	if(ideal > current) {if(move >= 180) move = move - 360;}
	else                {if(move <= -180) move = move + 360;}
	if(move > 0) {if(move > speed) move = speed;}
	else         {if(move < -speed) move = -speed;}
	ent->v.angles[1] = anglemod(current + move);
}

static sizebuf_t *WriteDest()
{
	s32 entnum;
	edict_t *ent;
	s32 dest = G_FLOAT(OFS_PARM0);
	switch(dest) {
		case MSG_BROADCAST: return &sv.datagram;
		case MSG_ONE:
			ent = PROG_TO_EDICT(pr_global_struct->msg_entity);
			entnum = NUM_FOR_EDICT(ent);
			if(entnum < 1 || entnum > svs.maxclients)
				PR_RunError("WriteDest: not a client");
			return &svs.clients[entnum-1].message;
		case MSG_ALL: return &sv.reliable_datagram;
		case MSG_INIT: return sv.signon;
		default: PR_RunError("WriteDest: bad destination"); break;
	}
	return NULL;
}

static void PF_WriteByte() { MSG_WriteByte(WriteDest(), G_FLOAT(OFS_PARM1)); }
static void PF_WriteChar() { MSG_WriteChar(WriteDest(), G_FLOAT(OFS_PARM1)); }
static void PF_WriteShort() { MSG_WriteShort(WriteDest(), G_FLOAT(OFS_PARM1)); }
static void PF_WriteLong() { MSG_WriteLong(WriteDest(), G_FLOAT(OFS_PARM1)); }
static void PF_WriteAngle()
{ MSG_WriteAngle(WriteDest(), G_FLOAT(OFS_PARM1), sv.protocolflags); }
static void PF_WriteCoord()
{ MSG_WriteCoord(WriteDest(), G_FLOAT(OFS_PARM1), sv.protocolflags); }
static void PF_WriteString()
{ MSG_WriteString(WriteDest(), LOC_GetString(G_STRING(OFS_PARM1))); }
static void PF_WriteEntity()
{ MSG_WriteShort(WriteDest(), G_EDICTNUM(OFS_PARM1)); }

static void PF_makestatic()
{
	s32 bits = 0;
	edict_t *ent = G_EDICT(OFS_PARM0);
	if(ent->alpha == ENTALPHA_ZERO){ ED_Free(ent); return; }
	if(sv.protocol == PROTOCOL_NETQUAKE) {
		if(SV_ModelIndex(PR_GetString(ent->v.model)) & 0xFF00
				|| (s32)(ent->v.frame) & 0xFF00) {
			ED_Free(ent); //can't display the correct model & frame
			return; // so don't show it at all
		}
	} else {
		if(SV_ModelIndex(PR_GetString(ent->v.model)) & 0xFF00)
			bits |= B_LARGEMODEL;
		if((s32)(ent->v.frame) & 0xFF00)
			bits |= B_LARGEFRAME;
		if(ent->alpha != ENTALPHA_DEFAULT)
			bits |= B_ALPHA;
		if(sv.protocol == PROTOCOL_RMQ) {
			eval_t* val = GetEdictFieldValue(ent, "scale");
			if(val) ent->scale = ENTSCALE_ENCODE(val->_float);
			else ent->scale = ENTSCALE_DEFAULT;
			if(ent->scale != ENTSCALE_DEFAULT) bits |= B_SCALE;
		}
	}
	SV_ReserveSignonSpace(34);
	if(bits) {
		MSG_WriteByte(sv.signon, svc_spawnstatic2);
		MSG_WriteByte(sv.signon, bits);
	}
	else MSG_WriteByte(sv.signon, svc_spawnstatic);
	if(bits & B_LARGEMODEL)
	   MSG_WriteShort(sv.signon, SV_ModelIndex(PR_GetString(ent->v.model)));
	else
	   MSG_WriteByte(sv.signon, SV_ModelIndex(PR_GetString(ent->v.model)));
	if(bits & B_LARGEFRAME) MSG_WriteShort(sv.signon, ent->v.frame);
	else MSG_WriteByte(sv.signon, ent->v.frame);
	MSG_WriteByte(sv.signon, ent->v.colormap);
	MSG_WriteByte(sv.signon, ent->v.skin);
	for(s32 i = 0; i < 3; i++) {
		MSG_WriteCoord(sv.signon, ent->v.origin[i], sv.protocolflags);
		MSG_WriteAngle(sv.signon, ent->v.angles[i], sv.protocolflags);
	}
	if(bits & B_ALPHA) MSG_WriteByte(sv.signon, ent->alpha);
	if(bits & B_SCALE) MSG_WriteByte(sv.signon, ent->scale);
	ED_Free(ent); // throw the entity away now
}

static void PF_setspawnparms()
{
	edict_t *ent = G_EDICT(OFS_PARM0);
	s32 i = NUM_FOR_EDICT(ent);
	if(i < 1 || i > svs.maxclients)
		PR_RunError("Entity is not a client");
	// copy spawn parms out of the client_t
	client_t *client = svs.clients + (i-1);
	for(i = 0; i < NUM_SPAWN_PARMS; i++)
		(&pr_global_struct->parm1)[i] = client->spawn_parms[i];
}

static void PF_changelevel()
{
	// make sure we don't issue two changelevels
	if(svs.changelevel_issued) return;
	svs.changelevel_issued = 1;
	const s8 *s = G_STRING(OFS_PARM0);
	Cbuf_AddText(va("changelevel %s\n",s));
}

static void PF_finalefinished() { G_FLOAT(OFS_RETURN) = 0; }
static void PF_CheckPlayerEXFlags() { G_FLOAT(OFS_RETURN) = 0; }
static void PF_walkpathtogoal() { G_FLOAT(OFS_RETURN) = 0; /* PATH_ERROR */ }

static void PF_localsound()
{
	s32 entnum = G_EDICTNUM(OFS_PARM0);
	const s8 *sample = G_STRING(OFS_PARM1);
	if(entnum < 1 || entnum > svs.maxclients){
		Con_Printf("tried to localsound to a non-client\n");
		return;
	}
	SV_LocalSound(&svs.clients[entnum-1], sample);
}

static void PF_Fixme() { PR_RunError("unimplemented builtin"); }

static builtin_t pr_builtin[] = {
PF_Fixme,PF_makevectors,PF_setorigin,PF_setmodel,PF_setsize,PF_Fixme,PF_break,
PF_random,PF_sound,PF_normalize,PF_error,PF_objerror,PF_vlen,PF_vectoyaw,
PF_Spawn,PF_Remove,PF_traceline,PF_checkclient,PF_Find,PF_precache_sound,
PF_precache_model,PF_stuffcmd,PF_findradius,PF_bprint,PF_sprint,PF_dprint,
PF_ftos,PF_vtos,PF_coredump,PF_traceon,PF_traceoff,PF_eprint,PF_walkmove,
PF_Fixme,PF_droptofloor,PF_lightstyle,PF_rint,PF_floor,PF_ceil,PF_Fixme,
PF_checkbottom,PF_pointcontents,PF_Fixme,PF_fabs,PF_aim,PF_cvar,PF_localcmd,
PF_nextent,PF_particle,PF_changeyaw,PF_Fixme,PF_vectoangles,PF_WriteByte,
PF_WriteChar,PF_WriteShort,PF_WriteLong,PF_WriteCoord,PF_WriteAngle,
PF_WriteString,PF_WriteEntity,PF_Fixme,PF_Fixme,PF_Fixme,PF_Fixme,PF_Fixme,
PF_Fixme,PF_Fixme,SV_MoveToGoal,PF_precache_file,PF_makestatic,PF_changelevel,
PF_Fixme,PF_cvar_set,PF_centerprint,PF_ambientsound,PF_precache_model,
PF_precache_sound,PF_precache_file,PF_setspawnparms,PF_finalefinished,
PF_localsound,PF_Fixme,PF_Fixme,PF_Fixme,PF_Fixme,PF_Fixme,PF_Fixme,PF_Fixme,
PF_Fixme,PF_Fixme,PF_CheckPlayerEXFlags,PF_walkpathtogoal,PF_Fixme, };
const builtin_t *pr_builtins = pr_builtin;
const s32 pr_numbuiltins = Q_COUNTOF(pr_builtin);
