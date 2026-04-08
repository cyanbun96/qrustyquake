// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// Copyright(C) 2016      Spike
// GPLv3 See LICENSE for details.
#include "quakedef.h"

static s8 pr_string_temp[STRINGTEMP_BUFFERS][STRINGTEMP_LENGTH];
static u8 pr_string_tempindex = 0;
static u8 *checkpvs; //ericw -- changed to malloc
static s32 checkpvs_capacity;
static s32 c_invis, c_notvis;
static struct
{
    char name[MAX_QPATH];
    unsigned int flags;
    qpic_t *pic;
} *qcpics;
static size_t numqcpics;
static size_t maxqcpics;

void PF_Fixme (void) { if(developer.value)PR_RunError ("unimplemented builtin"); }

static qpic_t *DrawQC_CachePic(const char *picname1, unsigned int flags)
{   //okay, so this is silly. we've ended up with 3 different cache levels. qcpics, pics, and images.
    char picname[MAX_OSPATH] = "gfx/";
    if(strncmp(picname, picname1, 4))
	    strncpy(picname+4, picname1, MAX_OSPATH-4-4);
    else
	    strncpy(picname, picname1, MAX_OSPATH-4);
    if(!strstr(picname, ".lmp"))
	    strncpy(picname+strlen(picname), ".lmp", 4);

    size_t i;
    unsigned int texflags;
    for (i = 0; i < numqcpics; i++)
    {   //binary search? something more sane?
        if (!strcmp(picname, qcpics[i].name) || !strcmp(picname1, qcpics[i].name))
        {
            if (qcpics[i].pic)
                return qcpics[i].pic;
            break;
        }
    }

    if (strlen(picname) >= MAX_QPATH)
        return NULL;    //too long. get lost.

    if (flags & PICFLAG_NOLOAD)
        return NULL;    //its a query, not actually needed.

    if (i+1 > maxqcpics)
    {
        maxqcpics = i + 32;
        qcpics = realloc(qcpics, maxqcpics * sizeof(*qcpics));
    }

    strcpy(qcpics[i].name, picname);
    qcpics[i].flags = flags;
    qcpics[i].pic = NULL;

    texflags = TEXPREF_ALPHA | TEXPREF_PAD | TEXPREF_NOPICMIP | TEXPREF_CLAMP | TEXPREF_UNCOMPRESSED;
    if (flags & PICFLAG_WRAP)
        texflags &= ~(TEXPREF_PAD | TEXPREF_CLAMP); //don't allow padding if its going to need to wrap (even if we don't enable clamp-to-edge normally). I just hope we have npot support.
    if (flags & PICFLAG_MIPMAP)
        texflags |= TEXPREF_MIPMAP;

    //try to load it from a wad if applicable.
    //the extra gfx/ crap is because DP insists on it for wad images. and its a nightmare to get things working in all engines if we don't accept that quirk too.
    qcpics[i].pic = Draw_TryCachePic(picname);
    if(!qcpics[i].pic){
	    qcpics[i].pic = Draw_PicFromWad(picname1);
	    strcpy(qcpics[i].name, picname1);
    }

    //okay, not a wad pic, try and load a lmp/tga/etc
    //TODOif (!qcpics[i].pic)
    //TODO    qcpics[i].pic = Draw_TryCachePic(picname, texflags);

    if (i == numqcpics)
        numqcpics++;

    return qcpics[i].pic;
}
void PR_ReloadPics(bool purge)
{
    numqcpics = 0;
    free(qcpics);
    qcpics = NULL;
    maxqcpics = 0;
}

static s8 *PR_GetTempString()
{ return pr_string_temp[(STRINGTEMP_BUFFERS-1) & ++pr_string_tempindex]; }

static const s8* PF_GetStringArg(s32 idx, void* userdata)
{
	if(userdata) idx += *(s32*)userdata;
	if(idx < 0 || idx >= pr_argc) return "";
	return LOC_GetString(G_STRING(OFS_PARM0 + idx * 3));
}

static char *PF_VarString (int  first)
{
    int     i;
    static char out[1024];
    const char *format;
    size_t s;

    out[0] = 0;
    s = 0;

    if (first >= qcvm->argc)
        return out;

    format = LOC_GetString(G_STRING((OFS_PARM0 + first * 3)));
    if (LOC_HasPlaceholders(format))
    {
        int offset = first + 1;
        s = LOC_Format(format, PF_GetStringArg, &offset, out, sizeof(out));
    }
    else
    {
        for (i = first; i < qcvm->argc; i++)
        {
            s = q_strlcat(out, LOC_GetString(G_STRING(OFS_PARM0+i*3)), sizeof(out));
            if (s >= sizeof(out))
            {
                Con_DPrintf("PF_VarString: overflow (string truncated)\n");
                return out;
            }
        }
    }
    if (s > 255)
    {
        /*TODOif (!dev_overflows.varstring || dev_overflows.varstring + CONSOLE_RESPAM_TIME < realtime)
        {
            Con_DPrintf("PF_VarString: %i characters exceeds standard limit of 255 (max = %d).\n",
                                (int) s, (int)(sizeof(out) - 1));
            dev_overflows.varstring = realtime;
        }*/
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
    edict_t *e;
    float   *minvec, *maxvec;

    e = G_EDICT(OFS_PARM0);
    minvec = G_VECTOR(OFS_PARM1);
    maxvec = G_VECTOR(OFS_PARM2);
    SetMinMaxSize (e, minvec, maxvec, false);
}

static void PF_setmodel (void)
{ // setmodel(entity, model)
    int     i;
    const char  *m, **check;
    model_t    *mod;
    edict_t     *e;

    e = G_EDICT(OFS_PARM0);
    m = G_STRING(OFS_PARM1);

// check to see if model was properly precached
    for (i = 0, check = sv.model_precache; *check; i++, check++)
    {
        if (!strcmp(*check, m))
            break;
    }

    if (!*check)
    {
        PR_RunError ("no precache: %s", m);
    }
    e->v.model = PR_SetEngineString(*check);
    e->v.modelindex = i; //SV_ModelIndex (m);

    mod = sv.models[ (int)e->v.modelindex];  // Mod_ForName (m, true);

    if (mod)
    //johnfitz -- correct physics cullboxes for bmodels
    {
        if (mod->type == mod_brush)
            SetMinMaxSize (e, mod->clipmins, mod->clipmaxs, true);
        else
            SetMinMaxSize (e, mod->mins, mod->maxs, true);
    }
    //johnfitz
    else
        SetMinMaxSize (e, vec3_origin, vec3_origin, true);
}

static void PF_bprint() // broadcast print to everyone on server
{ s8 *s = PF_VarString(0); SV_BroadcastPrintf("%s", s); }

/*KRIMZON_SV_PARSECLIENTCOMMAND added these two - note that for compatibility with DP, this tokenize builtin is veeery vauge and doesn't match the console*/
static void PF_Tokenize(void)
{
	puts("TODO PF_Tokenize");
    //G_FLOAT(OFS_RETURN) = tokenizeqc(G_STRING(OFS_PARM0), true);
}

static void PF_tokenize_console(void)
{
	puts("TODO PF_tokenize_console");
    //G_FLOAT(OFS_RETURN) = tokenizeqc(G_STRING(OFS_PARM0), false);
}

static void PF_ArgV(void)
{
	puts("TODO PF_ArgV");
    /*int idx = G_FLOAT(OFS_PARM0);

    //negative indexes are relative to the end
    if (idx < 0)
        idx += qctoken_count;

    if ((unsigned int)idx >= qctoken_count)
        G_INT(OFS_RETURN) = 0;
    else
    {
        char *ret = PR_GetTempString();
        q_strlcpy(ret, qctoken[idx].token, STRINGTEMP_LENGTH);
        G_INT(OFS_RETURN) = PR_SetEngineString(ret);
    }*/
}

static void PF_ArgC(void)
{
	puts("TODO PF_ArgC");
    //G_FLOAT(OFS_RETURN) = qctoken_count;
}

static void PF_sprintf()
{
	puts("TODO PF_sprintf");
}

static void PF_cl_registercommand(void)
{
    const s8 *cmdname = G_STRING(OFS_PARM0);
    Cmd_AddCommand((s8*)cmdname, NULL);
}

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
	else pr_global_struct->trace_ent = EDICT_TO_PROG(qcvm->edicts);
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
	if(qcvm->time - sv.lastchecktime >= 0.1) {
		sv.lastcheck = PF_newcheckclient(sv.lastcheck);
		sv.lastchecktime = qcvm->time;
	}
	edict_t *ent = EDICT_NUM(sv.lastcheck);
	if(ent->free || ent->v.health<=0){ //return check if it might be visible
		RETURN_EDICT(qcvm->edicts);
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
		RETURN_EDICT(qcvm->edicts);
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
	edict_t *chain = (edict_t *)qcvm->edicts;
	f32 *org = G_VECTOR(OFS_PARM0);
	f32 rad = G_FLOAT(OFS_PARM1);
	rad *= rad;
	edict_t *ent = NEXT_EDICT(qcvm->edicts);
	for(s32 i = 1; i < qcvm->num_edicts; i++, ent = NEXT_EDICT(ent)) {
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


static void PF_Find (void)
{ // entity (entity start, .string field, string match) find = #5;
    int     e;
    int     f;
    const char  *s, *t;
    edict_t *ed;

    e = G_EDICTNUM(OFS_PARM0);
    f = G_INT(OFS_PARM1);
    s = G_STRING(OFS_PARM2);
    if (!s)
        PR_RunError ("PF_Find: bad search string");

    for (e++ ; e < qcvm->num_edicts ; e++)
    {
        ed = EDICT_NUM(e);
        if (ed->free)
            continue;
        t = E_STRING(ed,f);
        if (!t)
            continue;
        if (!strcmp(t,s))
        {
            RETURN_EDICT(ed);
            return;
        }
    }

    RETURN_EDICT(qcvm->edicts);
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
		if(i == qcvm->num_edicts) {
			RETURN_EDICT(qcvm->edicts);
			return;
		}
		edict_t *ent = EDICT_NUM(i);
		if(!ent->free) {
			RETURN_EDICT(ent);
			return;
		}
	}
}

static void PF_checkcommand(void)
{
	puts("TODO: PF_checkcommand");
    /*const char *name = G_STRING(OFS_PARM0);
    if (Cmd_Exists(name))
        G_FLOAT(OFS_RETURN) = 1;
    else if (Cmd_AliasExists(name))
        G_FLOAT(OFS_RETURN) = 2;
    else if (Cvar_FindVar(name))
        G_FLOAT(OFS_RETURN) = 3;
    else
        G_FLOAT(OFS_RETURN) = 0;*/
}
static void PF_clientcommand(void)
{
    edict_t *ed             = G_EDICT(OFS_PARM0);
    const char *str         = G_STRING(OFS_PARM1);
    unsigned int i          = NUM_FOR_EDICT(ed)-1;
    if (i < (unsigned int)svs.maxclients && svs.clients[i].active)
    {
        client_t *ohc = host_client;
        host_client = &svs.clients[i];
        Cmd_ExecuteString (str, src_client);
        host_client = ohc;
    }
    else
        Con_Printf("PF_clientcommand: not a client\n");
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
	edict_t *check = NEXT_EDICT(qcvm->edicts);
	for(s32 i = 1; i < qcvm->num_edicts; i++, check = NEXT_EDICT(check)) {
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

static void PF_makestatic (void)
{
    edict_t *ent;
    int     i;
    int bits = 0; //johnfitz -- PROTOCOL_FITZQUAKE

    ent = G_EDICT(OFS_PARM0);

    //johnfitz -- don't send invisible static entities
    if (ent->alpha == ENTALPHA_ZERO) {
        ED_Free (ent);
        return;
    }
    //johnfitz

    //johnfitz -- PROTOCOL_FITZQUAKE
    if (sv.protocol == PROTOCOL_NETQUAKE)
    {
        if (SV_ModelIndex(PR_GetString(ent->v.model)) & 0xFF00 || (int)(ent->v.frame) & 0xFF00)
        {
            ED_Free (ent);
            return; //can't display the correct model & frame, so don't show it at all
        }
    }
    else
    {
        if (SV_ModelIndex(PR_GetString(ent->v.model)) & 0xFF00)
            bits |= B_LARGEMODEL;
        if ((int)(ent->v.frame) & 0xFF00)
            bits |= B_LARGEFRAME;
        if (ent->alpha != ENTALPHA_DEFAULT)
            bits |= B_ALPHA;

        if (sv.protocol == PROTOCOL_RMQ)
        {
            eval_t* val;
            val = GetEdictFieldValueByName(ent, "scale");
            if (val)
                ent->scale = ENTSCALE_ENCODE(val->_float);
            else
                ent->scale = ENTSCALE_DEFAULT;

            if (ent->scale != ENTSCALE_DEFAULT)
                bits |= B_SCALE;
        }
    }

    SV_ReserveSignonSpace (34);

    if (bits)
    {
        MSG_WriteByte (sv.signon, svc_spawnstatic2);
        MSG_WriteByte (sv.signon, bits);
    }
    else
        MSG_WriteByte (sv.signon, svc_spawnstatic);

    if (bits & B_LARGEMODEL)
        MSG_WriteShort (sv.signon, SV_ModelIndex(PR_GetString(ent->v.model)));
    else
        MSG_WriteByte (sv.signon, SV_ModelIndex(PR_GetString(ent->v.model)));

    if (bits & B_LARGEFRAME)
        MSG_WriteShort (sv.signon, ent->v.frame);
    else
        MSG_WriteByte (sv.signon, ent->v.frame);
    //johnfitz

    MSG_WriteByte (sv.signon, ent->v.colormap);
    MSG_WriteByte (sv.signon, ent->v.skin);
    for (i = 0; i < 3; i++)
    {
        MSG_WriteCoord(sv.signon, ent->v.origin[i], sv.protocolflags);
        MSG_WriteAngle(sv.signon, ent->v.angles[i], sv.protocolflags);
    }

    //johnfitz -- PROTOCOL_FITZQUAKE
    if (bits & B_ALPHA)
        MSG_WriteByte (sv.signon, ent->alpha);
    //johnfitz

    if (bits & B_SCALE)
        MSG_WriteByte (sv.signon, ent->scale);

// throw the entity away now
    ED_Free (ent);
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

static void PF_Sin(void)
{ G_FLOAT(OFS_RETURN) = sin(G_FLOAT(OFS_PARM0)); }
static void PF_asin(void)
{ G_FLOAT(OFS_RETURN) = asin(G_FLOAT(OFS_PARM0)); }
static void PF_Cos(void)
{ G_FLOAT(OFS_RETURN) = cos(G_FLOAT(OFS_PARM0)); }
static void PF_acos(void)
{ G_FLOAT(OFS_RETURN) = acos(G_FLOAT(OFS_PARM0)); }
static void PF_tan(void)
{ G_FLOAT(OFS_RETURN) = tan(G_FLOAT(OFS_PARM0)); }
static void PF_atan(void)
{ G_FLOAT(OFS_RETURN) = atan(G_FLOAT(OFS_PARM0)); }
static void PF_atan2(void)
{ G_FLOAT(OFS_RETURN) = atan2(G_FLOAT(OFS_PARM0), G_FLOAT(OFS_PARM1)); }
static void PF_Sqrt(void)
{ G_FLOAT(OFS_RETURN) = sqrt(G_FLOAT(OFS_PARM0)); }
static void PF_pow(void)
{ G_FLOAT(OFS_RETURN) = pow(G_FLOAT(OFS_PARM0), G_FLOAT(OFS_PARM1)); }

static void PF_stof(void)
{ G_FLOAT(OFS_RETURN) = atof(G_STRING(OFS_PARM0)); }
static void PF_stov(void)
{
    const char *s = G_STRING(OFS_PARM0);
    s = COM_Parse(s);
    G_VECTOR(OFS_RETURN)[0] = atof(com_token);
    s = COM_Parse(s);
    G_VECTOR(OFS_RETURN)[1] = atof(com_token);
    s = COM_Parse(s);
    G_VECTOR(OFS_RETURN)[2] = atof(com_token);
}
static void PF_etos(void)
{   //yes, this is lame
    char *result = PR_GetTempString();
    q_snprintf(result, STRINGTEMP_LENGTH, "entity %i", G_EDICTNUM(OFS_PARM0));
    G_INT(OFS_RETURN) = PR_SetEngineString(result);
}
static void PF_ftoi(void)
{ G_INT(OFS_RETURN) = G_FLOAT(OFS_PARM0); }
static void PF_itof(void)
{ G_FLOAT(OFS_RETURN) = G_INT(OFS_PARM0); }

static void PF_mod(void)
{
    float a = G_FLOAT(OFS_PARM0);
    float n = G_FLOAT(OFS_PARM1);

    if (n == 0)
    {
        Con_DPrintf("PF_mod: mod by zero\n");
        G_FLOAT(OFS_RETURN) = 0;
    }
    else
    {
        //because QC is inherantly floaty, lets use floats.
        G_FLOAT(OFS_RETURN) = a - (n * (int)(a/n));
    }
}
static void PF_min(void)
{
    float r = G_FLOAT(OFS_PARM0);
    int i;
    for (i = 1; i < qcvm->argc; i++)
    {
        if (r > G_FLOAT(OFS_PARM0 + i*3))
            r = G_FLOAT(OFS_PARM0 + i*3);
    }
    G_FLOAT(OFS_RETURN) = r;
}
static void PF_max(void)
{
    float r = G_FLOAT(OFS_PARM0);
    int i;
    for (i = 1; i < qcvm->argc; i++)
    {
        if (r < G_FLOAT(OFS_PARM0 + i*3))
            r = G_FLOAT(OFS_PARM0 + i*3);
    }
    G_FLOAT(OFS_RETURN) = r;
}
static void PF_bound(void)
{
    float minval = G_FLOAT(OFS_PARM0);
    float curval = G_FLOAT(OFS_PARM1);
    float maxval = G_FLOAT(OFS_PARM2);
    if (curval > maxval)
        curval = maxval;
    if (curval < minval)
        curval = minval;
    G_FLOAT(OFS_RETURN) = curval;
}

static void PF_vectorvectors(void)
{
    VectorCopy(G_VECTOR(OFS_PARM0), pr_global_struct->v_forward);
    VectorNormalize(pr_global_struct->v_forward);
    if (!pr_global_struct->v_forward[0] && !pr_global_struct->v_forward[1])
    {
        if (pr_global_struct->v_forward[2])
            pr_global_struct->v_right[1] = -1;
        else
            pr_global_struct->v_right[1] = 0;
        pr_global_struct->v_right[0] = pr_global_struct->v_right[2] = 0;
    }
    else
    {
        pr_global_struct->v_right[0] = pr_global_struct->v_forward[1];
        pr_global_struct->v_right[1] = -pr_global_struct->v_forward[0];
        pr_global_struct->v_right[2] = 0;
        VectorNormalize(pr_global_struct->v_right);
    }
    CrossProduct(pr_global_struct->v_right, pr_global_struct->v_forward, pr_global_struct->v_up);
}

static void PF_checkextension(void)
{
    const char *extname = G_STRING(OFS_PARM0);
    int i = PR_FindExtensionByName (extname);
    if (i)
        SetBit (qcvm->checked_ext, i);

    // Note: we expose FTE_QC_CHECKCOMMAND so that AD considers the engine
    // FTE-like instead of DP-like, in order to avoid a bug in the DP codepath
    // in older AD versions (e.g. 1.42, used in jam8)
    if (i == FTE_QC_CHECKCOMMAND)
    {
        G_FLOAT(OFS_RETURN) = true;
        SetBit (qcvm->advertised_ext, i);
    }
    else
        G_FLOAT(OFS_RETURN) = false;
}

static void PF_strlen(void)
{   //FIXME: doesn't try to handle utf-8
    const char *s = G_STRING(OFS_PARM0);
    G_FLOAT(OFS_RETURN) = strlen(s);
}
static void PF_strcat(void)
{
    int     i;
    char *out = PR_GetTempString();
    size_t s;

    out[0] = 0;
    s = 0;
    for (i = 0; i < qcvm->argc; i++)
    {
        s = q_strlcat(out, G_STRING((OFS_PARM0+i*3)), STRINGTEMP_LENGTH);
        if (s >= STRINGTEMP_LENGTH)
        {
            Con_DPrintf("PF_strcat: overflow (string truncated)\n");
            break;
        }
    }

    G_INT(OFS_RETURN) = PR_SetEngineString(out);
}
static void PF_substring(void)
{
    int start, length, slen;
    const char *s;
    char *string;

    s = G_STRING(OFS_PARM0);
    start = G_FLOAT(OFS_PARM1);
    length = G_FLOAT(OFS_PARM2);

    slen = strlen(s);   //utf-8 should use chars, not bytes.

    if (start < 0)
        start = slen+start;
    if (length < 0)
        length = slen-start+(length+1);
    if (start < 0)
    {
    //  length += start;
        start = 0;
    }

    if (start >= slen || length<=0)
    {
        G_INT(OFS_RETURN) = PR_SetEngineString("");
        return;
    }

    slen -= start;
    if (length > slen)
        length = slen;
    //utf-8 should switch to bytes now.
    s += start;

    if (length >= STRINGTEMP_LENGTH)
    {
        length = STRINGTEMP_LENGTH-1;
        Con_DPrintf("PF_substring: truncation\n");
    }

    string = PR_GetTempString();
    memcpy(string, s, length);
    string[length] = '\0';
    G_INT(OFS_RETURN) = PR_SetEngineString(string);
}
/*our zoned strings implementation is somewhat specific to quakespasm, so good luck porting*/
static void PF_strzone(void)
{
    char *buf;
    size_t len = 0;
    const char *s[8];
    size_t l[8];
    int i;
    size_t id;

    for (i = 0; i < qcvm->argc; i++)
    {
        s[i] = G_STRING(OFS_PARM0+i*3);
        l[i] = strlen(s[i]);
        len += l[i];
    }
    len++; /*for the null*/

    buf = Z_Malloc(len);
    G_INT(OFS_RETURN) = PR_SetEngineString(buf);
    id = -1-G_INT(OFS_RETURN);
    if (id >= qcvm->knownzonesize)
    {
        qcvm->knownzonesize = (id+32)&~7;
        qcvm->knownzone = Z_Realloc(qcvm->knownzone, (qcvm->knownzonesize+7)>>3);
    }
    qcvm->knownzone[id>>3] |= 1u<<(id&7);

    for (i = 0; i < qcvm->argc; i++)
    {
        memcpy(buf, s[i], l[i]);
        buf += l[i];
    }
    *buf = '\0';
}
static void PF_strunzone(void)
{
	puts("TODO PF_strunzone");
	/*
    size_t id;
    const char *foo = G_STRING(OFS_PARM0);

    if (!G_INT(OFS_PARM0))
        return; //don't bug out if they gave a null string
    id = -1-G_INT(OFS_PARM0);
    if (id < qcvm->knownzonesize && (qcvm->knownzone[id>>3] & (1u<<(id&7))))
    {
        qcvm->knownzone[id>>3] &= ~(1u<<(id&7));
        PR_ClearEngineString(G_INT(OFS_PARM0));
        Z_Free((void*)foo);
    }
    else
        Con_DPrintf("PF_strunzone: string wasn't strzoned\n");*/
}

static bool qc_isascii(unsigned int u)
{
    if (u < 256)    //should be just \n and 32-127, but we don't actually support any actual unicode and we don't really want to make things worse.
        return true;
    return false;
}
static void PF_str2chr(void)
{
    const char *instr = G_STRING(OFS_PARM0);
    int ofs = (qcvm->argc>1)?G_FLOAT(OFS_PARM1):0;

    if (ofs < 0)
        ofs = strlen(instr)+ofs;

    if (ofs && (ofs < 0 || ofs > (int)strlen(instr)))
        G_FLOAT(OFS_RETURN) = '\0';
    else
        G_FLOAT(OFS_RETURN) = (unsigned char)instr[ofs];
}
static void PF_chr2str(void)
{
    char *ret = PR_GetTempString(), *out;
    int i;
    for (i = 0, out=ret; out-ret < STRINGTEMP_LENGTH-6 && i < qcvm->argc; i++)
    {
        unsigned int u = G_FLOAT(OFS_PARM0 + i*3);
        if (u >= 0xe000 && u < 0xe100)
            *out++ = (unsigned char)u;  //quake chars.
        else if (qc_isascii(u))
            *out++ = u;
        else
            *out++ = '?';   //no unicode support
    }
    *out = 0;
    G_INT(OFS_RETURN) = PR_SetEngineString(ret);
}
static void PF_cl_drawfill(void)
{
    float *pos  = G_VECTOR(OFS_PARM0);
    float *size = G_VECTOR(OFS_PARM1);
    float *rgb  = G_VECTOR(OFS_PARM2);
    float alpha = G_FLOAT (OFS_PARM3);
//  int flags   = G_FLOAT (OFS_PARM4);

    puts("TODO PF_cl_drawfill");
    //Draw_FillEx (pos[0], pos[1], size[0], size[1], rgb, alpha);
}

static void PF_cl_drawpic(void)
{
    float *pos  = G_VECTOR(OFS_PARM0);
    qpic_t *pic = DrawQC_CachePic(G_STRING(OFS_PARM1), PICFLAG_AUTO);
    float *size = G_VECTOR(OFS_PARM2);
    float *rgb  = G_VECTOR(OFS_PARM3);
    float alpha = G_FLOAT (OFS_PARM4);
//  int flags   = G_FLOAT (OFS_PARM5);

    if (pic)
        Draw_TransPicScaled(pos[0], pos[1], pic, uiscale);
}

static void PF_cl_getimagesize(void)
{
    qpic_t *pic = DrawQC_CachePic(G_STRING(OFS_PARM0), PICFLAG_AUTO);
    if (pic)
        G_VECTORSET(OFS_RETURN, pic->width, pic->height, 0);
    else
        G_VECTORSET(OFS_RETURN, 0, 0, 0);
}
static void PF_cl_drawcharacter(void)
{
    float *pos  = G_VECTOR(OFS_PARM0);
    int charcode= (int)G_FLOAT (OFS_PARM1) & 0xff;
//  float *size = G_VECTOR(OFS_PARM2);
//  float *rgb  = G_VECTOR(OFS_PARM3);
//  float alpha = G_FLOAT (OFS_PARM4);
//  int flags   = G_FLOAT (OFS_PARM5);

    if (charcode == 32)
        return; //don't waste time on spaces

    Draw_CharacterScaled(pos[0], pos[1], charcode, uiscale);
}
static void PF_cl_drawrawstring(void)
{
    float *pos  = G_VECTOR(OFS_PARM0);
    const char *text = G_STRING (OFS_PARM1);
//  float *size = G_VECTOR(OFS_PARM2);
//  float *rgb  = G_VECTOR(OFS_PARM3);
//  float alpha = G_FLOAT (OFS_PARM4);
//  int flags   = G_FLOAT (OFS_PARM5);

    float x = pos[0];
    int c;

    if (!*text)
        return; //don't waste time on spaces

    while ((c = *text++))
    {
        Draw_CharacterScaled (x, pos[1], c, uiscale);
        x += 8*uiscale;//size[0];
    }
}
static void PF_cl_drawstring(void)
{
    float *pos  = G_VECTOR(OFS_PARM0);
    const char *text = G_STRING (OFS_PARM1);
    //float *size = G_VECTOR(OFS_PARM2);
    //float *rgb  = G_VECTOR(OFS_PARM3);
    //float alpha = G_FLOAT (OFS_PARM4);
//  int flags   = G_FLOAT (OFS_PARM5);

    //float x = pos[0];
    //struct markup_s mu;
    //int c;

    if (!*text)
        return; //don't waste time on spaces

    Draw_StringScaled(pos[0], pos[1], text, uiscale);
    //PR_Markup_Begin(&mu, text, rgb, alpha);

    //while ((c = PR_Markup_Parse(&mu)))
    //{
    //    GL_SetCanvasColor (mu.colour[0], mu.colour[1], mu.colour[2], mu.colour[3]);
    //    DrawQC_CharacterQuad (x, pos[1], c, size[0], size[1]);
    //    x += size[0];
    //}
    //GL_SetCanvasColor (1.f, 1.f, 1.f, 1.f);
}
static void PF_cl_precachepic(void)
{
    const char *name    = G_STRING(OFS_PARM0);
    unsigned int flags = G_FLOAT(OFS_PARM1);

    G_INT(OFS_RETURN) = G_INT(OFS_PARM0);   //return input string, for convienience

    if (!DrawQC_CachePic(name, flags) && (flags & PICFLAG_BLOCK))
        G_INT(OFS_RETURN) = 0;  //return input string, for convienience
}
static void PF_cl_iscachedpic(void)
{
    const char *name    = G_STRING(OFS_PARM0);
    if (DrawQC_CachePic(name, PICFLAG_NOLOAD))
        G_FLOAT(OFS_RETURN) = true;
    else
        G_FLOAT(OFS_RETURN) = false;
}
static void PF_cl_drawsetclip(void)
{
	puts("TODO PF_cl_drawsetclip");
    float x = G_FLOAT(OFS_PARM0);
    float y = G_FLOAT(OFS_PARM1);
    float w = G_FLOAT(OFS_PARM2);
    float h = G_FLOAT(OFS_PARM3);
    //Draw_SetClipRect (x, y, w, h);
}
static void PF_cl_drawresetclip(void)
{
	puts("TODO PF_cl_drawresetclip");
    //Draw_ResetClipping ();
}

static void PF_cl_stringwidth(void)
{
	puts("TODO PF_cl_stringwidth");
    /*static const float defaultfontsize[] = {8,8};
    const char *text = G_STRING (OFS_PARM0);
    bool usecolours = G_FLOAT(OFS_PARM1);
    const float *fontsize = (qcvm->argc>2)?G_VECTOR (OFS_PARM2):defaultfontsize;
    struct markup_s mu;
    int r = 0;

    if (!usecolours)
        r = strlen(text);
    else
    {
        PR_Markup_Begin(&mu, text, vec3_origin, 1);
        while (PR_Markup_Parse(&mu))
        {
            r += 1;
        }
    }

    //primitive and lame, but hey.
    G_FLOAT(OFS_RETURN) = fontsize[0] * r;*/
}
static void PF_cl_drawsubpic(void)
{
    float *pos  = G_VECTOR(OFS_PARM0);
    float *size = G_VECTOR(OFS_PARM1);
    qpic_t *pic = DrawQC_CachePic(G_STRING(OFS_PARM2), PICFLAG_AUTO);
    float *srcpos   = G_VECTOR(OFS_PARM3);
    float *srcsize  = G_VECTOR(OFS_PARM4);
    float *rgb  = G_VECTOR(OFS_PARM5);
    float alpha = G_FLOAT (OFS_PARM6);
//  int flags   = G_FLOAT (OFS_PARM7);

    if (pic)
        Draw_TransPicScaled(pos[0], pos[1], pic, uiscale);
}
static void PF_cl_getstat_int(void)
{
    int stnum = G_FLOAT(OFS_PARM0);
    if (stnum < 0 || stnum >= Q_COUNTOF(cl.stats))
        G_INT(OFS_RETURN) = 0;
    else
        G_INT(OFS_RETURN) = cl.stats[stnum];
}
static void PF_cl_getstat_float(void)
{
    int stnum = G_FLOAT(OFS_PARM0);
    if (stnum < 0 || stnum >= Q_COUNTOF(cl.stats))
        G_FLOAT(OFS_RETURN) = 0;
    /*else if (qcvm->argc > 1)
    {*/
        int firstbit = G_FLOAT(OFS_PARM1);
        int bitcount = G_FLOAT(OFS_PARM2);
        G_FLOAT(OFS_RETURN) = (cl.stats[stnum]>>firstbit) & ((1<<bitcount)-1);
    /*}
    else
        G_FLOAT(OFS_RETURN) = cl.statsf[stnum];*/
}
static void PF_cl_getstat_string(void)
{
	puts("PF_cl_getstat_string");
    /*int stnum = G_FLOAT(OFS_PARM0);
    if (stnum < 0 || stnum >= countof(cl.statss) || !cl.statss[stnum])
        G_INT(OFS_RETURN) = 0;
    else
    {
        char *result = PR_GetTempString();
        q_strlcpy(result, cl.statss[stnum], STRINGTEMP_LENGTH);
        G_INT(OFS_RETURN) = PR_SetEngineString(result);
    }*/
}
static void PF_cl_playerkey_s(void)
{
	puts("PF_cl_playerkey_s");
    int playernum = G_FLOAT(OFS_PARM0);
    const char *keyname = G_STRING(OFS_PARM1);
    //PF_cl_playerkey_internal(playernum, keyname, false);
}
static void PF_cl_playerkey_f(void)
{
	puts("PF_cl_playerkey_f");
    int playernum = G_FLOAT(OFS_PARM0);
    const char *keyname = G_STRING(OFS_PARM1);
    //PF_cl_playerkey_internal(playernum, keyname, true);
}

//part of PF_strconv
static int chrconv_number(int i, int base, int conv)
{
    i -= base;
    switch (conv)
    {
    default:
    case 5:
    case 6:
    case 0:
        break;
    case 1:
        base = '0';
        break;
    case 2:
        base = '0'+128;
        break;
    case 3:
        base = '0'-30;
        break;
    case 4:
        base = '0'+128-30;
        break;
    }
    return i + base;
}
//part of PF_strconv
static int chrconv_punct(int i, int base, int conv)
{
    i -= base;
    switch (conv)
    {
    default:
    case 0:
        break;
    case 1:
        base = 0;
        break;
    case 2:
        base = 128;
        break;
    }
    return i + base;
}
//part of PF_strconv
static int chrchar_alpha(int i, int basec, int baset, int convc, int convt, int charnum)
{
    //convert case and colour seperatly...

    i -= baset + basec;
    switch (convt)
    {
    default:
    case 0:
        break;
    case 1:
        baset = 0;
        break;
    case 2:
        baset = 128;
        break;

    case 5:
    case 6:
        baset = 128*((charnum&1) == (convt-5));
        break;
    }

    switch (convc)
    {
    default:
    case 0:
        break;
    case 1:
        basec = 'a';
        break;
    case 2:
        basec = 'A';
        break;
    }
    return i + basec + baset;
}
//FTE_STRINGS
//bulk convert a string. change case or colouring.
static void PF_strconv (void)
{
    int ccase = G_FLOAT(OFS_PARM0);     //0 same, 1 lower, 2 upper
    int redalpha = G_FLOAT(OFS_PARM1);  //0 same, 1 white, 2 red,  5 alternate, 6 alternate-alternate
    int rednum = G_FLOAT(OFS_PARM2);    //0 same, 1 white, 2 red, 3 redspecial, 4 whitespecial, 5 alternate, 6 alternate-alternate
    const unsigned char *string = (const unsigned char*)PF_VarString(3);
    int len = strlen((const char*)string);
    int i;
    unsigned char *resbuf = (unsigned char*)PR_GetTempString();
    unsigned char *result = resbuf;

    //UTF-8-FIXME: cope with utf+^U etc

    if (len >= STRINGTEMP_LENGTH)
        len = STRINGTEMP_LENGTH-1;

    for (i = 0; i < len; i++, string++, result++)   //should this be done backwards?
    {
        if (*string >= '0' && *string <= '9')   //normal numbers...
            *result = chrconv_number(*string, '0', rednum);
        else if (*string >= '0'+128 && *string <= '9'+128)
            *result = chrconv_number(*string, '0'+128, rednum);
        else if (*string >= '0'+128-30 && *string <= '9'+128-30)
            *result = chrconv_number(*string, '0'+128-30, rednum);
        else if (*string >= '0'-30 && *string <= '9'-30)
            *result = chrconv_number(*string, '0'-30, rednum);

        else if (*string >= 'a' && *string <= 'z')  //normal numbers...
            *result = chrchar_alpha(*string, 'a', 0, ccase, redalpha, i);
        else if (*string >= 'A' && *string <= 'Z')  //normal numbers...
            *result = chrchar_alpha(*string, 'A', 0, ccase, redalpha, i);
        else if (*string >= 'a'+128 && *string <= 'z'+128)  //normal numbers...
            *result = chrchar_alpha(*string, 'a', 128, ccase, redalpha, i);
        else if (*string >= 'A'+128 && *string <= 'Z'+128)  //normal numbers...
            *result = chrchar_alpha(*string, 'A', 128, ccase, redalpha, i);

        else if ((*string & 127) < 16 || !redalpha) //special chars..
            *result = *string;
        else if (*string < 128)
            *result = chrconv_punct(*string, 0, redalpha);
        else
            *result = chrconv_punct(*string, 128, redalpha);
    }
    *result = '\0';

    G_INT(OFS_RETURN) = PR_SetEngineString((char*)resbuf);
}
static struct svcustomstat_s *PR_CustomStat(int idx, int type)
{
    size_t i;
    if (idx < 0 || idx >= MAX_CL_STATS)
        return NULL;
    switch(type)
    {
    case ev_ext_integer:
    case ev_float:
    case ev_vector:
    case ev_entity:
        break;
    default:
        return NULL;
    }

    for (i = 0; i < sv.numcustomstats; i++)
    {
        if (sv.customstats[i].idx == idx && (sv.customstats[i].type==ev_string) == (type==ev_string))
            break;
    }
    if (i == sv.numcustomstats)
        sv.numcustomstats++;
    sv.customstats[i].idx = idx;
    sv.customstats[i].type = type;
    sv.customstats[i].fld = 0;
    sv.customstats[i].ptr = NULL;
    return &sv.customstats[i];
}
static void PF_clientstat(void)
{
    int idx = G_FLOAT(OFS_PARM0);
    int type = G_FLOAT(OFS_PARM1);
    int fldofs = G_INT(OFS_PARM2);
    struct svcustomstat_s *stat = PR_CustomStat(idx, type);
    if (!stat)
        return;
    stat->fld = fldofs;
}

#define PF_BOTH(x)  x,x
#define PF_CSQC(x)  NULL,x
#define PF_SSQC(x)  x,NULL

builtindef_t pr_builtindefs[] =
{
    {"makevectors",             PF_SSQC(PF_makevectors),        1, 0},     // void(entity e) makevectors       = #1
    {"setorigin",               PF_SSQC(PF_setorigin),          2, 0},     // void(entity e, vector o) setorigin   = #2
    {"setmodel",                PF_SSQC(PF_setmodel),           3, 0},     // void(entity e, string m) setmodel    = #3
    {"setsize",                 PF_SSQC(PF_setsize),            4, 0},     // void(entity e, vector min, vector max) setsize   = #4
    {"break",                   PF_BOTH(PF_break),              6, 0},     // void() break             = #6
    {"random",                  PF_BOTH(PF_random),             7, 0},     // float() random           = #7
    {"sound",                   PF_SSQC(PF_sound),              8, 0},     // void(entity e, float chan, string samp) sound    = #8
    {"normalize",               PF_BOTH(PF_normalize),          9, 0},     // vector(vector v) normalize       = #9
    {"error",                   PF_SSQC(PF_error),              10, 0},    // void(string e) error         = #10
    {"objerror",                PF_SSQC(PF_objerror),           11, 0},    // void(string e) objerror      = #11
    {"vlen",                    PF_BOTH(PF_vlen),               12, 0},    // float(vector v) vlen         = #12
    {"vectoyaw",                PF_BOTH(PF_vectoyaw),           13, 0},    // float(vector v) vectoyaw     = #13
    {"spawn",                   PF_SSQC(PF_Spawn),              14, 0},    // entity() spawn           = #14
    {"remove",                  PF_SSQC(PF_Remove),             15, 0},    // void(entity e) remove        = #15
    {"traceline",               PF_SSQC(PF_traceline),          16, 0},    // float(vector v1, vector v2, float tryents) traceline = #16
    {"checkclient",             PF_SSQC(PF_checkclient),        17, 0},    // entity() clientlist          = #17
    {"find",                    PF_SSQC(PF_Find),               18, 0},    // entity(entity start, .string fld, string match) find = #18
    {"precache_sound",          PF_SSQC(PF_precache_sound),     19, 0},    // void(string s) precache_sound    = #19
    {"precache_model",          PF_SSQC(PF_precache_model),     20, 0},    // void(string s) precache_model    = #20
    {"stuffcmd",                PF_SSQC(PF_stuffcmd),           21, 0},    // void(entity client, string s)stuffcmd    = #21
    {"findradius",              PF_SSQC(PF_findradius),         22, 0},    // entity(vector org, float rad) findradius = #22
    {"bprint",                  PF_SSQC(PF_bprint),             23, 0},    // void(string s) bprint        = #23
    {"sprint",                  PF_SSQC(PF_sprint),             24, 0},    // void(entity client, string s) sprint = #24
    {"dprint",                  PF_BOTH(PF_dprint),             25, 0},    // void(string s) dprint        = #25
    {"ftos",                    PF_BOTH(PF_ftos),               26, 0},    // void(string s) ftos          = #26
    {"vtos",                    PF_BOTH(PF_vtos),               27, 0},    // void(string s) vtos          = #27
    {"coredump",                PF_SSQC(PF_coredump),           28, 0},
    {"traceon",                 PF_BOTH(PF_traceon),            29, 0},
    {"traceoff",                PF_BOTH(PF_traceoff),           30, 0},
    {"eprint",                  PF_SSQC(PF_eprint),             31, 0},    // void(entity e) debug print an entire entity
    {"walkmove",                PF_SSQC(PF_walkmove),           32, 0},    // float(float yaw, float dist) walkmove
    {"droptofloor",             PF_SSQC(PF_droptofloor),        34, 0},
    {"lightstyle",              PF_SSQC(PF_lightstyle),         35, 0},
    {"rint",                    PF_BOTH(PF_rint),               36, 0},
    {"floor",                   PF_BOTH(PF_floor),              37, 0},
    {"ceil",                    PF_BOTH(PF_ceil),               38, 0},
    {"checkbottom",             PF_SSQC(PF_checkbottom),        40, 0},
    {"pointcontents",           PF_SSQC(PF_pointcontents),      41, 0},
    {"fabs",                    PF_BOTH(PF_fabs),               43, 0},
    {"aim",                     PF_SSQC(PF_aim),                44, 0},
    {"cvar",                    PF_BOTH(PF_cvar),               45, 0},
    {"localcmd",                PF_BOTH(PF_localcmd),           46, 0},
    {"nextent",                 PF_SSQC(PF_nextent),            47, 0},
    {"particle",                PF_SSQC(PF_particle),           48, 0},
    {"ChangeYaw",               PF_SSQC(PF_changeyaw),          49, 0},
    {"vectoangles",             PF_BOTH(PF_vectoangles),        51, 0},

    {"WriteByte",               PF_SSQC(PF_WriteByte),          52, 0},
    {"WriteChar",               PF_SSQC(PF_WriteChar),          53, 0},
    {"WriteShort",              PF_SSQC(PF_WriteShort),         54, 0},
    {"WriteLong",               PF_SSQC(PF_WriteLong),          55, 0},
    {"WriteCoord",              PF_SSQC(PF_WriteCoord),         56, 0},
    {"WriteAngle",              PF_SSQC(PF_WriteAngle),         57, 0},
    {"WriteString",             PF_SSQC(PF_WriteString),        58, 0},
    {"WriteEntity",             PF_SSQC(PF_WriteEntity),        59, 0},

    {"sin",                     PF_BOTH(PF_Sin),                60,     DP_QC_SINCOSSQRTPOW},   // float(float angle)
    {"cos",                     PF_BOTH(PF_Cos),                61,     DP_QC_SINCOSSQRTPOW},   // float(float angle)
    {"sqrt",                    PF_BOTH(PF_Sqrt),               62,     DP_QC_SINCOSSQRTPOW},   // float(float value)

    {"etos",                    PF_BOTH(PF_etos),               65,     DP_QC_ETOS},            // string(entity ent)

    {"movetogoal",              PF_SSQC(SV_MoveToGoal),         67, 0},
    {"precache_file",           PF_SSQC(PF_precache_file),      68, 0},
    {"makestatic",              PF_SSQC(PF_makestatic),         69, 0},

    {"changelevel",             PF_SSQC(PF_changelevel),        70, 0},

    {"cvar_set",                PF_BOTH(PF_cvar_set),           72, 0},
    {"centerprint",             PF_SSQC(PF_centerprint),        73, 0},

    {"ambientsound",            PF_SSQC(PF_ambientsound),       74, 0},

    {"precache_model2",         PF_SSQC(PF_precache_model),     75, 0},
    {"precache_sound2",         PF_SSQC(PF_precache_sound),     76, 0},    // precache_sound2 is different only for qcc
    {"precache_file2",          PF_SSQC(PF_precache_file),      77, 0},

    {"setspawnparms",           PF_SSQC(PF_setspawnparms),      78, 0},

    // 2021 re-release
    {"finaleFinished",          PF_SSQC(PF_finalefinished),     79, 0},    // float() finaleFinished = #79
    {"localsound",              PF_SSQC(PF_localsound),         80, 0},    // void localsound (entity client, string sample) = #80

    {"stof",                    PF_BOTH(PF_stof),               81,     FRIK_FILE},         // float(string)

    // 2021 re-release update 3
    {"ex_centerprint",          PF_SSQC(PF_centerprint), 0, 0},               // void(entity client, string s, ...)
    {"ex_bprint",               PF_SSQC(PF_bprint), 0, 0},                    // void(string s, ...)
    {"ex_sprint",               PF_SSQC(PF_sprint), 0, 0},                    // void(entity client, string s, ...)
    {"ex_finalefinished",       PF_SSQC(PF_finalefinished), 0, 0},            // float()
    {"ex_CheckPlayerEXFlags",   PF_SSQC(PF_CheckPlayerEXFlags), 0, 0},        // float(entity playerEnt)
    {"ex_walkpathtogoal",       PF_SSQC(PF_walkpathtogoal), 0, 0},            // float(float movedist, vector goal)
    {"ex_localsound",           PF_SSQC(PF_localsound), 0, 0},                // void(entity client, string sample)

    {"min",                     PF_BOTH(PF_min),                94,     DP_QC_MINMAXBOUND}, // float(float a, float b, ...)
    {"max",                     PF_BOTH(PF_max),                95,     DP_QC_MINMAXBOUND}, // float(float a, float b, ...)
    {"bound",                   PF_BOTH(PF_bound),              96,     DP_QC_MINMAXBOUND}, // float(float minimum, float val, float maximum)

    {"pow",                     PF_BOTH(PF_pow),                97,     DP_QC_SINCOSSQRTPOW},   // float(float value, float exp)

    {"checkextension",          PF_BOTH(PF_checkextension),     99, 0},    // float(string extname)

    {"strlen",                  PF_BOTH(PF_strlen),             114,    FRIK_FILE}, // float(string s)
    {"strcat",                  PF_BOTH(PF_strcat),             115,    FRIK_FILE}, // string(string s1, optional string s2, optional string s3, optional string s4, optional string s5, optional string s6, optional string s7, optional string s8)
    {"substring",               PF_BOTH(PF_substring),          116,    FRIK_FILE}, // string(string s, float start, float length)
    {"stov",                    PF_BOTH(PF_stov),               117,    FRIK_FILE}, // vector(string s)
    {"strzone",                 PF_BOTH(PF_strzone),            118,    FRIK_FILE}, // string(string s, ...)
    {"strunzone",               PF_BOTH(PF_strunzone),          119,    FRIK_FILE}, // void(string s)

    {"str2chr",                 PF_BOTH(PF_str2chr),            222,    FTE_STRINGS},   // float(string str, float index)
    {"chr2str",                 PF_BOTH(PF_chr2str),            223,    FTE_STRINGS},   // string(float chr, ...)
    {"strconv",                 PF_BOTH(PF_strconv),            224,    FTE_STRINGS},   // string(float ccase, float redalpha, float redchars, string str, ...)

    {"clientstat",              PF_SSQC(PF_clientstat),         232, 0},   // void(float num, float type, .__variant fld)

    {"mod",                     PF_BOTH(PF_mod),                245, 0},   // float(float a, float n)
									//     {"ftoi",                    PF_BOTH(PF_ftoi)},                      // int(float)
    {"itof",                    PF_BOTH(PF_itof), 0, 0},                      // float(int)

    {"checkcommand",            PF_BOTH(PF_checkcommand),       294,    FTE_QC_CHECKCOMMAND},   // float(string name)

    {"iscachedpic",             PF_CSQC(PF_cl_iscachedpic),     316, 0},   // float(string name)
    {"precache_pic",            PF_CSQC(PF_cl_precachepic),     317, 0},   // string(string name, optional float flags)
    {"drawgetimagesize",        PF_CSQC(PF_cl_getimagesize),    318, 0},   // #define draw_getimagesize drawgetimagesize\nvector(string picname)
    {"drawcharacter",           PF_CSQC(PF_cl_drawcharacter),   320, 0},   // float(vector position, float character, vector size, vector rgb, float alpha, optional float drawflag)
    {"drawrawstring",           PF_CSQC(PF_cl_drawrawstring),   321, 0},   // float(vector position, string text, vector size, vector rgb, float alpha, optional float drawflag)
    {"drawpic",                 PF_CSQC(PF_cl_drawpic),         322, 0},   // float(vector position, string pic, vector size, vector rgb, float alpha, optional float drawflag)
    {"drawfill",                PF_CSQC(PF_cl_drawfill),        323, 0},   // float(vector position, vector size, vector rgb, float alpha, optional float drawflag)
    {"drawsetcliparea",         PF_CSQC(PF_cl_drawsetclip),     324, 0},   // void(float x, float y, float width, float height)
    {"drawresetcliparea",       PF_CSQC(PF_cl_drawresetclip),   325, 0},   // void(void)
    {"drawstring",              PF_CSQC(PF_cl_drawstring),      326, 0},   // float(vector position, string text, vector size, vector rgb, float alpha, float drawflag)
    {"stringwidth",             PF_CSQC(PF_cl_stringwidth),     327, 0},   // float(string text, float usecolours, vector fontsize='8 8')
    {"drawsubpic",              PF_CSQC(PF_cl_drawsubpic),      328, 0},   // void(vector pos, vector sz, string pic, vector srcpos, vector srcsz, vector rgb, float alpha, optional float drawflag)

    {"getstati",                PF_CSQC(PF_cl_getstat_int),     330, 0},   // #define getstati_punf(stnum) (float)(__variant)getstati(stnum)\nint(float stnum)
    {"getstatf",                PF_CSQC(PF_cl_getstat_float),   331, 0},   // #define getstatbits getstatf\nfloat(float stnum, optional float firstbit, optional float bitcount)
    {"getstats",                PF_CSQC(PF_cl_getstat_string),  332, 0},   // string(float stnum)

    {"getplayerkeyvalue",       PF_CSQC(PF_cl_playerkey_s),     348, 0},   // string(float playernum, string keyname)
    {"getplayerkeyfloat",       PF_CSQC(PF_cl_playerkey_f), 0, 0},            // float(float playernum, string keyname, optional float assumevalue)

    {"registercommand",         PF_CSQC(PF_cl_registercommand), 352, 0},   // void(string cmdname)

    {"vectorvectors",           PF_BOTH(PF_vectorvectors),      432,    DP_QC_VECTORVECTORS},   // void(vector dir)

    {"clientcommand",           PF_SSQC(PF_clientcommand),      440,    KRIMZON_SV_PARSECLIENTCOMMAND}, // void(entity e, string s)
    {"tokenize",                PF_BOTH(PF_Tokenize),           441,    KRIMZON_SV_PARSECLIENTCOMMAND}, // float(string s)
    {"argv",                    PF_BOTH(PF_ArgV),               442,    KRIMZON_SV_PARSECLIENTCOMMAND}, // string(float n)
    {"argc",                    PF_BOTH(PF_ArgC), 0, 0},                      // float()

    {"asin",                    PF_BOTH(PF_asin),               471,    DP_QC_ASINACOSATANATAN2TAN},    // float(float s)
    {"acos",                    PF_BOTH(PF_acos),               472,    DP_QC_ASINACOSATANATAN2TAN},    // float(float c)
    {"atan",                    PF_BOTH(PF_atan),               473,    DP_QC_ASINACOSATANATAN2TAN},    // float(float t)
    {"atan2",                   PF_BOTH(PF_atan2),              474,    DP_QC_ASINACOSATANATAN2TAN},    // float(float c, float s)
    {"tan",                     PF_BOTH(PF_tan),                475,    DP_QC_ASINACOSATANATAN2TAN},    // float(float a)

    {"tokenize_console",        PF_BOTH(PF_tokenize_console),   514,    DP_QC_TOKENIZE_CONSOLE},        // float(string str)

    {"sprintf",                 PF_BOTH(PF_sprintf),            627,    DP_QC_SPRINTF},                 // string(string fmt, ...)
};

int pr_numbuiltindefs = Q_COUNTOF(pr_builtindefs);
