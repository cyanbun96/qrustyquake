// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
// sv_edict.c -- entity dictionary
#include "quakedef.h"

static s8 *pr_strings;
static s32 pr_stringssize;
static const s8 **pr_knownstrings;
static s32 pr_maxknownstrings;
static s32 pr_numknownstrings;
static ddef_t *pr_fielddefs;
static ddef_t *pr_globaldefs;
const s32 type_size[NUM_TYPE_SIZES] = {
	1, // ev_void
	1, // sizeof(string_t) / 4 // ev_string
	1, // ev_float
	3, // ev_vector
	1, // ev_entity
	1, // ev_field
	1, // sizeof(func_t) / 4 // ev_function
	1 // sizeof(void *) / 4 // ev_pointer
};
static ddef_t *ED_FieldAtOfs(s32 ofs);
static bool ED_ParseEpair(void *base, ddef_t *key, const s8 *s);
static gefv_cache gefvCache[GEFV_CACHESIZE] = { { NULL, "" }, { NULL, "" } };

void ED_ClearEdict(edict_t *e)
{
	memset(&e->v, 0, progs->entityfields * 4);
	e->free = 0;
}

// Either finds a free edict, or allocates a new one.
// Try to avoid reusing an entity that was recently freed, because it
// can cause the client to think the entity morphed into something else
// instead of being removed and recreated, which can cause interpolated
// angles and bad trails.
edict_t *ED_Alloc()
{
	s32 i = svs.maxclients + 1;
	for(; i < sv.num_edicts; i++) {
		edict_t *e = EDICT_NUM(i);
		// the first couple seconds of server time can involve a lot of
		// freeing and allocating, so relax the replacement policy
		if(e->free &&( e->freetime < 2 || sv.time - e->freetime > 0.5)){
			ED_ClearEdict(e);
			return e;
		}
	}
	if(i==sv.max_edicts)//johnfitz - use sv.max_edicts instead of MAX_EDICTS
	Host_Error("ED_Alloc: no free edicts(max_edicts is %i)", sv.max_edicts);
	sv.num_edicts++;
	edict_t *e = EDICT_NUM(i);
	memset(e, 0, pr_edict_size); // ericw -- switched sv.edicts to malloc()
	e->baseline.scale = ENTSCALE_DEFAULT;
	return e;
}

void ED_Free(edict_t *ed)
{
	SV_UnlinkEdict(ed); // unlink from world bsp
	ed->free = 1;
	ed->v.model = 0;
	ed->v.takedamage = 0;
	ed->v.modelindex = 0;
	ed->v.colormap = 0;
	ed->v.skin = 0;
	ed->v.frame = 0;
	VectorCopy(vec3_origin, ed->v.origin);
	VectorCopy(vec3_origin, ed->v.angles);
	ed->v.nextthink = -1;
	ed->v.solid = 0;
	ed->alpha = ENTALPHA_DEFAULT; //johnfitz -- reset alpha for next entity
	ed->scale = ENTSCALE_DEFAULT;
	ed->freetime = sv.time;
}

static ddef_t *ED_GlobalAtOfs(s32 ofs)
{
	for(s32 i = 0; i < progs->numglobaldefs; i++) {
		ddef_t *def = &pr_globaldefs[i];
		if(def->ofs == ofs) return def;
	}
	return NULL;
}

static ddef_t *ED_FieldAtOfs(s32 ofs)
{
	for(s32 i = 1; i < progs->numfielddefs; i++) {
		ddef_t *def = &pr_fielddefs[i];
		if(def->ofs == ofs) return def;
	}
	return NULL;
}

static ddef_t *ED_FindField(const s8 *name)
{
	for(s32 i = 0; i < progs->numfielddefs; i++) {
		ddef_t *def = &pr_fielddefs[i];
		if(!strcmp(PR_GetString(def->s_name), name))
			return def;
	}
	return NULL;
}

static ddef_t *ED_FindGlobal(const s8 *name)
{
	for(s32 i = 0; i < progs->numglobaldefs; i++) {
		ddef_t *def = &pr_globaldefs[i];
		if(!strcmp(PR_GetString(def->s_name), name))
			return def;
	}
	return NULL;
}

static dfunction_t *ED_FindFunction(const s8 *fn_name)
{
	for(s32 i = 0; i < progs->numfunctions; i++) {
		dfunction_t *func = &pr_functions[i];
		if(!strcmp(PR_GetString(func->s_name), fn_name))
			return func;
	}
	return NULL;
}

eval_t *GetEdictFieldValue(edict_t *ed, const s8 *field)
{
	ddef_t *def = NULL;
	static s32 rep = 0;
	for(s32 i = 0; i < GEFV_CACHESIZE; i++) {
		if(!strcmp(field, gefvCache[i].field)) {
			def = gefvCache[i].pcache;
			goto Done;
		}
	}
	def = ED_FindField(field);
	if(strlen(field) < MAX_FIELD_LEN) {
		gefvCache[rep].pcache = def;
		strcpy(gefvCache[rep].field, field);
		rep ^= 1;
	}
Done:
	if(!def) return NULL;
	return(eval_t *)((s8 *)&ed->v + def->ofs*4);
}

static const s8 *PR_ValueString(s32 type, eval_t *val)
{ // Returns a string describing *data in a type specific manner
	static s8 line[512];
	ddef_t *def;
	dfunction_t *f;
	type &= ~DEF_SAVEGLOBAL;
	switch(type) {
#define Q_ q_snprintf(line, sizeof(line),
	case ev_string: Q_"%s", PR_GetString(val->string)); break;
	case ev_void: Q_"void"); break;
	case ev_float: Q_"%5.1f", val->_float); break;
	case ev_pointer: Q_"pointer"); break;
	case ev_entity: Q_"entity %i",NUM_FOR_EDICT(PROG_TO_EDICT(val->edict)));
		break;
	case ev_function:
		f = pr_functions + val->function;
		Q_"%s()", PR_GetString(f->s_name));
		break;
	case ev_field:
		def = ED_FieldAtOfs(val->_int);
		Q_".%s", PR_GetString(def->s_name));
		break;
	case ev_vector:
	  Q_"'%5.1f %5.1f %5.1f'",val->vector[0],val->vector[1],val->vector[2]);
		break;
	default: Q_"bad type %i", type); break;
#undef Q_
	}
	return line;
}


static const s8 *PR_UglyValueString(s32 type, eval_t *val)
{ // Easier to parse than PR_ValueString
	static s8 line[1024];
	ddef_t *def;
	dfunction_t *f;
	type &= ~DEF_SAVEGLOBAL;
	switch(type)
	{
#define Q_ q_snprintf(line, sizeof(line),
	case ev_string: Q_"%s", PR_GetString(val->string)); break;
	case ev_entity: Q_"%i", NUM_FOR_EDICT(PROG_TO_EDICT(val->edict)));break;
	case ev_void: Q_"void"); break;
	case ev_float: Q_"%f", val->_float); break;
	case ev_function:
		f = pr_functions + val->function;
		Q_"%s", PR_GetString(f->s_name));
		break;
	case ev_field:
		def = ED_FieldAtOfs( val->_int );
		Q_"%s", PR_GetString(def->s_name));
		break;
	case ev_vector:
		Q_"%f %f %f", val->vector[0], val->vector[1], val->vector[2]);
		break;
	default: Q_"bad type %i", type); break;
#undef Q_
	}
	return line;
}

const s8 *PR_GlobalString(s32 ofs) // Returns a string with a description and
{ // the contents of a global, padded to 20 field width
	static s8 line[512];
	static const s32 lastchari = Q_COUNTOF(line) - 2;
	void *val = (void *)&pr_globals[ofs];
	ddef_t *def = ED_GlobalAtOfs(ofs);
	if(!def) q_snprintf(line, sizeof(line), "%i(?)", ofs);
	else {
		const s8 *s = PR_ValueString(def->type, (eval_t *)val);
		q_snprintf(line, sizeof(line), "%i(%s)%s", ofs,
				PR_GetString(def->s_name), s);
	}
	s32 i = strlen(line);
	for( ; i < 20; i++) strcat(line, " ");
	if(i < lastchari) strcat(line, " ");
	else line[lastchari] = ' ';
	return line;
}
const s8 *PR_GlobalStringNoContents(s32 ofs)
{
	static s8 line[512];
	static const s32 lastchari = Q_COUNTOF(line) - 2;
	ddef_t *def = ED_GlobalAtOfs(ofs);
	if(!def) q_snprintf(line, sizeof(line), "%i(?)", ofs);
	else q_snprintf(line, sizeof(line), "%i(%s)", ofs,
			PR_GetString(def->s_name));
	s32 i = strlen(line);
	for( ; i < 20; i++) strcat(line, " ");
	if(i < lastchari) strcat(line, " ");
	else line[lastchari] = ' ';
	return line;
}

void ED_Print(edict_t *ed)
{ // For debugging
	if(ed->free) { Con_Printf("FREE\n"); return; }
	Con_SafePrintf("\nEDICT %i:\n", NUM_FOR_EDICT(ed));
	for(s32 i = 1; i < progs->numfielddefs; i++) {
		ddef_t *d = &pr_fielddefs[i];
		const s8 *name = PR_GetString(d->s_name);
		s32 l = strlen(name);
		if(l > 1 && name[l - 2] == '_')
			continue; // skip _x, _y, _z vars
		s32 *v = (s32 *)((s8 *)&ed->v + d->ofs*4);
		// if the value is still all 0, skip the field
		s32 type = d->type & ~DEF_SAVEGLOBAL;
		if(type >= NUM_TYPE_SIZES) continue;
		s32 j = 0;
		for(; j < type_size[type]; j++)
			if(v[j]) break;
		if(j == type_size[type]) continue;
		Con_SafePrintf("%s", name);
		while(l++ < 15) Con_SafePrintf(" ");
		Con_SafePrintf("%s\n", PR_ValueString(d->type, (eval_t *)v));
	}
}


void ED_Write(FILE *f, edict_t *ed)
{ // For savegames
	fprintf(f, "{\n");
	if(ed->free) { fprintf(f, "}\n"); return; }
	for(s32 i = 1; i < progs->numfielddefs; i++) {
		ddef_t *d = &pr_fielddefs[i];
		const s8 *name = PR_GetString(d->s_name);
		s32 j = strlen(name);
		if(j > 1 && name[j - 2] == '_')
			continue; // skip _x, _y, _z vars
		s32 *v = (s32 *)((s8 *)&ed->v + d->ofs*4);
		// if the value is still all 0, skip the field
		s32 type = d->type & ~DEF_SAVEGLOBAL;
		if(type >= NUM_TYPE_SIZES) continue;
		for(j = 0; j < type_size[type]; j++)
			if(v[j]) break;
		if(j == type_size[type]) continue;
		fprintf(f, "\"%s\" ", name);
		fprintf(f, "\"%s\"\n", PR_UglyValueString(d->type,(eval_t *)v));
	}
//johnfitz -- save entity alpha manually when progs.dat doesn't know about alpha
	if(!pr_alpha_supported && ed->alpha != ENTALPHA_DEFAULT)
		fprintf(f, "\"alpha\" \"%f\"\n", ENTALPHA_TOSAVE(ed->alpha));
	fprintf(f, "}\n");
}

void ED_PrintNum(s32 ent) { ED_Print(EDICT_NUM(ent)); }

void ED_PrintEdicts()
{ // For debugging, prints all the entities in the current server
	if(!sv.active) return;
	Con_Printf("%i entities\n", sv.num_edicts);
	for(s32 i = 0; i < sv.num_edicts; i++)
		ED_PrintNum(i);
}

static void ED_PrintEdict_f()
{ // For debugging, prints a single edict
	if(!sv.active) return;
	s32 i = Q_atoi(Cmd_Argv(1));
	if(i < 0 || i >= sv.num_edicts) {
		Con_Printf("Bad edict number\n");
		return;
	}
	ED_PrintNum(i);
}

static void ED_Count()
{ // For debugging
	if(!sv.active) return;
	s32 active = 0, models = 0, solid = 0, step = 0;
	for(s32 i = 0; i < sv.num_edicts; i++) {
		edict_t *ent = EDICT_NUM(i);
		if(ent->free) continue;
		active++;
		if(ent->v.solid) solid++;
		if(ent->v.model) models++;
		if(ent->v.movetype == MOVETYPE_STEP) step++;
	}
	Con_Printf("num_edicts:%3i\n", sv.num_edicts);
	Con_Printf("active :%3i\n", active);
	Con_Printf("view :%3i\n", models);
	Con_Printf("touch :%3i\n", solid);
	Con_Printf("step :%3i\n", step);
}

void ED_WriteGlobals(FILE *f)
{
	fprintf(f, "{\n");
	for(s32 i = 0; i < progs->numglobaldefs; i++) {
		ddef_t *def = &pr_globaldefs[i];
		s32 type = def->type;
		if(!(def->type & DEF_SAVEGLOBAL)) continue;
		type &= ~DEF_SAVEGLOBAL;
		if(type != ev_string && type != ev_float && type != ev_entity)
			continue;
		const s8 *name = PR_GetString(def->s_name);
		fprintf(f, "\"%s\" ", name);
		fprintf(f, "\"%s\"\n", PR_UglyValueString(type,
					(eval_t *)&pr_globals[def->ofs]));
	}
	fprintf(f, "}\n");
}

const s8 *ED_ParseGlobals(const s8 *data)
{
	s8 keyname[64];
	while(1) {
		data = COM_Parse(data); // parse key
		if(com_token[0] == '}') break;
		if(!data)
			Host_Error("ED_ParseEntity: EOF without closing brace");
		q_strlcpy(keyname, com_token, sizeof(keyname));
		data = COM_Parse(data); // parse value
		if(!data)
			Host_Error("ED_ParseEntity: EOF without closing brace");
		if(com_token[0] == '}')
		       Host_Error("ED_ParseEntity: closing brace without data");
		ddef_t *key = ED_FindGlobal(keyname);
		if(!key) {
			Con_Printf("'%s' is not a global\n", keyname);
			continue;
		}
		if(!ED_ParseEpair((void *)pr_globals, key, com_token))
			Host_Error("ED_ParseGlobals: parse error");
	}
	return data;
}

static string_t ED_NewString(const s8 *string)
{
	s32 l = strlen(string) + 1;
	s8 *new_p;
	string_t num = PR_AllocString(l, &new_p);
	for(s32 i = 0; i < l; i++) {
		if(string[i] == '\\' && i < l-1) {
			i++;
			if(string[i] == 'n') *new_p++ = '\n';
			else *new_p++ = '\\';
		}
		else *new_p++ = string[i];
	}
	return num;
}


static bool ED_ParseEpair(void *base, ddef_t *key, const s8 *s)
{ // Can parse either fields or globals, returns 0 if error
	s8 string[128];
	void *d = (void *)((s32 *)base + key->ofs);
	dfunction_t *func; // keep here for OpenBSD
	ddef_t *def; // ditto
	switch(key->type & ~DEF_SAVEGLOBAL) {
	case ev_string: *(string_t *)d = ED_NewString(s); break;
	case ev_float: *(f32 *)d = atof(s); break;
	case ev_entity: *(s32 *)d = EDICT_TO_PROG(EDICT_NUM(atoi(s))); break;
	case ev_vector:
		q_strlcpy(string, s, sizeof(string));
		s8 *end = (s8 *)string + strlen(string);
		s8 *v = string;
		s8 *w = string;
		s32 i = 0;
		for(; i < 3 && (w <= end); i++) {
			// set v to the next space(or 0 u8)
			while(*v && *v != ' ') v++;
			*v = 0; // and change that s8 to a 0 u8
			((f32 *)d)[i] = atof(w);
			w = v = v+1;
		}
	// fill with 0 in case we hit the end of string before reading 3 floats
		if(i < 3) {
			printf("Avoided reading garbage for \"%s\" \"%s\"\n",
				PR_GetString(key->s_name), s);
			for(; i < 3; i++)
				((f32 *)d)[i] = 0.0f;
		}
		break;
	case ev_field:
		def = ED_FindField(s);
		if(!def) {
// suppress error becuase fog/sky fields might not be mentioned in defs.qc
			if(strncmp(s, "sky", 3) && strcmp(s, "fog"))
				Con_DPrintf("Can't find field %s\n", s);
			return 0;
		}
		*(s32 *)d = G_INT(def->ofs);
		break;
	case ev_function:
		func = ED_FindFunction(s);
		if(!func) {
			Con_Printf("Can't find function %s\n", s);
			return 0;
		}
		*(func_t *)d = func - pr_functions;
		break;
	default:
		break;
	}
	return 1;
}

// Parses an edict out of the given string, returning the new position
// ed should be a properly initialized empty edict.
// Used for initial level load and for savegames.
const s8 *ED_ParseEdict(const s8 *data, edict_t *ent)
{
	s8 keyname[256];
	bool init = 0;
	if(ent != sv.edicts) // clear it
		memset(&ent->v, 0, progs->entityfields * 4);
	// go through all the dictionary pairs
	while(1) { // parse key
		data = COM_Parse(data);
		if(com_token[0] == '}') break;
		if(!data)
			Host_Error("ED_ParseEntity: EOF without closing brace");
		// anglehack is to allow QuakeEd to write single scalar angles
		// and allow them to be turned into vectors.
		bool anglehack = 0;
		if(!strcmp(com_token, "angle")) {
			strcpy(com_token, "angles");
			anglehack = 1;
		}
		if(!strcmp(com_token, "light")) // hack for single light def
			strcpy(com_token, "light_lev");
		q_strlcpy(keyname, com_token, sizeof(keyname));
		// another hack to fix keynames with trailing spaces
		s32 n = strlen(keyname);
		while(n && keyname[n-1] == ' ') {
			keyname[n-1] = 0;
			n--;
		}
		// parse value
		// HACK: we allow truncation when reading the wad field,
		// otherwise maps using lots of wads with absolute paths
		// could cause a parse error
		data = COM_ParseEx(data, !strcmp(keyname, "wad") ?
				CPE_ALLOWTRUNC : CPE_NOTRUNC);
		if(!data)
			Host_Error("ED_ParseEntity: EOF without closing brace");
		if(com_token[0] == '}')
		       Host_Error("ED_ParseEntity: closing brace without data");
		init = 1;
		// keynames with a leading underscore are used for utility
		// comments, and are immediately discarded by quake
		if(keyname[0] == '_') continue;
		// support .alpha even when progs.dat doesn't know about it
		if(!strcmp(keyname, "alpha"))
			ent->alpha = ENTALPHA_ENCODE(Q_atof(com_token));
		ddef_t *key = ED_FindField(keyname);
		if(!key) {
// suppress error becuase fog/sky/alpha fields might not be mentioned in defs.qc
if(strncmp(keyname,"sky",3) && strcmp(keyname,"fog") && strcmp(keyname,"alpha"))
	Con_DPrintf("\"%s\" is not a field\n", keyname);
			continue;
		}
		if(anglehack) {
			s8 temp[32];
			strcpy(temp, com_token);
			sprintf(com_token, "0 %s 0", temp);
		}
		if(!ED_ParseEpair((void *)&ent->v, key, com_token))
			Host_Error("ED_ParseEdict: parse error");
	}
	if(!init) ent->free = 1;
	return data;
}

// The entities are directly placed in the array, rather than allocated with
// ED_Alloc, because otherwise an error loading the map would have entity
// number references out of order.
// Creates a server's entity / program execution context by
// parsing textual entity definitions out of an ent file.
// Used for both fresh maps and savegame loads. A fresh map would also need
// to call ED_CallSpawnFunctions() to let the objects initialize themselves.
void ED_LoadFromFile(const s8 *data)
{
	const s8 *classname;
	s8 spawnfunc[256];
	dfunction_t *func;
	edict_t *ent = NULL;
	s32 inhibit = 0;
	pr_global_struct->time = sv.time;
	// parse ents
	while(1) {
		// parse the opening brace
		data = COM_Parse(data);
		if(!data) break;
		if(com_token[0] != '{')
	Host_Error("ED_LoadFromFile: found %s when expecting {",com_token);
		if(!ent) ent = EDICT_NUM(0);
		else ent = ED_Alloc();
		data = ED_ParseEdict(data, ent);
		// remove things from different skill levels or deathmatch
		if(deathmatch.value) {
			if(((s32)ent->v.spawnflags & SPAWNFLAG_NOT_DEATHMATCH)){
				ED_Free(ent);
				inhibit++;
				continue;
			}
		}
else if((current_skill == 0 && ((s32)ent->v.spawnflags & SPAWNFLAG_NOT_EASY))
	|| (current_skill == 1 && ((s32)ent->v.spawnflags&SPAWNFLAG_NOT_MEDIUM))
	|| (current_skill >= 2 && ((s32)ent->v.spawnflags&SPAWNFLAG_NOT_HARD))){
			ED_Free(ent);
			inhibit++;
			continue;
		}
		// immediately call spawn function
		if(!ent->v.classname) {
			Con_SafePrintf("No classname for:\n");
			ED_Print(ent);
			ED_Free(ent);
			continue;
		}
		classname = PR_GetString(ent->v.classname);
		if(sv.nomonsters && !Q_strncmp(classname, "monster_", 8))
		{
			ED_Free(ent);
			inhibit++;
			continue;
		}
	// look for the spawn function
//erysdren: look for FTE/DP spawnfunc_* function first to support QuakeC classes
	    q_snprintf(spawnfunc, sizeof(spawnfunc), "spawnfunc_%s", classname);
		func = ED_FindFunction(spawnfunc);
		if(!func) {
			func = ED_FindFunction(classname);
			if(!func) {
				Con_SafePrintf("No spawn function for:\n");
				ED_Print(ent);
				ED_Free(ent);
				continue;
			}
		}
		SV_ReserveSignonSpace(512);
		pr_global_struct->self = EDICT_TO_PROG(ent);
		PR_ExecuteProgram(func - pr_functions);
	}
	Con_DPrintf("%i entities inhibited\n", inhibit);
}

static bool PR_HasGlobal(const s8 *name, f32 value)
{
	ddef_t *g = ED_FindGlobal(name);
	return g &&(g->type&~DEF_SAVEGLOBAL)==ev_float &&G_FLOAT(g->ofs)==value;
}

// Checks for the presence of Quake 2021 effects flags and returns a mask with
// the correspondings bits either on or off depending on the result, in order to
// avoid conflicts(e.g. Arcane Dimensions uses bit 32 for its explosions)
static s32 PR_FindSupportedEffects()
{
	bool isqex = PR_HasGlobal("EF_QUADLIGHT", EF_QEX_QUADLIGHT) &&
		(PR_HasGlobal("EF_PENTLIGHT", EF_QEX_PENTALIGHT) ||
		 PR_HasGlobal("EF_PENTALIGHT", EF_QEX_PENTALIGHT)) ;
	return isqex ? -1 : -1 & ~(EF_QEX_QUADLIGHT|EF_QEX_PENTALIGHT|
			EF_QEX_CANDLELIGHT);
}

static const exbuiltin_t exbuiltins[] = { // for 2021 re-release
//Upd-1 adds the following builtins with new ids. Patch them to use old indices.
//https://steamcommunity.com/games/2310/announcements/detail/2943653788150871156
	{ "centerprint", -90, -73 },
	{ "bprint", -91, -23 },
	{ "sprint", -92, -24 },
//Upd-3 changes its unique builtins to be looked up by name instead of builtin
//numbers, to avoid conflict with other engines. Patch them to use our indices.
//https://steamcommunity.com/games/2310/announcements/detail/3177861894960065435
	{ "ex_centerprint", 0, -73 },
	{ "ex_bprint", 0, -23 },
	{ "ex_sprint", 0, -24 },
	{ "ex_finaleFinished", 0, -79 },
	{ "ex_localsound", 0, -80 },
	{ "ex_draw_point", 0, -81 },
	{ "ex_draw_line", 0, -82 },
	{ "ex_draw_arrow", 0, -83 },
	{ "ex_draw_ray", 0, -84 },
	{ "ex_draw_circle", 0, -85 },
	{ "ex_draw_bounds", 0, -86 },
	{ "ex_draw_worldtext", 0, -87 },
	{ "ex_draw_sphere", 0, -88 },
	{ "ex_draw_cylinder", 0, -89 },
	{ "ex_CheckPlayerEXFlags", 0, -90 },
	{ "ex_walkpathtogoal", 0, -91 },
	{ "ex_bot_movetopoint", 0, -92 },
	{ "ex_bot_followentity", 0, -92 },
	{ NULL, 0, 0 } // end-of-list
};

static void PR_PatchRereleaseBuiltins()
{
	const exbuiltin_t *ex = exbuiltins;
	for( ; ex->name != NULL; ++ex) {
		dfunction_t *f = ED_FindFunction(ex->name);
		if(f && f->first_statement == ex->first_statement)
			f->first_statement = ex->patch_statement;
	}
}

void PR_LoadProgs()
{
	s32 i; // flush the non-C variable lookup cache
	for(i = 0; i < GEFV_CACHESIZE; i++)
		gefvCache[i].field[0] = 0;
	CRC_Init(&pr_crc);
	progs = (dprograms_t *)COM_LoadHunkFile("progs.dat", NULL);
	if(!progs)
		Host_Error("PR_LoadProgs: couldn't load progs.dat");
	Con_DPrintf("Programs occupy %iK.\n", com_filesize/1024);
	for(i = 0; i < com_filesize; i++)
		CRC_ProcessByte(&pr_crc, ((u8 *)progs)[i]);
	for(i = 0; i < (s32) sizeof(*progs) / 4; i++) // u8 swap the header
		((s32 *)progs)[i] = LittleLong( ((s32 *)progs)[i] );
	if(progs->version != PROG_VERSION)
	       Host_Error("progs.dat has wrong version number(%i should be %i)",
			progs->version, PROG_VERSION);
	if(progs->crc != PROGHEADER_CRC)
		Host_Error(
	 "progs.dat system vars have been modified, progdefs.h is out of date");
	pr_functions = (dfunction_t *)((u8 *)progs + progs->ofs_functions);
	pr_strings = (s8 *)progs + progs->ofs_strings;
	if(progs->ofs_strings + progs->numstrings >= com_filesize)
		Host_Error("progs.dat strings go past end of file\n");
	pr_numknownstrings = 0; // initialize the strings
	pr_maxknownstrings = 0;
	pr_stringssize = progs->numstrings;
	if(pr_knownstrings)
		Z_Free((void *)pr_knownstrings);
	pr_knownstrings = NULL;
	PR_SetEngineString("");
	pr_globaldefs = (ddef_t *)((u8 *)progs + progs->ofs_globaldefs);
	pr_fielddefs = (ddef_t *)((u8 *)progs + progs->ofs_fielddefs);
	pr_statements = (dstatement_t *)((u8 *)progs + progs->ofs_statements);
	pr_global_struct = (globalvars_t *)((u8 *)progs + progs->ofs_globals);
	pr_globals = (f32 *)pr_global_struct;
	// u8 swap the lumps
	for(i = 0; i < progs->numstatements; i++) {
		pr_statements[i].op = LittleShort(pr_statements[i].op);
		pr_statements[i].a = LittleShort(pr_statements[i].a);
		pr_statements[i].b = LittleShort(pr_statements[i].b);
		pr_statements[i].c = LittleShort(pr_statements[i].c);
	}
	for(i = 0; i < progs->numfunctions; i++) {
		pr_functions[i].first_statement =
			LittleLong(pr_functions[i].first_statement);
		pr_functions[i].parm_start =
			LittleLong(pr_functions[i].parm_start);
		pr_functions[i].s_name = LittleLong(pr_functions[i].s_name);
		pr_functions[i].s_file = LittleLong(pr_functions[i].s_file);
		pr_functions[i].numparms = LittleLong(pr_functions[i].numparms);
		pr_functions[i].locals = LittleLong(pr_functions[i].locals);
	}
	for(i = 0; i < progs->numglobaldefs; i++) {
		pr_globaldefs[i].type = LittleShort(pr_globaldefs[i].type);
		pr_globaldefs[i].ofs = LittleShort(pr_globaldefs[i].ofs);
		pr_globaldefs[i].s_name = LittleLong(pr_globaldefs[i].s_name);
	}
	pr_alpha_supported = 0; //johnfitz
	for(i = 0; i < progs->numfielddefs; i++) {
		pr_fielddefs[i].type = LittleShort(pr_fielddefs[i].type);
		if(pr_fielddefs[i].type & DEF_SAVEGLOBAL)
	Host_Error("PR_LoadProgs: pr_fielddefs[i].type & DEF_SAVEGLOBAL");
		pr_fielddefs[i].ofs = LittleShort(pr_fielddefs[i].ofs);
		pr_fielddefs[i].s_name = LittleLong(pr_fielddefs[i].s_name);
		if(!strcmp(pr_strings + pr_fielddefs[i].s_name,"alpha"))
			pr_alpha_supported=1;//detect alpha support in progs.dat
	}
	for(i = 0; i < progs->numglobals; i++)
		((s32 *)pr_globals)[i] = LittleLong(((s32 *)pr_globals)[i]);
	pr_edict_size=progs->entityfields*4+sizeof(edict_t)-sizeof(entvars_t);
	// round off to next highest whole word address(esp for Alpha)
	// this ensures that pointers in the engine data area are always
	// properly aligned
	pr_edict_size += sizeof(void *) - 1;
	pr_edict_size &= ~(sizeof(void *) - 1);
	PR_PatchRereleaseBuiltins();
	pr_effects_mask = PR_FindSupportedEffects();
}

static void ED_Nomonsters_f(cvar_t *cvar)
{ if(cvar->value) printf("\"%s\" can break gameplay.\n", cvar->name); }

void PR_Init()
{
	Cmd_AddCommand("edict", ED_PrintEdict_f);
	Cmd_AddCommand("edicts", ED_PrintEdicts);
	Cmd_AddCommand("edictcount", ED_Count);
	Cmd_AddCommand("profile", PR_Profile_f);
	Cvar_RegisterVariable(&nomonsters);
	Cvar_SetCallback(&nomonsters, ED_Nomonsters_f);
	Cvar_RegisterVariable(&gamecfg);
	Cvar_RegisterVariable(&scratch1);
	Cvar_RegisterVariable(&scratch2);
	Cvar_RegisterVariable(&scratch3);
	Cvar_RegisterVariable(&scratch4);
	Cvar_RegisterVariable(&savedgamecfg);
	Cvar_RegisterVariable(&saved1);
	Cvar_RegisterVariable(&saved2);
	Cvar_RegisterVariable(&saved3);
	Cvar_RegisterVariable(&saved4);
}

edict_t *EDICT_NUM(s32 n)
{
	if(n < 0 || n >= sv.max_edicts)
		Host_Error("EDICT_NUM: bad number %i", n);
	return(edict_t *)((u8 *)sv.edicts + (n)*pr_edict_size);
}

s32 NUM_FOR_EDICT(edict_t *e)
{
	s32 b = (u8 *)e - (u8 *)sv.edicts;
	b = b / pr_edict_size;
	if(b < 0 || b >= sv.num_edicts)
		Host_Error("NUM_FOR_EDICT: bad pointer");
	return b;
}

static void PR_AllocStringSlots()
{
	pr_maxknownstrings += PR_STRING_ALLOCSLOTS;
	pr_knownstrings = (const s8 **) Z_Realloc((void *)pr_knownstrings,
					pr_maxknownstrings * sizeof(s8 *));
}

const s8 *PR_GetString(s32 num)
{
	if(num >= 0 && num < pr_stringssize)
		return pr_strings + num;
	else if(num < 0 && num >= -pr_numknownstrings) {
		if(!pr_knownstrings[-1 - num]) {
     Host_Error("PR_GetString: attempt to get a non-existant string %d\n", num);
			return "";
		}
		return pr_knownstrings[-1 - num];
	} else {
		Host_Error("PR_GetString: invalid string offset %d\n", num);
		return "";
	}
}
s32 PR_SetEngineString(const s8 *s)
{
	if(!s) return 0;
	if(s >= pr_strings && s <= pr_strings + pr_stringssize - 2)
		return(s32)(s - pr_strings);
	s32 i = 0;
	for(; i < pr_numknownstrings; i++) {
		if(pr_knownstrings[i] == s)
			return -1 - i;
	}
	if(i >= pr_maxknownstrings)
		PR_AllocStringSlots();
	pr_numknownstrings++;
	pr_knownstrings[i] = s;
	return -1 - i;
}

s32 PR_AllocString(s32 size, s8 **ptr)
{
	if(!size) return 0;
	s32 i = 0;
	for(; i < pr_numknownstrings; i++)
		if(!pr_knownstrings[i]) break;
	if(i >= pr_maxknownstrings) PR_AllocStringSlots();
	pr_numknownstrings++;
	pr_knownstrings[i] = (s8 *)Hunk_AllocName(size, "string");
	if(ptr) *ptr = (s8 *) pr_knownstrings[i];
	return -1 - i;
}
