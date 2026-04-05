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
static bool ED_ParseEpair (void *base, ddef_t *key, const char *s, bool zoned);
static gefv_cache gefvCache[GEFV_CACHESIZE] = { { NULL, "" }, { NULL, "" } };

static void PR_HashInit (prhashtable_t *table, int capacity, const char *name)
{
    capacity *= 2; // 50% load factor
    table->capacity = capacity;
    table->strings = (const char **) Hunk_AllocName (sizeof(*table->strings) * capacity, name);
    table->indices = (int         *) Hunk_AllocName (sizeof(*table->indices) * capacity, name);
}

void PR_SwitchQCVM(qcvm_t *nvm)
{
    if (qcvm && nvm)
        Sys_Error("PR_SwitchQCVM: A qcvm was already active");
    qcvm = nvm;
    if (qcvm)
        pr_global_struct = (globalvars_t*)qcvm->globals;
    else
        pr_global_struct = NULL;
}

void PR_PushQCVM(qcvm_t *newvm, qcvm_t **oldvm)
{
    *oldvm = qcvm;
    PR_SwitchQCVM(NULL);
    PR_SwitchQCVM(newvm);
}

void PR_PopQCVM(qcvm_t *oldvm)
{
    PR_SwitchQCVM(NULL);
    PR_SwitchQCVM(oldvm);
}

static int PR_HashGet (prhashtable_t *table, const char *key)
{
    unsigned pos = COM_HashString (key) % table->capacity, end = pos;

    do
    {
        const char *s = table->strings[pos];
        if (!s)
            return -1;
        if (0 == strcmp(s, key))
            return table->indices[pos];

        ++pos;
        if (pos == table->capacity)
            pos = 0;
    }
    while (pos != end);

    return -1;
}

void PR_ClearEngineString (int num)
{
    if (num < 0 && num >= -qcvm->numknownstrings)
    {
        num = -1 - num;
        qcvm->knownstrings[num] = (const char*) qcvm->firstfreeknownstring;
        qcvm->firstfreeknownstring = &qcvm->knownstrings[num];
    }
}

// (etype_t type, eval_t *val)
// Returns a string describing *data in a type specific manner
// Easier to parse than PR_ValueString
static const char *PR_UglyValueString (int type, eval_t *val)
{
    static char line[1024];
    ddef_t      *def;
    dfunction_t *f;

    type &= ~DEF_SAVEGLOBAL;

    switch (type)
    {
    case ev_string:
        q_snprintf (line, sizeof(line), "%s", PR_GetString(val->string));
        break;
    case ev_entity:
        q_snprintf (line, sizeof(line), "%i", NUM_FOR_EDICT(PROG_TO_EDICT(val->edict)));
        break;
    case ev_function:
        f = qcvm->functions + val->function;
        q_snprintf (line, sizeof(line), "%s", PR_GetString(f->s_name));
        break;
    case ev_field:
        def = ED_FieldAtOfs ( val->_int );
        q_snprintf (line, sizeof(line), "%s", PR_GetString(def->s_name));
        break;
    case ev_void:
        q_snprintf (line, sizeof(line), "void");
        break;
    case ev_float:
        q_snprintf (line, sizeof(line), "%f", val->_float);
        break;
    case ev_vector:
        q_snprintf (line, sizeof(line), "%f %f %f", val->vector[0], val->vector[1], val->vector[2]);
        break;
    default:
        q_snprintf (line, sizeof(line), "bad type %i", type);
        break;
    }

    return line;
}

static ddef_t *ED_FindGlobal (const char *name)
{
    ddef_t      *def;
    int         i;

    if (qcvm->ht_globals.capacity > 0)
    {
        int index = PR_HashGet (&qcvm->ht_globals, name);
        if (index < 0)
            return NULL;
        return qcvm->globaldefs + index;
    }

    for (i = 0; i < qcvm->progs->numglobaldefs; i++)
    {
        def = &qcvm->globaldefs[i];
        if ( !strcmp(PR_GetString(def->s_name), name) )
            return def;
    }
    return NULL;
}

static dfunction_t *ED_FindFunction (const char *fn_name)
{
    dfunction_t     *func;
    int             i;

    if (qcvm->ht_functions.capacity > 0)
    {
        int index = PR_HashGet (&qcvm->ht_functions, fn_name);
        if (index < 0)
            return NULL;
        return qcvm->functions + index;
    }

    for (i = 0; i < qcvm->progs->numfunctions; i++)
    {
        func = &qcvm->functions[i];
        if ( !strcmp(PR_GetString(func->s_name), fn_name) )
            return func;
    }
    return NULL;
}

static func_t PR_FindExtFunction(const char *entryname)
{   //depends on 0 being an invalid function,
    dfunction_t *func = ED_FindFunction(entryname);
    if (func)
        return func - qcvm->functions;
    return 0;
}

static void *PR_FindExtGlobal(int type, const char *name)
{
    ddef_t *def = ED_FindGlobal(name);
    if (def && (def->type&~DEF_SAVEGLOBAL) == type && def->ofs < qcvm->progs->numglobals)
        return qcvm->globals + def->ofs;
    return NULL;
}

void PR_EnableExtensions (void)
{
    unsigned int i;
    unsigned int numautocvars = 0;

    if (!pr_checkextension.value && qcvm == &sv.qcvm)
    {
        Con_DPrintf("not enabling qc extensions\n");
        return;
    }

#define QCEXTFUNC(n,t) qcvm->extfuncs.n = PR_FindExtFunction(#n);
#define QCEXTGLOBAL_FLOAT(n) qcvm->extglobals.n = PR_FindExtGlobal(ev_float, #n);
#define QCEXTGLOBAL_INT(n) qcvm->extglobals.n = PR_FindExtGlobal(ev_ext_integer, #n);
#define QCEXTGLOBAL_VECTOR(n) qcvm->extglobals.n = PR_FindExtGlobal(ev_vector, #n);

    if (qcvm == &cl.qcvm)
    {   //csqc
        QCEXTFUNCS_CS
        QCEXTGLOBALS_CSQC
    }
    else
    {   //ssqc
        QCEXTFUNCS_SV
    }

#undef QCEXTGLOBAL_FLOAT
#undef QCEXTGLOBAL_INT
#undef QCEXTGLOBAL_VECTOR
#undef QCEXTFUNC

    //autocvars
    for (i = 0; i < (unsigned int)qcvm->progs->numglobaldefs; i++)
    {
        const char *n = PR_GetString(qcvm->globaldefs[i].s_name);
        if (!strncmp(n, "autocvar_", 9))
        {
            //really crappy approach
            cvar_t *var = Cvar_Create(n + 9, PR_UglyValueString (qcvm->globaldefs[i].type, (eval_t*)(qcvm->globals + qcvm->globaldefs[i].ofs)));
            numautocvars++;
            if (!var)
                continue;   //name conflicts with a command?

            if (!ED_ParseEpair ((void *)qcvm->globals, &qcvm->globaldefs[i], var->string, true))
                Con_DPrintf("EXT: Unable to configure %s\n", n);
            var->flags |= CVAR_AUTOCVAR;
        }
    }
    if (numautocvars)
        Con_DPrintf("Found %i autocvars\n", numautocvars);
}

static void ED_AddToFreeList (edict_t *ed)
{
    ed->free = true;
    if ((u8 *)ed <= (u8 *)qcvm->edicts + q_max (svs.maxclients, 1) * qcvm->edict_size)
        return;
    if (ed->freechain.prev)
        RemoveLink (&ed->freechain);
    InsertLinkBefore (&ed->freechain, &qcvm->free_edicts);
}

static void ED_RemoveFromFreeList (edict_t *ed)
{
    ed->free = false;
    if (ed->freechain.prev)
    {
        RemoveLink (&ed->freechain);
        ed->freechain.prev = ed->freechain.next = NULL;
    }
}

void ED_ClearEdict (edict_t *e)
{ // Sets everything to NULL
    if (!e->free)
        SV_UnlinkEdict (e);
    else
        ED_RemoveFromFreeList (e);
    memset (&e->v, 0, qcvm->progs->entityfields * 4);
}

// Either finds a free edict, or allocates a new one.
// Try to avoid reusing an entity that was recently freed, because it
// can cause the client to think the entity morphed into something else
// instead of being removed and recreated, which can cause interpolated
// angles and bad trails.
edict_t *ED_Alloc (void)
{
    edict_t     *e;

    if (qcvm->free_edicts.next != &qcvm->free_edicts)
    {
        e = STRUCT_FROM_LINK (qcvm->free_edicts.next, edict_t, freechain);
        if (!e->free)
            Host_Error ("ED_Alloc: free list entity still in use");
        // the first couple seconds of server time can involve a lot of
        // freeing and allocating, so relax the replacement policy
        if (e->freetime < 2 || qcvm->time - e->freetime > 0.5)
        {
            ED_ClearEdict (e);
            return e;
        }
    }

    if (qcvm->num_edicts == qcvm->max_edicts) //johnfitz -- use sv.max_edicts instead of MAX_EDICTS
        Host_Error ("ED_Alloc: no free edicts (max_edicts is %i)", qcvm->max_edicts);

    e = EDICT_NUM(qcvm->num_edicts++);
    memset(e, 0, qcvm->edict_size); // ericw -- switched sv.edicts to malloc(), so we are accessing uninitialized memory and must fully zero it, not just ED_ClearEdict
    e->baseline.scale = ENTSCALE_DEFAULT;

    return e;
}

void ED_Free(edict_t *ed)
{
	SV_UnlinkEdict(ed); // unlink from world bsp
	ED_AddToFreeList(ed);
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
	ed->freetime = qcvm->time;
}

static ddef_t *ED_FindField (const char *name)
{
    ddef_t      *def;
    int         i;

    if (qcvm->ht_fields.capacity > 0)
    {
        int index = PR_HashGet (&qcvm->ht_fields, name);
        if (index < 0)
            return NULL;
        return qcvm->fielddefs + index;
    }

    for (i = 0; i < qcvm->progs->numfielddefs; i++)
    {
        def = &qcvm->fielddefs[i];
        if ( !strcmp(PR_GetString(def->s_name), name) )
            return def;
    }
    return NULL;
}

static ddef_t *ED_GlobalAtOfs (int ofs)
{
    if (ofs < 0 || ofs > qcvm->maxglobalofs)
        return NULL;

    ofs = qcvm->ofstoglobal[ofs];
    if (ofs < 0)
        return NULL;

    return &qcvm->globaldefs[ofs];
}

/*
============
ED_FieldAtOfs
============
*/
static ddef_t *ED_FieldAtOfs (int ofs)
{
    if (ofs < 0 || ofs > qcvm->maxfieldofs)
        return NULL;

    ofs = qcvm->ofstofield[ofs];
    if (ofs < 0)
        return NULL;

    return &qcvm->fielddefs[ofs];
}

/*
============
ED_FindFieldOffset
============
*/
int ED_FindFieldOffset (const char *name)
{
    ddef_t      *def = ED_FindField(name);
    if (!def)
        return -1;
    return def->ofs;
}

eval_t *GetEdictFieldValue(edict_t *ed, int fldofs)
{
    if (fldofs < 0)
        return NULL;

    return (eval_t *)((char *)&ed->v + fldofs*4);
}

eval_t *GetEdictFieldValueByName(edict_t *ed, const char *name)
{ return GetEdictFieldValue(ed, ED_FindFieldOffset(name)); }

static const char *PR_FloatFormat (float f)
{ return fabs (f - round (f)) < 0.05f ? "% 5.0f  " : "% 7.1f"; }

static const s8 *PR_ValueString(s32 type, eval_t *val)
{ // Returns a string describing *data in a type specific manner
    static char line[512];
    char        fmt[64];
    const char  *str;
    ddef_t      *def;
    dfunction_t *f;
    edict_t     *ed;

    type &= ~DEF_SAVEGLOBAL;

    switch (type)
    {
    case ev_string:
        q_snprintf (line, sizeof(line), "%s", PR_GetString(val->string));
        break;
    case ev_entity:
        ed = PROG_TO_EDICT(val->edict);
        str = PR_GetString(ed->v.classname);
        q_snprintf (line, sizeof(line), *str ? "entity %i (%s)" : "entity %i", NUM_FOR_EDICT(ed), PR_GetString(ed->v.classname));
        break;
    case ev_function:
        f = qcvm->functions + val->function;
        q_snprintf (line, sizeof(line), "%s()", PR_GetString(f->s_name));
        break;
    case ev_field:
        def = ED_FieldAtOfs ( val->_int );
        q_snprintf (line, sizeof(line), ".%s", PR_GetString(def->s_name));
        break;
    case ev_void:
        q_snprintf (line, sizeof(line), "void");
        break;
    case ev_float:
        // Note: leading space, so that float fields are aligned with the first value in vector fields
        q_snprintf (fmt, sizeof(fmt), " %s", PR_FloatFormat (val->_float));
        q_snprintf (line, sizeof(line), fmt, val->_float);
        break;
    case ev_vector:
        q_snprintf (fmt, sizeof(fmt), "'%s %s %s'", PR_FloatFormat (val->vector[0]), PR_FloatFormat (val->vector[1]), PR_FloatFormat (val->vector[2]));
        q_snprintf (line, sizeof(line), fmt, val->vector[0], val->vector[1], val->vector[2]);
        break;
    case ev_pointer:
        q_snprintf (line, sizeof(line), "pointer");
        break;
    default:
        q_snprintf (line, sizeof(line), "bad type %i", type);
        break;
    }

    return line;
}

const char *PR_GlobalString (int ofs) // Returns a string with a description and
{ // the contents of a global, padded to 20 field width
    static char line[512];
    static const int lastchari = Q_COUNTOF(line) - 2;
    const char  *s;
    int     i;
    ddef_t      *def;
    void        *val;

    val = (void *)&qcvm->globals[ofs];
    def = ED_GlobalAtOfs(ofs);
    if (!def)
        q_snprintf (line, sizeof(line), "%i(?)", ofs);
    else
    {
        s = PR_ValueString (def->type, (eval_t *)val);
        q_snprintf (line, sizeof(line), "%i(%s)%s", ofs, PR_GetString(def->s_name), s);
    }

    i = strlen(line);
    for ( ; i < 20; i++)
        strcat (line, " ");

    if (i < lastchari)
        strcat (line, " ");
    else
        line[lastchari] = ' ';

    return line;
}

const char *PR_GlobalStringNoContents (int ofs)
{
    static char line[512];
    static const int lastchari = Q_COUNTOF(line) - 2;
    int     i;
    ddef_t      *def;

    def = ED_GlobalAtOfs(ofs);
    if (!def)
        q_snprintf (line, sizeof(line), "%i(?)", ofs);
    else
        q_snprintf (line, sizeof(line), "%i(%s)", ofs, PR_GetString(def->s_name));

    i = strlen(line);
    for ( ; i < 20; i++)
        strcat (line, " ");

    if (i < lastchari)
        strcat (line, " ");
    else
        line[lastchari] = ' ';

    return line;
}

static void ED_AppendFlagString (char *dst, size_t dstsize, const char *desc)
{
    if (*dst)
        q_strlcat (dst, " | ", dstsize);
    q_strlcat (dst, desc, dstsize);
}

const char *ED_FieldValueString (edict_t *ed, ddef_t *d)
{
    static char str[1024];
    int ofs = d->ofs*4;
    eval_t *val = (eval_t *)((char *)&ed->v + ofs);

    // .movetype
    if (ofs == offsetof (entvars_t, movetype) && val->_float == (int)val->_float)
    {
        switch ((int)val->_float)
        {
        #define MOVETYPE_CASE(x)    case x: return #x
        MOVETYPE_CASE (MOVETYPE_NONE);
        MOVETYPE_CASE (MOVETYPE_ANGLENOCLIP);
        MOVETYPE_CASE (MOVETYPE_ANGLECLIP);
        MOVETYPE_CASE (MOVETYPE_WALK);
        MOVETYPE_CASE (MOVETYPE_STEP);
        MOVETYPE_CASE (MOVETYPE_FLY);
        MOVETYPE_CASE (MOVETYPE_TOSS);
        MOVETYPE_CASE (MOVETYPE_PUSH);
        MOVETYPE_CASE (MOVETYPE_NOCLIP);
        MOVETYPE_CASE (MOVETYPE_FLYMISSILE);
        MOVETYPE_CASE (MOVETYPE_BOUNCE);
        MOVETYPE_CASE (MOVETYPE_GIB);
        #undef MOVETYPE_CASE
        default:
            break;
        }
    }

    // .solid
    if (ofs == offsetof (entvars_t, solid) && val->_float == (int)val->_float)
    {
        switch ((int)val->_float)
        {
        #define SOLID_CASE(x)   case x: return #x
        SOLID_CASE (SOLID_NOT);
        SOLID_CASE (SOLID_TRIGGER);
        SOLID_CASE (SOLID_BBOX);
        SOLID_CASE (SOLID_SLIDEBOX);
        SOLID_CASE (SOLID_BSP);
        #undef SOLID_CASE
        default:
            break;
        }
    }

    // .deadflag
    if (ofs == offsetof (entvars_t, deadflag) && val->_float == (int)val->_float)
    {
        switch ((int)val->_float)
        {
        #define DEAD_CASE(x)    case x: return #x
        DEAD_CASE (DEAD_NO);
        DEAD_CASE (DEAD_DYING);
        DEAD_CASE (DEAD_DEAD);
        DEAD_CASE (DEAD_RESPAWNABLE);
        #undef DEAD_CASE
        default:
            break;
        }
    }

    // .takedamage
    if (ofs == offsetof (entvars_t, takedamage) && val->_float == (int)val->_float)
    {
        switch ((int)val->_float)
        {
        #define TAKEDAMAGE_CASE(x)  case x: return #x
        TAKEDAMAGE_CASE (DAMAGE_NO);
        TAKEDAMAGE_CASE (DAMAGE_YES);
        TAKEDAMAGE_CASE (DAMAGE_AIM);
        #undef TAKEDAMAGE_CASE
        default:
            break;
        }
    }

    // bitfield: .flags, .spawnflags, .effects
    if ((ofs == offsetof (entvars_t, flags) || ofs == offsetof (entvars_t, spawnflags) || ofs == offsetof (entvars_t, effects)) && val->_float == (int)val->_float)
    {
        int bits = (int)val->_float;
        str[0] = '\0';

        #define BIT_CASE(f) do { if (bits & (int)f) { bits ^= (int)f; ED_AppendFlagString (str, sizeof (str), #f); } } while (0)

        if (ofs == offsetof (entvars_t, flags))
        {
            BIT_CASE (FL_FLY);
            BIT_CASE (FL_CONVEYOR);
            BIT_CASE (FL_CLIENT);
            BIT_CASE (FL_INWATER);
            BIT_CASE (FL_MONSTER);
            BIT_CASE (FL_GODMODE);
            BIT_CASE (FL_NOTARGET);
            BIT_CASE (FL_ITEM);
            BIT_CASE (FL_ONGROUND);
            BIT_CASE (FL_PARTIALGROUND);
            BIT_CASE (FL_WATERJUMP);
            BIT_CASE (FL_JUMPRELEASED);
        }
        else if (ofs == offsetof (entvars_t, spawnflags))
        {
            BIT_CASE (SPAWNFLAG_NOT_EASY);
            BIT_CASE (SPAWNFLAG_NOT_MEDIUM);
            BIT_CASE (SPAWNFLAG_NOT_HARD);
            BIT_CASE (SPAWNFLAG_NOT_DEATHMATCH);
        }
        else if (ofs == offsetof (entvars_t, effects))
        {
            BIT_CASE (EF_BRIGHTFIELD);
            BIT_CASE (EF_MUZZLEFLASH);
            BIT_CASE (EF_BRIGHTLIGHT);
            BIT_CASE (EF_DIMLIGHT);
        }

        #undef BIT_CASE

        while (bits)
        {
            int lowest = bits & -bits;
            bits ^= lowest;
            ED_AppendFlagString (str, sizeof (str), va ("%d", lowest));
        }

        return str;
    }

    // .nextthink
    if (ofs == offsetof (entvars_t, nextthink) && val->_float)
    {
        return va (" %7.1f (%+.2f)", val->_float, val->_float - qcvm->time);
    }

    // generic field
    return PR_ValueString (d->type, val);
}

bool ED_IsRelevantField (edict_t *ed, ddef_t *d)
{
    const char  *name;
    size_t      l;
    int         *v;
    int         type;
    int         i;

    name = PR_GetString (d->s_name);
    l = strlen (name);
    if (l > 1 && name[l - 2] == '_')
        return false;   // skip _x, _y, _z vars

    type = d->type & ~DEF_SAVEGLOBAL;
    if (type >= NUM_TYPE_SIZES)
        return false;

    // if the value is still all 0, skip the field
    v = (int *)((char *)&ed->v + d->ofs*4);
    for (i = 0; i < type_size[type]; i++)
        if (v[i])
            return true;

    return false;
}

void ED_Print (edict_t *ed)
{ // For debugging
    ddef_t  *d;
    int     i, l;
    char    field[4096], buf[4096], *p;

    if (ed->free)
    {
        Con_SafePrintf ("FREE\n");
        return;
    }

    q_snprintf (buf, sizeof (buf), "\nEDICT %i:\n", NUM_FOR_EDICT(ed)); //johnfitz -- was Con_Printf
    p = buf + strlen (buf);
    for (i = 1; i < qcvm->progs->numfielddefs; i++)
    {
        d = &qcvm->fielddefs[i];
        if (!ED_IsRelevantField (ed, d))
            continue;

        q_snprintf (field, sizeof (field), "%-14s %s\n", PR_GetString (d->s_name), ED_FieldValueString (ed, d)); // johnfitz -- was Con_Printf
        l = strlen (field);
        if (l + 1 > buf + sizeof (buf) - p)
        {
            Con_SafePrintf ("%s", buf);
            p = buf;
        }

        memcpy (p, field, l + 1);
        p += l;
    }

    Con_SafePrintf ("%s", buf);
}


void ED_Write(FILE *f, edict_t *ed)
{ // For savegames
	fprintf(f, "{\n");
	if(ed->free) { fprintf(f, "}\n"); return; }
	for(s32 i = 1; i < qcvm->progs->numfielddefs; i++) {
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
	int     i;
	qcvm_t  *oldqcvm;
	if (!sv.active)
		return;
	PR_PushQCVM(&sv.qcvm, &oldqcvm);
	Con_Printf ("%i entities\n", qcvm->num_edicts);
	for (i = 0; i < qcvm->num_edicts; i++)
		ED_PrintNum (i);
	PR_PopQCVM(oldqcvm);
}

static void ED_PrintEdict_f()
{ // For debugging, prints a single edict
	if(!sv.active) return;
	s32 i = Q_atoi(Cmd_Argv(1));
	if(i < 0 || i >= qcvm->num_edicts) {
		Con_Printf("Bad edict number\n");
		return;
	}
	ED_PrintNum(i);
}

static void ED_Count()
{ // For debugging
	if(!sv.active) return;
	s32 active = 0, models = 0, solid = 0, step = 0;
	for(s32 i = 0; i < qcvm->num_edicts; i++) {
		edict_t *ent = EDICT_NUM(i);
		if(ent->free) continue;
		active++;
		if(ent->v.solid) solid++;
		if(ent->v.model) models++;
		if(ent->v.movetype == MOVETYPE_STEP) step++;
	}
	Con_Printf("num_edicts:%3i\n", qcvm->num_edicts);
	Con_Printf("active :%3i\n", active);
	Con_Printf("view :%3i\n", models);
	Con_Printf("touch :%3i\n", solid);
	Con_Printf("step :%3i\n", step);
}

void ED_WriteGlobals(FILE *f)
{
	fprintf(f, "{\n");
	for(s32 i = 0; i < qcvm->progs->numglobaldefs; i++) {
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

const char *ED_ParseGlobals (const char *data)
{
    char    keyname[64];
    ddef_t  *key;

    while (1)
    {
    // parse key
        data = COM_Parse (data);
        if (com_token[0] == '}')
            break;
        if (!data)
            Host_Error ("ED_ParseEntity: EOF without closing brace");

        q_strlcpy (keyname, com_token, sizeof(keyname));                                                              
    // parse value
        data = COM_Parse (data);
        if (!data)
            Host_Error ("ED_ParseEntity: EOF without closing brace");

        if (com_token[0] == '}')
            Host_Error ("ED_ParseEntity: closing brace without data");

        key = ED_FindGlobal (keyname);
        if (!key)
        {
            Con_Printf ("'%s' is not a global\n", keyname);
            continue;
        }

        if (!ED_ParseEpair ((void *)qcvm->globals, key, com_token, false))
            Host_Error ("ED_ParseGlobals: parse error");
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


static void ED_RezoneString (string_t *ref, const char *str)
{
    char *buf;
    size_t len = strlen(str)+1;
    size_t id;

    if (*ref)
    {   //if the reference is already a zoned string then free it first.
        id = -1-*ref;
        if (id < qcvm->knownzonesize && (qcvm->knownzone[id>>3] & (1u<<(id&7))))
        {   //okay, it was zoned.
            qcvm->knownzone[id>>3] &= ~(1u<<(id&7));
            buf = (char*)PR_GetString(*ref);
            PR_ClearEngineString(*ref);
            Z_Free(buf);
        }
//      else
//          Con_Warning("ED_RezoneString: string wasn't strzoned\n");   //warnings would trigger from the default cvar value that autocvars are initialised with
    }

    buf = Z_Malloc(len);
    memcpy(buf, str, len);
    id = -1-(*ref = PR_SetEngineString(buf));
    //make sure its flagged as zoned so we can clean up properly after.
    if (id >= qcvm->knownzonesize)
    {
        qcvm->knownzonesize = (id+32)&~7;
        qcvm->knownzone = Z_Realloc(qcvm->knownzone, (qcvm->knownzonesize+7)>>3);
    }
    qcvm->knownzone[id>>3] |= 1u<<(id&7);
}

// Can parse either fields or globals
// returns false if error
static bool ED_ParseEpair (void *base, ddef_t *key, const char *s, bool zoned)
{
    int     i;
    char    string[128];
    ddef_t  *def;
    char    *v, *w;
    char    *end;
    void    *d;
    dfunction_t *func;

    d = (void *)((int *)base + key->ofs);

    switch (key->type & ~DEF_SAVEGLOBAL)
    {
    case ev_string:
        if (zoned)  //zoned version allows us to change the strings more freely
            ED_RezoneString((string_t *)d, s);
        else
            *(string_t *)d = ED_NewString(s);
        break;

    case ev_float:
        *(float *)d = atof (s);
        break;

    case ev_vector:
        q_strlcpy (string, s, sizeof(string));
        end = (char *)string + strlen(string);
        v = string;
        w = string;

        for (i = 0; i < 3 && (w <= end); i++) // ericw -- added (w <= end) check
        {
        // set v to the next space (or 0 byte), and change that char to a 0 byte
            while (*v && *v != ' ')
                v++;
            *v = 0;
            ((float *)d)[i] = atof (w);
            w = v = v+1;
        }
        // ericw -- fill remaining elements to 0 in case we hit the end of string
        // before reading 3 floats.
        if (i < 3)
        {
            Con_DPrintf ("Avoided reading garbage for \"%s\" \"%s\"\n", PR_GetString(key->s_name), s);
            for (; i < 3; i++)
                ((float *)d)[i] = 0.0f;
        }
        break;

    case ev_entity:
        *(int *)d = EDICT_TO_PROG(EDICT_NUM(atoi (s)));
        break;

    case ev_field:
        def = ED_FindField (s);
        if (!def)
        {
            //johnfitz -- HACK -- suppress error becuase fog/sky fields might not be mentioned in defs.qc
            if (strncmp(s, "sky", 3) && strcmp(s, "fog"))
                Con_DPrintf ("Can't find field %s\n", s);
            return false;
        }
        *(int *)d = G_INT(def->ofs);
        break;

    case ev_function:
        func = ED_FindFunction (s);
        if (!func)
        {
            Con_Printf ("Can't find function %s\n", s);
            return false;
        }
        *(func_t *)d = func - qcvm->functions;
        break;

    default:
        break;
    }
    return true;
}

// Parses an edict out of the given string, returning the new position
// ed should be a properly initialized empty edict.
// Used for initial level load and for savegames.
const char *ED_ParseEdict (const char *data, edict_t *ent)
{
    ddef_t      *key;
    char        keyname[256];
    bool    anglehack, init;
    int     n;

    init = false;

    // clear it
    if (ent != qcvm->edicts)    // hack
        memset (&ent->v, 0, qcvm->progs->entityfields * 4);

    // go through all the dictionary pairs
    while (1)
    {
        // parse key
        data = COM_Parse (data);
        if (com_token[0] == '}')
            break;
        if (!data)
            Host_Error ("ED_ParseEntity: EOF without closing brace");

        // anglehack is to allow QuakeEd to write single scalar angles
        // and allow them to be turned into vectors. (FIXME...)
        if (!strcmp(com_token, "angle"))
        {
            strcpy (com_token, "angles");
            anglehack = true;
        }
        else
            anglehack = false;

        // FIXME: change light to _light to get rid of this hack
        if (!strcmp(com_token, "light"))
            strcpy (com_token, "light_lev");    // hack for single light def

        q_strlcpy (keyname, com_token, sizeof(keyname));

        // another hack to fix keynames with trailing spaces
        n = strlen(keyname);
        while (n && keyname[n-1] == ' ')
        {
            keyname[n-1] = 0;
            n--;
        }

        // parse value
        // HACK: we allow truncation when reading the wad field,
        // otherwise maps using lots of wads with absolute paths
        // could cause a parse error
        data = COM_ParseEx (data, !strcmp (keyname, "wad") ? CPE_ALLOWTRUNC : CPE_NOTRUNC);
        if (!data)
            Host_Error ("ED_ParseEntity: EOF without closing brace");

        if (com_token[0] == '}')
            Host_Error ("ED_ParseEntity: closing brace without data");

        init = true;

        // keynames with a leading underscore are used for utility comments,
        // and are immediately discarded by quake
        if (keyname[0] == '_')
            continue;

        //johnfitz -- hack to support .alpha even when progs.dat doesn't know about it
        if (!strcmp(keyname, "alpha"))
            ent->alpha = ENTALPHA_ENCODE(Q_atof(com_token));
        //johnfitz

        key = ED_FindField (keyname);
        if (!key)
        {
            //johnfitz -- HACK -- suppress error becuase fog/sky/alpha fields might not be mentioned in defs.qc
            if (strncmp(keyname, "sky", 3) && strcmp(keyname, "fog") && strcmp(keyname, "alpha"))
                Con_DPrintf ("\"%s\" is not a field\n", keyname); //johnfitz -- was Con_Printf
            continue;
        }

        if (anglehack)
        {
            char    temp[32];
            strcpy (temp, com_token);
            sprintf (com_token, "0 %s 0", temp);
        }

        if (!ED_ParseEpair ((void *)&ent->v, key, com_token, qcvm != &sv.qcvm))
            Host_Error ("ED_ParseEdict: parse error");
    }

    if (!init)
        ED_Free (ent);

    return data;
}

static bool ED_IsSkillSelector (const edict_t *ent)
{
    int skill;
    const char *classname = PR_GetString (ent->v.classname);

    if (strcmp (classname, "trigger_setskill") == 0 || strcmp (classname, "target_setskill") == 0)
        return true;
    if (strcmp (classname, "info_command") == 0 && (int)ent->v.message != 0 && sscanf (PR_GetString (ent->v.message), "skill %d", &skill) == 1)
        return true;

    return false;
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
	const char  *classname;
	dfunction_t *func;
	edict_t     *ent = NULL;
	int     inhibit = 0;

	pr_global_struct->time = qcvm->time;

	// parse ents
	while (1)
	{
		// parse the opening brace
		data = COM_Parse (data);
		if (!data)
			break;
		if (com_token[0] != '{')
			Host_Error ("ED_LoadFromFile: found %s when expecting {",com_token);

		if (!ent)
			ent = EDICT_NUM(0);
		else
			ent = ED_Alloc ();
		data = ED_ParseEdict (data, ent);

		if (!ent->v.classname)
		{
			Con_DPrintf ("No classname for:\n"); //johnfitz -- was Con_Printf
			ED_Print (ent);
			ED_Free (ent);
			continue;
		}

		classname = PR_GetString (ent->v.classname);

		if (sv.mapchecks.active)
		{
			int skillflags = (int)ent->v.spawnflags & (SPAWNFLAG_NOT_EASY|SPAWNFLAG_NOT_MEDIUM|SPAWNFLAG_NOT_HARD);
			if (!(skillflags & SPAWNFLAG_NOT_EASY))
				sv.mapchecks.skill_ents[0]++;
			if (!(skillflags & SPAWNFLAG_NOT_MEDIUM))
				sv.mapchecks.skill_ents[1]++;
			if (!(skillflags & SPAWNFLAG_NOT_HARD))
				sv.mapchecks.skill_ents[2]++;

			if (strcmp (classname, "trigger_changelevel") == 0)
			{
				ddef_t *mapfield = ED_FindField ("map");
				sv.mapchecks.trigger_changelevel++;
				if (mapfield && (mapfield->type & ~DEF_SAVEGLOBAL) == ev_string)
				{
					eval_t      *val = GetEdictFieldValue (ent, mapfield->ofs);
					const char  *map = COM_SkipSpace (PR_GetString (val->string));
					if (*map)
					{
						sv.mapchecks.changelevel = map;
						sv.mapchecks.valid_changelevel++;
					}
				}
			}
			else if (ED_IsSkillSelector (ent))
				sv.mapchecks.skill_triggers++;
			else if (strcmp (classname, "info_intermission") == 0)
				sv.mapchecks.intermission++;
			else if (strcmp (classname, "info_player_coop") == 0)
				sv.mapchecks.coop_spawns++;
			else if (strcmp (classname, "info_player_deathmatch") == 0)
				sv.mapchecks.dm_spawns++;
		}

		// remove things from different skill levels or deathmatch
		if (deathmatch.value)
		{
			if (((int)ent->v.spawnflags & SPAWNFLAG_NOT_DEATHMATCH))
			{
				ED_Free (ent);
				inhibit++;
				continue;
			}
		}
		else if ((current_skill == 0 && ((int)ent->v.spawnflags & SPAWNFLAG_NOT_EASY))
				|| (current_skill == 1 && ((int)ent->v.spawnflags & SPAWNFLAG_NOT_MEDIUM))
				|| (current_skill >= 2 && ((int)ent->v.spawnflags & SPAWNFLAG_NOT_HARD)) )
		{
			ED_Free (ent);
			inhibit++;
			continue;
		}

		// remove monsters if nomonsters is set
		if (sv.nomonsters && !Q_strncmp (classname, "monster_", 8))
		{
			ED_Free (ent);
			inhibit++;
			continue;
		}


		//
		// immediately call spawn function
		//
		// look for the spawn function
		func = ED_FindFunction (classname);

		if (!func)
		{
			Con_DPrintf ("No spawn function for:\n"); //johnfitz -- was Con_Printf
			ED_Print (ent);
			ED_Free (ent);
			continue;
		}

		SV_ReserveSignonSpace (512);

		pr_global_struct->self = EDICT_TO_PROG(ent);
		PR_ExecuteProgram (func - qcvm->functions);
	}

	Con_DPrintf ("%i entities inhibited\n", inhibit);
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

void PR_UnzoneAll(void)
{   //called to clean up all zoned strings.
    while (qcvm->knownzonesize --> 0)
    {
        size_t id = qcvm->knownzonesize;
        if (qcvm->knownzone[id>>3] & (1u<<(id&7)))
        {
            string_t s = -1-(int)id;
            char *ptr = (char*)PR_GetString(s);
            PR_ClearEngineString(s);
            Z_Free(ptr);
        }
    }
    if (qcvm->knownzone)
        Z_Free(qcvm->knownzone);
    qcvm->knownzonesize = 0;
    qcvm->knownzone = NULL;
}

//makes sure extension fields are actually registered so they can be used for mappers without qc changes. eg so scale can be used.
static void PR_MergeEngineFieldDefs (void)
{
	struct {
		const char *fname;
		etype_t type;
		int newidx;
	} extrafields[] =
	{   //table of engine fields to add. we'll be using ED_FindFieldOffset for these later.
	    //this is useful for fields that should be defined for mappers which are not defined by the mod.
	    //future note: mutators will need to edit the mutator's globaldefs table too. remember to handle vectors and their 3 globals too.
		{"alpha",           ev_float},  //just because we can (though its already handled in a weird hacky way)
		{"scale",           ev_float},  //hurrah for being able to rescale entities.
		{"emiteffectnum",   ev_float},  //constantly emitting particles, even without moving.
		{"traileffectnum",  ev_float},  //custom effect for trails
						//{"glow_size",     ev_float},  //deprecated particle trail rubbish
						//{"glow_color",    ev_float},  //deprecated particle trail rubbish
		{"tag_entity",      ev_float},  //for setattachment to not bug out when omitted.
		{"tag_index",       ev_float},  //for setattachment to not bug out when omitted.
		{"modelflags",      ev_float},  //deprecated rubbish to fill the high 8 bits of effects.
						//{"vw_index",      ev_float},  //modelindex2
						//{"pflags",        ev_float},  //for rtlights
						//{"drawflags",     ev_float},  //hexen2 compat
						//{"abslight",      ev_float},  //hexen2 compat
		{"colormod",        ev_vector}, //lighting tints
						//{"glowmod",       ev_vector}, //fullbright tints
						//{"fatness",       ev_float},  //bloated rendering...
						//{"gravitydir",    ev_vector}, //says which direction gravity should act for this ent...

	};
	int maxofs = qcvm->progs->entityfields;
	int maxdefs = qcvm->progs->numfielddefs;
	unsigned int j, a;

	//figure out where stuff goes
	for (j = 0; j < Q_COUNTOF(extrafields); j++)
	{
		extrafields[j].newidx = ED_FindFieldOffset(extrafields[j].fname);
		if (extrafields[j].newidx < 0)
		{
			extrafields[j].newidx = maxofs;
			maxdefs++;
			if (extrafields[j].type == ev_vector)
				maxdefs+=3;
			maxofs+=type_size[extrafields[j].type];
		}
	}

	if (maxdefs != qcvm->progs->numfielddefs)
	{   //we now know how many entries we need to add...
		ddef_t *olddefs = qcvm->fielddefs;
		qcvm->fielddefs = malloc(maxdefs * sizeof(*qcvm->fielddefs));
		if (!qcvm->fielddefs)
			Sys_Error ("PR_MergeEngineFieldDefs: out of memory (%d defs)", maxdefs);
		memcpy(qcvm->fielddefs, olddefs, qcvm->progs->numfielddefs*sizeof(*qcvm->fielddefs));
		if (olddefs != (ddef_t *)((u8 *)qcvm->progs + qcvm->progs->ofs_fielddefs))
			free(olddefs);

		//allocate the extra defs
		for (j = 0; j < Q_COUNTOF(extrafields); j++)
		{
			if (extrafields[j].newidx >= qcvm->progs->entityfields && extrafields[j].newidx < maxofs)
			{   //looks like its new. make sure ED_FindField can find it.
				qcvm->fielddefs[qcvm->progs->numfielddefs].ofs = extrafields[j].newidx;
				qcvm->fielddefs[qcvm->progs->numfielddefs].type = extrafields[j].type;
				qcvm->fielddefs[qcvm->progs->numfielddefs].s_name = ED_NewString(extrafields[j].fname);
				qcvm->progs->numfielddefs++;

				if (extrafields[j].type == ev_vector)
				{   //vectors are weird and annoying.
					for (a = 0; a < 3; a++)
					{
						qcvm->fielddefs[qcvm->progs->numfielddefs].ofs = extrafields[j].newidx+a;
						qcvm->fielddefs[qcvm->progs->numfielddefs].type = ev_float;
						qcvm->fielddefs[qcvm->progs->numfielddefs].s_name = ED_NewString(va("%s_%c", extrafields[j].fname, 'x'+a));
						qcvm->progs->numfielddefs++;
					}
				}
			}
		}
		qcvm->progs->entityfields = maxofs;
	}
}

void PR_ShutdownExtensions (void)
{ // called at map end
    PR_UnzoneAll();
    if (qcvm == &cl.qcvm)
        PR_ReloadPics(true);
}

void PR_ClearProgs(qcvm_t *vm)
{
    qcvm_t *oldvm = qcvm;
    if (!vm->progs)
        return; //wasn't loaded.
    if (vm == &sv.qcvm)
        NULL;//Host_WaitForSaveThread ();
    qcvm = NULL;
    PR_SwitchQCVM(vm);
    PR_ShutdownExtensions();

    if (qcvm->knownstrings)
        Z_Free ((void *)qcvm->knownstrings);
    free(qcvm->edicts); // ericw -- sv.edicts switched to use malloc()
    if (qcvm->fielddefs != (ddef_t *)((u8 *)qcvm->progs + qcvm->progs->ofs_fielddefs))
        free(qcvm->fielddefs);
    memset(qcvm, 0, sizeof(*qcvm));

    qcvm = NULL;
    PR_SwitchQCVM(oldvm);
}
static void PR_HashAdd (prhashtable_t *table, int skey, int value)
{
    const char *name = PR_GetString (skey);
    unsigned pos = COM_HashString (name) % table->capacity, end = pos;

    do
    {
        if (!table->strings[pos])
        {
            table->strings[pos] = name;
            table->indices[pos] = value;
            return;
        }

        ++pos;
        if (pos == table->capacity)
            pos = 0;
    }
    while (pos != end);

    Sys_Error ("PR_HashAdd failed");
}

static void PR_InitHashTables (void)
{
    int i;

    PR_HashInit (&qcvm->ht_fields, qcvm->progs->numfielddefs, "ht_fields");
    for (i = 0; i < qcvm->progs->numfielddefs; i++)
        PR_HashAdd (&qcvm->ht_fields, qcvm->fielddefs[i].s_name, i);

    PR_HashInit (&qcvm->ht_functions, qcvm->progs->numfunctions, "ht_functions");
    for (i = 0; i < qcvm->progs->numfunctions; i++)
        PR_HashAdd (&qcvm->ht_functions, qcvm->functions[i].s_name, i);

    PR_HashInit (&qcvm->ht_globals, qcvm->progs->numglobaldefs, "ht_globals");
    for (i = 0; i < qcvm->progs->numglobaldefs; i++)
        PR_HashAdd (&qcvm->ht_globals, qcvm->globaldefs[i].s_name, i);
}

static void PR_InitBuiltins (void)
{
    dfunction_t *func;
    const char  *name;
    int         i, j;

    for (i = 0; i < MAX_BUILTINS; i++)
        qcvm->builtins[i] = PF_Fixme;

    for (i = MAX_BUILTINS - 2, j = 0; j < pr_numbuiltindefs; j++)
    {
        builtindef_t *def = &pr_builtindefs[j];
        builtin_t func = (qcvm == &sv.qcvm) ? def->ssqcfunc : def->csqcfunc;
        if (!def->number)
            def->number = i--;
        if (func)
        {
            qcvm->builtins[def->number] = func;
            qcvm->builtin_ext[def->number] = def->ext;
        }
    }

    qcvm->numbuiltins = MAX_BUILTINS;

    // remap progs functions with id 0
    for (i = 0; i < qcvm->progs->numfunctions; i++)
    {
        func = &qcvm->functions[i];
        if (func->first_statement || func->parm_start || func->locals)
            continue;

        name = PR_GetString (func->s_name);
        for (j = 0; j < pr_numbuiltindefs; j++)
        {
            builtindef_t *def = &pr_builtindefs[j];
            if (!strcmp (name, def->name))
            {
                func->first_statement = -def->number;
                break;
            }
        }
    }
}

static void PR_FindSavegameFields (void)
{ // Determines which fields should be stored in savefiles
    int i;
    for (i = 1; i < qcvm->progs->numfielddefs; i++)
    {
        ddef_t *field = &qcvm->fielddefs[i];
        const char *name = PR_GetString (field->s_name);
        size_t len = strlen (name);
        if (len < 2 || name[len - 2] != '_') // skip _x, _y, _z vars
            field->type |= DEF_SAVEGLOBAL;
    }
}

static int PR_CompareFunction (const void *pa, const void *pb)
{
    const dfunction_t *fa = &qcvm->functions[*(const int *)pa];
    const dfunction_t *fb = &qcvm->functions[*(const int *)pb];
    return fa->first_statement - fb->first_statement;
}

static void PR_FindFunctionRanges (void)
{
    int     i, mark;
    int     *order;

    qcvm->functionsizes = (int *) Hunk_AllocName (qcvm->progs->numfunctions * sizeof (*order), "func_sizes");
    mark = Hunk_LowMark ();

    order = (int *) Hunk_Alloc (qcvm->progs->numfunctions * sizeof (*order));
    for (i = 0; i < qcvm->progs->numfunctions; i++)
        order[i] = i;
    qsort (order, qcvm->progs->numfunctions, sizeof (*order), &PR_CompareFunction);

    for (i = 0; i < qcvm->progs->numfunctions; i++)
    {
        dfunction_t *f = &qcvm->functions[order[i]];
        if (f->first_statement <= 0)
            continue;
        if (i == qcvm->progs->numfunctions - 1)
            qcvm->functionsizes[order[i]] = qcvm->progs->numstatements - f->first_statement;
        else
            qcvm->functionsizes[order[i]] = qcvm->functions[order[i + 1]].first_statement - f->first_statement;
    }

    Hunk_FreeToLowMark (mark);
}

static void PR_FindEntityFields (void)
{ // Finds all the .entity fields (used by r_showbboxes & co for identifying entity links)
    int i, count;

    count = 0;
    for (i = 1; i < qcvm->progs->numfielddefs; i++)
    {
        ddef_t *field = &qcvm->fielddefs[i];
        if ((field->type & ~DEF_SAVEGLOBAL) == ev_entity)
            count++;
    }

    qcvm->numentityfields = count;
    qcvm->entityfieldofs = (int *) Hunk_AllocName (qcvm->numentityfields * sizeof (int), "entityfieldofs");
    qcvm->entityfields = (ddef_t **) Hunk_AllocName (qcvm->numentityfields * sizeof (ddef_t*), "entityfields");

    count = 0;
    for (i = 1; i < qcvm->progs->numfielddefs; i++)
    {
        ddef_t *field = &qcvm->fielddefs[i];
        if ((field->type & ~DEF_SAVEGLOBAL) == ev_entity)
        {
            qcvm->entityfieldofs[count] = field->ofs*4;
            qcvm->entityfields[count] = field;
            count++;
        }
    }
}

static void PR_FillOffsetTables (void)
{
    int     pass, i, maxofs, *data;
    struct
    {
        int         **offsets;
        int         *maxofs;
        int         numdefs;
        ddef_t      *defs;
        const char  *allocname;
    }
    passes[] =
    {
        { &qcvm->ofstofield,    &qcvm->maxfieldofs,     qcvm->progs->numfielddefs,  qcvm->fielddefs,    "ofstofield"    },
        { &qcvm->ofstoglobal,   &qcvm->maxglobalofs,    qcvm->progs->numglobaldefs, qcvm->globaldefs,   "ofstoglobal"   },
    };

    for (pass = 0; pass < (int) Q_COUNTOF (passes); pass++)
    {
        // find maximum offset
        for (i = 1, maxofs = 0; i < passes[pass].numdefs; i++)
            maxofs = q_max (maxofs, passes[pass].defs[i].ofs);
        *passes[pass].maxofs = maxofs;

        // alloc table and fill it with -1
        data = *passes[pass].offsets = (int *) Hunk_AllocName ((maxofs + 1) * sizeof (int), passes[pass].allocname);
        for (i = 0; i <= maxofs; i++)
            data[i] = -1;

        // fill actual offsets in descending order so that earlier defs are written last
        // this preserves the behavior of ED_FieldAtOfs/ED_GlobalAtOfs, which stopped at the first match
        for (i = passes[pass].numdefs - 1; i > 0; i--)
            data[passes[pass].defs[i].ofs] = i;
    }
}

bool PR_LoadProgs(const s8 *filename, bool fatal)
{
	int         i;

	PR_ClearProgs(qcvm);    //just in case.

	qcvm->progs = (dprograms_t *)COM_LoadHunkFile (filename, NULL);
	if (!qcvm->progs)
		return false;
	Con_DPrintf ("Programs occupy %" SDL_PRIs64 "K.\n", com_filesize/1024);

	qcvm->crc = CRC_Block ((u8*)qcvm->progs, com_filesize);

	// byte swap the header
	for (i = 0; i < (int) sizeof(*qcvm->progs) / 4; i++)
		((int *)qcvm->progs)[i] = LittleLong ( ((int *)qcvm->progs)[i] );

	if (qcvm->progs->version != PROG_VERSION) {
		if (fatal)
			Host_Error ("%s has wrong version number (%i should be %i)", filename, qcvm->progs->version, PROG_VERSION);
		else
		{
			Con_Printf("%s ABI set not supported\n", filename);
			qcvm->progs = NULL;
			return false;
		}
	}

	if (qcvm->progs->crc != PROGHEADER_CRC)
	{
		if (fatal)
			Host_Error ("%s system vars have been modified, progdefs.h is out of date", filename);
		else
		{
			switch(qcvm->progs->crc)
			{
				case 22390: //full csqc
					Con_Printf("%s - full csqc is not supported\n", filename);
					break;
				case 52195: //dp csqc
					Con_Printf("%s - obsolete csqc is not supported\n", filename);
					break;
				case 54730: //quakeworld
					Con_Printf("%s - quakeworld gamecode is not supported\n", filename);
					break;
				case 26940: //prerelease
					Con_Printf("%s - prerelease gamecode is not supported\n", filename);
					break;
				case 32401: //tenebrae
					Con_Printf("%s - tenebrae gamecode is not supported\n", filename);
					break;
				case 38488: //hexen2 release
				case 26905: //hexen2 mission pack
				case 14046: //hexen2 demo
					Con_Printf("%s - hexen2 gamecode is not supported\n", filename);
					break;
					//case 5927: //nq PROGHEADER_CRC as above. shouldn't happen, obviously.
				default:
					Con_Printf("%s system vars are not supported\n", filename);
					break;
			}
			qcvm->progs = NULL;
			return false;
		}
	}

	qcvm->functions = (dfunction_t *)((u8 *)qcvm->progs + qcvm->progs->ofs_functions);
	qcvm->strings = (char *)qcvm->progs + qcvm->progs->ofs_strings;
	if (qcvm->progs->ofs_strings + qcvm->progs->numstrings >= com_filesize)
		Host_Error ("progs.dat strings go past end of file\n");

	// initialize the strings
	qcvm->numknownstrings = 0;
	qcvm->maxknownstrings = 0;
	qcvm->stringssize = qcvm->progs->numstrings;
	if (qcvm->knownstrings)
		Z_Free ((void *)qcvm->knownstrings);
	qcvm->knownstrings = NULL;
	qcvm->firstfreeknownstring = NULL;
	PR_SetEngineString("");

	qcvm->globaldefs = (ddef_t *)((u8 *)qcvm->progs + qcvm->progs->ofs_globaldefs);
	qcvm->fielddefs = (ddef_t *)((u8 *)qcvm->progs + qcvm->progs->ofs_fielddefs);
	qcvm->statements = (dstatement_t *)((u8 *)qcvm->progs + qcvm->progs->ofs_statements);

	qcvm->globals = (float *)((u8 *)qcvm->progs + qcvm->progs->ofs_globals);
	pr_global_struct = (globalvars_t*)qcvm->globals;

	// byte swap the lumps
	for (i = 0; i < qcvm->progs->numstatements; i++)
	{
		qcvm->statements[i].op = LittleShort(qcvm->statements[i].op);
		qcvm->statements[i].a = LittleShort(qcvm->statements[i].a);
		qcvm->statements[i].b = LittleShort(qcvm->statements[i].b);
		qcvm->statements[i].c = LittleShort(qcvm->statements[i].c);
	}

	for (i = 0; i < qcvm->progs->numfunctions; i++)
	{
		qcvm->functions[i].first_statement = LittleLong (qcvm->functions[i].first_statement);
		qcvm->functions[i].parm_start = LittleLong (qcvm->functions[i].parm_start);
		qcvm->functions[i].s_name = LittleLong (qcvm->functions[i].s_name);
		qcvm->functions[i].s_file = LittleLong (qcvm->functions[i].s_file);
		qcvm->functions[i].numparms = LittleLong (qcvm->functions[i].numparms);
		qcvm->functions[i].locals = LittleLong (qcvm->functions[i].locals);
	}

	for (i = 0; i < qcvm->progs->numglobaldefs; i++)
	{
		qcvm->globaldefs[i].type = LittleShort (qcvm->globaldefs[i].type);
		qcvm->globaldefs[i].ofs = LittleShort (qcvm->globaldefs[i].ofs);
		qcvm->globaldefs[i].s_name = LittleLong (qcvm->globaldefs[i].s_name);
	}

	for (i = 0; i < qcvm->progs->numfielddefs; i++)
	{
		qcvm->fielddefs[i].type = LittleShort (qcvm->fielddefs[i].type);
		if (qcvm->fielddefs[i].type & DEF_SAVEGLOBAL)
			Host_Error ("PR_LoadProgs: pr_fielddefs[i].type & DEF_SAVEGLOBAL");
		qcvm->fielddefs[i].ofs = LittleShort (qcvm->fielddefs[i].ofs);
		qcvm->fielddefs[i].s_name = LittleLong (qcvm->fielddefs[i].s_name);
	}

	for (i = 0; i < qcvm->progs->numglobals; i++)
		((int *)qcvm->globals)[i] = LittleLong (((int *)qcvm->globals)[i]);

	//spike: detect extended fields from progs
	PR_MergeEngineFieldDefs ();
#define QCEXTFIELD(n,t) qcvm->extfields.n = ED_FindFieldOffset (#n);
	QCEXTFIELDS_ALL
		QCEXTFIELDS_GAME
		QCEXTFIELDS_SS
#undef QCEXTFIELD

		qcvm->edict_size = qcvm->progs->entityfields * 4 + sizeof(edict_t) - sizeof(entvars_t);
	// round off to next highest whole word address (esp for Alpha)
	// this ensures that pointers in the engine data area are always
	// properly aligned
	qcvm->edict_size += sizeof(void *) - 1;
	qcvm->edict_size &= ~(sizeof(void *) - 1);

	PR_InitHashTables ();
	PR_InitBuiltins ();
	PR_PatchRereleaseBuiltins ();
	PR_EnableExtensions ();
	PR_FindSavegameFields ();
	PR_FindEntityFields ();
	PR_FindFunctionRanges ();
	PR_FillOffsetTables ();

	qcvm->effects_mask = PR_FindSupportedEffects ();

	return true;
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
	if (n < 0 || n >= qcvm->max_edicts)
		Host_Error ("EDICT_NUM: bad number %i", n);
	return (edict_t *)((u8 *)qcvm->edicts + (n)*qcvm->edict_size);
}

int NUM_FOR_EDICT(edict_t *e)
{
    int     b;

    b = (u8 *)e - (u8 *)qcvm->edicts;
    b = b / qcvm->edict_size;

    if (b < 0 || b >= qcvm->num_edicts)
        Host_Error ("NUM_FOR_EDICT: bad pointer");
    return b;
}

static void PR_AllocStringSlots()
{
	pr_maxknownstrings += PR_STRING_ALLOCSLOTS;
	pr_knownstrings = (const s8 **) Z_Realloc((void *)pr_knownstrings,
					pr_maxknownstrings * sizeof(s8 *));
}

const char *PR_GetString (int num)
{
    if (num >= 0 && num < qcvm->stringssize)
        return qcvm->strings + num;
    else if (num < 0 && num >= -qcvm->numknownstrings)
    {
        if (!qcvm->knownstrings[-1 - num])
        {
            Host_Error ("PR_GetString: attempt to get a non-existant string %d\n", num);
            return "";
        }
        return qcvm->knownstrings[-1 - num];
    }
    else
    {
        Host_Error("PR_GetString: invalid string offset %d\n", num);
        return "";
    }
}
#define PR_STRING_ALLOCSLOTS    256

static int PR_AllocStringSlot (void)
{
    ptrdiff_t i;

    if (qcvm->firstfreeknownstring)
    {
        i = qcvm->firstfreeknownstring - qcvm->knownstrings;
        if (i < 0 || i >= qcvm->maxknownstrings)
            Sys_Error ("PR_AllocStringSlot failed: invalid free list index %" SDL_PRIs64 "/%i\n", (int64_t)i, qcvm->maxknownstrings);
        qcvm->firstfreeknownstring = (const char **) *qcvm->firstfreeknownstring;
    }
    else
    {
        i = qcvm->numknownstrings++;
        if (i >= qcvm->maxknownstrings)
        {
            qcvm->maxknownstrings += PR_STRING_ALLOCSLOTS;
            Con_DPrintf ("PR_AllocStringSlot: realloc'ing for %d slots\n", qcvm->maxknownstrings);
            qcvm->knownstrings = (const char **) Z_Realloc ((void *)qcvm->knownstrings, qcvm->maxknownstrings * sizeof(char *));
        }
    }

    return (int)i;
}
int PR_SetEngineString (const char *s)
{
    int     i;

    if (!s)
        return 0;
    if (s >= qcvm->strings && s <= qcvm->strings + qcvm->stringssize - 2)
        return (int)(s - qcvm->strings);
    for (i = 0; i < qcvm->numknownstrings; i++)
    {
        if (qcvm->knownstrings[i] == s)
            return -1 - i;
    }
    // new unknown engine string
    //Con_DPrintf ("PR_SetEngineString: new engine string %p\n", s);
    i = PR_AllocStringSlot ();
    qcvm->knownstrings[i] = s;
    return -1 - i;
}

int PR_AllocString (int size, char **ptr)
{
    int     i;

    if (!size)
        return 0;
    i = PR_AllocStringSlot ();
    qcvm->knownstrings[i] = (char *)Hunk_AllocName(size, "string");
    if (ptr)
        *ptr = (char *) qcvm->knownstrings[i];
    return -1 - i;
}
