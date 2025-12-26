// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
// cvar.c -- dynamic variable tracking
#include "quakedef.h"

static cvar_t *cvar_vars;
static s8 cvar_null_string[] = "";

void Cvar_List_f()
{
	cvar_t *cvar;
	s32 len = 0, count = 0;
	const s8 *partial = NULL;
	if(Cmd_Argc() > 1) {
		partial = Cmd_Argv(1);
		len = strlen(partial);
	}
	for(cvar = cvar_vars ; cvar ; cvar = cvar->next) {
		if(partial && strncmp(partial, cvar->name, len)) continue;
		Con_SafePrintf("%s%s %s \"%s\"\n",
				(cvar->flags & CVAR_ARCHIVE) ? "*" : " ",
				(cvar->flags & CVAR_NOTIFY) ? "s" : " ",
				cvar->name,
				cvar->string);
		count++;
	}
	Con_SafePrintf("%i cvars", count);
	if(partial) { Con_SafePrintf(" beginning with \"%s\"", partial); }
	Con_SafePrintf("\n");
}

void Cvar_Inc_f()
{
	switch(Cmd_Argc()) {
	default:
	case 1:
		Con_Printf("inc <cvar> [amount] : increment cvar\n");
		break;
	case 2:
		Cvar_SetValue(Cmd_Argv(1), Cvar_VariableValue(Cmd_Argv(1)) + 1);
		break;
	case 3:
Cvar_SetValue(Cmd_Argv(1),Cvar_VariableValue(Cmd_Argv(1))+Q_atof(Cmd_Argv(2)));
		break;
	}
}

void Cvar_Toggle_f()
{
	switch(Cmd_Argc()) {
	default:
	case 1:
		Con_Printf("toggle <cvar> : toggle cvar\n");
		break;
	case 2:
		if(Cvar_VariableValue(Cmd_Argv(1)))
			Cvar_Set(Cmd_Argv(1), "0");
		else    Cvar_Set(Cmd_Argv(1), "1");
		break;
	}
}

void Cvar_Cycle_f()
{

	if(Cmd_Argc() < 3) {
Con_Printf("cycle <cvar> <value list>: cycle cvar through a list of values\n");
		return; }
// loop through the args until you find one that matches the current cvar value.
// yes, this will get stuck on a list that contains the same value twice.
// it's not worth dealing with, and i'm not even sure it can be dealt with.
	s32 i = 2;
	for(; i < Cmd_Argc(); i++) {
// zero is assumed to be a string, even though it could actually be zero. The
// worst case is that the first time you call this command, it won't match on
// zero when it should, but after that, it will be comparing strings that all
// had the same source(the user) so it will work.
		if(Q_atof(Cmd_Argv(i)) == 0) {
		      if(!strcmp(Cmd_Argv(i), Cvar_VariableString(Cmd_Argv(1))))
				break;
		}else if(Q_atof(Cmd_Argv(i)) == Cvar_VariableValue(Cmd_Argv(1)))
				break;
	}
	if(i == Cmd_Argc())
		Cvar_Set(Cmd_Argv(1), Cmd_Argv(2)); // no match
	else if(i + 1 == Cmd_Argc())
		Cvar_Set(Cmd_Argv(1), Cmd_Argv(2)); //matched last value in list
	else
		Cvar_Set(Cmd_Argv(1), Cmd_Argv(i+1)); // matched earlier in list
}

void Cvar_Reset(const s8 *name)
{
	cvar_t *var = Cvar_FindVar(name);
	if(!var) Con_Printf("variable \"%s\" not found\n", name);
	else Cvar_SetQuick(var, var->default_string);
}

void Cvar_Reset_f()
{
	switch(Cmd_Argc()) {
	default:
	case 1: Con_Printf("reset <cvar> : reset cvar to default\n"); break;
	case 2: Cvar_Reset(Cmd_Argv(1)); break; }
}

void Cvar_ResetAll_f()
{ for(cvar_t *var = cvar_vars; var; var = var->next) Cvar_Reset(var->name); }

void Cvar_ResetCfg_f()
{
	for(cvar_t *var = cvar_vars; var; var = var->next)
		if(var->flags & CVAR_ARCHIVE) Cvar_Reset(var->name);
}

void Cvar_Init()
{
	Cmd_AddCommand("cvarlist", Cvar_List_f);
	Cmd_AddCommand("toggle", Cvar_Toggle_f);
	Cmd_AddCommand("cycle", Cvar_Cycle_f);
	Cmd_AddCommand("inc", Cvar_Inc_f);
	Cmd_AddCommand("reset", Cvar_Reset_f);
	Cmd_AddCommand("resetall", Cvar_ResetAll_f);
	Cmd_AddCommand("resetcfg", Cvar_ResetCfg_f);
}

cvar_t *Cvar_FindVar(const s8 *var_name)
{
	for(cvar_t *var = cvar_vars ; var ; var = var->next)
		if(!strcmp(var_name, var->name))
			return var;
	return NULL;
}

cvar_t *Cvar_FindVarAfter(const s8 *prev_name, u32 with_flags)
{
	cvar_t *var = cvar_vars;
	if(*prev_name) {
		var = Cvar_FindVar(prev_name);
		if(!var) return NULL;
		var = var->next;
	}
	// search for the next cvar matching the needed flags
	while(var) {
		if((var->flags & with_flags) || !with_flags) break;
		var = var->next;
	}
	return var;
}

void Cvar_LockVar(const s8 *var_name)
{
	cvar_t *var = Cvar_FindVar(var_name);
	if(var) var->flags |= CVAR_LOCKED;
}

void Cvar_UnlockVar(const s8 *var_name)
{
	cvar_t *var = Cvar_FindVar(var_name);
	if(var) var->flags &= ~CVAR_LOCKED;
}

void Cvar_UnlockAll()
{
	for(cvar_t *var = cvar_vars ; var ; var = var->next)
		var->flags &= ~CVAR_LOCKED;
}

f32 Cvar_VariableValue(const s8 *var_name)
{
	cvar_t *var;
	var = Cvar_FindVar(var_name);
	if(!var) return 0;
	return Q_atof(var->string);
}

const s8 *Cvar_VariableString(const s8 *var_name)
{
	cvar_t *var = Cvar_FindVar(var_name);
	if(!var) return cvar_null_string;
	return var->string;
}

const s8 *Cvar_CompleteVariable(const s8 *partial)
{
	s32 len = strlen(partial);
	if(!len) return NULL;
	for(cvar_t *cvar = cvar_vars; cvar; cvar = cvar->next) //check functions
		if(!strncmp(partial, cvar->name, len))
			return cvar->name;
	return NULL;
}

void Cvar_SetQuick(cvar_t *var, const s8 *value)
{
	if(var->flags & (CVAR_ROM|CVAR_LOCKED)) return;
	if(!(var->flags & CVAR_REGISTERED)) return;
	if(!var->string) var->string = Z_Strdup(value);
	else {
		if(!strcmp(var->string, value)) return; // no change
		var->flags |= CVAR_CHANGED;
		s32 len = strlen(value);
		if((u64)len != strlen(var->string)) {
			Z_Free((void *)var->string);
			var->string = (s8 *) Z_Malloc(len + 1);
		}
		memcpy((s8 *)var->string, value, len + 1);
	}
	var->value = Q_atof(var->string);
	//johnfitz -- save initial value for "reset" command
	if(!var->default_string) var->default_string = Z_Strdup(var->string);
	//johnfitz -- during initialization, update default too
	else if(!host_initialized) {
		Con_DPrintf("changing default of %s: %s -> %s\n", var->name,
				var->default_string, var->string);
		Z_Free((void *)var->default_string);
		var->default_string = Z_Strdup(var->string);
	}
	if(var->callback) var->callback(var);
}

void Cvar_SetValueQuick(cvar_t *var, const f32 value)
{
	s8 val[32], *ptr = val;
	if(value == (f32)((s32)value))
		q_snprintf(val, sizeof(val), "%i", (s32)value);
	else {
		q_snprintf(val, sizeof(val), "%f", value);
		while(*ptr) ptr++; // kill trailing zeroes
		while(--ptr > val && *ptr == '0' && ptr[-1] != '.') *ptr = '\0';
	}
	Cvar_SetQuick(var, val);
}

void Cvar_Set(const s8 *var_name, const s8 *value)
{
	cvar_t *var = Cvar_FindVar(var_name);
	if(!var) { // there is an error in C code if this happens
		Con_Printf("Cvar_Set: variable %s not found\n", var_name);
		return;
	}
	Cvar_SetQuick(var, value);
}

void Cvar_SetValue(const s8 *var_name, const f32 value)
{
	s8 val[32], *ptr = val;
	if(value == (f32)((s32)value))
		q_snprintf(val, sizeof(val), "%i", (s32)value);
	else {
		q_snprintf(val, sizeof(val), "%f", value);
		while(*ptr) ptr++; // kill trailing zeroes
		while(--ptr > val && *ptr == '0' && ptr[-1] != '.') *ptr = '\0';
	}
	Cvar_Set(var_name, val);
}

void Cvar_SetROM(const s8 *var_name, const s8 *value)
{
	cvar_t *var = Cvar_FindVar(var_name);
	if(var) {
		var->flags &= ~CVAR_ROM;
		Cvar_SetQuick(var, value);
		var->flags |= CVAR_ROM;
	}
}

void Cvar_SetValueROM(const s8 *var_name, const f32 value)
{
	cvar_t *var = Cvar_FindVar(var_name);
	if(var) {
		var->flags &= ~CVAR_ROM;
		Cvar_SetValueQuick(var, value);
		var->flags |= CVAR_ROM;
	}
}

void Cvar_RegisterVariable(cvar_t *variable)
{ // Adds a freestanding variable to the variable list.
	s8 value[512];
	if(Cvar_FindVar(variable->name)) { // first check if already defined
    Con_Printf("Can't register variable %s, already defined\n", variable->name);
		return; }
	if(Cmd_Exists(variable->name)) { // check for overlap with a command
	 Con_Printf("Cvar_RegisterVariable: %s is a command\n", variable->name);
		return; }
	// link the variable in
	//johnfitz -- insert each entry in alphabetical order
	if(cvar_vars == NULL || strcmp(variable->name, cvar_vars->name) < 0) {
		variable->next = cvar_vars; // insert at front
		cvar_vars = variable;
	}
	else { // insert later
		cvar_t *prev = cvar_vars;
		cvar_t *cursor = cvar_vars->next;
		while(cursor && (strcmp(variable->name, cursor->name) > 0)) {
			prev = cursor;
			cursor = cursor->next;
		}
		variable->next = prev->next;
		prev->next = variable;
	} //johnfitz
	variable->flags |= CVAR_REGISTERED;
	// copy the value off, because future sets will Z_Free it
	q_strlcpy(value, variable->string, sizeof(value));
	variable->string = NULL;
	variable->default_string = NULL;
	if(!(variable->flags & CVAR_CALLBACK)) variable->callback = NULL;
	// set it through the function to be consistent
	bool set_rom = (variable->flags & CVAR_ROM);
	variable->flags &= ~CVAR_ROM;
	Cvar_SetQuick(variable, value);
	if(set_rom) variable->flags |= CVAR_ROM;
}


void Cvar_SetCallback(cvar_t *var, cvarcallback_t func)
{ // Set a callback function to the var
	var->callback = func;
	if(func) var->flags |= CVAR_CALLBACK;
	else var->flags &= ~CVAR_CALLBACK;
}


bool Cvar_Command()
{ // Handles variable inspection and changing from the console
	cvar_t *v = Cvar_FindVar(Cmd_Argv(0)); // check variables
	if(!v) return 0;
	if(Cmd_Argc() == 1) { // perform a variable print or set
		Con_Printf("\"%s\" is \"%s\"\n", v->name, v->string);
		return 1;
	}
	Cvar_Set(v->name, Cmd_Argv(1));
	return 1;
}


void Cvar_WriteVariables(FILE *f) // Write lines containing "set variable value"
{ // for all variables with the archive flag set to 1.
	for(cvar_t *var = cvar_vars ; var ; var = var->next)
		if(var->flags & CVAR_ARCHIVE)
			fprintf(f, "%s \"%s\"\n", var->name, var->string);
}

