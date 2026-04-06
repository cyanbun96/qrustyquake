// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.

#include "quakedef.h"

static prstack_t pr_stack[MAX_STACK_DEPTH];
static s32 pr_depth;
static s32 localstack[LOCALSTACK_SIZE];
static s32 localstack_used;
static s32 pr_xstatement;

static const s8 *pr_opnames[] = {
	"DONE",
	"MUL_F", "MUL_V", "MUL_FV", "MUL_VF", "DIV",
	"ADD_F", "ADD_V", "SUB_F", "SUB_V",
	"EQ_F", "EQ_V", "EQ_S", "EQ_E", "EQ_FNC",
	"NE_F", "NE_V", "NE_S", "NE_E", "NE_FNC",
	"LE", "GE", "LT", "GT",
	"INDIRECT", "INDIRECT", "INDIRECT", "INDIRECT", "INDIRECT", "INDIRECT",
	"ADDRESS",
	"STORE_F", "STORE_V", "STORE_S", "STORE_ENT", "STORE_FLD", "STORE_FNC",
	"STOREP_F","STOREP_V","STOREP_S","STOREP_ENT","STOREP_FLD","STOREP_FNC",
	"RETURN",
	"NOT_F", "NOT_V", "NOT_S", "NOT_ENT", "NOT_FNC",
	"IF", "IFNOT",
	"CALL0", "CALL1", "CALL2", "CALL3", "CALL4",
	"CALL5", "CALL6", "CALL7", "CALL8",
	"STATE", "GOTO",
	"AND", "OR", "BITAND", "BITOR"
};

const s8 *PR_GlobalString(s32 ofs);
const s8 *PR_GlobalStringNoContents(s32 ofs);

static void PR_PrintStatement(dstatement_t *s)
{
	if((u32)s->op < Q_COUNTOF(pr_opnames)){
		Con_Printf("%s ", pr_opnames[s->op]);
		s32 i = strlen(pr_opnames[s->op]);
		for( ; i < 10; i++)
			Con_Printf(" ");
	}
	if(s->op == OP_IF || s->op == OP_IFNOT)
		Con_Printf("%sbranch %i", PR_GlobalString(s->a), s->b);
	else if(s->op == OP_GOTO) Con_Printf("branch %i", s->a);
	else if((u32)(s->op-OP_STORE_F) < 6) {
		Con_Printf("%s", PR_GlobalString(s->a));
		Con_Printf("%s", PR_GlobalStringNoContents(s->b));
	} else {
		if(s->a) Con_Printf("%s", PR_GlobalString(s->a));
		if(s->b) Con_Printf("%s", PR_GlobalString(s->b));
		if(s->c) Con_Printf("%s", PR_GlobalStringNoContents(s->c));
	}
	Con_Printf("\n");
}

static void PR_StackTrace()
{
	if(pr_depth == 0) { Con_Printf("<NO STACK>\n"); return; }
	pr_stack[pr_depth].f = pr_xfunction;
	for(s32 i = pr_depth; i >= 0; i--) {
		dfunction_t *f = pr_stack[i].f;
		if(!f) Con_Printf("<NO FUNCTION>\n");
		else Con_Printf("%12s : %s\n", PR_GetString(f->s_file),
						PR_GetString(f->s_name));
	}
}

void PR_Profile_f()
{
	if(!sv.active) return;
	s32 num = 0;
	dfunction_t *best = NULL;
	do {
		s32 pmax = 0;
		for(s32 i = 0; i < progs->numfunctions; i++) {
			dfunction_t *f = &pr_functions[i];
			if(f->profile > pmax) {
				pmax = f->profile;
				best = f;
			}
		} 
		if(best) {
			if(num < 10)
				Con_Printf("%7i %s\n", best->profile,
						PR_GetString(best->s_name));
			num++;
			best->profile = 0;
		}
	} while(best);
}


void PR_RunError(const s8 *error, ...)
{ // Aborts the currently executing function
	va_list argptr;
	s8 string[1024];
	va_start(argptr, error);
	q_vsnprintf(string, sizeof(string), error, argptr);
	va_end(argptr);
	PR_PrintStatement(pr_statements + pr_xstatement);
	PR_StackTrace();
	Con_Printf("%s\n", string);
	pr_depth = 0; // dump the stack so host_error can shutdown functions
	Host_Error("Program error");
}


static s32 PR_EnterFunction(dfunction_t *f)
{ // Returns the new program statement counter
	pr_stack[pr_depth].s = pr_xstatement;
	pr_stack[pr_depth].f = pr_xfunction;
	pr_depth++;
	if(pr_depth >= MAX_STACK_DEPTH)
		PR_RunError("stack overflow");
	s32 c = f->locals; // save off any locals that the new function steps on
	if(localstack_used + c > LOCALSTACK_SIZE)
		PR_RunError("PR_ExecuteProgram: locals stack overflow");
	for(s32 i = 0; i < c; i++)
	      localstack[localstack_used+i]=((s32*)pr_globals)[f->parm_start+i];
	localstack_used += c;
	s32 o = f->parm_start; // copy parameters
	for(s32 i = 0; i < f->numparms; i++) {
	for(s32 j = 0; j < f->parm_size[i]; j++) {
		((s32*)pr_globals)[o]=((s32*)pr_globals)[OFS_PARM0+i*3+j];
		o++;
	}}
	pr_xfunction = f;
	return f->first_statement - 1; // offset the s++
}

static s32 PR_LeaveFunction()
{
	if(pr_depth <= 0) Host_Error("prog stack underflow");
	s32 c = pr_xfunction->locals; // Restore locals from the stack
	localstack_used -= c;
	if(localstack_used < 0)
		PR_RunError("PR_ExecuteProgram: locals stack underflow");
	for(s32 i = 0; i < c; i++)
		((s32 *)pr_globals)[pr_xfunction->parm_start + i] =
			localstack[localstack_used + i];
	pr_depth--; // up stack
	pr_xfunction = pr_stack[pr_depth].f;
	return pr_stack[pr_depth].s;
}

#define OPA ((eval_t *)&pr_globals[(u16)st->a])
#define OPB ((eval_t *)&pr_globals[(u16)st->b])
#define OPC ((eval_t *)&pr_globals[(u16)st->c])
#define VEC vector
#define FP _float
void PR_ExecuteProgram(func_t fnum)
{ // The interpretation main loop
	eval_t *ptr;
	dfunction_t *newf;
	edict_t *ed;
	if(!fnum || fnum >= progs->numfunctions){
		if(pr_global_struct->self)
			ED_Print(PROG_TO_EDICT(pr_global_struct->self));
		Host_Error("PR_ExecuteProgram: NULL function");
	}
	dfunction_t *f = &pr_functions[fnum];
	pr_trace = 0;
	s32 exitdepth = pr_depth; // make a stack frame
	dstatement_t *st = &pr_statements[PR_EnterFunction(f)];
	s32 startprofile = 0, profile = 0;
	while(1) {
		st++; // next statement
		if(++profile > 0x1000000){ // was 100000
			pr_xstatement = st - pr_statements;
			PR_RunError("runaway loop error");
		}
		if(pr_trace) PR_PrintStatement(st);
		switch(st->op) {
case OP_ADD_F:   OPC->FP = OPA->FP + OPB->FP; break;
case OP_SUB_F:   OPC->FP = OPA->FP - OPB->FP; break;
case OP_MUL_F:   OPC->FP = OPA->FP * OPB->FP; break;
case OP_DIV_F:   OPC->FP = OPA->FP / OPB->FP; break;
case OP_BITAND:  OPC->FP = (s32)OPA->FP & (s32)OPB->FP; break;
case OP_BITOR:   OPC->FP = (s32)OPA->FP | (s32)OPB->FP; break;
case OP_GE:      OPC->FP = OPA->FP >= OPB->FP; break;
case OP_LE:      OPC->FP = OPA->FP <= OPB->FP; break;
case OP_GT:      OPC->FP = OPA->FP > OPB->FP; break;
case OP_LT:      OPC->FP = OPA->FP < OPB->FP; break;
case OP_AND:     OPC->FP = OPA->FP && OPB->FP; break;
case OP_OR:      OPC->FP = OPA->FP || OPB->FP; break;
case OP_NOT_F:   OPC->FP = !OPA->FP; break;
case OP_NOT_V:   OPC->FP = !OPA->VEC[0]&&!OPA->VEC[1]&&!OPA->VEC[2]; break;
case OP_NOT_S:   OPC->FP = !OPA->string || !*PR_GetString(OPA->string); break;
case OP_NOT_FNC: OPC->FP = !OPA->function; break;
case OP_NOT_ENT: OPC->FP = (PROG_TO_EDICT(OPA->edict) == sv.edicts); break;
case OP_EQ_F:    OPC->FP = OPA->FP == OPB->FP; break;
case OP_EQ_E:    OPC->FP = OPA->_int == OPB->_int; break;
case OP_EQ_FNC:  OPC->FP = OPA->function == OPB->function; break;
case OP_NE_F:    OPC->FP = OPA->FP != OPB->FP; break;
case OP_NE_E:    OPC->FP = OPA->_int != OPB->_int; break;
case OP_NE_FNC:  OPC->FP = OPA->function != OPB->function; break;
case OP_ADD_V:
	OPC->VEC[0] = OPA->VEC[0] + OPB->VEC[0];
	OPC->VEC[1] = OPA->VEC[1] + OPB->VEC[1];
	OPC->VEC[2] = OPA->VEC[2] + OPB->VEC[2];
	break;
case OP_SUB_V:
	OPC->VEC[0] = OPA->VEC[0] - OPB->VEC[0];
	OPC->VEC[1] = OPA->VEC[1] - OPB->VEC[1];
	OPC->VEC[2] = OPA->VEC[2] - OPB->VEC[2];
	break;
case OP_MUL_V:
	OPC->FP = OPA->VEC[0] * OPB->VEC[0] +
		OPA->VEC[1] * OPB->VEC[1] +
		OPA->VEC[2] * OPB->VEC[2];
	break;
case OP_MUL_FV:
	OPC->VEC[0] = OPA->FP * OPB->VEC[0];
	OPC->VEC[1] = OPA->FP * OPB->VEC[1];
	OPC->VEC[2] = OPA->FP * OPB->VEC[2];
	break;
case OP_MUL_VF:
	OPC->VEC[0] = OPB->FP * OPA->VEC[0];
	OPC->VEC[1] = OPB->FP * OPA->VEC[1];
	OPC->VEC[2] = OPB->FP * OPA->VEC[2];
	break;
case OP_EQ_V:
	OPC->FP = (OPA->VEC[0] == OPB->VEC[0]) &&
		(OPA->VEC[1] == OPB->VEC[1]) &&
		(OPA->VEC[2] == OPB->VEC[2]);
	break;
case OP_EQ_S: OPC->FP = !strcmp(PR_GetString(OPA->string),
			      PR_GetString(OPB->string));
	break;
case OP_NE_V:
	OPC->FP = (OPA->VEC[0] != OPB->VEC[0]) ||
		(OPA->VEC[1] != OPB->VEC[1]) ||
		(OPA->VEC[2] != OPB->VEC[2]);
	break;
case OP_NE_S:
	OPC->FP=strcmp(PR_GetString(OPA->string),PR_GetString(OPB->string));
	break;
case OP_STORE_F: case OP_STORE_ENT: case OP_STORE_FLD:
case OP_STORE_S: case OP_STORE_FNC:
	OPB->_int = OPA->_int; break;
case OP_STORE_V:
	OPB->VEC[0] = OPA->VEC[0];
	OPB->VEC[1] = OPA->VEC[1];
	OPB->VEC[2] = OPA->VEC[2];
	break;
case OP_STOREP_F: case OP_STOREP_ENT: case OP_STOREP_FLD:
case OP_STOREP_S: case OP_STOREP_FNC:
	ptr = (eval_t *)((u8 *)sv.edicts + OPB->_int);
	ptr->_int = OPA->_int;
	break;
case OP_STOREP_V:
	ptr = (eval_t *)((u8 *)sv.edicts + OPB->_int);
	ptr->VEC[0] = OPA->VEC[0];
	ptr->VEC[1] = OPA->VEC[1];
	ptr->VEC[2] = OPA->VEC[2];
	break;
case OP_ADDRESS:
	ed = PROG_TO_EDICT(OPA->edict);
	if(ed == (edict_t *)sv.edicts && sv.state == ss_active) {
		pr_xstatement = st - pr_statements;
		PR_RunError("assignment to world entity");
	}
	OPC->_int = (u8 *)((s32 *)&ed->v + OPB->_int) - (u8 *)sv.edicts;
	break;
case OP_LOAD_F: case OP_LOAD_FLD: case OP_LOAD_ENT:
case OP_LOAD_S: case OP_LOAD_FNC:
	ed = PROG_TO_EDICT(OPA->edict);
	OPC->_int = ((eval_t *)((s32 *)&ed->v + OPB->_int))->_int;
	break;
case OP_LOAD_V:
	ed = PROG_TO_EDICT(OPA->edict);
	ptr = (eval_t *)((s32 *)&ed->v + OPB->_int);
	OPC->VEC[0] = ptr->VEC[0];
	OPC->VEC[1] = ptr->VEC[1];
	OPC->VEC[2] = ptr->VEC[2];
	break;
case OP_IFNOT: if(!OPA->_int) st += st->b - 1; break; // -1 to offset the st++
case OP_IF: if(OPA->_int) st += st->b - 1; break; // -1 to offset the st++
case OP_GOTO: st += st->a - 1; break; // -1 to offset the st++
case OP_CALL0: case OP_CALL1: case OP_CALL2: case OP_CALL3: case OP_CALL4:
case OP_CALL5: case OP_CALL6: case OP_CALL7: case OP_CALL8:
	pr_xfunction->profile += profile - startprofile;
	startprofile = profile;
	pr_xstatement = st - pr_statements;
	pr_argc = st->op - OP_CALL0;
	if(!OPA->function) PR_RunError("NULL function");
	newf = &pr_functions[OPA->function];
	if(newf->first_statement < 0) { // Built-in function
		s32 i = -newf->first_statement;
		if(i >= pr_numbuiltins)
			PR_RunError("Bad builtin call number %d", i);
		pr_builtins[i]();
		break;
	}
	st = &pr_statements[PR_EnterFunction(newf)]; // Normal function
	break;
case OP_DONE: case OP_RETURN:
	pr_xfunction->profile += profile - startprofile;
	startprofile = profile;
	pr_xstatement = st - pr_statements;
	pr_globals[OFS_RETURN] = pr_globals[(u16)st->a];
	pr_globals[OFS_RETURN + 1] = pr_globals[(u16)st->a + 1];
	pr_globals[OFS_RETURN + 2] = pr_globals[(u16)st->a + 2];
	st = &pr_statements[PR_LeaveFunction()];
	if(pr_depth == exitdepth) return; // Done
	break;
case OP_STATE:
	ed = PROG_TO_EDICT(pr_global_struct->self);
	ed->v.nextthink = pr_global_struct->time + 0.1;
	ed->v.frame = OPA->FP;
	ed->v.think = OPB->function;
	break;
default:
	pr_xstatement = st - pr_statements;
	PR_RunError("Bad opcode %i", st->op);
}} // end of while(1) loop
}
#undef OPA
#undef OPB
#undef OPC
#undef VEC
#undef FP
