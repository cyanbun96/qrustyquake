// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.

#include "quakedef.h"

static prstack_t pr_stack[MAX_STACK_DEPTH];
static s32 pr_depth;

const s8 *PR_GlobalString(s32 ofs);
const s8 *PR_GlobalStringNoContents(s32 ofs);

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

static const s8 *const pr_extnames[QCEXT_COUNT] =
{
	"STD_QC",
	#define QCEXTENSION(name) #name,
	QCEXTENSIONS_ALL
	#undef QCEXTENSION
};

s32 PR_FindExtensionByName(const s8 *name)
{
	for(s32 i = 1; i < QCEXT_COUNT; i++)
		if(!strcmp(name, pr_extnames[i]))
			return i;
	return 0;
}

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
	else if((u32)(s->op-OP_STORE_F) < 6){
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
	if(pr_depth == 0){ Con_Printf("<NO STACK>\n"); return; }
	pr_stack[pr_depth].f = pr_xfunction;
	for(s32 i = pr_depth; i >= 0; i--){
		dfunction_t *f = pr_stack[i].f;
		if(!f) Con_Printf("<NO FUNCTION>\n");
		else Con_Printf("%12s : %s\n", PR_GetString(f->s_file),
						PR_GetString(f->s_name));
	}
}

void PR_Profile_f()
{
	if(!sv.active)
		return;
	PR_SwitchQCVM(&sv.qcvm);
	s32 num = 0;
	dfunction_t *best;
	do{
		s32 pmax = 0;
		best = NULL;
		for(s32 i = 0; i < qcvm->progs->numfunctions; i++){
			dfunction_t *f = &qcvm->functions[i];
			if(f->profile > pmax){
				pmax = f->profile;
				best = f;
			}
		}
		if(best){
			if(num < 10)
				Con_Printf("%7i %s\n", best->profile,
						PR_GetString(best->s_name));
			num++;
			best->profile = 0;
		}
	}while(best);
	PR_SwitchQCVM(NULL);
}

void PR_RunError(const s8 *error, ...)
{ // Aborts the currently executing function
	va_list argptr;
	s8 string[1024];
	va_start(argptr, error);
	q_vsnprintf(string, sizeof(string), error, argptr);
	va_end(argptr);
	PR_PrintStatement(qcvm->statements + qcvm->xstatement);
	PR_StackTrace();
	Con_Printf("%s\n", string);
	pr_depth = 0; // dump the stack so host_error can shutdown functions
	Host_Error("Program error");
}


static s32 PR_EnterFunction(dfunction_t *f)
{ // Returns the new program statement counter
	qcvm->stack[qcvm->depth].s = qcvm->xstatement;
	qcvm->stack[qcvm->depth].f = qcvm->xfunction;
	qcvm->depth++;
	if(qcvm->depth >= MAX_STACK_DEPTH)
		PR_RunError("stack overflow");
	s32 c = f->locals; // save off any locals that the new function steps on
	if(qcvm->localstack_used + c > LOCALSTACK_SIZE)
		PR_RunError("PR_ExecuteProgram: locals stack overflow");
	for(s32 i = 0; i < c ; i++)
		qcvm->localstack[qcvm->localstack_used + i] =
			((s32*)qcvm->globals)[f->parm_start + i];
	qcvm->localstack_used += c;
	s32 o = f->parm_start; // copy parameters
	for(s32 i = 0; i < f->numparms; i++){
		for(s32 j = 0; j < f->parm_size[i]; j++){
			((s32 *)qcvm->globals)[o] =
				((s32 *)qcvm->globals)[OFS_PARM0 + i*3 + j];
			o++;
		}
	}
	qcvm->xfunction = f;
	return f->first_statement - 1;  // offset the s++
}

static s32 PR_LeaveFunction()
{
	if(qcvm->depth <= 0)
		Host_Error("prog stack underflow");
	// Restore locals from the stack
	s32 c = qcvm->xfunction->locals;
	qcvm->localstack_used -= c;
	if(qcvm->localstack_used < 0)
		PR_RunError("PR_ExecuteProgram: locals stack underflow");
	for(s32 i = 0; i < c; i++)
		((s32 *)qcvm->globals)[qcvm->xfunction->parm_start + i] =
			qcvm->localstack[qcvm->localstack_used + i];
	// up stack
	qcvm->depth--;
	qcvm->xfunction = qcvm->stack[qcvm->depth].f;
	return qcvm->stack[qcvm->depth].s;
}

static void PR_CheckBuiltinExtension(dfunction_t *func)
{
	u32 builtin = -func->first_statement;
	u32 extnum = qcvm->builtin_ext[builtin];
	u32 checked, advertised;
	if(!extnum)
		return;
	checked = GetBit(qcvm->checked_ext, extnum);
	advertised = GetBit(qcvm->advertised_ext, extnum);
	if(checked && advertised)
		return;
	if(GetBit(qcvm->warned_builtin[checked], builtin))
		return;
	SetBit(qcvm->warned_builtin[checked], builtin);
	Con_DPrintf(checked ?
		"[%s] \"%s\" ignored when calling %s(%s: %s)\n" :
		"[%s] check \"%s\" before calling %s(%s: %s)\n",
		(qcvm == &cl.qcvm) ? "CL" : "SV",
		pr_extnames[extnum], PR_GetString(func->s_name),
		PR_GetString(qcvm->xfunction->s_file),
		PR_GetString(qcvm->xfunction->s_name));
}

#define OPA ((eval_t *)&qcvm->globals[(u16)st->a])
#define OPB ((eval_t *)&qcvm->globals[(u16)st->b])
#define OPC ((eval_t *)&qcvm->globals[(u16)st->c])

void PR_ExecuteProgram(func_t fnum)
{
	if(!fnum || fnum >= qcvm->progs->numfunctions){
		if(pr_global_struct->self)
			ED_Print(PROG_TO_EDICT(pr_global_struct->self));
		Host_Error("PR_ExecuteProgram: NULL function");
	}
	eval_t *ptr;
	dfunction_t *f, *newf;
	f = &qcvm->functions[fnum];
	qcvm->trace = false;
	s32 exitdepth = qcvm->depth;//make a stack frame
	dstatement_t *st = &qcvm->statements[PR_EnterFunction(f)];
	s32 startprofile = 0;
	s32 profile = 0;
	while(1){
		st++;//next statement
		if(++profile > 0x1000000){//was 100000
			qcvm->xstatement = st - qcvm->statements;
			PR_RunError("runaway loop error");
		}
		if(qcvm->trace)
			PR_PrintStatement(st);
		switch(st->op){
		case OP_ADD_F:
			OPC->_float = OPA->_float + OPB->_float;
			break;
		case OP_ADD_V:
			OPC->vector[0] = OPA->vector[0] + OPB->vector[0];
			OPC->vector[1] = OPA->vector[1] + OPB->vector[1];
			OPC->vector[2] = OPA->vector[2] + OPB->vector[2];
			break;
		case OP_SUB_F:
			OPC->_float = OPA->_float - OPB->_float;
			break;
		case OP_SUB_V:
			OPC->vector[0] = OPA->vector[0] - OPB->vector[0];
			OPC->vector[1] = OPA->vector[1] - OPB->vector[1];
			OPC->vector[2] = OPA->vector[2] - OPB->vector[2];
			break;
		case OP_MUL_F:
			OPC->_float = OPA->_float * OPB->_float;
			break;
		case OP_MUL_V:
			OPC->_float = OPA->vector[0] * OPB->vector[0] +
				OPA->vector[1] * OPB->vector[1] +
				OPA->vector[2] * OPB->vector[2];
			break;
		case OP_MUL_FV:
			OPC->vector[0] = OPA->_float * OPB->vector[0];
			OPC->vector[1] = OPA->_float * OPB->vector[1];
			OPC->vector[2] = OPA->_float * OPB->vector[2];
			break;
		case OP_MUL_VF:
			OPC->vector[0] = OPB->_float * OPA->vector[0];
			OPC->vector[1] = OPB->_float * OPA->vector[1];
			OPC->vector[2] = OPB->_float * OPA->vector[2];
			break;
		case OP_DIV_F:
			OPC->_float = OPA->_float / OPB->_float;
			break;
		case OP_BITAND:
			OPC->_float = (s32)OPA->_float & (s32)OPB->_float;
			break;
		case OP_BITOR:
			OPC->_float = (s32)OPA->_float | (s32)OPB->_float;
			break;
		case OP_GE:
			OPC->_float = OPA->_float >= OPB->_float;
			break;
		case OP_LE:
			OPC->_float = OPA->_float <= OPB->_float;
			break;
		case OP_GT:
			OPC->_float = OPA->_float > OPB->_float;
			break;
		case OP_LT:
			OPC->_float = OPA->_float < OPB->_float;
			break;
		case OP_AND:
			OPC->_float = OPA->_float && OPB->_float;
			break;
		case OP_OR:
			OPC->_float = OPA->_float || OPB->_float;
			break;
		case OP_NOT_F:
			OPC->_float = !OPA->_float;
			break;
		case OP_NOT_V:
			OPC->_float = !OPA->vector[0] && !OPA->vector[1]
							&& !OPA->vector[2];
			break;
		case OP_NOT_S:
			OPC->_float = !OPA->string||!*PR_GetString(OPA->string);
			break;
		case OP_NOT_FNC:
			OPC->_float = !OPA->function;
			break;
		case OP_NOT_ENT:
			OPC->_float = (PROG_TO_EDICT(OPA->edict)==qcvm->edicts);
			break;
		case OP_EQ_F:
			OPC->_float = OPA->_float == OPB->_float;
			break;
		case OP_EQ_V:
			OPC->_float = (OPA->vector[0] == OPB->vector[0]) &&
				(OPA->vector[1] == OPB->vector[1]) &&
				(OPA->vector[2] == OPB->vector[2]);
			break;
		case OP_EQ_S:
			OPC->_float = !strcmp(PR_GetString(OPA->string),
						PR_GetString(OPB->string));
			break;
		case OP_EQ_E:
			OPC->_float = OPA->_int == OPB->_int;
			break;
		case OP_EQ_FNC:
			OPC->_float = OPA->function == OPB->function;
			break;

		case OP_NE_F:
			OPC->_float = OPA->_float != OPB->_float;
			break;
		case OP_NE_V:
			OPC->_float = (OPA->vector[0] != OPB->vector[0]) ||
				(OPA->vector[1] != OPB->vector[1]) ||
				(OPA->vector[2] != OPB->vector[2]);
			break;
		case OP_NE_S:
			OPC->_float = strcmp(PR_GetString(OPA->string),
						PR_GetString(OPB->string));
			break;
		case OP_NE_E:
			OPC->_float = OPA->_int != OPB->_int;
			break;
		case OP_NE_FNC:
			OPC->_float = OPA->function != OPB->function;
			break;
		case OP_STORE_F:
		case OP_STORE_ENT:
		case OP_STORE_FLD: // integers
		case OP_STORE_S:
		case OP_STORE_FNC: // pointers
			OPB->_int = OPA->_int;
			break;
		case OP_STORE_V:
			OPB->vector[0] = OPA->vector[0];
			OPB->vector[1] = OPA->vector[1];
			OPB->vector[2] = OPA->vector[2];
			break;
		case OP_STOREP_F:
		case OP_STOREP_ENT:
		case OP_STOREP_FLD: // integers
		case OP_STOREP_S:
		case OP_STOREP_FNC: // pointers
			ptr = (eval_t *)((u8 *)qcvm->edicts + OPB->_int);
			ptr->_int = OPA->_int;
			break;
		case OP_STOREP_V:
			ptr = (eval_t *)((u8 *)qcvm->edicts + OPB->_int);
			ptr->vector[0] = OPA->vector[0];
			ptr->vector[1] = OPA->vector[1];
			ptr->vector[2] = OPA->vector[2];
			break;
		case OP_ADDRESS:
			edict_t *ed = PROG_TO_EDICT(OPA->edict);
			if(ed==(edict_t *)qcvm->edicts && sv.state==ss_active){
				qcvm->xstatement = st - qcvm->statements;
				PR_RunError("assignment to world entity");
			}
			OPC->_int = (u8 *)((s32 *)&ed->v + OPB->_int)
					- (u8 *)qcvm->edicts;
			break;
		case OP_LOAD_F:
		case OP_LOAD_FLD:
		case OP_LOAD_ENT:
		case OP_LOAD_S:
		case OP_LOAD_FNC:
			ed = PROG_TO_EDICT(OPA->edict);
			OPC->_int = ((eval_t*)((s32*)&ed->v + OPB->_int))->_int;
			break;
		case OP_LOAD_V:
			ed = PROG_TO_EDICT(OPA->edict);
			ptr = (eval_t *)((s32 *)&ed->v + OPB->_int);
			OPC->vector[0] = ptr->vector[0];
			OPC->vector[1] = ptr->vector[1];
			OPC->vector[2] = ptr->vector[2];
			break;
		case OP_IFNOT:
			if(!OPA->_int)
				st += st->b - 1; // -1 to offset the st++
			break;
		case OP_IF:
			if(OPA->_int)
				st += st->b - 1; // -1 to offset the st++
			break;
		case OP_GOTO:
			st += st->a - 1; // -1 to offset the st++
			break;
		case OP_CALL0:
		case OP_CALL1:
		case OP_CALL2:
		case OP_CALL3:
		case OP_CALL4:
		case OP_CALL5:
		case OP_CALL6:
		case OP_CALL7:
		case OP_CALL8:
			qcvm->xfunction->profile += profile - startprofile;
			startprofile = profile;
			qcvm->xstatement = st - qcvm->statements;
			qcvm->argc = st->op - OP_CALL0;
			if(!OPA->function)
				PR_RunError("NULL function");
			newf = &qcvm->functions[OPA->function];
			if(newf->first_statement < 0){ // Built-in function
				int i = -newf->first_statement;
				if(i >= qcvm->numbuiltins)
				    PR_RunError("Bad builtin call number %d",i);
				PR_CheckBuiltinExtension(newf);
				qcvm->builtins[i]();
				break;
			}
			// Normal function
			st = &qcvm->statements[PR_EnterFunction(newf)];
			break;
		case OP_DONE:
		case OP_RETURN:
			qcvm->xfunction->profile += profile - startprofile;
			startprofile = profile;
			qcvm->xstatement = st - qcvm->statements;
			qcvm->globals[OFS_RETURN] = qcvm->globals[(u16)st->a];
			qcvm->globals[OFS_RETURN+1]=qcvm->globals[(u16)st->a+1];
			qcvm->globals[OFS_RETURN+2]=qcvm->globals[(u16)st->a+2];
			st = &qcvm->statements[PR_LeaveFunction()];
			if(qcvm->depth == exitdepth) // Done
				return;
			break;
		case OP_STATE:
			ed = PROG_TO_EDICT(pr_global_struct->self);
			ed->v.nextthink = pr_global_struct->time + 0.1;
			ed->v.frame = OPA->_float;
			ed->v.think = OPB->function;
			break;
		default:
			qcvm->xstatement = st - qcvm->statements;
			PR_RunError("Bad opcode %i", st->op);
		}
	}//end of while(1) loop
}

#undef OPA
#undef OPB
#undef OPC
