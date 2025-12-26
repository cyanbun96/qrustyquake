// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// cl_parse.c -- parse a message received from the server
#include "quakedef.h"

static entity_t cl_static_entities[MAX_STATIC_ENTITIES];
static u8 net_olddata[NET_MAXMESSAGE];
static s8 model_precache[MAX_MODELS][MAX_QPATH];
static s8 sound_precache[MAX_SOUNDS][MAX_QPATH];

static s8 *svc_strings[] = {
	"svc_bad",
	"svc_nop",
	"svc_disconnect",
	"svc_updatestat",
	"svc_version", // [s64] server version
	"svc_setview", // [s16] entity number
	"svc_sound", // <see code>
	"svc_time", // [f32] server time
	"svc_print", // [string] null terminated string
	"svc_stufftext", // [string] stuffed into client's console buffer
			 // the string should be \n terminated
	"svc_setangle", // [vec3] set the view angle to this absolute value
	"svc_serverinfo", // [s64] version
			  // [string] signon string
			  // [string]..[0]model cache [string]...[0]sounds cache
			  // [string]..[0]item cache
	"svc_lightstyle", // [u8] [string]
	"svc_updatename", // [u8] [string]
	"svc_updatefrags", // [u8] [s16]
	"svc_clientdata", // <shortbits + data>
	"svc_stopsound", // <see code>
	"svc_updatecolors", // [u8] [u8]
	"svc_particle", // [vec3] <variable>
	"svc_damage", // [u8] impact [u8] blood [vec3] from
	"svc_spawnstatic",
	"OBSOLETE svc_spawnbinary",
	"svc_spawnbaseline",
	"svc_temp_entity", // <variable>
	"svc_setpause",
	"svc_signonnum",
	"svc_centerprint",
	"svc_killedmonster",
	"svc_foundsecret",
	"svc_spawnstaticsound",
	"svc_intermission",
	"svc_finale", // [string] music [string] text
	"svc_cdtrack", // [u8] track [u8] looptrack
	"svc_sellscreen",
	"svc_cutscene"
};

void V_RestoreAngles()
{
	entity_t *ent = &cl_entities[cl.viewentity];
	VectorCopy(ent->msg_angles[0], ent->angles);
}

void R_CheckEfrags()
{
	if(cls.signon < 2) return; //don't spam when still parsing signon
	if(cl.num_efrags > 640) // packet full of static ents
		Con_DPrintf("%i efrags exceeds standard limit of 640.\n",
				cl.num_efrags);
}

void Fog_ParseServerMessage()
{
	f32 density, red, green, blue, time;
	density = MSG_ReadByte() / 255.0;
	red = MSG_ReadByte() / 255.0;
	green = MSG_ReadByte() / 255.0;
	blue = MSG_ReadByte() / 255.0;
	time = MSG_ReadShort() / 100.0;
	if(time < 0.0f) time = 0.0f;
	Fog_Update(density, red, green, blue);
}

entity_t *CL_EntityNum(s32 num)
{
	if(num >= cl.num_entities) {
		if(num >= MAX_EDICTS)
			Host_Error("CL_EntityNum: %i is an invalid number",num);
		while(cl.num_entities<=num) {
			cl_entities[cl.num_entities].colormap = CURWORLDCMAP;
			cl_entities[cl.num_entities].baseline.scale =
				ENTSCALE_DEFAULT;
			cl.num_entities++;
		}
	}
	return &cl_entities[num];
}

void CL_ParseStartSoundPacket()
{
	s32 field_mask = MSG_ReadByte();
	s32 vol;
	if(field_mask & SND_VOLUME) vol = MSG_ReadByte();
	else vol = DEFAULT_SOUND_PACKET_VOLUME;
	f32 attn;
	if(field_mask & SND_ATTENUATION) attn = MSG_ReadByte() / 64.0;
	else attn = DEFAULT_SOUND_PACKET_ATTENUATION;
	s32 ch, ent;
	if(field_mask & SND_LARGEENTITY) { //johnfitz -- PROTOCOL_FITZQUAKE
		ent =(u16) MSG_ReadShort ();
		ch = MSG_ReadByte();
	} else {
		ch =(u16) MSG_ReadShort ();
		ent = ch >> 3;
		ch &= 7;
	}
	s32 snum;
	if(field_mask & SND_LARGESOUND) snum = (u16) MSG_ReadShort ();
	else snum = MSG_ReadByte();
	if(snum >= MAX_SOUNDS) //johnfitz -- check soundnum
		Host_Error("CL_ParseStartSoundPacket: %i > MAX_SOUNDS", snum);
	vec3_t pos;
	for(s32 i = 0; i < 3; i++)
		pos[i] = MSG_ReadCoord(cl.protocolflags);
	S_StartSound(ent, ch, cl.sound_precache[snum], pos, vol/255.0, attn);
}

void CL_ParseLocalSound()
{
	s32 field_mask = MSG_ReadByte();
	s32 sn = (field_mask&SND_LARGESOUND) ? MSG_ReadShort() : MSG_ReadByte();
	if(sn >= MAX_SOUNDS)
		Host_Error("CL_ParseLocalSound: %i > MAX_SOUNDS", sn);
	S_LocalSound(cl.sound_precache[sn]->name);
}

void CL_KeepaliveMessage()
{
	static f32 lastmsg;
	if(sv.active) return; // no need if server is local
	if(cls.demoplayback) return; // read from server, should just be nops
	u8 *olddata = net_olddata;
	sizebuf_t old = net_message;
	memcpy(olddata, net_message.data, net_message.cursize);
	s32 ret;
	do {
		ret = CL_GetMessage();
		switch(ret)
		{
		default:Host_Error("CL_KeepaliveMessage: CL_GetMessage failed");
		case 0: break; // nothing waiting
		case 1:Host_Error("CL_KeepaliveMessage: received a message");
			break;
		case 2: if(MSG_ReadByte() != svc_nop)
				Host_Error(
				"CL_KeepaliveMessage: datagram wasn't a nop");
			break;
		}
	} while(ret);
	net_message = old;
	memcpy(net_message.data, olddata, net_message.cursize);
	f32 time = Sys_DoubleTime(); // check time
	if(time - lastmsg < 5) return;
	lastmsg = time;
	Con_Printf("--> client to server keepalive\n"); // write out a nop
	MSG_WriteByte(&cls.message, clc_nop);
	NET_SendMessage(cls.netcon, &cls.message);
	SZ_Clear(&cls.message);
}

void CL_ParseServerInfo()
{
	Con_DPrintf("Serverinfo packet received.\n");
	if(cls.demoplayback) // ericw - loading plaque for map changes in a demo
		SCR_BeginLoadingPlaque(); // it will be hidden in CL_SignonReply
	CL_ClearState(); // wipe the client_state_t struct
	s32 i = MSG_ReadLong(); // parse protocol version number
	// johnfitz -- support multiple protocols
	if(i!=PROTOCOL_NETQUAKE && i!=PROTOCOL_FITZQUAKE && i!=PROTOCOL_RMQ){
		Con_Printf("\n"); // because there's no newline after serverinfo
		Host_Error("Server returned version %i, not %i or %i or %i", i,
			PROTOCOL_NETQUAKE, PROTOCOL_FITZQUAKE, PROTOCOL_RMQ);
	}
	cl.protocol = i; //johnfitz
	if(cl.protocol == PROTOCOL_RMQ) {
		const u32 supportedflags = (PRFL_SHORTANGLE | PRFL_FLOATANGLE |
			PRFL_24BITCOORD | PRFL_FLOATCOORD | PRFL_EDICTSCALE |
			PRFL_INT32COORD);
		// mh - read protocol flags from server so that we know what
		cl.protocolflags =(u32) MSG_ReadLong (); // features to expect
		if(0 != (cl.protocolflags & (~supportedflags))) Con_Printf(
		   "PROTOCOL_RMQ protocolflags %i contains unsupported flags\n",
		   cl.protocolflags);
	}
	else cl.protocolflags = 0;
	cl.maxclients = MSG_ReadByte(); // parse maxclients
	if(cl.maxclients < 1 || cl.maxclients > MAX_SCOREBOARD)
		Host_Error("Bad maxclients (%u) from server", cl.maxclients);
	cl.scores =(scoreboard_t *) Hunk_AllocName
		(cl.maxclients*sizeof(*cl.scores), "scores");
	cl.gametype = MSG_ReadByte(); // parse gametype
	const s8 *str = MSG_ReadString(); // parse signon message
	q_strlcpy(cl.levelname, str, sizeof(cl.levelname));
	// seperate the printfs so the server message can have a color
	Con_Printf("%c%s\n", 2, str);
	Con_Printf("Using protocol %i\n", i); //johnfitz
	// first we go through and touch all of the precache data that still
	// happens to be in the cache, so precaching something else doesn't
	// needlessly purge it
	memset(cl.model_precache, 0, sizeof(cl.model_precache));
	s32 nummodels = 1;
	for(;; nummodels++) { // precache models
		str = MSG_ReadString();
		if(!str[0]) break;
		if(nummodels == MAX_MODELS)
			Host_Error("Server sent too many model precaches");
		q_strlcpy(model_precache[nummodels], str, MAX_QPATH);
		Mod_TouchModel(str);
	}
	if(nummodels >= 256) //johnfitz -- check for excessive models
	    Con_DPrintf("%i models exceeds standard limit of 256 (max = %d).\n",
			nummodels, MAX_MODELS);
	memset(cl.sound_precache, 0, sizeof(cl.sound_precache));//precache sound
	s32 numsounds = 1;
	for(;; numsounds++) {
		str = MSG_ReadString();
		if(!str[0]) break;
		if(numsounds == MAX_SOUNDS)
			Host_Error("Server sent too many sound precaches");
		q_strlcpy(sound_precache[numsounds], str, MAX_QPATH);
		S_TouchSound(str);
	}
	if(numsounds >= 256) //johnfitz -- check for excessive sounds
	    Con_DPrintf("%i sounds exceeds standard limit of 256 (max = %d).\n",
			numsounds, MAX_SOUNDS);
	// now we try to load everything else until a cache allocation fails
	// copy the naked name of the map file to the cl structure -- O.S
	COM_StripExtension(COM_SkipPath(model_precache[1]), cl.mapname, sizeof(cl.mapname));
	for(i = 1; i < nummodels; i++) {
		cl.model_precache[i] = Mod_ForName(model_precache[i], 0);
		if(cl.model_precache[i] == NULL)
			Host_Error("Model %s not found", model_precache[i]);
		CL_KeepaliveMessage();
	}
	for(i = 1; i < numsounds; i++) {
		cl.sound_precache[i] = S_PrecacheSound(sound_precache[i]);
		CL_KeepaliveMessage();
	}
	cl_entities[0].model = cl.worldmodel = cl.model_precache[1];
	R_NewMap();
	Hunk_Check(); // make sure nothing is hurt
	noclip_anglehack = 0; // noclip is turned off at start
}

void CL_ParseBaseline(entity_t *ent, s32 version) //johnfitz -- added argument
{
	s32 bits = (version == 2) ? MSG_ReadByte() : 0; // PROTOCOL_FITZQUAKE
	ent->baseline.modelindex = (bits & B_LARGEMODEL)
		? MSG_ReadShort() : MSG_ReadByte();
	ent->baseline.frame = (bits & B_LARGEFRAME)
		? MSG_ReadShort() : MSG_ReadByte();
	ent->baseline.colormap = MSG_ReadByte();
	ent->baseline.skin = MSG_ReadByte();
	for(s32 i = 0; i < 3; i++) {
		ent->baseline.origin[i] = MSG_ReadCoord(cl.protocolflags);
		ent->baseline.angles[i] = MSG_ReadAngle(cl.protocolflags);
	} //johnfitz -- PROTOCOL_FITZQUAKE
	ent->baseline.alpha = (bits & B_ALPHA)?MSG_ReadByte():ENTALPHA_DEFAULT;
	ent->baseline.scale = (bits & B_SCALE)?MSG_ReadByte():ENTSCALE_DEFAULT;
}


void CL_ParseUpdate(s32 bits) // If an entities model or origin changes from
{ // frame to frame, it must be relinked. Others can change without relinking.
	if(cls.signon == SIGNONS - 1){ // first update is the final signon stage
		cls.signon = SIGNONS;
		CL_SignonReply();
	}
	if(bits & U_MOREBITS) bits |= (MSG_ReadByte()<<8);
	if(cl.protocol == PROTOCOL_FITZQUAKE || cl.protocol == PROTOCOL_RMQ) {
		if(bits & U_EXTEND1) bits |= MSG_ReadByte() << 16;
		if(bits & U_EXTEND2) bits |= MSG_ReadByte() << 24;
	}
	s32 num = bits & U_LONGENTITY ? MSG_ReadShort() : MSG_ReadByte();
	entity_t *ent = CL_EntityNum(num);
	bool forcelink = ent->msgtime != cl.mtime[1];
	ent->msgtime = cl.mtime[0];
	s32 modnum = bits & U_MODEL ? MSG_ReadByte() : ent->baseline.modelindex;
	if(modnum >= MAX_MODELS) Host_Error("CL_ParseModel: bad modnum");
	ent->frame = bits & U_FRAME ? MSG_ReadByte() : ent->baseline.frame;
	s32 i = bits & U_COLORMAP ? MSG_ReadByte() : ent->baseline.colormap;
	if(!i) ent->colormap = CURWORLDCMAP;
	else { if(i > cl.maxclients)
			Sys_Error("i >= cl.maxclients");
		ent->colormap = cl.scores[i-1].translations;
	}
	s32 skin = bits & U_SKIN ? MSG_ReadByte() : ent->baseline.skin;
	if(skin != ent->skinnum) ent->skinnum = skin;
	ent->effects = bits&U_EFFECTS ? MSG_ReadByte() : ent->baseline.effects;
	// shift the known values for interpolation
	VectorCopy(ent->msg_origins[0], ent->msg_origins[1]);
	VectorCopy(ent->msg_angles[0], ent->msg_angles[1]);
	ent->msg_origins[0][0] = bits & U_ORIGIN1 ?
		MSG_ReadCoord(cl.protocolflags) : ent->baseline.origin[0];
	ent->msg_angles[0][0] = bits & U_ANGLE1 ?
		 MSG_ReadAngle(cl.protocolflags) : ent->baseline.angles[0];
	ent->msg_origins[0][1] = bits & U_ORIGIN2 ?
		MSG_ReadCoord(cl.protocolflags) : ent->baseline.origin[1];
	ent->msg_angles[0][1] = bits & U_ANGLE2 ?
		 MSG_ReadAngle(cl.protocolflags) : ent->baseline.angles[1];
	ent->msg_origins[0][2] = bits & U_ORIGIN3 ?
		 MSG_ReadCoord(cl.protocolflags) : ent->baseline.origin[2];
	ent->msg_angles[0][2] = bits & U_ANGLE3 ?
		 MSG_ReadAngle(cl.protocolflags) : ent->baseline.angles[2];
	if(cl.protocol == PROTOCOL_FITZQUAKE || cl.protocol == PROTOCOL_RMQ) {
		ent->alpha = bits&U_ALPHA ? MSG_ReadByte():ent->baseline.alpha;
		ent->scale = bits&U_SCALE ? MSG_ReadByte():ent->baseline.scale;
		if(bits & U_FRAME2)
			ent->frame = (ent->frame & 0x00FF)|(MSG_ReadByte()<<8);
		if(bits & U_MODEL2)
			modnum = (modnum & 0x00FF) | (MSG_ReadByte() << 8);
		if(bits & U_LERPFINISH) {
			ent->lerpfinish=ent->msgtime+(f32)(MSG_ReadByte())/255;
			ent->lerpflags |= LERP_FINISH;
		} else ent->lerpflags &= ~LERP_FINISH;
	}
	else if(cl.protocol == PROTOCOL_NETQUAKE) {
		if(bits & U_TRANS) { //HACK: if this bit is set, assume NEHAHRA
			f32 a, b;
			a = MSG_ReadFloat();
			b = MSG_ReadFloat(); //alpha
			if(a == 2) MSG_ReadFloat();
			ent->alpha = ENTALPHA_ENCODE(b);
		}
		else ent->alpha = ent->baseline.alpha;
		ent->scale = ent->baseline.scale;
	}
	model_t *model = cl.model_precache[modnum];
	if(model != ent->model) { // automatic animation(torches, etc) can be
		ent->model = model; // either all together or randomized
		if(model)ent->syncbase = model->synctype == ST_RAND ?
				(f32)(rand()&0x7fff) / 0x7fff : 0.0;
		else forcelink = 1; // hack to make null model players work
	}
	if(forcelink) { // didn't have an update last message
		VectorCopy(ent->msg_origins[0], ent->msg_origins[1]);
		VectorCopy(ent->msg_origins[0], ent->origin);
		VectorCopy(ent->msg_angles[0], ent->msg_angles[1]);
		VectorCopy(ent->msg_angles[0], ent->angles);
		ent->forcelink = 1;
	}
}

void CL_ParseClientdata()
{
	s32 bits = (u16)MSG_ReadShort ();//johnfitz -- read bits here isntead of
	if(bits & SU_EXTEND1) // in CL_ParseServerMessage()
		bits |= (MSG_ReadByte() << 16);
	if(bits & SU_EXTEND2)
		bits |= (MSG_ReadByte() << 24);
	cl.viewheight = bits & SU_VIEWHEIGHT?MSG_ReadChar():DEFAULT_VIEWHEIGHT;
	cl.idealpitch = bits & SU_IDEALPITCH?MSG_ReadChar():0;
	VectorCopy(cl.mvelocity[0], cl.mvelocity[1]);
	s32 i, j;
	for(i = 0; i < 3; i++) {
		cl.punchangle[i] = bits & (SU_PUNCH1<<i) ? MSG_ReadChar() : 0;
		cl.mvelocity[0][i] = bits&(SU_VELOCITY1<<i)?MSG_ReadChar()*16:0;
	}
	i = MSG_ReadLong();
	if(cl.items != i) { // set flash times
		Sbar_Changed();
		for(j = 0; j < 32; j++)
			if((i & (1<<j)) && !(cl.items & (1<<j)))
				cl.item_gettime[j] = cl.time;
		cl.items = i;
	}
	cl.onground = (bits & SU_ONGROUND) != 0;
	cl.inwater = (bits & SU_INWATER) != 0;
	cl.stats[STAT_WEAPONFRAME] = (bits&SU_WEAPONFRAME) ? MSG_ReadByte() : 0;
	i = (bits & SU_ARMOR) ? MSG_ReadByte() : 0;
	if (cl.stats[STAT_ARMOR] != i) {
		cl.stats[STAT_ARMOR] = i;
		Sbar_Changed();
	}
	i = (bits & SU_WEAPON) ? MSG_ReadByte() : 0;
	if (cl.stats[STAT_WEAPON] != i) {
		cl.stats[STAT_WEAPON] = i;
		Sbar_Changed();
	}
	i = MSG_ReadShort();
	if (cl.stats[STAT_HEALTH] != i) {
		cl.stats[STAT_HEALTH] = i;
		Sbar_Changed();
	}
	i = MSG_ReadByte();
	if (cl.stats[STAT_AMMO] != i) {
		cl.stats[STAT_AMMO] = i;
		Sbar_Changed();
	}
	for (i = 0; i < 4; i++) {
		j = MSG_ReadByte();
		if (cl.stats[STAT_SHELLS + i] != j) {
			cl.stats[STAT_SHELLS + i] = j;
			Sbar_Changed();
		}
	}
	i = MSG_ReadByte();
	if (standard_quake) {
		if (cl.stats[STAT_ACTIVEWEAPON] != i) {
			cl.stats[STAT_ACTIVEWEAPON] = i;
			Sbar_Changed();
		}
	} else {
		if (cl.stats[STAT_ACTIVEWEAPON] != (1 << i)) {
			cl.stats[STAT_ACTIVEWEAPON] = (1 << i);
			Sbar_Changed();
		}
	}
	if(bits & SU_WEAPON2) //johnfitz -- PROTOCOL_FITZQUAKE
		cl.stats[STAT_WEAPON] |=(MSG_ReadByte() << 8);
	if(bits & SU_ARMOR2)
		cl.stats[STAT_ARMOR] |=(MSG_ReadByte() << 8);
	if(bits & SU_AMMO2)
		cl.stats[STAT_AMMO] |=(MSG_ReadByte() << 8);
	if(bits & SU_SHELLS2)
		cl.stats[STAT_SHELLS] |=(MSG_ReadByte() << 8);
	if(bits & SU_NAILS2)
		cl.stats[STAT_NAILS] |=(MSG_ReadByte() << 8);
	if(bits & SU_ROCKETS2)
		cl.stats[STAT_ROCKETS] |=(MSG_ReadByte() << 8);
	if(bits & SU_CELLS2)
		cl.stats[STAT_CELLS] |=(MSG_ReadByte() << 8);
	if(bits & SU_WEAPONFRAME2)
		cl.stats[STAT_WEAPONFRAME] |=(MSG_ReadByte() << 8);
	cl.viewent.alpha = bits&SU_WEAPONALPHA?MSG_ReadByte():ENTALPHA_DEFAULT;
}

void CL_NewTranslation(s32 slot)
{
	if(slot > cl.maxclients)
		Sys_Error("CL_NewTranslation: slot > cl.maxclients");
	u8 *dest = cl.scores[slot].translations;
	u8 *source = CURWORLDCMAP;
	memcpy(dest, CURWORLDCMAP, sizeof(cl.scores[slot].translations));
	s32 top = cl.scores[slot].colors & 0xf0;
	s32 bottom =(cl.scores[slot].colors &15)<<4;
	for(s32 i = 0; i < VID_GRADES; i++, dest += 256, source+=256) {
		if(top < 128) // the artists made some backwards ranges. sigh.
			memcpy(dest + TOP_RANGE, source + top, 16);
		else for(s32 j = 0; j < 16; j++)
				dest[TOP_RANGE+j] = source[top+15-j];
		if(bottom < 128)
			memcpy(dest + BOTTOM_RANGE, source + bottom, 16);
		else for(s32 j = 0; j < 16; j++)
				dest[BOTTOM_RANGE+j] = source[bottom+15-j];
	}
}

void CL_ParseStatic(s32 version) //johnfitz -- added a parameter
{
	s32 i = cl.num_statics;
	if(i >= MAX_STATIC_ENTITIES) Host_Error("Too many static entities");
	entity_t *ent = &cl_static_entities[i];
	cl.num_statics++;
	CL_ParseBaseline(ent, version); //johnfitz -- added second parameter
	ent->model = cl.model_precache[ent->baseline.modelindex];
	ent->frame = ent->baseline.frame;
	ent->colormap = CURWORLDCMAP;
	ent->skinnum = ent->baseline.skin;
	ent->effects = ent->baseline.effects;
	ent->alpha = ent->baseline.alpha; //johnfitz -- alpha
	ent->scale = ent->baseline.scale;
	VectorCopy(ent->baseline.origin, ent->origin);
	VectorCopy(ent->baseline.angles, ent->angles);
	R_AddEfrags(ent);
}

void CL_ParseStaticSound(s32 version) //johnfitz -- added argument
{
	vec3_t org;
	for(s32 i = 0; i < 3; i++)
		org[i] = MSG_ReadCoord(cl.protocolflags);
	s32 sound_num = version == 2 ? MSG_ReadShort() : MSG_ReadByte(); //FITZQ
	s32 vol = MSG_ReadByte();
	s32 atten = MSG_ReadByte();
	S_StaticSound(cl.sound_precache[sound_num], org, vol, atten);
}

void CL_ParseServerMessage()
{ // if recording demos, copy the message out
	if(cl_shownet.value == 1) Con_Printf("%i ",net_message.cursize);
	else if(cl_shownet.value == 2) Con_Printf("------------------\n");
	cl.onground = 0; // unless the server says otherwise
	MSG_BeginReading(); // parse the message
	s32 lastcmd = 0;
	while(1){
	const s8 *str;
	if(msg_badread)
		Host_Error("CL_ParseServerMessage: Bad server message");
	s32 cmd = MSG_ReadByte();
	if(cmd == -1){
		if(cl_shownet.value==2)
			Con_Printf("%3i:%s\n",msg_readcount-1,
					"END OF MESSAGE");
		return; // end of message
	}
	// if the high bit of the command u8 is set, it is a fast update
	if(cmd & U_SIGNAL){ //johnfitz -- was 128, changed for clarity
		if(cl_shownet.value==2)
			Con_Printf("%3i:%s\n",msg_readcount-1,
					"fast update");
		CL_ParseUpdate(cmd&127);
		continue;
	}
	if(cmd < (s32)Q_COUNTOF(svc_strings)){
		if(cl_shownet.value==2)
			Con_Printf("%3i:%s\n",msg_readcount-1,
					svc_strings[cmd]);
	}
	s32 i, j; // keep here for OpenBSD
	switch(cmd) { // other commands
	default:
		Host_Error(
		"Illegible server message %d (previous was %s)",
		cmd, svc_strings[lastcmd]);
		break;
	case svc_nop:
		Con_Printf("svc_nop\n");
		break;
	case svc_time:
		cl.mtime[1] = cl.mtime[0];
		cl.mtime[0] = MSG_ReadFloat();
		break;
	case svc_clientdata:
		CL_ParseClientdata();
		break;
	case svc_version:
		i = MSG_ReadLong();
		if(i!=PROTOCOL_NETQUAKE&&i!=PROTOCOL_FITZQUAKE
				&&i!=PROTOCOL_RMQ) Host_Error(
		   "Server returned version %i, not %i or %i or %i", i,
		   PROTOCOL_NETQUAKE, PROTOCOL_FITZQUAKE, PROTOCOL_RMQ);
		cl.protocol = i;
		break;
	case svc_disconnect:
		Host_EndGame("Server disconnected\n");
		break;
	case svc_print:
		Con_Printf("%s", MSG_ReadString ());
		break;
	case svc_centerprint:
		str = MSG_ReadString();
		SCR_CenterPrint(str);
		break;
	case svc_stufftext:
		Cbuf_AddText(MSG_ReadString ());
		break;
	case svc_damage:
		V_ParseDamage();
		break;
	case svc_serverinfo:
		CL_ParseServerInfo();
		vid.recalc_refdef = 1; // leave intermission full screen
		break;
	case svc_setangle:
		for(i = 0; i < 3; i++)
		       cl.viewangles[i]=MSG_ReadAngle(cl.protocolflags);
		break;
	case svc_setview:
		cl.viewentity = MSG_ReadShort();
		break;
	case svc_lightstyle:
		i = MSG_ReadByte();
		if(i >= MAX_LIGHTSTYLES)
			Sys_Error("svc_lightstyle > MAX_LIGHTSTYLES");
		q_strlcpy(cl_lightstyle[i].map, MSG_ReadString(),
				MAX_STYLESTRING);
		cl_lightstyle[i].length=Q_strlen(cl_lightstyle[i].map);
		if(cl_lightstyle[i].length){ //save extra info
			s32 total = 0;
			cl_lightstyle[i].peak = 'a';
			for(j = 0; j<cl_lightstyle[i].length; j++) {
				total += cl_lightstyle[i].map[j] - 'a';
				cl_lightstyle[i].peak = q_max(
					cl_lightstyle[i].peak,
					cl_lightstyle[i].map[j]);
			}
			cl_lightstyle[i].average =
				total / cl_lightstyle[i].length + 'a';
		} else cl_lightstyle[i].average =
				cl_lightstyle[i].peak = 'm';
		break;
	case svc_sound:
		CL_ParseStartSoundPacket();
		break;
	case svc_stopsound:
		i = MSG_ReadShort();
		S_StopSound(i>>3, i&7);
		break;
	case svc_updatename:
		Sbar_Changed();
		i = MSG_ReadByte();
		if(i >= cl.maxclients) Host_Error(
		      "CL_ParseServerMessage: svc_updatename > MAX_SCOREBOARD");
		q_strlcpy(cl.scores[i].name, MSG_ReadString(),
				MAX_SCOREBOARDNAME);
		break;
	case svc_updatefrags:
		Sbar_Changed();
		i = MSG_ReadByte();
		if(i >= cl.maxclients) Host_Error(
		     "CL_ParseServerMessage: svc_updatefrags > MAX_SCOREBOARD");
		cl.scores[i].frags = MSG_ReadShort();
		break;
	case svc_updatecolors:
		Sbar_Changed();
		i = MSG_ReadByte();
		if(i >= cl.maxclients) Host_Error(
		    "CL_ParseServerMessage: svc_updatecolors > MAX_SCOREBOARD");
		cl.scores[i].colors = MSG_ReadByte();
		CL_NewTranslation(i);
		break;
	case svc_particle:
		R_ParseParticleEffect();
		break;
	case svc_spawnbaseline:
		i = MSG_ReadShort(); // must use CL_EntityNum() to force
		CL_ParseBaseline(CL_EntityNum(i), 1); // cl.num_entities up
		break;
	case svc_spawnstatic:
		CL_ParseStatic(1);
		break;
	case svc_temp_entity:
		CL_ParseTEnt();
		break;
	case svc_setpause:
		cl.paused = MSG_ReadByte();
		if(cl.paused)
			CDAudio_Pause();
		else
			CDAudio_Resume();
		break;
	case svc_signonnum:
		i = MSG_ReadByte();
		if(i <= cls.signon)
		     Host_Error("Received signon %i when at %i", i, cls.signon);
		cls.signon = i;//if signonnum==2, signon packet has been fully
		if(i==2){//parsed, so check for excessive static ents and efrags
			if(cl.num_statics > 128) Con_DPrintf(
		"%i static entities exceeds standard limit of 128 (max = %d).\n"
					, cl.num_statics, MAX_STATIC_ENTITIES);
			R_CheckEfrags();
		}
		CL_SignonReply();
		break;
	case svc_killedmonster:
		cl.stats[STAT_MONSTERS]++;
		break;
	case svc_foundsecret:
		cl.stats[STAT_SECRETS]++;
		break;
	case svc_updatestat:
		i = MSG_ReadByte();
		if(i < 0 || i >= MAX_CL_STATS)
			Sys_Error("svc_updatestat: %i is invalid", i);
		cl.stats[i] = MSG_ReadLong();;
		break;
	case svc_spawnstaticsound:
		CL_ParseStaticSound(1); //johnfitz -- added parameter
		break;
	case svc_cdtrack:
		cl.cdtrack = MSG_ReadByte();
		cl.looptrack = MSG_ReadByte();
		if((cls.demoplayback || cls.demorecording) && (cls.forcetrack != -1))
			CDAudio_Play((u8)cls.forcetrack, true);
		else
			CDAudio_Play((u8)cl.cdtrack, true);
		break;
	case svc_intermission:
		cl.intermission = 1;
		cl.completed_time = cl.time;
		vid.recalc_refdef = 1; // go to full screen
		V_RestoreAngles();
		break;
	case svc_finale:
		cl.intermission = 2;
		cl.completed_time = cl.time;
		vid.recalc_refdef = 1; // go to full screen
		str = MSG_ReadString();
		SCR_CenterPrint(str);
		V_RestoreAngles();
		break;
	case svc_cutscene:
		cl.intermission = 3;
		cl.completed_time = cl.time;
		vid.recalc_refdef = 1; // go to full screen
		str = MSG_ReadString();
		SCR_CenterPrint(str);
		V_RestoreAngles();
		break;
	case svc_sellscreen:
		Cmd_ExecuteString("help", src_command);
		break;
		//johnfitz -- new svc types
	case svc_skybox:
		Sky_LoadSkyBox(MSG_ReadString());
		break;
	case svc_bf:
		Cmd_ExecuteString("bf", src_command);
		break;
	case svc_fog:
		Fog_ParseServerMessage();
		break;
	case svc_spawnbaseline2: //PROTOCOL_FITZQUAKE
		i = MSG_ReadShort();
		CL_ParseBaseline(CL_EntityNum(i), 2);
		break;
	case svc_spawnstatic2: //PROTOCOL_FITZQUAKE
		CL_ParseStatic(2);
		break;
	case svc_spawnstaticsound2: //PROTOCOL_FITZQUAKE
		CL_ParseStaticSound(2);
		break; // used by the 2021 rerelease
	case svc_achievement:
		str = MSG_ReadString();
		Con_DPrintf("Ignoring svc_achievement(%s)\n", str);
		break;
	case svc_localsound:
		CL_ParseLocalSound();
		break;
	}
	lastcmd = cmd; //johnfitz
	}
}
