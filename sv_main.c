// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
// sv_main.c -- server main program
#include "quakedef.h"

static s8 localmodels[MAX_MODELS][8]; // inline model names for precache
static s32 sv_protocol = PROTOCOL_FITZQUAKE; //johnfitz
static s32 fatbytes;
static u8 *fatpvs;
static s32 fatpvs_capacity;

void SV_Protocol_f()
{
	s32 i; // keep here for OpenBSD
	switch(Cmd_Argc()) {
	case 1:
		Con_Printf("\"sv_protocol\" is \"%i\"\n", sv_protocol);
		break;
	case 2:
		i = atoi(Cmd_Argv(1));
		if(i != PROTOCOL_NETQUAKE && i != PROTOCOL_FITZQUAKE
				&& i != PROTOCOL_RMQ)
			Con_Printf("sv_protocol must be %i or %i or %i\n",
			PROTOCOL_NETQUAKE, PROTOCOL_FITZQUAKE, PROTOCOL_RMQ);
		else {
			sv_protocol = i;
			if(sv.active)
	Con_Printf("changes will not take effect until the next level load.\n");
		}
		break;
	default:
		Con_SafePrintf("usage: sv_protocol <protocol>\n");
		break;
	}
}

void SV_Init()
{
	sv.edicts = NULL; // ericw -- sv.edicts switched to use malloc()
	Cvar_RegisterVariable(&sv_maxvelocity);
	Cvar_RegisterVariable(&sv_gravity);
	Cvar_RegisterVariable(&sv_friction);
	Cvar_RegisterVariable(&sv_edgefriction);
	Cvar_RegisterVariable(&sv_stopspeed);
	Cvar_RegisterVariable(&sv_maxspeed);
	Cvar_RegisterVariable(&sv_accelerate);
	Cvar_RegisterVariable(&sv_idealpitchscale);
	Cvar_RegisterVariable(&sv_aim);
	Cvar_RegisterVariable(&sv_nostep);
	Cvar_RegisterVariable(&sv_freezenonclients);
	Cvar_RegisterVariable(&sv_altnoclip); //johnfitz
	Cmd_AddCommand("sv_protocol", &SV_Protocol_f); //johnfitz
	for(s32 i = 0; i < MAX_MODELS; i++)
		sprintf(localmodels[i], "*%i", i);
	s32 i = COM_CheckParm("-protocol");
	if(i && i < com_argc - 1) sv_protocol = atoi(com_argv[i + 1]);
	const s8 *p;
	switch(sv_protocol) {
	case PROTOCOL_NETQUAKE:
		p = "NetQuake";
		break;
	case PROTOCOL_FITZQUAKE:
		p = "FitzQuake";
		break;
	case PROTOCOL_RMQ:
		p = "RMQ";
		break;
	default:
Sys_Error("Bad protocol version request %i. Accepted values: %i, %i, %i.",
	sv_protocol, PROTOCOL_NETQUAKE, PROTOCOL_FITZQUAKE, PROTOCOL_RMQ);
		return; // silence compiler
	}
	Sys_Printf("Server using protocol %i(%s)\n", sv_protocol, p);
}

void SV_StartParticle(vec3_t org, vec3_t dir, s32 color, s32 count)
{ // Make sure the event gets sent to all clients
	if(sv.datagram.cursize > MAX_DATAGRAM-18) return;
	MSG_WriteByte(&sv.datagram, svc_particle);
	MSG_WriteCoord(&sv.datagram, org[0], sv.protocolflags);
	MSG_WriteCoord(&sv.datagram, org[1], sv.protocolflags);
	MSG_WriteCoord(&sv.datagram, org[2], sv.protocolflags);
	for(s32 i = 0; i < 3; i++) {
		s32 v = dir[i]*16;
		if(v > 127) v = 127;
		else if(v < -128) v = -128;
		MSG_WriteChar(&sv.datagram, v);
	}
	MSG_WriteByte(&sv.datagram, count);
	MSG_WriteByte(&sv.datagram, color);
}

// Each entity can have eight independant sound sources, like voice,
// weapon, feet, etc.
// Channel 0 is an auto-allocate channel, the others override anything
// already running on that entity/channel pair.
// An attenuation of 0 will play full volume everywhere in the level.
// Larger attenuations will drop off. (max 4 attenuation)
void SV_StartSound(edict_t *entity, s32 channel, const s8 *sample, s32 volume, f32 attenuation)
{
	if(volume < 0 || volume > 255)
		Host_Error("SV_StartSound: volume = %i", volume);
	if(attenuation < 0 || attenuation > 4)
		Host_Error("SV_StartSound: attenuation = %f", attenuation);
	if(channel < 0 || channel > 7)
		Host_Error("SV_StartSound: channel = %i", channel);
	if(sv.datagram.cursize > MAX_DATAGRAM-21)
		return;
	s32 sound_num = 1; // find precache number for sound
	for(;sound_num<MAX_SOUNDS && sv.sound_precache[sound_num]; sound_num++){
		if(!strcmp(sample, sv.sound_precache[sound_num]))
			break;
	}
	if(sound_num == MAX_SOUNDS || !sv.sound_precache[sound_num]) {
		Con_Printf("SV_StartSound: %s not precached\n", sample);
		return;
	}
	s32 ent = NUM_FOR_EDICT(entity);
	s32 field_mask = 0;
	if(volume != DEFAULT_SOUND_PACKET_VOLUME)
		field_mask |= SND_VOLUME;
	if(attenuation != DEFAULT_SOUND_PACKET_ATTENUATION)
		field_mask |= SND_ATTENUATION;
	if(ent >= 8192) {
		if(sv.protocol == PROTOCOL_NETQUAKE)
			return; //don't send any info protocol can't support
		field_mask |= SND_LARGEENTITY;
	}
	if(sound_num >= 256 || channel >= 8) {
		if(sv.protocol == PROTOCOL_NETQUAKE)
			return; //don't send any info protocol can't support
		field_mask |= SND_LARGESOUND;
	}
	if(sv.datagram.cursize > MAX_DATAGRAM-21) return;
	// directed messages go only to the entity the are targeted on
	MSG_WriteByte(&sv.datagram, svc_sound);
	MSG_WriteByte(&sv.datagram, field_mask);
	if(field_mask & SND_VOLUME)
		MSG_WriteByte(&sv.datagram, volume);
	if(field_mask & SND_ATTENUATION)
		MSG_WriteByte(&sv.datagram, attenuation*64);
	if(field_mask & SND_LARGEENTITY) {
		MSG_WriteShort(&sv.datagram, ent);
		MSG_WriteByte(&sv.datagram, channel);
	} else MSG_WriteShort(&sv.datagram, (ent<<3) | channel);
	if(field_mask & SND_LARGESOUND) MSG_WriteShort(&sv.datagram, sound_num);
	else MSG_WriteByte(&sv.datagram, sound_num);
	for(s32 i = 0; i < 3; i++)
		MSG_WriteCoord(&sv.datagram, entity->v.origin[i]+0.5*
			(entity->v.mins[i]+entity->v.maxs[i]),sv.protocolflags);
}

void SV_LocalSound(client_t *client, const s8 *sample)
{ // for 2021 rerelease
	s32 sound_num = 1;
	for(;sound_num<MAX_SOUNDS && sv.sound_precache[sound_num]; sound_num++){
		if(!strcmp(sample, sv.sound_precache[sound_num]))
			break;
	}
	if(sound_num == MAX_SOUNDS || !sv.sound_precache[sound_num]) {
		Con_Printf("SV_LocalSound: %s not precached\n", sample);
		return;
	}
	s32 field_mask = 0;
	if(sound_num >= 256) {
		if(sv.protocol == PROTOCOL_NETQUAKE) return;
		field_mask = SND_LARGESOUND;
	}
	if(client->message.cursize > client->message.maxsize-4) return;
	MSG_WriteByte(&client->message, svc_localsound);
	MSG_WriteByte(&client->message, field_mask);
	if(field_mask & SND_LARGESOUND)
		MSG_WriteShort(&client->message, sound_num);
	else MSG_WriteByte(&client->message, sound_num);
}

static bool SV_IsLocalClient(client_t *client)
{ return Q_strcmp(client->netconnection->address, "LOCAL") == 0; }

// Sends the first message from the server to a connected client.
// This will be sent on the initial connection and upon each server load.
void SV_SendServerinfo(client_t *client)
{
	s8 message[2048];
	MSG_WriteByte(&client->message, svc_print);
	sprintf(message, "%c\nFITZQUAKE %1.2f SERVER(%i CRC)\n",
			2, FITZQUAKE_VERSION, pr_crc);
	MSG_WriteString(&client->message,message);
	MSG_WriteByte(&client->message, svc_serverinfo);
	MSG_WriteLong(&client->message, sv.protocol);
// mh - send protocol flags so the client knows the protocol features to expect
	if(sv.protocol == PROTOCOL_RMQ)
		MSG_WriteLong(&client->message, sv.protocolflags);
	MSG_WriteByte(&client->message, svs.maxclients);
	if(!coop.value && deathmatch.value)
		MSG_WriteByte(&client->message, GAME_DEATHMATCH);
	else MSG_WriteByte(&client->message, GAME_COOP);
	MSG_WriteString(&client->message, PR_GetString(sv.edicts->v.message));
//johnfitz - only send the first 256 model and sound precaches if protocol is 15
	s32 i;
	const s8 **s;
	for(i = 1, s = sv.model_precache+1; *s; s++,i++)
		if(sv.protocol != PROTOCOL_NETQUAKE || i < 256)
			MSG_WriteString(&client->message, *s);
	MSG_WriteByte(&client->message, 0);
	for(i = 1, s = sv.sound_precache+1; *s; s++, i++)
		if(sv.protocol != PROTOCOL_NETQUAKE || i < 256)
			MSG_WriteString(&client->message, *s);
	MSG_WriteByte(&client->message, 0);
	MSG_WriteByte(&client->message, svc_cdtrack);
	MSG_WriteByte(&client->message, sv.edicts->v.sounds);
	MSG_WriteByte(&client->message, sv.edicts->v.sounds);
	MSG_WriteByte(&client->message, svc_setview);
	MSG_WriteShort(&client->message, NUM_FOR_EDICT(client->edict));
	MSG_WriteByte(&client->message, svc_signonnum);
	MSG_WriteByte(&client->message, 1);
	client->sendsignon = PRESPAWN_FLUSH;
	client->spawned = 0; // need prespawn, spawn, etc
}

// Initializes a client_t for a new net connection. This will only be called
// once for a player each game, not once for each level change.
void SV_ConnectClient(s32 clientnum)
{
	f32 spawn_parms[NUM_SPAWN_PARMS];
	client_t *client = svs.clients + clientnum;
	Con_DPrintf("Client %s connected\n", client->netconnection->address);
	s32 edictnum = clientnum+1;
	edict_t *ent = EDICT_NUM(edictnum);
	struct qsocket_s *netconnection = client->netconnection;
	if(sv.loadgame) // set up the client_t
		memcpy(spawn_parms, client->spawn_parms, sizeof(spawn_parms));
	memset(client, 0, sizeof(*client));
	client->netconnection = netconnection;
	strcpy(client->name, "unconnected");
	client->active = 1;
	client->spawned = 0;
	client->edict = ent;
	client->message.data = client->msgbuf;
	client->message.maxsize = sizeof(client->msgbuf);
	client->message.allowoverflow = 1; // we can catch it
	if(sv.loadgame)
		memcpy(client->spawn_parms, spawn_parms, sizeof(spawn_parms));
	else {
		// call the progs to get default spawn parms for the new client
		PR_ExecuteProgram(pr_global_struct->SetNewParms);
		for(s32 i = 0; i < NUM_SPAWN_PARMS; i++)
			client->spawn_parms[i] = (&pr_global_struct->parm1)[i];
	}
	SV_SendServerinfo(client);
}

void SV_CheckForNewClients()
{
	while(1) { // check for new connections
		struct qsocket_s *ret = NET_CheckNewConnections();
		if(!ret) break;
		s32 i = 0; // init a new client structure
		for(; i < svs.maxclients; i++)
			if(!svs.clients[i].active) break;
		if(i == svs.maxclients)
			Sys_Error("Host_CheckForNewClients: no free clients");
		svs.clients[i].netconnection = ret;
		SV_ConnectClient(i);
		net_activeconnections++;
	}
}

void SV_ClearDatagram() { SZ_Clear(&sv.datagram); }

// The PVS must include a small area around the client to allow head bobbing
// or other small motion on the client side. Otherwise, a bob might cause an
// entity that should be visible to not show up, especially when the bob
// crosses a waterline.
void SV_AddToFatPVS(vec3_t org, mnode_t *node, model_t *worldmodel)
{
	while(1) {
		if(node->contents < 0) { // if leaf, accumulate the pvs bits
			if(node->contents != CONTENTS_SOLID) {
				u8 *pvs=Mod_LeafPVS((mleaf_t*)node,worldmodel);
				for(s32 i = 0; i < fatbytes; i++)
					fatpvs[i] |= pvs[i];
			}
			return;
		}
		mplane_t *plane = node->plane;
		f32 d = DotProduct(org, plane->normal) - plane->dist;
		if(d > 8) node = node->children[0];
		else if(d < -8) node = node->children[1];
		else { // go down both
			SV_AddToFatPVS(org, node->children[0], worldmodel);
			node = node->children[1];
		}
	}
}


u8 *SV_FatPVS(vec3_t org, model_t *worldmodel) // Calculates a PVS that is the
{ // inclusive or of all leafs within 8 pixels of the given point.
	fatbytes = (worldmodel->numleafs+7)>>3; // ericw - was +31, assumed typo
	if(fatpvs == NULL || fatbytes > fatpvs_capacity) {
		fatpvs_capacity = fatbytes;
		fatpvs = (u8 *) realloc(fatpvs, fatpvs_capacity);
		if(!fatpvs)
	  Sys_Error("SV_FatPVS: realloc() failed on %d bytes", fatpvs_capacity);
	}
	Q_memset(fatpvs, 0, fatbytes);
	SV_AddToFatPVS(org, worldmodel->nodes, worldmodel);
	return fatpvs;
}

bool SV_VisibleToClient(edict_t *client, edict_t *test, model_t *worldmodel)
{ // PVS test encapsulated in a nice function
	vec3_t org;
	VectorAdd(client->v.origin, client->v.view_ofs, org);
	u8 *pvs = SV_FatPVS(org, worldmodel);
	for(s32 i = 0; i < test->num_leafs; i++)
		if(pvs[test->leafnums[i] >> 3] & (1 << (test->leafnums[i]&7)))
			return 1;
	return 0;
}

void SV_WriteEntitiesToClient(edict_t *clent, sizebuf_t *msg)
{
	vec3_t org; // find the client's PVS
	VectorAdd(clent->v.origin, clent->v.view_ofs, org);
	u8 *pvs = SV_FatPVS(org, sv.worldmodel);
	// send over all entities(excpet the client) that touch the pvs
	edict_t *ent = NEXT_EDICT(sv.edicts);
	for(s32 e = 1; e < sv.num_edicts; e++, ent = NEXT_EDICT(ent)) {
		s32 i;
		if(ent != clent){ // clent is ALLWAYS sent
			// ignore ents without visible models
			if(!ent->v.modelindex || !PR_GetString(ent->v.model)[0])
				continue;
			if(sv.protocol == PROTOCOL_NETQUAKE &&
					(s32)ent->v.modelindex & 0xFF00)
				continue;
			// ignore if not touching a PV leaf
			for(i = 0; i < ent->num_leafs; i++)
				if(pvs[ent->leafnums[i] >> 3] &
					(1 << (ent->leafnums[i]&7)))
						break;
	// ericw -- added ent->num_leafs < MAX_ENT_LEAFS condition.
	// if ent->num_leafs == MAX_ENT_LEAFS, the ent is visible from too many
	// leafs for us to say whether it's in the PVS, so don't try to vis cull
	// it. this commonly happens with rotators, because they often have huge
	// bboxes spanning the entire map, or really tall lifts, etc.
			if(i == ent->num_leafs && ent->num_leafs < MAX_ENT_LEAFS)
				continue; // not visible
		}
	// johnfitz - max size for protocol 15 is 18 bytes, not 16 as originally
	// assumed here. And, for protocol 85 the max size is actually 24 bytes.
	// For f32 coords and angles the limit is 40.
	// FIXME: Use tighter limit according to protocol flags and send bits.
		if(msg->cursize + 40 > msg->maxsize)
			Con_DPrintf("Packet overflow!\n");
		s32 bits = 0; // send an update
		for(i = 0; i < 3; i++) {
			f32 miss = ent->v.origin[i] - ent->baseline.origin[i];
			if(miss < -0.1 || miss > 0.1) bits |= U_ORIGIN1<<i;
		}
		if(ent->v.angles[0] != ent->baseline.angles[0])bits |= U_ANGLE1;
		if(ent->v.angles[1] != ent->baseline.angles[1])bits |= U_ANGLE2;
		if(ent->v.angles[2] != ent->baseline.angles[2])bits |= U_ANGLE3;
		if(ent->v.movetype == MOVETYPE_STEP) bits |= U_STEP;
		if(ent->baseline.colormap != ent->v.colormap)bits |= U_COLORMAP;
		if(ent->baseline.skin != ent->v.skin) bits |= U_SKIN;
		if(ent->baseline.frame != ent->v.frame) bits |= U_FRAME;
		if((ent->baseline.effects^(s32)ent->v.effects)&pr_effects_mask)
			bits |= U_EFFECTS;
		if(ent->baseline.modelindex!=ent->v.modelindex)bits |= U_MODEL;
		eval_t *val;
		if(pr_alpha_supported) {
			val = GetEdictFieldValue(ent, "alpha");
			if(val) ent->alpha = ENTALPHA_ENCODE(val->_float);
		}
		//don't send invisible entities unless they have effects
		if(ent->alpha == ENTALPHA_ZERO &&
			!((s32)ent->v.effects & pr_effects_mask))
				continue;
		val = GetEdictFieldValue(ent, "scale");
		if(val) ent->scale = ENTSCALE_ENCODE(val->_float);
		else ent->scale = ENTSCALE_DEFAULT;
		if(sv.protocol != PROTOCOL_NETQUAKE) {
			if(ent->baseline.alpha != ent->alpha) bits |= U_ALPHA;
			if(ent->baseline.scale != ent->scale) bits |= U_SCALE;
			if(bits & U_FRAME && (s32)ent->v.frame & 0xFF00)
				bits |= U_FRAME2;
			if(bits & U_MODEL && (s32)ent->v.modelindex & 0xFF00)
				bits |= U_MODEL2;
			if(ent->sendinterval) bits |= U_LERPFINISH;
			if(bits >= 65536) bits |= U_EXTEND1;
			if(bits >= 16777216) bits |= U_EXTEND2;
		}
		if(e >= 256) bits |= U_LONGENTITY;
		if(bits >= 256) bits |= U_MOREBITS;
		MSG_WriteByte(msg, bits | U_SIGNAL); // write the message
		if(bits & U_MOREBITS) MSG_WriteByte(msg, bits>>8);
		if(bits & U_EXTEND1) MSG_WriteByte(msg, bits>>16);
		if(bits & U_EXTEND2) MSG_WriteByte(msg, bits>>24);
		if(bits & U_LONGENTITY) MSG_WriteShort(msg,e);
		else MSG_WriteByte(msg,e);
		if(bits & U_MODEL) MSG_WriteByte(msg, ent->v.modelindex);
		if(bits & U_FRAME) MSG_WriteByte(msg, ent->v.frame);
		if(bits & U_COLORMAP) MSG_WriteByte(msg, ent->v.colormap);
		if(bits & U_SKIN) MSG_WriteByte(msg, ent->v.skin);
		if(bits & U_EFFECTS)
			MSG_WriteByte(msg, (s32)ent->v.effects&pr_effects_mask);
		if(bits & U_ORIGIN1)
			MSG_WriteCoord(msg, ent->v.origin[0], sv.protocolflags);
		if(bits & U_ANGLE1)
			MSG_WriteAngle(msg, ent->v.angles[0], sv.protocolflags);
		if(bits & U_ORIGIN2)
			MSG_WriteCoord(msg, ent->v.origin[1], sv.protocolflags);
		if(bits & U_ANGLE2)
			MSG_WriteAngle(msg, ent->v.angles[1], sv.protocolflags);
		if(bits & U_ORIGIN3)
			MSG_WriteCoord(msg, ent->v.origin[2], sv.protocolflags);
		if(bits & U_ANGLE3)
			MSG_WriteAngle(msg, ent->v.angles[2], sv.protocolflags);
		if(bits & U_ALPHA) MSG_WriteByte(msg, ent->alpha);
		if(bits & U_SCALE) MSG_WriteByte(msg, ent->scale);
		if(bits & U_FRAME2) MSG_WriteByte(msg, (s32)ent->v.frame >> 8);
		if(bits & U_MODEL2)MSG_WriteByte(msg,(s32)ent->v.modelindex>>8);
		if(bits & U_LERPFINISH)
			MSG_WriteByte(msg,
				(u8)(Q_rint((ent->v.nextthink-sv.time)*255)));
	}
}

void SV_CleanupEnts()
{
	edict_t *ent = NEXT_EDICT(sv.edicts);
	for(s32 e = 1; e < sv.num_edicts; e++, ent = NEXT_EDICT(ent))
		ent->v.effects = (s32)ent->v.effects & ~EF_MUZZLEFLASH;
}

void SV_WriteClientdataToMessage(edict_t *ent, sizebuf_t *msg)
{
	if(ent->v.dmg_take || ent->v.dmg_save) { // send a damage message
		edict_t *other = PROG_TO_EDICT(ent->v.dmg_inflictor);
		MSG_WriteByte(msg, svc_damage);
		MSG_WriteByte(msg, ent->v.dmg_save);
		MSG_WriteByte(msg, ent->v.dmg_take);
		for(s32 i = 0; i < 3; i++)
			MSG_WriteCoord(msg, other->v.origin[i] + 
		   0.5*(other->v.mins[i] + other->v.maxs[i]), sv.protocolflags);
		ent->v.dmg_take = 0;
		ent->v.dmg_save = 0;
	}
	// send the current viewpos offset from the view entity
	SV_SetIdealPitch(); // how much to look up / down ideally a fixangle
			    // might get lost in a dropped packet. Oh well.
	if(ent->v.fixangle) {
		MSG_WriteByte(msg, svc_setangle);
		for(s32 i = 0; i < 3; i++)
			MSG_WriteAngle(msg, ent->v.angles[i], sv.protocolflags);
		ent->v.fixangle = 0;
	}
	s32 bits = 0;
	if(ent->v.view_ofs[2] != DEFAULT_VIEWHEIGHT) bits |= SU_VIEWHEIGHT;
	if(ent->v.idealpitch) bits |= SU_IDEALPITCH;
// stuff sigil bits into high bits of items for sbar, or else mix in items2
	eval_t *val = GetEdictFieldValue(ent, "items2");
	s32 items;
	if(val) items = (s32)ent->v.items | ((s32)val->_float << 23);
	else items = (s32)ent->v.items|((s32)pr_global_struct->serverflags<<28);
	bits |= SU_ITEMS;
	if((s32)ent->v.flags & FL_ONGROUND) bits |= SU_ONGROUND;
	if(ent->v.waterlevel >= 2) bits |= SU_INWATER;
	for(s32 i = 0; i < 3; i++) {
		if(ent->v.punchangle[i]) bits |= (SU_PUNCH1<<i);
		if(ent->v.velocity[i]) bits |= (SU_VELOCITY1<<i);
	}
	if(ent->v.weaponframe) bits |= SU_WEAPONFRAME;
	if(ent->v.armorvalue) bits |= SU_ARMOR;
	bits |= SU_WEAPON;
	if(sv.protocol != PROTOCOL_NETQUAKE) {
		if(bits & SU_WEAPON && SV_ModelIndex(PR_GetString
			(ent->v.weaponmodel)) & 0xFF00) bits |= SU_WEAPON2;
		if((s32)ent->v.armorvalue & 0xFF00) bits |= SU_ARMOR2;
		if((s32)ent->v.currentammo & 0xFF00) bits |= SU_AMMO2;
		if((s32)ent->v.ammo_shells & 0xFF00) bits |= SU_SHELLS2;
		if((s32)ent->v.ammo_nails & 0xFF00) bits |= SU_NAILS2;
		if((s32)ent->v.ammo_rockets & 0xFF00) bits |= SU_ROCKETS2;
		if((s32)ent->v.ammo_cells & 0xFF00) bits |= SU_CELLS2;
		if(bits & SU_WEAPONFRAME && (s32)ent->v.weaponframe & 0xFF00)
						bits |= SU_WEAPONFRAME2;
		if(bits & SU_WEAPON && ent->alpha != ENTALPHA_DEFAULT)
						bits |= SU_WEAPONALPHA;
		if(bits >= 65536) bits |= SU_EXTEND1;
		if(bits >= 16777216) bits |= SU_EXTEND2;
	}
	MSG_WriteByte(msg, svc_clientdata); // send the data
	MSG_WriteShort(msg, bits);
	if(bits & SU_EXTEND1) MSG_WriteByte(msg, bits>>16);
	if(bits & SU_EXTEND2) MSG_WriteByte(msg, bits>>24);
	if(bits & SU_VIEWHEIGHT) MSG_WriteChar(msg, ent->v.view_ofs[2]);
	if(bits & SU_IDEALPITCH) MSG_WriteChar(msg, ent->v.idealpitch);
	for(s32 i = 0; i < 3; i++) {
		if(bits & (SU_PUNCH1<<i))
			MSG_WriteChar(msg, ent->v.punchangle[i]);
		if(bits & (SU_VELOCITY1<<i))
			MSG_WriteChar(msg, ent->v.velocity[i]/16);
	}
	MSG_WriteLong(msg, items); // [always sent] if(bits & SU_ITEMS)
	if(bits & SU_WEAPONFRAME) MSG_WriteByte(msg, ent->v.weaponframe);
	if(bits & SU_ARMOR) MSG_WriteByte(msg, ent->v.armorvalue);
	if(bits & SU_WEAPON)
	    MSG_WriteByte(msg, SV_ModelIndex(PR_GetString(ent->v.weaponmodel)));
	MSG_WriteShort(msg, ent->v.health);
	MSG_WriteByte(msg, ent->v.currentammo);
	MSG_WriteByte(msg, ent->v.ammo_shells);
	MSG_WriteByte(msg, ent->v.ammo_nails);
	MSG_WriteByte(msg, ent->v.ammo_rockets);
	MSG_WriteByte(msg, ent->v.ammo_cells);
	if(standard_quake) { MSG_WriteByte(msg, ent->v.weapon); }
	else for(s32 i = 0; i < 32; i++)
			if(((s32)ent->v.weapon) & (1<<i)) {
				MSG_WriteByte(msg, i);
				break;
			}
	if(bits & SU_WEAPON2) MSG_WriteByte(msg,
			SV_ModelIndex(PR_GetString(ent->v.weaponmodel)) >> 8);
	if(bits & SU_ARMOR2) MSG_WriteByte(msg, (s32)ent->v.armorvalue >> 8);
	if(bits & SU_AMMO2) MSG_WriteByte(msg, (s32)ent->v.currentammo >> 8);
	if(bits & SU_SHELLS2) MSG_WriteByte(msg, (s32)ent->v.ammo_shells >> 8);
	if(bits & SU_NAILS2) MSG_WriteByte(msg, (s32)ent->v.ammo_nails >> 8);
	if(bits & SU_ROCKETS2)MSG_WriteByte(msg, (s32)ent->v.ammo_rockets >> 8);
	if(bits & SU_CELLS2) MSG_WriteByte(msg, (s32)ent->v.ammo_cells >> 8);
	if(bits & SU_WEAPONFRAME2)MSG_WriteByte(msg,(s32)ent->v.weaponframe>>8);
	if(bits & SU_WEAPONALPHA) MSG_WriteByte(msg, ent->alpha);
}

bool SV_SendClientDatagram(client_t *client)
{
	u8 buf[MAX_DATAGRAM];
	sizebuf_t msg;
	msg.data = buf;
	msg.maxsize = sizeof(buf);
	msg.cursize = 0;
// if client is nonlocal, use smaller max size so packets aren't fragmented
	if(Q_strcmp(client->netconnection->address, "LOCAL") != 0) //johnfitz
		msg.maxsize = DATAGRAM_MTU;
	MSG_WriteByte(&msg, svc_time);
	MSG_WriteFloat(&msg, sv.time);
	// add the client specific data to the datagram
	SV_WriteClientdataToMessage(client->edict, &msg);
	SV_WriteEntitiesToClient(client->edict, &msg);
	// copy the server datagram if there is space
	if(msg.cursize + sv.datagram.cursize < msg.maxsize)
		SZ_Write(&msg, sv.datagram.data, sv.datagram.cursize);
	// send the datagram
	if(NET_SendUnreliableMessage(client->netconnection, &msg) == -1) {
		SV_DropClient(1);// if the message couldn't send, kick off
		return 0;
	}
	return 1;
}

void SV_UpdateToReliableMessages()
{
	s32 i, j;
	client_t *cl;
	// check for changes to be sent over the reliable streams
	for(i=0,host_client=svs.clients; i<svs.maxclients; i++,host_client++){
		if(host_client->old_frags != host_client->edict->v.frags) {
			for(j=0, cl = svs.clients; j<svs.maxclients; j++, cl++){
				if(!cl->active) continue;
				MSG_WriteByte(&cl->message, svc_updatefrags);
				MSG_WriteByte(&cl->message, i);
				MSG_WriteShort(&cl->message,
						host_client->edict->v.frags);
			}
			host_client->old_frags = host_client->edict->v.frags;
		}
	}
	for(j=0, cl = svs.clients ; j<svs.maxclients ; j++, cl++) {
		if(!cl->active) continue;
		SZ_Write(&cl->message, sv.reliable_datagram.data,
				sv.reliable_datagram.cursize);
	}
	SZ_Clear(&sv.reliable_datagram);
}

void SV_SendNop(client_t *client) // Send a nop message without trashing or
{ // sending the accumulated client message buffer
	sizebuf_t msg;
	u8 buf[4];
	msg.data = buf;
	msg.maxsize = sizeof(buf);
	msg.cursize = 0;
	MSG_WriteChar(&msg, svc_nop);
	if(NET_SendUnreliableMessage(client->netconnection, &msg) == -1)
		SV_DropClient(1); // if the message couldn't send, kick off
	client->last_message = realtime;
}

void SV_SendClientMessages()
{
	SV_UpdateToReliableMessages(); // update frags, names, etc
	// build individual updates
	s32 i = 0;
	for(host_client = svs.clients; i < svs.maxclients; i++, host_client++){
		if(!host_client->active) continue;
		if(host_client->spawned) {
			if(!SV_SendClientDatagram(host_client))
				continue;
		} else {
			// the player isn't totally in the game yet send small
			// keepalive messages if too much time has passed send a
			// full message when the next signon stage has been
			// requested some other message data(name changes, etc)
			// may accumulate between signon stages
			if(!host_client->sendsignon) {
				if(realtime - host_client->last_message > 5)
					SV_SendNop(host_client);
				continue; // don't send out non-signon messages
			}
			if(host_client->sendsignon == PRESPAWN_SIGNONBUFS) {
				bool local = SV_IsLocalClient(host_client);
				while(host_client->signonidx < sv.num_signon_buffers) {
					sizebuf_t *signon = sv.signon_buffers[host_client->signonidx];
					if(host_client->message.cursize + signon->cursize > host_client->message.maxsize)
						break;
					SZ_Write(&host_client->message, signon->data, signon->cursize);
					host_client->signonidx++;
					// only send multiple buffers at once when playing locally,
					// otherwise we send one signon at a time to avoid overflowing
					// the datagram buffer for clients using a lower limit(e.g. 32000 in QS)
					if(!local) break;
				}
				if(host_client->signonidx == sv.num_signon_buffers)
					host_client->sendsignon = PRESPAWN_SIGNONMSG;
			}
			if(host_client->sendsignon == PRESPAWN_SIGNONMSG) {
				if(host_client->message.cursize + 2 < host_client->message.maxsize) {
					MSG_WriteByte(&host_client->message, svc_signonnum);
					MSG_WriteByte(&host_client->message, 2);
					host_client->sendsignon=PRESPAWN_FLUSH;
				}
			}
		}
		// check for an overflowed message. Should only happen
		// on a very fucked up connection that backs up a lot, then
		// changes level
		if(host_client->message.overflowed) {
			SV_DropClient(1);
			host_client->message.overflowed = 0;
			continue;
		}
		if(host_client->message.cursize || host_client->dropasap) {
			if(!NET_CanSendMessage(host_client->netconnection)) {
				continue;
			}
			if(host_client->dropasap)
				SV_DropClient(0); // went to another level
			else {
				if(NET_SendMessage(host_client->netconnection,
						&host_client->message) == -1)
					SV_DropClient(1); // if the message couldn't send, kick off
				SZ_Clear(&host_client->message);
				host_client->last_message = realtime;
				if(host_client->sendsignon == PRESPAWN_FLUSH)
					host_client->sendsignon = PRESPAWN_DONE;
			}
		}
	}
	SV_CleanupEnts(); // clear muzzle flashes
}

static void SV_AddSignonBuffer()
{
	if(sv.num_signon_buffers >= MAX_SIGNON_BUFFERS)
		Host_Error("SV_AddSignonBuffer overflow\n");
	sizebuf_t *sb = Hunk_AllocName(sizeof(sizebuf_t)+SIGNON_SIZE,"signon");
	sb->data = (u8 *)(sb + 1);
	sb->maxsize = SIGNON_SIZE;
	sv.signon_buffers[sv.num_signon_buffers++] = sb;
	sv.signon = sb;
}

void SV_ReserveSignonSpace(s32 numbytes)
{ if(sv.signon->cursize + numbytes > sv.signon->maxsize) SV_AddSignonBuffer(); }

s32 SV_ModelIndex(const s8 *name)
{
	if(!name || !name[0]) return 0;
	s32 i = 0;
	for(;i < MAX_MODELS && sv.model_precache[i]; i++)
		if(!strcmp(sv.model_precache[i], name)) return i;
	if(i==MAX_MODELS || !sv.model_precache[i])
		Sys_Error("SV_ModelIndex: model %s not precached", name);
	return i;
}

void SV_CreateBaseline()
{
	edict_t *sve;
	for(s32 entnum = 0; entnum < sv.num_edicts ; entnum++) {
		sve = EDICT_NUM(entnum); // get the current server version
		if(sve->free) continue;
		if(entnum > svs.maxclients && !sve->v.modelindex) continue;
		VectorCopy(sve->v.origin, sve->baseline.origin);//create ent
		VectorCopy(sve->v.angles, sve->baseline.angles);//baseline
		sve->baseline.frame = sve->v.frame;
		sve->baseline.skin = sve->v.skin;
		if(entnum > 0 && entnum <= svs.maxclients) {
			sve->baseline.colormap = entnum;
			sve->baseline.modelindex =
				SV_ModelIndex("progs/player.mdl");
			sve->baseline.alpha = ENTALPHA_DEFAULT;
			sve->baseline.scale = ENTSCALE_DEFAULT;
		} else {
			sve->baseline.colormap = 0;
			sve->baseline.modelindex =
				SV_ModelIndex(PR_GetString(sve->v.model));
			sve->baseline.alpha = sve->alpha;
			sve->baseline.scale = ENTSCALE_DEFAULT;
			if(sv.protocol == PROTOCOL_RMQ) {
				eval_t* val = GetEdictFieldValue(sve, "scale");
				if(val) sve->baseline.scale =
					ENTSCALE_ENCODE(val->_float);
			}
		}
		s32 bits = 0;
		if(sv.protocol == PROTOCOL_NETQUAKE) {
			if(sve->baseline.modelindex & 0xFF00)
				sve->baseline.modelindex = 0;
			if(sve->baseline.frame & 0xFF00)
				sve->baseline.frame = 0;
			sve->baseline.alpha = ENTALPHA_DEFAULT;
			sve->baseline.scale = ENTSCALE_DEFAULT;
		} else { //decide which extra data needs to be sent
			if(sve->baseline.modelindex&0xFF00)bits|=B_LARGEMODEL;
			if(sve->baseline.frame & 0xFF00) bits |= B_LARGEFRAME;
			if(sve->baseline.alpha!=ENTALPHA_DEFAULT)bits|=B_ALPHA;
			if(sve->baseline.scale!=ENTSCALE_DEFAULT)bits|=B_SCALE;
		}
		// add to the message
		SV_ReserveSignonSpace(35);
		if(bits)
			MSG_WriteByte(sv.signon, svc_spawnbaseline2);
		else MSG_WriteByte(sv.signon, svc_spawnbaseline);
		MSG_WriteShort(sv.signon,entnum);
		if(bits) MSG_WriteByte(sv.signon, bits);
		if(bits & B_LARGEMODEL)
			MSG_WriteShort(sv.signon, sve->baseline.modelindex);
		else MSG_WriteByte(sv.signon, sve->baseline.modelindex);
		if(bits & B_LARGEFRAME)
			MSG_WriteShort(sv.signon, sve->baseline.frame);
		else MSG_WriteByte(sv.signon, sve->baseline.frame);
		MSG_WriteByte(sv.signon, sve->baseline.colormap);
		MSG_WriteByte(sv.signon, sve->baseline.skin);
		for(s32 i = 0; i < 3; i++) {
			MSG_WriteCoord(sv.signon, sve->baseline.origin[i],
					sv.protocolflags);
			MSG_WriteAngle(sv.signon, sve->baseline.angles[i],
					sv.protocolflags);
		}
		if(bits&B_ALPHA)MSG_WriteByte(sv.signon, sve->baseline.alpha);
		if(bits&B_SCALE)MSG_WriteByte(sv.signon, sve->baseline.scale);
	}
}

void SV_SendReconnect()
{ // Tell all the clients that the server is changing levels
	u8 data[128];
	sizebuf_t msg;
	msg.data = data;
	msg.cursize = 0;
	msg.maxsize = sizeof(data);
	MSG_WriteChar(&msg, svc_stufftext);
	MSG_WriteString(&msg, "reconnect\n");
	NET_SendToAll(&msg, 5.0);
	if(!isDedicated) Cmd_ExecuteString("reconnect\n", src_command);
}


void SV_SaveSpawnparms() // Grabs the current state of each client for saving
{ // across the transition to another level
	svs.serverflags = pr_global_struct->serverflags;
	s32 i = 0;
	for(host_client = svs.clients; i < svs.maxclients; i++, host_client++){
		if(!host_client->active) continue;
		// call the progs to get default spawn parms for the new client
		pr_global_struct->self = EDICT_TO_PROG(host_client->edict);
		PR_ExecuteProgram(pr_global_struct->SetChangeParms);
		for(s32 j = 0; j < NUM_SPAWN_PARMS; j++)
		    host_client->spawn_parms[j] = (&pr_global_struct->parm1)[j];
	}
}

void SV_SpawnServer(const s8 *server)
{ // This is called at the start of each level
	static s8 dummy[8] = { 0,0,0,0,0,0,0,0 };
	if(hostname.string[0] == 0) Cvar_Set("hostname", "UNNAMED");
	scr_centertime_off = 0;
	lit_loaded = 0;
	Con_DPrintf("SpawnServer: %s\n",server);
	svs.changelevel_issued = 0; // now safe to issue another
	// tell all connected clients that we are going to a new level
	if(sv.active) SV_SendReconnect();
	if(coop.value) Cvar_Set("deathmatch", "0"); // make cvars consistant
	current_skill = (s32)(skill.value + 0.5);
	if(current_skill < 0) current_skill = 0;
	if(current_skill > 3) current_skill = 3;
	Cvar_SetValue("skill", (f32)current_skill);
	Host_ClearMemory(); // set up the new server
	q_strlcpy(sv.name, server, sizeof(sv.name));
	sv.protocol = sv_protocol;
	if(sv.protocol == PROTOCOL_RMQ)
		sv.protocolflags = PRFL_INT32COORD | PRFL_SHORTANGLE;
	else sv.protocolflags = 0;
	PR_LoadProgs(); // load progs to get entity field count
	// allocate server memory
	sv.max_edicts = CLAMP(MIN_EDICTS,(s32)max_edicts.value,MAX_EDICTS);
	sv.edicts = (edict_t *) malloc(sv.max_edicts*pr_edict_size);
	sv.datagram.maxsize = sizeof(sv.datagram_buf);
	sv.datagram.cursize = 0;
	sv.datagram.data = sv.datagram_buf;
	sv.reliable_datagram.maxsize = sizeof(sv.reliable_datagram_buf);
	sv.reliable_datagram.cursize = 0;
	sv.reliable_datagram.data = sv.reliable_datagram_buf;
	SV_AddSignonBuffer();
	// leave slots at start for clients only
	sv.num_edicts = svs.maxclients+1;
	memset(sv.edicts, 0, sv.num_edicts*pr_edict_size);
	for(s32 i = 0; i < svs.maxclients; i++) {
		edict_t *ent = EDICT_NUM(i+1);
		svs.clients[i].edict = ent;
	}
	sv.state = ss_loading;
	sv.paused = 0;
	sv.nomonsters = (nomonsters.value != 0.f);
	sv.time = 1.0;
	q_strlcpy(sv.name, server, sizeof(sv.name));
	snprintf(sv.modelname, sizeof(sv.modelname), "maps/%s.bsp", server);
	sv.worldmodel = Mod_ForName(sv.modelname, 0);
	if(!sv.worldmodel) {
		Con_Printf("Couldn't spawn server %s\n", sv.modelname);
		sv.active = 0;
		return;
	}
	sv.models[1] = sv.worldmodel;
	SV_ClearWorld(); // clear world interaction links
	sv.sound_precache[0] = dummy;
	sv.model_precache[0] = dummy;
	sv.model_precache[1] = sv.modelname;
	for(s32 i = 1; i < sv.worldmodel->numsubmodels; i++) {
		sv.model_precache[1+i] = localmodels[i];
		sv.models[i+1] = Mod_ForName(localmodels[i], 0);
	}
	edict_t *ent = EDICT_NUM(0); // load the rest of the entities
	memset(&ent->v, 0, progs->entityfields * 4);
	ent->free = 0;
	ent->v.model = PR_SetEngineString(sv.worldmodel->name);
	ent->v.modelindex = 1; // world model
	ent->v.solid = SOLID_BSP;
	ent->v.movetype = MOVETYPE_PUSH;
	if(coop.value) pr_global_struct->coop = coop.value;
	else pr_global_struct->deathmatch = deathmatch.value;
	pr_global_struct->mapname = PR_SetEngineString(sv.name);
	// serverflags are for cross level information(sigils)
	pr_global_struct->serverflags = svs.serverflags;
	ED_LoadFromFile(sv.worldmodel->entities);
	sv.active = 1;
	// all setup is completed, any further precache statements are errors
	sv.state = ss_active;
	// run two frames to allow everything to settle
	host_frametime = 0.1;
	SV_Physics();
	SV_Physics();
	// create a baseline for more efficient communications
	SV_CreateBaseline();
	s32 i = 0, signonsize = 0;
	for(; i < sv.num_signon_buffers; i++)
		signonsize += sv.signon_buffers[i]->cursize;
	if(signonsize > 64000-2)
printf("%i u8 signon buffer exceeds QS limit of 63998.\n", signonsize);
	else if(signonsize > 8000-2)
printf("%i u8 signon buffer exceeds standard limit of 7998.\n", signonsize);
	// send serverinfo to all connected clients
	for(i=0,host_client = svs.clients; i<svs.maxclients; i++, host_client++)
		if(host_client->active) SV_SendServerinfo(host_client);
	Con_DPrintf("Server spawned.\n");
}
