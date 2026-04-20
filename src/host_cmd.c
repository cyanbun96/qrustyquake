// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2007-2008 Kristian Duske
// Copyright(C) 2010-2014 QuakeSpasm developers
#include "quakedef.h"

static s32 maxlevelnamelen = 0;
static s32 maxmodnamelen = 0;
static const char *const knownmods[][2] = {
	{"id1",         "Quake"},
	{"hipnotic",    "Scourge of Armagon"},
	{"rogue",       "Dissolution of Eternity"},
	{"dopa",        "Dimension of the Past"},
	{"mg1",         "Dimension of the Machine"},
	{"q64",         "Quake(Nintendo 64)"},
	{"ctf",         "Capture The Flag"},
	{"udob",        "Underdark Overbright"},
	{"ad",          "Arcane Dimensions"},
	{"quoth",       "Quoth"},
	{"impel",       "Abyss of Pandemonium"},
	{"alk1.2",      "Alkaline 1.2"},
	{"rubicon2",    "Rubicon 2"},
	{"zer",         "Zerstorer - Testament of the Destroyer"},
};


void Host_Quit_f()
{
	if(key_dest != key_console && cls.state != ca_dedicated){
		M_Menu_Quit_f();
		return;
	}
	CL_Disconnect();
	Host_ShutdownServer(0);
	Sys_Quit();
}

void Host_Status_f()
{
	void(*print)(const s8 *fmt, ...);
	if(!sv.active){
		Cmd_ForwardToServer();
		return;
	}
	print = Con_Printf;
	print("host: %s\n", Cvar_VariableString("hostname"));
	print("version: %4.2f\n", VERSION);
	if(tcpipAvailable) print("tcp/ip: %s\n", my_tcpip_address);
	print("map: %s\n", sv.name);
	print("players: %i active(%i max)\n\n", net_activeconnections,
			svs.maxclients);
	s32 j = 0;
	for(client_t *client = svs.clients; j < svs.maxclients; j++, client++){
		if(!client->active) continue;
		s32 seconds=(s32)(net_time-client->netconnection->connecttime);
		s32 minutes = seconds / 60;
		s32 hours = 0;
		if(minutes){
			seconds -= (minutes * 60);
			hours = minutes / 60;
			if(hours) minutes -= (hours * 60);
		}
		print("#%-2u %-16.16s %3i %2i:%02i:%02i\n", j + 1, client->name,
			(s32)client->edict->v.frags, hours, minutes, seconds);
		print(" %s\n", client->netconnection->address);
	}
}

void Host_God_f() // Sets client to godmode
{
	if(cmd_source == src_command){
		Cmd_ForwardToServer();
		return;
	}
	if(pr_global_struct->deathmatch) return;
	sv_player->v.flags = (s32)sv_player->v.flags ^ FL_GODMODE;
	if(!((s32)sv_player->v.flags & FL_GODMODE))
		SV_ClientPrintf("godmode OFF\n");
	else SV_ClientPrintf("godmode ON\n");
}

void Host_Notarget_f()
{
	if(cmd_source == src_command){
		Cmd_ForwardToServer();
		return;
	}
	if(pr_global_struct->deathmatch) return;
	sv_player->v.flags = (s32)sv_player->v.flags ^ FL_NOTARGET;
	if(!((s32)sv_player->v.flags & FL_NOTARGET))
		SV_ClientPrintf("notarget OFF\n");
	else SV_ClientPrintf("notarget ON\n");
}
void Host_Noclip_f()
{
	if(cmd_source == src_command){
		Cmd_ForwardToServer();
		return;
	}
	if(pr_global_struct->deathmatch) return;
	if(sv_player->v.movetype != MOVETYPE_NOCLIP){
		noclip_anglehack = 1;
		sv_player->v.movetype = MOVETYPE_NOCLIP;
		SV_ClientPrintf("noclip ON\n");
	} else {
		noclip_anglehack = 0;
		sv_player->v.movetype = MOVETYPE_WALK;
		SV_ClientPrintf("noclip OFF\n");
	}
}

void Host_Fly_f()
{ // Sets client to flymode
	if(cmd_source == src_command){
		Cmd_ForwardToServer();
		return;
	}
	if(pr_global_struct->deathmatch) return;
	if(sv_player->v.movetype != MOVETYPE_FLY){
		sv_player->v.movetype = MOVETYPE_FLY;
		SV_ClientPrintf("flymode ON\n");
	} else {
		sv_player->v.movetype = MOVETYPE_WALK;
		SV_ClientPrintf("flymode OFF\n");
	}
}

void Host_Ping_f()
{
	if(cmd_source == src_command){
		Cmd_ForwardToServer();
		return;
	}
	SV_ClientPrintf("Client ping times:\n");
	s32 i = 0;
	for(client_t *client = svs.clients; i < svs.maxclients; i++, client++){
		if(!client->active) continue;
		f32 total = 0;
		for(s32 j = 0; j < NUM_PING_TIMES; j++)
			total += client->ping_times[j];
		total /= NUM_PING_TIMES;
		SV_ClientPrintf("%4i %s\n", (s32)(total * 1000), client->name);
	}
}


static void Host_Map_f()
{ // map <servername> command from console. Active clients are kicked off.
	if(Cmd_Argc() < 2){//no map name given
		if(cls.state == ca_dedicated){
			if(sv.active)
				Con_Printf("Current map: %s\n", sv.name);
			else
				Con_Printf("Server not active\n");
		}else if(cls.state == ca_connected)
			Con_Printf("Current map: %s( %s )\n",
					cl.levelname, cl.mapname);
		else
			Con_Printf("map <levelname>: start a new server\n");
		return;
	}
	if(cmd_source != src_command)
		return;
	cls.demonum = -1;// stop demo loop in case this fails
	CL_Disconnect();
	Host_ShutdownServer(false);
	key_dest = key_game;// remove console or menu
	SCR_BeginLoadingPlaque();
	svs.serverflags = 0;// haven't completed an episode yet
	s8 name[MAX_QPATH];
	q_strlcpy(name, Cmd_Argv(1), sizeof(name));
	// remove(any) trailing ".bsp" from mapname -- S.A.
	s8 *p = strstr(name, ".bsp");
	if(p && p[4] == '\0')
		*p = '\0';
	PR_SwitchQCVM(&sv.qcvm);
	SV_SpawnServer(name);
	PR_SwitchQCVM(NULL);
	if(!sv.active)
		return;
	if(cls.state != ca_dedicated){
		memset(cls.spawnparms, 0, MAX_MAPSTRING);
		for(s32 i = 2; i < Cmd_Argc(); i++){
			q_strlcat(cls.spawnparms, Cmd_Argv(i), MAX_MAPSTRING);
			q_strlcat(cls.spawnparms, " ", MAX_MAPSTRING);
		}
		Cmd_ExecuteString("connect local", src_command);
	}
}

void Host_Changelevel_f()
{ // Goes to a new map, taking all clients along
	s8 level[MAX_QPATH];
	if(Cmd_Argc() != 2){
	Con_Printf("changelevel <levelname> : continue game on a new level\n");
		return;
	}
	if(!sv.active || cls.demoplayback){
		Con_Printf("Only the server may changelevel\n");
		return;
	}
	strcpy(level, Cmd_Argv(1));
	key_dest = key_game;
	PR_SwitchQCVM(&sv.qcvm);
	SV_SaveSpawnparms();
	SV_SpawnServer(level);
	PR_SwitchQCVM(NULL);
}

void Host_Restart_f()
{ // Restarts the current server for a dead player
	s8 mapname[MAX_QPATH];
	if(cls.demoplayback || !sv.active) return;
	if(cmd_source != src_command) return;
	strcpy(mapname, sv.name); // must copy out, because it gets
	PR_SwitchQCVM(&sv.qcvm);
	SV_SpawnServer(mapname); // cleared in sv_spawnserver
	PR_SwitchQCVM(NULL);
}


void Host_Reconnect_f() // This command causes the client to wait for the signon
{ // messages again. This is sent just before a server changes levels
	SCR_BeginLoadingPlaque();
	cls.signon = 0; // need new connection messages
}


void Host_Connect_f()
{ // User command to connect to server
	s8 name[MAX_QPATH];
	cls.demonum = -1; // stop demo loop in case this fails
	if(cls.demoplayback){
		CL_StopPlayback();
		CL_Disconnect();
	}
	strcpy(name, Cmd_Argv(1));
	CL_EstablishConnection(name);
	Host_Reconnect_f();
}

void Host_SavegameComment(s8 *text)
{ // Writes a SAVEGAME_COMMENT_LENGTH character comment describing the current
	s8 kills[20];
	for(s32 i = 0; i < SAVEGAME_COMMENT_LENGTH; i++)
		text[i] = ' ';
	memcpy(text, cl.levelname, strlen(cl.levelname));
	sprintf(kills, "kills:%3i/%3i", cl.stats[STAT_MONSTERS],
			cl.stats[STAT_TOTALMONSTERS]);
	memcpy(text + 22, kills, strlen(kills));
	for(s32 i = 0; i < SAVEGAME_COMMENT_LENGTH; i++)
		if(text[i] == ' ') // convert space to _ to make stdio happy
			text[i] = '_';
	text[SAVEGAME_COMMENT_LENGTH] = '\0';
}

void Host_Savegame_f()
{
	s8 name[256];
	s8 comment[SAVEGAME_COMMENT_LENGTH + 1];
	if(cmd_source != src_command) return;
if(!sv.active){ Con_Printf("Not playing a local game.\n"); return; }
if(cl.intermission){ Con_Printf("Can't save in intermission.\n"); return; }
if(svs.maxclients != 1){ Con_Printf("Can't save multiplayer games.\n"); return;}
if(Cmd_Argc() != 2){ Con_Printf("save <savename> : save a game\n"); return; }
	if(strstr(Cmd_Argv(1), "..")){
		Con_Printf("Relative pathnames are not allowed.\n");
		return;
	}
	for(s32 i = 0; i < svs.maxclients; i++)
		if(svs.clients[i].active&&(svs.clients[i].edict->v.health<=0)){
			Con_Printf("Can't savegame with a dead player\n");
			return;
		} 
	snprintf(name, sizeof(name), "%s/%s", com_gamedir, Cmd_Argv(1));
	COM_AddExtension(name, ".sav", sizeof(name));
	Con_Printf("Saving game to %s...\n", name);
	PR_SwitchQCVM(&sv.qcvm);
	FILE *f = fopen(name, "w");
	if(!f){ Con_Printf("ERROR: couldn't open.\n"); return; }
	fprintf(f, "%i\n", SAVEGAME_VERSION);
	Host_SavegameComment(comment);
	fprintf(f, "%s\n", comment);
	for(s32 i = 0; i < NUM_SPAWN_PARMS; i++)
		fprintf(f, "%f\n", svs.clients->spawn_parms[i]);
	fprintf(f, "%d\n", current_skill);
	fprintf(f, "%s\n", sv.name);
	fprintf(f, "%f\n", qcvm->time);
	for(s32 i = 0; i < MAX_LIGHTSTYLES; i++){ // write the light styles
		if(sv.lightstyles[i]) fprintf(f, "%s\n", sv.lightstyles[i]);
		else fprintf(f, "m\n");
	}
	ED_WriteGlobals(f);
	for(s32 i = 0; i < qcvm->num_edicts; i++){
		ED_Write(f, EDICT_NUM(i));
		fflush(f);
	}
	fclose(f);
	PR_SwitchQCVM(NULL);
	Con_Printf("done.\n");
}

void Host_Loadgame_f()
{
	s8 name[MAX_OSPATH+2];
	s8 mapname[MAX_QPATH];
	s8 str[32768];
	f32 spawn_parms[NUM_SPAWN_PARMS];
	if(cmd_source != src_command) return;
if(Cmd_Argc() != 2){ Con_Printf("load <savename> : load a game\n"); return; }
	cls.demonum = -1; // stop demo loop in case this fails
	snprintf(name, sizeof(name), "%s/%s", com_gamedir, Cmd_Argv(1));
	COM_AddExtension(name, ".sav", sizeof(name));
	// we can't call SCR_BeginLoadingPlaque because too much stack space has
	// been used. The menu calls it before stuffing loadgame command
	// SCR_BeginLoadingPlaque();
	Con_Printf("Loading game from %s...\n", name);
	FILE *f = fopen(name, "r");
	if(!f){ Con_Printf("ERROR: couldn't open.\n"); return; }
	s32 version;
	fscanf(f, "%i\n", &version);
	if(version != SAVEGAME_VERSION){
		fclose(f);
      Con_Printf("Savegame is version %i, not %i\n", version, SAVEGAME_VERSION);
		return;
	}
	fscanf(f, "%s\n", str);
	for(s32 i = 0; i < NUM_SPAWN_PARMS; i++)
		fscanf(f, "%f\n", &spawn_parms[i]);
// this silliness is so we can load 1.06 save files, which have f32 skill values
	f32 tfloat;
	fscanf(f, "%f\n", &tfloat);
	current_skill = (s32)(tfloat + 0.1);
	Cvar_SetValue("skill", (f32)current_skill);
	fscanf(f, "%s\n", mapname);
	f32 time;
	fscanf(f, "%f\n", &time);
	CL_Disconnect_f();
	PR_SwitchQCVM(&sv.qcvm);
	SV_SpawnServer(mapname);
	if(!sv.active){
		PR_SwitchQCVM(NULL);
		Con_Printf("Couldn't load map\n");
		return; 
	}
	sv.paused = 1; // pause until all clients connect
	sv.loadgame = 1;
	refresh_palette = 10;
	s32 i = 0;
	for(; i < MAX_LIGHTSTYLES; i++){ // load the light styles
		fscanf(f, "%s\n", str);
		s8 *dst = Hunk_Alloc(strlen(str) + 1);
		strcpy(dst, str);
		sv.lightstyles[i] = dst;
	}
	// load the edicts out of the savegame file
	s32 entnum = -1; // -1 is the globals
	while(!feof(f)){
		for(i = 0; i < (s32)sizeof(str) - 1; i++){
			s32 r = fgetc(f);
			if(r == EOF || !r) break;
			str[i] = r;
			if(r == '}'){ i++; break; }
		}
		if(i == sizeof(str) - 1) Sys_Error("Loadgame buffer overflow");
		str[i] = 0;
		const s8 *start = str;
		start = COM_Parse(str);
		if(!com_token[0]) break; // end of file
		if(strcmp(com_token,"{"))Sys_Error("First token isn't a brace");
		if(entnum == -1){ ED_ParseGlobals(start);
		} else { // parse an edict
			edict_t *ent = EDICT_NUM(entnum);
			if(entnum < qcvm->num_edicts){ ED_ClearEdict(ent); }
			else {
				memset(ent, 0, qcvm->edict_size);
				ent->baseline.scale = ENTSCALE_DEFAULT;
			}
			start = ED_ParseEdict(start, ent);
			if(!ent->free) // link it into the bsp tree
				SV_LinkEdict(ent, false);
		}
		entnum++;
	}
	// Free edicts allocated during map loading but no longer used after restoring saved game state
	// Note: we use ED_ClearEdict instead of ED_Free to avoid placing entities >= num_edicts in the free list
	// This is different from QuakeSpasm, which doesn't use a free list
	for(i = entnum; i < qcvm->num_edicts; i++)
		ED_ClearEdict(EDICT_NUM(i));
	qcvm->num_edicts = entnum;
	qcvm->time = time;
	fclose(f);
	for(s32 i = 0; i < NUM_SPAWN_PARMS; i++)
		svs.clients->spawn_parms[i] = spawn_parms[i];
	PR_SwitchQCVM(NULL);
	if(cls.state != ca_dedicated){
		CL_EstablishConnection("local");
		Host_Reconnect_f();
	}
}

static void Host_Name_f()
{
	s8 newName[32];
	if(Cmd_Argc() == 1){
		Con_Printf("\"name\" is \"%s\"\n", cl_name.string);
		return;
	}
	if(Cmd_Argc() == 2) q_strlcpy(newName, Cmd_Argv(1), sizeof(newName));
	else q_strlcpy(newName, Cmd_Args(), sizeof(newName));
	newName[15] = 0; // client_t structure actually says name[32].
	if(cmd_source == src_command){
		if(Q_strcmp(cl_name.string, newName) == 0) return;
		Cvar_Set("_cl_name", newName);
		if(cls.state == ca_connected) Cmd_ForwardToServer();
		return;
	}
	if(host_client->name[0] && strcmp(host_client->name, "unconnected")
		&& Q_strcmp(host_client->name, newName) != 0)
		Con_Printf("%s renamed to %s\n", host_client->name, newName);
	Q_strcpy(host_client->name, newName);
	host_client->edict->v.netname = PR_SetEngineString(host_client->name);
	// send notification to all clients
	MSG_WriteByte(&sv.reliable_datagram, svc_updatename);
	MSG_WriteByte(&sv.reliable_datagram, host_client - svs.clients);
	MSG_WriteString(&sv.reliable_datagram, host_client->name);
}

void Host_Version_f()
{
	Con_Printf("Version %4.2f\n", VERSION);
	Con_Printf("Exe: " __TIME__ " " __DATE__ "\n");
}

void Host_Say(bool teamonly)
{
	s8 text[64];
	bool fromServer = 0;
	if(cmd_source == src_command){
		if(cls.state == ca_dedicated){ fromServer = 1; teamonly = 0; }
		else { Cmd_ForwardToServer(); return; }
	}
	if(Cmd_Argc() < 2) return;
	client_t *save = host_client;
	s8 *p = Cmd_Args();
	if(*p == '"'){ // remove quotes if present
		p++;
		p[Q_strlen(p) - 1] = 0;
	}
	if(!fromServer)sprintf(text,"%c%s: ",1,save->name);//turn on color set 1
	else sprintf(text, "%c<%s> ", 1, hostname.string);
	s32 j = sizeof(text)-2-Q_strlen(text); // -2 for /n and null terminator
	if(Q_strlen(p) > j) p[j] = 0;
	strcat(text, p);
	strcat(text, "\n");
	j = 0;
	for(client_t *client = svs.clients; j < svs.maxclients; j++, client++){
		if(!client || !client->active || !client->spawned) continue;
		if(teamplay.value && teamonly
			&& client->edict->v.team != save->edict->v.team)
			continue;
		host_client = client;
		SV_ClientPrintf("%s", text);
	}
	host_client = save;
	Sys_Printf("%s", &text[1]);
}

void Host_Say_f(){ Host_Say(0); }
void Host_Say_Team_f(){ Host_Say(1); }
void Host_Tell_f()
{
	s8 text[64];
	if(cmd_source == src_command){ Cmd_ForwardToServer(); return; }
	if(Cmd_Argc() < 3) return;
	Q_strcpy(text, host_client->name);
	Q_strcat(text, ": ");
	s8 *p = Cmd_Args();
	if(*p == '"'){ // remove quotes if present
		p++;
		p[Q_strlen(p) - 1] = 0;
	}
	// check length & truncate if necessary
	s32 j = sizeof(text)-2-Q_strlen(text); // -2 for /n and null terminator
	if(Q_strlen(p) > j)
		p[j] = 0;
	strcat(text, p);
	strcat(text, "\n");
	client_t *save = host_client;
	j = 0;
	for(client_t *client = svs.clients; j < svs.maxclients; j++, client++){
		if(!client->active || !client->spawned)
			continue;
		if(q_strcasecmp(client->name, Cmd_Argv(1)))
			continue;
		host_client = client;
		SV_ClientPrintf("%s", text);
		break;
	}
	host_client = save;
}

void Host_Color_f()
{
	if(Cmd_Argc() == 1){
		Con_Printf("\"color\" is \"%i %i\"\n",
			((s32)cl_color.value)>>4, ((s32)cl_color.value) & 0x0f);
		Con_Printf("color <0-13> [0-13]\n");
		return;
	}
	s32 top, bottom;
	if(Cmd_Argc() == 2) top = bottom = atoi(Cmd_Argv(1));
	else {
		top = atoi(Cmd_Argv(1));
		bottom = atoi(Cmd_Argv(2));
	}
	top &= 15;
	if(top > 13) top = 13;
	bottom &= 15;
	if(bottom > 13) bottom = 13;
	s32 playercolor = top * 16 + bottom;
	if(cmd_source == src_command){
		Cvar_SetValue("_cl_color", playercolor);
		if(cls.state == ca_connected) Cmd_ForwardToServer();
		return;
	}
	host_client->colors = playercolor;
	host_client->edict->v.team = bottom + 1;
	// send notification to all clients
	MSG_WriteByte(&sv.reliable_datagram, svc_updatecolors);
	MSG_WriteByte(&sv.reliable_datagram, host_client - svs.clients);
	MSG_WriteByte(&sv.reliable_datagram, host_client->colors);
}

void Host_Kill_f()
{
	if(cmd_source == src_command){ Cmd_ForwardToServer(); return; }
	if(sv_player->v.health <= 0){
		SV_ClientPrintf("Can't suicide -- allready dead!\n");
		return;
	}
	pr_global_struct->time = qcvm->time;
	pr_global_struct->self = EDICT_TO_PROG(sv_player);
	PR_ExecuteProgram(pr_global_struct->ClientKill);
}

void Host_Pause_f()
{
	if(cmd_source == src_command){ Cmd_ForwardToServer(); return; }
	if(!pausable.value) SV_ClientPrintf("Pause not allowed.\n");
	else {
		sv.paused ^= 1;
		if(sv.paused){
			SV_BroadcastPrintf("%s paused the game\n",
					PR_GetString(sv_player->v.netname));
		} else {
			SV_BroadcastPrintf("%s unpaused the game\n",
					PR_GetString(sv_player->v.netname));
		}
		// send notification to all clients
		MSG_WriteByte(&sv.reliable_datagram, svc_setpause);
		MSG_WriteByte(&sv.reliable_datagram, sv.paused);
	}
}

void Host_PreSpawn_f()
{
	if(cmd_source == src_command){
		Con_Printf("prespawn is not valid from the console\n");
		return;
	}
	if(host_client->spawned){
		Con_Printf("prespawn not valid -- allready spawned\n");
		return;
	}
	host_client->sendsignon = PRESPAWN_SIGNONBUFS;
	host_client->signonidx = 0;
}

static void Host_Spawn_f()
{
	if(cmd_source == src_command) {
		Con_Printf("spawn is not valid from the console\n");
		return;
	}
	if(host_client->spawned) {
		Con_Printf("Spawn not valid -- already spawned\n");
		return;
	}
	// run the entrance script
	if(sv.loadgame){// loaded games are fully inited already
	    // if this is the last client to be connected, unpause
		sv.paused = false;
	} else {
		edict_t *ent = host_client->edict; // set up the edict
		memset(&ent->v, 0, qcvm->progs->entityfields * 4);
		ent->v.colormap = NUM_FOR_EDICT(ent);
		ent->v.team = (host_client->colors & 15) + 1;
		ent->v.netname = PR_SetEngineString(host_client->name);
		// copy spawn parms out of the client_t
		for(s32 i = 0; i < NUM_SPAWN_PARMS; i++)
		    (&pr_global_struct->parm1)[i]=host_client->spawn_parms[i];
		// call the spawn function
		pr_global_struct->time = qcvm->time;
		pr_global_struct->self = EDICT_TO_PROG(sv_player);
		PR_ExecuteProgram(pr_global_struct->ClientConnect);
		if((Sys_DoubleTime() - host_client->netconnection->connecttime)
				<= qcvm->time)
			Sys_Printf("%s entered the game\n", host_client->name);
		PR_ExecuteProgram(pr_global_struct->PutClientInServer);
	}
	// send all current names, colors, and frag counts
	SZ_Clear(&host_client->message);
	MSG_WriteByte(&host_client->message, svc_time); // send time of update
	MSG_WriteFloat(&host_client->message, qcvm->time);
	s32 i;
	client_t *client;
	for(i = 0, client = svs.clients; i < svs.maxclients; i++, client++){
		MSG_WriteByte(&host_client->message, svc_updatename);
		MSG_WriteByte(&host_client->message, i);
		MSG_WriteString(&host_client->message, client->name);
		MSG_WriteByte(&host_client->message, svc_updatefrags);
		MSG_WriteByte(&host_client->message, i);
		MSG_WriteShort(&host_client->message, client->old_frags);
		MSG_WriteByte(&host_client->message, svc_updatecolors);
		MSG_WriteByte(&host_client->message, i);
		MSG_WriteByte(&host_client->message, client->colors);
	}
	// send all current light styles
	for(i = 0; i < MAX_LIGHTSTYLES; i++){
		MSG_WriteByte(&host_client->message, svc_lightstyle);
		MSG_WriteByte(&host_client->message, (s8)i);
		MSG_WriteString(&host_client->message, sv.lightstyles[i]);
	}
	// send some stats
	MSG_WriteByte(&host_client->message, svc_updatestat);
	MSG_WriteByte(&host_client->message, STAT_TOTALSECRETS);
	MSG_WriteLong(&host_client->message, pr_global_struct->total_secrets);
	MSG_WriteByte(&host_client->message, svc_updatestat);
	MSG_WriteByte(&host_client->message, STAT_TOTALMONSTERS);
	MSG_WriteLong(&host_client->message, pr_global_struct->total_monsters);
	MSG_WriteByte(&host_client->message, svc_updatestat);
	MSG_WriteByte(&host_client->message, STAT_SECRETS);
	MSG_WriteLong(&host_client->message, pr_global_struct->found_secrets);
	MSG_WriteByte(&host_client->message, svc_updatestat);
	MSG_WriteByte(&host_client->message, STAT_MONSTERS);
	MSG_WriteLong(&host_client->message, pr_global_struct->killed_monsters);
	// send a fixangle
	// Never send a roll angle, because savegames can catch the server
	// in a state where it is expecting the client to correct the angle
	// and it won't happen if the game was just loaded, so you wind up
	// with a permanent head tilt
	edict_t *ent = EDICT_NUM(1 + (host_client - svs.clients));
	MSG_WriteByte(&host_client->message, svc_setangle);
	for(i = 0; i < 2; i++)
		if(sv.loadgame)
			MSG_WriteAngle(&host_client->message,
					ent->v.v_angle[i], sv.protocolflags );
		else
			MSG_WriteAngle(&host_client->message,
					ent->v.angles[i], sv.protocolflags );
	MSG_WriteAngle(&host_client->message, 0, sv.protocolflags );
	SV_WriteClientdataToMessage(sv_player, &host_client->message);
	MSG_WriteByte(&host_client->message, svc_signonnum);
	MSG_WriteByte(&host_client->message, 3);
	host_client->sendsignon = PRESPAWN_FLUSH;
}

void Host_Begin_f()
{
	if(cmd_source == src_command){
		Con_Printf("begin is not valid from the console\n");
		return;
	}
	host_client->spawned = 1;
}

void Host_Kick_f() // Kicks a user off of the server
{
	if(cmd_source == src_command){
		if(!sv.active){
			Cmd_ForwardToServer();
			return;
		}
	} else if(pr_global_struct->deathmatch) return;
	client_t *save = host_client;
	s32 i;
	bool byNumber = 0;
	if(Cmd_Argc() > 2 && Q_strcmp(Cmd_Argv(1), "#") == 0){
		i = Q_atof(Cmd_Argv(2)) - 1;
		if(i < 0 || i >= svs.maxclients) return;
		if(!svs.clients[i].active) return;
		host_client = &svs.clients[i];
		byNumber = 1;
	} else {
		for(i = 0, host_client = svs.clients; i < svs.maxclients;
							i++, host_client++){
			if(!host_client->active) continue;
			if(q_strcasecmp(host_client->name,Cmd_Argv(1))==0)break;
		}
	}
	if(i < svs.maxclients){
		const s8 *who;
		if(cmd_source == src_command)
		      who=cls.state==ca_dedicated?(s8*)"Console":cl_name.string;
		else who = save->name;
		if(host_client == save) return; // can't kick yourself!
		const s8 *message = NULL;
		if(Cmd_Argc() > 2){
			message = COM_Parse(Cmd_Args());
			if(byNumber){
				message++; // skip the #
				while(*message == ' ') // skip white space
					message++;
				message += Q_strlen(Cmd_Argv(2)); // skip number
			}
			while(*message && *message == ' ') message++;
		}
		if(message) SV_ClientPrintf("Kicked by %s: %s\n", who, message);
		else SV_ClientPrintf("Kicked by %s\n", who);
		SV_DropClient(0);
	}
	host_client = save;
}

void Host_Give_f()
{
	eval_t *val;
	if(cmd_source == src_command){ Cmd_ForwardToServer(); return; }
	if(pr_global_struct->deathmatch) return;
	s8 *t = Cmd_Argv(1);
	s32 v = atoi(Cmd_Argv(2));
	switch(t[0]){
	case '0': case '1': case '2': case '3': case '4':
	case '5': case '6': case '7': case '8': case '9':
	if(hipnotic){
	  if(t[0]=='6')
		sv_player->v.items=(s32)sv_player-> v.items|(t[1]=='a'?
				HIT_PROXIMITY_GUN :IT_GRENADE_LAUNCHER);
	  else if(t[0]=='9')
		sv_player->v.items=(s32)sv_player-> v.items|HIT_LASER_CANNON;
	  else if(t[0]=='0')
		sv_player->v.items=(s32)sv_player-> v.items|HIT_MJOLNIR;
	  else if(t[0]>='2')
	   sv_player->v.items=(s32)sv_player-> v.items|(IT_SHOTGUN<<(t[0]-'2'));
	} else if(t[0]>='2')
	sv_player->v.items=(s32)sv_player-> v.items|(IT_SHOTGUN<<(t[0]-'2'));
	break;
	case 's':
		if(rogue){
			val = GetEdictFieldValueByName(sv_player,"ammo_shells1");
			if(val) val->_float = v; }
		sv_player->v.ammo_shells = v; break;
	case 'n':
		if(rogue){
			val = GetEdictFieldValueByName(sv_player,"ammo_nails1");
			if(val){ val->_float = v;
				if(sv_player->v.weapon <= IT_LIGHTNING)
					sv_player->v.ammo_nails = v; }
		} else sv_player->v.ammo_nails = v; break;
	case 'l':
		if(rogue){
		      val=GetEdictFieldValueByName(sv_player,"ammo_lava_nails");
			if(val){ val->_float = v;
				if(sv_player->v.weapon > IT_LIGHTNING)
					sv_player->v.ammo_nails = v; } } break;
	case 'r':
		if(rogue){
		      val = GetEdictFieldValueByName(sv_player,"ammo_rockets1");
			if(val){ val->_float = v;
				if(sv_player->v.weapon <= IT_LIGHTNING)
					sv_player->v.ammo_rockets = v; }
		} else sv_player->v.ammo_rockets = v; break;
	case 'm':
	       if(rogue){
		   val=GetEdictFieldValueByName(sv_player,"ammo_multi_rockets");
			if(val){ val->_float = v;
				if(sv_player->v.weapon > IT_LIGHTNING)
					sv_player->v.ammo_rockets = v; }} break;
	case 'h': sv_player->v.health = v; break;
	case 'c':
		if(rogue){
			val = GetEdictFieldValueByName(sv_player,"ammo_cells1");
			if(val){ val->_float = v;
				if(sv_player->v.weapon <= IT_LIGHTNING)
					sv_player->v.ammo_cells = v; }
		} else sv_player->v.ammo_cells = v; break;
	case 'p':
		if(rogue){
			val = GetEdictFieldValueByName(sv_player,"ammo_plasma");
			if(val){ val->_float = v;
				if(sv_player->v.weapon > IT_LIGHTNING)
					sv_player->v.ammo_cells = v;
			} } break;
	}
}

edict_t *FindViewthing()
{
	PR_SwitchQCVM(&sv.qcvm);
	for(s32 i = 0; i < qcvm->num_edicts; i++){
		edict_t *e = EDICT_NUM(i);
		if(!strcmp(PR_GetString(e->v.classname), "viewthing"))
			return e;
	}
	Con_Printf("No viewthing on map\n");
	PR_SwitchQCVM(NULL);
	return NULL;
}

void Host_Viewmodel_f()
{
	edict_t *e = FindViewthing();
	if(!e) return;
	model_t *m = Mod_ForName(Cmd_Argv(1), 0);
	if(!m){ Con_Printf("Can't load %s\n", Cmd_Argv(1)); return; }
	PR_SwitchQCVM(&sv.qcvm);
	e->v.frame = 0;
	cl.model_precache[(s32)e->v.modelindex] = m;
	PR_SwitchQCVM(NULL);
}

void Host_Viewframe_f()
{
	edict_t *e = FindViewthing();
	if(!e) return;
	model_t *m = cl.model_precache[(s32)e->v.modelindex];
	s32 f = atoi(Cmd_Argv(1));
	if(f >= m->numframes) f = m->numframes - 1;
	e->v.frame = f;
}

void PrintFrameName(model_t *m, s32 frame)
{
	aliashdr_t *hdr = (aliashdr_t *) Mod_Extradata(m);
	if(!hdr) return;
	maliasframedesc_t *pframedesc = &hdr->frames[frame];
	Con_Printf("frame %i: %s\n", frame, pframedesc->name);
}

void Host_Viewnext_f()
{
	edict_t *e = FindViewthing();
	if(!e) return;
	model_t *m = cl.model_precache[(s32)e->v.modelindex];
	e->v.frame = e->v.frame + 1;
	if(e->v.frame >= m->numframes) e->v.frame = m->numframes - 1;
	PrintFrameName(m, e->v.frame);
}

void Host_Viewprev_f()
{
	edict_t *e = FindViewthing();
	if(!e) return;
	model_t *m = cl.model_precache[(s32)e->v.modelindex];
	e->v.frame = e->v.frame - 1;
	if(e->v.frame < 0) e->v.frame = 0;
	PrintFrameName(m, e->v.frame);
}

void Host_Startdemos_f()
{
	if(cls.state==ca_dedicated&&!sv.active){
		Cbuf_AddText("map start\n");
		return; }
	s32 c = Cmd_Argc() - 1;
	if(c > MAX_DEMOS){
		Con_Printf("Max %i demos in demoloop\n", MAX_DEMOS);
		c = MAX_DEMOS;
	}
	Con_Printf("%i demo(s) in loop\n", c);
	for(s32 i = 1; i < c + 1; i++)
		strncpy(cls.demos[i - 1], Cmd_Argv(i), sizeof(cls.demos[0])-1);
	if(!sv.active && cls.demonum != -1 && !cls.demoplayback){
		cls.demonum = 0;
		CL_NextDemo();
	} else cls.demonum = -1;
}

void Host_Demos_f()
{ // Return to looping demos
	if(cls.state == ca_dedicated) return;
	if(cls.demonum == -1) cls.demonum = 1;
	CL_Disconnect_f();
	CL_NextDemo();
}

void Host_Stopdemo_f()
{ // Return to looping demos
	if(cls.state == ca_dedicated || !cls.demoplayback) return;
	CL_StopPlayback();
	CL_Disconnect();
}

time_t Mod_GetMapDate(const s8 *map)
{
	s8 path[MAX_QPATH];
	if((size_t) q_snprintf(path, sizeof(path), "%s/maps/%s.bsp",
				com_gamedir, map) >= sizeof(path))
		return 0;
	SDL_PathInfo info;
	SDL_GetPathInfo(path, &info);
	u64 ns = info.modify_time;
	time_t seconds = ns / 1000000000ULL;
	return seconds;
}

void FileList_Add(const char *name, const char *desc, filelist_item_t **list)
{
	filelist_item_t *item,*cursor,*prev;
	for(item = *list; item; item = item->next) // ignore duplicate
		if(!Q_strcmp(name, item->name)) return;
	item = (filelist_item_t *) Z_Malloc(sizeof(filelist_item_t));
	q_strlcpy(item->name, name, sizeof(item->name));
	if(desc) q_strlcpy(item->desc, desc, sizeof(item->desc));
	// insert each entry in alphabetical order
	if(*list == NULL || q_strcasecmp(item->name, (*list)->name) < 0){
		item->next = *list; //insert at front
		*list = item;
	} else { //insert later
		prev = *list;
		cursor = (*list)->next;
		while(cursor && (q_strcasecmp(item->name, cursor->name) > 0)){
			prev = cursor;
			cursor = cursor->next;
		}
		item->next = prev->next;
		prev->next = item;
	}
}

void FileList_AddMap(const char *name, const char *desc, filelist_item_t **list)
{
	filelist_item_t *item,*cursor,*prev;
	for(item = *list; item; item = item->next) // ignore duplicate
		if(!Q_strcmp(name, item->name)) return;
	item = (filelist_item_t *) Z_Malloc(sizeof(filelist_item_t));
	q_strlcpy(item->name, name, sizeof(item->name));
	if(desc) q_strlcpy(item->desc, desc, sizeof(item->desc));
	item->data1 = Mod_CountMonsters(name);
	item->data2 = Mod_CountSecrets(name);
	item->date = Mod_GetMapDate(name);
	// insert each entry in alphabetical order
	if(*list == NULL || q_strcasecmp(item->name, (*list)->name) < 0){
		item->next = *list; //insert at front
		*list = item;
	} else { //insert later
		prev = *list;
		cursor = (*list)->next;
		while(cursor && (q_strcasecmp(item->name, cursor->name) > 0)){
			prev = cursor;
			cursor = cursor->next;
		}
		item->next = prev->next;
		prev->next = item;
	}
}

void ExtraMaps_Add(const s8 *name, const s8 *game)
{
	s8 buf[128];
	filelist_item_t **list = &extralevels;
	if(!Mod_LoadMapDescription(buf, sizeof(buf), name))
		return;
	if(Q_strncmp(COM_SkipPath(game), "id1", 4))
		list = &extralevels_mod;
	FileList_AddMap(name, buf, list);
}

static void ExtraMaps_Init_SearchDir(searchpath_t *search)
{
	s8 filestring[MAX_OSPATH];
	s8 mapname[32];
#ifdef _WIN32
	WIN32_FIND_DATA fdat;
	HANDLE fhnd;
	q_snprintf(filestring, sizeof(filestring), "%s/maps/*.bsp",
							search->filename);
	fhnd = FindFirstFile(filestring, &fdat);
	if(fhnd == INVALID_HANDLE_VALUE) continue;
	do {
		COM_StripExtension(fdat.cFileName, mapname, sizeof(mapname));
		if(maxlevelnamelen < Q_strlen(mapname))
			maxlevelnamelen = Q_strlen(mapname);
		ExtraMaps_Add(mapname, search->filename);
	} while(FindNextFile(fhnd, &fdat));
	FindClose(fhnd);
#else
	DIR *dir_p;
	struct dirent *dir_t;
	q_snprintf(filestring, sizeof(filestring), "%s/maps/",search->filename);
	dir_p = opendir(filestring);
	if(dir_p == NULL) return;;
	while((dir_t = readdir(dir_p)) != NULL){
		if(q_strcasecmp(COM_FileGetExtension(dir_t->d_name),"bsp") != 0)
			continue;
		COM_StripExtension(dir_t->d_name, mapname, sizeof(mapname));
		if(maxlevelnamelen < Q_strlen(mapname))
			maxlevelnamelen = Q_strlen(mapname);
		ExtraMaps_Add(mapname, search->filename);
	}
	closedir(dir_p);
#endif
}

static void ExtraMaps_Init_SearchPak(searchpath_t *search)
{
	s8 mapname[32];
	s8 ignorepakdir[32];
	// we don't want to list the maps in id1 pakfiles,
	// because these are not "add-on" levels
	q_snprintf(ignorepakdir, sizeof(ignorepakdir), "/%s/", GAMENAME);
	if(strstr(search->pack->filename, ignorepakdir))
		return;
	s32 i = 0;
	for(pack_t *pak = search->pack; i < pak->numfiles; i++){
		if(!strcmp(COM_FileGetExtension(pak->files[i].name), "bsp") &&
		    pak->files[i].filelen > 32*1024){
			// don't list files under 32k(ammo boxes etc)
			COM_StripExtension(pak->files[i].name + 5, mapname,
							sizeof(mapname));
			ExtraMaps_Add(mapname, com_gamedir);
		} 
	}
}

void ExtraMaps_Init()
{
	maxlevelnamelen = 0;
	for(searchpath_t *search=com_searchpaths; search; search=search->next){
		if(*search->filename){
			ExtraMaps_Init_SearchDir(search);
		} else {
			ExtraMaps_Init_SearchPak(search);
		}
	}
}

static void FileList_Clear(filelist_item_t **list)
{
	filelist_item_t *blah;
	while(*list){
		blah = (*list)->next;
		Z_Free(*list);
		*list = blah;
	}
}

static void ExtraMaps_Clear()
{
	FileList_Clear(&extralevels);
	FileList_Clear(&extralevels_mod);
}

void ExtraMaps_NewGame()
{
	ExtraMaps_Clear();
	ExtraMaps_Init();
}

static const s8 *RightPad(const s8 *str, size_t minlen, s8 c)
{
	static s8 buf[1024];
	size_t len = strlen(str);
	minlen = q_min(minlen, sizeof(buf) - 1);
	if(len >= minlen)
		return str;
	memcpy(buf, str, len);
	for(; len < minlen; len++)
		buf[len] = c;
	buf[len] = '\0';
	return buf;
}

void Host_Maps_f()
{
	s32 i;
	s32 tot = 0;
	s8 padchar = '.' | 0x80;
	filelist_item_t *level;
	Con_Printf("%s\n", RightPad("id1", 32, '-'|0x80));
	for(level = extralevels, i = 0; level; level = level->next, i++){
		Con_Printf("   %s%c%s\n", RightPad(level->name,
			maxlevelnamelen, padchar), padchar, level->desc);
		++tot;
	}
	if(!Q_strncmp("id1", COM_SkipPath(com_gamedir), 4))
		goto host_maps_f_fin;
	Con_Printf("%s\n", RightPad(COM_SkipPath(com_gamedir), 32, '-'|0x80));
	for(level = extralevels_mod, i = 0; level; level = level->next, i++){
		Con_Printf("   %s%c%s\n", RightPad(level->name,
			maxlevelnamelen, padchar), padchar, level->desc);
		++tot;
	}
host_maps_f_fin:
	if(tot>1) Con_Printf("%i maps\n", i);
	else if(tot) Con_Printf("1 map\n");
	else Con_Printf("no maps found\n");
}

void Modlist_Add(const char *name, const char *desc){
	FileList_Add(name, desc, &modlist);
}

static const char *Modlist_KnownDescription(const char *modname)
{
	for(u32 i = 0; i < sizeof(knownmods) / sizeof(knownmods[0]); i++){
		if(!q_strcasecmp(modname, knownmods[i][0]))
			return knownmods[i][1];
	}
	return NULL;
}

char *Modlist_ReadDescription(const char *mod_path)
{
	static s8 desc[128];
	FILE *f;
	s8 path[MAX_OSPATH];
	q_snprintf(path, sizeof(path), "%s/descript.ion", mod_path);
	f = fopen(path, "rb");
	if(!f)
		return NULL;
	size_t len = fread(desc, 1, sizeof(desc) - 1, f);
	fclose(f);
	if(len == 0)
		return NULL;
	desc[len] = '\0';
	// strip trailing newline(s)
	while(len > 0 && (desc[len-1] == '\n' || desc[len-1] == '\r')){
		desc[len-1] = '\0';
		len--;
	}
	return desc;
}

#ifdef _WIN32
void Modlist_Init()
{
	WIN32_FIND_DATA fdat;
	HANDLE fhnd;
	DWORD attribs;
	s8 dir_string[MAX_OSPATH], mod_string[MAX_OSPATH];
	q_snprintf(dir_string, sizeof(dir_string), "%s/*", com_basedir);
	fhnd = FindFirstFile(dir_string, &fdat);
	maxmodnamelen = 0;
	if(fhnd == INVALID_HANDLE_VALUE) return;
	do {
		if(!strcmp(fdat.cFileName,".") || !strcmp(fdat.cFileName,".."))
			continue;
		q_snprintf(mod_string, sizeof(mod_string), "%s/%s",
					com_basedir, fdat.cFileName);
		attribs = GetFileAttributes(mod_string);
		if(attribs != INVALID_FILE_ATTRIBUTES &&
				(attribs & FILE_ATTRIBUTE_DIRECTORY)){
			// don't bother testing for pak files / progs.dat
			s8 *file_desc = Modlist_ReadDescription(mod_string);
			s8 *desc = file_desc ? file_desc :
				Modlist_KnownDescription(fdat.cFileName);
			if(maxmodnamelen < Q_strlen(fdat.cFileName))
				maxmodnamelen = Q_strlen(fdat.cFileName);
			Modlist_Add(fdat.cFileName, desc);

		}
	} while(FindNextFile(fhnd, &fdat));
	FindClose(fhnd);
}
#else
void Modlist_Init()
{
	DIR *dir_p, *mod_dir_p;
	struct dirent *dir_t;
	s8 dir_string[MAX_OSPATH], mod_string[MAX_OSPATH];
	q_snprintf(dir_string, sizeof(dir_string), "%s/", com_basedir);
	dir_p = opendir(dir_string);
	maxmodnamelen = 0;
	if(dir_p == NULL) return;
	while((dir_t = readdir(dir_p)) != NULL){
		if(!strcmp(dir_t->d_name, ".") || !strcmp(dir_t->d_name, ".."))
			continue;
		if(!q_strcasecmp(COM_FileGetExtension(dir_t->d_name), "app"))
			continue; // skip .app bundles on macOS
		q_snprintf(mod_string, sizeof(mod_string), "%s%s/",
				dir_string, dir_t->d_name);
		mod_dir_p = opendir(mod_string);
		if(mod_dir_p == NULL)
			continue;
		// don't bother testing for pak files / progs.dat
		char *file_desc = Modlist_ReadDescription(mod_string);
		const char *desc = file_desc ? file_desc
				: Modlist_KnownDescription(dir_t->d_name);
		if(maxmodnamelen < Q_strlen(dir_t->d_name))
			maxmodnamelen = Q_strlen(dir_t->d_name);
		Modlist_Add(dir_t->d_name, desc);
		closedir(mod_dir_p);
	}
	closedir(dir_p);
}
#endif

void Host_Mods_f()
{//list all potential mod directories(contain either a pak file or a progs.dat)
	s32 i;
	filelist_item_t *mod;
	s8 padchar = '.' | 0x80;
	for(mod = modlist, i=0; mod; mod = mod->next, i++)
		Con_Printf("   %s%c%s\n", RightPad(mod->name,
			maxmodnamelen, padchar), padchar, mod->desc);
	if(i == 1) Con_Printf("1 mod\n");
	if(i) Con_Printf("%i mods\n", i);
	else Con_Printf("no mods found\n");
}

void Cmd_Resurrect_f()
{
	PR_SwitchQCVM(&sv.qcvm);
	if(!sv.active || cls.demoplayback){
		Con_Printf("You must be in a game to resurrect.\n");
		goto resurrect_ret;
	}
	if(pr_global_struct->deathmatch || sv.qcvm.max_edicts <= 1)
		goto resurrect_ret;
	edict_t *ent = EDICT_NUM(1);
	if(!ent || ent->free){
		Con_Printf("Error: PLAYER edict missing/freed.\n");
		goto resurrect_ret;
	}
	if(ent->v.health > 0){
		Con_Printf("You are still alive!\n");
		goto resurrect_ret;
	}
	if(ent->v.nextthink != -1){
		Con_Printf("Not fully dead yet.\n");
		goto resurrect_ret;
	}
	// 2. Chained assignments compress the Physics & Camera resets
	if(ent->v.health < -40){ // gibbed
		ent->v.origin[2] += 24;
		ent->v.flags = (u64)ent->v.flags & ~FL_ONGROUND;
		ent->v.mins[2] = -24;
		ent->v.maxs[2] = 32;
	}
	ent->v.health = 100;
	ent->v.deadflag = DEAD_NO;
	ent->v.takedamage = DAMAGE_AIM;
	ent->v.solid = SOLID_SLIDEBOX;
	ent->v.movetype = MOVETYPE_WALK;
	ent->v.modelindex = 1;
	ent->v.punchangle[0] = ent->v.punchangle[1] = ent->v.punchangle[2] = 0;
	cl.punchangle[0] = cl.punchangle[1] = cl.punchangle[2] = 0;
	ent->v.view_ofs[0] = ent->v.view_ofs[1] = ent->v.angles[2] = cl.viewangles[2] = 0;
	ent->v.view_ofs[2] = cl.viewheight = 22;
	// 3. Reset powerups
	ent->v.effects = 0;
	ent->v.items = (u64)ent->v.items & ~IT_QUAD;
	ent->v.items = (u64)ent->v.items & ~IT_SUIT;
	ent->v.items = (u64)ent->v.items & ~IT_INVISIBILITY;
	ent->v.items = (u64)ent->v.items & ~IT_INVULNERABILITY;
	// 4. Compact Weapon Logic (Ternary Chain)
	s32 w = ent->v.weapon;
	ent->v.weapon = 0;
	ent->v.impulse=(w==IT_AXE)?1:(w==IT_SHOTGUN)?2:(w==IT_SUPER_SHOTGUN)?3:
	  (w==IT_NAILGUN)?4:(w==IT_SUPER_NAILGUN)?5:(w==IT_GRENADE_LAUNCHER)?6:
	  (w==IT_ROCKET_LAUNCHER)?7:(w==IT_LIGHTNING)?8:10;
	// 5. Synchronous VM Execution
	s32 old_self = pr_global_struct->self;
	pr_global_struct->self = EDICT_TO_PROG(ent);
	PR_ExecuteProgram(pr_global_struct->PlayerPostThink);
	pr_global_struct->self = old_self;
	Con_Printf("Resurrection successful.\n");
resurrect_ret:
	PR_SwitchQCVM(NULL);
}

void Host_InitCommands()
{
	Cmd_AddCommand("maps", Host_Maps_f);
	Cmd_AddCommand("mods", Host_Mods_f);
	Cmd_AddCommand("games", Host_Mods_f);
	Cmd_AddCommand("status", Host_Status_f);
	Cmd_AddCommand("quit", Host_Quit_f);
	Cmd_AddCommand("god", Host_God_f);
	Cmd_AddCommand("resurrect", Cmd_Resurrect_f);
	Cmd_AddCommand("notarget", Host_Notarget_f);
	Cmd_AddCommand("fly", Host_Fly_f);
	Cmd_AddCommand("map", Host_Map_f);
	Cmd_AddCommand("restart", Host_Restart_f);
	Cmd_AddCommand("changelevel", Host_Changelevel_f);
	Cmd_AddCommand("connect", Host_Connect_f);
	Cmd_AddCommand("reconnect", Host_Reconnect_f);
	Cmd_AddCommand("name", Host_Name_f);
	Cmd_AddCommand("noclip", Host_Noclip_f);
	Cmd_AddCommand("version", Host_Version_f);
	Cmd_AddCommand("say", Host_Say_f);
	Cmd_AddCommand("say_team", Host_Say_Team_f);
	Cmd_AddCommand("tell", Host_Tell_f);
	Cmd_AddCommand("color", Host_Color_f);
	Cmd_AddCommand("kill", Host_Kill_f);
	Cmd_AddCommand("pause", Host_Pause_f);
	Cmd_AddCommand("spawn", Host_Spawn_f);
	Cmd_AddCommand("begin", Host_Begin_f);
	Cmd_AddCommand("prespawn", Host_PreSpawn_f);
	Cmd_AddCommand("kick", Host_Kick_f);
	Cmd_AddCommand("ping", Host_Ping_f);
	Cmd_AddCommand("load", Host_Loadgame_f);
	Cmd_AddCommand("save", Host_Savegame_f);
	Cmd_AddCommand("give", Host_Give_f);
	Cmd_AddCommand("startdemos", Host_Startdemos_f);
	Cmd_AddCommand("demos", Host_Demos_f);
	Cmd_AddCommand("stopdemo", Host_Stopdemo_f);
	Cmd_AddCommand("viewmodel", Host_Viewmodel_f);
	Cmd_AddCommand("viewframe", Host_Viewframe_f);
	Cmd_AddCommand("viewnext", Host_Viewnext_f);
	Cmd_AddCommand("viewprev", Host_Viewprev_f);
}
