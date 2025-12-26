// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

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
	s32 i;
	s8 name[MAX_QPATH], *p;
	if(Cmd_Argc() < 2){ //no map name given
		if(cls.state == ca_dedicated) {
			if(sv.active) Con_Printf("Current map: %s\n", sv.name);
			else Con_Printf("Server not active\n");
		}
		else if(cls.state == ca_connected)
		Con_Printf("Current map: %s( %s )\n", cl.levelname, cl.mapname);
		else Con_Printf("map <levelname>: start a new server\n");
		return;
	}
	if(cmd_source != src_command) return;
	cls.demonum = -1; // stop demo loop in case this fails
	CL_Disconnect();
	Host_ShutdownServer(0);
	key_dest = key_game; // remove console or menu
	SCR_BeginLoadingPlaque();
	svs.serverflags = 0; // haven't completed an episode yet
	q_strlcpy(name, Cmd_Argv(1), sizeof(name));
	p = strstr(name, ".bsp"); // remove trailing ".bsp" from mapname -- S.A.
	if(p && p[4] == '\0') *p = '\0';
	SV_SpawnServer(name);
	if(!sv.active) return;
	if(cls.state != ca_dedicated) {
		memset(cls.spawnparms, 0, MAX_MAPSTRING);
		for(i = 2; i < Cmd_Argc(); i++) {
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
	SV_SaveSpawnparms();
	strcpy(level, Cmd_Argv(1));
	SV_SpawnServer(level);
}

void Host_Restart_f()
{ // Restarts the current server for a dead player
	s8 mapname[MAX_QPATH];
	if(cls.demoplayback || !sv.active) return;
	if(cmd_source != src_command) return;
	strcpy(mapname, sv.name); // must copy out, because it gets
	SV_SpawnServer(mapname); // cleared in sv_spawnserver
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
	FILE *f = fopen(name, "w");
	if(!f){ Con_Printf("ERROR: couldn't open.\n"); return; }
	fprintf(f, "%i\n", SAVEGAME_VERSION);
	Host_SavegameComment(comment);
	fprintf(f, "%s\n", comment);
	for(s32 i = 0; i < NUM_SPAWN_PARMS; i++)
		fprintf(f, "%f\n", svs.clients->spawn_parms[i]);
	fprintf(f, "%d\n", current_skill);
	fprintf(f, "%s\n", sv.name);
	fprintf(f, "%f\n", sv.time);
	for(s32 i = 0; i < MAX_LIGHTSTYLES; i++){ // write the light styles
		if(sv.lightstyles[i]) fprintf(f, "%s\n", sv.lightstyles[i]);
		else fprintf(f, "m\n");
	}
	ED_WriteGlobals(f);
	for(s32 i = 0; i < sv.num_edicts; i++){
		ED_Write(f, EDICT_NUM(i));
		fflush(f);
	}
	fclose(f);
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
	SV_SpawnServer(mapname);
	if(!sv.active){ Con_Printf("Couldn't load map\n"); return; }
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
			memset(&ent->v, 0, progs->entityfields * 4);
			ent->free = 0;
			ED_ParseEdict(start, ent);
			if(!ent->free) // link it into the bsp tree
				SV_LinkEdict(ent, 0);
		}
		entnum++;
	}
	sv.num_edicts = entnum;
	sv.time = time;
	fclose(f);
	for(s32 i = 0; i < NUM_SPAWN_PARMS; i++)
		svs.clients->spawn_parms[i] = spawn_parms[i];
	if(cls.state != ca_dedicated){
		CL_EstablishConnection("local");
		Host_Reconnect_f();
	}
}

static void Host_Name_f()
{
	s8 newName[32];
	if(Cmd_Argc() == 1) {
		Con_Printf("\"name\" is \"%s\"\n", cl_name.string);
		return;
	}
	if(Cmd_Argc() == 2) q_strlcpy(newName, Cmd_Argv(1), sizeof(newName));
	else q_strlcpy(newName, Cmd_Args(), sizeof(newName));
	newName[15] = 0; // client_t structure actually says name[32].
	if(cmd_source == src_command) {
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
	pr_global_struct->time = sv.time;
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
	VID_SetPalette(CURWORLDPAL, screen);
	if(cmd_source == src_command) {
		Con_Printf("spawn is not valid from the console\n");
		return;
	}
	if(host_client->spawned) {
		Con_Printf("Spawn not valid -- already spawned\n");
		return;
	}
	// run the entrance script
	if(sv.loadgame) // loaded games are fully inited already
		sv.paused = 0; //if this is last client to be connected, unpause
	else {
		edict_t *ent = host_client->edict; // set up the edict
		memset(&ent->v, 0, progs->entityfields * 4);
		ent->v.colormap = NUM_FOR_EDICT(ent);
		ent->v.team = (host_client->colors & 15) + 1;
		ent->v.netname = PR_SetEngineString(host_client->name);
		for(s32 i = 0; i < NUM_SPAWN_PARMS; i++) // copy out of client_t
		    (&pr_global_struct->parm1)[i] = host_client->spawn_parms[i];
		// call the spawn function
		pr_global_struct->time = sv.time;
		pr_global_struct->self = EDICT_TO_PROG(sv_player);
		PR_ExecuteProgram(pr_global_struct->ClientConnect);
		if((Sys_DoubleTime() - host_client->netconnection->connecttime)
				<= sv.time)
			Sys_Printf("%s entered the game\n", host_client->name);
		PR_ExecuteProgram(pr_global_struct->PutClientInServer);
	}

	SZ_Clear(&host_client->message);//send all current names colors frags
	MSG_WriteByte(&host_client->message, svc_time); // send time of update
	MSG_WriteFloat(&host_client->message, sv.time);
	client_t *client = svs.clients;
	for(s32 i = 0; i < svs.maxclients; i++, client++) {
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
	for(s32 i = 0; i < MAX_LIGHTSTYLES; i++){//send all current light styles
		MSG_WriteByte(&host_client->message, svc_lightstyle);
		MSG_WriteByte(&host_client->message, (s8)i);
		MSG_WriteString(&host_client->message, sv.lightstyles[i]);
	}
	MSG_WriteByte(&host_client->message, svc_updatestat); // send some stats
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
	edict_t *ent = EDICT_NUM( 1 + (host_client - svs.clients) );
	MSG_WriteByte(&host_client->message, svc_setangle);
	for(s32 i = 0; i < 2; i++)
		MSG_WriteAngle(&host_client->message,ent->v.angles[i],
				sv.protocolflags);
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
	const s8 *who;
	const s8 *message = NULL;
	client_t *save;
	s32 i;
	bool byNumber = 0;
	if(cmd_source == src_command){
		if(!sv.active){
			Cmd_ForwardToServer();
			return;
		}
	} else if(pr_global_struct->deathmatch) return;
	save = host_client;
	if(Cmd_Argc() > 2 && Q_strcmp(Cmd_Argv(1), "#") == 0){
		i = Q_atof(Cmd_Argv(2)) - 1;
		if(i < 0 || i >= svs.maxclients) return;
		if(!svs.clients[i].active) return;
		host_client = &svs.clients[i];
		byNumber = 1;
	} else {
		for(i = 0, host_client = svs.clients; i < svs.maxclients; i++, host_client++){
			if(!host_client->active) continue;
			if(q_strcasecmp(host_client->name,Cmd_Argv(1))==0)break;
		}
	}
	if(i < svs.maxclients){
		if(cmd_source == src_command)
		      who=cls.state==ca_dedicated?(s8*)"Console":cl_name.string;
		else who = save->name;
		if(host_client == save) return; // can't kick yourself!
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
	if (hipnotic) {
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
		if(rogue){ val = GetEdictFieldValue(sv_player, "ammo_shells1");
			if(val) val->_float = v; }
		sv_player->v.ammo_shells = v; break;
	case 'n':
		if(rogue){ val = GetEdictFieldValue(sv_player, "ammo_nails1");
			if(val){ val->_float = v;
				if(sv_player->v.weapon <= IT_LIGHTNING)
					sv_player->v.ammo_nails = v; }
		} else sv_player->v.ammo_nails = v; break;
	case 'l':
		if(rogue){val=GetEdictFieldValue(sv_player, "ammo_lava_nails");
			if(val){ val->_float = v;
				if(sv_player->v.weapon > IT_LIGHTNING)
					sv_player->v.ammo_nails = v; } } break;
	case 'r':
		if(rogue){ val = GetEdictFieldValue(sv_player, "ammo_rockets1");
			if(val){ val->_float = v;
				if(sv_player->v.weapon <= IT_LIGHTNING)
					sv_player->v.ammo_rockets = v; }
		} else sv_player->v.ammo_rockets = v; break;
	case 'm':
	       if(rogue){val=GetEdictFieldValue(sv_player,"ammo_multi_rockets");
			if(val){ val->_float = v;
				if(sv_player->v.weapon > IT_LIGHTNING)
					sv_player->v.ammo_rockets = v; }} break;
	case 'h': sv_player->v.health = v; break;
	case 'c':
		if(rogue){ val = GetEdictFieldValue(sv_player, "ammo_cells1");
			if(val){ val->_float = v;
				if(sv_player->v.weapon <= IT_LIGHTNING)
					sv_player->v.ammo_cells = v; }
		} else sv_player->v.ammo_cells = v; break;
	case 'p':
		if(rogue){ val = GetEdictFieldValue(sv_player, "ammo_plasma");
			if(val){ val->_float = v;
				if(sv_player->v.weapon > IT_LIGHTNING)
					sv_player->v.ammo_cells = v;
			} } break;
	}
}

edict_t *FindViewthing()
{
	for(s32 i = 0; i < sv.num_edicts; i++){
		edict_t *e = EDICT_NUM(i);
		if(!strcmp(PR_GetString(e->v.classname), "viewthing"))
			return e;
	}
	Con_Printf("No viewthing on map\n");
	return NULL;
}

void Host_Viewmodel_f()
{
	edict_t *e = FindViewthing();
	if(!e) return;
	model_t *m = Mod_ForName(Cmd_Argv(1), 0);
	if(!m){ Con_Printf("Can't load %s\n", Cmd_Argv(1)); return; }
	e->v.frame = 0;
	cl.model_precache[(s32)e->v.modelindex] = m;
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

void Host_InitCommands()
{
	Cmd_AddCommand("status", Host_Status_f);
	Cmd_AddCommand("quit", Host_Quit_f);
	Cmd_AddCommand("god", Host_God_f);
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
