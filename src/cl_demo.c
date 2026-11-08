// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2005 John Fitzgibbons and others
// Copyright(C) 2007-2008 Kristian Duske
// GPLv3 See LICENSE for details.
// When a demo is playing back, all NET_SendMessages are skipped, and
// NET_GetMessages are read from the demo file.
// Whenever cl.time gets past the last received message, another message is
// read from the demo file.
#include "quakedef.h"

// timedemo_threads: the demo with each r_threads from 0 up to what -1 would
// try, then with -1, and a table of them to compare computers with
static struct {
	c8 demo[MAX_QPATH];
	s32 run, last; // the r_threads of this run, -1 the last one
	f32 fps[R_MAXWORKERS + 1];
	c8 saved[32];
} td_threads;

static void CL_TimeDemoThreadsRun(s32 n)
{
	td_threads.run = n;
	Cvar_SetValue("r_threads", n);
	Cbuf_AddText(va("timedemo %s\n", td_threads.demo));
}

static void CL_TimeDemoThreadsStop()
{
	Cvar_Set("r_threads", td_threads.saved);
	td_threads.demo[0] = 0;
}

static void CL_TimeDemoThreadsDone(f32 fps)
{
	if(td_threads.run >= 0){
		td_threads.fps[td_threads.run] = fps;
		CL_TimeDemoThreadsRun(td_threads.run < td_threads.last
				? td_threads.run + 1 : -1);
		return;
	}
	Con_Printf("r_threads   fps\n");
	for(s32 i = 0; i <= td_threads.last; i++)
		Con_Printf("%9d %6.1f\n", i, td_threads.fps[i]);
	Con_Printf("%9d %6.1f, with %d at the end\n", -1, fps,
			R_ThreadsActive());
	CL_TimeDemoThreadsStop();
}

void CL_TimeDemoThreads_f()
{ // timedemo_threads [demoname [most]]
	if(cmd_source != src_command) return;
	if(td_threads.demo[0]){
		CL_TimeDemoThreadsStop();
		Con_Printf("timedemo_threads stopped\n");
		return;
	}
	if(Cmd_Argc() < 2){
		Con_Printf("timedemo_threads <demoname> [most] : demo speeds "
			"with each r_threads up to most and with -1, run it "
			"again to stop\n");
		return;
	}
	SDL_strlcpy(td_threads.saved, r_threads.string,
			sizeof(td_threads.saved));
	SDL_strlcpy(td_threads.demo, Cmd_Argv(1), sizeof(td_threads.demo));
	td_threads.last = Cmd_Argc() > 2
		? CLAMP(0, atoi(Cmd_Argv(2)), R_MAXWORKERS)
		: R_ThreadsAutoMax();
	CL_TimeDemoThreadsRun(0);
}

void CL_FinishTimeDemo()
{
	cls.timedemo = 0; // the first frame didn't count
	s32 frames = (host_framecount - cls.td_startframe) - 1;
	f32 t = realtime - cls.td_starttime;
	if(!t) t = 1;
	Con_Printf("%i frames %5.1f seconds %5.1f fps\n", frames, t, frames/t);
	if(td_threads.demo[0])
		CL_TimeDemoThreadsDone(frames/t);
}

void CL_StopPlayback()
{ // Called when a demo file runs out, or the user starts a game
	if(!cls.demoplayback) return;
	fclose(cls.demofile);
	cls.demoplayback = 0;
	cls.demofile = NULL;
	cls.state = ca_disconnected;
	if(cls.timedemo)
		CL_FinishTimeDemo();
}

void CL_WriteDemoMessage()
{ // Dumps the current net message, prefixed by the length and view angles
	s32 len = net_message.cursize;
	fwrite(&len, 4, 1, cls.demofile);
	for(s32 i = 0; i < 3; i++){
		f32 f = cl.viewangles[i];
		fwrite(&f, 4, 1, cls.demofile);
	}
	fwrite(net_message.data, net_message.cursize, 1, cls.demofile);
	fflush(cls.demofile);
}

s32 CL_GetMessage()
{ // Handles recording and playback of demos, on top of NET_ code
	if(cls.demoplayback){ // decide if it is time to grab the next message
		if(cls.signon == SIGNONS){ // allways grab until fully connected
			if(cls.timedemo){
				if(host_framecount == cls.td_lastframe)
					return 0; // allready read this frame
				cls.td_lastframe = host_framecount;
				// if second frame, grab the real td_starttime
				// so the bogus time on the first doesn't count
				if(host_framecount == cls.td_startframe + 1)
					cls.td_starttime = realtime;
			}else if(cl.time <= cl.mtime[0])
				return 0; // don't need another message yet
		} // get the next message
		fread(&net_message.cursize, 4, 1, cls.demofile);
		VectorCopy(cl.mviewangles[0], cl.mviewangles[1]);
		f32 f;
		for(s32 i = 0; i < 3; i++){
			fread(&f, 4, 1, cls.demofile);
			cl.mviewangles[0][i] = f;
		}
		if(net_message.cursize > MAX_MSGLEN)
			Sys_Error("Demo message > MAX_MSGLEN");
		if(fread(net_message.data, net_message.cursize,
				1, cls.demofile) != 1){
			CL_StopPlayback();
			return 0;
		}
		return 1;
	}
	s32 r;
	while(1){
		r = NET_GetMessage(cls.netcon);
		if(r != 1 && r != 2) return r;
		if(net_message.cursize == 1 && net_message.data[0] == svc_nop)
			Con_Printf("<-- server to client keepalive\n");
		else break;
	}
	if(cls.demorecording) CL_WriteDemoMessage();
	return r;
}

void CL_Stop_f()
{ // stop recording a demo
	if(cmd_source != src_command) return;
	if(!cls.demorecording){
		Con_Printf("Not recording a demo.\n");
		return;
	}
	SZ_Clear(&net_message); // write a disconnect message to the demo file
	MSG_WriteByte(&net_message, svc_disconnect);
	CL_WriteDemoMessage();
	fclose(cls.demofile); // finish up
	cls.demofile = NULL;
	cls.demorecording = 0;
	Con_Printf("Completed demo\n");
}


void CL_Record_f()
{ // record <demoname> <map> [cd track]
	if(cmd_source != src_command) return;
	s32 c = Cmd_Argc();
	if(c != 2 && c != 3 && c != 4){
		Con_Printf("record <demoname> [<map> [cd track]]\n");
		return;
	}
	if(strstr(Cmd_Argv(1), "..")){
		Con_Printf("Relative pathnames are not allowed.\n");
		return;
	}
	if(c == 2 && cls.state == ca_connected){
		Con_Printf("Can not record - already connected to server\nClient demo recording must be started before connecting\n");
		return;
	}
	s32 track = -1; // write the forced cd track number, or -1
	if(c == 4){
		track = atoi(Cmd_Argv(3));
		Con_Printf("Forcing CD track to %i\n", cls.forcetrack);
	}
	c8 name[MAX_OSPATH+2];
	sprintf(name, "%s/%s", com_gamedir, Cmd_Argv(1));
	if(c > 2) // start the map up
		Cmd_ExecuteString(va("map %s", Cmd_Argv(2)), src_command);
	COM_AddExtension(name, ".dem", sizeof(name));
	Con_Printf("recording to %s.\n", name);
	cls.demofile = fopen(name, "wb");
	if(!cls.demofile){
		Con_Printf("ERROR: couldn't open.\n");
		return;
	}
	cls.forcetrack = track;
	fprintf(cls.demofile, "%i\n", cls.forcetrack);
	cls.demorecording = 1;
}


void CL_PlayDemo_f()
{ // play [demoname]
	if(cmd_source != src_command) return;
	if(Cmd_Argc() != 2){
		Con_Printf("play <demoname> : plays a demo\n");
		return;
	}
	if(key_dest == key_console) key_dest = key_game; //withdraw console/menu
	CL_Disconnect(); // disconnect from server
	c8 name[256];
	strcpy(name, Cmd_Argv(1)); // open the demo file
	COM_AddExtension(name, ".dem", sizeof(name));
	Con_Printf("Playing demo from %s.\n", name);
	COM_FOpenFile(name, &cls.demofile, NULL);
	if(!cls.demofile){
		Con_Printf("ERROR: couldn't open.\n");
		cls.demonum = -1; // stop demo loop
		return;
	}
	cls.demoplayback = 1;
	cls.state = ca_connected;
	r_numskies = 0;
	fscanf(cls.demofile, "%i\n", &cls.forcetrack);
}

void CL_TimeDemo_f()
{ // timedemo [demoname]
	if(cmd_source != src_command) return;
	if(Cmd_Argc() != 2){
		Con_Printf("timedemo <demoname> : gets demo speeds\n");
		return;
	}
	CL_PlayDemo_f();
	// cls.td_starttime will be grabbed at the second frame of the demo so
	// all the loading time doesn't get counted
	cls.timedemo = 1;
	cls.td_startframe = host_framecount;
	cls.td_lastframe = -1;// get a new message this frame
}
