// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.

#include "quakedef.h"

static s32 con_current; // where next message will be printed
static s8 *con_text = 0;
static s32 con_linewidth;
static f32 con_cursorspeed = 4;
static s32 con_vislines;
static bool con_debuglog;
static f32 con_times[NUM_CON_TIMES];
static s32 con_x; // offset in current line for next print

void Con_ToggleConsole_f()
{
	if(key_dest == key_console){
		if(cls.state == ca_connected){
			key_dest = key_game;
			key_lines[edit_line][1] = 0; // clear any typing
			key_linepos = 1;
		} else { M_Menu_Main_f(); }
	} else key_dest = key_console;
	SCR_EndLoadingPlaque();
	memset(con_times, 0, sizeof(con_times));
}

void Con_Clear_f()
{ if(con_text) memset(con_text, ' ', CON_TEXTSIZE); }
void Con_ClearNotify()
{ memset(con_times, 0, sizeof(con_times)); }
void Con_MessageMode_f()
{ key_dest = key_message; team_message = 0; }
void Con_MessageMode2_f()
{ key_dest = key_message; team_message = 1; }

void Con_CheckResize()
{ // If the line width has changed, reformat the buffer.
	s32 i, j, width, oldwidth, oldtotallines, numlines, numchars;
	s8 tbuf[CON_TEXTSIZE];
	width = (vid.width >> 3) - 2;
	if(uiscale) width /= uiscale;
	if(width == con_linewidth) return;
	if(width < 1){	// video hasn't been initialized yet
		width = 38;
		con_linewidth = width;
		con_totallines = CON_TEXTSIZE / con_linewidth;
		Q_memset(con_text, ' ', CON_TEXTSIZE);
	} else {
		oldwidth = con_linewidth;
		con_linewidth = width;
		oldtotallines = con_totallines;
		con_totallines = CON_TEXTSIZE / con_linewidth;
		numlines = oldtotallines;
		if(con_totallines < numlines) numlines = con_totallines;
		numchars = oldwidth;
		if(con_linewidth < numchars) numchars = con_linewidth;
		Q_memcpy(tbuf, con_text, CON_TEXTSIZE);
		Q_memset(con_text, ' ', CON_TEXTSIZE);
		for(i = 0; i < numlines; i++){
		for(j = 0; j < numchars; j++){
			con_text[(con_totallines - 1 - i) * con_linewidth + j] =
				tbuf[((con_current - i + oldtotallines) %
						oldtotallines) * oldwidth + j];
		}}
		Con_ClearNotify();
	}
	con_backscroll = 0;
	con_current = con_totallines - 1;
}

void Con_Init()
{
	s8 temp[MAXGAMEDIRLEN + 1];
	s8 *t2 = "/qconsole.log";
	con_debuglog = COM_CheckParm("-condebug");
	if(con_debuglog){
		if(strlen(com_gamedir) < (MAXGAMEDIRLEN - strlen(t2))){
			sprintf(temp, "%s%s", com_gamedir, t2);
#ifdef _WIN32
			_unlink(temp);
#else
			unlink(temp);
#endif
		}
	}
	con_text = Hunk_AllocName(CON_TEXTSIZE, "context");
	Q_memset(con_text, ' ', CON_TEXTSIZE);
	con_linewidth = -1;
	Con_CheckResize();
	Con_Printf("Console initialized.\n");
	Cvar_RegisterVariable(&con_notifytime);
	Cvar_RegisterVariable(&con_logcenterprint);
	Cmd_AddCommand("toggleconsole", Con_ToggleConsole_f);
	Cmd_AddCommand("messagemode", Con_MessageMode_f);
	Cmd_AddCommand("messagemode2", Con_MessageMode2_f);
	Cmd_AddCommand("clear", Con_Clear_f);
	con_initialized = 1;
}

void Con_Linefeed()
{
	if(!con_initialized) return;
	con_x = 0;
	con_current++;
	Q_memset(&con_text[(con_current % con_totallines) * con_linewidth]
		 , ' ', con_linewidth);
}

// Handles cursor positioning, line wrapping, etc
// All console printing must go through this in order to be logged to disk
// If no console is visible, the notify window will pop up (if notify != 0).
void Con_Print(s8 *txt, s32 notify)
{
	static s32 cr;
	if(!con_initialized) return;
	con_backscroll = 0;
	s32 mask = 0;
	if(txt[0] == 1){
		mask = 128; // go to colored text
		S_LocalSound("misc/talk.wav"); // play talk wav
		txt++;
	} else if(txt[0] == 2){
		mask = 128; // go to colored text
		txt++;
	}
	s32 c, l, y; // keep here for OpenBSD compiler
	while((c = *txt)){
		for(l = 0; l < con_linewidth; l++) // count word length
			if(txt[l] <= ' ') break;
		if(l != con_linewidth && (con_x+l>con_linewidth)) // word wrap
			con_x = 0;
		txt++;
		if(cr){ con_current--; cr = 0; }
		if(!con_x){
			Con_Linefeed();
			if(con_current >= 0 && notify) // top-left overlay
				con_times[con_current%NUM_CON_TIMES] = realtime;
		}
		switch(c){
		case '\n': con_x = 0; break;
		case '\r': con_x = 0; cr = 1; break;
		default: // display character and advance
			y = con_current % con_totallines;
			con_text[y * con_linewidth + con_x] = c | mask;
			con_x++;
			con_x = con_x >= con_linewidth ? 0 : con_x;
			break;
		}
	}
}

void Con_DebugLog(s8 *file, s8 *fmt, ...)
{
	va_list argptr;
	va_start(argptr, fmt);
	static s8 data[1024];
	vsprintf(data, fmt, argptr);
	va_end(argptr);
	s32 fd = open(file, O_WRONLY | O_CREAT | O_APPEND, 0666);
	write(fd, data, strlen(data));
	close(fd);
}

void Con_Printf(const s8 *fmt, ...)
{ // Handles cursor positioning, line wrapping, etc
	va_list argptr;
	va_start(argptr, fmt);
	s8 msg[MAXPRINTMSG];
	vsprintf(msg, fmt, argptr);
	va_end(argptr);
	Sys_Printf("%s", msg);// also echo to debugging console
	if(con_debuglog) // log all messages to file
		Con_DebugLog(va("%s/qconsole.log", com_gamedir), "%s", msg);
	if(!con_initialized || cls.state == ca_dedicated)
		return;
	Con_Print(msg, 1); // write it to the scrollable buffer
}

void Con_DPrintf(s8 *fmt, ...)
{ // A Con_Printf that only shows up if the "developer" cvar is set
	va_list argptr;
	s8 msg[MAXPRINTMSG];
	if(!developer.value)
		return; // don't confuse non-developers with techie stuff...
	va_start(argptr, fmt);
	vsprintf(msg, fmt, argptr);
	va_end(argptr);
	Con_Printf("%s", msg);
}

void Con_SafePrintf(s8 *fmt, ...)
{ // Okay to call even when the screen can't be updated
	s8 msg[1024];
	va_list argptr;
	va_start(argptr, fmt);
	vsprintf(msg, fmt, argptr);
	va_end(argptr);
	s32 temp = scr_disabled_for_loading;
	scr_disabled_for_loading = 1;
	Con_Printf("%s", msg);
	scr_disabled_for_loading = temp;
}


void Con_DrawInput()
{ // The input line scrolls horizontally if typing goes beyond the right edge
	if(key_dest != key_console && !con_forcedup)
		return; // don't draw anything
	s8 *text = key_lines[edit_line];
	// add the cursor frame
	text[key_linepos] = 10 + ((s32)(realtime * con_cursorspeed) & 1);
	for(s32 i = key_linepos + 1; i < con_linewidth; i++)
		text[i] = ' '; // fill out remainder with spaces
	if(key_linepos >= con_linewidth) // prestep if horizontally scrolling
		text += 1 + key_linepos - con_linewidth;
	for(s32 i = 0; i < con_linewidth; i++) // draw it
		Draw_CharacterScaled(((i + 1) << 3) * uiscale,
			con_vislines - 16 * uiscale, text[i], uiscale);
	key_lines[edit_line][key_linepos] = 0; // remove cursor
}


void Con_DrawNotify()
{ // Draws the last few lines of output transparently over the game top
	extern s8 chat_buffer[];
	s32 x = 0, v = 0;
	drawlayer = lyr_notify.value;
	for(s32 i = con_current - NUM_CON_TIMES + 1; i <= con_current; i++){
		if(i < 0) continue;
		f32 time = con_times[i % NUM_CON_TIMES];
		if(time == 0) continue;
		time = realtime - time;
		if(time > con_notifytime.value) continue;
		s8 *text = con_text + (i % con_totallines) * con_linewidth;
		clearnotify = 0;
		for(x = 0; x < con_linewidth; x++)
			Draw_CharacterScaled(((x + 1) << 3) * uiscale,
					     v * uiscale, text[x], uiscale);
		v += 8;
	}
	if(key_dest == key_message){
		clearnotify = 0;
		x = 0;
		Draw_StringScaled(8 * uiscale, v * uiscale, "say: ", uiscale);
		while(chat_buffer[x]){
			Draw_CharacterScaled(((x + 5) << 3) * uiscale,
					v * uiscale, chat_buffer[x], uiscale);
			x++;
		}
		Draw_CharacterScaled(((x + 5) << 3) * uiscale, v * uiscale,
			10 + ((s32)(realtime * con_cursorspeed) & 1), uiscale);
		v += 8 * uiscale;
	}
	if(v > con_notifylines) con_notifylines = v;
	drawlayer = lyr_main.value;
}

void Con_DrawConsole(s32 lines, bool drawinput) // Draws console with solid bg
{ // Typing input line at the bottom should only be drawn if typing is allowed
	s8 *text;
	if(lines <= 0) return;
	drawlayer = lyr_console.value;
	Draw_ConsoleBackground(lines); // draw the background
	con_vislines = lines; // draw the text
	s32 rows = ((lines - 16) >> 3) * uiscale; // rows of text to draw
	s32 y = lines-16*uiscale-(rows<<3)*uiscale; // may start a bit negative
	for(s32 i = con_current-rows+1; i <= con_current; i++, y += 8*uiscale){
		s32 j = i - con_backscroll;
		j = j < 0 ? 0 : j;
		text = con_text + (j % con_totallines) * con_linewidth;
		for(s32 x = 0; x < con_linewidth; x++)
			Draw_CharacterScaled(((x + 1) << 3) * uiscale, y,
					     text[x], uiscale);
	} // draw the input prompt, user text and cursor if desired
	if(drawinput) Con_DrawInput();
	drawlayer = lyr_main.value;
}
