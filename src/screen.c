// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// screen.c -- master for refresh, status bar, console, chat, notify, etc
#include "quakedef.h"

static bool scr_initialized; // ready to draw
static qpic_t *scr_ram;
static qpic_t *scr_net;
static qpic_t *scr_turtle;
static u32 clearconsole;
static bool scr_drawloading;
static f32 scr_disabled_time;
static s8 scr_centerstring[1024]; // center printing
static f32 scr_centertime_start; // for slow victory printing
static s32 scr_center_lines;
static s32 scr_erase_lines;
static u32 scr_erase_center;
static f32 oldscreensize, oldfov;
static f32 scr_conlines; // lines of console to display
static s8 *scr_notifystring;
static bool scr_drawdialog;

void SCR_ScreenShot_f();
void SCR_HUDStyle_f (cvar_t *cvar);

void SCR_CenterPrint(const s8 *str) // Called for important messages
{ // that should stay in the center of the screen for a few moments
	if (con_logcenterprint.value) {
		Con_Print((s8*)str, 0);
		Con_Print("\n", 0);
	}
	strncpy(scr_centerstring, str, sizeof(scr_centerstring) - 1);
	scr_centertime_off = scr_centertime.value;
	scr_centertime_start = cl.time;
	scr_center_lines = 1; // count the number of lines for centering
	while (*str) {
		if (*str == '\n')
			scr_center_lines++;
		str++;
	}
}

void SCR_EraseCenterString()
{
	if (scr_erase_center++ > vid.numpages) {
		scr_erase_lines = 0;
		return;
	}
	s32 y = 48;
	if (scr_center_lines <= 4)
		y = vid.height * 0.35;
	Draw_TileClear(0, y, vid.width, 8 * scr_erase_lines);
}

void SCR_DrawCenterString()
{
	s32 remaining = cl.intermission ? scr_printspeed.value // the finale
		* (cl.time - scr_centertime_start) : 9999; // s8-by-s8 print
	scr_erase_center = 0;
	s8 *start = scr_centerstring;
	s32 y = scr_center_lines <= 4 ? vid.height * 0.35 : 48 * uiscale;
	drawlayer = lyr_centerprint.value;
	do {
		s32 l = 0;
		for (; l < 40; l++) // scan the width of the line
			if (start[l] == '\n' || !start[l])
				break;
		s32 x = (vid.width - l * 8 * uiscale) / 2;
		for (s32 j = 0; j < l; j++, x += 8 * uiscale) {
			Draw_CharacterScaled(x, y, start[j], uiscale);
			if (!remaining--) { drawlayer = lyr_main.value; return; }
		}
		y += 8 * uiscale;
		while (*start && *start != '\n')
			start++;
		if (!*start)
			break;
		start++; // skip the \n
	} while (1);
	drawlayer = lyr_main.value;
}

void SCR_CheckDrawCenterString()
{
	if (scr_center_lines > scr_erase_lines)
		scr_erase_lines = scr_center_lines;
	scr_centertime_off -= host_frametime;
	if ((scr_centertime_off <= 0 && !cl.intermission)
		|| key_dest != key_game)
		return;
	SCR_DrawCenterString();
}

f32 CalcFov(f32 fov_x, f32 width, f32 height)
{
	if (fov_x < 1 || fov_x > 179)
		Sys_Error("Bad fov: %f", fov_x);
	return atan(height / (width / tan(fov_x / 360 * M_PI))) * 360 / M_PI;
}

static void SCR_CalcRefdef() // Must be called whenever vid changes
{ // Internal use only
	fog_initialized = 0;
	scr_fullupdate = 0; // force a background redraw
	vid.recalc_refdef = 0;
	Sbar_CalcPos();
	Sbar_Changed(); // force the status bar to redraw
	if (scr_viewsize.value < 30) // bound viewsize
		Cvar_Set("viewsize", "30");
	if (scr_viewsize.value > 120)
		Cvar_Set("viewsize", "120");
	if (scr_fov.value < 10) // bound field of view
		Cvar_Set("fov", "10");
	if (scr_fov.value > 170)
		Cvar_Set("fov", "170");
	r_refdef.fov_x = scr_fov.value;
	r_refdef.fov_y = CalcFov(r_refdef.fov_x,
		r_refdef.vrect.width, r_refdef.vrect.height);
	f32 size = cl.intermission ? 120 : scr_viewsize.value;
	if (size >= 120 || hudstyle != HUD_CLASSIC)
		sb_lines = 0; // no status bar at all
	else if (size >= 110)
		sb_lines = 24 * uiscale; // no inventory
	else
		sb_lines = (24 + 16 + 8) * uiscale;
	vrect_t vrect;
	vrect.x = 0; // these calculations mirror those in R_Init()
	vrect.y = 0; // for r_refdef, but take no account of water warping
	vrect.width = vid.width;
	vrect.height = vid.height;
	R_SetVrect(&vrect, &scr_vrect, sb_lines);
	// guard against going from one mode to another that's
	// less than half the vertical resolution
	if (scr_con_current > vid.height)
		scr_con_current = vid.height;
	// notify the refresh of the change
	R_ViewChanged(&vrect, sb_lines, vid.aspect);
}

void SCR_SizeUp_f() // Keybinding command
{
	Cvar_SetValue("viewsize", scr_viewsize.value + 10);
	vid.recalc_refdef = 1;
}

void SCR_SizeDown_f() // Keybinding command
{
	Cvar_SetValue("viewsize", scr_viewsize.value - 10);
	vid.recalc_refdef = 1;
}

void SCR_Init()
{
	scr_ram = Draw_PicFromWad("ram");
	scr_net = Draw_PicFromWad("net");
	scr_turtle = Draw_PicFromWad("turtle");
	if (scr_initialized) return;
	Cvar_RegisterVariable(&scr_fov);
	Cvar_RegisterVariable(&scr_viewsize);
	Cvar_RegisterVariable(&scr_conspeed);
	Cvar_RegisterVariable(&scr_showram);
	Cvar_RegisterVariable(&scr_showturtle);
	Cvar_RegisterVariable(&scr_showpause);
	Cvar_RegisterVariable(&scr_centertime);
	Cvar_RegisterVariable(&scr_printspeed);
	Cvar_RegisterVariable(&scr_showfps);
	Cvar_RegisterVariable (&scr_hudstyle);
	Cvar_SetCallback (&scr_hudstyle, SCR_HUDStyle_f);
	Cmd_AddCommand("screenshot", SCR_ScreenShot_f);
	Cmd_AddCommand("sizeup", SCR_SizeUp_f);
	Cmd_AddCommand("sizedown", SCR_SizeDown_f);
	scr_initialized = 1;
}

void SCR_DrawFPS()
{ // -- johnfitz
	static f64 oldtime = 0;
	static f64 lastfps = 0;
	static s32 oldframecount = 0;
	if (con_forcedup) {
		oldtime = realtime;
		oldframecount = r_framecount;
		lastfps = 0;
		return;
	}
	f64 elapsed_time = realtime - oldtime;
	s32 frames = r_framecount - oldframecount;
	if (elapsed_time < 0 || frames < 0) {
		oldtime = realtime;
		oldframecount = r_framecount;
		return;
	}
	if (elapsed_time > 0.75) { // update value every 3/4 second
		lastfps = frames / elapsed_time;
		oldtime = realtime;
		oldframecount = r_framecount;
	}
	if (scr_showfps.value && scr_viewsize.value < 130 && lastfps) {
		s8 st[16];
		if (scr_showfps.value > 0.f)
			sprintf(st, "%4.0f fps", lastfps);
		else
			sprintf(st, "%.2f ms", 1000.f / lastfps);
		s32 x = vid.width - (strlen(st) * 8 * uiscale);
		Draw_StringScaled(x, 0, st, uiscale);
	}
}

void SCR_DrawRam()
{
	if (!scr_showram.value || !r_cache_thrash)
		return;
	Draw_PicScaled(scr_vrect.x + 32*uiscale, scr_vrect.y, scr_ram, uiscale);
}

void SCR_DrawTurtle()
{
	static s32 count;
	if (!scr_showturtle.value)
		return;
	if (host_frametime < 0.1) {
		count = 0;
		return;
	}
	count++;
	if (count < 3)
		return;
	Draw_PicScaled(scr_vrect.x, scr_vrect.y, scr_turtle, uiscale);
}

void SCR_DrawNet()
{
	if (realtime - cl.last_received_message < 0.3 || cls.demoplayback)
		return;
	Draw_PicScaled(scr_vrect.x + 64*uiscale, scr_vrect.y, scr_net, uiscale);
}

void SCR_DrawPause()
{
	if (!scr_showpause.value || !cl.paused) // turn off for screenshots
		return;
	qpic_t *pic = Draw_CachePic("gfx/pause.lmp");
	Draw_TransPicScaled((vid.width - pic->width * uiscale) / 2,
		       (vid.height - 48 * uiscale - pic->height * uiscale) / 2,
		       pic, uiscale);
}

void SCR_DrawLoading()
{
	if (!scr_drawloading)
		return;
	qpic_t *pic = Draw_CachePic("gfx/loading.lmp");
	Draw_TransPicScaled((vid.width - pic->width * uiscale) / 2,
		       (vid.height - 48 * uiscale - pic->height * uiscale) / 2,
		       pic, uiscale);
}

void SCR_SetUpToDrawConsole()
{
	Con_CheckResize();
	if (scr_drawloading)
		return;	// never a console with loading plaque
	con_forcedup = !cl.worldmodel || cls.signon != SIGNONS; // height
	if (con_forcedup) {
		scr_conlines = vid.height; // full screen
		scr_con_current = scr_conlines;
	} else if (key_dest == key_console)
		scr_conlines = vid.height / 2; // half screen
	else
		scr_conlines = 0; // none visible
	if (scr_conlines < scr_con_current) {
		scr_con_current -= scr_conspeed.value * host_frametime *uiscale;
		if (scr_conlines > scr_con_current)
			scr_con_current = scr_conlines;

	} else if (scr_conlines > scr_con_current) {
		scr_con_current += scr_conspeed.value * host_frametime *uiscale;
		if (scr_conlines < scr_con_current)
			scr_con_current = scr_conlines;
	}
	if (clearconsole++ < vid.numpages) {
		Draw_TileClear(0, (s32)scr_con_current, vid.width,
			       vid.height - (s32)scr_con_current);
		Sbar_Changed();
	} else if (clearnotify++ < vid.numpages) {
		Draw_TileClear(0, 0, vid.width, con_notifylines);
	} else
		con_notifylines = 0;
}

void SCR_DrawConsole()
{
	if (scr_con_current) {
		Con_DrawConsole(scr_con_current, 1);
		clearconsole = 0;
	} else
		if (key_dest == key_game || key_dest == key_message)
			Con_DrawNotify(); // only draw notify in game
}

void WritePCXfile(s8 *filename, u8 *data, s32 width, s32 height,
		  s32 rowbytes, u8 *palette)
{
	pcx_t *pcx = Hunk_TempAlloc(width * height * 2 + 1000);
	if (pcx == NULL) {
		Con_Printf("SCR_ScreenShot_f: not enough memory\n");
		return;
	}
	pcx->manufacturer = 0x0a; // PCX id
	pcx->version = 5; // 256 color
	pcx->encoding = 1; // uncompressed
	pcx->bits_per_pixel = 8; // 256 color
	pcx->xmin = 0;
	pcx->ymin = 0;
	pcx->xmax = LittleShort((s16)(width - 1));
	pcx->ymax = LittleShort((s16)(height - 1));
	pcx->hres = LittleShort((s16)width);
	pcx->vres = LittleShort((s16)height);
	Q_memset(pcx->palette, 0, sizeof(pcx->palette));
	pcx->color_planes = 1; // chunky image
	pcx->bytes_per_line = LittleShort((s16)width);
	pcx->palette_type = LittleShort(2); // not a grey scale
	Q_memset(pcx->filler, 0, sizeof(pcx->filler));
	u8 *pack = &pcx->data; // pack the image
	for (s32 i = 0; i < height; i++) {
		for (s32 j = 0; j < width; j++) {
			if ((*data & 0xc0) != 0xc0)
				*pack++ = *data++;
			else {
				*pack++ = 0xc1;
				*pack++ = *data++;
			}
		}
		data += rowbytes - width;
	}
	*pack++ = 0x0c; // palette ID u8
	for (s32 i = 0; i < 768; i++) // write the palette
		*pack++ = *palette++;
	s32 length = pack - (u8 *) pcx; // write output file 
	COM_WriteFile(filename, pcx, length);
}


void WriteBMPfile(s8 *filename, u8 *data, s32 width, s32 height,
				  s32 rowbytes, u8 *palette)
{
	SDL_Surface *surface = SDL_CreateSurfaceFrom(
		width,
		height,
		SDL_PIXELFORMAT_INDEX8,
		data,
		rowbytes
	);

	SDL_Palette *surfacepalette = SDL_CreateSurfacePalette(surface);

	for (s32 i = 0; i < surfacepalette->ncolors; i++)
	{
		surfacepalette->colors[i].r = palette[i * 3 + 0];
		surfacepalette->colors[i].g = palette[i * 3 + 1];
		surfacepalette->colors[i].b = palette[i * 3 + 2];
		surfacepalette->colors[i].a = 255;
	}

	size_t len_temp = width * height * 2 + 1000;
	void *temp = Hunk_TempAlloc(len_temp);

	SDL_IOStream *io = SDL_IOFromMem(temp, len_temp);

	SDL_SaveBMP_IO(surface, io, false);

	Sint64 length = SDL_TellIO(io);

	SDL_CloseIO(io);

	COM_WriteFile(filename, temp, length);

	SDL_DestroySurface(surface);
}

void SCR_ScreenShot_f()
{
	s8 bmpname[80];
	s8 checkname[MAX_OSPATH*2];
	strcpy(bmpname, "quake00.bmp"); // find a file name to save it to
	s32 i = 0;
	for (; i <= 99; i++) {
		bmpname[5] = i / 10 + '0';
		bmpname[6] = i % 10 + '0';
		sprintf(checkname, "%s/%s", com_gamedir, bmpname);
		if (Sys_FileTime(checkname) == -1)
			break; // file doesn't exist
	}
	if (i == 100) {
		Con_Printf("SCR_ScreenShot_f: Couldn't create a BMP file\n");
		return;
	}
	WriteBMPfile(bmpname, vid.buffer, vid.width, vid.height, vid.width,
		     vid_curpal); // save the bmp file
	snprintf(checkname, 95, "echo Wrote %s", bmpname);
	Cbuf_AddText(checkname);
}

void SCR_BeginLoadingPlaque()
{
	S_StopAllSounds(1);
	if (cls.state != ca_connected || cls.signon != SIGNONS)
		return;
	Con_ClearNotify(); // redraw with no console and the loading plaque
	scr_centertime_off = 0;
	scr_con_current = 0;
	scr_drawloading = 1;
	scr_fullupdate = 0;
	Sbar_Changed();
	SCR_UpdateScreen();
	scr_drawloading = 0;
	scr_disabled_for_loading = 1;
	scr_disabled_time = realtime;
	scr_fullupdate = 0;
}

void SCR_EndLoadingPlaque()
{
	scr_disabled_for_loading = 0;
	scr_fullupdate = 0;
	Con_ClearNotify();
}

void SCR_DrawNotifyString()
{
	s8 *start = scr_notifystring;
	s32 y = vid.height * 0.35;
	do {
		s32 l = 0;
		for (; l < 40; l++) // scan the width of the line
			if (start[l] == '\n' || !start[l])
				break;
		s32 x = (vid.width - l * 8 * uiscale) / 2;
		for (s32 j = 0; j < l; j++, x += 8 * uiscale)
			Draw_CharacterScaled(x, y, start[j], uiscale);
		y += 8 * uiscale;
		while (*start && *start != '\n')
			start++;
		if (!*start)
			break;
		start++; // skip the \n
	} while (1);
}

s32 SCR_ModalMessage(s8 *text) // Displays a text string in the center
{ // of the screen and waits for a Y or N keypress.  
	if (cls.state == ca_dedicated)
		return 1;
	scr_notifystring = text;
	scr_fullupdate = 0; // draw a fresh screen
	scr_drawdialog = 1;
	SCR_UpdateScreen();
	scr_drawdialog = 0;
	S_ClearBuffer(); // so dma doesn't loop current sound
	do {
		key_count = -1; // wait for a key down and up
		Sys_SendKeyEvents();
	} while (key_lastpress != 'y' && key_lastpress != 'n'
		 && key_lastpress != K_ESCAPE);
	scr_fullupdate = 0;
	SCR_UpdateScreen();
	return key_lastpress == 'y';
}

void SCR_UpdateScreen() // This is called every frame,
{ // and can also be called explicitly to flush text to the screen.
	static f32 oldscr_viewsize;
	if (scr_disabled_for_loading) {
		if (realtime - scr_disabled_time > 60) {
			scr_disabled_for_loading = 0;
			Con_Printf("load failed.\n");
		} else
			return;
	}
	if (cls.state == ca_dedicated || !scr_initialized || !con_initialized)
		return;
	if (scr_viewsize.value != oldscr_viewsize) {
		oldscr_viewsize = scr_viewsize.value;
		vid.recalc_refdef = 1;
	}
	if (oldfov != scr_fov.value) { // check for vid changes
		oldfov = scr_fov.value;
		vid.recalc_refdef = 1;
	}
	if (oldscreensize != scr_viewsize.value) {
		oldscreensize = scr_viewsize.value;
		vid.recalc_refdef = 1;
	}
	if (vid.recalc_refdef) { // something changed, so reorder the screen
		SCR_CalcRefdef();
	}
	// do 3D refresh drawing, and then update the screen
	Draw_TileClear(0, 0, vid.width, vid.height); // clear the entire screen
	Sbar_Changed();
	SCR_SetUpToDrawConsole();
	SCR_EraseCenterString();
	V_RenderView(); // stay mapped in for linear writes all the time
	if (scr_drawdialog) {
		Sbar_Draw();
		fadescreen = 1;
		drawlayer = lyr_centerprint.value;
		SCR_DrawNotifyString();
		drawlayer = lyr_main.value;
	} else if (scr_drawloading) {
		SCR_DrawLoading();
		Sbar_Draw();
	} else if (cl.intermission == 1 && key_dest == key_game) {
		Sbar_IntermissionOverlay();
	} else if (cl.intermission == 2 && key_dest == key_game) {
		Sbar_FinaleOverlay();
		SCR_CheckDrawCenterString();
	} else if (cl.intermission == 3 && key_dest == key_game) {
		SCR_CheckDrawCenterString();
	} else {
		SCR_DrawRam();
		SCR_DrawNet();
		SCR_DrawTurtle();
		SCR_DrawPause();
		SCR_CheckDrawCenterString();
		Sbar_Draw();
		SCR_DrawConsole();
		M_Draw();
		SCR_DrawFPS();
	}
	V_UpdatePalette();
	VID_Update();
}

void SCR_HUDStyle_f (cvar_t *cvar)
{ // Updates hudstyle variable and invalidates refdef when scr_hudstyle changes
    s32 val = (s32) cvar->value;
    hudstyle = (hudstyle_t) CLAMP (0, val, (s32) HUD_COUNT - 1);
    vid.recalc_refdef = 1;
}
