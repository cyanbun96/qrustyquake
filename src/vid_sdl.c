#include "quakedef.h"

static SDL_Renderer *renderer;
static SDL_Surface *argbbuffer;
static SDL_Texture *texture;
static SDL_Rect blitRect;
static u8 *screenpixels;
static u8 *toppixels;
static u8 *uipixels;
static u8 *sbarpixels;
static SDL_PixelFormat window_format;
static SDL_Palette *sdlworldpal;
static SDL_Palette *sdltoppal;
static SDL_Palette *sdluipal;
static SDL_Palette *sdlsbarpal;
static s32 renderscale;
static s32 scalemode;

void VID_CalcScreenDimensions(cvar_t *cvar);
void VID_AllocBuffers();
void VID_VidFullscreenCommand_f();

s32 VID_GetConfigCvar(FILE *file, const c8 *cvname)
{
	// CyanBun96: _vid_default_mode_win gets read from config.cfg only after
	// video is initialized. To avoid creating a window and textures just to
	// destroy them right away and only then replace them with valid ones,
	// this function reads the default mode independently of the cvar status
	if(!file) return -2;
	c8 line[256];
	s32 ret = -1, len = strlen(cvname);
	rewind(file);
	while(fgets(line, sizeof(line), file)){
		c8 *p = line;
		while(*p == ' ' || *p == '\t') ++p;
		if(strncmp(p, cvname, len) || (p[len] != ' ' && p[len] != '\t'))
			continue; // not this cvar, or a bind/alias mentioning it
		c8 *start = strchr(p + len, '"');
		if(start){
			c8 *end = strchr(start + 1, '"');
			if(end){
				*end = '\0';
				ret = atoi(start + 1);
				break;
			}
		}
	}
	if(ret != -1) printf("%s: %d\n", cvname, ret);
	else printf("%s not found in config.cfg\n", cvname);
	return ret;
}

s32 VID_DetermineMode(s32 win)
{
	for(s32 i = 0; i < NUM_OLDMODES; ++i){
		if((i > 2 ? !win : win) && vid.width == oldmodes[i * 2]
				&& vid.height == oldmodes[i * 2 + 1])
			return i;
	}
	return -1;
}

void VID_SetPalette(u8 *palette, SDL_Surface *dest)
{
	SDL_Color colors[256];
	for(s32 i = 0; i < 256; ++i){
		colors[i].r = *palette++;
		colors[i].g = *palette++;
		colors[i].b = *palette++;
		colors[i].a = 255;
	}
	SDL_Palette *pal = sdlworldpal;
	if(dest == screentop) pal = sdltoppal;
	else if(dest == screenui) pal = sdluipal;
	else if(dest == screensbar) pal = sdlsbarpal;
	SDL_SetPaletteColors(pal, colors, 0, 256);
	SDL_SetSurfacePalette(dest, pal);
}

void VID_SetWindowed()
{
	SDL_SetWindowFullscreen(window, 0);
	SDL_SetWindowSize(window,vid.width*renderscale,vid.height*renderscale);
	SDL_SetWindowPosition(window, SDL_WINDOWPOS_UNDEFINED,
			SDL_WINDOWPOS_UNDEFINED);
}

void VID_SetFullscreen(s32 w, s32 h)
{
	s32 moden = -1;
	SDL_DisplayMode *mode = 0;
	SDL_DisplayMode **modes = SDL_GetFullscreenDisplayModes(
			SDL_GetDisplayForWindow(window), &moden);
	for(s32 i = 0; i < moden; ++i){
		if(modes[i]->w == w && modes[i]->h == h){
			mode = modes[i];
			break;
		}}
	if(!mode)Con_Printf("Couldn't find an exclusive fullscreen mode with specified dimensions(%dx%d), using borderless fullscreen\n", w, h);
	else Con_Printf("Setting exclusive fullscreen mode(%dx%d)\n", w, h);
	SDL_SetWindowFullscreenMode(window, mode);
	SDL_SetWindowFullscreen(window, 1);
	SDL_free(modes);
}

void VID_SetBorderless()
{
	SDL_SetWindowFullscreenMode(window, 0);
	SDL_SetWindowFullscreen(window, 1);
}

void VID_UpdateUIScale(SDL_UNUSED cvar_t *cvar)
{
	if(vid.width / 320 >= scr_uiscale.value && scr_uiscale.value >= 1){
		uiscale = scr_uiscale.value;
		vid.recalc_refdef = 1;
	}else{
		Con_Printf("Invalid UI scale value %d\n",
				(s32)scr_uiscale.value);
		Con_Printf("Max %d for current resolution\n", vid.width/320);
		Cvar_SetValue("scr_uiscale", uiscale);
	}
}

void VID_Init(SDL_UNUSED u8 *palette)
{
	c8 caption[50];
	Cvar_RegisterVariable(&_windowed_mouse);
	Cvar_RegisterVariable(&_vid_default_mode_win);
	Cvar_RegisterVariable(&vid_mode);
	Cvar_RegisterVariable(&scr_uiscale);
	Cvar_RegisterVariable(&sensitivityyscale);
	Cvar_RegisterVariable(&newoptions);
	Cvar_RegisterVariable(&aspectr);
	Cvar_RegisterVariable(&aspectr_lock);
	Cvar_RegisterVariable(&realwidth);
	Cvar_RegisterVariable(&realheight);
	void(*vid_callback)(cvar_t *) =
		(void(*)(cvar_t *))VID_CalcScreenDimensions;
	Cvar_SetCallback(&aspectr, vid_callback);
	Cvar_SetCallback(&realwidth, vid_callback);
	Cvar_SetCallback(&realheight, vid_callback);
	Cvar_SetCallback(&scr_uiscale, VID_UpdateUIScale);
	Cmd_AddCommand("vid_fullscreen", VID_VidFullscreenCommand_f);
	// Set up display mode(width and height)
	vid.width = 320;
	vid.height = 240;
	s32 winmode = 0;
#ifdef __EMSCRIPTEN__
	// erysdren: emscripten requires windowed mode
	winmode = 1;
#endif
	c8 cfgpath[MAX_OSPATH];
	snprintf(cfgpath, sizeof(cfgpath), "%s/config.cfg", com_gamedir);
	FILE *cfg = fopen(cfgpath, "r");
	if(!cfg) printf("VID_Init: Failed to open %s\n", cfgpath);
	s32 defmode = -1;
	if(!(COM_CheckParm("-width") || COM_CheckParm("-height")
	   || COM_CheckParm("-window") || COM_CheckParm("-fullscreen")
	   || COM_CheckParm("-borderless") || COM_CheckParm("-winsize"))){
		defmode = VID_GetConfigCvar(cfg, "_vid_default_mode_win");
		if(defmode >= 0 && defmode < NUM_OLDMODES){
			vid.width = oldmodes[defmode * 2];
			vid.height = oldmodes[defmode * 2 + 1];
			winmode = defmode < 3;
		}
	}
	s32 confwidth = VID_GetConfigCvar(cfg, "vid_cwidth");
	s32 confheight = VID_GetConfigCvar(cfg, "vid_cheight");
	s32 confwmode = VID_GetConfigCvar(cfg, "vid_cwmode");
	renderscale = VID_GetConfigCvar(cfg, "r_renderscale");
	if(cfg) fclose(cfg);
	if(renderscale < 1) renderscale = 1;
	if(confwidth >= 320 && confheight >= 200 &&
		confwidth <= MAXWIDTH && confheight <= MAXHEIGHT){
		vid.width = confwidth;
		vid.height = confheight;
		Sys_Printf("Using vid_cwidth x vid_cheight to set mode\n");
	}
	else if(defmode != -1)
		Sys_Printf("Using _vid_default_mode_win to set mode\n");
	s32 pnum;
	if((pnum = COM_CheckParm("-winsize"))){
		if(pnum >= com_argc - 2)
			Sys_Error("VID: -winsize <width> <height>\n");
		vid.width = atoi(com_argv[pnum + 1]);
		vid.height = atoi(com_argv[pnum + 2]);
		if(!vid.width || !vid.height)
			Sys_Error("VID: Bad window width/height\n");
	}
	if((pnum = COM_CheckParm("-width"))){
		if(pnum >= com_argc - 1)
			Sys_Error("VID: -width <width>\n");
		vid.width = atoi(com_argv[pnum + 1]);
		if(!vid.width)
			Sys_Error("VID: Bad window width\n");
	}
	if((pnum = COM_CheckParm("-height"))){
		if(pnum >= com_argc - 1)
			Sys_Error("VID: -height <height>\n");
		vid.height = atoi(com_argv[pnum + 1]);
		if(!vid.height)
			Sys_Error("VID: Bad window height\n");
	}
	vimmode = COM_CheckParm("-vimmode");
	if(vid.width > 1280 || vid.height > 1024)
		Sys_Printf("vanilla maximum resolution is 1280x1024\n");
	aspectr.value = 1.333333;
	realwidth.value = vid.width*renderscale;
	realheight.value = (s32)(vid.width / aspectr.value + 0.5) * renderscale;
	window = SDL_CreateWindow("QrustyQuake",realwidth.value,
			realheight.value,SDL_WINDOW_RESIZABLE);
	if(!window) Sys_Error("Couldn't create window: %s", SDL_GetError());
	SDL_SetWindowMinimumSize(window, 320, 200);
	screen = SDL_CreateSurfaceFrom(vid.width, vid.height,
			SDL_PIXELFORMAT_INDEX8, NULL, vid.width);
		sdlworldpal = SDL_CreateSurfacePalette(screen);
	screensbar = SDL_CreateSurfaceFrom(vid.width, vid.height,
			SDL_PIXELFORMAT_INDEX8, NULL, vid.width);
		sdlsbarpal = SDL_CreateSurfacePalette(screensbar);
	screenui = SDL_CreateSurfaceFrom(vid.width, vid.height,
			SDL_PIXELFORMAT_INDEX8, NULL, vid.width);
		sdluipal = SDL_CreateSurfacePalette(screenui);
	screentop = SDL_CreateSurfaceFrom(vid.width, vid.height,
			SDL_PIXELFORMAT_INDEX8, NULL, vid.width);
		sdltoppal = SDL_CreateSurfacePalette(screentop);
	screen->pixels = screenpixels;
	screenui->pixels = uipixels;
	screensbar->pixels = sbarpixels;
	screentop->pixels = toppixels;
	screen->pitch = vid.width;
	screenui->pitch = vid.width;
	screensbar->pitch = vid.width;
	screentop->pitch = vid.width;
	vid.buffer = screenpixels;
	scrbuffs[0] = screen; scrbuffs[1] = screentop;
	scrbuffs[2] = screenui; scrbuffs[3] = screensbar;
	VID_SetPalette(host_basepal, screen);
	VID_SetPalette(host_basepal, screentop);
	SDL_SetSurfaceColorKey(screentop, 1, 255);
	VID_SetPalette(host_basepal, screenui);
	SDL_SetSurfaceColorKey(screenui, 1, 255);
	VID_SetPalette(host_basepal, screensbar);
	SDL_SetSurfaceColorKey(screensbar, 1, 255);
	renderer = SDL_CreateRenderer(window, NULL);
	if(!renderer)
		Sys_Error("Couldn't create renderer: %s", SDL_GetError());
	window_format = SDL_GetWindowPixelFormat(window);
	// pixels and pitch are filled in by SDL_LockTexture every frame
	argbbuffer = SDL_CreateSurfaceFrom(vid.width*renderscale,
			vid.height*renderscale, window_format, 0, 0);
	texture = SDL_CreateTexture(renderer, window_format,
			SDL_TEXTUREACCESS_STREAMING, vid.width*renderscale,
			vid.height * renderscale);
	if(!texture) Sys_Error("Couldn't create texture: %s", SDL_GetError());
	SDL_SetTextureScaleMode(texture, scalemode);
	sprintf(caption, "QrustyQuake - Version %4.2f", VERSION);
	SDL_SetWindowTitle(window, (const c8 *)&caption);
	vid.aspect = ((f32)vid.height / (f32)vid.width) * (320.0 / 240.0);
	vid.numpages = 1;
	vid.colormap = host_colormap;
	VID_AllocBuffers(); // allocate z buffer, surface cache and the fbuffer
	s32 win = 0;
	if(COM_CheckParm("-fullscreen") || confwmode == 1)
		VID_SetFullscreen(vid.width*renderscale,vid.height*renderscale);
	else if(COM_CheckParm("-borderless")||confwmode==2)VID_SetBorderless();
	else if(COM_CheckParm("-window") || confwmode == 0 || winmode){
		VID_SetWindowed();
		win = 1;
	}else VID_SetBorderless();
	if(defmode >= 0) vid_modenum = defmode;
	else vid_modenum = VID_DetermineMode(win);
	if(vid_modenum < 0) Con_Printf("WARNING: non-standard video mode\n");
	else Con_Printf("Detected video mode %d\n", vid_modenum);
	SDLWindowFlags = SDL_GetWindowFlags(window);
	realwidth.value = 0;
	VID_CalcScreenDimensions(0);
	memset(scrshot_name, 0, sizeof(scrshot_name));
}

void VID_Shutdown(){ SDL_QuitSubSystem(SDL_INIT_VIDEO); }

void VID_CalcScreenDimensions(SDL_UNUSED cvar_t *cvar)
{
	if(scr_lockuiscale.value == 0) uiscale = vid.width / 320;
	if(uiscale * 200 > vid.height)
		uiscale = 1; // For weird resolutions like 800x200
	Cvar_SetValue("scr_uiscale", uiscale);
	blitRect.x = 0;
	blitRect.y = 0;
	blitRect.w = vid.width;
	blitRect.h = vid.height;
	s32 winW, winH;
	if(realwidth.value == 0 || realheight.value == 0){
		SDL_GetWindowSize(window, &winW, &winH);
	}else{
		winW = realwidth.value;
		winH = realheight.value;
	}
	s32 bufW = vid.width * renderscale;
	s32 bufH = vid.height * renderscale;
	f32 bufAspect;
	if(aspectr.value == 0)
		bufAspect = (f32)bufW / bufH;
	else
		bufAspect = aspectr.value;
	Cvar_SetValue("aspectr", bufAspect);
	s32 destW, destH;
	if((f32)winW / winH > bufAspect){
		// Window is wider than buffer, black bars on sides
		destH = winH;
		destW = (s32)(winH * bufAspect + 0.5);
	}else{
		// Window is taller than buffer, black bars on top/bottom
		destW = winW;
		destH = (s32)(winW / bufAspect + 0.5);
	}
	destRect.x = (winW - destW) / 2; // Center the destination rectangle
	destRect.y = (winH - destH) / 2;
	destRect.w = destW;
	destRect.h = destH;
	Sbar_CalcPos();
}

void VID_Update()
{
	bool blitscreentop = false;
	bool blitscreenui = false;
	bool blitscreensbar = false;
	if(lyr_main.value==1||lyr_sbar.value==1||lyr_menu.value==1||
		lyr_centerprint.value==1||lyr_console.value==1||
		lyr_notify.value==1||lyr_crosshair.value==1)blitscreentop=true;
	if(lyr_main.value==2||lyr_sbar.value==2||lyr_menu.value==2||
		lyr_centerprint.value==2||lyr_console.value==2||
		lyr_notify.value==2||lyr_crosshair.value==2)blitscreenui=true;
	if(lyr_main.value==3||lyr_sbar.value==3||lyr_menu.value==3||
		lyr_centerprint.value==3||lyr_console.value==3||
		lyr_notify.value==3||lyr_crosshair.value==3)blitscreensbar=true;
	SDL_Rect dst = {0, 0, vid.width*renderscale, vid.height*renderscale};
	static SDL_Rect scRect;
	scRect.w = vid.width * (scr_uixscale.value>0 ? scr_uixscale.value : 1);
	scRect.h = vid.height * (scr_uiyscale.value>0 ? scr_uiyscale.value : 1);
	scRect.x = (vid.width - scRect.w) / 2;
	scRect.y = (vid.height - scRect.h) / 2;
	static SDL_Rect scRect2;
	scRect2.w = vid.width * renderscale * (scr_uixscale.value>0 ?
			scr_uixscale.value : 1);
	scRect2.h = vid.height * renderscale * (scr_uiyscale.value>0 ?
			scr_uiyscale.value : 1);
	scRect2.x = (vid.width*renderscale - scRect2.w) / 2;
	scRect2.y = (vid.height*renderscale - scRect2.h) / 2;
	if(SDL_LockTexture(texture,0, &argbbuffer->pixels, &argbbuffer->pitch)){
		if(blitscreensbar)
			SDL_BlitSurfaceScaled(screensbar, &blitRect, screen,
				&scRect, SDL_SCALEMODE_NEAREST); // 8-bit to 8-bit
		if(fadescreen == 1 && lyr_menu.value != 0) Draw_FadeScreen();
		SDL_BlitSurfaceScaled(screen, &blitRect, argbbuffer, &dst,
					scalemode);
		if(blitscreenui)
			SDL_BlitSurfaceScaled(screenui, &blitRect, argbbuffer,
					&scRect2, scalemode);
		if(blitscreentop)
			SDL_BlitSurfaceScaled(screentop, &blitRect, argbbuffer,
					&dst, scalemode);
		if(scrshot_name[0]){ // the pixels are only valid while locked
			if(!SDL_SaveBMP(argbbuffer, scrshot_name))
				Con_Printf("Screenshot: couldn't save %s: %s\n",
						scrshot_name, SDL_GetError());
			else
				Con_Printf("Wrote %s\n", scrshot_name);
			memset(scrshot_name, 0, sizeof(scrshot_name));
		}
		SDL_UnlockTexture(texture);
	}else{ printf("Couldn't lock texture %s\n", SDL_GetError()); }
	SDL_RenderClear(renderer);
	SDL_RenderTexture(renderer, texture, NULL, &destRect);
	SDL_RenderPresent(renderer);
	if(vid_testingmode){
		if(realtime >= vid_testendtime){
			VID_SetMode(vid_realmode, 0, 0, 0);
			vid_testingmode = 0;
		}
	}else if((s32)vid_mode.value != vid_realmode && vid_realmode != -1){
		VID_SetMode((s32)vid_mode.value, 0, 0, 0);
		Cvar_SetValue("vid_mode", (f32)vid_modenum);
		vid_realmode = vid_modenum;
	}
	// unused layers aren't drawn to, so they stay clear from AllocBuffers
	if(blitscreentop) memset(screentop->pixels, 255, vid.width*vid.height);
	if(blitscreensbar) memset(screensbar->pixels,255,vid.width*vid.height);
	if(blitscreenui) memset(screenui->pixels, 255, vid.width*vid.height);
}

c8 *VID_GetModeDescription(s32 mode)
{
	static c8 pinfo[40];
	if((mode < 0) || (mode >= NUM_OLDMODES))
		sprintf(pinfo, "Custom fullscreen");
	else if(mode >= 3)
		sprintf(pinfo, "%s fullscreen", modelist[mode]);
	else
		sprintf(pinfo, "%s windowed", modelist[mode]);
	return pinfo;
}

void VID_AllocBuffers()
{
	s32 area = vid.width * (vid.height + 1); // +safety padding line
	screenpixels = realloc(screenpixels, area);
	toppixels = realloc(toppixels, area);
	uipixels = realloc(uipixels, area);
	sbarpixels = realloc(sbarpixels, area);
	if(!screenpixels||!toppixels||!uipixels||!sbarpixels)
		Sys_Error("Not enough memory for video mode");
	memset(toppixels, 255, area); // 255 is transparent in the overlays
	memset(uipixels, 255, area);
	memset(sbarpixels, 255, area);
	screen->pixels = vid.buffer = screenpixels;
	screentop->pixels = toppixels;
	screenui->pixels = uipixels;
	screensbar->pixels = sbarpixels;
	d_pzbuffer_size = 0; // reallocate caches
	D_AllocCaches();
}

void VID_SetMode(s32 modenum, s32 custw, s32 custh, s32 custwinm)
{
	Con_DPrintf("SetMode: %d %dx%d %d\n", modenum, custw, custh, custwinm);
	if(modenum == vid_modenum && (!custw || !custh))
		return;
	if(custw && custh){
		vid.width = custw;
		vid.height = custh;
		vid_modenum = -1;
	}else{
		vid.width = oldmodes[modenum * 2];
		vid.height = oldmodes[modenum * 2 + 1];
		vid_modenum = modenum;
	}
	SDL_DestroySurface(screen);
	SDL_DestroySurface(screentop);
	SDL_DestroySurface(screenui);
	SDL_DestroySurface(screensbar);
	screen = SDL_CreateSurfaceFrom(vid.width, vid.height,
			SDL_PIXELFORMAT_INDEX8, NULL, vid.width);
		sdlworldpal = SDL_CreateSurfacePalette(screen);
	screensbar = SDL_CreateSurfaceFrom(vid.width, vid.height,
			SDL_PIXELFORMAT_INDEX8, NULL, vid.width);
		sdlsbarpal = SDL_CreateSurfacePalette(screensbar);
	screenui = SDL_CreateSurfaceFrom(vid.width, vid.height,
			SDL_PIXELFORMAT_INDEX8, NULL, vid.width);
		sdluipal = SDL_CreateSurfacePalette(screenui);
	screentop = SDL_CreateSurfaceFrom(vid.width, vid.height,
			SDL_PIXELFORMAT_INDEX8, NULL, vid.width);
		sdltoppal = SDL_CreateSurfacePalette(screentop);
	scrbuffs[0] = screen; scrbuffs[1] = screentop;
	scrbuffs[2] = screenui; scrbuffs[3] = screensbar;
	SDL_DestroySurface(argbbuffer);
	// pixels and pitch are filled in by SDL_LockTexture every frame
	argbbuffer = SDL_CreateSurfaceFrom(vid.width*renderscale,
			vid.height*renderscale, window_format, 0, 0);
	SDL_DestroyTexture(texture);
	texture = SDL_CreateTexture(renderer, window_format,
		SDL_TEXTUREACCESS_STREAMING, vid.width*renderscale,
		vid.height * renderscale);
	if(!texture) Sys_Error("Couldn't create texture: %s", SDL_GetError());
	SDL_SetTextureScaleMode(texture, scalemode);
	vid.aspect = ((f32)vid.height / (f32)vid.width) * (320.0 / 240.0);
	screen->pixels = screenpixels;
	screenui->pixels = uipixels;
	screensbar->pixels = sbarpixels;
	screentop->pixels = toppixels;
	screen->pitch = vid.width;
	screenui->pitch = vid.width;
	screensbar->pitch = vid.width;
	screentop->pitch = vid.width;
	vid.buffer = screen->pixels;
	VID_AllocBuffers();
	vid.recalc_refdef = 1;
	VID_SetPalette(worldpalname[0]?worldpal:host_basepal, screen);
	VID_SetPalette(uipalname[0]?uipal:host_basepal, screenui);
	SDL_SetSurfaceColorKey(screenui, 1, 255);
	VID_SetPalette(uipalname[0]?uipal:host_basepal, screensbar);
	SDL_SetSurfaceColorKey(screensbar, 1, 255);
	VID_SetPalette(host_basepal, screentop);
	SDL_SetSurfaceColorKey(screentop, 1, 255);
	if(!custw || !custh){
		if(modenum <= 2) VID_SetWindowed();
		else VID_SetBorderless();
	}else{
		if(custwinm == 0) VID_SetWindowed();
		else if(custwinm== 1) VID_SetFullscreen(custw*renderscale,
							custh*renderscale);
		else VID_SetBorderless();
	}
	realwidth.value = 0;
	VID_CalcScreenDimensions(0);
	Cvar_SetValue("vid_mode", (f32)vid_modenum);
	Cvar_SetValue("vid_cwidth", (f32)custw);
	Cvar_SetValue("vid_cheight", (f32)custh);
	Cvar_SetValue("vid_cwmode", (f32)custwinm);
}

void VID_VidSetModeCommand_f()
{
	s32 custw, custh; // keep here for OpenBSD
	switch(Cmd_Argc()){
	default:
	case 1:
		Con_Printf("usage:\n");
		Con_Printf("   vid_setmode <width> <height> <mode>\n");
		Con_Printf("   modes: 0 - windowed\n");
		Con_Printf("          1 - fullscreen\n");
		Con_Printf("          2 - borderless\n");
		return;
	case 4:
		custw = atoi(Cmd_Argv(1));
		custh = atoi(Cmd_Argv(2));
		if(custw<320||custw>MAXWIDTH||custh<200||custh>MAXHEIGHT){
			Con_Printf("320 <= width <= %d\n", MAXWIDTH);
			Con_Printf("200 <= height <= %d\n", MAXHEIGHT);
			return;
		}
		vid_realmode = -1;
		VID_SetMode(-1, custw, custh, atoi(Cmd_Argv(3)));
		break;
	}
}

void VID_VidFullscreenCommand_f()
{
	s32 m = 0;
	switch(Cmd_Argc()){
	default:
	case 1:
		Con_Printf("usage:\n");
		Con_Printf("   vid_fullscreen <mode>\n");
		Con_Printf("   modes: 0 - windowed\n");
		Con_Printf("          1 - fullscreen\n");
		Con_Printf("          2 - borderless\n");
		return;
	case 2:
		m = atoi(Cmd_Argv(1));
		if(!m)
			VID_SetWindowed();
		else if(m == 1)
			VID_SetFullscreen(vid.width * renderscale,
					vid.height * renderscale);
		else{
			VID_SetBorderless();
			m = 2;
		}
		Cvar_SetValue("vid_cwmode", m);
		break;
	}
}

s32 VID_GetWinMode()
{ // 0 - windowed, 1 - exclusive fullscreen, 2 - borderless
	if(!(SDL_GetWindowFlags(window) & SDL_WINDOW_FULLSCREEN)) return 0;
	return SDL_GetWindowFullscreenMode(window) ? 1 : 2;
}

void VID_SetRenderScaleCommand_f(SDL_UNUSED cvar_t *cvar)
{
	if((s32)r_renderscale.value < 1)
		r_renderscale.value = 1;
	if(renderscale == (s32)r_renderscale.value)
		return; // e.g. config.cfg, VID_Init already used it
	renderscale = r_renderscale.value;
	s32 winm = VID_GetWinMode(); // keep the current window mode
	VID_SetMode(-1, vid.width, vid.height, winm);
	vid_modenum = VID_DetermineMode(!winm);
	Cvar_SetValue("vid_mode", (f32)vid_modenum);
}

void VID_SetScaleModeCommand_f()
{
	s32 i;
	switch(Cmd_Argc()){
	default:
	case 1:
		Con_Printf("valid values:\n");
		Con_Printf("   0 - nearest(default, pixelated)\n");
		Con_Printf("   1 - linear(smooth)\n");
		Con_Printf("   2 - pixelart(pixelated, smoother)\n");
		Con_Printf("   current: %d\n", scalemode);
		return;
	case 2:
		i = atoi(Cmd_Argv(1));
		if(i < 0 || i > 2){
			Con_Printf("valid values:\n");
			Con_Printf("   0 - nearest(default, pixelated)\n");
			Con_Printf("   1 - linear(smooth)\n");
			Con_Printf("   2 - pixelart(pixelated, smoother)\n");
			Con_Printf("   current: %d\n", scalemode);
			return;
		}
		else{
			scalemode = i;
			SDL_SetTextureScaleMode(texture, scalemode);
		}
		break;
	}
}
