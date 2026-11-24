#include "quakedef.h"
static FILE *sys_handles[MAX_HANDLES];
#ifdef SDL_PLATFORM_ANDROID
#include <SDL3/SDL_main.h> // main becomes SDL_main, run by SDLActivity
static jmp_buf sys_exitjmp; // Android wants main to return, not exit()
static s32 sys_exitcode;
static SDL_ThreadID sys_mainthread;
static void Sys_Exit(s32 code)
{
	if(SDL_GetCurrentThreadID() != sys_mainthread) exit(code);
	sys_exitcode = code;
	longjmp(sys_exitjmp, 1);
}
#else
#define Sys_Exit exit
#endif

void Sys_Printf(const c8 *fmt, ...)
{
	va_list argptr;
	c8 qtext[MAXPRINTMSG];
	c8 u8text[MAXPRINTMSG*4];
	va_start(argptr, fmt);
	vsnprintf(qtext, sizeof(qtext), fmt, argptr);
	va_end(argptr);
	UTF8_FromQuake(u8text, sizeof(u8text), qtext);
#ifdef SDL_PLATFORM_ANDROID
	SDL_Log("%s", u8text); // to logcat, stdout goes nowhere
#else
	printf("%s", u8text);
#endif
}

void Sys_Quit()
{
	CDAudio_Stop();
	bool endscreen = !exitstyle.value;
#ifdef SDL_PLATFORM_ANDROID // in the background or being closed, unseen
	if(SDL_GetWindowFlags(window) & SDL_WINDOW_MINIMIZED) endscreen = 0;
#endif
	if(endscreen){
		u16 *screen; // erysdren(it/its)
		if(registered.value)screen=(u16*)COM_LoadHunkFile("end2.bin",0);
		else screen = (u16*)COM_LoadHunkFile("end1.bin", 0);
		SDL_SetWindowSize(window, 640, 400);
		if(screen) vgatext_main(window, screen);
	}
	Host_Shutdown();
#ifdef __EMSCRIPTEN__
	emscripten_cancel_main_loop();
#else
	Sys_Exit(0);
#endif
}

void Sys_Error(const c8 *error, ...)
{
	va_list argptr;
	c8 str[1024];
	CDAudio_Stop();
	va_start(argptr, error);
	vsnprintf(str, sizeof(str), error, argptr);
	va_end(argptr);
	fprintf(stderr, "Error: %s\n", str);
	SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "QrustyQuake", str, 0);
	Host_Shutdown();
#ifdef SDL_PLATFORM_DOS // again, as going back to text mode cleared the screen
	fprintf(stderr, "Error: %s\n", str);
#endif
#ifdef __EMSCRIPTEN__
	emscripten_cancel_main_loop();
#else
	Sys_Exit(1);
#endif
}

s32 findhandle()
{
	for(s32 i = 1; i < MAX_HANDLES; i++)
		if(!sys_handles[i]) return i;
	Sys_Error("out of handles");
	return -1;
}

static s32 Qfilelength(FILE *f)
{
	s32 pos = ftell(f);
	fseek(f, 0, SEEK_END);
	s32 end = ftell(f);
	fseek(f, pos, SEEK_SET);
	return end;
}

static f64 Sys_WaitUntil(f64 endtime)
{
	static f64 estimate = 1e-3;
	static f64 mean = 1e-3;
	static f64 m2 = 0.0;
	static f64 count = 1.0;
	f64 now = Sys_DoubleTime();
	f64 before, observed, delta, stddev;
	endtime -= 1e-6; // allow finishing 1 microsecond earlier than requested
	while(now + estimate < endtime){
		before = now;
#ifdef __EMSCRIPTEN__ // SDL_Delay spins without SDL's sleeps, see main
		emscripten_sleep(1);
#else
		SDL_Delay(1);
#endif
		now = Sys_DoubleTime();
	// Determine Sleep(1) mean duration & variance using Welford's algorithm
	// https://blog.bearcats.nl/accurate-sleep-function/
		if(count >= 1e6) continue;
		// skip this if we already have more than enough samples
		++count;
		observed = now - before;
		delta = observed - mean;
		mean += delta / count;
		m2 += delta * (observed - mean);
		stddev = sqrt(m2 / (count - 1.0));
		estimate = mean + 1.5 * stddev;
		estimate = q_min(estimate, 2e-3);
	}
	while(now < endtime){ now = Sys_DoubleTime(); }
	return now;
}

static f64 Sys_Throttle(f64 oldtime)
{ return Sys_WaitUntil(oldtime + Host_GetFrameInterval()); }

s32 Sys_FileOpenRead(const c8 *path, s32 *hndl)
{
	s32 i = findhandle();
	FILE *f = fopen(path, "rb");
	if(!f){ *hndl = -1; return -1; }
	sys_handles[i] = f;
	*hndl = i;
	return Qfilelength(f);
}


void Sys_FileClose(s32 handle)
{ fclose(sys_handles[handle]); sys_handles[handle] = NULL; }

void Sys_FileSeek(s32 handle, s32 position)
{ fseek(sys_handles[handle], position, SEEK_SET); }

s32 Sys_FileRead(s32 handle, void *dst, s32 count)
{ return fread(dst, 1, count, sys_handles[handle]); }

s32 Sys_FileTime(const c8 *path)
{
	FILE *f = fopen(path, "rb");
	if(f){ fclose(f); return 1; }
	return -1;
}

f64 Sys_DoubleTime()
{
	static u64 starttime = 0;
	if(!starttime) starttime = SDL_GetTicksNS();
	return((f64)(SDL_GetTicksNS() - starttime)) / 1000000000;
}

s32 Sys_FileType(const c8* path)
{
	SDL_PathInfo info;
	if(SDL_GetPathInfo(path, &info)){
		switch(info.type){
		case SDL_PATHTYPE_FILE: return FS_ENT_FILE;
		case SDL_PATHTYPE_DIRECTORY: return FS_ENT_DIRECTORY;
		default: break;
		}
	}
	return 0;
}

#ifdef __EMSCRIPTEN__
static void main_loop()
{ // the browser calls it once a refresh, a timedemo runs most of one
	static f64 oldtime = 0;
	f64 start = Sys_DoubleTime();
	do{
		if(oldtime == 0) oldtime = Sys_DoubleTime() - 0.1;
		f64 newtime = Sys_Throttle(oldtime);
		f64 deltatime = newtime - oldtime;
		if(cls.state == ca_dedicated) deltatime = sys_ticrate.value;
		if(deltatime > sys_ticrate.value * 2) oldtime = newtime;
		else oldtime += deltatime;
		Host_Frame(deltatime);
	}while(cls.timedemo && Sys_DoubleTime() - start < 0.012);
}
#endif

#ifdef SDL_PLATFORM_DOS // SDL's main runs this one, after locking memory for
#include <SDL3/SDL_main.h> // its interrupt handlers and enabling near pointers
#endif
int main(int c, char **v)
{ // the identifier names the window to the desktop (Wayland app_id, X11
  // class), the same as AppImage/qrustyquake.desktop, whatever the binary
#ifdef SDL_PLATFORM_ANDROID
	sys_mainthread = SDL_GetCurrentThreadID();
	if(setjmp(sys_exitjmp)) return sys_exitcode;
#endif
	c8 version[16];
	snprintf(version, sizeof(version), "%.2f", VERSION);
	SDL_SetAppMetadata("QrustyQuake", version, "qrustyquake");
	SDL_SetAppMetadataProperty(SDL_PROP_APP_METADATA_TYPE_STRING, "game");
#ifdef __EMSCRIPTEN__ // SDL would also sleep at every present, but main_loop
	// already hands each frame back to the browser. Asyncify stays for the
	// sleeps that have to wait in place (SCR_ModalMessage, Sys_WaitUntil)
	SDL_SetHint(SDL_HINT_EMSCRIPTEN_ASYNCIFY, "0");
#endif
	if(!SDL_Init(SDL_INIT_VIDEO))
		Sys_Error("SDL_Init failed: %s", SDL_GetError());
	host_parms.argc = c;
	host_parms.argv = (c8**)v;
	COM_InitArgv(host_parms.argc, (c8**)host_parms.argv);
	// the first hunk segment, as in Ironwail; more are added when it fills
#ifdef SDL_PLATFORM_DOS // CWSDPMI backs a big one with its swap file, and
	host_parms.memsize = 24*1024*1024; // DOSBox died on 162 MB of that
#else
	host_parms.memsize = 384*1024*1024;
#endif
	if(COM_CheckParm("-heapsize")){
		s32 t = COM_CheckParm("-heapsize") + 1;
		if(t < c) host_parms.memsize = q_max(atoi(v[t]), 1024) * 1024;
	}
	host_parms.membase = malloc(host_parms.memsize);
	if(host_parms.membase == NULL){
		printf("WARNING: Allocation of requested %d bytes failed.\n",
				host_parms.memsize);
		do{ // the biggest heap a small machine can give, down to 16 MB
			host_parms.memsize = host_parms.memsize / 4 * 3;
			if(host_parms.memsize < 16*1024*1024)
				host_parms.memsize = 16*1024*1024;
			host_parms.membase = malloc(host_parms.memsize);
		}while(!host_parms.membase
				&& host_parms.memsize > 16*1024*1024);
		printf("WARNING: Using fallback %d byte heap size.\n",
				host_parms.memsize);
	}
	host_parms.basedir = ".";
	host_parms.userdir = ".";
	Host_Init();
#ifdef __EMSCRIPTEN__
	emscripten_set_main_loop(main_loop, 0, 1);
#else
	f64 oldtime = Sys_DoubleTime() - 0.1;
	while(1){
		f64 newtime = Sys_Throttle(oldtime);
		f64 deltatime = newtime - oldtime;
		if(cls.state == ca_dedicated) deltatime = sys_ticrate.value;
		if(deltatime > sys_ticrate.value * 2) oldtime = newtime;
		else oldtime += deltatime;
		Host_Frame(deltatime);
	}
#endif
}
