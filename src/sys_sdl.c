#include "quakedef.h"
static FILE *sys_handles[MAX_HANDLES];

void Sys_Printf(const s8 *fmt, ...)
{
	va_list argptr;
	s8 text[1024];
	va_start(argptr, fmt);
	vsprintf(text, fmt, argptr);
	va_end(argptr);
	fprintf(stderr, "%s", text);
}

void Sys_Quit()
{
	CDAudio_Stop();
        Uint16 *screen; // erysdren (it/its)
        if(registered.value) screen = (Uint16 *)COM_LoadHunkFile("end2.bin", 0);
        else screen = (Uint16 *)COM_LoadHunkFile("end1.bin", 0);
	SDL_SetWindowSize(window, 640, 400);
        if(screen) vgatext_main(window, screen);
	Host_Shutdown();
#ifdef __EMSCRIPTEN__
	emscripten_cancel_main_loop();
#else
	exit(0);
#endif
}

void Sys_Error(const s8 *error, ...)
{
	va_list argptr;
	s8 str[1024];
	CDAudio_Stop();
	va_start(argptr, error);
	vsprintf(str, error, argptr);
	va_end(argptr);
	fprintf(stderr, "Error: %s\n", str);
#ifdef DEBUG
	__asm__("int3");
#endif
	SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "QrustyQuake", str, 0);
	Host_Shutdown();
#ifdef __EMSCRIPTEN__
	emscripten_cancel_main_loop();
#else
	exit(1);
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

static f64 Sys_WaitUntil (f64 endtime)
{
	static f64 estimate = 1e-3;
	static f64 mean = 1e-3;
	static f64 m2 = 0.0;
	static f64 count = 1.0;
	f64 now = Sys_DoubleTime();
	f64 before, observed, delta, stddev;
	endtime -= 1e-6; // allow finishing 1 microsecond earlier than requested
	while(now + estimate < endtime) {
		before = now;
		SDL_Delay(1);
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
	while(now < endtime) { now = Sys_DoubleTime (); }
	return now;
}

static f64 Sys_Throttle (f64 oldtime)
{ return Sys_WaitUntil (oldtime + Host_GetFrameInterval ()); }

s32 Sys_FileOpenRead(const s8 *path, s32 *hndl)
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

s32 Sys_FileWrite(s32 handle, const void *src, s32 count)
{ return fwrite(src, 1, count, sys_handles[handle]); }

s32 Sys_FileTime(const s8 *path)
{
	FILE *f = fopen(path, "rb");
	if(f){ fclose(f); return 1; }
	return -1;
}

s32 Sys_FileOpenWrite(const s8 *path)
{
	s32 i = findhandle();
	FILE *f = fopen(path, "wb");
	if(!f) Sys_Error("Error opening %s: %s", path, strerror(errno));
	sys_handles[i] = f;
	return i;
}

#ifndef _WIN32
s32 Sys_FileType (const s8 *path)
{
        struct stat st;
        if(stat(path, &st) != 0) return FS_ENT_NONE;
        if(S_ISDIR(st.st_mode)) return FS_ENT_DIRECTORY;
        if(S_ISREG(st.st_mode)) return FS_ENT_FILE;
        return FS_ENT_NONE;
}

f64 Sys_DoubleTime()
{
	struct timeval tp;
	struct timezone tzp;
	static s32 secbase;
	gettimeofday(&tp, &tzp);
	if(!secbase){
		secbase = tp.tv_sec;
		return tp.tv_usec / 1000000.0;
	}
	return (tp.tv_sec - secbase) + tp.tv_usec / 1000000.0;
}

void Sys_mkdir(const s8 *path) { mkdir(path, 0777); }
#else

s32 Sys_FileType (const s8 *path)
{
        s32 result = GetFileAttributes(path);
        if(result == -1) return FS_ENT_NONE;
        if(result & FILE_ATTRIBUTE_DIRECTORY) return FS_ENT_DIRECTORY;
        return FS_ENT_FILE;
}

f64 Sys_DoubleTime()
{
	static s32 starttime = 0;
	if(!starttime) starttime = clock();
	return (clock() - starttime) * 1.0 / 1024;
}

void Sys_mkdir(const s8 *path) { _mkdir(path); }
#endif

#ifdef __EMSCRIPTEN__
static void main_loop(void)
{
	static f64 oldtime = 0;
	if (oldtime == 0) oldtime = Sys_DoubleTime() - 0.1;
	f64 newtime = Sys_Throttle(oldtime);
	f64 deltatime = newtime - oldtime;
	if(cls.state == ca_dedicated) deltatime = sys_ticrate.value;
	if(deltatime > sys_ticrate.value * 2) oldtime = newtime;
	else oldtime += deltatime;
	Host_Frame(deltatime);
}
#endif

int main(int c, char **v)
{
	if(!SDL_Init(SDL_INIT_VIDEO))
		Sys_Error("SDL_Init failed: %s", SDL_GetError());
	host_parms.argc = c;
	host_parms.argv = (s8**)v;
	COM_InitArgv(host_parms.argc, (s8**)host_parms.argv);
	host_parms.memsize = DEFAULT_MEMORY;
	if(COM_CheckParm("-heapsize")){
		s32 t = COM_CheckParm("-heapsize") + 1;
		if(t < c) host_parms.memsize = Q_atoi(v[t]) * 1024;
	}
	host_parms.membase = malloc(host_parms.memsize);
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
