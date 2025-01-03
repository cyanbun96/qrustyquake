#include "quakedef.h"

char *basedir = ".";
FILE *sys_handles[MAX_HANDLES];

// =============================================================================
// General routines
// =============================================================================

void Sys_Printf(char *fmt, ...)
{
	va_list argptr;
	char text[1024];
	va_start(argptr, fmt);
	vsprintf(text, fmt, argptr);
	va_end(argptr);
	fprintf(stderr, "%s", text);
}

void Sys_Quit()
{
	Host_Shutdown();
	exit(0);
}

void Sys_Error(char *error, ...)
{
	va_list argptr;
	char string[1024];
	va_start(argptr, error);
	vsprintf(string, error, argptr);
	va_end(argptr);
	fprintf(stderr, "Error: %s\n", string);
	Host_Shutdown();
	exit(1);
}

void Sys_DebugLog(char *file, char *fmt, ...)
{
	va_list argptr;
	static char data[1024];
	va_start(argptr, fmt);
	vsprintf(data, fmt, argptr);
	va_end(argptr);
	FILE *fp = fopen(file, "a");
	fwrite(data, strlen(data), 1, fp);
	fclose(fp);
}

// =============================================================================
// File IO
// =============================================================================

int findhandle()
{
	for (int i = 1; i < MAX_HANDLES; i++)
		if (!sys_handles[i])
			return i;
	Sys_Error("out of handles");
	return -1;
}

static int Qfilelength(FILE *f)
{
	int pos = ftell(f);
	fseek(f, 0, SEEK_END);
	int end = ftell(f);
	fseek(f, pos, SEEK_SET);
	return end;
}

int Sys_FileOpenRead(char *path, int *hndl)
{
	int i = findhandle();
	FILE *f = fopen(path, "rb");
	if (!f) {
		*hndl = -1;
		return -1;
	}
	sys_handles[i] = f;
	*hndl = i;
	return Qfilelength(f);
}

int Sys_FileOpenWrite(char *path)
{
	int i = findhandle();
	FILE *f = fopen(path, "wb");
	if (!f)
		Sys_Error("Error opening %s: %s", path, strerror(errno));
	sys_handles[i] = f;
	return i;
}

void Sys_FileClose(int handle)
{
	fclose(sys_handles[handle]);
	sys_handles[handle] = NULL;
}

void Sys_FileSeek(int handle, int position)
{
	fseek(sys_handles[handle], position, SEEK_SET);
}

int Sys_FileRead(int handle, void *dst, int count)
{
	return fread(dst, 1, count, sys_handles[handle]);
}

int Sys_FileWrite(int handle, void *src, int count)
{
	return fwrite(src, 1, count, sys_handles[handle]);
}

int Sys_FileTime(char *path)
{
	FILE *f = fopen(path, "rb");
	if (f) {
		fclose(f);
		return 1;
	}
	return -1;
}

void Sys_mkdir(char *path)
{
#ifdef _WIN32
	_mkdir(path);
#else
	mkdir(path, 0777);
#endif
}

double Sys_FloatTime()
{
#ifdef _WIN32
	static int starttime = 0;
	if (!starttime)
		starttime = clock();
	return (clock() - starttime) * 1.0 / 1024;
#else
	struct timeval tp;
	struct timezone tzp;
	static int secbase;
	gettimeofday(&tp, &tzp);
	if (!secbase) {
		secbase = tp.tv_sec;
		return tp.tv_usec / 1000000.0;
	}
	return (tp.tv_sec - secbase) + tp.tv_usec / 1000000.0;
#endif
}

byte *Sys_ZoneBase(int *size)
{
	char *QUAKEOPT = getenv("QUAKEOPT");
	*size = 0xc00000;
	if (QUAKEOPT) {
		while (*QUAKEOPT)
			if (tolower(*QUAKEOPT++) == 'm') {
				*size = atof(QUAKEOPT) * 1024 * 1024;
				break;
			}
	}
	return malloc(*size);
}

#ifdef _WIN32
// Function to parse lpCmdLine into argc and argv
void CommandLineToArgv(const char *lpCmdLine, int *argc, char ***argv)
{
	*argc = 0;
	char *cmdLine = strdup(lpCmdLine);
	char *token = strtok(cmdLine, " ");
	while (token) {
		(*argc)++;
		token = strtok(NULL, " ");
	}
	*argv = (char **)malloc((*argc + 1) * sizeof(char *));
	strcpy(cmdLine, lpCmdLine);
	token = strtok(cmdLine, " ");
	for (int i = 0; i < *argc; i++) {
		(*argv)[i] = strdup(token);
		token = strtok(NULL, " ");
	}
	(*argv)[*argc] = NULL;
	free(cmdLine);
}

int WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine,
	    int nCmdShow)
{
	int argc;
	char **argv;
	CommandLineToArgv(lpCmdLine, &argc, &argv);
	int retCode = main(argc, argv);
	for (int i = 0; i < argc; i++) {
		free(argv[i]);
	}
	free(argv);
	return retCode;
}
#endif

int main(int c, char **v)
{
	double time, oldtime, newtime;
	quakeparms_t parms;
	extern int vcrFile;
	extern int recording;

	if (SDL_Init(SDL_INIT_EVERYTHING) < 0)
		Sys_Error("SDL_Init failed: %s", SDL_GetError());

	signal(SIGFPE, SIG_IGN);
	parms.memsize = DEFAULT_MEMORY;
	if (COM_CheckParm("-heapsize")) {
		int t = COM_CheckParm("-heapsize") + 1;
		if (t < c)
			parms.memsize = Q_atoi(v[t]) * 1024;
	}
	parms.membase = malloc(parms.memsize);
	parms.basedir = basedir;
	// Disable cache, else it looks in the cache for config.cfg.
	parms.cachedir = NULL;

	COM_InitArgv(c, v);
	parms.argc = com_argc;
	parms.argv = com_argv;

	Host_Init(&parms);

	oldtime = Sys_FloatTime() - 0.1;
	while (1) {
		newtime = Sys_FloatTime();
		// find time spent rendering last frame
		time = newtime - oldtime;

		if (cls.state == ca_dedicated) {
			// play vcrfiles at max speed
			if (time < sys_ticrate.value
			    && (vcrFile == -1 || recording)) {
				SDL_Delay(1);
				// not time to run a server only tic yet
				continue;
			}
			time = sys_ticrate.value;
		}

		if (time > sys_ticrate.value * 2)
			oldtime = newtime;
		else
			oldtime += time;

		Host_Frame(time);
	}
}
