#include <SDL3/SDL_platform_defines.h> // for SDL_PLATFORM_* macros
#include <errno.h>
#include <stddef.h>
#include <limits.h> // for PATH_MAX, usually 4096
#include <stdbool.h>
#include <stdlib.h>
#include <fcntl.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <setjmp.h>
#include <float.h>
#include <time.h>
#include <sys/types.h>
#ifndef SDL_PLATFORM_WINDOWS
#include <unistd.h>
#include <sys/param.h>
#include <sys/ioctl.h>
#ifndef SDL_PLATFORM_DOS // DJGPP has no sockets, loopback only (net_bsd.c)
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#endif
#else // SDL_PLATFORM_WINDOWS
#define WIN32_LEAN_AND_MEAN
#ifndef _USE_WINSOCK2
#define _USE_WINSOCK2
#endif
#include <winsock2.h>
#include <ws2tcpip.h>
#include <direct.h>
#include <io.h>
#endif
#ifdef SDL_PLATFORM_DOS // DJGPP has no fenv.h; the rounding r_main.c and
#define FE_TONEAREST RC_NEAR // r_thread.c set is the x87 control word's
#define FE_TOWARDZERO RC_CHOP
typedef unsigned int fenv_t;
#define fesetround(r) _control87((r), MCW_RC)
#define fegetenv(e) (*(e) = _control87(0, 0))
#define fesetenv(e) _control87(*(e), 0xffff)
// nor C99 fmin and fmax (model.c)
static inline double fmin(double a, double b){return b!=b || a<b ? a : b;}
static inline double fmax(double a, double b){return b!=b || a>b ? a : b;}
#endif
#include <SDL3/SDL.h>
#ifdef AVAIL_SDL3MIXER
#include <SDL3_mixer/SDL_mixer.h>
#endif
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif
