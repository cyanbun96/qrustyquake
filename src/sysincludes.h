#include <SDL3/SDL_platform_defines.h> // for SDL_PLATFORM_* macros
#ifndef SDL_PLATFORM_WINDOWS
#include <stdatomic.h>
#include <limits.h> // for PATH_MAX, usually 4096
#endif
#include <errno.h>
#include <stddef.h>
#include <limits.h>
#include <stdbool.h>
#include <signal.h>
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
#include <assert.h>
#ifndef SDL_PLATFORM_WINDOWS
#include <dirent.h>
#include <unistd.h>
#include <sys/param.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <sys/ipc.h>
#ifndef SDL_PLATFORM_HAIKU
#include <sys/shm.h>
#endif
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <sys/mman.h>
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
#include <SDL3/SDL.h>
#ifdef AVAIL_SDL3MIXER
#include <SDL3_mixer/SDL_mixer.h>
#endif
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif
