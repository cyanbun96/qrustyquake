#ifndef _WIN32
#include <stdatomic.h>
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
#ifndef _WIN32
#include <unistd.h>
#include <sys/param.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <sys/ipc.h>
#ifndef __HAIKU__
#include <sys/shm.h>
#endif
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <sys/mman.h>
#else
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
