// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.

// quakedef.h -- primary header for client

#ifndef __QUAKEDEF__
#define __QUAKEDEF__
// CyanBun96: imagine if we could keep all the includes in one place instead
// of all over the .c and .h files that'd be cool right
#include <SDL2/SDL.h>
#ifdef __WIN32__
#include <windows.h>
#include <time.h>
#include <direct.h>
#include <io.h>
#else
#include <unistd.h>
#include <sys/ipc.h>
#ifndef __HAIKU__
#include <sys/shm.h>
#endif
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <sys/mman.h>
#endif
#include <stdbool.h>
#include <signal.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/types.h>
#include <fcntl.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <setjmp.h>

#define VERSION 1.09
#define FITZQUAKE_VERSION 0.85
#define GAMENAME "id1"
#define CMDLINE_LENGTH 256
#define DEFAULT_MEMORY (128 * 1024 * 1024)
#define CACHE_SIZE 32 // used to align key data structures
#define MINIMUM_MEMORY 0x550000
#define MINIMUM_MEMORY_LEVELPAK (MINIMUM_MEMORY + 0x100000)
#define MAX_NUM_ARGVS 50
#define PITCH 0 // up / down
#define YAW 1 // left / right
#define ROLL 2 // fall over
#define MAX_QPATH 64 // max length of a quake game pathname
#define MAX_OSPATH 128 // max length of a filesystem pathname
#define ON_EPSILON 0.1 // point on plane side epsilon
#define MAX_MSGLEN 65536 // max length of a reliable message, from Mark V
#define MAX_DATAGRAM 65527 // max length of unreliable message, from Mark V
#define MAX_EDICTS 32000 // FIXME: ouch! ouch! ouch!
#define MIN_EDICTS 256 // johnfitz -- lowest allowed value for max_edicts cvar
#define MAX_LIGHTSTYLES 64
#define MAX_MODELS 2048 // these are sent over the net as bytes
#define MAX_SOUNDS 2048 // so they cannot be blindly increased
#define SAVEGAME_COMMENT_LENGTH 39
#define MAX_STYLESTRING 64
#define MAX_CL_STATS 32 // stats are integers communicated to
#define STAT_HEALTH 0 // the client by the server
#define STAT_FRAGS 1
#define STAT_WEAPON 2
#define STAT_AMMO 3
#define STAT_ARMOR 4
#define STAT_WEAPONFRAME 5
#define STAT_SHELLS 6
#define STAT_NAILS 7
#define STAT_ROCKETS 8
#define STAT_CELLS 9
#define STAT_ACTIVEWEAPON 10
#define STAT_TOTALSECRETS 11
#define STAT_TOTALMONSTERS 12
#define STAT_SECRETS 13 // bumped on client side by svc_foundsecret
#define STAT_MONSTERS 14 // bumped by svc_killedmonster
#define IT_SHOTGUN 1 // stock defines
#define IT_SUPER_SHOTGUN 2
#define IT_NAILGUN 4
#define IT_SUPER_NAILGUN 8
#define IT_GRENADE_LAUNCHER 16
#define IT_ROCKET_LAUNCHER 32
#define IT_LIGHTNING 64
#define IT_SUPER_LIGHTNING 128
#define IT_SHELLS 256
#define IT_NAILS 512
#define IT_ROCKETS 1024
#define IT_CELLS 2048
#define IT_AXE 4096
#define IT_ARMOR1 8192
#define IT_ARMOR2 16384
#define IT_ARMOR3 32768
#define IT_SUPERHEALTH 65536
#define IT_KEY1 131072
#define IT_KEY2 262144
#define IT_INVISIBILITY 524288
#define IT_INVULNERABILITY 1048576
#define IT_SUIT 2097152
#define IT_QUAD 4194304
#define IT_SIGIL1 (1<<28)
#define IT_SIGIL2 (1<<29)
#define IT_SIGIL3 (1<<30)
#define IT_SIGIL4 (1<<31)
#define RIT_SHELLS 128 // rogue changed and added defines
#define RIT_NAILS 256
#define RIT_ROCKETS 512
#define RIT_CELLS 1024
#define RIT_AXE 2048
#define RIT_LAVA_NAILGUN 4096
#define RIT_LAVA_SUPER_NAILGUN 8192
#define RIT_MULTI_GRENADE 16384
#define RIT_MULTI_ROCKET 32768
#define RIT_PLASMA_GUN 65536
#define RIT_ARMOR1 8388608
#define RIT_ARMOR2 16777216
#define RIT_ARMOR3 33554432
#define RIT_LAVA_NAILS 67108864
#define RIT_PLASMA_AMMO 134217728
#define RIT_MULTI_ROCKETS 268435456
#define RIT_SHIELD 536870912
#define RIT_ANTIGRAV 1073741824
#define RIT_SUPERHEALTH 2147483648
#define HIT_PROXIMITY_GUN_BIT 16 // MED 01/04/97 added hipnotic defines
#define HIT_MJOLNIR_BIT 7 // hipnotic added defines
#define HIT_LASER_CANNON_BIT 23
#define HIT_PROXIMITY_GUN (1<<HIT_PROXIMITY_GUN_BIT)
#define HIT_MJOLNIR (1<<HIT_MJOLNIR_BIT)
#define HIT_LASER_CANNON (1<<HIT_LASER_CANNON_BIT)
#define HIT_WETSUIT (1<<(23+2))
#define HIT_EMPATHY_SHIELDS (1<<(23+3))
#define MAX_SCOREBOARD 16
#define MAX_SCOREBOARDNAME 32
#define SOUND_CHANNELS 8
#define DATAGRAM_MTU 1400 // johnfitz -- actual limit for unreliable messages to nonlocal clients
#define DIST_EPSILON (0.03125) // 1/32 epsilon to keep floating point happy (moved from world.c)
// This makes anyone on id's net privileged
// Use for multiplayer testing only - VERY dangerous!!!
// #define IDGODS

#include "mathlib.h"

// the host system specifies the base of the directory tree, the
// command line parms passed to the program, and the amount of memory
// available for the program to use
typedef struct {
	char *basedir;
	char *userdir;
	char *cachedir; // for development over ISDN lines
	int argc;
	char **argv;
	void *membase;
	int memsize;
} quakeparms_t;

#include "q_stdinc.h"
#include "common.h"
#include "bspfile.h"
#include "vid.h"
#include "sys.h"
#include "zone.h"
#include "wad.h"
#include "draw.h"
#include "cvar.h"
#include "screen.h"
#include "net.h"
#include "protocol.h"
#include "cmd.h"
#include "sbar.h"
#include "q_sound.h"
#include "render.h"
#include "client.h"
#include "progs.h"
#include "server.h"
#include "model.h"
#include "d_iface.h"
#include "input.h"
#include "world.h"
#include "keys.h"
#include "console.h"
#include "view.h"
#include "menu.h"
#include "crc.h"
#include "vid.h"
#include "vgatext.h"

extern double host_time;
extern SDL_Window *window; // global for checking windowed state in options
extern Uint32 SDLWindowFlags;
extern qboolean noclip_anglehack;
extern quakeparms_t host_parms;
extern cvar_t sys_ticrate;
extern cvar_t sys_nostdout;
extern cvar_t developer;
extern qboolean host_initialized; // true if into command execution
extern double host_frametime;
extern byte *host_basepal;
extern byte *host_colormap;
extern qboolean isDedicated;
extern int minimum_memory;
extern cvar_t chase_active;
extern int host_framecount; // incremented every frame, never reset
extern double realtime; // not bounded in any way, changed at
			// start of every frame, never reset
extern qboolean msg_suppress_1; // suppresses resolution and cache size console
				// output an fullscreen DIB focus gain/loss
extern int current_skill; // skill level for currently loaded level (in case
			  // the user changes the cvar while the level is
			  // running, this reflects the level actually in use)
extern byte r_foundtranswater, r_wateralphapass; // Manoel Kasimier - translucent water
extern cvar_t max_edicts;

void Host_ClearMemory();
void Host_ServerFrame();
void Host_InitCommands();
void Host_Init();
void Host_Shutdown();
void Host_Error(char *error, ...);
void Host_EndGame(char *message, ...);
void Host_Frame(float time);
void Host_Quit_f();
void Host_ClientCommands(char *fmt, ...);
void Host_ShutdownServer(qboolean crash);
void Chase_Init();
void Chase_Update();
void Cvar_SetCallback(cvar_t *var, cvarcallback_t func);
#endif
