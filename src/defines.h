#ifndef GLOBDEFINES_
#define GLOBDEFINES_

#define q_min(a, b) (((a) < (b)) ? (a) : (b)) // common.h
#define q_max(a, b) (((a) > (b)) ? (a) : (b))
#define CLAMP(_minval, x, _maxval) \
 ((x) < (_minval) ? (_minval) : \
 (x) > (_maxval) ? (_maxval) : (x))
#define STRUCT_FROM_LINK(l,t,m) ((t *)((u8 *)l - offsetof(t,m)))
#define FS_ENT_NONE (0)
#define FS_ENT_FILE (1 << 0)
#define FS_ENT_DIRECTORY (1 << 1)
#define CURWORLDPAL ((u8*)(worldpalname[0]?worldpal:host_basepal))
#define CURWORLDCMAP ((u8*)(worldcmapname[0]?worldcmap:host_colormap))
#define Q_rint(x) ((x) > 0 ? (s32)((x) + 0.5) : (s32)((x) - 0.5)) // johnfitz -- from joequake
#define Q_COUNTOF(x) (sizeof(x)/sizeof((x)[0])) // for array size
#define NUM_TYPE_SIZES 8 // progs.h
#define MAX_ENT_LEAFS 32
#define EDICT_FROM_AREA(l) STRUCT_FROM_LINK(l,edict_t,area)
#define NEXT_EDICT(e) ((edict_t*)((u8*)e + pr_edict_size))
#define EDICT_TO_PROG(e) (s32)((u8*)e - (u8*)sv.edicts)
#define PROG_TO_EDICT(e) ((edict_t*)((u8*)sv.edicts + e))
#define G_FLOAT(o) (pr_globals[o])
#define G_INT(o) (*(s32*)&pr_globals[o])
#define G_EDICT(o) ((edict_t*)((u8*)sv.edicts+ *(s32*)&pr_globals[o]))
#define G_EDICTNUM(o) NUM_FOR_EDICT(G_EDICT(o))
#define G_VECTOR(o) (&pr_globals[o])
#define G_STRING(o) (PR_GetString(*(string_t*)&pr_globals[o]))
#define G_FUNCTION(o) (*(func_t*)&pr_globals[o])
#define E_FLOAT(e,o) (((f32*)&e->v)[o])
#define E_INT(e,o) (*(s32*)&((f32*)&e->v)[o])
#define E_VECTOR(e,o) (&((f32*)&e->v)[o])
#define E_STRING(e,o) (PR_GetString(*(string_t*)&((f32*)&e->v)[o]))
#define PROGHEADER_CRC 5927
#define WAV_FORMAT_PCM 1 // q_sound.h
#define MAX_CHANNELS 1024
#define MAX_DYNAMIC_CHANNELS 128
#define MAX_RAW_SAMPLES 8192
#define NET_HEADERSIZE (2 * sizeof(u32)) // net_defs.h
#define NET_DATAGRAMSIZE (MAX_DATAGRAM + NET_HEADERSIZE)
#define NETFLAG_LENGTH_MASK 0x0000ffff // Must fit NET_MAXMESSAGE
#define NETFLAG_DATA 0x00010000
#define NETFLAG_ACK 0x00020000
#define NETFLAG_NAK 0x00040000
#define NETFLAG_EOM 0x00080000
#define NETFLAG_UNRELIABLE 0x00100000
#define NETFLAG_CTL 0x80000000
#define NET_PROTOCOL_VERSION 3
#define CCREQ_CONNECT 0x01
#define CCREQ_SERVER_INFO 0x02
#define CCREQ_PLAYER_INFO 0x03
#define CCREQ_RULE_INFO 0x04
#define CCREP_ACCEPT 0x81
#define CCREP_REJECT 0x82
#define CCREP_SERVER_INFO 0x83
#define CCREP_PLAYER_INFO 0x84
#define CCREP_RULE_INFO 0x85
#define MAX_NET_DRIVERS 8
#define HOSTCACHESIZE 8
#define IS_LOOP_DRIVER(p) ((p) == 0) // Loop driver must always be registered the first
#define PROTOCOL_NETQUAKE 15 //johnfitz -- standard quake protocol // protocol.h
#define PROTOCOL_FITZQUAKE 666 //johnfitz -- added new protocol for fitzquake 0.85
#define PROTOCOL_RMQ 999
#define PRFL_SHORTANGLE (1 << 1)
#define PRFL_FLOATANGLE (1 << 2)
#define PRFL_24BITCOORD (1 << 3)
#define PRFL_FLOATCOORD (1 << 4)
#define PRFL_EDICTSCALE (1 << 5)
#define PRFL_ALPHASANITY (1 << 6) // cleanup insanity with alpha
#define PRFL_INT32COORD (1 << 7)
#define PRFL_MOREFLAGS (1 << 31) // not supported
#define U_MOREBITS (1<<0)
#define U_ORIGIN1 (1<<1)
#define U_ORIGIN2 (1<<2)
#define U_ORIGIN3 (1<<3)
#define U_ANGLE2 (1<<4)
#define U_STEP (1<<5) //johnfitz -- was U_NOLERP, renamed since it's only used for MOVETYPE_STEP
#define U_FRAME (1<<6)
#define U_SIGNAL (1<<7) // just differentiates from other updates
#define U_ANGLE1 (1<<8) // svc_update can pass all of the fast update bits, plus more
#define U_ANGLE3 (1<<9)
#define U_MODEL (1<<10)
#define U_COLORMAP (1<<11)
#define U_SKIN (1<<12)
#define U_EFFECTS (1<<13)
#define U_LONGENTITY (1<<14)
#define U_EXTEND1 (1<<15)
#define U_ALPHA (1<<16) // 1 u8, uses ENTALPHA_ENCODE, not sent if equal to baseline
#define U_FRAME2 (1<<17) // 1 u8, this is .frame & 0xFF00 (second u8)
#define U_MODEL2 (1<<18) // 1 u8, this is .modelindex & 0xFF00 (second u8)
#define U_LERPFINISH (1<<19) // 1 u8, 0.0-1.0 maps to 0-255, not sent if exactly 0.1, this is ent->v.nextthink - sv.time, used for lerping
#define U_SCALE (1<<20) // 1 u8, for PROTOCOL_RMQ PRFL_EDICTSCALE
#define U_UNUSED21 (1<<21)
#define U_UNUSED22 (1<<22)
#define U_EXTEND2 (1<<23) // another u8 to follow, future expansion
#define U_TRANS (1<<15)
#define SU_VIEWHEIGHT (1<<0)
#define SU_IDEALPITCH (1<<1)
#define SU_PUNCH1 (1<<2)
#define SU_PUNCH2 (1<<3)
#define SU_PUNCH3 (1<<4)
#define SU_VELOCITY1 (1<<5)
#define SU_VELOCITY2 (1<<6)
#define SU_VELOCITY3 (1<<7)
#define SU_UNUSED8 (1<<8) // AVAILABLE BIT
#define SU_ITEMS (1<<9)
#define SU_ONGROUND (1<<10) // no data follows, the bit is it
#define SU_INWATER (1<<11) // no data follows, the bit is it
#define SU_WEAPONFRAME (1<<12)
#define SU_ARMOR (1<<13)
#define SU_WEAPON (1<<14)
#define SU_EXTEND1 (1<<15) // another u8 to follow
#define SU_WEAPON2 (1<<16) // 1 u8, this is .weaponmodel & 0xFF00 (second u8)
#define SU_ARMOR2 (1<<17) // 1 u8, this is .armorvalue & 0xFF00 (second u8)
#define SU_AMMO2 (1<<18) // 1 u8, this is .currentammo & 0xFF00 (second u8)
#define SU_SHELLS2 (1<<19) // 1 u8, this is .ammo_shells & 0xFF00 (second u8)
#define SU_NAILS2 (1<<20) // 1 u8, this is .ammo_nails & 0xFF00 (second u8)
#define SU_ROCKETS2 (1<<21) // 1 u8, this is .ammo_rockets & 0xFF00 (second u8)
#define SU_CELLS2 (1<<22) // 1 u8, this is .ammo_cells & 0xFF00 (second u8)
#define SU_EXTEND2 (1<<23) // another u8 to follow
#define SU_WEAPONFRAME2 (1<<24) // 1 u8, this is .weaponframe & 0xFF00 (second u8)
#define SU_WEAPONALPHA (1<<25) // 1 u8, this is alpha for weaponmodel, uses ENTALPHA_ENCODE, not sent if ENTALPHA_DEFAULT
#define SU_UNUSED26 (1<<26)
#define SU_UNUSED27 (1<<27)
#define SU_UNUSED28 (1<<28)
#define SU_UNUSED29 (1<<29)
#define SU_UNUSED30 (1<<30)
#define SU_EXTEND3 (1<<31) // another u8 to follow, future expansion
#define SND_VOLUME (1<<0) // a u8 // a sound with no channel is a local only sound
#define SND_ATTENUATION (1<<1) // a u8
#define SND_LOOPING (1<<2) // a s64
#define DEFAULT_SOUND_PACKET_VOLUME 255
#define DEFAULT_SOUND_PACKET_ATTENUATION 1.0
#define SND_LARGEENTITY (1<<3) // a s16 + u8 (instead of just a s16)
#define SND_LARGESOUND (1<<4) // a s16 soundindex (instead of a u8)
#define B_LARGEMODEL (1<<0) // modelindex is s16 instead of u8
#define B_LARGEFRAME (1<<1) // frame is s16 instead of u8
#define B_ALPHA (1<<2) // 1 u8, uses ENTALPHA_ENCODE, not sent if ENTALPHA_DEFAULT
#define B_SCALE (1<<3)
#define ENTALPHA_DEFAULT 0 //entity's alpha is "default" (i.e. water obeys r_wateralpha) -- must be zero so zeroed out memory works
#define ENTALPHA_ZERO 1 //entity is invisible (lowest possible alpha)
#define ENTALPHA_ONE 255 //entity is fully opaque (highest possible alpha)
#define ENTSCALE_DEFAULT 16 // Equivalent to f32 1.0f due to u8 packing.
#define DEFAULT_VIEWHEIGHT 22 // defaults for clientinfo messages
#define GAME_COOP 0 // game types sent by serverinfo
#define GAME_DEATHMATCH 1 // these determine which intermission screen plays
// note that there are some defs.qc that mirror to these numbers
// also related to svc_strings[] in cl_parse
#define svc_bad 0 // server to client
#define svc_nop 1
#define svc_disconnect 2
#define svc_updatestat 3 // [u8] [s64]
#define svc_version 4 // [s64] server version
#define svc_setview 5 // [s16] entity number
#define svc_sound 6 // <see code>
#define svc_time 7 // [f32] server time
#define svc_print 8 // [string] null terminated string
#define svc_stufftext 9 // [string] stuffed into client's console buffer, \n terminated
#define svc_setangle 10 // [angle3] set the view angle to this absolute value
#define svc_serverinfo 11 // [s64] version
 // [string] signon string
 // [string]..[0]model cache
 // [string]...[0]sounds cache
#define svc_lightstyle 12 // [u8] [string]
#define svc_updatename 13 // [u8] [string]
#define svc_updatefrags 14 // [u8] [s16]
#define svc_clientdata 15 // <shortbits + data>
#define svc_stopsound 16 // <see code>
#define svc_updatecolors 17 // [u8] [u8]
#define svc_particle 18 // [vec3] <variable>
#define svc_damage 19
#define svc_spawnstatic 20
#define svc_spawnbinary 21
#define svc_spawnbaseline 22
#define svc_temp_entity 23
#define svc_setpause 24 // [u8] on / off
#define svc_signonnum 25 // [u8] used for the signon sequence
#define svc_centerprint 26 // [string] to put in center of the screen
#define svc_killedmonster 27
#define svc_foundsecret 28
#define svc_spawnstaticsound 29 // [coord3] [u8] samp [u8] vol [u8] aten
#define svc_intermission 30 // [string] music
#define svc_finale 31 // [string] music [string] text
#define svc_cdtrack 32 // [u8] track [u8] looptrack
#define svc_sellscreen 33
#define svc_cutscene 34
#define svc_skybox 37 // [string] name // PROTOCOL_FITZQUAKE -- new server messages
#define svc_bf 40
#define svc_fog 41 // [u8] density [u8] red [u8] green [u8] blue [f32] time
#define svc_spawnbaseline2 42 // support for large modelindex, large framenum, alpha, using flags
#define svc_spawnstatic2 43 // support for large modelindex, large framenum, alpha, using flags
#define svc_spawnstaticsound2 44 // [coord3] [s16] samp [u8] vol [u8] aten
#define svc_botchat 38 // 2021 re-release server messages - see:
#define svc_setviews 45 // https://steamcommunity.com/sharedfiles/filedetails/?id=2679459726
#define svc_updateping 46
#define svc_updatesocial 47
#define svc_updateplinfo 48
#define svc_rawprint 49
#define svc_servervars 50
#define svc_seq 51
#define svc_achievement 52 // [string] id // Note: svc_achievement has same value as svcdp_effect!
#define svc_chat 53
#define svc_levelcompleted 54
#define svc_backtolobby 55
#define svc_localsound 56
#define clc_bad 0 // client to server
#define clc_nop 1
#define clc_disconnect 2
#define clc_move 3 // [usercmd_t]
#define clc_stringcmd 4 // [string] message
#define TE_SPIKE 0 // temp entity events
#define TE_SUPERSPIKE 1
#define TE_GUNSHOT 2
#define TE_EXPLOSION 3
#define TE_TAREXPLOSION 4
#define TE_LIGHTNING1 5
#define TE_LIGHTNING2 6
#define TE_WIZSPIKE 7
#define TE_KNIGHTSPIKE 8
#define TE_LIGHTNING3 9
#define TE_LAVASPLASH 10
#define TE_TELEPORT 11
#define TE_EXPLOSION2 12
#define TE_BEAM 13
#define ENTALPHA_ENCODE(a) (((a)==0)?ENTALPHA_DEFAULT:Q_rint(CLAMP(1.0f,(a)*254.0f+1,255.0f))) // server convert to u8 to send to client
#define ENTALPHA_DECODE(a) (((a)==ENTALPHA_DEFAULT)?1.0f:((f32)(a)-1)/(254)) // client convert to f32 for rendering
#define ENTALPHA_TOSAVE(a) (((a)==ENTALPHA_DEFAULT)?0.0f:(((a)==ENTALPHA_ZERO)?-1.0f:((f32)(a)-1)/(254))) // server convert to f32 for savegame
#define ENTSCALE_ENCODE(a) ((a) ? ((a) * ENTSCALE_DEFAULT) : ENTSCALE_DEFAULT) // Convert to u8
#define ENTSCALE_DECODE(a) ((f32)(a) / ENTSCALE_DEFAULT) // Convert to f32 for rendering
#define CVAR_NONE 0 // cvar.h
#define CVAR_ARCHIVE (1U << 0) // if set, causes it to be saved to config
#define CVAR_NOTIFY (1U << 1) // changes will be broadcasted to all players (q1)
#define CVAR_SERVERINFO (1U << 2) // added to serverinfo will be sent to clients (q1/net_dgrm.c and qwsv)
#define CVAR_USERINFO (1U << 3) // added to userinfo, will be sent to server (qwcl)
#define CVAR_CHANGED (1U << 4)
#define CVAR_ROM (1U << 6)
#define CVAR_LOCKED (1U << 8) // locked temporarily
#define CVAR_REGISTERED (1U << 10) // the var is added to the list of variables
#define CVAR_CALLBACK (1U << 16) // var has a callback
#define NET_NAMELEN 64 // net.h
#define NET_MAXMESSAGE 65535 // ericw -- was 32000
#define R_SKY_SMASK 0x007F0000 // d_local.h
#define R_SKY_TMASK 0x007F0000
#define DS_SPAN_LIST_END -128
#define SKYBOX_MAX_SIZE 1024
#define SURFCACHE_SIZE_AT_320X200 600*1024
#define FOG_LUT_LEVELS 32
#define MAX_ALIAS_NAME 32 // cmd.h
#define MAX_ARGS 80
#define MOVE_NORMAL 0 // world.h
#define MOVE_NOMONSTERS 1
#define MOVE_MISSILE 2
#define MAX_HANDLES 32 // sys.h
#define MAX_SIGNON_BUFFERS 256 // server.h
#define NUM_PING_TIMES 16
#define NUM_SPAWN_PARMS 16
#define MOVETYPE_NONE 0 // never moves // edict->movetype values
#define MOVETYPE_ANGLENOCLIP 1
#define MOVETYPE_ANGLECLIP 2
#define MOVETYPE_WALK 3 // gravity
#define MOVETYPE_STEP 4 // gravity, special edge handling
#define MOVETYPE_FLY 5
#define MOVETYPE_TOSS 6 // gravity
#define MOVETYPE_PUSH 7 // no clip to world, push and crush
#define MOVETYPE_NOCLIP 8
#define MOVETYPE_FLYMISSILE 9 // extra size to monsters
#define MOVETYPE_BOUNCE 10
#define MOVETYPE_GIB 11 // 2021 rerelease gibs
#define SOLID_NOT 0 // no interaction with other objects // edict->solid values
#define SOLID_TRIGGER 1 // touch on edge, but not blocking
#define SOLID_BBOX 2 // touch on edge, block
#define SOLID_SLIDEBOX 3 // touch on edge, but not an onground
#define SOLID_BSP 4 // bsp clip, touch on edge, block
#define DEAD_NO 0 // edict->deadflag values
#define DEAD_DYING 1
#define DEAD_DEAD 2
#define DAMAGE_NO 0
#define DAMAGE_YES 1
#define DAMAGE_AIM 2
#define FL_FLY 1 // edict->flags
#define FL_SWIM 2
#define FL_CONVEYOR 4
#define FL_CLIENT 8
#define FL_INWATER 16
#define FL_MONSTER 32
#define FL_GODMODE 64
#define FL_NOTARGET 128
#define FL_ITEM 256
#define FL_ONGROUND 512
#define FL_PARTIALGROUND 1024 // not all corners are valid
#define FL_WATERJUMP 2048 // player jumping out of water
#define FL_JUMPRELEASED 4096 // for jump debouncing
#define EF_BRIGHTFIELD 1 // entity effects
#define EF_MUZZLEFLASH 2
#define EF_BRIGHTLIGHT 4
#define EF_DIMLIGHT 8
#define SPAWNFLAG_NOT_EASY 256
#define SPAWNFLAG_NOT_MEDIUM 512
#define SPAWNFLAG_NOT_HARD 1024
#define SPAWNFLAG_NOT_DEATHMATCH 2048
#define OFS_NULL 0 // pr_comp.h
#define OFS_RETURN 1
#define OFS_PARM0 4 // leave 3 ofs for each parm to hold vectors
#define OFS_PARM1 7
#define OFS_PARM2 10
#define OFS_PARM3 13
#define OFS_PARM4 16
#define OFS_PARM5 19
#define OFS_PARM6 22
#define OFS_PARM7 25
#define RESERVED_OFS 28
#define DEF_SAVEGLOBAL (1<<15)
#define MAX_PARMS 8
#define PROG_VERSION 6
#define NUM_OLDMODES 9 // vid_sdl.h
#define VID_GRADES (1 << VID_CBITS)
#define VID_CBITS 6
#define MAX_MODE_LIST 30
#define VID_ROW_SIZE 3
#define MAX_COLUMN_SIZE 5
#define MODE_AREA_HEIGHT (MAX_COLUMN_SIZE + 6)
#define RANDARR_SIZE 65536 // d_fog.h
#define FOG_LUT_LEVELS 32 // rgbtoi.h
#define LIT_LUT_RES 64
#define RGB_LUT_SIZE 512 // 8×8×8 color space
#define CSHIFT_CONTENTS 0 // client.h
#define CSHIFT_DAMAGE 1
#define CSHIFT_BONUS 2
#define CSHIFT_POWERUP 3
#define NUM_CSHIFTS 4
#define NAME_LENGTH 64
#define SIGNONS 4 // signon messages to receive before connected
#define MAX_DLIGHTS 64
#define MAX_BEAMS 32
#define MAX_EFRAGS 2048
#define MAX_MAPSTRING 2048
#define MAX_DEMOS 8
#define MAX_DEMONAME 16
#define MAX_TEMP_ENTITIES 1024 // lightning bolts, etc
#define MAX_STATIC_ENTITIES 2048 // torches, etc
#define MAX_VISEDICTS 16384 // larger, now we support BSP2
#define K_TAB 9 // keys.h
#define K_ENTER 13 // these are key numbers that should be passed to Key_Event
#define K_ESCAPE 27 // normal keys should be passed as lowercase ascii
#define K_SPACE 32
#define K_BACKSPACE 127
#define K_UPARROW 128
#define K_DOWNARROW 129
#define K_LEFTARROW 130
#define K_RIGHTARROW 131
#define K_ALT 132
#define K_CTRL 133
#define K_SHIFT 134
#define K_F1 135
#define K_F2 136
#define K_F3 137
#define K_F4 138
#define K_F5 139
#define K_F6 140
#define K_F7 141
#define K_F8 142
#define K_F9 143
#define K_F10 144
#define K_F11 145
#define K_F12 146
#define K_INS 147
#define K_DEL 148
#define K_PGDN 149
#define K_PGUP 150
#define K_HOME 151
#define K_END 152
#define K_PAUSE 255
#define K_MOUSE1 200 // mouse buttons generate virtual keys
#define K_MOUSE2 201
#define K_MOUSE3 202
#define K_JOY1 203 // joystick buttons
#define K_JOY2 204
#define K_JOY3 205
#define K_JOY4 206
#define K_AUX1 207 // aux keys are for multi-buttoned joysticks to generate
#define K_AUX2 208 // so they can use the normal binding process
#define K_AUX3 209
#define K_AUX4 210
#define K_AUX5 211
#define K_AUX6 212
#define K_AUX7 213
#define K_AUX8 214
#define K_AUX9 215
#define K_AUX10 216
#define K_AUX11 217
#define K_AUX12 218
#define K_AUX13 219
#define K_AUX14 220
#define K_AUX15 221
#define K_AUX16 222
#define K_AUX17 223
#define K_AUX18 224
#define K_AUX19 225
#define K_AUX20 226
#define K_AUX21 227
#define K_AUX22 228
#define K_AUX23 229
#define K_AUX24 230
#define K_AUX25 231
#define K_AUX26 232
#define K_AUX27 233
#define K_AUX28 234
#define K_AUX29 235
#define K_AUX30 236
#define K_AUX31 237
#define K_AUX32 238
#define K_MWHEELUP 239
#define K_MWHEELDOWN 240
#define K_MOUSE4 241
#define K_MOUSE5 242
#define MAXCLIPPLANES 11 // render.h
#define TOP_RANGE 16 // soldier uniform colors
#define BOTTOM_RANGE 96
#define LERP_FINISH (1<<4) //use lerpfinish time from server update instead of assuming interval of 0.1
#undef HAVE_SA_LEN // net_sys.h
#define SA_FAM_OFFSET 0
#ifdef _WIN32
#define MAXHOSTNAMELEN 1024
#define selectsocket select
#define SOCKETERRNO WSAGetLastError()
#define NET_EWOULDBLOCK WSAEWOULDBLOCK
#define NET_ECONNREFUSED WSAECONNREFUSED
#define socketerror(x) __WSAE_StrError((x)) // must include wsaerror.h for this
#else
#define INVALID_SOCKET (-1)
#define SOCKET_ERROR (-1)
#define SOCKETERRNO errno
#define ioctlsocket ioctl
#define closesocket close
#define selectsocket select
#define NET_EWOULDBLOCK EWOULDBLOCK
#define NET_ECONNREFUSED ECONNREFUSED
#define socketerror(x) strerror((x))
#endif
#define ALIAS_BASE_SIZE_RATIO (1.0 / 11.0) // r_local.h // normalizing factor so player model works out to about 1 pixel per triangle
#define BMODEL_FULLY_CLIPPED 0x10 // value returned by R_BmodelCheckBBox () if bbox is trivially rejected
#define CLIP_EPSILON 0.001
#define BACKFACE_EPSILON 0.01
#define DIST_NOT_SET 98765
#define NEAR_CLIP 0.01
#define ALIAS_Z_CLIP_PLANE 0.1 // CyanBun96: was 5, lowered for high FOV guns
#define MAX_BTOFPOLYS 5000 // FIXME: tune this
#define AMP 8*0x10000 // turbulence stuff
#define AMP2 3
#define SPEED 20
#define LIT_LUT_RES 64
#define QUANT(x) (((x) * (LIT_LUT_RES - 1)) / 255)
#define CMP_NONE 0 // wad.h
#define CMP_LZSS 1
#define TYP_NONE 0
#define TYP_LABEL 1
#define TYP_LUMPY 64 // 64 + grab command number
#define TYP_PALETTE 64
#define TYP_QTEX 65
#define TYP_QPIC 66
#define TYP_SOUND 67
#define TYP_MIPTEX_PALETTE 67
#define TYP_MIPTEX 68
#define WADID ('W' | ('A' << 8) | ('D' << 16) | ('2' << 24))
#define WADID_VALVE ('W' | ('A' << 8) | ('D' << 16) | ('3' << 24))
#define WADFILENAME "gfx.wad" //johnfitz -- filename is now hard-coded for honesty
#define VERSION 0.62 // quakedef.h
#define FITZQUAKE_VERSION 0.85
#define GAMENAME "id1"
#define CMDLINE_LENGTH 256
#define DEFAULT_MEMORY (2 * 1000 * 1000 * 1000)
#define CACHE_SIZE 32 // used to align key data structures
#define MAX_NUM_ARGVS 50
#define PITCH 0 // up / down
#define YAW 1 // left / right
#define ROLL 2 // fall over
#define MAX_QPATH 64 // max length of a quake game pathname
#define MAX_OSPATH 128 // max length of a filesystem pathname
#define ON_EPSILON 0.1 // point on plane side epsilon
#define MAX_MSGLEN 65536 // max length of a reliable message, from Mark V
#define MAX_DATAGRAM 65527 // max length of unreliable message, from Mark V
#define MAX_EDICTS 32000
#define MIN_EDICTS 256 // johnfitz -- lowest allowed value for max_edicts cvar
#define MAX_LIGHTSTYLES 64
#define MAX_MODELS 4096 // these are sent over the net as bytes
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
#define WARP_WIDTH MAXWIDTH // d_iface.h
#define WARP_HEIGHT MAXHEIGHT
#define MAX_LBM_HEIGHT 480
#define PARTICLE_Z_CLIP 8.0
#define DR_SOLID 0 // transparency types for D_DrawRect()
#define DR_TRANSPARENT 1
#define TRANSPARENT_COLOR 0xFF
#define TILE_SIZE 128 // size of textures generated by R_GenTiledSurf
#define SKYSHIFT 7
#define SKYSIZE (1 << SKYSHIFT)
#define SKYMASK (SKYSIZE - 1)
#define CYCLE 128 // turbulent cycle size
#define MAXVERTS 16 // max points in a surface polygon // r_shared.h
#define MAXWORKINGVERTS (MAXVERTS+4) // max points in an intermediate polygon (while processing)
#ifdef __EMSCRIPTEN__
#define MAXHEIGHT 600 // erysdren
#define MAXWIDTH 800
#else
#define MAXHEIGHT 8640 // CyanBun96: 16k resolution. futureproofing.
#define MAXWIDTH 15360
#endif
#define MAXDIMENSION ((MAXHEIGHT > MAXWIDTH) ? MAXHEIGHT : MAXWIDTH)
#define SIN_BUFFER_SIZE (MAXDIMENSION+CYCLE)
#define INFINITE_DISTANCE 0x10000 // distance that's always guaranteed to be farther away than anything in the scene
#define NUMSTACKEDGES 262144 // CyanBun96: expanding limits
#define MINEDGES NUMSTACKEDGES
#define NUMSTACKSURFACES 32768 // CyanBun96: expanding limits
#define MINSURFACES NUMSTACKSURFACES
#define MAXSPANS 32768 // CyanBun96: expanding limits
#define ALIAS_LEFT_CLIP 0x0001 // flags in finalvert_t.flags
#define ALIAS_TOP_CLIP 0x0002
#define ALIAS_RIGHT_CLIP 0x0004
#define ALIAS_BOTTOM_CLIP 0x0008
#define ALIAS_Z_CLIP 0x0010
#define ALIAS_XY_CLIP_MASK 0x000F
#define TEX_SPECIAL 1 // sky or slime, no lightmap or 256 subd // bspfile.h
#define TEX_MISSING 2 // johnfitz -- this texinfo does not have a texture
#define MAX_MAP_HULLS 4 // upper design bounds
#define MAX_KEY 32 // key / value pair sizes
#define MAX_VALUE 1024
#define BSPVERSION 29
#ifdef BSP29_VALVE
#define BSPVERSION_VALVE 30
#endif
// RMQ support (2PSB). 32bits instead of shorts for all but bbox sizes (which still use shorts)
#define BSP2VERSION_2PSB (('B' << 24) | ('S' << 16) | ('P' << 8) | '2')
// BSP2 support. 32bits instead of shorts for everything (bboxes use floats)
#define BSP2VERSION_BSP2 (('B' << 0) | ('S' << 8) | ('P' << 16) | ('2'<<24))
#define BSPVERSION_QUAKE64 (('Q' << 24) | ('6' << 16) | ('4' << 8) | ' ')
#define TOOLVERSION 2
#define LUMP_ENTITIES 0
#define LUMP_PLANES 1
#define LUMP_TEXTURES 2
#define LUMP_VERTEXES 3
#define LUMP_VISIBILITY 4
#define LUMP_NODES 5
#define LUMP_TEXINFO 6
#define LUMP_FACES 7
#define LUMP_LIGHTING 8
#define LUMP_CLIPNODES 9
#define LUMP_LEAFS 10
#define LUMP_MARKSURFACES 11
#define LUMP_EDGES 12
#define LUMP_SURFEDGES 13
#define LUMP_MODELS 14
#define HEADER_LUMPS 15
#define MIPLEVELS 4
#define PLANE_X 0 // 0-2 are axial planes
#define PLANE_Y 1
#define PLANE_Z 2
#define PLANE_ANYX 3 // 3-5 are non-axial planes snapped to the nearest
#define PLANE_ANYY 4
#define PLANE_ANYZ 5
#define CONTENTS_EMPTY -1
#define CONTENTS_SOLID -2
#define CONTENTS_WATER -3
#define CONTENTS_SLIME -4
#define CONTENTS_LAVA -5
#define CONTENTS_SKY -6
#define CONTENTS_ORIGIN -7 // removed at csg time
#define CONTENTS_CLIP -8 // changed to contents_solid
#define CONTENTS_CURRENT_0 -9
#define CONTENTS_CURRENT_90 -10
#define CONTENTS_CURRENT_180 -11
#define CONTENTS_CURRENT_270 -12
#define CONTENTS_CURRENT_UP -13
#define CONTENTS_CURRENT_DOWN -14
#define TEX_SPECIAL 1 // sky or slime, no lightmap or 256 subdivision
#define MAXLIGHTMAPS 4
#define AMBIENT_WATER 0
#define AMBIENT_SKY 1
#define AMBIENT_SLIME 2
#define AMBIENT_LAVA 3
#define NUM_AMBIENTS 4 // automatic ambient sounds
#define SPRITE_VERSION 1 // spritegn.h
#define SPR_VP_PARALLEL_UPRIGHT 0
#define SPR_FACING_UPRIGHT 1
#define SPR_VP_PARALLEL 2
#define SPR_ORIENTED 3
#define SPR_VP_PARALLEL_ORIENTED 4
#define IDSPRITEHEADER (('P'<<24)+('S'<<16)+('D'<<8)+'I') // little-endian IDSP
#define ALIAS_VERSION 6 // modelgen.h
#define ALIAS_ONSEAM 0x0020
#define DT_FACES_FRONT 0x0010
#define IDPOLYHEADER (('O'<<24)+('P'<<16)+('D'<<8)+'I') // little-endian "IDPO"
#define SIDE_FRONT 0 // model.h
#define SIDE_BACK 1
#define SIDE_ON 2
#define SURF_PLANEBACK 2
#define SURF_DRAWSKY 4
#define SURF_DRAWSPRITE 8
#define SURF_DRAWTURB 0x10
#define SURF_DRAWTILED 0x20
#define SURF_DRAWBACKGROUND 0x40
#define SURF_DRAWCUTOUT 0x80
#define SURF_DRAWLAVA 0x400
#define SURF_DRAWSLIME 0x800
#define SURF_DRAWTELE 0x1000
#define SURF_DRAWWATER 0x2000
#define SURF_DRAWSKYBOX 0x80000 // Manoel Kasimier - skyboxes
#define SURF_WINQUAKE_DRAWTRANSLUCENT (SURF_DRAWLAVA | SURF_DRAWSLIME | SURF_DRAWWATER | SURF_DRAWTELE)
#define EF_ROCKET 1 // leave a trail
#define EF_GRENADE 2 // leave a trail
#define EF_GIB 4 // leave a trail
#define EF_ROTATE 8 // rotate(bonus items)
#define EF_TRACER 16 // green split trail
#define EF_ZOMGIB 32 // small blood trail
#define EF_TRACER2 64 // orange split trail + rotate
#define EF_TRACER3 128 // purple trail
#define EF_BRIGHTFIELD 1 // entity effects
#define EF_MUZZLEFLASH 2
#define EF_BRIGHTLIGHT 4
#define EF_DIMLIGHT 8
#define EF_QEX_QUADLIGHT 16 // 2021 rerelease
#define EF_QEX_PENTALIGHT 32 // 2021 rerelease
#define EF_QEX_CANDLELIGHT 64 // 2021 rerelease
#define SIDE_FRONT 0
#define SIDE_BACK 1
#define SIDE_ON 2
#define SURF_PLANEBACK 2
#define SURF_DRAWSKY 4
#define SURF_DRAWSPRITE 8
#define SURF_DRAWTURB 0x10
#define SURF_DRAWTILED 0x20
#define SURF_DRAWBACKGROUND 0x40
#define SURF_UNDERWATER 0x80
#define SURF_NOTEXTURE 0x100 //johnfitz
#define SURF_DRAWFENCE 0x200
#define SURF_DRAWLAVA 0x400
#define SURF_DRAWSLIME 0x800
#define SURF_DRAWTELE 0x1000
#define SURF_DRAWWATER 0x2000
#define VERTEXSIZE 7
#define MAX_SKINS 32
#define MAXALIASVERTS 3984 //Baker: 1024 is GLQuake original limit. 2000 is WinQuake original limit. Baker 2000->3984
#define MAXALIASFRAMES 1024 //spike -- was 256
#define MAXALIASTRIS 4096 //ericw -- was 2048
#define EF_ROCKET 1 // leave a trail
#define EF_GRENADE 2 // leave a trail
#define EF_GIB 4 // leave a trail
#define EF_ROTATE 8 // rotate (bonus items)
#define EF_TRACER 16 // green split trail
#define EF_ZOMGIB 32 // small blood trail
#define EF_TRACER2 64 // orange split trail + rotate
#define EF_TRACER3 128 // purple trail
#define MF_HOLEY (1u<<14) // MarkV/QSS -- make index 255 transparent on mdl's
#define MOD_NOSHADOW 512 //don't cast a shadow
#define MOD_FBRIGHTHACK 1024 //when fullbrights are disabled, use a hack to render this model brighter
#define TEXPREF_NONE 0x0000
#define TEXPREF_MIPMAP 0x0001 // generate mipmaps
 // TEXPREF_NEAREST and TEXPREF_LINEAR aren't supposed to be ORed with TEX_MIPMAP
#define TEXPREF_LINEAR 0x0002 // force linear
#define TEXPREF_NEAREST 0x0004 // force nearest
#define TEXPREF_ALPHA 0x0008 // allow alpha
#define TEXPREF_PAD 0x0010 // allow padding
#define TEXPREF_PERSIST 0x0020 // never free
#define TEXPREF_OVERWRITE 0x0040 // overwrite existing same-name texture
#define TEXPREF_NOPICMIP 0x0080 // always load full-sized
#define TEXPREF_FULLBRIGHT 0x0100 // use fullbright mask palette
#define TEXPREF_NOBRIGHT 0x0200 // use nobright mask palette
#define TEXPREF_CONCHARS 0x0400 // use conchars palette
#define TEXPREF_WARPIMAGE 0x0800 // resize this texture when warpimagesize changes
//johnfitz -- extra flags for rendering
#define MOD_NOLERP 256 // don't lerp when animating
#define MOD_NOSHADOW 512 // don't cast a shadow
#define MOD_FBRIGHTHACK 1024 // when fullbrights are disabled, use a hack to render this model brighter

#define MAX_CACHED_PICS 128 // draw.c

#define MAXLEFTCLIPEDGES 100 // r_draw.c
#define FULLY_CLIPPED_CACHED 0x80000000
#define FRAMECOUNT_MASK 0x7FFFFFFF

#define BLINK_HZ ((1000 / 70) * 16) // vgatext.c

#define MAX_BMODEL_VERTS 2500 // CyanBun96: was 500 // r_bsp.c
#define MAX_BMODEL_EDGES 5000 // was 1000

#define NUMVERTEXNORMALS 162 // r_alias.c
#define LIGHT_MIN 5

#define GUARDSIZE 4 // d_surf.c

#define ANIM_CYCLE 2 // model.c
#define MAX_MOD_KNOWN 4096 // johnfitz -- was 512
#define NL_PRESENT 0 // values for model_t's needload
#define NL_NEEDS_LOADED 1
#define NL_UNREFERENCED 2
#define VISPATCH_HEADER_LEN 36

#define CON_TEXTSIZE 65536 // console.c
#define NUM_CON_TIMES 4
#define MAXCMDLINE 256
#define MAXGAMEDIRLEN 1000
#define MAXPRINTMSG 4096

#define MAX_STACK_DEPTH 64 // pr_exec.c
#define LOCALSTACK_SIZE 2048

#define CRC_INIT_VALUE 0xffff // crc.c
#define CRC_XOR_VALUE 0x0000

#define MAIN_ITEMS 5 // menu.c
#define SINGLEPLAYER_ITEMS 3
#define MAX_SAVEGAMES 12
#define MULTIPLAYER_ITEMS 3
#define NUM_SETUP_CMDS 5
#define SLIDER_RANGE 10
#define NUMCOMMANDS (s32)(sizeof(bindnames)/sizeof(bindnames[0]))
#define StartingGame (m_multiplayer_cursor == 1)
#define JoiningGame (m_multiplayer_cursor == 0)
#define SerialConfig (m_net_cursor == 0)
#define DirectConfig (m_net_cursor == 1)
#define IPXConfig (m_net_cursor == 2)
#define TCPIPConfig (m_net_cursor == 3)
#define NUM_HELP_PAGES 6
#define NUM_LANCONFIG_CMDS 3
#define NUM_GAMEOPTIONS 9

#define MOVE_EPSILON 0.01 // sv_phys.c
#define STOP_EPSILON 0.1
#define MAX_CLIP_PLANES 5
#define STEPSIZE 18

#define SKY_SPAN_SHIFT 5 // d_sky.c
#define SKY_SPAN_MAX (1 << SKY_SPAN_SHIFT)

#define MAX_FORWARD 6 // sv_user.c

#define PAK0_COUNT 339 // id1/pak0.pak - v1.0x // common.c
#define PAK0_CRC_V100 13900 // id1/pak0.pak - v1.00
#define PAK0_CRC_V101 62751 // id1/pak0.pak - v1.01
#define PAK0_CRC_V106 32981 // id1/pak0.pak - v1.06
#define PAK0_CRC (PAK0_CRC_V106)
#define PAK0_COUNT_V091 308 // id1/pak0.pak - v0.91/0.92, not supported
#define PAK0_CRC_V091 28804 // id1/pak0.pak - v0.91/0.92, not supported
#define VA_NUM_BUFFS 4
#define VA_BUFFERLEN 1024
#define MAX_FILES_IN_PACK 2048
#define LOADFILE_ZONE 0
#define LOADFILE_HUNK 1
#define LOADFILE_TEMPHUNK 2
#define LOADFILE_CACHE 3
#define LOADFILE_STACK 4
#define LOADFILE_MALLOC 5

#define MAX_FIELD_LEN 64 // pr_edict.c
#define GEFV_CACHESIZE 2
#define PR_STRING_ALLOCSLOTS 256

#define STEPSIZE 18 // sv_move.c
#define DI_NODIR -1

#define STRINGTEMP_BUFFERS 16 // pr_cmds.c
#define STRINGTEMP_LENGTH 1024
#define RETURN_EDICT(e) (((s32 *)pr_globals)[OFS_RETURN] = EDICT_TO_PROG(e))
#define MSG_BROADCAST 0 // unreliable to all
#define MSG_ONE 1 // reliable to one (msg_entity)
#define MSG_ALL 2 // reliable to all
#define MSG_INIT 3 // write to the init string
#define MAX_CHECK 16

#ifdef _WIN32 // platform dependant (v)snprintf function names: // common.c
#define snprintf_func _snprintf
#define vsnprintf_func _vsnprintf
#else
#define snprintf_func snprintf
#define vsnprintf_func vsnprintf
#endif

#define NUMVERTEXNORMALS 162 // r_part.c
#define MAX_PARTICLES 2048
#define ABSOLUTE_MIN_PARTICLES 512
#define DPS_MAXSPANS MAXHEIGHT+1 // +1 for spanpackage marking end // d_polyse.c
#define SIGNON_SIZE 31500 // sv_main.c
#define MAX_TIMINGS 100 // r_misc.c
#define DYNAMIC_SIZE (4 * 1024 * 1024) // ericw -- was 512KB // zone.c
#define ZONEID 0x1d4a11
#define MINFRAGMENT 64
#define HUNK_SENTINEL 0x1df001ed
#define HUNKNAME_LEN 24
#define CACHENAME_LEN 32
#define NUM_MIPS 4 // d_init.c
#define PAINTBUFFER_SIZE 2048 // snd_mix.c
#define AREA_DEPTH 4 // world.c
#define AREA_NODES 32
#define SAVEGAME_VERSION 5 // host_cmd.c
#define EXTRA_EFRAGS 128 // r_efrag.c
#define SND_FILTERQUALITY_DEFAULT "1" // snd_dma.c
#define sound_nominal_clip_dist 1000.0
#define MAX_SFX 1024
#define SPAN_NORMAL 1 // d_edge.c
#define SPAN_SKYBOX 2
#define SPAN_TRANS 4
#define SPAN_CUTOUT 8
#endif
