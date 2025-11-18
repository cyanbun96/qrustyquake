#include "quakedef.h"
u8 color_mix_lut[256][256][FOG_LUT_LEVELS];                          // rgbtoi.c
s32 fog_lut_built = 0;
u8 lit_lut[LIT_LUT_RES*LIT_LUT_RES*LIT_LUT_RES];
u8 lit_lut_initialized = 0;
s32 fog_initialized = 0;                                              // d_fog.c
f32 fog_density;
u8 fog_pal_index;
u32 oldmodes[NUM_OLDMODES * 2] = {                                  // vid_sdl.c
        320, 240,       640, 480,       800, 600,
        320, 200,       320, 240,       640, 350,
        640, 400,       640, 480,       800, 600 };
SDL_Window *window;
SDL_Surface *windowSurface;
// scrbuffs[0] = screen; 1 = screentop; 2 = screenui; 3 = screensbar;
SDL_Surface *screen; // custom world palette, not scaled
SDL_Surface *screentop; // default palette, not scaled
SDL_Surface *screenui; // custom ui palette, scaled
SDL_Surface *screensbar; // custom ui palette, scaled
SDL_Surface *scrbuffs[4]; // contains the above four
s8 modelist[NUM_OLDMODES][8]; // "320x240" etc. for menus
u32 SDLWindowFlags;
u32 uiscale;
u32 vimmode;
s32 vid_modenum;
s32 vid_testingmode;
s32 vid_realmode;
f64 vid_testendtime;
u8 vid_curpal[256 * 3];
viddef_t vid; // global video state
s32 fadescreen;
u64 d_pzbuffer_size;
u8 r_foundtranswater, r_alphapass;                                   // r_main.c
void *colormap;
s32 r_outofsurfaces;
s32 r_outofedges;
bool r_dowarp, r_dowarpold, r_viewchanged;
mvertex_t *r_pcurrentvertbase;
s32 c_surf;
s32 r_maxsurfsseen, r_maxedgesseen;
s32 r_clipflags;
u8 *r_warpbuffer;
vec3_t vup, base_vup; // view origin
vec3_t vpn, base_vpn;
vec3_t vright, base_vright;
vec3_t r_origin;
refdef_t r_refdef; // screen size info
f32 xcenter, ycenter;
f32 xscale, yscale;
f32 xscaleinv, yscaleinv;
f32 xscaleshrink, yscaleshrink;
f32 aliasxscale, aliasyscale, aliasxcenter, aliasycenter;
s32 screenwidth;
f32 pixelAspect;
mplane_t screenedge[4];
s64 r_visframecount; // refresh flags
s64 r_framecount = 1; // so frame counts initialized to 0 don't match
s32 r_polycount;
s32 r_drawnpolycount;
s32 *pfrustum_indexes[4];
mleaf_t *r_viewleaf, *r_oldviewleaf;
texture_t *r_notexture_mip;
f32 r_aliastransition, r_resfudge;
s32 d_lightstylevalue[256]; // 8.8 fraction of base light value
f64 d_times[15];
s32 colored_aliaslight;
f32 d_sdivzstepu, d_tdivzstepu, d_zistepu;                           // d_vars.c
f32 d_sdivzstepv, d_tdivzstepv, d_zistepv;
f32 d_sdivzorigin, d_tdivzorigin, d_ziorigin;
s32 sadjust, tadjust, bbextents, bbextentt;
u8 *cacheblock;
s32 cachewidth, cacheheight;
u8 *d_viewbuffer;
s16 *d_pzbuffer;
u32 d_zwidth;
kbutton_t in_mlook, in_strafe;                                     // cl_input.c
s32 jaxis_move_x = 0;                                                // in_sdl.c
s32 jaxis_move_y = 0;
s32 jaxis_look_x = 0;
s32 jaxis_look_y = 0;
s32 jaxis_trig_1 = -32768;
s32 jaxis_trig_2 = -32768;
SDL_Joystick *joystick = 0;
edict_t *sv_player;                                                 // sv_user.c
f32 *angles;
f32 *origin;
f32 *velocity;
bool onground;
usercmd_t cmd;
client_static_t cls;                                                // cl_main.c
client_state_t cl;
entity_t cl_entities[MAX_EDICTS];
lightstyle_t cl_lightstyle[MAX_LIGHTSTYLES];
dlight_t cl_dlights[MAX_DLIGHTS];
s32 cl_numvisedicts;
entity_t *cl_visedicts[MAX_VISEDICTS];
entity_t cl_temp_entities[MAX_TEMP_ENTITIES];                       // cl_tent.c
beam_t cl_beams[MAX_BEAMS];
s32 safemode;                                                        // common.c
s8 com_token[1024];
s32 com_argc;
s8 **com_argv;
bool standard_quake = 1, rogue, hipnotic;
s32 msg_readcount;
bool msg_badread;
s32 com_filesize;
s8 com_gamedir[MAX_OSPATH];
s8 com_basedir[MAX_OSPATH];
s32 file_from_pak;
s16 (*BigShort) (s16 l);
s16 (*LittleShort) (s16 l);
s32 (*BigLong) (s32 l);
s32 (*LittleLong) (s32 l);
f32 (*BigFloat) (f32 l);
f32 (*LittleFloat) (f32 l);
bool con_forcedup; // because no entities to refresh                // console.c
s32 con_totallines; // total lines in console scrollback
s32 con_backscroll; // lines up from bottom to display
bool con_initialized;
s32 con_notifylines; // scan lines to clear for notify lines
surfcache_t *d_initial_rover;                                        // d_init.c
bool d_roverwrapped;
s32 d_minmip;
f32 d_scalemip[NUM_MIPS - 1];
s32 d_y_aspect_shift, d_pix_min, d_pix_max, d_pix_shift;           // d_modech.c
s32 d_vrectx, d_vrecty, d_vrectright_particle, d_vrectbottom_particle;
s32 d_scantable[MAXHEIGHT];
s16 *zspantable[MAXHEIGHT];
u32 sb_updates; // if >= vid.numpages, no update needed                // sbar.c
s32 sb_lines; // scan lines to draw
bool r_cache_thrash; // set if surface cache is thrashing            // d_surf.c
surfcache_t *sc_rover;
s32 lmonly; // render lightmap only, for lit water
s32 drawlayer = 0;                                                     // draw.c
qpic_t *draw_disc;
quakeparms_t host_parms;                                               // host.c
bool host_initialized; // 1 if into command execution
bool isDedicated;
f64 host_frametime;
f64 realtime; // without any filtering or bounding
s32 host_framecount;
client_t *host_client; // current client
u8 *host_basepal;
u8 *host_colormap;
s32 current_skill;                                                 // host_cmd.c
bool noclip_anglehack;
s8 key_lines[32][MAXCMDLINE];                                          // keys.c
s32 key_linepos;
s32 key_lastpress;
s32 edit_line = 0;
keydest_t key_dest;
s32 key_count; // incremented every key event
s8 *keybindings[256];
s8 chat_buffer[32];
bool team_message = 0;
s32 lwmark = 0;                                                      // d_scan.c
u8 *litwater_base;
bool pr_trace;                                                      // pr_exec.c
dfunction_t *pr_xfunction;
s32 pr_argc;
bool insubmodel; // current entity info                               // r_bsp.c
entity_t *currententity;
vec3_t modelorg, base_modelorg; // viewpoint reletive to currently rendering ent
vec3_t r_entorigin; // the currently rendering entity in world coordinates
s32 r_currentbkey;
s32 c_faceclip; // number of faces clipped                           // r_draw.c
clipplane_t view_clipplanes[4];
s32 sintable[SIN_BUFFER_SIZE];
s32 intsintable[SIN_BUFFER_SIZE];
f32 winquake_surface_liquid_alpha;
edge_t *last_pcheck[MAXHEIGHT]; // indexed by scanline v
edge_t *r_edges, *edge_p, *edge_max;                                 // r_edge.c
surf_t *surfaces, *surface_p, *surf_max;
edge_t *newedges[MAXHEIGHT];
edge_t *removeedges[MAXHEIGHT];
s32 r_currentkey;
u32 scr_fullupdate;                                                  // screen.c
u32 clearnotify;
vrect_t scr_vrect;
bool scr_disabled_for_loading;
f32 scr_centertime_off;
f32 scr_con_current;
hudstyle_t hudstyle;
mnode_t *r_pefragtopnode;                                           // r_efrag.c
vec3_t r_emins, r_emaxs;
dma_t *shm = NULL;                                                  // snd_dma.c
channel_t snd_channels[MAX_CHANNELS];
s32 total_channels;
s32 paintedtime; // sample PAIRS
s32 s_rawend;
portable_samplepair_t s_rawsamples[MAX_RAW_SAMPLES];
server_t sv;                                                        // sv_main.c
server_static_t svs;
bool pr_alpha_supported; //johnfitz                                // pr_edict.c
s32 pr_effects_mask; // only enable rerelease quad/penta dlights when applicable
dprograms_t *progs;
dfunction_t *pr_functions;
dstatement_t *pr_statements;
globalvars_t *pr_global_struct;
f32 *pr_globals; // same as pr_global_struct
s32 pr_edict_size; // in bytes
u16 pr_crc;
u8 r_skypixels[6][SKYBOX_MAX_SIZE*SKYBOX_MAX_SIZE];                   // r_sky.c
f32 skyspeed, skyspeed2;
f32 skytime;
u8 *r_skysource;
s32 r_skymade;
s8 skybox_name[1024]; // name of current skybox, or "" if no skybox
drawsurf_t r_drawsurf;                                               // r_surf.c
u8 lit_loaded = 0;
u8 worldpal[768]; // custom world palette, set with "worldpal"       // common.c
u8 worldcmap[64*256];
u8 uipal[768];
s8 worldpalname[MAX_OSPATH];
s8 worldcmapname[MAX_OSPATH];
s8 uipalname[MAX_OSPATH];
s32 refresh_palette = 0;                                               // view.c
vec3_t vec3_origin = {0,0,0};                                       // mathlib.h
sspan_t spans[MAXHEIGHT + 1];                                      // d_sprite.c
s32 r_dlightframecount;                                             // r_light.c
vec3_t r_pright, r_pup, r_ppn;                                       // r_part.c
spritedesc_t r_spritedesc;                                         // r_sprite.c
cmd_source_t cmd_source;                                                // cmd.c
f32 scale_for_mip;                                                   // d_edge.c
s32 miplevel;
f32 cur_ent_alpha = 1;                                             // d_polyse.c
vec3_t lightcolor; //johnfitz -- lit support via lordhavoc          // r_light.c
const s32 dither_s[4] = { // Unreal-style 2×2 dither matrix scaled   // d_scan.c
	(s32)(-0.25 * 65536.0f), // negative values move the sampling the other
	(s32)(-0.50 * 65536.0f), // way, for the sole purpose of making skill
	(s32)(-0.75 * 65536.0f), // gates not sample the texels from that weird
	(s32)(-0.00 * 65536.0f)  // unused right side of the texture
};
const s32 dither_t[4] = {
	(s32)(-0.00 * 65536.0f),
	(s32)(-0.75 * 65536.0f),
	(s32)(-0.50 * 65536.0f),
	(s32)(-0.25 * 65536.0f)
};
f32 r_avertexnormals[NUMVERTEXNORMALS][3] = {                       // r_alias.c
{-0.525731,0.0,0.850651},{-0.442863,0.238856,0.864188},
{-0.295242,0.0,0.955423},{-0.309017,0.5,0.809017},
{-0.162460,0.262866,0.951056},{0.0,0.0,1.0},
{0.0,0.850651,0.525731},{-0.147621,0.716567,0.681718},
{0.147621,0.716567,0.681718},{0.0,0.525731,0.850651},
{0.309017,0.5,0.809017},{0.525731,0.0,0.850651},
{0.295242,0.0,0.955423},{0.442863,0.238856,0.864188},
{0.162460,0.262866,0.951056},{-0.681718,0.147621,0.716567},
{-0.809017,0.309017,0.5},{-0.587785,0.425325,0.688191},
{-0.850651,0.525731,0.0},{-0.864188,0.442863,0.238856},
{-0.716567,0.681718,0.147621},{-0.688191,0.587785,0.425325},
{-0.5,0.809017,0.309017},{-0.238856,0.864188,0.442863},
{-0.425325,0.688191,0.587785},{-0.716567,0.681718,-0.147621},
{-0.5,0.809017,-0.309017},{-0.525731,0.850651,0.0},
{0.0,0.850651,-0.525731},{-0.238856,0.864188,-0.442863},
{0.0,0.955423,-0.295242},{-0.262866,0.951056,-0.162460},
{0.0,1.0,0.0},{0.0,0.955423,0.295242},
{-0.262866,0.951056,0.162460},{0.238856,0.864188,0.442863},
{0.262866,0.951056,0.162460},{0.5,0.809017,0.309017},
{0.238856,0.864188,-0.442863},{0.262866,0.951056,-0.162460},
{0.5,0.809017,-0.309017},{0.850651,0.525731,0.0},
{0.716567,0.681718,0.147621},{0.716567,0.681718,-0.147621},
{0.525731,0.850651,0.0},{0.425325,0.688191,0.587785},
{0.864188,0.442863,0.238856},{0.688191,0.587785,0.425325},
{0.809017,0.309017,0.5},{0.681718,0.147621,0.716567},
{0.587785,0.425325,0.688191},{0.955423,0.295242,0.0},
{1.0,0.0,0.0},{0.951056,0.162460,0.262866},
{0.850651,-0.525731,0.0},{0.955423,-0.295242,0.0},
{0.864188,-0.442863,0.238856},{0.951056,-0.162460,0.262866},
{0.809017,-0.309017,0.5},{0.681718,-0.147621,0.716567},
{0.850651,0.0,0.525731},{0.864188,0.442863,-0.238856},
{0.809017,0.309017,-0.5},{0.951056,0.162460,-0.262866},
{0.525731,0.0,-0.850651},{0.681718,0.147621,-0.716567},
{0.681718,-0.147621,-0.716567},{0.850651,0.0,-0.525731},
{0.809017,-0.309017,-0.5},{0.864188,-0.442863,-0.238856},
{0.951056,-0.162460,-0.262866},{0.147621,0.716567,-0.681718},
{0.309017,0.5,-0.809017},{0.425325,0.688191,-0.587785},
{0.442863,0.238856,-0.864188},{0.587785,0.425325,-0.688191},
{0.688191,0.587785,-0.425325},{-0.147621,0.716567,-0.681718},
{-0.309017,0.5,-0.809017},{0.0,0.525731,-0.850651},
{-0.525731,0.0,-0.850651},{-0.442863,0.238856,-0.864188},
{-0.295242,0.0,-0.955423},{-0.162460,0.262866,-0.951056},
{0.0,0.0,-1.0},{0.295242,0.0,-0.955423},
{0.162460,0.262866,-0.951056},{-0.442863,-0.238856,-0.864188},
{-0.309017,-0.5,-0.809017},{-0.162460,-0.262866,-0.951056},
{0.0,-0.850651,-0.525731},{-0.147621,-0.716567,-0.681718},
{0.147621,-0.716567,-0.681718},{0.0,-0.525731,-0.850651},
{0.309017,-0.5,-0.809017},{0.442863,-0.238856,-0.864188},
{0.162460,-0.262866,-0.951056},{0.238856,-0.864188,-0.442863},
{0.5,-0.809017,-0.309017},{0.425325,-0.688191,-0.587785},
{0.716567,-0.681718,-0.147621},{0.688191,-0.587785,-0.425325},
{0.587785,-0.425325,-0.688191},{0.0,-0.955423,-0.295242},
{0.0,-1.0,0.0},{0.262866,-0.951056,-0.162460},
{0.0,-0.850651,0.525731},{0.0,-0.955423,0.295242},
{0.238856,-0.864188,0.442863},{0.262866,-0.951056,0.162460},
{0.5,-0.809017,0.309017},{0.716567,-0.681718,0.147621},
{0.525731,-0.850651,0.0},{-0.238856,-0.864188,-0.442863},
{-0.5,-0.809017,-0.309017},{-0.262866,-0.951056,-0.162460},
{-0.850651,-0.525731,0.0},{-0.716567,-0.681718,-0.147621},
{-0.716567,-0.681718,0.147621},{-0.525731,-0.850651,0.0},
{-0.5,-0.809017,0.309017},{-0.238856,-0.864188,0.442863},
{-0.262866,-0.951056,0.162460},{-0.864188,-0.442863,0.238856},
{-0.809017,-0.309017,0.5},{-0.688191,-0.587785,0.425325},
{-0.681718,-0.147621,0.716567},{-0.442863,-0.238856,0.864188},
{-0.587785,-0.425325,0.688191},{-0.309017,-0.5,0.809017},
{-0.147621,-0.716567,0.681718},{-0.425325,-0.688191,0.587785},
{-0.162460,-0.262866,0.951056},{0.442863,-0.238856,0.864188},
{0.162460,-0.262866,0.951056},{0.309017,-0.5,0.809017},
{0.147621,-0.716567,0.681718},{0.0,-0.525731,0.850651},
{0.425325,-0.688191,0.587785},{0.587785,-0.425325,0.688191},
{0.688191,-0.587785,0.425325},{-0.955423,0.295242,0.0},
{-0.951056,0.162460,0.262866},{-1.0,0.0,0.0},
{-0.850651,0.0,0.525731},{-0.955423,-0.295242,0.0},
{-0.951056,-0.162460,0.262866},{-0.864188,0.442863,-0.238856},
{-0.951056,0.162460,-0.262866},{-0.809017,0.309017,-0.5},
{-0.864188,-0.442863,-0.238856},{-0.951056,-0.162460,-0.262866},
{-0.809017,-0.309017,-0.5},{-0.681718,0.147621,-0.716567},
{-0.681718,-0.147621,-0.716567},{-0.850651,0.0,-0.525731},
{-0.688191,0.587785,-0.425325},{-0.587785,0.425325,-0.688191},
{-0.425325,0.688191,-0.587785},{-0.425325,-0.688191,-0.587785},
{-0.587785,-0.425325,-0.688191},{-0.688191,-0.587785,-0.425325} };
