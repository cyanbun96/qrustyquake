#ifndef QTYPEDEFS_
#define QTYPEDEFS_
typedef unsigned char      u8;
typedef char               s8;
typedef unsigned short     u16;
typedef short              s16;
typedef unsigned int       u32;
typedef int                s32;
typedef unsigned long      u64;
typedef unsigned long long ul64;
typedef long               s64;
typedef float              f32;
typedef double             f64;

typedef struct { // specified by the host system                   // quakedef.h
	s8 *basedir;
	s8 *userdir;
	s8 *cachedir; // for development over ISDN lines
	s32 argc;
	s8 **argv;
	void *membase;
	s32 memsize;
} quakeparms_t;

typedef f32 vec_t;                                                 // q_stdinc.h
typedef vec_t vec3_t[3];
typedef vec_t vec4_t[4];
typedef vec_t vec5_t[5];

typedef struct sizebuf_s {                                           // common.h
	bool allowoverflow; // if 0, do a Sys_Error
	bool overflowed; // set to 1 if the buffer size failed
	u8 *data;
	s32 maxsize;
	s32 cursize;
} sizebuf_t;
typedef enum {
	CPE_NOTRUNC, // return parse error in case of overflow
	CPE_ALLOWTRUNC // truncate com_token in case of overflow
} cpe_mode;
typedef struct link_s {
	struct link_s *prev, *next;
} link_t;
typedef struct vec_header_t {
	size_t capacity;
	size_t size;
} vec_header_t;
typedef struct {
	s8 name[MAX_QPATH];
	s32 filepos, filelen;
} packfile_t;
typedef struct pack_s {
	s8 filename[MAX_OSPATH];
	s32 handle;
	s32 numfiles;
	packfile_t *files;
} pack_t;
typedef struct searchpath_s {
	u32 path_id; // identifier assigned to the game directory
	s8 filename[MAX_OSPATH];
	pack_t *pack; // only one of filename / pack will be used
	struct searchpath_s *next;
} searchpath_t;
typedef struct _fshandle_t {
	FILE *file;
	bool pak; /* is the file read from a pak */
	s64 start; /* file or data start position */
	s64 length; /* file or data size */
	s64 pos; /* current position relative to start */
} fshandle_t;

typedef struct {                                                    // bspfile.h
	s32 fileofs, filelen;
} lump_t;
typedef struct {
	f32 mins[3], maxs[3];
	f32 origin[3];
	s32 headnode[MAX_MAP_HULLS];
	s32 visleafs; // not including the solid leaf 0
	s32 firstface, numfaces;
} dmodel_t;
typedef struct {
	s32 version;
	lump_t lumps[HEADER_LUMPS];
} dheader_t;
typedef struct {
	s32 nummiptex;
	s32 dataofs[4]; // [nummiptex]
} dmiptexlump_t;
typedef struct miptex_s {
	s8 name[16];
	unsigned width, height;
	unsigned offsets[MIPLEVELS]; // four mip maps stored
} miptex_t;
typedef struct {
	f32 point[3];
} dvertex_t;
typedef struct {
	f32 normal[3];
	f32 dist;
	s32 type; // PLANE_X - PLANE_ANYZ ?remove? trivial to regenerate
} dplane_t;
typedef struct {
	s32 planenum;
	s16 children[2]; // negative numbers are -(leafs+1), not nodes
	s16 mins[3]; // for sphere culling
	s16 maxs[3];
	u16 firstface;
	u16 numfaces; // counting both sides
} dnode_t;
typedef struct {
	s32 planenum;
	s32 children[2]; // negative numbers are -(leafs+1), not nodes
	s16 mins[3]; // for sphere culling
	s16 maxs[3];
	u32 firstface;
	u32 numfaces; // counting both sides
} dl1node_t;
typedef struct {
	s32 planenum;
	s32 children[2]; // negative numbers are -(leafs+1), not nodes
	f32 mins[3]; // for sphere culling
	f32 maxs[3];
	u32 firstface;
	u32 numfaces; // counting both sides
} dl2node_t;
typedef struct {
	s32 planenum;
	s16 children[2]; // negative numbers are contents
} dclipnode_t;
typedef struct texinfo_s {
	f32 vecs[2][4]; // [s/t][xyz offset]
	s32 miptex;
	s32 flags;
} texinfo_t;
typedef struct { // note that edge 0 is never used, because negative edge nums
		 // are used for counterclockwise use of the edge in a face
	u16 v[2]; // vertex numbers
} dedge_t;
typedef struct {
	u32 v[2]; // vertex numbers
} dledge_t;
typedef struct {
	s16 planenum;
	s16 side;
	s32 firstedge; // we must support > 64k edges
	s16 numedges;
	s16 texinfo;
	u8 styles[MAXLIGHTMAPS]; // lighting info
	s32 lightofs; // start of [numstyles*surfsize] samples
} dface_t;
typedef struct {
	s32 planenum;
	s32 side;
	s32 firstedge; // we must support > 64k edges
	s32 numedges;
	s32 texinfo;
	// lighting info
	u8 styles[MAXLIGHTMAPS];
	s32 lightofs; // start of [numstyles*surfsize] samples
} dlface_t;
typedef struct { // leaf 0 is the generic CONTENTS_SOLID leaf, used for all
		 // solid areas, all other leafs need visibility info
	s32 contents;
	s32 visofs; // -1 = no visibility info
	s16 mins[3]; // for frustum culling
	s16 maxs[3];
	u16 firstmarksurface;
	u16 nummarksurfaces;
	u8 ambient_level[NUM_AMBIENTS];
} dleaf_t;
typedef struct {
	s32 contents;
	s32 visofs; // -1 = no visibility info
	s16 mins[3]; // for frustum culling
	s16 maxs[3];
	u32 firstmarksurface;
	u32 nummarksurfaces;
	u8 ambient_level[NUM_AMBIENTS];
} dl1leaf_t;
typedef struct {
	s32 contents;
	s32 visofs; // -1 = no visibility info
	f32 mins[3]; // for frustum culling
	f32 maxs[3];
	u32 firstmarksurface;
	u32 nummarksurfaces;
	u8 ambient_level[NUM_AMBIENTS];
} dl2leaf_t;

typedef struct cvar_s cvar_t;                                          // cvar.h
typedef void (*cvarcallback_t)(cvar_t *);
struct cvar_s {
	const s8 *name;
	const s8 *string;
	u32 flags;
	f32  value;
	const s8 *default_string; // remember defaults for reset function
	cvarcallback_t callback;
	cvar_t  *next;
};

// d*_t structures are on-disk representations
// m*_t structures are in-memory

typedef struct { f32 l, a, b; } lab_t;                               // rgbtoi.h

typedef struct vrect_s {                                                // vid.h
	s32 x,y,width,height; struct vrect_s *pnext;
} vrect_t;
typedef struct {
	u8 *buffer; // invisible buffer
	u8 *colormap; // 256 * VID_GRADES size
	u32 width;
	u32 height;
	f32 aspect; // width / height -- < 0 is taller than wide
	u32 numpages;
	u32 recalc_refdef; // if 1, recalc vid-based stuff
} viddef_t;

typedef struct {                                                      // model.h
	vec3_t position;
} mvertex_t;
typedef struct texture_s {
	s8 name[16];
	unsigned width, height;
	unsigned shift; // Q64
	bool update_warp; //johnfitz -- update warp this frame
	struct msurface_s *texturechains[2]; // for texture chains
	s32 anim_total; // total tenths in sequence ( 0 = no)
	s32 anim_min, anim_max; // time for this frame min <=time< max
	struct texture_s *anim_next; // in the animation sequence
	struct texture_s *alternate_anims; // bmodels in frmae 1 use these
	unsigned offsets[MIPLEVELS]; // four mip maps stored
} texture_t;
typedef struct {
	f32 vecs[2][4];
	f32 mipadjust;
	texture_t *texture;
	s32 flags;
} mtexinfo_t;
typedef struct mplane_s {
	vec3_t normal;
	f32 dist;
	u8 type; // for texture axis selection and fast side tests
	u8 signbits; // signx + signy<<1 + signz<<1
	u8 pad[2];
} mplane_t;
typedef struct msurface_s {
	s32 visframe; // should be drawn when node is crossed
	f32 mins[3]; // johnfitz -- for frustum culling
	f32 maxs[3]; // johnfitz -- for frustum culling
	mplane_t *plane;
	s32 flags;
	s32 firstedge; // look up in model->surfedges[], negative numbers
	s32 numedges; // are backwards edges
	s16 texturemins[2];
	s16 extents[2];
	s32 light_s, light_t; // gl lightmap coordinates
	struct msurface_s *texturechain;
	mtexinfo_t *texinfo;
	s32 vbo_firstvert; // index of surface's first vert in VBO lighting info
	s32 dlightframe;
	s32 dlightbits;
	s32 lightmaptexturenum;
	u8 styles[MAXLIGHTMAPS];
	s32 cached_light[MAXLIGHTMAPS]; // values currently used in lightmap
	bool cached_dlight; // 1 if dynamic light in cache
	u8 *samples; // [numstyles*surfsize]
	struct surfcache_s *cachespots[MIPLEVELS]; // surface generation data
} msurface_t;

typedef struct { // viewmodel lighting                              // r_local.h
	s32 ambientlight;
	s32 shadelight;
	f32 *plightvec;
} alight_t;
typedef struct bedge_s // clipped bmodel edges
{
	mvertex_t *v[2];
	struct bedge_s *pnext;
} bedge_t;
typedef struct {
	f32 fv[3]; // viewspace x, y
} auxvert_t;
typedef struct clipplane_s {
	vec3_t normal;
	f32 dist;
	struct clipplane_s *next;
	u8 leftedge;
	u8 rightedge;
	u8 reserved[2];
} clipplane_t;
typedef struct btofpoly_s {
	s32 clipflags;
	msurface_t *psurf;
} btofpoly_t;
typedef enum {
	TEXTYPE_DEFAULT,
	TEXTYPE_CUTOUT,
	TEXTYPE_SKY,
	TEXTYPE_LAVA,
	TEXTYPE_SLIME,
	TEXTYPE_TELE,
	TEXTYPE_WATER,
	TEXTYPE_COUNT,
	TEXTYPE_FIRSTLIQUID = TEXTYPE_LAVA,
	TEXTYPE_LASTLIQUID = TEXTYPE_WATER,
	TEXTYPE_NUMLIQUIDS = TEXTYPE_LASTLIQUID + 1 - TEXTYPE_FIRSTLIQUID,
} textype_t;

typedef enum { ST_SYNC=0, ST_RAND } synctype_t;                    // modelgen.h
typedef enum { ALIAS_SINGLE=0, ALIAS_GROUP } aliasframetype_t;
typedef enum { ALIAS_SKIN_SINGLE=0, ALIAS_SKIN_GROUP } aliasskintype_t;
typedef struct {
	s32 ident;
	s32 version;
	vec3_t scale;
	vec3_t scale_origin;
	f32 boundingradius;
	vec3_t eyeposition;
	s32 numskins;
	s32 skinwidth;
	s32 skinheight;
	s32 numverts;
	s32 numtris;
	s32 numframes;
	synctype_t synctype;
	s32 flags;
	f32 size;
} mdl_t;
typedef struct {
	s32 onseam;
	s32 s;
	s32 t;
} stvert_t;
typedef struct dtriangle_s {
	s32 facesfront;
	s32 vertindex[3];
} dtriangle_t;
typedef struct {
	u8 v[3];
	u8 lightnormalindex;
} trivertx_t;
typedef struct {
	trivertx_t bboxmin; // lightnormal isn't used
	trivertx_t bboxmax; // lightnormal isn't used
	s8 name[16]; // frame name from grabbing
} daliasframe_t;
typedef struct {
	s32 numframes;
	trivertx_t bboxmin; // lightnormal isn't used
	trivertx_t bboxmax; // lightnormal isn't used
} daliasgroup_t;
typedef struct {
	s32 numskins;
} daliasskingroup_t;
typedef struct {
	f32 interval;
} daliasinterval_t;
typedef struct {
	f32 interval;
} daliasskininterval_t;
typedef struct {
	aliasframetype_t type;
} daliasframetype_t;
typedef struct {
	aliasskintype_t type;
} daliasskintype_t;

typedef enum {SPR_SINGLE=0,SPR_GROUP,SPR_ANGLED}spriteframetype_t; // spritegn.h
typedef struct {
	s32 ident;
	s32 version;
	s32 type;
	f32 boundingradius;
	s32 width;
	s32 height;
	s32 numframes;
	f32 beamlength;
	synctype_t synctype;
} dsprite_t;
typedef struct {
	s32 origin[2];
	s32 width;
	s32 height;
} dspriteframe_t;
typedef struct { s32 numframes; } dspritegroup_t;
typedef struct { f32 interval; } dspriteinterval_t;
typedef struct { spriteframetype_t type; } dspriteframetype_t;

typedef struct cache_user_s                                            // zone.h
{
	void *data;
} cache_user_t;

typedef struct {                                                   // protocol.h
	vec3_t  origin;
	vec3_t  angles;
	u16 modelindex; //johnfitz -- was s32
	u16 frame;  //johnfitz -- was s32
	u8 colormap; //johnfitz -- was s32
	u8 skin;  //johnfitz -- was s32
	u8 alpha;  //johnfitz -- added
	u8 scale;  //Quakespasm: for model scale support.
	s32  effects;
} entity_state_t;
typedef struct {
	vec3_t viewangles;
	f32 forwardmove; // intended velocities
	f32 sidemove;
	f32 upmove;
} usercmd_t;

typedef struct entity_s {                                            // render.h
	bool forcelink; // model changed
	s32 update_type;
	entity_state_t baseline; // to fill in defaults in updates
	f64 msgtime; // time of last update
	vec3_t msg_origins[2]; // last two updates (0 is newest)
	vec3_t origin;
	vec3_t msg_angles[2]; // last two updates (0 is newest)
	vec3_t angles;
	struct model_s *model; // NULL = no model
	struct efrag_s *efrag; // linked list of efrags
	s32 frame;
	f32 syncbase; // for client-side animations
	u8 *colormap;
	s32 effects; // light, particles, etc
	s32 skinnum; // for Alias models
	s32 visframe; // last frame this entity was found in an active leaf
	s32 dlightframe; // dynamic lighting
	s32 dlightbits;
	s32 trivial_accept;
	struct mnode_s *topnode; // for bmodels, first world node that splits
				 // bmodel, or NULL if not split
	u8 alpha; //johnfitz -- alpha
	u8 scale;
	u8 lerpflags; //johnfitz -- lerping
	f32 lerpstart; //johnfitz -- animation lerping
	f32 lerptime; //johnfitz -- animation lerping
	f32 lerpfinish; //johnfitz -- lerping
	s16 previouspose; //johnfitz -- animation lerping
	s16 currentpose; //johnfitz -- animation lerping
	s16 futurepose; //johnfitz -- animation lerping
	f32 movelerpstart; //johnfitz -- transform lerping
	vec3_t previousorigin; //johnfitz -- transform lerping
	vec3_t currentorigin; //johnfitz -- transform lerping
	vec3_t previousangles; //johnfitz -- transform lerping
	f32 traildelay; // time left until next particle trail update
	vec3_t trailorg; // previous particle trail point
} entity_t;
typedef struct efrag_s {
	struct mleaf_s *leaf;
	struct efrag_s *leafnext;
	struct entity_s *entity;
	struct efrag_s *entnext;
} efrag_t;
typedef struct {
	vrect_t vrect; // subwindow in video for refresh
		       // FIXME: not need vrect next field here?
	vrect_t aliasvrect; // scaled Alias version
	s64 vrectright, vrectbottom; // right & bottom screen coords
	s64 aliasvrectright, aliasvrectbottom; // scaled Alias versions
	f32 vrectrightedge; // rightmost right edge we care about,
			    // for use in edge list
	f32 fvrectx, fvrecty; // for floating-point compares
	f32 fvrectx_adj, fvrecty_adj; // left and top edges, for clamping
	s64 vrect_x_adj_shift20; //(vrect.x + 0.5 - epsilon) << 20
	s64 vrectright_adj_shift20; //(vrectright + 0.5 - epsilon) << 20
	f32 fvrectright_adj, fvrectbottom_adj; // right and bottom edges,
					       // for clamping
	f32 fvrectright; // rightmost edge, for Alias clamping
	f32 fvrectbottom; // bottommost edge, for Alias clamping
	f32 horizontalFieldOfView; // at Z = 1.0, this many X is visible 
				   // 2.0 = 90 degrees
	f32 xOrigin; // should probably allways be 0.5
	f32 yOrigin; // between be around 0.3 to 0.5

	vec3_t vieworg;
	vec3_t viewangles;

	f32 fov_x, fov_y;

	s32 ambientlight;
} refdef_t;

typedef struct {                                                      // model.h
	aliasskintype_t type;
	void *pcachespot;
	s32 skin;
} maliasskindesc_t;
typedef struct {
	s32 numskins;
	s32 intervals;
	maliasskindesc_t skindescs[1];
} maliasskingroup_t;
typedef enum { // ericw -- each texture has two chains, so we can clear the
	chain_world = 0, // model chains without affecting the world
	chain_model = 1
} texchain_t;
typedef struct {
	u32 v[2];
	u32 cachededgeoffset;
} medge_t;
typedef struct mnode_s {
	s32 contents; // common with leaf. 0, to differentiate from leafs
	s32 visframe; // node needs to be traversed if current
	f32 minmaxs[6]; // for bounding box culling
	struct mnode_s *parent;
	mplane_t *plane; // node specific
	struct mnode_s *children[2];
	u32 firstsurface;
	u32 numsurfaces;
} mnode_t;
typedef struct mleaf_s {
	s32 contents; // common with node. Will be a negative contents number
	s32 visframe; // node needs to be traversed if current
	f32 minmaxs[6]; // for bounding box culling
	struct mnode_s *parent;
	u8 *compressed_vis; // leaf specific
	efrag_t *efrags;
	msurface_t **firstmarksurface;
	s32 nummarksurfaces;
	s32 key; // BSP sequence number for leaf's contents
	u8 ambient_sound_level[NUM_AMBIENTS];
} mleaf_t;
typedef struct mclipnode_s { //johnfitz -- for clipnodes>32k
	s32 planenum;
	s32 children[2]; // negative numbers are contents
} mclipnode_t;
typedef struct {
	s32 planenum;
	s16 children[2]; // negative numbers are contents
} dsclipnode_t;
typedef struct {
	s32 planenum;
	s32 children[2]; // negative numbers are contents
} dlclipnode_t;
typedef struct {
	mclipnode_t *clipnodes; //johnfitz -- was dclipnode_t
	mplane_t *planes;
	s32 firstclipnode;
	s32 lastclipnode;
	vec3_t clip_mins;
	vec3_t clip_maxs;
} hull_t;
typedef struct mspriteframe_s {
	s32 width, height;
	f32 up, down, left, right;
	f32 smax, tmax; //johnfitz -- image might be padded
	u8 pixels[4];
} mspriteframe_t;
typedef struct {
	s32 numframes;
	f32 *intervals;
	mspriteframe_t *frames[1];
} mspritegroup_t;
typedef struct {
	spriteframetype_t type;
	mspriteframe_t *frameptr;
} mspriteframedesc_t;
typedef struct {
	s32 type;
	s32 maxwidth;
	s32 maxheight;
	s32 numframes;
	f32 beamlength;
	void *cachespot;
	mspriteframedesc_t frames[1];
} msprite_t;
// Alias models are position independent, so the cache manager can move them.
typedef struct aliasmesh_s{//from RMQEngine, split out to keep vertex sizes down
	f32 st[2];
	u16 vertindex;
} aliasmesh_t;
typedef struct meshxyz_s {
	u8 xyz[4];
	s8 normal[4];
} meshxyz_t;
typedef struct meshst_s {
	f32 st[2];
} meshst_t; // RMQEngine end
typedef struct {
	aliasframetype_t type;
	s32 firstpose;
	s32 numposes;
	f32 interval;
	trivertx_t bboxmin;
	trivertx_t bboxmax;
	s32 frame;
	s8 name[16];
} maliasframedesc_t;
typedef struct {
	trivertx_t bboxmin;
	trivertx_t bboxmax;
	s32 frame;
} maliasgroupframedesc_t;
typedef struct {
	s32 numframes;
	s32 intervals;
	maliasgroupframedesc_t frames[1];
} maliasgroup_t;
typedef struct mtriangle_s {
	s32 facesfront;
	s32 vertindex[3];
} mtriangle_t;
typedef struct {
	s32 ident;
	s32 version;
	vec3_t scale;
	vec3_t scale_origin;
	f32 boundingradius;
	vec3_t eyeposition;
	s32 numskins;
	s32 skinwidth;
	s32 skinheight;
	s32 numverts;
	s32 numtris;
	s32 numframes;
	synctype_t synctype;
	s32 flags;
	f32 size;
	//ericw -- used to populate vbo
	s32 numverts_vbo; // number of verts with unique x,y,z,s,t
	intptr_t meshdesc; // offset into extradata: numverts_vbo aliasmesh_t
	s32 numindexes;
	intptr_t indexes; // offset into extradata: numindexes unsigned shorts
	intptr_t vertexes; // offset into extradata: numposes*vertsperframe trivertx_t
	s32 numposes;
	s32 poseverts;
	s32 posedata; // numposes*poseverts trivert_t
	s32 commands; // gl command list with embedded s/t
	s32 texels[MAX_SKINS]; // only for player skins
	s32 model;
	s32 stverts;
	s32 skindesc;
	s32 triangles;
	maliasframedesc_t frames[1]; // variable sized
} aliashdr_t;
typedef enum {mod_brush, mod_sprite, mod_alias} modtype_t;
typedef struct model_s {
	f32 radius;
	s8 name[MAX_QPATH];
	u32 path_id; // path id of game dir that this model came from
	bool needload; // bmodels and sprites don't cache normally
	modtype_t type;
	s32 numframes;
	synctype_t synctype;
	s32 flags;
	vec3_t mins, maxs; // volume occupied by the model graphics
	vec3_t ymins, ymaxs; //johnfitz -- bounds for entities with nonzero yaw
	vec3_t rmins, rmaxs; //johnfitz -- bounds for entities with nonzero pitch or roll
	bool clipbox; // solid volume for clipping
	vec3_t clipmins, clipmaxs;
	s32 firstmodelsurface, nummodelsurfaces; // brush model
	s32 numsubmodels;
	dmodel_t *submodels;
	s32 numplanes;
	mplane_t *planes;
	s32 numleafs; // number of visible leafs, not counting 0
	mleaf_t *leafs;
	s32 numvertexes;
	mvertex_t *vertexes;
	s32 numedges;
	medge_t *edges;
	s32 numnodes;
	mnode_t *nodes;
	s32 numtexinfo;
	mtexinfo_t *texinfo;
	s32 numsurfaces;
	msurface_t *surfaces;
	s32 numsurfedges;
	s32 *surfedges;
	s32 numclipnodes;
	mclipnode_t *clipnodes; //johnfitz -- was dclipnode_t
	s32 nummarksurfaces;
	msurface_t **marksurfaces;
	hull_t hulls[MAX_MAP_HULLS];
	s32 numtextures;
	texture_t **textures;
	u8 *visdata;
	u8 *lightdata;
	s8 *entities;
	bool viswarn; // for Mod_DecompressVis()
	s32 bspversion;
	s32 contentstransparent; //spike -- added this so we can disable glitchy
				 //wateralpha where its not supported.
	bool haslitwater;
	// alias model
	s32 vboindexofs; // offset in vbo of the hdr->numindexes unsigned shorts
	s32 vboxyzofs; // offset of hdr->numposes*hdr->numverts_vbo meshxyz_t
	s32 vbostofs; // offset in vbo of hdr->numverts_vbo meshst_t
	cache_user_t cache;//extra model data, only access through Mod_Extradata
} model_t;

typedef struct espan_s {                                           // r_shared.h
	s64 u, v, count;
	struct espan_s *pnext;
} espan_t;
typedef struct surf_s {
	struct surf_s *next; // active surface stack in r_edge.c
	struct surf_s *prev; // used in r_edge.c for active surf stack
	struct espan_s *spans; // pointer to linked list of spans to draw
	s32 key; // sorting key (BSP order)
	s64 last_u; // set during tracing
	s32 spanstate; // 0 = not in span
		       // 1 = in span
		       // -1 = in inverted span (end before start)
	s32 flags; // currentface flags
	void *data; // associated data like msurface_t
	entity_t *entity;
	f32 nearzi; // nearest 1/z on surface, for mipmapping
	bool insubmodel;
	f32 d_ziorigin, d_zistepu, d_zistepv;
	s32 pad[2]; // to 64 bytes
} surf_t;
typedef struct edge_s {
	s64 u;
	s64 u_step;
	struct edge_s *prev, *next;
	u16 surfs[2];
	struct edge_s *nextremove;
	f32 nearzi;
	medge_t *owner;
} edge_t;

typedef struct {                                                    // d_iface.h
	f32 u, v;
	f32 s, t;
	f32 zi;
} emitpoint_t;
typedef enum {
	pt_static, pt_grav, pt_slowgrav, pt_fire,
	pt_explode, pt_explode2, pt_blob, pt_blob2
} ptype_t;
typedef struct particle_s {
	// driver-usable fields
	vec3_t org;
	f32 color;
	// drivers never touch the following fields
	struct particle_s *next;
	vec3_t vel;
	f32 ramp;
	f32 die;
	ptype_t type;
} particle_t;
typedef struct polyvert_s {
	f32 u, v, zi, s, t;
} polyvert_t;
typedef struct polydesc_s {
	s32 numverts;
	f32 nearzi;
	msurface_t *pcurrentface;
	polyvert_t *pverts;
} polydesc_t;
typedef struct finalvert_s {
	s32 v[6]; // u, v, s, t, l, 1/z
	s32 flags;
	f32 reserved;
} finalvert_t;
typedef struct {
	void *pskin;
	maliasskindesc_t *pskindesc;
	s32 skinwidth;
	s32 skinheight;
	mtriangle_t *ptriangles;
	finalvert_t *pfinalverts;
	s32 numtriangles;
	s32 drawtype;
	s32 seamfixupX16;
} affinetridesc_t;
typedef struct {
	f32 u, v, zi, color;
} screenpart_t;
typedef struct {
	s32 nump;
	emitpoint_t *pverts; // there's room for an extra element at [nump], if the driver wants to duplicate element [0] at element [nump] to avoid dealing with wrapping
	mspriteframe_t *pspriteframe;
	vec3_t vup, vright, vpn; // in worldspace
	f32 nearzi;
} spritedesc_t;
typedef struct {
	s32 u, v;
	f32 zi;
	s32 color;
} zpointdesc_t;
typedef struct {
	u8 *surfdat; // destination for generated surface
	s32 rowbytes; // destination logical width in bytes
	msurface_t *surf; // description for surface to generate
	s32 lightadj[MAXLIGHTMAPS];
	// adjust for lightmap levels for dynamic lighting
	texture_t *texture; // corrected for animating textures
	s32 surfmip; // mipmapped ratio of surface texels / world pixels
	s32 surfwidth; // in mipmapped texels
	s32 surfheight; // in mipmapped texels
} drawsurf_t;

typedef enum {key_game, key_console, key_message, key_menu}keydest_t;  // keys.h
typedef struct { s8 *name; s32 keynum; } keyname_t;

#ifndef _WIN32                                                      // net_sys.h
typedef s32 sys_socket_t;
#else
typedef u_long in_addr_t;
typedef SOCKET sys_socket_t;
#endif

struct qsockaddr {                                                 // net_defs.h
	s16 qsa_family;
	u8 qsa_data[14];
};
typedef struct qsocket_s {
	struct qsocket_s *next;
	f64 connecttime;
	f64 lastMessageTime;
	f64 lastSendTime;
	bool disconnected;
	bool canSend;
	bool sendNext;
	s32 driver;
	s32 landriver;
	sys_socket_t socket;
	void *driverdata;
	u32 ackSequence;
	u32 sendSequence;
	u32 unreliableSendSequence;
	s32 sendMessageLength;
	u8 sendMessage [NET_MAXMESSAGE];
	u32 receiveSequence;
	u32 unreliableReceiveSequence;
	s32 receiveMessageLength;
	u8 receiveMessage [NET_MAXMESSAGE];
	struct qsockaddr addr;
	s8 address[NET_NAMELEN];

} qsocket_t;
typedef struct {
	const s8 *name;
	bool initialized;
	sys_socket_t controlSock;
	sys_socket_t (*Init) ();
	void (*Shutdown) ();
	void (*Listen) (bool state);
	sys_socket_t (*Open_Socket) (s32 port);
	s32 (*Close_Socket) (sys_socket_t socketid);
	s32 (*Connect) (sys_socket_t socketid, struct qsockaddr *addr);
	sys_socket_t (*CheckNewConnections) ();
	s32 (*Read) (sys_socket_t socketid, u8 *buf, s32 len, struct qsockaddr *addr);
	s32 (*Write) (sys_socket_t socketid, u8 *buf, s32 len, struct qsockaddr *addr);
	s32 (*Broadcast) (sys_socket_t socketid, u8 *buf, s32 len);
	const s8 * (*AddrToString) (struct qsockaddr *addr);
	s32 (*StringToAddr) (const s8 *string, struct qsockaddr *addr);
	s32 (*GetSocketAddr) (sys_socket_t socketid, struct qsockaddr *addr);
	s32 (*GetNameFromAddr) (struct qsockaddr *addr, s8 *name);
	s32 (*GetAddrFromName) (const s8 *name, struct qsockaddr *addr);
	s32 (*AddrCompare) (struct qsockaddr *addr1, struct qsockaddr *addr2);
	s32 (*GetSocketPort) (struct qsockaddr *addr);
	s32 (*SetSocketPort) (struct qsockaddr *addr, s32 port);
} net_landriver_t;
typedef struct {
	const s8 *name;
	bool initialized;
	s32 (*Init) ();
	void (*Listen) (bool state);
	void (*SearchForHosts) (bool xmit);
	qsocket_t *(*Connect) (const s8 *host);
	qsocket_t *(*CheckNewConnections) ();
	s32 (*QGetMessage) (qsocket_t *sock);
	s32 (*QSendMessage) (qsocket_t *sock, sizebuf_t *data);
	s32 (*SendUnreliableMessage) (qsocket_t *sock, sizebuf_t *data);
	bool (*CanSendMessage) (qsocket_t *sock);
	bool (*CanSendUnreliableMessage) (qsocket_t *sock);
	void (*Close) (qsocket_t *sock);
	void (*Shutdown) ();
} net_driver_t;
typedef struct {
	s8 name[16];
	s8 map[16];
	s8 cname[32];
	s32 users;
	s32 maxusers;
	s32 driver;
	s32 ldriver;
	struct qsockaddr addr;
} hostcache_t;
typedef struct _PollProcedure {
	struct _PollProcedure *next;
	f64 nextTime;
	void (*procedure)(void *);
	void *arg;
} PollProcedure;

typedef enum hudstyle_t {                                            // screen.h
	HUD_CLASSIC,
	HUD_MODERN_CENTERAMMO, // Modern 1
	HUD_MODERN_SIDEAMMO, // Modern 2
	HUD_QUAKEWORLD,
	HUD_COUNT,
} hudstyle_t;

typedef struct {                                                        // wad.h
	s32 width, height;
	u8 data[]; // variably sized
} qpic_t;
typedef struct {
	s8 identification[4]; // should be WAD2 or 2DAW
	s32 numlumps;
	s32 infotableofs;
} wadinfo_t;
typedef struct {
	s32 filepos;
	s32 disksize;
	s32 size; // uncompressed
	s8 type;
	s8 compression;
	s8 pad1, pad2;
	s8 name[16]; // must be null terminated
} lumpinfo_t;
typedef struct wad_s {
	s8 name[MAX_QPATH];
	s32 id;
	fshandle_t fh;
	s32 numlumps;
	lumpinfo_t *lumps;
	struct wad_s *next;
} wad_t;

typedef struct surfcache_s {                                        // d_local.h
	struct surfcache_s *next;
	struct surfcache_s **owner; // NULL is an empty chunk of memory
	s32 lightadj[MAXLIGHTMAPS]; // checked for strobe flush
	s32 dlight;
	s32 size; // including header
	unsigned width;
	unsigned height; // DEBUG only needed for debug
	f32 mipscale;
	struct texture_s *texture; // checked for animating textures
	u8 data[4]; // width*height elements
} surfcache_t;
typedef struct sspan_s {
	s32 u, v, count;
} sspan_t;

typedef s32 func_t;                                                 // pr_comp.h
typedef s32 string_t;
typedef enum {
	ev_bad = -1, ev_void = 0, ev_string, ev_float,
	ev_vector, ev_entity, ev_field, ev_function,
	ev_pointer
} etype_t;
enum { OP_DONE, OP_MUL_F, OP_MUL_V, OP_MUL_FV, OP_MUL_VF,
	OP_DIV_F, OP_ADD_F, OP_ADD_V, OP_SUB_F, OP_SUB_V,
	OP_EQ_F, OP_EQ_V, OP_EQ_S, OP_EQ_E, OP_EQ_FNC,
	OP_NE_F, OP_NE_V, OP_NE_S, OP_NE_E, OP_NE_FNC,
	OP_LE, OP_GE, OP_LT, OP_GT,
	OP_LOAD_F, OP_LOAD_V, OP_LOAD_S,
	OP_LOAD_ENT, OP_LOAD_FLD, OP_LOAD_FNC,
	OP_ADDRESS,
	OP_STORE_F, OP_STORE_V, OP_STORE_S,
	OP_STORE_ENT, OP_STORE_FLD, OP_STORE_FNC,
	OP_STOREP_F, OP_STOREP_V, OP_STOREP_S,
	OP_STOREP_ENT, OP_STOREP_FLD, OP_STOREP_FNC,
	OP_RETURN, OP_NOT_F, OP_NOT_V, OP_NOT_S, OP_NOT_ENT,
	OP_NOT_FNC, OP_IF, OP_IFNOT,
	OP_CALL0, OP_CALL1, OP_CALL2, OP_CALL3, OP_CALL4,
	OP_CALL5, OP_CALL6, OP_CALL7, OP_CALL8,
	OP_STATE, OP_GOTO, OP_AND, OP_OR,
	OP_BITAND, OP_BITOR
};
typedef struct statement_s {
	u16 op;
	s16 a, b, c;
} dstatement_t;
typedef struct {
	u16 type;
	u16 ofs;
	s32 s_name;
} ddef_t;
typedef struct {
	s32 first_statement; // negative numbers are builtins
	s32 parm_start;
	s32 locals; // total ints of parms + locals
	s32 profile; // runtime
	s32 s_name;
	s32 s_file; // source file defined in
	s32 numparms;
	u8 parm_size[MAX_PARMS];
} dfunction_t;
typedef struct {
	s32 version;
	s32 crc; // check of header file
	s32 ofs_statements;
	s32 numstatements; // statement 0 is an error
	s32 ofs_globaldefs;
	s32 numglobaldefs;
	s32 ofs_fielddefs;
	s32 numfielddefs;
	s32 ofs_functions;
	s32 numfunctions; // function 0 is an empty
	s32 ofs_strings;
	s32 numstrings; // first string is a null string
	s32 ofs_globals;
	s32 numglobals;
	s32 entityfields;
} dprograms_t;

typedef struct {                                                   // progdefs.h
	s32 pad[28];
	s32 self;
	s32 other;
	s32 world;
	f32 time;
	f32 frametime;
	f32 force_retouch;
	string_t mapname;
	f32 deathmatch;
	f32 coop;
	f32 teamplay;
	f32 serverflags;
	f32 total_secrets;
	f32 total_monsters;
	f32 found_secrets;
	f32 killed_monsters;
	f32 parm1;
	f32 parm2;
	f32 parm3;
	f32 parm4;
	f32 parm5;
	f32 parm6;
	f32 parm7;
	f32 parm8;
	f32 parm9;
	f32 parm10;
	f32 parm11;
	f32 parm12;
	f32 parm13;
	f32 parm14;
	f32 parm15;
	f32 parm16;
	vec3_t v_forward;
	vec3_t v_up;
	vec3_t v_right;
	f32 trace_allsolid;
	f32 trace_startsolid;
	f32 trace_fraction;
	vec3_t trace_endpos;
	vec3_t trace_plane_normal;
	f32 trace_plane_dist;
	s32 trace_ent;
	f32 trace_inopen;
	f32 trace_inwater;
	s32 msg_entity;
	func_t main;
	func_t StartFrame;
	func_t PlayerPreThink;
	func_t PlayerPostThink;
	func_t ClientKill;
	func_t ClientConnect;
	func_t PutClientInServer;
	func_t ClientDisconnect;
	func_t SetNewParms;
	func_t SetChangeParms;
} globalvars_t;
typedef struct {
	f32 modelindex;
	vec3_t absmin;
	vec3_t absmax;
	f32 ltime;
	f32 movetype;
	f32 solid;
	vec3_t origin;
	vec3_t oldorigin;
	vec3_t velocity;
	vec3_t angles;
	vec3_t avelocity;
	vec3_t punchangle;
	string_t classname;
	string_t model;
	f32 frame;
	f32 skin;
	f32 effects;
	vec3_t mins;
	vec3_t maxs;
	vec3_t size;
	func_t touch;
	func_t use;
	func_t think;
	func_t blocked;
	f32 nextthink;
	s32 groundentity;
	f32 health;
	f32 frags;
	f32 weapon;
	string_t weaponmodel;
	f32 weaponframe;
	f32 currentammo;
	f32 ammo_shells;
	f32 ammo_nails;
	f32 ammo_rockets;
	f32 ammo_cells;
	f32 items;
	f32 takedamage;
	s32 chain;
	f32 deadflag;
	vec3_t view_ofs;
	f32 button0;
	f32 button1;
	f32 button2;
	f32 impulse;
	f32 fixangle;
	vec3_t v_angle;
	f32 idealpitch;
	string_t netname;
	s32 enemy;
	f32 flags;
	f32 colormap;
	f32 team;
	f32 max_health;
	f32 teleport_time;
	f32 armortype;
	f32 armorvalue;
	f32 waterlevel;
	f32 watertype;
	f32 ideal_yaw;
	f32 yaw_speed;
	s32 aiment;
	s32 goalentity;
	f32 spawnflags;
	string_t target;
	string_t targetname;
	f32 dmg_take;
	f32 dmg_save;
	s32 dmg_inflictor;
	s32 owner;
	vec3_t movedir;
	string_t message;
	f32 sounds;
	string_t noise;
	string_t noise1;
	string_t noise2;
	string_t noise3;
} entvars_t;

typedef union eval_s {                                                // progs.h
	string_t string;
	f32 _float;
	f32 vector[3];
	func_t function;
	s32 _int;
	s32 edict;
} eval_t;
typedef struct edict_s {
	bool free;
	link_t area;
	s32 num_leafs;
	s32 leafnums[MAX_ENT_LEAFS];
	entity_state_t baseline;
	u8 alpha;
	u8 scale;
	bool sendinterval;
	f32 oldframe;
	f32 oldthinktime;
	f32 freetime;
	entvars_t v;
} edict_t;
typedef struct {
	const s8 *name;
	s32 first_statement;
	s32 patch_statement;
} exbuiltin_t;
typedef void (*builtin_t) ();

typedef struct {                                                      // world.h
	vec3_t normal;
	f32 dist;
} plane_t;
typedef struct {
	bool allsolid; // if 1, plane is not valid
	bool startsolid; // if 1, the initial point was in a solid area
	bool inopen, inwater;
	f32 fraction; // time completed, 1.0 = didn't hit anything
	vec3_t endpos; // final position
	plane_t plane; // surface normal at impact
	edict_t *ent; // entity the surface is on
} trace_t;

typedef void (*xcommand_t) ();                                          // cmd.h
typedef struct cmdalias_s {
	struct cmdalias_s *next;
	s8 name[MAX_ALIAS_NAME];
	s8 *value;
} cmdalias_t;
typedef struct cmd_function_s {
	struct cmd_function_s *next;
	s8 *name;
	xcommand_t function;
} cmd_function_t;
typedef enum { src_client, src_command } cmd_source_t;

typedef struct { s32 left; s32 right; } portable_samplepair_t;      // q_sound.h
typedef struct sfx_s { s8 name[MAX_QPATH]; cache_user_t cache; } sfx_t;
typedef struct {
	s32 length;
	s32 loopstart;
	s32 speed;
	s32 width;
	s32 stereo;
	u8 data[1]; /* variable sized */
} sfxcache_t;
typedef struct {
	s32 channels;
	s32 samples; /* mono samples in buffer */
	s32 submission_chunk; /* don't mix less than this # */
	s32 samplepos; /* in mono samples */
	s32 samplebits;
	s32 signed8; /* device opened for S8 format? (e.g. Amiga AHI) */
	s32 speed;
	u8 *buffer;
} dma_t;
typedef struct {
	sfx_t *sfx; /* sfx number */
	s32 leftvol; /* 0-255 volume */
	s32 rightvol; /* 0-255 volume */
	s32 end; /* end time in global paintsamples */
	s32 pos; /* sample position in sfx */
	s32 looping; /* where to loop, -1 = no looping */
	s32 entnum; /* to allow overriding a specific sound */
	s32 entchannel;
	vec3_t origin; /* origin of sound effect */
	vec_t dist_mult; /* distance multiplier (attenuation/clipK) */
	s32 master_vol; /* 0-255 master volume */
} channel_t;
typedef struct {
	s32 rate;
	s32 width;
	s32 channels;
	s32 loopstart;
	s32 samples;
	s32 dataofs; /* chunk starts this many bytes from file start */
} wavinfo_t;

typedef struct {                                                     // client.h
	s32 length;
	s8 map[MAX_STYLESTRING];
	s8 average; //johnfitz
	s8 peak; //johnfitz
} lightstyle_t;
typedef struct {
	s8 name[MAX_SCOREBOARDNAME];
	f32 entertime;
	s32 frags;
	s32 colors; // two 4 bit fields
	u8 translations[VID_GRADES*256];
} scoreboard_t;
typedef struct {
	s32 destcolor[3];
	s32 percent; // 0-256
} cshift_t;
typedef struct {
	vec3_t origin;
	f32 radius;
	f32 die; // stop lighting after this time
	f32 decay; // drop this each second
	f32 minlight; // don't add when contributing less
	s32 key;
	vec3_t color; //johnfitz -- lit support via lordhavoc
} dlight_t;
typedef struct {
	s32 entity;
	struct model_s *model;
	f32 endtime;
	vec3_t start, end;
} beam_t;
typedef enum {
	ca_dedicated, // dedicated server with no ability to start a client
	ca_disconnected, // full screen console with no connection
	ca_connected // valid netcon, talking to a server
} cactive_t;
typedef struct {
	cactive_t state; // personalization data sent to server
	s8 spawnparms[MAX_MAPSTRING]; // to restart a level
	s32 demonum; // -1 = don't play demos
	s8 demos[MAX_DEMOS][MAX_DEMONAME];
	bool demorecording;
	bool demoplayback;
	bool demopaused;
	bool timedemo;
	s32 forcetrack; // -1 = use normal cd track
	FILE *demofile;
	s32 td_lastframe; // to meter out one message a frame
	s32 td_startframe; // host_framecount at start
	f32 td_starttime; // realtime at second frame of timedemo
	s32 signon; // 0 to SIGNONS
	struct qsocket_s *netcon;
	sizebuf_t message; // writing buffer to send to server
} client_static_t;
extern client_static_t cls;
typedef struct {
	s32 movemessages;
	usercmd_t cmd; // last command sent to the server
	usercmd_t pendingcmd;
	s32 stats[MAX_CL_STATS]; // health, etc
	s32 items; // inventory bit flags
	f32 item_gettime[32]; // cl.time of aquiring item, for blinking
	f32 faceanimtime; // use anim frame if cl.time < this
	cshift_t cshifts[NUM_CSHIFTS]; // color shifts for damage, powerups
	cshift_t prev_cshifts[NUM_CSHIFTS];
	vec3_t mviewangles[2];
	vec3_t viewangles;
	vec3_t mvelocity[2];
	vec3_t velocity; // lerped between mvelocity[0] and [1]
	vec3_t punchangle; // temporary offset
	f32 idealpitch;
	f32 pitchvel;
	bool nodrift;
	f32 driftmove;
	f64 laststop;
	f32 viewheight;
	f32 crouch; // local amount for smoothing stepups
	bool paused; // send over by server
	bool onground;
	bool inwater;
	s32 intermission; // don't change view angle, full screen, etc
	s32 completed_time; // latched at intermission start
	f64 mtime[2]; // the timestamp of last two messages
	f64 time;
	f64 oldtime;
	f32 last_received_message;
	struct model_s *model_precache[MAX_MODELS];
	struct sfx_s *sound_precache[MAX_SOUNDS];
	s8 mapname[128];
	s8 levelname[128];
	s32 viewentity; // cl_entitites[cl.viewentity] = player
	s32 maxclients;
	s32 gametype;
	struct model_s *worldmodel; // cl_entitites[0].model
	struct efrag_s *free_efrags;
	s32 num_efrags;
	s32 num_entities; // held in cl_entities array
	s32 num_statics; // held in cl_staticentities array
	entity_t viewent; // the gun model
	s32 cdtrack, looptrack; // cd audio
	scoreboard_t *scores; // [cl.maxclients]
	unsigned protocol; //johnfitz
	unsigned protocolflags;
} client_state_t;
typedef struct {
	s32 down[2]; // key nums holding it down
	s32 state; // low bit is down state
} kbutton_t;

typedef struct {                                                     // client.h
	s32 maxclients;
	s32 maxclientslimit;
	struct client_s *clients; // [maxclients]
	s32 serverflags; // episode completion information
	bool changelevel_issued; // cleared when at SV_SpawnServer
} server_static_t;
typedef enum {ss_loading, ss_active} server_state_t;
typedef struct {
	bool active; // 0 if only a net client
	bool paused;
	bool loadgame; // handle connections specially
	bool nomonsters; // server started with 'nomonsters' cvar active
	f64 time;
	s32 lastcheck; // used by PF_checkclient
	f64 lastchecktime;
	s8 name[64]; // map name
	s8 modelname[64]; // maps/<name>.bsp, for model_precache[0]
	struct model_s *worldmodel;
	const s8 *model_precache[MAX_MODELS]; // NULL terminated
	struct model_s *models[MAX_MODELS];
	const s8 *sound_precache[MAX_SOUNDS]; // NULL terminated
	const s8 *lightstyles[MAX_LIGHTSTYLES];
	s32 num_edicts;
	s32 max_edicts;
	edict_t *edicts;
	server_state_t state; // some actions are only valid during load
	sizebuf_t datagram;
	u8 datagram_buf[MAX_DATAGRAM];
	sizebuf_t reliable_datagram;
	u8 reliable_datagram_buf[MAX_DATAGRAM];
	sizebuf_t *signon;
	s32 num_signon_buffers;
	sizebuf_t *signon_buffers[MAX_SIGNON_BUFFERS];
	unsigned protocol; //johnfitz
	unsigned protocolflags;
} server_t;
enum sendsignon_e {
	PRESPAWN_DONE,
	PRESPAWN_FLUSH=1,
	PRESPAWN_SIGNONBUFS,
	PRESPAWN_SIGNONMSG,
};
typedef struct client_s {
	bool active; // 0 = client is free
	bool spawned; // 0 = don't send datagrams
	bool dropasap; // has been told to go to another level
	enum sendsignon_e sendsignon; // only valid before spawned
	s32 signonidx;
	f64 last_message; // reliable messages must be sent periodically
	struct qsocket_s *netconnection; // communications handle
	usercmd_t cmd; // movement
	vec3_t wishdir; // intended motion calced from cmd
	sizebuf_t message;
	u8 msgbuf[MAX_MSGLEN];
	edict_t *edict; // EDICT_NUM(clientnum+1)
	s8 name[32]; // for printing to other people
	s32 colors;
	f32 ping_times[NUM_PING_TIMES];
	s32 num_pings; // ping_times[num_pings%NUM_PING_TIMES]
	f32 spawn_parms[NUM_SPAWN_PARMS];
	s32 old_frags;
} client_t;

typedef struct {                                                       // draw.c
	vrect_t rect;
	s32 width;
	s32 height;
	u8 *ptexbytes;
	s32 rowbytes;
} rectdesc_t;
typedef struct cachepic_s {
	s8 name[MAX_QPATH];
	cache_user_t cache;
} cachepic_t;

typedef struct {                                                     // r_draw.c
	f32 u, v;
	s32 ceilv;
} evert_t;

typedef enum { touchessolid, drawnode, nodrawnode } solidstate_t;     // r_bsp.c

typedef struct { s32 index0; s32 index1; } aedge_t;                 // r_alias.c

typedef struct vispatch_s { // External VIS file support              // model.c
	s8 mapname[32];
	s32 filelen; // length of data after header (VIS+Leafs)
} vispatch_t;

typedef struct { s32 s; dfunction_t *f; } prstack_t;                // pr_exec.c

typedef struct { s8 *name; s8 *description; } level_t;                 // menu.c
typedef struct { s8 *description; s32 firstLevel; s32 levels; } episode_t;

typedef struct {                                                     // screen.c
	s8 manufacturer;
	s8 version;
	s8 encoding;
	s8 bits_per_pixel;
	u16 xmin, ymin, xmax, ymax;
	u16 hres, vres;
	u8 palette[48];
	s8 reserved;
	s8 color_planes;
	u16 bytes_per_line;
	u16 palette_type;
	s8 filler[58];
	u8 data; // unbounded
} pcx_t;

typedef struct{ddef_t *pcache;s8 field[MAX_FIELD_LEN];}gefv_cache; // pr_edict.c

typedef struct {                                                   // d_polyse.c
	s32 quotient;
	s32 remainder;
} adivtab_t;
typedef struct {
	void *pdest;
	s16 *pz;
	s32 count;
	u8 *ptex;
	s32 sfrac, tfrac, light, zi;
} spanpackage_t;
typedef struct {
	s32 isflattop;
	s32 numleftedges;
	s32 *pleftedgevert0;
	s32 *pleftedgevert1;
	s32 *pleftedgevert2;
	s32 numrightedges;
	s32 *prightedgevert0;
	s32 *prightedgevert1;
	s32 *prightedgevert2;
} edgetable;

typedef struct cache_system_s {                                        // zone.c
	s32 size; // including this header
	cache_user_t *user;
	s8 name[CACHENAME_LEN];
	struct cache_system_s *prev, *next;
	struct cache_system_s *lru_prev, *lru_next; // for LRU flushing
} cache_system_t;
typedef struct memblock_s {
	s32 size; // including the header and possibly tiny fragments
	s32 tag; // a tag of 0 is a free block
	s32 id; // should be ZONEID
	s32 pad; // pad to 64 bit boundary
	struct memblock_s *next, *prev;
} memblock_t;
typedef struct {
	s32 size; // total bytes malloced, including header
	memblock_t blocklist; // start / end cap for linked list
	memblock_t *rover;
} memzone_t;
typedef struct {
	s32 sentinel;
	s32 size; // including sizeof(hunk_t), -1 = not allocated
	s8 name[HUNKNAME_LEN];
} hunk_t;

typedef struct { // on-disk pakfile                                  // common.c
	s8 name[56];
	s32 filepos, filelen;
} dpackfile_t;
typedef struct {
	s8 id[4];
	s32 dirofs;
	s32 dirlen;
} dpackheader_t;
typedef struct {
	s8 *key;
	s8 *value;
} locentry_t;
typedef struct {
	s32 numentries;
	s32 maxnumentries;
	s32 numindices;
	unsigned *indices;
	locentry_t *entries;
	s8 *text;
} localization_t;

typedef struct stdio_buffer_s {                                       // image.c
	FILE *f;
	u8 buffer[1024];
	s32 size;
	s32 pos;
} stdio_buffer_t;
typedef struct targaheader_s {
	u8 id_length, colormap_type, image_type;
	u16 colormap_index, colormap_length;
	u8 colormap_size;
	u16 x_origin, y_origin, width, height;
	u8 pixel_size, attributes;
} targaheader_t;

typedef struct {                                                      // world.c
	vec3_t boxmins, boxmaxs; // enclose the test object along entire move
	f32 *mins, *maxs; // size of the moving object
	vec3_t mins2, maxs2; // size when clipping against mosnters
	f32 *start, *end;
	trace_t trace;
	s32 type;
	edict_t *passedict;
} moveclip_t;
typedef struct areanode_s {
	s32 axis; // -1 = leaf node
	f32 dist;
	struct areanode_s *children[2];
	link_t trigger_edicts;
	link_t solid_edicts;
} areanode_t;

typedef struct {                                                    // snd_mix.c
	f32 *memory; // kernelsize floats
	f32 *kernel; // kernelsize floats
	s32 kernelsize; // M+1, rounded up to be a multiple of 16
	s32 M; // M value used to make kernel, even
	s32 parity; // 0-3
	f32 f_c; // cutoff frequency, [0..1], fraction of sample rate
} filter_t;
#endif
