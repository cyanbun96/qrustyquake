// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.

#include "modelgen.h"
#include "spritegn.h"

#define SIDE_FRONT 0
#define SIDE_BACK 1
#define SIDE_ON 2
#define SURF_PLANEBACK 2
#define SURF_DRAWSKY 4
#define SURF_DRAWSPRITE 8
#define SURF_DRAWTURB 0x10
#define SURF_DRAWTILED 0x20
#define SURF_DRAWBACKGROUND 0x40
#define SURF_DRAWCUTOUT 0x80
#define EF_ROCKET 1 // leave a trail
#define EF_GRENADE 2 // leave a trail
#define EF_GIB 4 // leave a trail
#define EF_ROTATE 8 // rotate(bonus items)
#define EF_TRACER 16 // green split trail
#define EF_ZOMGIB 32 // small blood trail
#define EF_TRACER2 64 // orange split trail + rotate
#define EF_TRACER3 128 // purple trail

// d*_t structures are on-disk representations
// m*_t structures are in-memory

typedef struct
{
	vec3_t position;
} mvertex_t;

typedef struct mplane_s
{
	vec3_t normal;
	float dist;
	byte type; // for texture axis selection and fast side tests
	byte signbits; // signx + signy<<1 + signz<<1
	byte pad[2];
} mplane_t;

typedef struct texture_s
{
	char name[16];
	unsigned width, height;
	int anim_total; // total tenths in sequence( 0 = no)
	int anim_min, anim_max; // time for this frame min <=time< max
	struct texture_s *anim_next; // in the animation sequence
	struct texture_s *alternate_anims; // bmodels in frmae 1 use these
	unsigned offsets[MIPLEVELS]; // four mip maps stored
} texture_t;

typedef struct
{
	unsigned short v[2];
	unsigned int cachededgeoffset;
} medge_t;

typedef struct
{
	float vecs[2][4];
	float mipadjust;
	texture_t *texture;
	int flags;
} mtexinfo_t;

typedef struct msurface_s
{
	unsigned int visframe; // should be drawn when node is crossed
	unsigned int dlightframe;
	int dlightbits;
	mplane_t *plane;
	int flags;
	int firstedge; // look up in model->surfedges[], negative numbers
	int numedges; // are backwards edges
	struct surfcache_s *cachespots[MIPLEVELS]; // surface generation data
	short texturemins[2];
	short extents[2];
	mtexinfo_t *texinfo;
	byte styles[MAXLIGHTMAPS]; // lighting info
	byte *samples; // [numstyles*surfsize]
} msurface_t;

typedef struct mnode_s
{
	// common with leaf
	int contents; // 0, to differentiate from leafs
	unsigned int visframe; // node needs to be traversed if current
	short minmaxs[6]; // for bounding box culling
	struct mnode_s *parent;
	// node specific
	mplane_t *plane;
	struct mnode_s *children[2]; 
	unsigned short firstsurface;
	unsigned short numsurfaces;
} mnode_t;

typedef struct mleaf_s
{
	// common with node
	int contents; // wil be a negative contents number
	unsigned int visframe; // node needs to be traversed if current
	short minmaxs[6]; // for bounding box culling
	struct mnode_s *parent;
	// leaf specific
	byte *compressed_vis;
	efrag_t *efrags;
	msurface_t **firstmarksurface;
	int nummarksurfaces;
	int key; // BSP sequence number for leaf's contents
	byte ambient_sound_level[NUM_AMBIENTS];
} mleaf_t;

typedef struct
{
	dclipnode_t *clipnodes;
	mplane_t *planes;
	int firstclipnode;
	int lastclipnode;
	vec3_t clip_mins;
	vec3_t clip_maxs;
} hull_t;

typedef struct mspriteframe_s // FIXME: shorten these?
{
	int width;
	int height;
	void *pcachespot; // remove?
	float up, down, left, right;
	byte pixels[4];
} mspriteframe_t;

typedef struct
{
	int numframes;
	float *intervals;
	mspriteframe_t *frames[1];
} mspritegroup_t;

typedef struct
{
	spriteframetype_t type;
	mspriteframe_t *frameptr;
} mspriteframedesc_t;

typedef struct
{
	int type;
	int maxwidth;
	int maxheight;
	int numframes;
	float beamlength; // remove?
	void *cachespot; // remove?
	mspriteframedesc_t frames[1];
} msprite_t;

typedef struct
{ // Alias models are position independent, so the cache manager can move them.
	aliasframetype_t type;
	trivertx_t bboxmin;
	trivertx_t bboxmax;
	int frame;
	char name[16];
} maliasframedesc_t;

typedef struct
{
	aliasskintype_t type;
	void *pcachespot;
	int skin;
} maliasskindesc_t;

typedef struct
{
	trivertx_t bboxmin;
	trivertx_t bboxmax;
	int frame;
} maliasgroupframedesc_t;

typedef struct
{
	int numframes;
	int intervals;
	maliasgroupframedesc_t frames[1];
} maliasgroup_t;

typedef struct
{
	int numskins;
	int intervals;
	maliasskindesc_t skindescs[1];
} maliasskingroup_t;

typedef struct mtriangle_s {
	int facesfront;
	int vertindex[3];
} mtriangle_t;

typedef struct {
	int model;
	int stverts;
	int skindesc;
	int triangles;
	maliasframedesc_t frames[1];
} aliashdr_t;

typedef enum {mod_brush, mod_sprite, mod_alias} modtype_t;

typedef struct model_s
{
	char name[MAX_QPATH];
	qboolean needload; // bmodels and sprites don't cache normally
	modtype_t type;
	int numframes;
	synctype_t synctype;
	int flags;
	vec3_t mins, maxs; // volume occupied by the model
	float radius;
	int firstmodelsurface, nummodelsurfaces; // brush model
	int numsubmodels;
	dmodel_t *submodels;
	int numplanes;
	mplane_t *planes;
	int numleafs; // number of visible leafs, not counting 0
	mleaf_t *leafs;
	int numvertexes;
	mvertex_t *vertexes;
	int numedges;
	medge_t *edges;
	int numnodes;
	mnode_t *nodes;
	int numtexinfo;
	mtexinfo_t *texinfo;
	int numsurfaces;
	msurface_t *surfaces;
	int numsurfedges;
	int *surfedges;
	int numclipnodes;
	dclipnode_t *clipnodes;
	int nummarksurfaces;
	msurface_t **marksurfaces;
	hull_t hulls[MAX_MAP_HULLS];
	int numtextures;
	texture_t **textures;
	byte *visdata;
	byte *lightdata;
	char *entities;
	cache_user_t cache; // only access through Mod_Extradata
} model_t;

void Mod_Init();
void Mod_ClearAll();
model_t *Mod_ForName(char *name, qboolean crash);
void *Mod_Extradata(model_t *mod); // handles caching
void Mod_TouchModel(char *name);
mleaf_t *Mod_PointInLeaf(float *p, model_t *model);
byte *Mod_LeafPVS(mleaf_t *leaf, model_t *model);
