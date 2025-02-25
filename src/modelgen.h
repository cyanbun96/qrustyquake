// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.

// modelgen.h: header file for model generation program

// *********************************************************
// * This file must be identical in the modelgen directory *
// * and in the Quake directory, because it's used to *
// * pass data from one to the other via model files. *
// *********************************************************

#define ALIAS_VERSION 6
#define ALIAS_ONSEAM 0x0020
#define DT_FACES_FRONT 0x0010
#define IDPOLYHEADER (('O'<<24)+('P'<<16)+('D'<<8)+'I') // little-endian "IDPO"
#ifndef SYNCTYPE_T
#define SYNCTYPE_T // must match definition in spritegn.h
typedef enum { ST_SYNC=0, ST_RAND } synctype_t;
#endif
typedef enum { ALIAS_SINGLE=0, ALIAS_GROUP } aliasframetype_t;
typedef enum { ALIAS_SKIN_SINGLE=0, ALIAS_SKIN_GROUP } aliasskintype_t;

typedef struct {
	int ident;
	int version;
	vec3_t scale;
	vec3_t scale_origin;
	float boundingradius;
	vec3_t eyeposition;
	int numskins;
	int skinwidth;
	int skinheight;
	int numverts;
	int numtris;
	int numframes;
	synctype_t synctype;
	int flags;
	float size;
} mdl_t;

typedef struct {
	int onseam;
	int s;
	int t;
} stvert_t;

typedef struct dtriangle_s {
	int facesfront;
	int vertindex[3];
} dtriangle_t;

typedef struct { // This mirrors trivert_t in trilib.h, is present so Quake
	byte v[3]; // knows how to load this data
	byte lightnormalindex;
} trivertx_t;

typedef struct {
	trivertx_t bboxmin; // lightnormal isn't used
	trivertx_t bboxmax; // lightnormal isn't used
	char name[16]; // frame name from grabbing
} daliasframe_t;

typedef struct {
	int numframes;
	trivertx_t bboxmin; // lightnormal isn't used
	trivertx_t bboxmax; // lightnormal isn't used
} daliasgroup_t;

typedef struct {
	int numskins;
} daliasskingroup_t;

typedef struct {
	float interval;
} daliasinterval_t;

typedef struct {
	float interval;
} daliasskininterval_t;

typedef struct {
	aliasframetype_t type;
} daliasframetype_t;

typedef struct {
	aliasskintype_t type;
} daliasskintype_t;
