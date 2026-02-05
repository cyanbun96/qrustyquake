// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// models are the only shared resource between a client and server running
// on the same machine.
#include "quakedef.h"

static aliashdr_t *pheader;
// pose is a single set of vertexes. frame may be an animating sequence of poses
static trivertx_t *poseverts[MAXALIASFRAMES];
static s32 posenum;
static model_t* loadmodel;
static s8 loadname[32]; // for hunk tags
static u8 *mod_base;
static void Mod_LoadSpriteModel(model_t *mod, void *buffer);
static void Mod_LoadBrushModel(model_t *mod, void *buffer);
static void Mod_LoadAliasModel(model_t *mod, void *buffer);
static model_t *Mod_LoadModel(model_t *mod, bool crash);
static void Mod_Print();
static u8 *mod_novis;
static s32 mod_novis_capacity;
static u8 *mod_decompressed;
static s32 mod_decompressed_capacity;
static model_t mod_known[MAX_MOD_KNOWN];
static s32 mod_numknown;
static texture_t *r_notexture_mip2; //non-lightmapped surfs without a texture
extern texture_t *r_notexture_mip;

u8 *Mod_NoVisPVS(model_t *model);

void Mod_Init()
{
	Cvar_RegisterVariable(&external_vis);
	Cvar_RegisterVariable(&external_ents);
	Cvar_RegisterVariable(&external_textures);
	Cmd_AddCommand("mcache", Mod_Print);
	r_notexture_mip = Hunk_AllocName(sizeof(texture_t), "r_notexture_mip");
	strcpy(r_notexture_mip->name, "notexture");
	r_notexture_mip->height = r_notexture_mip->width = 32;
	r_notexture_mip2 = Hunk_AllocName(sizeof(texture_t),"r_notexture_mip2");
	strcpy(r_notexture_mip2->name, "notexture2");
	r_notexture_mip2->height = r_notexture_mip2->width = 32;
}

void *Mod_Extradata(model_t *mod)
{ // Caches the data if needed
	void *r = Cache_Check(&mod->cache);
	if(r) return r;
	Mod_LoadModel(mod, 1);
	if(!mod->cache.data) Sys_Error("Mod_Extradata: caching failed");
	return mod->cache.data;
}

mleaf_t *Mod_PointInLeaf(vec3_t p, model_t *model)
{
	if(!model || !model->nodes) Sys_Error("Mod_PointInLeaf: bad model");
	mnode_t *node = model->nodes;
	while(1) {
		if(node->contents < 0) return(mleaf_t *)node;
		mplane_t *plane = node->plane;
		f32 d = DotProduct(p,plane->normal) - plane->dist;
		node = node->children[d > 0 ? 0 : 1];
	}
	return NULL; // never reached
}

static u8 *Mod_DecompressVis(u8 *in, model_t *model)
{
	s32 row = (model->numleafs+7)>>3;
	if(mod_decompressed == NULL || row > mod_decompressed_capacity){
		mod_decompressed_capacity = row;
		mod_decompressed = (u8 *) realloc(mod_decompressed,
				mod_decompressed_capacity);
		if(!mod_decompressed)
		    Sys_Error("Mod_DecompressVis: realloc() failed on %d bytes",
				    mod_decompressed_capacity);
	}
	u8 *out = mod_decompressed;
	u8 *outend = mod_decompressed + row;
	if(!in){ // no vis info, so make all visible
		while(row){
			*out++ = 0xff;
			row--;
		}
		return mod_decompressed;
	} do {
		if(*in){
			*out++ = *in++;
			continue;
		}
		s32 c = in[1];
		in += 2;
		if (c > row - (out - mod_decompressed))
			c = row - (out - mod_decompressed);	//now that we're dynamically allocating pvs buffers, we have to be more careful to avoid heap overflows with buggy maps.
		while(c){
			if(out == outend){
				if(!model->viswarn){
					model->viswarn = 1;
     printf("Mod_DecompressVis: output overrun on model \"%s\"\n", model->name);
				}
				return mod_decompressed;
			}
			*out++ = 0;
			c--;
		}
	} while(out - mod_decompressed < row);
	return mod_decompressed;
}

u8 *Mod_LeafPVS(mleaf_t *leaf, model_t *model)
{
	if(leaf == model->leafs) return Mod_NoVisPVS(model);
	return Mod_DecompressVis(leaf->compressed_vis, model);
}


u8 *Mod_NoVisPVS(model_t *model)
{
	s32 pvsbytes = (model->numleafs+7)>>3;
	if(mod_novis == NULL || pvsbytes > mod_novis_capacity){
		mod_novis_capacity = pvsbytes;
		mod_novis = (u8 *) realloc(mod_novis, mod_novis_capacity);
		if(!mod_novis)
    Sys_Error("Mod_NoVisPVS: realloc() failed on %d bytes", mod_novis_capacity);
		memset(mod_novis, 0xff, mod_novis_capacity);
	}
	return mod_novis;
}

void Mod_ClearAll()
{
	s32 i = 0;
	model_t *mod = mod_known;
	for(; i < mod_numknown; i++, mod++){
		mod->needload = NL_UNREFERENCED;
		if(mod->type == mod_sprite) mod->cache.data = NULL;
	}
}

void Mod_ResetAll()
{
	s32 i = 0;
	model_t *mod = mod_known;
	for(; i < mod_numknown; i++, mod++)
		memset(mod, 0, sizeof(model_t));
	mod_numknown = 0;
}

static model_t *Mod_FindName(const s8 *name)
{
	if(!name[0]) Sys_Error("Mod_FindName: NULL name");
	s32 i = 0; // search the currently loaded models
	model_t *mod = mod_known;
	for(; i<mod_numknown ; i++, mod++)
		if(!strcmp(mod->name, name)) break;
	if(i == mod_numknown){
		if(mod_numknown == MAX_MOD_KNOWN)
			Sys_Error("mod_numknown == MAX_MOD_KNOWN");
		q_strlcpy(mod->name, name, MAX_QPATH);
		mod->needload = 1;
		mod_numknown++;
	}
	return mod;
}

void Mod_TouchModel(const s8 *name)
{
	model_t *mod = Mod_FindName(name);
	if(!mod->needload == NL_PRESENT && mod->type == mod_alias)
		Cache_Check(&mod->cache);
}


static model_t *Mod_LoadModel(model_t *mod, bool crash)
{ // Loads a model into the cache
	u8 stackbuf[1024]; // avoid dirtying the cache heap
	if(mod->type == mod_alias){
		if(Cache_Check(&mod->cache)){
			mod->needload = NL_PRESENT;
			return mod;
		}
	}
	else if(mod->needload == NL_PRESENT) return mod;
	// because the world is so huge, load it one piece at a time
	u8 *buf = 
	COM_LoadStackFile(mod->name, stackbuf, sizeof(stackbuf), &mod->path_id);
	if(!buf){
		if(crash) Host_Error("Mod_LoadModel: %s not found", mod->name);
		return NULL;
	}
	COM_FileBase(mod->name,loadname,sizeof(loadname));//allocate a new model
	loadmodel = mod;
	mod->needload = 0; // fill it in, call the apropriate loader
	s32 mod_type = (buf[0] | (buf[1]<<8) | (buf[2]<<16) | (buf[3]<<24));
	switch(mod_type){
		case IDPOLYHEADER: Mod_LoadAliasModel(mod, buf); break;
		case IDSPRITEHEADER: Mod_LoadSpriteModel(mod, buf); break;
		default: Mod_LoadBrushModel(mod, buf); break;
	}
	return mod;
}

model_t *Mod_ForName(const s8 *name, bool crash)
{ // Loads in a model for the given name
	model_t *mod = Mod_FindName(name);
	return Mod_LoadModel(mod, crash);
}

static u8 palette_search(u8 r, u8 g, u8 b)
{ // qremip by erysdren (it/its), originally cc0 or public domain
	s32 pen, dist = 9999999;
	for (s32 i = 256 - 32; i--;) { // CyanBun96: don't check fullbrights
		s32 rr = host_basepal[i * 3] - r;
		s32 gg = host_basepal[i * 3 + 1] - g;
		s32 bb = host_basepal[i * 3 + 2] - b;
		s32 rank = rr*rr + gg*gg + bb*bb;
		if (rank < dist) {
			pen = i;
			dist = rank;
		}
	}
	return pen;
}

static u8 blend4_u8(f32 alpha, f32 beta, u8 a, u8 b, u8 c, u8 d)
{ // qremip by erysdren (it/its), originally cc0 or public domain
	u8 *a_rgb, *b_rgb, *c_rgb, *d_rgb;
	s32 result[3];
	a_rgb = &host_basepal[a * 3];
	b_rgb = &host_basepal[b * 3];
	c_rgb = &host_basepal[c * 3];
	d_rgb = &host_basepal[d * 3];
	for (s32 i = 0; i < 3; i++) {
		result[i] = (1.0 - alpha) * (1.0 - beta) * a_rgb[i] +
				    alpha * (1.0 - beta) * b_rgb[i] +
				    (1.0 - alpha) * beta * c_rgb[i] +
					    alpha * beta * d_rgb[i];
	}
	return palette_search(result[0], result[1], result[2]);
}

static void bilinear_u8(u8 *in, u8 *out, s32 w, s32 h, s32 new_w, s32 new_h)
{ // qremip by erysdren (it/its), originally cc0 or public domain
	f32 xscale = (f32)w / (f32)new_w; // get x and y scale
	f32 yscale = (f32)h / (f32)new_h;
	for (s32 y = 0; y < new_h; y++) { // loop over all pixels in the new mip
	for (s32 x = 0; x < new_w; x++) {
		f32 src_x = (f32)x * xscale; // get weights
		f32 src_y = (f32)y * yscale;
		s32 xlow = (s32)src_x; // get high and low pixels
		s32 ylow = (s32)src_y;
		s32 xhigh = (xlow + 1 < w) ? xlow + 1 : w - 1;
		s32 yhigh = (ylow + 1 < h) ? ylow + 1 : h - 1;
		u8 a = in[ylow * w + xlow]; // high and low pixels values
		u8 b = in[ylow * w + xhigh];
		u8 c = in[yhigh * w + xlow];
		u8 d = in[yhigh * w + xhigh];
		if (a==0xFF || b==0xFF || c==0xFF || d==0xFF) {
			out[y * new_w + x] = 0xFF; // don't blend transparents
		} else if (a>=0xE0 || b>=0xE0 || c>=0xE0 || d>=0xE0) {
			out[y * new_w + x] = a; // don't blend fullbrights
		} else { // blend color
			f32 alpha = src_x - (f32)xlow;
			f32 beta = src_y - (f32)ylow;
			out[y * new_w + x] = blend4_u8(alpha, beta, a, b, c, d);
		}
	}
	}
}

void Mod_LoadTextures(lump_t *l)
{
	texture_t *anims[10];
	texture_t *altanims[10];
	if(!l->filelen){
		loadmodel->textures = NULL;
		return;
	}
	dmiptexlump_t *m = (dmiptexlump_t *) (mod_base + l->fileofs);
	m->nummiptex = LittleLong(m->nummiptex);
	loadmodel->numtextures = m->nummiptex;
	loadmodel->textures = Hunk_AllocName(m->nummiptex *
			sizeof(*loadmodel->textures), loadname);
	for(s32 i = 0; i < m->nummiptex; i++){
		m->dataofs[i] = LittleLong(m->dataofs[i]);
		if(m->dataofs[i] == -1) continue;
		miptex_t *mt = (miptex_t *) ((u8 *) m + m->dataofs[i]);
		mt->width = LittleLong(mt->width);
		mt->height = LittleLong(mt->height);
		for(s32 j = 0; j < MIPLEVELS; j++)
			mt->offsets[j] = LittleLong(mt->offsets[j]);
		if((mt->width & 15) || (mt->height & 15))
			Sys_Error("Texture %s is not 16 aligned", mt->name);
		s32 pixels = mt->width * mt->height / 64 * 85;
		texture_t *tx=Hunk_AllocName(sizeof(texture_t)+pixels,loadname);
		loadmodel->textures[i] = tx;
		memcpy(tx->name, mt->name, sizeof(tx->name));
		tx->width = mt->width;
		tx->height = mt->height;
		for(s32 j = 0; j < MIPLEVELS; j++)
			tx->offsets[j] = mt->offsets[j] + sizeof(texture_t) -
				sizeof(miptex_t);
		// the pixels immediately follow the structures
		memcpy(tx + 1, mt + 1, pixels);
		if(!Q_strncmp(mt->name, "sky", 3))
			R_InitSky(tx);
		if (!r_rebuildmips.value) continue;
		u8 *base = (u8 *)tx + LittleLong(tx->offsets[0]);
		if (r_rebuildmips.value == 2) { // only cutouts and fullbrights
			for (s32 h = 0; h < tx->height; ++h)
			for (s32 w = 0; w < tx->width; ++w)
				if (base[w + h*tx->width] >= 0xE0) goto rebuild;
			continue;
		}
rebuild:
		u8 *mip1 = (u8 *)tx + LittleLong(tx->offsets[1]);
		u8 *mip2 = (u8 *)tx + LittleLong(tx->offsets[2]);
		u8 *mip3 = (u8 *)tx + LittleLong(tx->offsets[3]);
		bilinear_u8(base, mip1, tx->width, tx->height, tx->width/2, tx->height/2);
		bilinear_u8(base, mip2, tx->width, tx->height, tx->width/4, tx->height/4);
		bilinear_u8(base, mip3, tx->width, tx->height, tx->width/8, tx->height/8);
	}
	for(s32 i = 0; i < m->nummiptex; i++){ // sequence the animations
		texture_t *tx = loadmodel->textures[i];
		if(!tx || tx->name[0] != '+') continue;
		if(tx->anim_next) continue; // allready sequenced
		// find the number of frames in the animation
		memset(anims, 0, sizeof(anims));
		memset(altanims, 0, sizeof(altanims));
		s32 max = tx->name[1];
		s32 altmax = 0;
		if(max >= 'a' && max <= 'z')
			max -= 'a' - 'A';
		if(max >= '0' && max <= '9'){
			max -= '0';
			altmax = 0;
			anims[max] = tx;
			max++;
		} else if(max >= 'A' && max <= 'J'){
			altmax = max - 'A';
			max = 0;
			altanims[altmax] = tx;
			altmax++;
		} else Sys_Error("Bad animating texture %s", tx->name);
		for(s32 j = i + 1; j < m->nummiptex; j++){
			texture_t *tx2 = loadmodel->textures[j];
			if(!tx2 || tx2->name[0] != '+') continue;
			if(strcmp(tx2->name + 2, tx->name + 2)) continue;
			s32 num = tx2->name[1];
			if(num >= 'a' && num <= 'z')
				num -= 'a' - 'A';
			if(num >= '0' && num <= '9'){
				num -= '0';
				anims[num] = tx2;
				if(num + 1 > max)
					max = num + 1;
			} else if(num >= 'A' && num <= 'J'){
				num = num - 'A';
				altanims[num] = tx2;
				if(num + 1 > altmax)
					altmax = num + 1;
			} else
				Sys_Error("Bad animating texture %s", tx->name);
		}
		for(s32 j = 0; j < max; j++){ // link them all together
			texture_t *tx2 = anims[j];
			if(!tx2)Sys_Error("Missing frame %i of %s", j,tx->name);
			tx2->anim_total = max * ANIM_CYCLE;
			tx2->anim_min = j * ANIM_CYCLE;
			tx2->anim_max = (j + 1) * ANIM_CYCLE;
			tx2->anim_next = anims[(j + 1) % max];
			if(altmax)
				tx2->alternate_anims = altanims[0];
		}
		for(s32 j = 0; j < altmax; j++){
			texture_t *tx2 = altanims[j];
			if(!tx2)Sys_Error("Missing frame %i of %s", j,tx->name);
			tx2->anim_total = altmax * ANIM_CYCLE;
			tx2->anim_min = j * ANIM_CYCLE;
			tx2->anim_max = (j + 1) * ANIM_CYCLE;
			tx2->anim_next = altanims[(j + 1) % altmax];
			if(max) tx2->alternate_anims = anims[0];
		}
	}
}

static void Mod_LoadLighting(lump_t *l)
{ // johnfitz -- replaced with lit support code via lordhavoc
	s8 litfilename[MAX_OSPATH];
	loadmodel->lightdata = NULL;
	// LordHavoc: check for a .lit file
	q_strlcpy(litfilename, loadmodel->name, sizeof(litfilename));
	COM_StripExtension(litfilename, litfilename, sizeof(litfilename));
	q_strlcat(litfilename, ".lit", sizeof(litfilename));
	s32 mark = Hunk_LowMark();
	u32 path_id;
	u8 *data = COM_LoadHunkFile(litfilename, &path_id);
	if(data){
		// use lit file only from the same gamedir as the map
		// itself or from a searchpath with higher priority.
		if(path_id < loadmodel->path_id){
			Hunk_FreeToLowMark(mark);
		  Con_DPrintf("ignored %s from a gamedir with lower priority\n",
				litfilename);
		}
		else if(data[0]=='Q'&&data[1]=='L'&&data[2]=='I'&&data[3]=='T'){
			s32 i = LittleLong(((s32 *)data)[1]);
			if(i == 1){
				if(8+l->filelen*3 == com_filesize){
					lit_loaded = 1;
					printf("%s loaded\n", litfilename);
					loadmodel->lightdata = data + 8;
					return;
				}
				Hunk_FreeToLowMark(mark);
		Con_Printf("Outdated .lit file(%s should be %u bytes, not %u)\n"
				, litfilename, 8+l->filelen*3, com_filesize);
			}
			else { Hunk_FreeToLowMark(mark);
			Con_Printf("Unknown .lit file version(%d)\n",i); }
		} else { Hunk_FreeToLowMark(mark);
		Con_Printf("Corrupt .lit file(old version?), ignoring\n");}
	} // LordHavoc: no .lit found, expand the white lighting data to color
	if(!l->filelen) return;
	loadmodel->lightdata = (u8 *)Hunk_AllocName( l->filelen*3, litfilename);
	u8 *in = loadmodel->lightdata + l->filelen*2; // place the file at the
		// end, so it will not be overwritten until the very last write
	u8 *out = loadmodel->lightdata;
	memcpy(in, mod_base + l->fileofs, l->filelen);
	for(s32 i = 0; i < l->filelen; i++) {
		u8 d = *in++;
		*out++ = d; *out++ = d; *out++ = d;
	}
}

static void Mod_LoadVisibility(lump_t *l)
{
	loadmodel->viswarn = 0;
	if(!l->filelen)
	{
		loadmodel->visdata = NULL;
		return;
	}
	loadmodel->visdata = (u8 *) Hunk_AllocName(l->filelen, loadname);
	memcpy(loadmodel->visdata, mod_base + l->fileofs, l->filelen);
}

static void Mod_LoadEdges(lump_t *l, s32 bsp2)
{
	if(bsp2){
		dledge_t *in = (dledge_t *)(mod_base + l->fileofs);
		if(l->filelen % sizeof(*in))
	Sys_Error("MOD_LoadEdges: funny lump size in %s",loadmodel->name);
		s32 count = l->filelen / sizeof(*in);
		medge_t *out = (medge_t *)Hunk_AllocName((count + 13)
				* sizeof(*out), loadname);
		loadmodel->edges = out;
		loadmodel->numedges = count;
		for(s32 i = 0; i < count; i++, in++, out++){
			out->v[0] = LittleLong(in->v[0]);
			out->v[1] = LittleLong(in->v[1]);
		}
	} else {
		dedge_t *in = (dedge_t *)(mod_base + l->fileofs);
		if(l->filelen % sizeof(*in))
	Sys_Error("MOD_LoadEdges: funny lump size in %s",loadmodel->name);
		s32 count = l->filelen / sizeof(*in);
		medge_t *out = (medge_t *) Hunk_AllocName
			((count + 13) * sizeof(*out), loadname);
		loadmodel->edges = out;
		loadmodel->numedges = count;
		for(s32 i = 0; i<count; i++, in++, out++){
			out->v[0] = (u16)LittleShort(in->v[0]);
			out->v[1] = (u16)LittleShort(in->v[1]);
		}
	}
}

static void Mod_LoadVertexes(lump_t *l)
{
	dvertex_t *in = (dvertex_t *)(mod_base + l->fileofs);
	if(l->filelen % sizeof(*in))
	   Sys_Error("MOD_LoadVertexes: funny lump size in %s",loadmodel->name);
	s32 count = l->filelen / sizeof(*in);
	// Manoel Kasimier - skyboxes - extra for skybox 
	// Code taken from the ToChriS engine - Author: Vic
	mvertex_t *out = 
		(mvertex_t*)Hunk_AllocName((count+8)*sizeof(*out),loadname);
	loadmodel->vertexes = out;
	loadmodel->numvertexes = count;
	for(s32 i = 0; i < count; i++, in++, out++){
		out->position[0] = LittleFloat(in->point[0]);
		out->position[1] = LittleFloat(in->point[1]);
		out->position[2] = LittleFloat(in->point[2]);
	}
}

static void Mod_LoadEntities(lump_t *l)
{
	s8 basemapname[MAX_QPATH];
	s8 entfilename[MAX_QPATH+16];
	if(!external_ents.value) goto _load_embedded;
	s32 mark = Hunk_LowMark();
	u32 crc = 0;
	if(l->filelen > 0) crc = CRC_Block(mod_base+l->fileofs, l->filelen-1);
	q_strlcpy(basemapname, loadmodel->name, sizeof(basemapname));
	COM_StripExtension(basemapname, basemapname, sizeof(basemapname));
	snprintf(entfilename,sizeof(entfilename),"%s@%04x.ent",basemapname,crc);
	u32 path_id;
	s8 *ents = (s8*)COM_LoadHunkFile(entfilename, &path_id);
	if(!ents) {
		snprintf(entfilename,sizeof(entfilename),"%s.ent", basemapname);
		ents = (s8*)COM_LoadHunkFile(entfilename, &path_id);
	}
	if(ents) {
		// use ent file only from the same gamedir as the map
		// itself or from a searchpath with higher priority.
		if(path_id < loadmodel->path_id) {
			Hunk_FreeToLowMark(mark);
    Con_DPrintf("ignored %s from a gamedir with lower priority\n", entfilename);
		} else {
			loadmodel->entities = ents;
		   Con_DPrintf("Loaded external entity file %s\n", entfilename);
			return;
		}
	}
_load_embedded:
	if(!l->filelen) {
		loadmodel->entities = NULL;
		return;
	}
	loadmodel->entities = (s8 *) Hunk_AllocName( l->filelen, loadname);
	memcpy(loadmodel->entities, mod_base + l->fileofs, l->filelen);
}

static void Mod_LoadTexinfo(lump_t *l)
{
	texinfo_t *in = (void *)(mod_base + l->fileofs);
	if(l->filelen % sizeof(*in))
	   Sys_Error("MOD_LoadTexinfo: funny lump size in %s", loadmodel->name);
	s32 count = l->filelen / sizeof(*in);
	mtexinfo_t *out = Hunk_AllocName((count+6) * sizeof(*out), loadname);
	loadmodel->texinfo = out;
	loadmodel->numtexinfo = count;
	for(s32 i = 0; i < count; i++, in++, out++){
		for(s32 k = 0; k < 2; k++)
			for(s32 j = 0; j < 4; j++)
				out->vecs[k][j] = LittleFloat(in->vecs[k][j]);
		f32 len1 = VectorLength(out->vecs[0]);
		f32 len2 = VectorLength(out->vecs[1]);
		len1 = (len1 + len2) / 2;
		if(len1 < 0.32)      out->mipadjust = 4;
		else if(len1 < 0.49) out->mipadjust = 3;
		else if(len1 < 0.99) out->mipadjust = 2;
		else                 out->mipadjust = 1;
		s32 miptex = LittleLong(in->miptex);
		out->flags = LittleLong(in->flags);
		if(!loadmodel->textures){
			out->texture = r_notexture_mip; // checkerboard texture
			out->flags = 0;
		} else {
			if(miptex >= loadmodel->numtextures)
				Sys_Error("miptex >= loadmodel->numtextures");
			out->texture = loadmodel->textures[miptex];
			if(!out->texture){
				out->texture = r_notexture_mip; // tex not found
				out->flags = 0;
			}
		}
	}
}

static void CalcSurfaceExtents(msurface_t *s)
{ // Fills in s->texturemins[] and s->extents[]
	f32 mins[2], maxs[2], val;
	mvertex_t *v;
	s32 bmins[2], bmaxs[2];
	mins[0] = mins[1] = FLT_MAX;
	maxs[0] = maxs[1] = -FLT_MAX;
	mtexinfo_t *tex = s->texinfo;
	for(s32 i = 0; i < s->numedges; i++) {
		s32 e = loadmodel->surfedges[s->firstedge+i];
		if(e >= 0) v = &loadmodel->vertexes[loadmodel->edges[e].v[0]];
		else v = &loadmodel->vertexes[loadmodel->edges[-e].v[1]];
		for(s32 j = 0; j < 2; j++) {
		    // The following calculation is sensitive to floating-point
		    // precision. It needs to produce the same result that the
		    // light compiler does, because R_BuildLightMap uses surf->
		    // extents to know the width/height of a surface's lightmap,
		    // and incorrect rounding here manifests itself as patches
		    // of "corrupted" looking lightmaps.
		    // Most light compilers are win32 executables, so they use
		    // x87 floating point. This means the multiplies and adds
		    // are done at 80-bit precision, and the result is rounded
		    // down to 32-bits and stored in val.
		    // Adding the casts to f64 seems to be good enough to fix
		    // lighting glitches when Quakespasm is compiled as x86_64
		    // and using SSE2 floating-point. A potential trouble spot
		    // is the hallway at the beginning of mfxsp17. -- ericw
			val = ((f64)v->position[0] * (f64)tex->vecs[j][0]) +
				((f64)v->position[1] * (f64)tex->vecs[j][1]) +
				((f64)v->position[2] * (f64)tex->vecs[j][2]) +
				(f64)tex->vecs[j][3];
			if(val < mins[j]) mins[j] = val;
			if(val > maxs[j]) maxs[j] = val;
		}
	}
	for(s32 i = 0; i < 2; i++) {
		bmins[i] = floor(mins[i]/16);
		bmaxs[i] = ceil(maxs[i]/16);
		s->texturemins[i] = bmins[i] * 16;
		s->extents[i] = (bmaxs[i] - bmins[i]) * 16;
		if(!(tex->flags & TEX_SPECIAL) && s->extents[i] > 2000)
			Sys_Error("Bad surface extents");
	}
}

static void Mod_CalcSurfaceBounds(msurface_t *s) // -- johnfitz
{ // calculate bounding box for per-surface frustum culling
	s->mins[0] = s->mins[1] = s->mins[2] = FLT_MAX;
	s->maxs[0] = s->maxs[1] = s->maxs[2] = -FLT_MAX;
	for(s32 i = 0; i < s->numedges; i++){
		s32 e = loadmodel->surfedges[s->firstedge+i];
		mvertex_t *v = &loadmodel->vertexes[loadmodel->edges[-e].v[1]];
		if(e >= 0) v = &loadmodel->vertexes[loadmodel->edges[e].v[0]];
		if(s->mins[0] > v->position[0]) s->mins[0] = v->position[0];
		if(s->mins[1] > v->position[1]) s->mins[1] = v->position[1];
		if(s->mins[2] > v->position[2]) s->mins[2] = v->position[2];
		if(s->maxs[0] < v->position[0]) s->maxs[0] = v->position[0];
		if(s->maxs[1] < v->position[1]) s->maxs[1] = v->position[1];
		if(s->maxs[2] < v->position[2]) s->maxs[2] = v->position[2];
	}
}

static void Mod_LoadFaces(lump_t *l, bool bsp2)
{
	dface_t *ins;
	dlface_t *inl;
	s32 i, count, lofs;
	s32 planenum, side, texinfon;
	if(bsp2) {
		ins = NULL;
		inl = (dlface_t *)(mod_base + l->fileofs);
		if(l->filelen % sizeof(*inl))
	     Sys_Error("MOD_LoadBmodel: funny lump size in %s",loadmodel->name);
		count = l->filelen / sizeof(*inl);
	} else {
		ins = (dface_t *)(mod_base + l->fileofs);
		inl = NULL;
		if(l->filelen % sizeof(*ins))
	     Sys_Error("MOD_LoadBmodel: funny lump size in %s",loadmodel->name);
		count = l->filelen / sizeof(*ins);
	}
	msurface_t *out = Hunk_AllocName((count+6)*sizeof(*out), loadname);
	if(count > 32767 && !bsp2)
	      Con_DPrintf("%i faces exceeds standard limit of 32767.\n", count);
	loadmodel->surfaces = out;
	loadmodel->numsurfaces = count;
	for(s32 surfnum = 0; surfnum < count; surfnum++, out++) {
		if(bsp2) {
			out->firstedge = LittleLong(inl->firstedge);
			out->numedges = LittleLong(inl->numedges);
			planenum = LittleLong(inl->planenum);
			side = LittleLong(inl->side);
			texinfon = LittleLong(inl->texinfo);
			for(i=0 ; i<MAXLIGHTMAPS ; i++)
				out->styles[i] = inl->styles[i];
			lofs = LittleLong(inl->lightofs);
			inl++;
		} else {
			out->firstedge = LittleLong(ins->firstedge);
			out->numedges = LittleShort(ins->numedges);
			planenum = LittleShort(ins->planenum);
			side = LittleShort(ins->side);
			texinfon = LittleShort(ins->texinfo);
			for(i=0 ; i<MAXLIGHTMAPS ; i++)
				out->styles[i] = ins->styles[i];
			lofs = LittleLong(ins->lightofs);
			ins++;
		}
		out->flags = 0;
		if(out->numedges < 3)
		  printf("surfnum %d: bad numedges %d\n",surfnum,out->numedges);
		if(side) out->flags |= SURF_PLANEBACK;
		out->plane = loadmodel->planes + planenum;
		out->texinfo = loadmodel->texinfo + texinfon;
		CalcSurfaceExtents(out);
		Mod_CalcSurfaceBounds(out); //jfitz-per-surface frustum culling
		if(lofs == -1) out->samples = NULL; // lighting info
		else out->samples = loadmodel->lightdata + (lofs * 3);
		if(!q_strncasecmp(out->texinfo->texture->name,"sky",3))
			out->flags |= (SURF_DRAWSKY | SURF_DRAWTILED);
		else if(out->texinfo->texture->name[0] == '*' ||
			out->texinfo->texture->name[0] == '!') { // warp surface
			out->flags |= SURF_DRAWTURB;
			if(out->texinfo->flags & TEX_SPECIAL)
				out->flags |= SURF_DRAWTILED;
			else if(out->samples && !loadmodel->haslitwater) {
				Con_DPrintf("Map has lit water\n");
				loadmodel->haslitwater = 1;
			} // detect special liquid types
			if(!strncmp(out->texinfo->texture->name, "*lava", 5))
				out->flags |= SURF_DRAWLAVA;
		       else if(!strncmp(out->texinfo->texture->name,"*slime",6))
				out->flags |= SURF_DRAWSLIME;
			else if(!strncmp(out->texinfo->texture->name,"*tele",5))
				out->flags |= SURF_DRAWTELE;
			else out->flags |= SURF_DRAWWATER;
		} else if(out->texinfo->texture->name[0] == '{'){//ericw--fence
			out->flags |= SURF_DRAWCUTOUT;
		} else if(out->texinfo->flags & TEX_MISSING){// missing from bsp
			if(out->samples)out->flags|=SURF_NOTEXTURE;//lightmapped
			else out->flags|=(SURF_NOTEXTURE|SURF_DRAWTILED);//not.
		}
	}

}

static void Mod_SetParent(mnode_t *node, mnode_t *parent)
{
	node->parent = parent;
	if(node->contents < 0) return;
	Mod_SetParent(node->children[0], node);
	Mod_SetParent(node->children[1], node);
}

static void Mod_LoadNodes_S(lump_t *l)
{
	dnode_t *in = (dnode_t *)(mod_base + l->fileofs);
	if(l->filelen % sizeof(*in))
	    Sys_Error("MOD_LoadNodes_S: funny lump size in %s",loadmodel->name);
	s32 count = l->filelen / sizeof(*in);
	mnode_t *out = Hunk_AllocName( count*sizeof(*out), loadname);
	if(count > 32767) //johnfitz -- warn mappers about exceeding old limits
	      Con_DPrintf("%i nodes exceeds standard limit of 32767.\n", count);
	loadmodel->nodes = out;
	loadmodel->numnodes = count;
	for(s32 i = 0; i < count; i++, in++, out++){
		for(s32 j = 0; j < 3; j++){
			out->minmaxs[j] = LittleShort(in->mins[j]);
			out->minmaxs[3+j] = LittleShort(in->maxs[j]);
		}
		s32 p = LittleLong(in->planenum);
		out->plane = loadmodel->planes + p;
		out->firstsurface = (u16)LittleShort(in->firstface);
		out->numsurfaces = (u16)LittleShort(in->numfaces);
		for(s32 j = 0; j < 2; j++){
			p = (u16)LittleShort(in->children[j]);
			if(p < count)
				out->children[j] = loadmodel->nodes + p;
			else {//note this uses 65535 intentionally, -1 is leaf 0
				p = 65535 - p;
				if(p < loadmodel->numleafs)
			   out->children[j] = (mnode_t *)(loadmodel->leafs + p);
				else {
    Con_Printf("Mod_LoadNodes: invalid leaf index %i(file has only %i leafs)\n",
							p, loadmodel->numleafs);
					out->children[j] =
			(mnode_t *)(loadmodel->leafs);//map it to the solid leaf
				}
			}
		}
	}
}

static void Mod_LoadNodes_L1(lump_t *l)
{
	dl1node_t *in = (dl1node_t *)(mod_base + l->fileofs);
	if(l->filelen % sizeof(*in))
	      Sys_Error("Mod_LoadNodes: funny lump size in %s",loadmodel->name);
	s32 count = l->filelen / sizeof(*in);
	mnode_t *out = (mnode_t *)Hunk_AllocName( count*sizeof(*out), loadname);
	loadmodel->nodes = out;
	loadmodel->numnodes = count;
	for(s32 i = 0; i < count; i++, in++, out++){
		for(s32 j = 0; j < 3; j++){
			out->minmaxs[j] = LittleShort(in->mins[j]);
			out->minmaxs[3+j] = LittleShort(in->maxs[j]);
		}
		s32 p = LittleLong(in->planenum);
		out->plane = loadmodel->planes + p;
		out->firstsurface = LittleLong(in->firstface);
		out->numsurfaces = LittleLong(in->numfaces);
		for(s32 j = 0; j < 2; j++){
			p = LittleLong(in->children[j]);
			if(p >= 0 && p < count)
				out->children[j] = loadmodel->nodes + p;
			else {//note this uses 65535 intentionally, -1 is leaf 0
				p = 0xffffffff - p;
				if(p >= 0 && p < loadmodel->numleafs)
					out->children[j] =
					    (mnode_t *)(loadmodel->leafs + p);
				else {
    Con_Printf("Mod_LoadNodes: invalid leaf index %i(file has only %i leafs)\n",
						        p, loadmodel->numleafs);
					out->children[j] =
			(mnode_t *)(loadmodel->leafs);//map it to the solid leaf
				}
			}
		}
	}
}

static void Mod_LoadNodes_L2(lump_t *l)
{
	dl2node_t *in = (dl2node_t *)(mod_base + l->fileofs);
	if(l->filelen % sizeof(*in))
	      Sys_Error("Mod_LoadNodes: funny lump size in %s",loadmodel->name);
	s32 count = l->filelen / sizeof(*in);
	mnode_t *out = (mnode_t *)Hunk_AllocName( count*sizeof(*out), loadname);
	loadmodel->nodes = out;
	loadmodel->numnodes = count;
	for(s32 i = 0; i < count; i++, in++, out++){
		for(s32 j = 0; j < 3; j++){
			out->minmaxs[j] = LittleFloat(in->mins[j]);
			out->minmaxs[3+j] = LittleFloat(in->maxs[j]);
		}
		s32 p = LittleLong(in->planenum);
		out->plane = loadmodel->planes + p;
		out->firstsurface = LittleLong(in->firstface);
		out->numsurfaces = LittleLong(in->numfaces);
		for(s32 j = 0; j < 2; j++){
			p = LittleLong(in->children[j]);
			if(p > 0 && p < count)
				out->children[j] = loadmodel->nodes + p;
			else {//note this uses 65535 intentionally, -1 is leaf 0
				p = 0xffffffff - p;
				if(p >= 0 && p < loadmodel->numleafs)
					out->children[j] =
					    (mnode_t *)(loadmodel->leafs + p);
				else {
    Con_Printf("Mod_LoadNodes: invalid leaf index %i(file has only %i leafs)\n",
		    					p, loadmodel->numleafs);
					out->children[j] =
			(mnode_t *)(loadmodel->leafs);//map it to the solid leaf
				}
			}
		}
	}
}

static void Mod_LoadNodes(lump_t *l, s32 bsp2)
{
	if(bsp2 == 2) Mod_LoadNodes_L2(l);
	else if(bsp2) Mod_LoadNodes_L1(l);
	else Mod_LoadNodes_S(l);
	Mod_SetParent(loadmodel->nodes, NULL); // sets nodes and leafs
}

static void Mod_ProcessLeafs_S(dleaf_t *in, s32 filelen)
{
	if(filelen % sizeof(*in))
	  Sys_Error("Mod_ProcessLeafs: funny lump size in %s", loadmodel->name);
	s32 count = filelen / sizeof(*in);
	mleaf_t *out = (mleaf_t *)Hunk_AllocName( count*sizeof(*out), loadname);
	if(count > 32767) // johnfitz
	   Host_Error("Mod_LoadLeafs: %i leafs exceeds limit of 32767.", count);
	loadmodel->leafs = out;
	loadmodel->numleafs = count;
	for(s32 i = 0; i < count; i++, in++, out++){
		for(s32 j = 0; j < 3; j++){
			out->minmaxs[j] = LittleShort(in->mins[j]);
			out->minmaxs[3+j] = LittleShort(in->maxs[j]);
		}
		s32 p = LittleLong(in->contents);
		out->contents = p;
		out->firstmarksurface = loadmodel->marksurfaces 
				+ (u16)LittleShort(in->firstmarksurface);
		out->nummarksurfaces = (u16)LittleShort(in->nummarksurfaces);
		p = LittleLong(in->visofs);
		if(p == -1) out->compressed_vis = NULL;
		else out->compressed_vis = loadmodel->visdata + p;
		out->efrags = NULL;
		for(s32 j = 0; j < 4; j++)
			out->ambient_sound_level[j] = in->ambient_level[j];
	}
}

static void Mod_ProcessLeafs_L1(dl1leaf_t *in, s32 filelen)
{
	if(filelen % sizeof(*in))
	   Sys_Error("Mod_ProcessLeafs: funny lump size in %s",loadmodel->name);
	s32 count = filelen / sizeof(*in);
	mleaf_t *out = (mleaf_t *) Hunk_AllocName(count*sizeof(*out), loadname);
	loadmodel->leafs = out;
	loadmodel->numleafs = count;
	for(s32 i = 0; i < count; i++, in++, out++){
		for(s32 j = 0; j < 3; j++){
			out->minmaxs[j] = LittleShort(in->mins[j]);
			out->minmaxs[3+j] = LittleShort(in->maxs[j]);
		}
		s32 p = LittleLong(in->contents);
		out->contents = p;
		out->firstmarksurface = loadmodel->marksurfaces 
			+ LittleLong(in->firstmarksurface);
		out->nummarksurfaces = LittleLong(in->nummarksurfaces);
		p = LittleLong(in->visofs);
		if(p == -1) out->compressed_vis = NULL;
		else out->compressed_vis = loadmodel->visdata + p;
		out->efrags = NULL;
		for(s32 j = 0; j < 4; j++)
			out->ambient_sound_level[j] = in->ambient_level[j];
	}
}

static void Mod_ProcessLeafs_L2(dl2leaf_t *in, s32 filelen)
{
	if(filelen % sizeof(*in))
	  Sys_Error("Mod_ProcessLeafs: funny lump size in %s", loadmodel->name);
	s32 count = filelen / sizeof(*in);
	mleaf_t *out = (mleaf_t *) Hunk_AllocName(count*sizeof(*out), loadname);
	loadmodel->leafs = out;
	loadmodel->numleafs = count;
	for(s32 i = 0; i < count; i++, in++, out++){
		for(s32 j = 0; j < 3; j++){
			out->minmaxs[j] = LittleFloat(in->mins[j]);
			out->minmaxs[3+j] = LittleFloat(in->maxs[j]);
		}
		s32 p = LittleLong(in->contents);
		out->contents = p;
		out->firstmarksurface = loadmodel->marksurfaces
			+ LittleLong(in->firstmarksurface);
		out->nummarksurfaces = LittleLong(in->nummarksurfaces);
		p = LittleLong(in->visofs);
		if(p == -1) out->compressed_vis = NULL;
		else out->compressed_vis = loadmodel->visdata + p;
		out->efrags = NULL;
		for(s32 j = 0; j < 4; j++)
			out->ambient_sound_level[j] = in->ambient_level[j];
	}
}

static void Mod_LoadLeafs(lump_t *l, s32 bsp2)
{
	void *in = (void *)(mod_base + l->fileofs);
	if(bsp2 == 2) Mod_ProcessLeafs_L2((dl2leaf_t *)in, l->filelen);
	else if(bsp2) Mod_ProcessLeafs_L1((dl1leaf_t *)in, l->filelen);
	else Mod_ProcessLeafs_S((dleaf_t *) in, l->filelen);
}

static void Mod_CheckWaterVis()
{
	mleaf_t *leaf, *other;
	msurface_t * surf;
	s32 i, j, k;
	s32 numclusters = loadmodel->submodels[0].visleafs;
	s32 contentfound = 0;
	s32 contenttransparent = 0;
	s32 contenttype;
	u32 hascontents = 0;
	if (r_novis.value) { // all can be
		loadmodel->contentstransparent = (SURF_DRAWWATER|SURF_DRAWTELE|
				SURF_DRAWSLIME|SURF_DRAWLAVA);
		return;
	}
	//pvs is 1-based. leaf 0 sees all (the solid leaf).
	//leaf 0 has no pvs, and does not appear in other leafs either, so watch out for the biases.
	for (i=0,leaf=loadmodel->leafs+1 ; i<numclusters ; i++, leaf++) {
		u8 *vis;
		if (leaf->contents < 0) //err... wtf?
			hascontents |= 1u<<-leaf->contents;
		if (leaf->contents == CONTENTS_WATER) {
			if ((contenttransparent & (SURF_DRAWWATER|SURF_DRAWTELE))==(SURF_DRAWWATER|SURF_DRAWTELE))
				continue;
			//this check is somewhat risky, but we should be able to get away with it.
			for (contenttype = 0, j = 0; j < leaf->nummarksurfaces; j++) {
				surf = leaf->firstmarksurface[j];
				if (surf->flags & (SURF_DRAWWATER|SURF_DRAWTELE)) {
					contenttype = surf->flags & (SURF_DRAWWATER|SURF_DRAWTELE);
					break;
				}
			}
			//its possible that this leaf has absolutely no surfaces in it, turb or otherwise.
			if (contenttype == 0) continue;
		}
		else if (leaf->contents == CONTENTS_SLIME)
			contenttype = SURF_DRAWSLIME;
		else if (leaf->contents == CONTENTS_LAVA)
			contenttype = SURF_DRAWLAVA;
		//fixme: tele
		else continue;
		if (contenttransparent & contenttype) {
nextleaf:
			continue;   //found one of this type already
		}
		contentfound |= contenttype;
		vis = Mod_DecompressVis(leaf->compressed_vis, loadmodel);
		for (j = 0; j < (numclusters+7)/8; j++) {
			if (!vis[j]) continue;
			for (k = 0; k < 8; k++) {
				if (!(vis[j] & (1u<<k))) continue; 
				other = &loadmodel->leafs[(j<<3)+k+1];
				if (leaf->contents != other->contents) {
					contenttransparent |= contenttype;
					goto nextleaf;
				}
			}
		}
	}
	if (!contenttransparent) {   //no water leaf saw a non-water leaf
	    //but only warn when there's actually water somewhere there...
		if (hascontents & ((1<<-CONTENTS_WATER)
				|  (1<<-CONTENTS_SLIME)
				|  (1<<-CONTENTS_LAVA)))
			Con_DPrintf("%s is not watervised\n", loadmodel->name);
	} else {
		Con_DPrintf("%s is vised for transparent", loadmodel->name);
		if (contenttransparent & SURF_DRAWWATER)
			Con_DPrintf(" water");
		if (contenttransparent & SURF_DRAWTELE)
			Con_DPrintf(" tele");
		if (contenttransparent & SURF_DRAWLAVA)
			Con_DPrintf(" lava");
		if (contenttransparent & SURF_DRAWSLIME)
			Con_DPrintf(" slime");
		Con_DPrintf("\n");
		if (((r_wateralpha.value > 0 && r_wateralpha.value < 1) || 
			(r_wateralpha.value > 0 && r_wateralpha.value < 1) ||
			(r_wateralpha.value > 0 && r_wateralpha.value < 1) ||
			(r_wateralpha.value > 0 && r_wateralpha.value < 1))
			&& r_alphastyle.value != 1) R_BuildColorMixLUT(0);
	}
	// any types that we didn't find are assumed to be transparent.
	// this allows submodels to work okay
	// (eg: ad uses func_illusionary teleporters for some reason).
	loadmodel->contentstransparent = contenttransparent | (~contentfound &
		(SURF_DRAWWATER|SURF_DRAWTELE|SURF_DRAWSLIME|SURF_DRAWLAVA));
}

static void Mod_LoadClipnodes(lump_t *l, bool bsp2)
{
	dsclipnode_t *ins;
	dlclipnode_t *inl;
	s32 count;
	if(bsp2){
		ins = NULL;
		inl = (dlclipnode_t *)(mod_base + l->fileofs);
		if(l->filelen % sizeof(*inl))
	  Sys_Error("Mod_LoadClipnodes: funny lump size in %s",loadmodel->name);
		count = l->filelen / sizeof(*inl);
	} else {
		ins = (dsclipnode_t *)(mod_base + l->fileofs);
		inl = NULL;
		if(l->filelen % sizeof(*ins))
	  Sys_Error("Mod_LoadClipnodes: funny lump size in %s",loadmodel->name);
		count = l->filelen / sizeof(*ins);
	}
	mclipnode_t *out;
	out = (mclipnode_t *) Hunk_AllocName( count*sizeof(*out), loadname);
	if(count > 32767 && !bsp2) //johnfitz -- warn about exceeding old limits
	  Con_DPrintf("%i clipnodes exceeds standard limit of 32767.\n", count);
	loadmodel->clipnodes = out;
	loadmodel->numclipnodes = count;
	hull_t *hull = &loadmodel->hulls[1];
	hull->clipnodes = out;
	hull->firstclipnode = 0;
	hull->lastclipnode = count-1;
	hull->planes = loadmodel->planes;
	hull->clip_mins[0] = -16;
	hull->clip_mins[1] = -16;
	hull->clip_mins[2] = -24;
	hull->clip_maxs[0] = 16;
	hull->clip_maxs[1] = 16;
	hull->clip_maxs[2] = 32;
	hull = &loadmodel->hulls[2];
	hull->clipnodes = out;
	hull->firstclipnode = 0;
	hull->lastclipnode = count-1;
	hull->planes = loadmodel->planes;
	hull->clip_mins[0] = -32;
	hull->clip_mins[1] = -32;
	hull->clip_mins[2] = -24;
	hull->clip_maxs[0] = 32;
	hull->clip_maxs[1] = 32;
	hull->clip_maxs[2] = 64;
	if(bsp2){ for(s32 i = 0; i < count; i++, out++, inl++){ // johnfitz
		out->planenum = LittleLong(inl->planenum);//bounds check
		if(out->planenum < 0 || out->planenum >= loadmodel->numplanes)
			Host_Error("Mod_LoadClipnodes: planenum out of bounds");
		out->children[0] = LittleLong(inl->children[0]);
		out->children[1] = LittleLong(inl->children[1]);
	} } else { for(s32 i = 0; i < count; i++, out++, ins++){
		out->planenum = LittleLong(ins->planenum);
		if(out->planenum < 0 || out->planenum >= loadmodel->numplanes)
			Host_Error("Mod_LoadClipnodes: planenum out of bounds");
		out->children[0] = (u16)LittleShort(ins->children[0]);
		out->children[1] = (u16)LittleShort(ins->children[1]);
		if(out->children[0] >= count) out->children[0] -= 65536;
		if(out->children[1] >= count) out->children[1] -= 65536;
	} }
}


static void Mod_MakeHull0()
{ // Duplicate the drawing hull structure as a clipping hull
	hull_t *hull = &loadmodel->hulls[0];
	mnode_t *in = loadmodel->nodes;
	s32 count = loadmodel->numnodes;
	mclipnode_t *out = Hunk_AllocName(count*sizeof(*out), loadname);
	hull->clipnodes = out;
	hull->firstclipnode = 0;
	hull->lastclipnode = count-1;
	hull->planes = loadmodel->planes;
	for(s32 i = 0; i < count; i++, out++, in++){
		out->planenum = in->plane - loadmodel->planes;
		for(s32 j = 0; j < 2; j++){
			mnode_t *child = in->children[j];
			if(child->contents<0) out->children[j]=child->contents;
			else out->children[j] = child - loadmodel->nodes;
		}
	}
}

static void Mod_LoadMarksurfaces(lump_t *l, s32 bsp2)
{
	if(bsp2){
		u32 *in = (u32 *)(mod_base + l->fileofs);
		if(l->filelen % sizeof(*in))
      Host_Error("Mod_LoadMarksurfaces: funny lump size in %s",loadmodel->name);
		s32 count = l->filelen / sizeof(*in);
		msurface_t **out = Hunk_AllocName(count*sizeof(*out), loadname);
		loadmodel->marksurfaces = out;
		loadmodel->nummarksurfaces = count;
		for(s32 i = 0; i < count; i++){
			u64 j = LittleLong(in[i]);
			if(j >= (u64)loadmodel->numsurfaces)
			 Host_Error("Mod_LoadMarksurfaces: bad surface number");
			out[i] = loadmodel->surfaces + j;
		}
	} else {
		s16 *in = (s16 *)(mod_base + l->fileofs);
		if(l->filelen % sizeof(*in))
      Host_Error("Mod_LoadMarksurfaces: funny lump size in %s",loadmodel->name);
		s32 count = l->filelen / sizeof(*in);
		msurface_t **out = Hunk_AllocName(count*sizeof(*out), loadname);
		loadmodel->marksurfaces = out;
		loadmodel->nummarksurfaces = count;
		if(count > 32767)
	Con_DPrintf("%i marksurfaces exceeds standard limit of 32767.\n",count);
		for(s32 i = 0; i < count; i++){
			u16 j = (u16)LittleShort(in[i]);
			if(j >= loadmodel->numsurfaces)
			  Sys_Error("Mod_LoadMarksurfaces: bad surface number");
			out[i] = loadmodel->surfaces + j;
		}
	}
}

static void Mod_LoadSurfedges(lump_t *l)
{
	s32 *in = (s32 *)(mod_base + l->fileofs);
	if(l->filelen % sizeof(*in))
	  Sys_Error("MOD_LoadSurfedges: funny lump size in %s",loadmodel->name);
	s32 count = l->filelen / sizeof(*in);
	s32 *out = (s32 *) Hunk_AllocName( (count+24)*sizeof(*out), loadname);
	loadmodel->surfedges = out;
	loadmodel->numsurfedges = count;
	for(s32 i = 0; i < count; i++)
		out[i] = LittleLong(in[i]);
}


static void Mod_LoadPlanes(lump_t *l)
{
	dplane_t *in = (dplane_t *)(mod_base + l->fileofs);
	if(l->filelen % sizeof(*in))
	     Sys_Error("MOD_LoadPlanes: funny lump size in %s",loadmodel->name);
	s32 count = l->filelen / sizeof(*in);
	mplane_t *out = Hunk_AllocName( (count+6)*2*sizeof(*out), loadname);
	loadmodel->planes = out;
	loadmodel->numplanes = count;
	for(s32 i = 0; i < count; i++, in++, out++){
		s32 bits = 0;
		for(s32 j = 0; j < 3; j++){
			out->normal[j] = LittleFloat(in->normal[j]);
			if(out->normal[j] < 0)
				bits |= 1<<j;
		}
		out->dist = LittleFloat(in->dist);
		out->type = LittleLong(in->type);
		out->signbits = bits;
	}
}

static f32 RadiusFromBounds(vec3_t min, vec3_t max)
{
	vec3_t corner;
	for(s32 i = 0; i < 3; i++)
		corner[i]=fabs(min[i])>fabs(max[i])?fabs(min[i]):fabs(max[i]);
	return VectorLength(corner);
}

static void Mod_LoadSubmodels(lump_t *l)
{
	dmodel_t *in = (dmodel_t *)(mod_base + l->fileofs);
	if(l->filelen % sizeof(*in))
	  Sys_Error("MOD_LoadSubmodels: funny lump size in %s",loadmodel->name);
	s32 count = l->filelen / sizeof(*in);
	dmodel_t *out = Hunk_AllocName(count*sizeof(*out), loadname);
	loadmodel->submodels = out;
	loadmodel->numsubmodels = count;
	for(s32 i = 0; i < count; i++, in++, out++){
		for(s32 j = 0; j < 3; j++){ // spread the mins / maxs by a pixel
			out->mins[j] = LittleFloat(in->mins[j]) - 1;
			out->maxs[j] = LittleFloat(in->maxs[j]) + 1;
			out->origin[j] = LittleFloat(in->origin[j]);
		}
		for(s32 j = 0; j < MAX_MAP_HULLS; j++)
			out->headnode[j] = LittleLong(in->headnode[j]);
		out->visleafs = LittleLong(in->visleafs);
		out->firstface = LittleLong(in->firstface);
		out->numfaces = LittleLong(in->numfaces);
	}
	out = loadmodel->submodels;
	if(out->visleafs > 8192) // johnfitz -- check world visleafs, from bjp
     Con_DPrintf("%i visleafs exceeds standard limit of 8192.\n",out->visleafs);
}

static FILE *Mod_FindVisibilityExternal()
{
	vispatch_t header;
	s8 visfilename[MAX_QPATH];
	u32 path_id;
	FILE *f;
	size_t r;
	snprintf(visfilename, sizeof(visfilename), "maps/%s.vis", loadname);
	if(COM_FOpenFile(visfilename, &f, &path_id) < 0) {
		Con_DPrintf("%s not found, trying ", visfilename);
		snprintf(visfilename, sizeof(visfilename), "%s.vis",
					COM_SkipPath(com_gamedir));
		Con_DPrintf("%s\n", visfilename);
		if(COM_FOpenFile(visfilename, &f, &path_id) < 0) {
			Con_DPrintf("external vis not found\n");
			return NULL;
		}
	}
	if(path_id < loadmodel->path_id) {
		fclose(f);
     Con_DPrintf("ignored %s from a gamedir with lower priority\n",visfilename);
		return NULL;
	}
	Con_DPrintf("Found external VIS %s\n", visfilename);
	const s8 *shortname = COM_SkipPath(loadmodel->name);
	s64 pos = 0;
	while((r=fread(&header,1,VISPATCH_HEADER_LEN,f))==VISPATCH_HEADER_LEN){
		header.filelen = LittleLong(header.filelen);
		if(header.filelen <= 0){ // bad entry -- don't trust the rest.
			fclose(f);
			return NULL;
		}
		if(!q_strcasecmp(header.mapname, shortname)) break;
		pos += header.filelen + VISPATCH_HEADER_LEN;
		fseek(f, pos, SEEK_SET);
	}
	if(r != VISPATCH_HEADER_LEN){
		fclose(f);
		Con_DPrintf("%s not found in %s\n", shortname, visfilename);
		return NULL;
	}
	return f;
}

static u8 *Mod_LoadVisibilityExternal(FILE *f)
{
	s32 filelen = 0;
	if(!fread(&filelen, 4, 1, f)) return NULL;
	filelen = LittleLong(filelen);
	if(filelen <= 0) return NULL;
	Con_DPrintf("...%d bytes visibility data\n", filelen);
	u8 *visdata = (u8 *) Hunk_AllocName(filelen, "EXT_VIS");
	if(!fread(visdata, filelen, 1, f)) return NULL;
	return visdata;
}

static void Mod_LoadLeafsExternal(FILE *f)
{
	s32 filelen = 0;
	if(!fread(&filelen, 4, 1, f)) return;
	filelen = LittleLong(filelen);
	if(filelen <= 0) return;
	Con_DPrintf("...%d bytes leaf data\n", filelen);
	void *in = Hunk_AllocName(filelen, "EXT_LEAF");
	if(!fread(in, filelen, 1, f)) return;
	Mod_ProcessLeafs_S((dleaf_t *)in, filelen);
}

static void Mod_LoadBrushModel(model_t *mod, void *buffer)
{
	loadmodel->type = mod_brush;
	dheader_t *header = (dheader_t *)buffer;
	mod->bspversion = LittleLong(header->version);
	s32 bsp2 = 0;
	switch(mod->bspversion){
		case BSPVERSION: bsp2 = 0; break;
		case BSP2VERSION_2PSB: bsp2 = 1; break; //first iteration
		case BSP2VERSION_BSP2: bsp2 = 2; break; //sanitised revision
		case BSPVERSION_QUAKE64: bsp2 = 0; break;
		default:
	  Sys_Error("Mod_LoadBrushModel: %s has unsupported version number(%i)",
			mod->name, mod->bspversion);
			break;
	}
	mod_base = (u8 *)header; // swap all the lumps
	for(s32 i = 0; i < (s32)sizeof(dheader_t) / 4; i++)
		((s32 *)header)[i] = LittleLong( ((s32 *)header)[i]);
	Mod_LoadVertexes(&header->lumps[LUMP_VERTEXES]); // load into heap
	Mod_LoadEdges(&header->lumps[LUMP_EDGES], bsp2);
	Mod_LoadSurfedges(&header->lumps[LUMP_SURFEDGES]);
	Mod_LoadEntities(&header->lumps[LUMP_ENTITIES]);
	Mod_LoadTextures(&header->lumps[LUMP_TEXTURES]);
	Mod_LoadLighting(&header->lumps[LUMP_LIGHTING]);
	Mod_LoadPlanes(&header->lumps[LUMP_PLANES]);
	Mod_LoadTexinfo(&header->lumps[LUMP_TEXINFO]);
	Mod_LoadFaces(&header->lumps[LUMP_FACES], bsp2);
	Mod_LoadMarksurfaces(&header->lumps[LUMP_MARKSURFACES], bsp2);
	if(mod->bspversion == BSPVERSION && external_vis.value &&
			sv.modelname[0] && !q_strcasecmp(loadname, sv.name)){
		Con_DPrintf("trying to open external vis file\n");
		FILE *fvis = Mod_FindVisibilityExternal();
		if(fvis){
			s32 mark = Hunk_LowMark();
			loadmodel->leafs = NULL;
			loadmodel->numleafs = 0;
			Con_DPrintf("found valid external .vis file for map\n");
			loadmodel->visdata = Mod_LoadVisibilityExternal(fvis);
			if(loadmodel->visdata){
				Mod_LoadLeafsExternal(fvis);
			}
			fclose(fvis);
			if(loadmodel->visdata && loadmodel->leafs
				&& loadmodel->numleafs) goto visdone;
			Hunk_FreeToLowMark(mark);
		 Con_DPrintf("External VIS data failed, using standard vis.\n");
		}
	}
	Mod_LoadVisibility(&header->lumps[LUMP_VISIBILITY]);
	Mod_LoadLeafs(&header->lumps[LUMP_LEAFS], bsp2);
visdone:
	Mod_LoadNodes(&header->lumps[LUMP_NODES], bsp2);
	Mod_LoadClipnodes(&header->lumps[LUMP_CLIPNODES], bsp2);
	Mod_LoadSubmodels(&header->lumps[LUMP_MODELS]);
	Mod_MakeHull0();
	mod->numframes = 2; // regular and alternate animation
	Mod_CheckWaterVis();
	// set up the submodels (FIXME: this is confusing)
// johnfitz -- okay, so that i stop getting confused every time i look at this 
// loop, here's how it works: we're looping through the submodels starting at 0.
// Submodel 0 is the main model, so we don't have to worry about clobbering data
// the first time through, since it's the same data. At the end of the loop, we
// create a new copy of the data to use the next time through.
	for(s32 i = 0; i < mod->numsubmodels; i++){
		dmodel_t *bm = &mod->submodels[i];
		mod->hulls[0].firstclipnode = bm->headnode[0];
		for(s32 j = 1; j < MAX_MAP_HULLS; j++){
			mod->hulls[j].firstclipnode = bm->headnode[j];
			mod->hulls[j].lastclipnode = mod->numclipnodes-1;
		}
		mod->firstmodelsurface = bm->firstface;
		mod->nummodelsurfaces = bm->numfaces;
		VectorCopy(bm->maxs, mod->maxs);
		VectorCopy(bm->mins, mod->mins);
		//johnfitz -- calculate rotate bounds and yaw bounds
		f32 radius = RadiusFromBounds(mod->mins, mod->maxs);
		mod->rmaxs[0] = mod->rmaxs[1] = mod->rmaxs[2] = mod->ymaxs[0] =
			mod->ymaxs[1] = mod->ymaxs[2] = radius;
		mod->rmins[0] = mod->rmins[1] = mod->rmins[2] = mod->ymins[0] =
			mod->ymins[1] = mod->ymins[2] = -radius;
//johnfitz -- correct physics cullboxes so that outlying clip brushes on doors
//and stuff are handled right. skip submodel 0 of sv.worldmodel-the actual world
		if(i > 0 || strcmp(mod->name, sv.modelname) != 0){
			VectorCopy(mod->maxs, mod->clipmaxs);
			VectorCopy(mod->mins, mod->clipmins);
		}
		mod->numleafs = bm->visleafs;
		if(i < mod->numsubmodels-1){ // duplicate the basic information
			s8 name[12];
			sprintf(name, "*%i", i+1);
			loadmodel = Mod_FindName(name);
			*loadmodel = *mod;
			strcpy(loadmodel->name, name);
			mod = loadmodel;
		}
	}
}

void *Mod_LoadAliasFrame(void *pin, s32 *pframeindex, s32 numv,
	SDL_UNUSED trivertx_t *pbboxmin, SDL_UNUSED trivertx_t *pbboxmax, aliashdr_t *pheader,
	s8 *name, maliasframedesc_t *frame, s32 recursed) {
	daliasframe_t *pdaliasframe = (daliasframe_t *) pin;
	q_strlcpy(name, pdaliasframe->name, 16);
	if(!recursed){
		frame->firstpose = posenum;
		frame->numposes = 1;
	}
	for(s32 i = 0; i < 3; i++){
		frame->bboxmin.v[i] = pdaliasframe->bboxmin.v[i];
		frame->bboxmax.v[i] = pdaliasframe->bboxmax.v[i];
	}
	trivertx_t *pinframe = (trivertx_t *) (pdaliasframe + 1);
	trivertx_t *pframe = Hunk_AllocName(numv * sizeof(*pframe), loadname);
	*pframeindex = (u8 *) pframe - (u8 *) pheader;
	poseverts[posenum] = pinframe;
	posenum++;
	for(s32 j = 0; j < numv; j++){
		pframe[j].lightnormalindex = pinframe[j].lightnormalindex;
		for(s32 k = 0; k < 3; k++)
			pframe[j].v[k] = pinframe[j].v[k];
	}
	pinframe += numv;
	return(void *)pinframe;
}

void *Mod_LoadAliasGroup(void *pin, s32 *pframeindex, s32 numv,
		SDL_UNUSED trivertx_t *pbboxmin, SDL_UNUSED trivertx_t *pbboxmax,
		aliashdr_t *pheader, s8 *name, maliasframedesc_t *frame)
{
	daliasgroup_t *pingroup = (daliasgroup_t *) pin;
	s32 numframes = LittleLong(pingroup->numframes);
	frame->firstpose = posenum;
	maliasgroup_t *paliasgroup = Hunk_AllocName(sizeof(maliasgroup_t) +
		(numframes - 1) * sizeof(paliasgroup->frames[0]), loadname);
	paliasgroup->numframes = frame->numposes = numframes;
	for(s32 i = 0; i < 3; i++){
		frame->bboxmin.v[i] = pingroup->bboxmin.v[i];
		frame->bboxmax.v[i] = pingroup->bboxmax.v[i];
	}
	*pframeindex = (u8 *) paliasgroup - (u8 *) pheader;
	daliasinterval_t *pin_intervals = (daliasinterval_t *) (pingroup + 1);
	frame->interval = LittleFloat(pin_intervals->interval);
	f32 *poutintervals = Hunk_AllocName(numframes*sizeof(f32),loadname);
	paliasgroup->intervals = (u8 *) poutintervals - (u8 *) pheader;
	for(s32 i = 0; i < numframes; i++){
		*poutintervals = LittleFloat(pin_intervals->interval);
		if(*poutintervals <= 0.0)
			Sys_Error("Mod_LoadAliasGroup: interval<=0");
		poutintervals++;
		pin_intervals++;
	}
	void *ptemp = (void *)pin_intervals;
	for(s32 i = 0; i < numframes; i++)
		ptemp = Mod_LoadAliasFrame(ptemp, &paliasgroup->frames[i].frame,
			numv, &paliasgroup->frames[i].bboxmin,
			&paliasgroup->frames[i].bboxmax, pheader, name,frame,1);
	return ptemp;
}

static void Mod_CalcAliasBounds(aliashdr_t *a) // johnfitz -- calculate bounds
{ // of alias model for nonrotated, yawrotated, and fullrotated cases
	for(s32 i = 0; i < 3; i++){ //clear out all data
	    loadmodel->mins[i]=loadmodel->ymins[i]=loadmodel->rmins[i]=999999;
	    loadmodel->maxs[i]=loadmodel->ymaxs[i]=loadmodel->rmaxs[i]=-999999;
	}
	f32 radius = 0;
	f32 yawradius = 0;
	for(s32 i = 0; i < a->numposes; i++) //process verts
		for(s32 j = 0; j < a->numverts; j++){
			vec3_t v;
			for(s32 k = 0; k < 3; k++)
				v[k] = poseverts[i][j].v[k]*pheader->scale[k]
					+pheader->scale_origin[k];
			for(s32 k = 0; k < 3; k++){
			    loadmodel->mins[k] = fmin(loadmodel->mins[k], v[k]);
			    loadmodel->maxs[k] = fmax(loadmodel->maxs[k], v[k]);
			}
			f32 dist = v[0] * v[0] + v[1] * v[1];
			if(yawradius < dist)
				yawradius = dist;
			dist += v[2] * v[2];
			if(radius < dist)
				radius = dist;
		}
}

bool nameInList(const s8 *list, const s8 *name)
{
	s8 tmp[MAX_QPATH];
	const s8 *s = list;
	while(*s){
		s32 i = 0; // make a copy until the next comma or end of string
		while(*s && *s != ','){
			if(i < MAX_QPATH - 1) tmp[i++] = *s;
			s++;
		}
		tmp[i] = '\0';
		if(!strcmp(name, tmp)) //compare it to the model name
			return 1;
		while(*s && *s == ',')
			s++; // search forwards to next comma or end of string
	}
	return 0;
}

void Mod_SetExtraFlags(model_t *mod)
{ // johnfitz -- set up extra flags that aren't in the mdl
	if(!mod || mod->type != mod_alias)
		return;
	mod->flags &= (0xFF | MF_HOLEY); //only preserve first u8, plus MF_HOLEY
	if(nameInList(r_nolerp_list.string, mod->name))
		mod->flags |= MOD_NOLERP;
	if(nameInList(r_noshadow_list.string, mod->name))
		mod->flags |= MOD_NOSHADOW;
	// fullbright hack(TODO: make this a cvar list)
	if(!strcmp(mod->name, "progs/flame2.mdl") ||
			!strcmp(mod->name, "progs/flame.mdl") ||
			!strcmp(mod->name, "progs/boss.mdl"))
		mod->flags |= MOD_FBRIGHTHACK;
}

void *Mod_LoadAliasSkin(void *pin, s32 *pskinindex, s32 skinsize,
		aliashdr_t *pheader)
{
	u8 *pskin = Hunk_AllocName(skinsize, loadname);
	u8 *pinskin = (u8 *)pin;
	*pskinindex = (u8 *)pskin - (u8 *)pheader;
	memcpy(pskin, pinskin, skinsize);
	pinskin += skinsize;
	return((void *)pinskin);
}

void *Mod_LoadAliasSkinGroup(void *pin, s32 *pskinindex, s32 skinsize,
		aliashdr_t *pheader)
{
	daliasskingroup_t *pinskingroup = (daliasskingroup_t *) pin;
	s32 numskins = LittleLong(pinskingroup->numskins);
	maliasskingroup_t *paliasskingroup =
		Hunk_AllocName(sizeof(maliasskingroup_t) + (numskins - 1) *
				sizeof(paliasskingroup->skindescs[0]),loadname);
	paliasskingroup->numskins = numskins;
	*pskinindex = (u8 *) paliasskingroup - (u8 *) pheader;
	daliasskininterval_t *pinskinintervals =
		(daliasskininterval_t *) (pinskingroup + 1);
	f32 *poutskinintervals = Hunk_AllocName(numskins*sizeof(f32), loadname);
	paliasskingroup->intervals = (u8*)poutskinintervals - (u8*)pheader;
	for(s32 i = 0; i < numskins; i++){
		*poutskinintervals = LittleFloat(pinskinintervals->interval);
		if(*poutskinintervals <= 0)
			Sys_Error("Mod_LoadAliasSkinGroup: interval<=0");
		poutskinintervals++;
		pinskinintervals++;
	}
	void *ptemp = (void *)pinskinintervals;
	for(s32 i = 0; i < numskins; i++)
		ptemp = Mod_LoadAliasSkin(ptemp,
			&paliasskingroup->skindescs[i].skin, skinsize, pheader);
	return ptemp;
}

void Mod_LoadAliasModel(model_t *mod, void *buffer)
{
	stvert_t *pstverts, *pinstverts;
	mtriangle_t *ptri;
	dtriangle_t *pintriangles;
	daliasframetype_t *pframetype;
	daliasskintype_t *pskintype;
	maliasskindesc_t *pskindesc;
	s32 start = Hunk_LowMark();
	mdl_t *pinmodel = buffer;
	mod_base = (u8 *)buffer; //johnfitz
	s32 version = LittleLong(pinmodel->version);
	if(version != ALIAS_VERSION) Host_Error(
	     "Mod_LoadAliasModel: %s has wrong version number(%d should be %d)",
				mod->name, version, ALIAS_VERSION);
	// allocate space for a working header, plus all the data except the
	// frames, skin and group info
	s32 size = sizeof(aliashdr_t) +
		(LittleLong(pinmodel->numframes)-1) * sizeof(pheader->frames[0])
		+ sizeof(mdl_t) + LittleLong(pinmodel->numverts)
		* sizeof(stvert_t)
		+ LittleLong(pinmodel->numtris) * sizeof(mtriangle_t) ;
	pheader = (aliashdr_t *)Hunk_AllocName(size, loadname);
	mdl_t *pmodel = (mdl_t*)((u8*)&pheader[1]+(LittleLong(
			pinmodel->numframes)-1)*sizeof(pheader->frames[0]));
	pheader->flags = mod->flags = LittleLong(pinmodel->flags);
	// endian-adjust and copy the data, starting with the alias model header
	pmodel->boundingradius = pheader->boundingradius
					= LittleFloat(pinmodel->boundingradius);
	pmodel->numskins = pheader->numskins = LittleLong(pinmodel->numskins);
	pmodel->skinwidth = pheader->skinwidth =LittleLong(pinmodel->skinwidth);
	pmodel->skinheight=pheader->skinheight=LittleLong(pinmodel->skinheight);
	if(pheader->skinheight > MAX_LBM_HEIGHT) Host_Error(
		"Mod_LoadAliasModel: model %s has a skin taller than %d",
			mod->name, MAX_LBM_HEIGHT);
	pmodel->numverts = pheader->numverts = LittleLong(pinmodel->numverts);
	if(pheader->numverts <= 0) Host_Error(
	    "Mod_LoadAliasModel: model %s has no vertices", mod->name);
	if(pheader->numverts > MAXALIASVERTS) Host_Error(
	    "Mod_LoadAliasModel: model %s has too many vertices(%d, max = %d)",
	    mod->name, pheader->numverts, MAXALIASVERTS);
	if(pheader->numtris > MAXALIASTRIS) Host_Error(
	    "Mod_LoadAliasModel: model %s has too many triangles(%d, max = %d)",
	    mod->name, pheader->numtris, MAXALIASTRIS);
	if(pmodel->numverts > MAXALIASVERTS) Host_Error(
	    "Mod_LoadAliasModel: model %s has too many vertexes(%d, max = %d)",
	    mod->name, pmodel->numverts, MAXALIASVERTS);
	    pmodel->numtris = pheader->numtris = LittleLong(pinmodel->numtris);
	if(pheader->numtris <= 0) Host_Error(
		"Mod_LoadAliasModel: model %s has no triangles", mod->name);
	if(pheader->numtris > MAXALIASTRIS) Host_Error(
		"model %s has too many triangles(%d > %d)",
		mod->name, pheader->numtris, MAXALIASTRIS);
	pmodel->numframes = pheader->numframes =LittleLong(pinmodel->numframes);
	s32 numframes = pheader->numframes;
	if(numframes < 1)
	    Host_Error("Mod_LoadAliasModel: Invalid # of frames: %d",numframes);
	pmodel->size = pheader->size = LittleFloat(pinmodel->size)
						* ALIAS_BASE_SIZE_RATIO;
	mod->synctype = (synctype_t)LittleLong(pinmodel->synctype);
	mod->numframes = pheader->numframes;
	for(s32 i = 0; i < 3; i++) {
		pmodel->scale[i] = pheader->scale[i] =
			LittleFloat(pinmodel->scale[i]);
		pmodel->scale_origin[i] = pheader->scale_origin[i] =
			LittleFloat(pinmodel->scale_origin[i]);
		pmodel->eyeposition[i] = pheader->eyeposition[i] =
			LittleFloat(pinmodel->eyeposition[i]);
	}
	s32 numskins = pmodel->numskins;
	if(pmodel->skinwidth & 0x03)
	    Con_DPrintf("Mod_LoadAliasModel: skinwidth not multiple of 4\n");
	pheader->model = (u8 *)pmodel - (u8 *)pheader;
	s32 skinsize = pmodel->skinheight * pmodel->skinwidth; // load the skins
	if(numskins < 1)
	     Host_Error("Mod_LoadAliasModel: Invalid # of skins: %d", numskins);
	pskintype = (daliasskintype_t *)&pinmodel[1];
	pskindesc = Hunk_AllocName(numskins*sizeof(maliasskindesc_t), loadname);
	pheader->skindesc = (u8 *)pskindesc - (u8 *)pheader;
	for(s32 i = 0; i < numskins; i++) {
		aliasskintype_t skintype;
		skintype = LittleLong(pskintype->type);
		pskindesc[i].type = skintype;
		if(skintype == ALIAS_SKIN_SINGLE) {
			pskintype = (daliasskintype_t *) Mod_LoadAliasSkin(
			  pskintype + 1, &pskindesc[i].skin, skinsize, pheader);
		} else {
			pskintype = (daliasskintype_t *) Mod_LoadAliasSkinGroup(
			  pskintype + 1, &pskindesc[i].skin, skinsize, pheader);
		}
	}
	pstverts = (stvert_t *)&pmodel[1]; // set base s and t vertices
	pinstverts = (stvert_t *)pskintype;
	pheader->stverts = (u8 *)pstverts - (u8 *)pheader;
	for(s32 i = 0; i < pheader->numverts; i++) {
		pstverts[i].onseam = LittleLong(pinstverts[i].onseam);
		// put s and t in 16.16 format
		pstverts[i].s = LittleLong(pinstverts[i].s) << 16;
		pstverts[i].t = LittleLong(pinstverts[i].t) << 16;
	}
	ptri = (mtriangle_t *)&pstverts[pmodel->numverts];//set up the triangles
	pintriangles = (dtriangle_t *)&pinstverts[pmodel->numverts];
	pheader->triangles = (u8 *)ptri - (u8 *)pheader;
	for(s32 i = 0; i<pheader->numtris ; i++) {
		ptri[i].facesfront = LittleLong(pintriangles[i].facesfront);
		for(s32 j = 0; j < 3; j++)
		  ptri[i].vertindex[j]=LittleLong(pintriangles[i].vertindex[j]);
	}
	posenum = 0; // load the frames
	pframetype = (daliasframetype_t *)&pintriangles[pmodel->numtris];
	for(s32 i = 0; i < numframes; i++) {
		aliasframetype_t frametype;
		frametype = (aliasframetype_t)LittleLong(pframetype->type);
		pheader->frames[i].type = frametype;
		if(frametype == ALIAS_SINGLE)
			pframetype = (daliasframetype_t *) Mod_LoadAliasFrame(
				pframetype + 1, &pheader->frames[i].frame,
				pmodel->numverts, &pheader->frames[i].bboxmin,
				&pheader->frames[i].bboxmax, pheader,
				pheader->frames[i].name, &pheader->frames[i],0);
		else pframetype = (daliasframetype_t *) Mod_LoadAliasGroup(
				pframetype + 1, &pheader->frames[i].frame,
				pmodel->numverts, &pheader->frames[i].bboxmin,
				&pheader->frames[i].bboxmax, pheader,
				pheader->frames[i].name, &pheader->frames[i]);
	}
	pheader->numposes = posenum; // Baker: Watch this and compare vs. Mark V
	mod->type = mod_alias;
	Mod_SetExtraFlags(mod); //johnfitz
	Mod_CalcAliasBounds(pheader); //johnfitz
	// move the complete, relocatable alias model to the cache
	s32 end = Hunk_LowMark();
	s32 total = end - start;
	Cache_Alloc(&mod->cache, total, loadname);
	if(!mod->cache.data) return;
	memcpy(mod->cache.data, pheader, total);
	Hunk_FreeToLowMark(start);
}

void *Mod_LoadSpriteFrame(void *pin, mspriteframe_t **ppframe)
{
	s32 width, height, size, origin[2];
	dspriteframe_t *pinframe = (dspriteframe_t *) pin;
	width = LittleLong(pinframe->width);
	height = LittleLong(pinframe->height);
	size = width * height;
	mspriteframe_t *pspriteframe = Hunk_AllocName(sizeof(mspriteframe_t) +
			size, loadname);
	Q_memset(pspriteframe, 0, sizeof(mspriteframe_t) + size);
	*ppframe = pspriteframe;
	pspriteframe->width = width;
	pspriteframe->height = height;
	origin[0] = LittleLong(pinframe->origin[0]);
	origin[1] = LittleLong(pinframe->origin[1]);
	pspriteframe->up = origin[1];
	pspriteframe->down = origin[1] - height;
	pspriteframe->left = origin[0];
	pspriteframe->right = width + origin[0];
	Q_memcpy(&pspriteframe->pixels[0], (u8 *) (pinframe+1), size);
	return(void *)((u8 *) pinframe + sizeof(dspriteframe_t) + size);
}

void *Mod_LoadSpriteGroup(void *pin, mspriteframe_t **ppframe)
{
	dspritegroup_t *pingroup = (dspritegroup_t *) pin;
	s32 numframes = LittleLong(pingroup->numframes);
	mspritegroup_t *pspritegroup = Hunk_AllocName(sizeof(mspritegroup_t) +
		(numframes - 1) * sizeof(pspritegroup->frames[0]), loadname);
	pspritegroup->numframes = numframes;
	*ppframe = (mspriteframe_t *) pspritegroup;
	dspriteinterval_t *pin_intervals = (dspriteinterval_t *) (pingroup + 1);
	f32 *poutintervals = Hunk_AllocName(numframes*sizeof(f32),loadname);
	pspritegroup->intervals = poutintervals;
	for(s32 i = 0; i < numframes; i++){
		*poutintervals = LittleFloat(pin_intervals->interval);
		if(*poutintervals <= 0.0)
			Sys_Error("Mod_LoadSpriteGroup: interval<=0");
		poutintervals++;
		pin_intervals++;
	}
	void *ptemp = (void *)pin_intervals;
	for(s32 i = 0; i < numframes; i++)
		ptemp = Mod_LoadSpriteFrame(ptemp, &pspritegroup->frames[i]);
	return ptemp;
}

void Mod_LoadSpriteModel(model_t *mod, void *buffer)
{
	dspriteframetype_t *pframetype;
	dsprite_t *pin = (dsprite_t *) buffer;
	s32 version = LittleLong(pin->version);
	if(version != SPRITE_VERSION)
		Sys_Error("%s has wrong version number " "(%i should be %i)",
				mod->name, version, SPRITE_VERSION);
	s32 numframes = LittleLong(pin->numframes);
	msprite_t *psprite;
	s32 size = sizeof(msprite_t) + (numframes - 1)*sizeof(psprite->frames);
	psprite = Hunk_AllocName(size, loadname);
	mod->cache.data = psprite;
	psprite->type = LittleLong(pin->type);
	psprite->maxwidth = LittleLong(pin->width);
	psprite->maxheight = LittleLong(pin->height);
	psprite->beamlength = LittleFloat(pin->beamlength);
	mod->synctype = LittleLong(pin->synctype);
	psprite->numframes = numframes;
	mod->mins[0] = mod->mins[1] = -psprite->maxwidth / 2;
	mod->maxs[0] = mod->maxs[1] = psprite->maxwidth / 2;
	mod->mins[2] = -psprite->maxheight / 2;
	mod->maxs[2] = psprite->maxheight / 2;
	// load the frames
	if(numframes < 1)
		Sys_Error("Mod_LoadSpriteModel: Invalid # of frames: %d\n",
				numframes);
	mod->numframes = numframes;
	mod->flags = 0;
	pframetype = (dspriteframetype_t *) (pin + 1);
	for(s32 i = 0; i < numframes; i++){
		spriteframetype_t frametype;
		frametype = LittleLong(pframetype->type);
		psprite->frames[i].type = frametype;
		if(frametype == SPR_SINGLE){
			pframetype = (dspriteframetype_t *)
				Mod_LoadSpriteFrame(pframetype + 1,
						&psprite->frames[i].frameptr);
		} else {
			pframetype = (dspriteframetype_t *)
				Mod_LoadSpriteGroup(pframetype + 1,
						&psprite->frames[i].frameptr);
		}
	}
	mod->type = mod_sprite;
}

static void Mod_Print()
{
	Con_SafePrintf("Cached models:\n");
	s32 i = 0;
	model_t *mod = mod_known;
	for(; i < mod_numknown; i++, mod++)
		Con_SafePrintf("%8p : %s\n", mod->cache.data, mod->name);
	Con_Printf("%i models\n",mod_numknown);
}
