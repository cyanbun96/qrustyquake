// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2007-2008 Kristian Duske
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
#include "quakedef.h"

static msurface_t *r_skyfaces;
static mplane_t r_skyplanes[6]; // Manoel Kasimier - edited
static mtexinfo_t r_skytexinfo[6];
static mvertex_t *r_skyverts;
static u8 rgb_lut[RGB_LUT_SIZE];
static s32 rgb_lut_built = 0;
static s32 skybox_planes[12] = {2,-128, 0,-128, 2,128, 1,128, 0,128, 1,-128};
static vec3_t box_vecs[6][2] = { {{0,-1,0}, {-1,0,0} }, { {0,1,0}, {0,0,-1} },
				{ {0,-1,0}, {1,0,0} },  { {1,0,0}, {0,0,-1} },
				{ {0,-1,0}, {0,0,-1} }, { {-1,0,0}, {0,0,-1}}};
static vec3_t box_bigvecs[6][2] = { {{0,-2,0}, {-2,0,0} }, {{0,2,0}, {0,0,-2} },
				   { {0,-2,0}, {2,0,0} }, { {2,0,0}, {0,0,-2} },
				   { {0,-2,0}, {0,0,-2}}, {{-2,0,0}, {0,0,-2}}};
static vec3_t box_bigbigvecs[6][2] = {{{0,-4,0},{-4,0,0} },{{0,4,0}, {0,0,-4}},
				     { {0,-4,0},{4,0,0} }, {{4,0,0}, {0,0,-4}},
				     { {0,-4,0},{0,0,-4} },{{-4,0,0},{0,0,-4}}};

s32 R_LoadSkybox(const c8 *name);

void Sky_LoadSkyBox(const c8 *name)
{
	if(r_enableskybox.value == 0) return;
	if(strcmp(skybox_name, name) == 0)
		return; //no change
	if(name[0] == 0){ // turn off skybox if sky is set to ""
		skybox_name[0] = 0;
		return;
	}
	static c8 last_skybox_name[1024];
	// Baker: If name matches, we already have the pixels loaded and we
	// don't actually need to reload
	if(strcmp(name, last_skybox_name)){
		if(cl.worldmodel && !R_LoadSkybox(name)){
			skybox_name[0] = 0;
			return;
		}
	}
	SDL_strlcpy(skybox_name, name, sizeof(skybox_name));
	SDL_strlcpy(last_skybox_name, skybox_name, sizeof(last_skybox_name));
}

void R_BuildRGBLUT(){
	for(s32 r = 0; r < 256; r += 32){
		for(s32 g = 0; g < 256; g += 32){
			for(s32 b = 0; b < 256; b += 32){
				s32 i = ((r/32) << 6) | ((g/32) << 3) | (b/32);
				rgb_lut[i] = rgbtoi(r, g, b);
			}
		}
	}
	rgb_lut_built = 1;
}

u8 fast_rgbtoi(u8 r, u8 g, u8 b)
{
	s32 i = ((r / 32) << 6) | ((g / 32) << 3) | (b / 32);
	return rgb_lut[i];
}

u8 *Upscale_NearestNeighbor(u8 *src, s32 width, s32 height, s32 *new_width,
	s32 *new_height){
	*new_width = (width < 256) ? 256 : width;
	*new_height = (height < 256) ? 256 : height;
	if(*new_width == width && *new_height == height)
		return src;
	u8 *dest = malloc((*new_width) * (*new_height));
	if(!dest) return 0;
	for(s32 y = 0; y < *new_height; ++y){ // sizes that don't divide 256
		s32 sy = y * height / *new_height; // need a fractional step
		for(s32 x = 0; x < *new_width; ++x){
			s32 sx = x * width / *new_width;
			dest[y*(*new_width)+x] = src[sy*width+sx];
		}
	}
	return dest;
}

u8 *Resize_Nearest(const u8 *src, s32 sw, s32 sh, s32 dw, s32 dh)
{
	u8 *dst = malloc(dw * dh);
	if(!dst) return NULL;
	for(s32 y = 0; y < dh; y++){
		s32 sy = (y * sh) / dh;
		for(s32 x = 0; x < dw; x++){
			s32 sx = (x * sw) / dw;
			dst[y * dw + x] = src[sy * sw + sx];
		}
	}
	return dst;
}

s32 R_LoadSkybox(const c8 *name)
{
	if(!name || !name[0]){
		skybox_name[0] = 0;
		return 0;
	}
	if(!strcmp(name, skybox_name)) // the same skybox we are using now
		return 1;
	SDL_strlcpy(skybox_name, name, sizeof(skybox_name));
	s32 mark = Hunk_LowMark();
	for(s32 i = 0 ; i < 6 ; i++){
		c8 pathname[1024];
		c8 *suf[6] = {"rt", "bk", "lf", "ft", "up", "dn"};
		s32 r_skysideimage[6] = {5, 2, 4, 1, 0, 3};
		snprintf(pathname, sizeof(pathname), "gfx/env/%s%s", name,
				suf[r_skysideimage[i]]);
		s32 width, height;
		u8 *pic = Image_LoadImage(pathname, &width, &height);
		if(!pic){
			Con_Printf("Couldn't load %s", pathname);
			return 0;
		}
		u8 *pdest = pic; // CyanBun96: palettize in place
		u8 *psrc = pic; // with some noise dithering
		if(!rgb_lut_built)
			R_BuildRGBLUT();
		for(s32 j = 0; j < width*height; j++, psrc+=4){
			s32 r = psrc[0] + ((f32)(((s32)lfsr_random() % 16) - 7)
					* r_skynoise.value);
			s32 g = psrc[1] + ((f32)(((s32)lfsr_random() % 16) - 7)
					* r_skynoise.value);
			s32 b = psrc[2] + ((f32)(((s32)lfsr_random() % 16) - 7)
					* r_skynoise.value);
			pdest[j] = fast_rgbtoi(CLAMP(0,r,255),
						CLAMP(0,g,255),
						CLAMP(0,b,255));
		}
		u8 *final_pic = Upscale_NearestNeighbor(pic,
				width, height, &width, &height);
		u8 *origpic = pic;
		pic = final_pic;
		s32 target = 1024;
		s32 m = width > height ? width : height;
		bool freescaled = false;
		if(m <= 512) target = 512;
		if(m <= 256) target = 256;
		u8 *scaled;
		if(width != target || height != target){
			scaled = Resize_Nearest(pic,width,height,target,target);
			freescaled = true;
			if(!scaled){
			    Con_Printf("skybox resize failed for %s", pathname);
				Hunk_FreeToLowMark(mark);
				return 0;
			}
			width = height = target;
			pic = scaled;
		}
		switch(width){ // Manoel Kasimier - hi-res skyboxes - begin
			case 1024: // falls through
			case 512: // falls through
			case 256: // We're good!
				break;
			default: // We aren't good
Con_Printf("skybox width(%d) for %s must be 256, 512, 1024", width, pathname);
				pic = origpic;
				Hunk_FreeToLowMark(mark);
				if(final_pic != origpic)
					free(final_pic);
				return 0;
		}
		switch(height){
			case 1024: // falls through
			case 512: // falls through
			case 256: // We're good!
				break;
			default:
Con_Printf("skybox height(%d) for %s must be 256, 512, 1024", height, pathname);
				pic = origpic;
				Hunk_FreeToLowMark(mark);
				if(final_pic != origpic)
					free(final_pic);
				return 0;
		}
		static texture_t r_skytextures[6];
		r_skytexinfo[i].texture = &r_skytextures[i];
		r_skytexinfo[i].texture->width = width;
		r_skytexinfo[i].texture->height = height;
		r_skytexinfo[i].texture->offsets[0] = i;
		switch(width){
		case 1024:
		    VectorCopy(box_bigbigvecs[i][0], r_skytexinfo[i].vecs[0]);
		    break;
		case 512:
			VectorCopy(box_bigvecs[i][0], r_skytexinfo[i].vecs[0]);
			break;
		default:
			VectorCopy(box_vecs[i][0], r_skytexinfo[i].vecs[0]);
			break;
		}
		switch(height){
		case 1024:
		    VectorCopy(box_bigbigvecs[i][1], r_skytexinfo[i].vecs[1]);
		    break;
		case 512:
			VectorCopy(box_bigvecs[i][1], r_skytexinfo[i].vecs[1]);
			break;
		default:
			VectorCopy(box_vecs[i][1], r_skytexinfo[i].vecs[1]);
			break;
		}
		// This is if one is already loaded and the size changed
		if(r_skyfaces){
			r_skyfaces[i].texturemins[0] = -(width/2);
			r_skyfaces[i].texturemins[1] = -(height/2);
			r_skyfaces[i].extents[0] = width;
			r_skyfaces[i].extents[1] = height;
		}
		else
	Con_DPrintf("Warning: No surface to load yet for WinQuake skybox");
		memcpy(r_skypixels[i], pic, width*height);
		pic = origpic;
		if(freescaled) free(scaled);
		Hunk_FreeToLowMark(mark);
		if(final_pic != origpic)
			free(final_pic);
	}
	return 1;
}

// Manoel Kasimier - skyboxes - begin
void R_InitSkyBox()
{
	model_t *loadmodel = cl.worldmodel; // Manoel Kasimier - edited
	r_skyfaces = loadmodel->surfaces + loadmodel->numsurfaces;
	loadmodel->numsurfaces += 6;
	r_skyverts = loadmodel->vertexes + loadmodel->numvertexes;
	loadmodel->numvertexes += 8;
	static medge_t *r_skyedges;
	r_skyedges = loadmodel->edges + loadmodel->numedges;
	loadmodel->numedges += 12;
	static s32 *r_skysurfedges;
	r_skysurfedges = loadmodel->surfedges + loadmodel->numsurfedges;
	loadmodel->numsurfedges += 24;
	memset(r_skyfaces, 0, 6 * sizeof(*r_skyfaces));
	for(s32 i = 0; i < 6; i++){
		r_skyplanes[i].normal[skybox_planes[i*2]] = 1;
		r_skyplanes[i].dist = skybox_planes[i*2+1];
		r_skyfaces[i].plane = &r_skyplanes[i];
		r_skyfaces[i].numedges = 4;
		static s32 box_faces[6] = {0,0,2,2,2,0};
		r_skyfaces[i].flags = box_faces[i] | SURF_DRAWSKYBOX;
		r_skyfaces[i].firstedge = loadmodel->numsurfedges-24+i*4;
		r_skyfaces[i].texinfo = &r_skytexinfo[i];
		// Manoel Kasimier - hi-res skyboxes - begin
		s32 width, height;
		if(r_skytexinfo[i].texture){
			width = r_skytexinfo[i].texture->width;
			height = r_skytexinfo[i].texture->height;
		}
		else width = height = 256;
		switch(width){
		case 1024:
		    VectorCopy(box_bigbigvecs[i][0], r_skytexinfo[i].vecs[0]);
		    break;
		case 512:
			VectorCopy(box_bigvecs[i][0], r_skytexinfo[i].vecs[0]);
			break;
		default:
			VectorCopy(box_vecs[i][0], r_skytexinfo[i].vecs[0]);
			break;
		}
		switch(height){
		case 1024:
		    VectorCopy(box_bigbigvecs[i][1], r_skytexinfo[i].vecs[1]);
		    break;
		case 512:
			VectorCopy(box_bigvecs[i][1], r_skytexinfo[i].vecs[1]);
			break;
		default:
			VectorCopy(box_vecs[i][1], r_skytexinfo[i].vecs[1]);
			break;
		}
		r_skyfaces[i].texturemins[0] = -(width/2);
		r_skyfaces[i].texturemins[1] = -(height/2);
		r_skyfaces[i].extents[0] = width;
		r_skyfaces[i].extents[1] = height;
		// Manoel Kasimier - hi-res skyboxes - end
	}
	static s32 box_surfedges[24] = {
		1,2,3,4,  -1,5,6,7,  8,9,-6,10,  -2,-7,-9,11,
		12,-3,-11,-8,     -12,-10,-5,-4};
	for(s32 i = 0; i < 24; i++)
		if(box_surfedges[i] > 0)
			r_skysurfedges[i] = loadmodel->numedges-13
						+ box_surfedges[i];
		else
			r_skysurfedges[i] = - (loadmodel->numedges-13
						+ -box_surfedges[i]);
	for(s32 i = 0; i < 12; i++){
		static s32 box_edges[24] =
			{1,2,2,3,3,4,4,1,1,5,5,6,6,2,7,8,8,6,5,7,8,3,7,4};
		r_skyedges[i].v[0] = loadmodel->numvertexes-9+box_edges[i*2+0];
		r_skyedges[i].v[1] = loadmodel->numvertexes-9+box_edges[i*2+1];
		r_skyedges[i].cachededgeoffset = 0;
	}
}

void R_EmitSkyBox()
{
	static f32 box_verts[8][3] = {
		{-1,-1,-1}, {-1,1,-1}, {1,1,-1}, {1,-1,-1},
		{-1,-1,1}, {-1,1,1}, {1,-1,1}, {1,1,1} };
	for(s32 i = 0 ; i < 8 ; i++) // set the eight fake vertexes
		for(s32 j = 0 ; j < 3 ; j++)
		  r_skyverts[i].position[j] = r_origin[j] + box_verts[i][j]*128;
	for(s32 i = 0 ; i < 6 ; i++) // set the six fake planes
		if(skybox_planes[i*2+1] > 0)
			r_skyplanes[i].dist = r_origin[skybox_planes[i*2]]+128;
		else
			r_skyplanes[i].dist = r_origin[skybox_planes[i*2]]-128;
	for(s32 i = 0 ; i < 6; i++){ // fix texture offsets
		r_skytexinfo[i].vecs[0][3] =
			-DotProduct(r_origin, r_skytexinfo[i].vecs[0]);
		r_skytexinfo[i].vecs[1][3] =
			-DotProduct(r_origin, r_skytexinfo[i].vecs[1]);
	}
	s32 oldkey = r_currentkey; // emit the six faces
	r_currentkey = 0x7ffffff0;
	skybox_surf_p = surface_p;
	for(s32 i = 0; i < 6; i++)
		R_RenderFace(r_skyfaces + i, 15);
	r_currentkey = oldkey; // bsp sorting order
} // Manoel Kasimier - skyboxes - end

void Sky_NewMap()
{ // read worldspawn(this is so ugly, and shouldn't it be done on the server?)
	c8 key[128], value[4096];
	const c8 *data = cl.worldmodel->entities;
	data = COM_Parse(data);
	if(!data || com_token[0] != '{') // should never happen
		return; // error
	while(1){
		data = COM_Parse(data);
		if(!data)
			return; // error
		if(com_token[0] == '}')
			break; // end of worldspawn
		if(com_token[0] == '_')
			SDL_strlcpy(key, com_token + 1, sizeof(key));
		else
			SDL_strlcpy(key, com_token, sizeof(key));
		while(key[0] && key[strlen(key)-1] == ' ')
			key[strlen(key)-1] = 0; // remove trailing spaces
		data = COM_ParseEx(data, CPE_ALLOWTRUNC);
		if(!data)
			return; // error
		SDL_strlcpy(value, com_token, sizeof(value));
		if(!strcmp("sky", key))
			Sky_LoadSkyBox(value);
		if(!strcmp("skyfog", key))
			Cvar_SetValue("r_skyfog", atof(value));
		else if(!strcmp("skyname", key)) // half-life
			Sky_LoadSkyBox(value);
		else if(!strcmp("qlsky", key)) // quake lives
			Sky_LoadSkyBox(value);
	}
}

void Sky_SkyCommand_f()
{
	switch(Cmd_Argc()){
	case 1:
		Con_Printf("\"sky\" is \"%s\"\n", skybox_name);
		break;
	case 2:
		Sky_LoadSkyBox(Cmd_Argv(1));
		break;
	default:
		Con_Printf("usage: sky <skyname>\n");
	}
}

void Sky_Init()
{
	Cmd_AddCommand("sky", Sky_SkyCommand_f);
	Cvar_RegisterVariable(&r_enableskybox);
	Cvar_RegisterVariable(&r_skyfog);
	Cvar_RegisterVariable(&r_skynoise);
	skybox_name[0] = 0;
}

s32 R_SkyIndexForTexture(texture_t *mt)
{
	for(s32 i = 0; i < r_numskies; ++i)
		if(!strncmp(mt->name, r_skies[i].name, 16))
			return i;
	Sys_Error("Sky not found: %s\n", mt->name);
	return 0;
}

void R_InitSky(texture_t *mt)
{ // A sky texture is 256*128, with the right side being a masked overlay
	static s32 maxskies;
	s32 newi = 0;
	while(newi < r_numskies && strncmp(r_skies[newi].name, mt->name, 16))
		newi++; // a reloaded model reuses its slot
	if(newi == r_numskies){
		if(r_numskies == maxskies){
			maxskies += 8;
			r_skies = realloc(r_skies, maxskies*sizeof(*r_skies));
			if(!r_skies) Sys_Error(
				"R_InitSky: realloc() failed on %d skies",
				maxskies);
		}
		r_numskies++;
		strncpy(r_skies[newi].name, mt->name, 16);
	}
	sky_t *sky = &r_skies[newi];
	sky->made = 0;
	sky->xlast = sky->ylast = -1;
	// skies of any other size get sampled down or up to it
	u8 *src = (u8 *) mt + mt->offsets[0];
	s32 w = mt->width, h = mt->height, half = w / 2;
	if(half < 1){
		Con_Printf("Sky %.16s is too small (%dx%d)\n", mt->name, w, h);
		memset(sky->source, 0, sizeof(sky->source));
		memset(sky->bottomsky, 0, sizeof(sky->bottomsky));
		memset(sky->bottommask, 0xff, sizeof(sky->bottommask));
		return;
	}
	// row i, column j of a 128x128 layer; the left half is the front
#define SKYSRC(i, j, left) src[(i) * h / 128 * w + (left) + (j) * half / 128]
	for(s32 i = 0; i < 128; i++)
		for(s32 j = 0; j < 128; j++)
			sky->source[(i*256)+j+128] = SKYSRC(i, j, half);
	for(s32 i = 0; i < 128; i++)
		for(s32 j = 0; j < 131; j++){
			if(SKYSRC(i, j & 0x7F, 0)){
				sky->bottomsky[(i * 131) + j] =
				    SKYSRC(i, j & 0x7F, 0);
				sky->bottommask[(i * 131) + j] = 0;
			}else{
				sky->bottomsky[(i * 131) + j] = 0;
				sky->bottommask[(i * 131) + j] = 0xff;
			}
		}
#undef SKYSRC
}

void R_MakeSky(texture_t *mt)
{
	sky_t *sky = &r_skies[R_SkyIndexForTexture(mt)];
	s32 xshift = skytime * skyspeed;
	s32 yshift = skytime * skyspeed;
	if((xshift == sky->xlast) && (yshift == sky->ylast))
		return;
	sky->xlast = xshift;
	sky->ylast = yshift;
	u32 *pnewsky = (u32 *)sky->source;
	for(s32 y = 0; y < SKYSIZE; y++){
		s32 baseofs = ((y + yshift) & SKYMASK) * 131;
		for(s32 x = 0; x < SKYSIZE; x++){
			s32 ofs = baseofs + ((x + xshift) & SKYMASK);
			*(u8 *) pnewsky = (*((u8 *) pnewsky + 128) &
					sky->bottommask[ofs]) |
					sky->bottomsky[ofs];
			pnewsky = (unsigned *)((u8 *) pnewsky + 1);
		}
		pnewsky += 128 / sizeof(unsigned);
	}
	sky->made = 1;
}

void R_SetSkyFrame()
{
	static s32 iskyspeed = 8;
	skyspeed = iskyspeed;
	static s32 iskyspeed2 = 2;
	s32 g = GreatestCommonDivisor(iskyspeed, iskyspeed2);
	s32 s1 = iskyspeed / g;
	s32 s2 = iskyspeed2 / g;
	f32 temp = SKYSIZE * s1 * s2;
	skytime = cl.time - ((s32)(cl.time / temp) * temp);
	for(s32 skyi = 0; skyi < r_numskies; ++skyi)
		r_skies[skyi].made = 0;
}
