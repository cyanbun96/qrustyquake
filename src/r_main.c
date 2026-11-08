// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"
#include <fenv.h>

#ifndef FE_TOWARDZERO // not every target can change it
#define FE_TOWARDZERO FE_TONEAREST
#endif

typedef struct { s32 x0, y0, x1, y1; } debugline_t;
typedef struct { vec3_t origin; } debugpoint_t;
static debugline_t *r_debuglines = 0;
static debugpoint_t *r_debugpoints = 0;
static vec3_t viewlightvec;
static alight_t r_viewlighting = { 128, 192, viewlightvec };
static f32 aliaslightvec[3] = { -1, 0, 0 };

static void R_AddDebugLine(f32 x0, f32 y0, f32 x1, f32 y1)
{
	if(r_numdebuglines >= MAX_DEBUG_LINES)
		return;
	debugline_t *l = &r_debuglines[r_numdebuglines++];
	l->x0 = x0;
	l->y0 = y0;
	l->x1 = x1;
	l->y1 = y1;
}

bool R_ProjectPointToScreen(vec3_t world, s32 *screenX, s32 *screenY)
{ // Helper to project a 3D point(Used by point entities)
	vec3_t local, transformed;
	VectorSubtract(world, r_origin, local);
	TransformVector(local, transformed);
	if(transformed[2] < NEAR_CLIP)return false;
	f32 lzi = 1.0 / transformed[2];
	*screenX = (int)(xcenter + (xscale * lzi) * transformed[0]);
	*screenY = (int)(ycenter - (yscale * lzi) * transformed[1]);
	return true;
}

void R_DrawDebugLine3D(vec3_t p1, vec3_t p2)
{
	vec3_t t1, t2, local1, local2;
	VectorSubtract(p1, modelorg, local1);
	TransformVector(local1, t1);
	VectorSubtract(p2, modelorg, local2);
	TransformVector(local2, t2);
	if(t1[2] < NEAR_CLIP && t2[2] < NEAR_CLIP)return; // Near Plane Clipping
	if(t1[2] < NEAR_CLIP){
		f32 frac = (NEAR_CLIP - t1[2]) / (t2[2] - t1[2]);
		t1[0] += frac * (t2[0] - t1[0]);
		t1[1] += frac * (t2[1] - t1[1]);
		t1[2] = NEAR_CLIP;
	}else if(t2[2] < NEAR_CLIP){
		f32 frac = (NEAR_CLIP - t2[2]) / (t1[2] - t2[2]);
		t2[0] += frac * (t1[0] - t2[0]);
		t2[1] += frac * (t1[1] - t2[1]);
		t2[2] = NEAR_CLIP;
	}
	// 3D Frustum Boundary Clipping Math
	// We dynamically derive the 4 frustum planes from the projection matrix
	// A point is inside the screen if its calculated distance is >= 0.
	f32 w = vid.width;
	f32 h = vid.height;
	#define DIST_LEFT(t)   ( xscale * t[0] + xcenter * t[2])
	#define DIST_RIGHT(t)  (-xscale * t[0] + (w - xcenter) * t[2])
	#define DIST_TOP(t)    (-yscale * t[1] + ycenter * t[2])
	#define DIST_BOTTOM(t) ( yscale * t[1] + (h - ycenter) * t[2])
	f32 d1[4] = {DIST_LEFT(t1),DIST_RIGHT(t1),DIST_TOP(t1),DIST_BOTTOM(t1)};
	f32 d2[4] = {DIST_LEFT(t2),DIST_RIGHT(t2),DIST_TOP(t2),DIST_BOTTOM(t2)};
	for(s32 i = 0; i < 4; i++){
		if(d1[i]<0 && d2[i]<0)return; //Fully outside this frustum plane
		if(d1[i] < 0){
			f32 frac = d1[i] / (d1[i] - d2[i]);
			t1[0] += frac * (t2[0] - t1[0]);
			t1[1] += frac * (t2[1] - t1[1]);
			t1[2] += frac * (t2[2] - t1[2]);
			d1[0] = DIST_LEFT(t1); d1[1] = DIST_RIGHT(t1);
			d1[2] = DIST_TOP(t1);  d1[3] = DIST_BOTTOM(t1);
		}else if(d2[i] < 0){
			f32 frac = d2[i] / (d2[i] - d1[i]);
			t2[0] += frac * (t1[0] - t2[0]);
			t2[1] += frac * (t1[1] - t2[1]);
			t2[2] += frac * (t1[2] - t2[2]);
			d2[0] = DIST_LEFT(t2); d2[1] = DIST_RIGHT(t2);
			d2[2] = DIST_TOP(t2);  d2[3] = DIST_BOTTOM(t2);
		}
	}
	f32 lzi1 = 1.0 / t1[2]; // 2D Projection
	f32 lzi2 = 1.0 / t2[2];
	R_AddDebugLine(xcenter + (xscale * lzi1) * t1[0],
		ycenter - (yscale * lzi1) * t1[1],
		xcenter + (xscale * lzi2) * t2[0],
		ycenter - (yscale * lzi2) * t2[1]);
}

static void R_DebugDrawBBox(vec3_t origin, vec3_t mins, vec3_t maxs)
{
	vec3_t old_modelorg;//save current state,force modelorg to global spaces
	VectorCopy(modelorg, old_modelorg);
	VectorCopy(r_origin, modelorg);
	vec3_t corners[8];
	for(s32 i = 0; i < 8; i++)
		for(s32 k = 0; k < 3; k++)
			corners[i][k] = origin[k] +
				((i >> k) & 1 ? maxs[k] : mins[k]);
	for(s32 i = 0; i < 8; i++) // the edges along each axis
		for(s32 k = 1; k < 8; k <<= 1)
			if(!(i & k))
				R_DrawDebugLine3D(corners[i], corners[i | k]);
	VectorCopy(old_modelorg, modelorg); // restore prev modelorg state
}

static s32 ComputeOutCode(s32 x, s32 y, s32 w, s32 h)
{
	s32 code = 0;
	if(x < 0)code |= CLIP_LEFT;
	else if(x >= w) code |= CLIP_RIGHT;
	if(y < 0)code |= CLIP_TOP;
	else if(y >= h)code |= CLIP_BOTTOM;
	return code;
}

void R_DrawDebugLine(s32 x0, s32 y0, s32 x1, s32 y1, u8 color)
{ // Highly optimized Bresenham algorithm with Cohen-Sutherland pre-clipping
	s32 w = vid.width;
	s32 h = vid.height;
	s32 outcode0 = ComputeOutCode(x0, y0, w, h); // clipping
	s32 outcode1 = ComputeOutCode(x1, y1, w, h);
	bool accept = false;
	while(1){
		if(!(outcode0 | outcode1)){ // completely inside bounds
			accept = true;
			break;
		}else if(outcode0 & outcode1){ // completely outside bounds
			break;
		}else{
			// partially inside, calculate the intersection point
			s32 x = 0, y = 0;
			s32 outcodeOut = outcode0 ? outcode0 : outcode1;
			if(outcodeOut & CLIP_BOTTOM){
			    x=x0+(s32)((((s64)x1-x0)*(h-1-y0))/((s64)y1-y0));
				y = h - 1;
			}else if(outcodeOut & CLIP_TOP){
				x=x0+(s32)((((s64)x1-x0)*(0-y0))/((s64)y1-y0));
				y = 0;
			}else if(outcodeOut & CLIP_RIGHT){
			    y=y0+(s32)((((s64)y1-y0)*(w-1-x0))/((s64)x1-x0));
				x = w - 1;
			}else if(outcodeOut & CLIP_LEFT){
				y=y0+(s32)((((s64)y1-y0)*(0-x0))/((s64)x1-x0));
				x = 0;
			}
			// move the outside point to the intersection point
			if(outcodeOut == outcode0){
				x0 = x; y0 = y;
				outcode0 = ComputeOutCode(x0, y0, w, h);
			}else{
				x1 = x; y1 = y;
				outcode1 = ComputeOutCode(x1, y1, w, h);
			}
		}
	}
	if(!accept)return; // completely off-screen
	s32 dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1; // rasterization
	s32 dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
	s32 err = dx + dy, e2;
	for(;;){
		// keep a minor bounds check here strictly as a failsafe against
		// integer division rounding errors from the clipping math above
		if(x0 >= 0 && x0 < w && y0 >= 0 && y0 < h)
			vid.buffer[y0 * w + x0] = color;
		if(x0 == x1 && y0 == y1)break;
		e2 = 2 * err;
		if(e2 >= dy){err += dy; x0 += sx;}
		if(e2 <= dx){err += dx; y0 += sy;}
	}
}

static void R_AllocShowTris()
{
	if(!r_debuglines)r_debuglines =
		malloc(MAX_DEBUG_LINES * sizeof(debugline_t));
	if(!r_debugpoints)r_debugpoints =
		malloc(MAX_DEBUG_POINTS * sizeof(debugpoint_t));
	if(!r_debuglines || !r_debugpoints)
		Sys_Error("Failed to allocate memory for r_showtris");
}

static void R_ParseDebugEntities()
{
	R_AllocShowTris();
	r_numdebugpoints = 0;
	if(!cl.worldmodel || !cl.worldmodel->entities) return;
	const c8 *data = cl.worldmodel->entities;
	c8 key[128], value[4096];
	while(1){
		data = COM_Parse(data);
		if(!data)break;
		if(com_token[0] != '{')continue;
		bool is_point = true;
		vec3_t parsed_origin = {0, 0, 0};
		while(1){
			data = COM_Parse(data);
			if(!data || com_token[0] == '}')break;
			strcpy(key, com_token);
			data = COM_ParseEx(data, CPE_ALLOWTRUNC);
			if(!data)break;
			strcpy(value, com_token);
			if(!strcmp(key, "origin"))
				sscanf(value, "%f %f %f", &parsed_origin[0],
					&parsed_origin[1], &parsed_origin[2]);
			else if(!strcmp(key, "model")) // Has a brush model
				is_point = false; // not a point entity
		}
		if((is_point && (parsed_origin[0] != 0 || parsed_origin[1] != 0
				|| parsed_origin[2] != 0))
				&& r_numdebugpoints < MAX_DEBUG_POINTS){
			VectorCopy(parsed_origin,
					r_debugpoints[r_numdebugpoints].origin);
			r_numdebugpoints++;
		}
	}
}

static void R_DebugDrawPoint(vec3_t origin)
{ // a cross
	s32 sx, sy;
	if(!R_ProjectPointToScreen(origin, &sx, &sy))
		return;
	R_AddDebugLine(sx - 4, sy, sx + 4, sy);
	R_AddDebugLine(sx, sy - 4, sx, sy + 4);
}

void R_InitTextures()
{ // create a simple checkerboard texture for the default
	r_notexture_mip = Hunk_AllocName(sizeof(texture_t)
			+ 16 * 16 + 8 * 8 + 4 * 4 + 2 * 2, "notexture");
	r_notexture_mip->width = r_notexture_mip->height = 16;
	r_notexture_mip->offsets[0] = sizeof(texture_t);
	r_notexture_mip->offsets[1] = r_notexture_mip->offsets[0] + 16 * 16;
	r_notexture_mip->offsets[2] = r_notexture_mip->offsets[1] + 8 * 8;
	r_notexture_mip->offsets[3] = r_notexture_mip->offsets[2] + 4 * 4;
	for(s32 m = 0; m < 4; m++){
		u8 *dest = (u8 *) r_notexture_mip
			+ r_notexture_mip->offsets[m];
		for(s32 y = 0; y < (16 >> m); y++)
			for(s32 x = 0; x < (16 >> m); x++){
				if((y < (8 >> m)) ^ (x < (8 >> m)))
					*dest++ = 0;
				else
					*dest++ = 0xff;
			}
	}
}

static void R_InitTurb()
{
	s32 amp = 8*0x10000;
	for(s32 i = 0; i < (SIN_BUFFER_SIZE); i++){
		sintable[i] = amp + sin(i * 3.14159 * 2 / TURB_CYCLE) * amp;
		intsintable[i]=TURB_AMP2+sin(i*3.14159*2/TURB_CYCLE)*TURB_AMP2;
	}
}

static void R_ViewChangedCallback(SDL_UNUSED cvar_t *cvar)
{ vid.recalc_refdef = 1; }

void R_Init()
{
	Cmd_AddCommand("timerefresh", R_TimeRefresh_f);
	Cmd_AddCommand("pointfile", R_ReadPointFile_f);
	Cmd_AddCommand("fog", Fog_FogCommand_f);
	Cmd_AddCommand("vid_scalemode", VID_SetScaleModeCommand_f);
	Cmd_AddCommand("vid_setmode", VID_VidSetModeCommand_f);
	Cvar_RegisterVariable(&r_draworder);
	Cvar_RegisterVariable(&r_speeds);
	Cvar_RegisterVariable(&r_drawflat);
	Cvar_RegisterVariable(&r_ambient);
	Cvar_RegisterVariable(&r_clearcolor);
	Cvar_RegisterVariable(&r_waterwarp);
	Cvar_RegisterVariable(&r_fullbright);
	Cvar_RegisterVariable(&r_drawentities);
	Cvar_RegisterVariable(&r_drawviewmodel);
	Cvar_RegisterVariable(&r_reportsurfout);
	Cvar_RegisterVariable(&r_numsurfs);
	Cvar_RegisterVariable(&r_reportedgeout);
	Cvar_RegisterVariable(&r_numedges);
	Cvar_RegisterVariable(&r_aliastransbase);
	Cvar_RegisterVariable(&r_aliastransadj);
	Cvar_RegisterVariable(&r_dithertex);
	Cvar_RegisterVariable(&r_wateralpha);
	Cvar_RegisterVariable(&r_slimealpha);
	Cvar_RegisterVariable(&r_lavaalpha);
	Cvar_RegisterVariable(&r_telealpha);
	Cvar_RegisterVariable(&r_fogstyle);
	Cvar_RegisterVariable(&r_nofog);
	Cvar_RegisterVariable(&r_alphastyle);
	Cvar_RegisterVariable(&r_entalpha);
	Cvar_RegisterVariable(&r_entocclusion);
	Cvar_RegisterVariable(&r_threads);
	Cvar_RegisterVariable(&r_bmodelocclusion);
	Cvar_RegisterVariable(&r_labmixpal);
	Cvar_RegisterVariable(&r_rgblighting);
	Cvar_RegisterVariable(&r_fogbrightness);
	Cvar_RegisterVariable(&r_fogfactor);
	Cvar_RegisterVariable(&r_fogscale);
	Cvar_RegisterVariable(&r_fognoise);
	Cvar_RegisterVariable(&r_lockfog);
	Cvar_RegisterVariable(&r_lockfogd);
	Cvar_RegisterVariable(&r_lockfogr);
	Cvar_RegisterVariable(&r_lockfogg);
	Cvar_RegisterVariable(&r_lockfogb);
	Cvar_RegisterVariable(&r_rebuildmips);
	Cvar_RegisterVariable(&r_fogdepthcorrection);
	Cvar_RegisterVariable(&r_fullbright_list);
	Cvar_RegisterVariable(&r_nolerp_list);
	Cvar_RegisterVariable(&r_litwater);
	Cvar_RegisterVariable(&r_novis);
	Cvar_RegisterVariable(&r_particlescale);
	Cvar_RegisterVariable(&r_particlestyle);
	Cvar_RegisterVariable(&r_particlesize);
	Cvar_RegisterVariable(&r_particlealpha);
	Cvar_RegisterVariable(&r_fovmode);
	Cvar_RegisterVariable(&r_hlwater);
	Cvar_RegisterVariable(&r_hlripplescale);
	Cvar_RegisterVariable(&r_hlwavescale);
	Cvar_RegisterVariable(&r_hlwaterquality);
	Cvar_RegisterVariable(&vid_cwidth);
	Cvar_RegisterVariable(&vid_cheight);
	Cvar_RegisterVariable(&vid_cwmode);
	Cvar_RegisterVariable(&scr_uixscale);
	Cvar_RegisterVariable(&scr_uiyscale);
	Cvar_RegisterVariable(&yaspectscale);
	Cvar_RegisterVariable(&scr_lockuiscale);
	Cvar_RegisterVariable(&r_mipscale);
	Cvar_RegisterVariable(&r_cutoutmipscale);
	Cvar_RegisterVariable(&scr_menubgstyle);
	Cvar_RegisterVariable(&lyr_main);
	Cvar_RegisterVariable(&lyr_sbar);
	Cvar_RegisterVariable(&lyr_menu);
	Cvar_RegisterVariable(&lyr_centerprint);
	Cvar_RegisterVariable(&lyr_console);
	Cvar_RegisterVariable(&lyr_notify);
	Cvar_RegisterVariable(&lyr_crosshair);
	Cvar_RegisterVariable(&r_renderscale);
	Cvar_RegisterVariable(&cl_gun_fovscale);
	Cvar_RegisterVariable(&r_showtris);
	Cvar_RegisterVariable(&r_showtris_color);
	Cvar_RegisterVariable(&r_coarseocclusion);
	Cvar_RegisterVariable(&r_lerpmodels);
	Cvar_RegisterVariable(&r_lerpmove);
	Cvar_RegisterVariable(&r_lerplightstyles);
	Cvar_SetCallback(&r_labmixpal, R_LabMixPalChanged);
	Cvar_SetCallback(&r_rgblighting, D_FlushCaches);
	Cvar_SetCallback(&r_fogbrightness, Fog_SetPalIndex);
	Cvar_SetCallback(&r_wateralpha, R_SetWateralpha_f);
	Cvar_SetCallback(&r_lavaalpha, R_SetLavaalpha_f);
	Cvar_SetCallback(&r_telealpha, R_SetTelealpha_f);
	Cvar_SetCallback(&r_slimealpha, R_SetSlimealpha_f);
	Cvar_SetCallback(&r_fovmode, R_ViewChangedCallback);
	Cvar_SetCallback(&r_renderscale, VID_SetRenderScaleCommand_f);
	Cvar_SetCallback(&r_ambient, D_FlushCaches);
	Cvar_SetCallback(&yaspectscale, R_ViewChangedCallback);
	view_clipplanes[0].leftedge = 1;
	view_clipplanes[1].rightedge = 1;
	view_clipplanes[1].leftedge = view_clipplanes[2].leftedge =
		view_clipplanes[3].leftedge = 0;
	view_clipplanes[0].rightedge = view_clipplanes[2].rightedge =
		view_clipplanes[3].rightedge = 0;
	r_refdef.xOrigin = 0.5;
	r_refdef.yOrigin = 0.5;
	R_InitTurb();
	R_InitParticles();
	D_Init();
	Sky_Init();
}

static void Pal_ParseWorldspawn()
{ // called at map load
	c8 key[128], value[4096];
	bool found = 0;
	const c8 *data = COM_Parse(cl.worldmodel->entities);
	if(!data || com_token[0] != '{'){ // error
		ResetWorldPal();
		return;
	}
	while(COM_ParseEntityPair(&data, key, value)){
		if(!strcmp("palette", key)){
			c8 pal[MAX_OSPATH] = "";
			c8 cmap[MAX_OSPATH];
			sscanf(value, "%259s", pal); // fits any MAX_OSPATH
			SDL_strlcpy(cmap, pal, MAX_OSPATH);
			SDL_strlcat(cmap, "_colormap", MAX_OSPATH);
			Con_DPrintf("Requested palettes %s %s\n", pal, cmap);
			found = SetWorldPal(pal, cmap);
		}
	}
	if(!found) ResetWorldPal();
	Fog_SetPalIndex(0);
}

void R_NewMap()
{
	skybox_name[0] = 0; // before a new palette can try to reload it
	lit_loaded = cl.worldmodel->haslit;
	Pal_ParseWorldspawn();
	R_InitSkyBox(); // Manoel Kasimier - skyboxes
	// clear out efrags in case the level hasn't been reloaded
	for(s32 i = 0; i < cl.worldmodel->numleafs; i++)
		cl.worldmodel->leafs[i].efrags = NULL;
	r_viewleaf = NULL;
	R_ClearParticles();
	r_maxedgesseen = 0;
	r_maxsurfsseen = 0;
	r_dowarpold = 0;
	r_viewchanged = 0;
	for(s32 i = 0; i < r_numskies; ++i)
		r_skies[i].made = 0;
	Sky_NewMap();
	Fog_ParseWorldspawn();
	R_ParseWorldspawn();
	R_ParseDebugEntities(); //init point-ent wireframe
}

void R_SetVrect(vrect_t *pvrectin, vrect_t *pvrect, s32 lineadj)
{
	f32 size = scr_viewsize.value > 100 ? 100 : scr_viewsize.value;
	if(cl.intermission){
		size = 100;
		lineadj = 0;
	}
	size /= 100;
	s32 h = pvrectin->height - lineadj;
	pvrect->width = pvrectin->width * size;
	if(pvrect->width < 96){
		size = 96.0 / pvrectin->width;
		pvrect->width = 96; // min for icons
	}
	if(pvrect->width < pvrectin->width) // a full view fills any screen
		pvrect->width &= ~7;
	pvrect->height = pvrectin->height * size;
	if(pvrect->height >= h)
		pvrect->height = h;
	else
		pvrect->height &= ~1;
	pvrect->x = (pvrectin->width - pvrect->width) / 2;
	pvrect->y = (h - pvrect->height) / 2;
}

void R_ViewChanged(vrect_t *pvrect, s32 lineadj, f32 aspect)
{ // Called every time the vid structure or r_refdef changes.
 // Guaranteed to be called before the first refresh
	if(r_fovmode.value == 1) aspect *= (f32)vid.width/(f32)vid.height*0.75;
	else if(r_fovmode.value == 2) aspect *= yaspectscale.value;
	r_viewchanged = 1;
	R_SetVrect(pvrect, &r_refdef.vrect, lineadj);
	r_refdef.horizontalFieldOfView = 2.0 * tan(r_refdef.fov_x / 360 * M_PI);
	r_refdef.fvrectx = (f32)r_refdef.vrect.x;
	r_refdef.fvrectx_adj = (f32)r_refdef.vrect.x - 0.5;
	r_refdef.vrect_x_adj_shift20 =
		((s64)r_refdef.vrect.x << 20) + (1 << 19) - 1;
	r_refdef.fvrecty = (f32)r_refdef.vrect.y;
	r_refdef.fvrecty_adj = (f32)r_refdef.vrect.y - 0.5;
	r_refdef.vrectright = r_refdef.vrect.x + r_refdef.vrect.width;
	r_refdef.vrectright_adj_shift20 =
		(r_refdef.vrectright << 20) + (1 << 19) - 1;
	r_refdef.fvrectright = (f32)r_refdef.vrectright;
	r_refdef.fvrectright_adj = (f32)r_refdef.vrectright - 0.5;
	r_refdef.vrectbottom = r_refdef.vrect.y + r_refdef.vrect.height;
	r_refdef.fvrectbottom = (f32)r_refdef.vrectbottom;
	r_refdef.fvrectbottom_adj = (f32)r_refdef.vrectbottom - 0.5;
	r_refdef.aliasvrect.x = r_refdef.vrect.x;
	r_refdef.aliasvrect.y = r_refdef.vrect.y;
	r_refdef.aliasvrect.width = r_refdef.vrect.width;
	r_refdef.aliasvrect.height = r_refdef.vrect.height;
	r_refdef.aliasvrectright =
		r_refdef.aliasvrect.x + r_refdef.aliasvrect.width;
	r_refdef.aliasvrectbottom =
		r_refdef.aliasvrect.y + r_refdef.aliasvrect.height;
	pixelAspect = aspect;
	static f32 xOrigin, yOrigin;
	xOrigin = r_refdef.xOrigin;
	yOrigin = r_refdef.yOrigin;
	f32 screenAspect = r_refdef.vrect.width * pixelAspect /
		r_refdef.vrect.height;
	static f32 verticalFieldOfView;
	// 320*200 1.0 pixelAspect = 1.6 screenAspect
	// 320*240 1.0 pixelAspect = 1.3333 screenAspect
	// proper 320*200 pixelAspect = 0.8333333
	verticalFieldOfView = r_refdef.horizontalFieldOfView / screenAspect;
	// values for perspective projection
	// if math were exact, the values would range from 0.5 to to range+0.5
	// hopefully they wll be in the 0.000001 to range+.999999 and truncate
	// the polygon rasterization will never render in the first row or
	// column but will definately render in the [range] row and column, so
	// adjust the buffer origin to get an exact edge to edge fill
	xcenter = ((f32)r_refdef.vrect.width*0.5)+r_refdef.vrect.x-0.5;
	aliasxcenter = xcenter;
	ycenter = ((f32)r_refdef.vrect.height*0.5)+r_refdef.vrect.y-0.5;
	aliasycenter = ycenter;
	xscale = r_refdef.vrect.width / r_refdef.horizontalFieldOfView;
	aliasxscale = xscale;
	xscaleinv = 1.0 / xscale;
	yscale = xscale * pixelAspect;
	aliasyscale = yscale;
	yscaleinv = 1.0 / yscale;
	xscaleshrink = (r_refdef.vrect.width-6)/r_refdef.horizontalFieldOfView;
	yscaleshrink = xscaleshrink * pixelAspect;
	screenedge[0].normal[0] = // left side clip
		-1.0 / (xOrigin * r_refdef.horizontalFieldOfView);
	screenedge[0].normal[1] = 0;
	screenedge[0].normal[2] = 1;
	screenedge[0].type = PLANE_ANYZ;
	screenedge[1].normal[0] = // right side clip
		1.0 / ((1.0 - xOrigin) * r_refdef.horizontalFieldOfView);
	screenedge[1].normal[1] = 0;
	screenedge[1].normal[2] = 1;
	screenedge[1].type = PLANE_ANYZ;
	screenedge[2].normal[0] = 0; // top side clip
	screenedge[2].normal[1] = -1.0 / (yOrigin * verticalFieldOfView);
	screenedge[2].normal[2] = 1;
	screenedge[2].type = PLANE_ANYZ;
	screenedge[3].normal[0] = 0; // bottom side clip
	screenedge[3].normal[1] = 1.0 / ((1.0 - yOrigin) * verticalFieldOfView);
	screenedge[3].normal[2] = 1;
	screenedge[3].type = PLANE_ANYZ;
	for(s32 i = 0; i < 4; i++)
		VectorNormalize(screenedge[i].normal);
	f32 res_scale = sqrt((f64)(r_refdef.vrect.width * r_refdef.vrect.height)
			/ (320.0*152.0)) * (2.0/r_refdef.horizontalFieldOfView);
	r_aliastransition = r_aliastransbase.value * res_scale;
	r_resfudge = r_aliastransadj.value * res_scale;
	D_ViewChanged();
}

static void R_MarkLeaves()
{
	if(r_oldviewleaf == r_viewleaf) return;
	r_visframecount++;
	r_oldviewleaf = r_viewleaf;
	u8 *vis = r_novis.value ? Mod_NoVisPVS(cl.worldmodel)
		: Mod_LeafPVS(r_viewleaf, cl.worldmodel);
	for(s32 i = 0; i < cl.worldmodel->numleafs; i++){
		if(vis[i >> 3] & (1 << (i & 7))){
			mnode_t *nd = (mnode_t *) & cl.worldmodel->leafs[i + 1];
			do{
				if((u32)nd->visframe == r_visframecount)
					break;
				nd->visframe = r_visframecount;
				nd = nd->parent;
			}while(nd);
		}
	}
}

static void R_SetupEntityTransform(entity_t *e)
{ // johnfitz -- set up transform part of lerpdata
	if(e->lerpflags & LERP_RESETMOVE){ // kill any lerps in progress
		e->movelerpstart = 0;
		VectorCopy(e->origin, e->previousorigin);
		VectorCopy(e->origin, e->currentorigin);
		VectorCopy(e->angles, e->previousangles);
		VectorCopy(e->angles, e->currentangles);
		e->lerpflags -= LERP_RESETMOVE;
	}else if(!VectorCompare(e->origin, e->currentorigin) ||
		!VectorCompare(e->angles, e->currentangles)){ // origin/angles
		e->movelerpstart = cl.time; // changed, start new lerp
		VectorCopy(e->currentorigin, e->previousorigin);
		VectorCopy(e->origin, e->currentorigin);
		VectorCopy(e->currentangles, e->previousangles);
		VectorCopy(e->angles, e->currentangles);
	}
	if(!r_lerpmove.value || e == &cl.viewent ||
			!(e->lerpflags & LERP_MOVESTEP))
		return;
	// Baker: Software uses these fields in too many places, so the lerp
	// goes straight into them
	f32 blend = CLAMP(0.f, (cl.time - e->movelerpstart) /
			(e->lerpflags & LERP_FINISH ?
			 e->lerpfinish - e->movelerpstart : 0.1f), 1.f);
	vec3_t d;
	VectorSubtract(e->currentorigin, e->previousorigin, d);
	for(s32 i = 0; i < 3; i++)
		e->origin[i] = e->previousorigin[i] + d[i] * blend;
	VectorSubtract(e->currentangles, e->previousangles, d);
	for(s32 i = 0; i < 3; i++){
		if(d[i] > 180) d[i] -= 360;
		if(d[i] < -180) d[i] += 360;
		e->angles[i] = e->previousangles[i] + d[i] * blend;
	}
}

static void R_SetEntityModelOrg(entity_t *ent)
{
	VectorCopy(ent->origin, r_entorigin);
	VectorSubtract(r_origin, r_entorigin, modelorg);
}
 
static f32 R_EntityAlpha(entity_t *ent)
{ return ent->alpha ? (f32)ent->alpha/255 : 1; }
 
static void R_ClampAliasLighting(alight_t *light)
{
	if(light->ambientlight > 128)
		light->ambientlight = 128;
	if(light->ambientlight + light->shadelight > 192)
		light->shadelight = 192 - light->ambientlight;
}
 
static void R_AddDlightsToAliasLight(alight_t *light, vec3_t origin)
{
	for(s32 lnum = 0; lnum < MAX_DLIGHTS; lnum++){
		dlight_t *dl = &cl_dlights[lnum];
		if(!dl->radius || dl->die < cl.time)
			continue;
		vec3_t dist;
		VectorSubtract(origin, dl->origin, dist);
		f32 len = VectorLength(dist);
		f32 add = dl->radius - len;
		if(add > 0)
			light->ambientlight += add;
	}
}

// Hue of the coloured light at unit brightness. The brightness itself comes
// from the colormap like with uncoloured lighting, so ambient/shade/clamps/
// dlights apply. lit is the R_LightPoint() level the light started from,
// anything added on top of it is white.
static void R_SetupAliasTint(alight_t *light, s32 lit)
{
	f32 white = light->ambientlight - lit;
	if(white < 0)
		white = 0;
	f32 sum = lightcolor[0] + lightcolor[1] + lightcolor[2] + 3 * white;
	if(!r_rgblighting.value || sum <= 0){
		colored_aliaslight = 0;
		return;
	}
	for(s32 i = 0; i < 3; i++)
		alias_tint[i] = (lightcolor[i] + white) * (3 * 256) / sum;
}
 
// r_ambient, which R_LightPoint() applied in WinQuake. Here it leaves the zero
// test for a second sample working, and R_SetupAliasTint() takes it as white
static inline s32 R_WithAmbient(s32 light)
{ return light < r_refdef.ambientlight ? r_refdef.ambientlight : light; }

static void R_SetupAliasEntityLighting(alight_t *lighting)
{
	s32 lit = R_LightPoint(currententity->origin, 0.f);
	if(!lit)
		lit = R_LightPoint(currententity->origin,
				currententity->model->maxs[2]*0.5f);
	s32 j = R_WithAmbient(lit);
	lighting->ambientlight = j;
	lighting->shadelight = j;
	lighting->plightvec = aliaslightvec;
	R_AddDlightsToAliasLight(lighting, currententity->origin);
	R_SetupAliasTint(lighting, lit);
	R_ClampAliasLighting(lighting);
}
 
static void R_DrawAliasEntity()
{
	// always track the move state, so turning r_lerpmove on doesn't start
	// a lerp from wherever the entity was when it was turned off
	if(!(cls.demoplayback && (cl.paused & 2)))
		R_SetupEntityTransform(currententity); // move lerp
	R_SetEntityModelOrg(currententity);
	R_AliasValidateFrame();
	// see if the bounding box lets us trivially reject
	// also sets trivial accept status
	if(!R_AliasCheckBBox())
		return;
	// hidden behind what's already in the z-buffer: skip the lighting,
	// transform and triangle setup, but keep its animation state current.
	// The workers can only test what's in their bands when they get to it.
	bool cull = r_ent_occ_active && r_threaded;
	if(r_ent_occ_active && !cull && R_AliasOccluded()){
		R_AliasSyncLerpState();
		r_amodels_occluded++;
		return;
	}
	alight_t lighting;
	R_SetupAliasEntityLighting(&lighting);
	cur_ent_alpha = r_entalpha.value == 1 ? R_EntityAlpha(currententity):1;
	if(colored_aliaslight &&
		nameInList(r_fullbright_list.string,currententity->model->name))
		colored_aliaslight = 0;
	R_AliasDrawModel(&lighting, cull);
}
 
static void R_DrawEntitiesOnList()
{
	if(!r_drawentities.value)
		return;
	for(s32 i = 0; i < cl_numvisedicts; i++){
		currententity = cl_visedicts[i];
		if(currententity == &cl_entities[cl.viewentity]
			&& !chase_active.value) continue;//don't draw the player
		if(r_showtris.value && currententity->model)
			R_DebugDrawBBox(currententity->origin,currententity->
					model->mins,currententity->model->maxs);
		switch(currententity->model->type){
		case mod_sprite:
			R_SetEntityModelOrg(currententity);
			R_DrawSprite();
			break;
		case mod_alias:
			R_DrawAliasEntity();
			break;
		default:
			break;
		}
	}
	cur_ent_alpha = 1;
}
 
static void R_SetupViewModelLighting()
{
	s32 lit = R_LightPoint(currententity->origin, 0.f);
	if(lit < 24)
		lit = 24; // always give some light on the gun
	s32 j = R_WithAmbient(lit);
	r_viewlighting.ambientlight = j;
	r_viewlighting.shadelight = j;
	r_viewlighting.plightvec = aliaslightvec;
	R_AddDlightsToAliasLight(&r_viewlighting, currententity->origin);
	R_SetupAliasTint(&r_viewlighting, lit);
	R_ClampAliasLighting(&r_viewlighting);
}
 
static void R_DrawViewModel()
{
	if(!r_drawviewmodel.value || cl.items & IT_INVISIBILITY
		|| cl.stats[STAT_HEALTH] <= 0 || !cl.viewent.model
		|| cl.viewent.model->type != mod_alias || chase_active.value)
                return;
	currententity = &cl.viewent;
	R_SetEntityModelOrg(currententity);
	R_AliasValidateFrame();
	VectorCopy(vup, viewlightvec);
	VectorInverse(viewlightvec);
	R_SetupViewModelLighting();
	R_AliasDrawModel(&r_viewlighting, 0);
}
 
static f64 R_ClipPlaneDist(vec3_t point, s32 i)
{
	f64 d = DotProduct(point, view_clipplanes[i].normal);
	d -= view_clipplanes[i].dist;
	return d;
}
 
static s32 R_BmodelCheckRotated(model_t *clmodel)
{
	s32 clipflags = 0;
	// the model turns around its origin, so measure from there, not from
	// the middle of its box, which can be far away from it
	f32 radius = clmodel->rmaxs[0];
	for(s32 i = 0; i < 4; i++){
		f64 d = R_ClipPlaneDist(currententity->origin, i);
		if(d <= -radius)
			return BMODEL_FULLY_CLIPPED;
		if(d <= radius)
			clipflags |= (1 << i);
	}
	return clipflags;
}
 
static s32 R_BmodelCheckAxisAligned(f32 *minmaxs)
{
	s32 clipflags = 0;
	for(s32 i = 0; i < 4; i++){
		s32 *pindex = frustum_indexes[i];
		vec3_t acceptpt, rejectpt;
		for(s32 k = 0; k < 3; k++){
			rejectpt[k] = minmaxs[pindex[k]];
			acceptpt[k] = minmaxs[pindex[3 + k]];
		}
		if(R_ClipPlaneDist(rejectpt, i) <= 0)
			return BMODEL_FULLY_CLIPPED;
		if(R_ClipPlaneDist(acceptpt, i) <= 0)
			clipflags |= (1 << i);
	}
	return clipflags;
}
 
static s32 R_BmodelCheckBBox(model_t *clmodel, f32 *minmaxs)
{
	if(currententity->angles[0] || currententity->angles[1]
			|| currententity->angles[2])
		return R_BmodelCheckRotated(clmodel);
	return R_BmodelCheckAxisAligned(minmaxs);
}
 
static s32 R_SetupBmodelAlpha(entity_t *ent)
{
	if(r_entalpha.value != 1){
		cur_ent_alpha = 1;
		return 1;
	}
	if(!r_alphapass){
		if(ent->alpha && ent->alpha != 255){
			r_foundtranswater = 1;
			return 0;
		}
		return 1;
	}
	cur_ent_alpha = R_EntityAlpha(ent);
	return 1;
}
 
static void R_BmodelWorldBounds(entity_t *ent, f32 *minmaxs)
{ // a turned model can be anywhere within its distance from its origin
	model_t *clmodel = ent->model;
	bool rotated = ent->angles[0] || ent->angles[1] || ent->angles[2];
	const f32 *mn = rotated ? clmodel->rmins : clmodel->mins;
	const f32 *mx = rotated ? clmodel->rmaxs : clmodel->maxs;
	for(s32 j = 0; j < 3; j++){
		minmaxs[j] = ent->origin[j] + mn[j];
		minmaxs[3 + j] = ent->origin[j] + mx[j];
	}
}
 
static void R_MarkBmodelDlights(model_t *clmodel)
{
	if(clmodel->firstmodelsurface == 0)
		return;
	for(s32 k = 0; k < MAX_DLIGHTS; k++){
		if(cl_dlights[k].die < cl.time || !cl_dlights[k].radius)
			continue;
		R_MarkLights(&cl_dlights[k], 1ULL << k,
			clmodel->nodes + clmodel->hulls[0].firstclipnode);
	}
}
 

static void R_DrawBmodelPolygons(model_t *clmodel, f32 *minmaxs, s32 clipflags)
{
	r_pefragtopnode = NULL;
	for(s32 j = 0; j < 3; j++){
		r_emins[j] = minmaxs[j];
		r_emaxs[j] = minmaxs[3 + j];
	}
	R_SplitEntityOnNode2(cl.worldmodel->nodes);
	if(!r_pefragtopnode)
		return;
	// the world traversal didn't get to it: it's off screen or hidden
	if(r_pefragtopnode->keypass != r_worldpass)
		return;
	currententity->topnode = r_pefragtopnode;
	if(r_pefragtopnode->contents >= 0){
		// not a leaf; has to be clipped to the world BSP
		r_clipflags = clipflags;
		R_DrawSolidClippedSubmodelPolygons(clmodel);
	}else{
// falls entirely in one leaf, so we just put all the edges in the edge list and
// let 1/z sorting handle drawing order
		R_DrawSubmodelPolygons(clmodel, clipflags);
	}
	currententity->topnode = NULL;
}

static void R_RestoreWorldView(vec3_t oldorigin)
{
	VectorCopy(base_vpn, vpn);
	VectorCopy(base_vup, vup);
	VectorCopy(base_vright, vright);
	VectorCopy(oldorigin, modelorg);
	R_TransformFrustum();
}
 

static void R_DrawBrushEntity(vec3_t oldorigin)
{
	model_t *clmodel = currententity->model;
	f32 minmaxs[6];
	if(!R_SetupBmodelAlpha(currententity))
		return;
	// none of its faces get drawn in the alpha pass, see R_SurfaceInPass
	if(r_alphapass && !(cur_ent_alpha < 1 && r_entalpha.value == 1) &&
			(currententity->alpha || !clmodel->translucentfaces))
		return;
	// see if the bounding box lets us trivially reject
	// also sets trivial accept status
	R_BmodelWorldBounds(currententity, minmaxs);
	s32 clipflags = R_BmodelCheckBBox(clmodel, minmaxs);
	if(clipflags == BMODEL_FULLY_CLIPPED)
		return;
	r_bmodel_occ_faces = r_bmodel_occ_active &&
		(r_alphapass || cur_ent_alpha >= 1);
	// can nothing of it show up, because opaque world faces nearer than
	// all of it cover everything it could be drawn on?
	if(r_bmodel_occ_faces && R_CoarseOcclusionTestBox(minmaxs, minmaxs+3)){
		r_bmodels_occluded++; // skips everything below
		r_bmodel_occ_faces = 0;
		return;
	}
	R_SetEntityModelOrg(currententity);
	r_pcurrentvertbase = clmodel->vertexes;
	R_RotateBmodel();
	R_MarkBmodelDlights(clmodel);
	R_DrawBmodelPolygons(clmodel, minmaxs, clipflags);
	r_bmodel_occ_faces = 0;
	R_RestoreWorldView(oldorigin);
}
 
static void R_DrawBEntitiesOnList()
{
	if(!r_drawentities.value)
		return;
	vec3_t oldorigin;
	VectorCopy(modelorg, oldorigin);
	insubmodel = 1;
	r_dlightframecount = r_framecount;
	for(s32 i = 0; i < cl_numvisedicts; i++){
		currententity = cl_visedicts[i];
		if(currententity->model->type == mod_brush)
			R_DrawBrushEntity(oldorigin);
	}
	insubmodel = 0;
	cur_ent_alpha = 1;
}

static void R_TimeStamp(s32 i)
{ if(r_speeds.value) d_times[i] = Sys_DoubleTime(); }

// The world and the brush entities into the edge list, then the spans to the
// screen. The main pass draws everything opaque, the alpha pass what's left.
// Timed from d_times[t] to t + 3.
static void R_EdgeDrawing(bool alpha, s32 t)
{
	r_alphapass = alpha;
	R_BeginEdgeFrame();
	R_TimeStamp(t);
	R_RenderWorld();
	R_TimeStamp(t + 1);
	R_DrawBEntitiesOnList();
	R_TimeStamp(t + 2);
	R_ScanEdges();
	R_TimeStamp(t + 3);
}

void R_RenderView()
{
	R_ThreadsBeginFrame();
	R_AllocShowTris();
	R_TimeStamp(0);
	R_SetupFrame();
	R_MarkLeaves(); // done here so we know if we're in water
	fesetround(FE_TOWARDZERO); // Sys_LowFPPrecision set this in WinQuake
	R_ThreadsFPEnv(); // and the workers too
	R_CoarseOcclusionBeginFrame();
	r_foundtranswater = 0;
	R_EdgeDrawing(0, 1);
	R_TimeStamp(5);
	R_DrawEntitiesOnList();
	R_TimeStamp(6);
	if(r_foundtranswater)
		R_EdgeDrawing(1, 7);
	else for(s32 i = 7; i <= 10; i++)
		R_TimeStamp(i);
	R_DrawViewModel();
	R_TimeStamp(11);
	R_DrawParticles();
	R_TimeStamp(12);
	if(r_dowarp) D_WarpScreen();
	R_TimeStamp(13);
	R_DrawFog();
	R_ThreadsEndFrame();
	R_TimeStamp(14);
	V_SetContentsColor(r_viewleaf->contents);
	if(r_showtris.value){
		for(s32 i = 0; i < r_numdebugpoints; i++)
			R_DebugDrawPoint(r_debugpoints[i].origin);
		for(s32 i = 0; i < r_numdebuglines; i++){
			R_DrawDebugLine(r_debuglines[i].x0, r_debuglines[i].y0,
					r_debuglines[i].x1, r_debuglines[i].y1,
					r_showtris_color.value);
		}
		r_numdebuglines = 0;
	}
	if(coccl_enable && r_speeds.value)
		Con_Printf("Coarse rejected: %d (nodes %d)\n",
				r_coarse_rejected, r_nodes_occluded);
	if(r_bmodel_occ_active && r_speeds.value)
		Con_Printf("Brush entities occluded: %d (faces %d)\n",
				r_bmodels_occluded, r_bmodel_faces_occluded);
	if(r_ent_occ_active && r_speeds.value)
		Con_Printf("Entities occluded: %d (drawn %d)\n",
				r_amodels_occluded, r_amodels_drawn);
	fesetround(FE_TONEAREST);
}
