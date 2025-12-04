// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

static vec3_t viewlightvec;
static alight_t r_viewlighting = { 128, 192, viewlightvec };
static f32 verticalFieldOfView;
static f32 xOrigin, yOrigin;
static u8 warpbuffer[WARP_WIDTH * WARP_HEIGHT];

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

void R_InitTurb()
{
	for(s32 i = 0; i < (SIN_BUFFER_SIZE); i++){
		sintable[i] = AMP + sin(i * 3.14159 * 2 / CYCLE) * AMP;
		intsintable[i] = AMP2 + sin(i * 3.14159 * 2 / CYCLE) * AMP2;
	}
	r_warpbuffer = warpbuffer;
}

void R_ViewChangedCallback(SDL_UNUSED cvar_t*cvar)
{ R_ViewChanged(&r_refdef.vrect, sb_lines, pixelAspect); }

void R_Init()
{
	Cmd_AddCommand("timerefresh", R_TimeRefresh_f);
	Cmd_AddCommand("pointfile", R_ReadPointFile_f);
	Cmd_AddCommand("fog", Fog_FogCommand_f);
	Cmd_AddCommand("vid_setmode", VID_VidSetModeCommand_f);
	Cvar_RegisterVariable(&r_draworder);
	Cvar_RegisterVariable(&r_speeds);
	Cvar_RegisterVariable(&r_timegraph);
	Cvar_RegisterVariable(&r_graphheight);
	Cvar_RegisterVariable(&r_drawflat);
	Cvar_RegisterVariable(&r_ambient);
	Cvar_RegisterVariable(&r_clearcolor);
	Cvar_RegisterVariable(&r_waterwarp);
	Cvar_RegisterVariable(&r_fullbright);
	Cvar_RegisterVariable(&r_drawentities);
	Cvar_RegisterVariable(&r_drawviewmodel);
	Cvar_RegisterVariable(&r_aliasstats);
	Cvar_RegisterVariable(&r_dspeeds);
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
	Cvar_RegisterVariable(&r_litwater);
	Cvar_RegisterVariable(&r_novis);
	Cvar_RegisterVariable(&r_particlescale);
	Cvar_RegisterVariable(&r_fovmode);
	Cvar_RegisterVariable(&r_hlwater);
	Cvar_RegisterVariable(&vid_cwidth);
	Cvar_RegisterVariable(&vid_cheight);
	Cvar_RegisterVariable(&vid_cwmode);
	Cvar_RegisterVariable(&scr_uixscale);
	Cvar_RegisterVariable(&scr_uiyscale);
	Cvar_RegisterVariable(&yaspectscale);
	Cvar_RegisterVariable(&scr_lockuiscale);
	Cvar_RegisterVariable(&r_mipscale);
	Cvar_RegisterVariable(&scr_menubgstyle);
	Cvar_RegisterVariable(&lyr_main);
	Cvar_RegisterVariable(&lyr_sbar);
	Cvar_RegisterVariable(&lyr_menu);
	Cvar_RegisterVariable(&lyr_centerprint);
	Cvar_RegisterVariable(&lyr_console);
	Cvar_RegisterVariable(&lyr_notify);
	Cvar_RegisterVariable(&lyr_crosshair);
	Cvar_SetCallback(&r_labmixpal, build_color_mix_lut);
	Cvar_SetCallback(&r_rgblighting, D_FlushCaches);
	Cvar_SetCallback(&r_fogbrightness, Fog_SetPalIndex);
	Cvar_SetCallback(&r_wateralpha, R_SetWateralpha_f);
	Cvar_SetCallback(&r_lavaalpha, R_SetLavaalpha_f);
	Cvar_SetCallback(&r_telealpha, R_SetTelealpha_f);
	Cvar_SetCallback(&r_slimealpha, R_SetSlimealpha_f);
	Cvar_SetCallback(&r_fovmode, R_ViewChangedCallback);
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

void Pal_ParseWorldspawn ()
{ // called at map load
	s8 key[128], value[4096];
	const s8 *data = COM_Parse(cl.worldmodel->entities);
	if(!data || com_token[0] != '{') return; // error
	while(1){
		if(!data) return; // error
		if(com_token[0] == '}') break; // end of worldspawn
		if(com_token[0] == '_')q_strlcpy(key, com_token+1, sizeof(key));
		else q_strlcpy(key, com_token, sizeof(key));
		while(key[0] && key[strlen(key)-1] == ' ') // no trailing spaces
			key[strlen(key)-1] = 0;
		data = COM_ParseEx(data, CPE_ALLOWTRUNC);
		if(!data) return; // error
		q_strlcpy(value, com_token, sizeof(value));
		if(!strcmp("palette", key)){
			s8 pal[MAX_OSPATH];
			s8 cmap[MAX_OSPATH];
			sscanf(value, "%s", pal);
			q_strlcpy(cmap, pal, MAX_OSPATH);
			q_strlcat(cmap, "_colormap", MAX_OSPATH);
			Con_DPrintf("Requested palettes %s %s\n", pal, cmap);
			SetWorldPal(pal, cmap);
		}
	}
	Fog_SetPalIndex(0);
}

void R_NewMap()
{
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
	skybox_name[0] = 0;
	r_skymade = 0;
	Sky_NewMap();
	Fog_ParseWorldspawn();
	R_ParseWorldspawn();
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
	pvrect->width &= ~7;
	pvrect->height = pvrectin->height * size;
	if(pvrect->height > pvrectin->height - lineadj)
		pvrect->height = pvrectin->height - lineadj;
	pvrect->height &= ~1;
	pvrect->x = (pvrectin->width - pvrect->width) / 2;
	pvrect->y = (h - pvrect->height) / 2;
}

void R_ViewChanged(vrect_t *pvrect, s32 lineadj, f32 aspect)
{ // Called every time the vid structure or r_refdef changes.
 // Guaranteed to be called before the first refresh
	if (r_fovmode.value == 1) aspect *= (f32)vid.width/(f32)vid.height*0.75;
	else if (r_fovmode.value == 2) aspect *= yaspectscale.value;
	r_viewchanged = 1;
	R_SetVrect(pvrect, &r_refdef.vrect, lineadj);
	r_refdef.horizontalFieldOfView = 2.0 * tan(r_refdef.fov_x / 360 * M_PI);
	r_refdef.fvrectx = (f32)r_refdef.vrect.x;
	r_refdef.fvrectx_adj = (f32)r_refdef.vrect.x - 0.5;
	r_refdef.vrect_x_adj_shift20 = (r_refdef.vrect.x << 20) + (1 << 19) - 1;
	r_refdef.fvrecty = (f32)r_refdef.vrect.y;
	r_refdef.fvrecty_adj = (f32)r_refdef.vrect.y - 0.5;
	r_refdef.vrectright = r_refdef.vrect.x + r_refdef.vrect.width;
	r_refdef.vrectright_adj_shift20 =
		(r_refdef.vrectright << 20) + (1 << 19) - 1;
	r_refdef.fvrectright = (f32)r_refdef.vrectright;
	r_refdef.fvrectright_adj = (f32)r_refdef.vrectright - 0.5;
	r_refdef.vrectrightedge = (f32)r_refdef.vrectright - 0.99;
	r_refdef.vrectbottom = r_refdef.vrect.y + r_refdef.vrect.height;
	r_refdef.fvrectbottom = (f32)r_refdef.vrectbottom;
	r_refdef.fvrectbottom_adj = (f32)r_refdef.vrectbottom - 0.5;
	r_refdef.aliasvrect.x = (s32)(r_refdef.vrect.x);
	r_refdef.aliasvrect.y = (s32)(r_refdef.vrect.y);
	r_refdef.aliasvrect.width = (s32)(r_refdef.vrect.width);
	r_refdef.aliasvrect.height = (s32)(r_refdef.vrect.height);
	r_refdef.aliasvrectright =
		r_refdef.aliasvrect.x + r_refdef.aliasvrect.width;
	r_refdef.aliasvrectbottom =
		r_refdef.aliasvrect.y + r_refdef.aliasvrect.height;
	pixelAspect = aspect;
	xOrigin = r_refdef.xOrigin;
	yOrigin = r_refdef.yOrigin;
	f32 screenAspect = r_refdef.vrect.width * pixelAspect /
		r_refdef.vrect.height;
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
		/ (320.0 * 152.0)) * (2.0 / r_refdef.horizontalFieldOfView);
	r_aliastransition = r_aliastransbase.value * res_scale;
	r_resfudge = r_aliastransadj.value * res_scale;
	D_ViewChanged();
}

void R_MarkLeaves()
{
	if(r_oldviewleaf == r_viewleaf) return;
	r_visframecount++;
	r_oldviewleaf = r_viewleaf;
	u8 *vis = Mod_LeafPVS(r_viewleaf, cl.worldmodel);
	for(s32 i = 0; i < cl.worldmodel->numleafs; i++){
		if(vis[i >> 3] & (1 << (i & 7))){
			mnode_t *nd = (mnode_t *) & cl.worldmodel->leafs[i + 1];
			do {
				if((u32)nd->visframe == r_visframecount)
					break;
				nd->visframe = r_visframecount;
				nd = nd->parent;
			} while(nd);
		}
	}
}

void R_DrawEntitiesOnList()
{
	f32 lightvec[3] = { -1, 0, 0 }; // FIXME: remove and do real lighting
	if(!r_drawentities.value)
		return;
	for(s32 i = 0; i < cl_numvisedicts; i++){
		currententity = cl_visedicts[i];
		if(currententity == &cl_entities[cl.viewentity])
			continue; // don't draw the player
		switch(currententity->model->type){
		case mod_sprite:
			VectorCopy(currententity->origin, r_entorigin);
			VectorSubtract(r_origin, r_entorigin, modelorg);
			R_DrawSprite();
			break;
		case mod_alias:
			VectorCopy(currententity->origin, r_entorigin);
			VectorSubtract(r_origin, r_entorigin, modelorg);
			// see if the bounding box lets us trivially reject, also sets
			// trivial accept status
			if(R_AliasCheckBBox()){
				s32 j = R_LightPoint(currententity->origin);
				alight_t lighting;
				lighting.ambientlight = j;
				lighting.shadelight = j;
				lighting.plightvec = lightvec;
				for(s32 lnum = 0; lnum < MAX_DLIGHTS; lnum++)
					if(cl_dlights[lnum].die >= cl.time){
						vec3_t dist;
						VectorSubtract(currententity->
							       origin,
							       cl_dlights[lnum].
							       origin, dist);
						f32 add =
							cl_dlights[lnum].radius
							- VectorLength(dist);
						if (add > 0)
							lighting.ambientlight +=
								add;
					}
				// clamp lighting so it doesn't overbright as much
				if(lighting.ambientlight > 128)
					lighting.ambientlight = 128;
				if(lighting.ambientlight +
				    lighting.shadelight > 192)
					lighting.shadelight =
					    192 - lighting.ambientlight;
				cur_ent_alpha = currententity->alpha && r_entalpha.value == 1 ?
					(f32)currententity->alpha/255 : 1;
				if(colored_aliaslight &&
					nameInList(r_fullbright_list.string, currententity->model->name))
					colored_aliaslight = 0;
				R_AliasDrawModel(&lighting);
			}
			break;
		default:
			break;
		}
	}
	cur_ent_alpha = 1;
}

void R_DrawViewModel()
{
	f32 lightvec[3] = { -1, 0, 0 }; // FIXME: remove and do real lighting
	if(!r_drawviewmodel.value || cl.items & IT_INVISIBILITY
		|| cl.stats[STAT_HEALTH] <= 0 || !cl.viewent.model)
		return;
	currententity = &cl.viewent;
	VectorCopy(currententity->origin, r_entorigin);
	VectorSubtract(r_origin, r_entorigin, modelorg);
	VectorCopy(vup, viewlightvec);
	VectorInverse(viewlightvec);
	s32 j = R_LightPoint(currententity->origin);
	if(j < 24)
		j = 24; // allways give some light on the gun
	r_viewlighting.ambientlight = j;
	r_viewlighting.shadelight = j;
	colored_aliaslight = r_rgblighting.value && lit_loaded ? 1 : 0;
	for(s32 lnum = 0; lnum < MAX_DLIGHTS; lnum++){ // add dynamic lights
		dlight_t *dl = &cl_dlights[lnum];
		if(!dl->radius || dl->die < cl.time)
			continue;
		vec3_t dist;
		VectorSubtract(currententity->origin, dl->origin, dist);
		f32 add = dl->radius - VectorLength(dist);
		if(add > 150 && VectorLength(dist) < 50) // hack in the muzzleflash
			colored_aliaslight = 0; // FIXME and do it properly
		if(add > 0)
			r_viewlighting.ambientlight += add;
	}
	// clamp lighting so it doesn't overbright as much
	if(r_viewlighting.ambientlight > 128)
		r_viewlighting.ambientlight = 128;
	if(r_viewlighting.ambientlight + r_viewlighting.shadelight > 192)
		r_viewlighting.shadelight = 192 - r_viewlighting.ambientlight;
	r_viewlighting.plightvec = lightvec;
	R_AliasDrawModel(&r_viewlighting);
}

s32 R_BmodelCheckBBox(model_t *clmodel, f32 *minmaxs)
{
	s32 clipflags = 0;
	if(currententity->angles[0] || currententity->angles[1]
	    || currententity->angles[2]){
		for(s32 i = 0; i < 4; i++){
			f64 d = DotProduct(currententity->origin,
				       view_clipplanes[i].normal);
			d -= view_clipplanes[i].dist;
			if(d <= -clmodel->radius)
				return BMODEL_FULLY_CLIPPED;
			if(d <= clmodel->radius)
				clipflags |= (1 << i);
		}
	} else {
		for(s32 i = 0; i < 4; i++){
			// generate accept and reject points
			// FIXME: do with fast look-ups or integer tests
			// based on the sign bit of the floating point values
			s32 *pindex = pfrustum_indexes[i];
			vec3_t acceptpt, rejectpt;
			rejectpt[0] = minmaxs[pindex[0]];
			rejectpt[1] = minmaxs[pindex[1]];
			rejectpt[2] = minmaxs[pindex[2]];
			f64 d = DotProduct(rejectpt, view_clipplanes[i].normal);
			d -= view_clipplanes[i].dist;
			if(d <= 0)
				return BMODEL_FULLY_CLIPPED;
			acceptpt[0] = minmaxs[pindex[3 + 0]];
			acceptpt[1] = minmaxs[pindex[3 + 1]];
			acceptpt[2] = minmaxs[pindex[3 + 2]];
			d = DotProduct(acceptpt, view_clipplanes[i].normal);
			d -= view_clipplanes[i].dist;
			if(d <= 0)
				clipflags |= (1 << i);
		}
	}
	return clipflags;
}

void R_DrawBEntitiesOnList()
{
	f32 minmaxs[6];
	if(!r_drawentities.value)
		return;
	vec3_t oldorigin;
	VectorCopy(modelorg, oldorigin);
	insubmodel = 1;
	r_dlightframecount = r_framecount;
	model_t *clmodel; // keep here for the OpenBSD compiler
	for(s32 i = 0; i < cl_numvisedicts; i++){
		currententity = cl_visedicts[i];
		switch(currententity->model->type){
		case mod_brush:
			if(r_entalpha.value == 1){
				if(!r_alphapass && currententity->alpha &&
						currententity->alpha != 255){
					r_foundtranswater = 1;
					continue;
				}
				if(r_alphapass)
					cur_ent_alpha = currententity->alpha ?
						(f32)currententity->alpha/255 : 1;
			}
			else
				cur_ent_alpha = 1;
			clmodel = currententity->model;
			// see if the bounding box lets us trivially reject
			// also sets trivial accept status
			for(s32 j = 0; j < 3; j++){
				minmaxs[j] = currententity->origin[j] +
				    clmodel->mins[j];
				minmaxs[3 + j] = currententity->origin[j] +
				    clmodel->maxs[j];
			}
			s32 clipflags = R_BmodelCheckBBox(clmodel, minmaxs);
			if(clipflags != BMODEL_FULLY_CLIPPED){
				VectorCopy(currententity->origin, r_entorigin);
				VectorSubtract(r_origin, r_entorigin, modelorg);
				r_pcurrentvertbase = clmodel->vertexes;
				// FIXME: stop transforming twice
				R_RotateBmodel();
				// calculate dynamic lighting for bmodel
				// if it's not an instanced model
				if(clmodel->firstmodelsurface != 0){
					for(s32 k = 0; k < MAX_DLIGHTS; k++){
						if((cl_dlights[k].die <cl.time)
						    || (!cl_dlights[k].radius))
							continue;
						R_MarkLights(&cl_dlights[k],
							     1 << k,
							     clmodel->nodes +
							     clmodel->hulls[0].
							     firstclipnode);
					}
				}
				r_pefragtopnode = NULL;
				for(s32 j = 0; j < 3; j++){
					r_emins[j] = minmaxs[j];
					r_emaxs[j] = minmaxs[3 + j];
				}
				R_SplitEntityOnNode2(cl.worldmodel-> nodes);
				if(r_pefragtopnode){
					currententity->topnode =r_pefragtopnode;
					if(r_pefragtopnode->contents >= 0){
						// not a leaf; has to be clipped to the world BSP
						r_clipflags = clipflags;
						R_DrawSolidClippedSubmodelPolygons(clmodel);
					} else {
						// falls entirely in one leaf, so we just put all the
						// edges in the edge list and let 1/z sorting handle
						// drawing order
						R_DrawSubmodelPolygons(clmodel, clipflags);
					}
					currententity->topnode = NULL;
				}
				// put back world rotation and frustum clipping         
				// FIXME: R_RotateBmodel should just work off base_vxx
				VectorCopy(base_vpn, vpn);
				VectorCopy(base_vup, vup);
				VectorCopy(base_vright, vright);
				VectorCopy(base_modelorg, modelorg);
				VectorCopy(oldorigin, modelorg);
				R_TransformFrustum();
			}
			break;

		default:
			break;
		}
	}
	insubmodel = 0;
	cur_ent_alpha = 1;
}

void R_EdgeDrawing()
{
	r_foundtranswater =  r_alphapass = 0;
	R_BeginEdgeFrame();
	if(r_dspeeds.value) d_times[1] = Sys_DoubleTime();
	R_RenderWorld();
	if(r_dspeeds.value) d_times[2] = Sys_DoubleTime();
	R_DrawBEntitiesOnList();
	if(r_dspeeds.value) d_times[3] = Sys_DoubleTime();
	R_ScanEdges();
	if(r_dspeeds.value) d_times[4] = Sys_DoubleTime();
}

void R_EdgeDrawingAlpha()
{
	if(!r_foundtranswater || !r_entalpha.value){
		if(r_dspeeds.value)
			d_times[7]=d_times[8]=d_times[9]=Sys_DoubleTime();
		return;
	}
	r_alphapass = 1;
	R_BeginEdgeFrame();
	if(r_dspeeds.value) d_times[7] = Sys_DoubleTime();
	R_RenderWorld();
	if(r_dspeeds.value) d_times[8] = Sys_DoubleTime();
	R_DrawBEntitiesOnList();
	if(r_dspeeds.value) d_times[9] = Sys_DoubleTime();
	R_ScanEdges();
}

void R_RenderView()
{
	if(r_timegraph.value || r_speeds.value || r_dspeeds.value)
		d_times[0] = Sys_DoubleTime();
	R_SetupFrame();
	R_MarkLeaves(); // done here so we know if we're in water
	R_EdgeDrawing();
	if(r_dspeeds.value) d_times[5] = Sys_DoubleTime();
	R_DrawEntitiesOnList();
	if(r_dspeeds.value) d_times[6] = Sys_DoubleTime();
	R_EdgeDrawingAlpha();
	if(r_dspeeds.value) d_times[10] = Sys_DoubleTime();
	R_DrawViewModel();
	if(r_dspeeds.value) d_times[11] = Sys_DoubleTime();
	R_DrawParticles();
	if(r_dspeeds.value) d_times[12] = Sys_DoubleTime();
	if(r_dowarp) D_WarpScreen();
	if(r_dspeeds.value) d_times[13] = Sys_DoubleTime();
        if(fog_density < 1) R_DrawFog();
	if(r_dspeeds.value) d_times[14] = Sys_DoubleTime();
	V_SetContentsColor(r_viewleaf->contents);
}
