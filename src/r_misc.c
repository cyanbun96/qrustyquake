// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

static f32 map_fallbackalpha;
static f32 map_wateralpha;
static f32 map_lavaalpha;
static f32 map_telealpha;
static f32 map_slimealpha;
static s32 r_frustum_indexes[4 * 6];

void R_ParseWorldspawn()
{
	s8 key[128], value[4096];
	map_fallbackalpha = r_wateralpha.value;
	s32 ct = cl.worldmodel->contentstransparent;
	map_wateralpha = (ct&SURF_DRAWWATER)?r_wateralpha.value:1;
	map_lavaalpha =  (ct&SURF_DRAWLAVA)? r_lavaalpha.value:1;
	map_telealpha =  (ct&SURF_DRAWTELE)? r_telealpha.value:1;
	map_slimealpha = (ct&SURF_DRAWSLIME)?r_slimealpha.value:1;
	const s8 *data = COM_Parse(cl.worldmodel->entities);
	if(!data)return; // error
	if(com_token[0] != '{')return; // error
	while(1){
		data = COM_Parse(data);
		if(!data)return; // error
		if(com_token[0]=='}')break; // end of worldspawn
		if(com_token[0]=='_')q_strlcpy(key, com_token + 1, sizeof(key));
		else q_strlcpy(key, com_token, sizeof(key));
		while(key[0]&&key[strlen(key)-1]==' ') // remove trailing spaces
			key[strlen(key)-1] = 0;
		data = COM_ParseEx(data, CPE_ALLOWTRUNC);
		if(!data)return; // error
		q_strlcpy(value, com_token, sizeof(value));
		if(!strcmp("wateralpha", key))map_wateralpha = atof(value);
		if(!strcmp("lavaalpha", key)) map_lavaalpha  = atof(value);
		if(!strcmp("telealpha", key)) map_telealpha  = atof(value);
		if(!strcmp("slimealpha", key))map_slimealpha = atof(value);
	}
}

void R_SetWateralpha_f(cvar_t *var) 
{ // -- ericw
	if (cls.signon == SIGNONS && cl.worldmodel &&
		!(cl.worldmodel->contentstransparent&SURF_DRAWWATER)
		&& var->value < 1)
		Con_Printf("Map does not appear to be water-vised\n");
	map_wateralpha = var->value;
	map_fallbackalpha = var->value;
}

void R_SetLavaalpha_f(cvar_t *var)
{
	if (cls.signon == SIGNONS && cl.worldmodel &&
		!(cl.worldmodel->contentstransparent&SURF_DRAWLAVA)
		&& var->value && var->value < 1)
		Con_Printf("Map does not appear to be lava-vised\n");
	map_lavaalpha = var->value;
}

void R_SetTelealpha_f(cvar_t *var)
{
	if (cls.signon == SIGNONS && cl.worldmodel &&
		!(cl.worldmodel->contentstransparent&SURF_DRAWTELE)
		&& var->value && var->value < 1)
		Con_Printf("Map does not appear to be tele-vised\n");
	map_telealpha = var->value;
}

void R_SetSlimealpha_f(cvar_t *var)
{
	if (cls.signon == SIGNONS && cl.worldmodel &&
		!(cl.worldmodel->contentstransparent&SURF_DRAWSLIME)
		&& var->value && var->value < 1)
		Con_Printf("Map does not appear to be slime-vised\n");
	map_slimealpha = var->value;
}

f32 R_LiquidAlphaForFlags(s32 flags)
{
	if (flags & SURF_DRAWLAVA)
		return map_lavaalpha > 0 ? map_lavaalpha : map_fallbackalpha;
	else if (flags & SURF_DRAWSLIME)
		return map_slimealpha > 0 ? map_slimealpha : map_fallbackalpha;
	else if (flags & SURF_DRAWTELE)
		return map_telealpha > 0 ? map_telealpha : map_fallbackalpha;
	else
		return map_wateralpha > 0 ? map_wateralpha : map_fallbackalpha;
}

void R_CheckVariables()
{
	static f32 oldbright;
	if (r_fullbright.value != oldbright) {
		oldbright = r_fullbright.value;
		D_FlushCaches(0); // so all lighting changes
	}
}

void R_TimeRefresh_f()
{ // For program optimization
	s32 startangle = r_refdef.viewangles[1];
	f32 start = Sys_DoubleTime();
	for (s32 i = 0; i < 128; i++) {
		r_refdef.viewangles[1] = i / 128.0 * 360.0;
		R_RenderView();
		VID_Update();
	}
	f32 stop = Sys_DoubleTime();
	f32 time = stop - start;
	Con_Printf("%f seconds (%f fps)\n", time, 128 / time);
	r_refdef.viewangles[1] = startangle;
}

void R_LineGraph(s32 x, s32 y, s32 h)
{ // Only called by R_DisplayTime
// FIXME: should be disabled on no-buffer adapters, or should be in the driver
	x += r_refdef.vrect.x;
	y += r_refdef.vrect.y;
	u8 *dest = vid.buffer + vid.width * y + x;
	s32 s = r_graphheight.value;
	if (h > s)
		h = s;
	s32 i = 0;
	for (; i < h; i++, dest -= vid.width * 2) {
		dest[0] = 0xff;
		*(dest - vid.width) = 0x30;
	}
	for (; i < s; i++, dest -= vid.width * 2) {
		dest[0] = 0x30;
		*(dest - vid.width) = 0x30;
	}
}

void R_TimeGraph()
{ // Performance monitoring tool
	static s32 timex;
	static u8 r_timings[MAX_TIMINGS];
	f32 r_time2 = Sys_DoubleTime();
	s32 a = (r_time2 - d_times[0]) / 0.01;
	r_timings[timex] = a;
	a = timex;
	s32 x;
	if (r_refdef.vrect.width <= MAX_TIMINGS)
		x = r_refdef.vrect.width - 1;
	else
		x = r_refdef.vrect.width -
		    (r_refdef.vrect.width - MAX_TIMINGS) / 2;
	do {
		R_LineGraph(x, r_refdef.vrect.height - 2, r_timings[a]);
		if (x == 0)
			break; // screen too small to hold entire thing
		x--;
		a--;
		if (a == -1)
			a = MAX_TIMINGS - 1;
	} while (a != timex);
	timex = (timex + 1) % MAX_TIMINGS;
}

void R_PrintTimes()
{
	f32 r_time2 = Sys_DoubleTime();
	f32 ms = 1000 * (r_time2 - d_times[0]);
	Con_Printf("%5.1f ms %3i/%3i/%3i poly %3i surf\n",
		   ms, c_faceclip, r_polycount, r_drawnpolycount, c_surf);
	c_surf = 0;
}

void R_PrintDSpeeds()
{
	f64 t[17] = {0};
	f64 r_time2 = Sys_DoubleTime();
	for(s32 i = 0; i < 14; i++)
		t[i] = (d_times[i+1] - d_times[i]) * 1000000;
	f64 tot = (r_time2 - d_times[0]) * 1000000;
	Con_Printf("%.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\n"
		,tot,t[0],t[1],t[2],t[3], t[4], t[5], t[6], t[7]);
	Con_Printf("%.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\n",
		t[8], t[9], t[10], t[11], t[12], t[13], t[14]);
}

void R_PrintAliasStats()
{
	Con_Printf("%3i polygon model drawn\n", r_amodels_drawn);
}

void R_TransformFrustum()
{
	for (s32 i = 0; i < 4; i++) {
		vec3_t v, v2;
		v[0] = screenedge[i].normal[2];
		v[1] = -screenedge[i].normal[0];
		v[2] = screenedge[i].normal[1];
		v2[0] = v[1] * vright[0] + v[2] * vup[0] + v[0] * vpn[0];
		v2[1] = v[1] * vright[1] + v[2] * vup[1] + v[0] * vpn[1];
		v2[2] = v[1] * vright[2] + v[2] * vup[2] + v[0] * vpn[2];
		VectorCopy(v2, view_clipplanes[i].normal);
		view_clipplanes[i].dist = DotProduct(modelorg, v2);
	}
}

void R_SetUpFrustumIndexes()
{
	s32 *pindex = r_frustum_indexes;
	for (s32 i = 0; i < 4; i++) {
		for (s32 j = 0; j < 3; j++) {
			if (view_clipplanes[i].normal[j] < 0) {
				pindex[j] = j;
				pindex[j + 3] = j + 3;
			} else {
				pindex[j] = j + 3;
				pindex[j + 3] = j;
			}
		}
		pfrustum_indexes[i] = pindex; // FIXME: do just once at start
		pindex += 6;
	}
}

void R_SetupFrame()
{
	if (cl.maxclients > 1) { // don't allow cheats in multiplayer
		Cvar_Set("r_draworder", "0");
		Cvar_Set("r_fullbright", "0");
		Cvar_Set("r_ambient", "0");
		Cvar_Set("r_drawflat", "0");
	}
	if (r_numsurfs.value) {
		if (surface_p - surfaces > r_maxsurfsseen)
			r_maxsurfsseen = surface_p - surfaces;
		Con_Printf("Used %d of %d surfs; %d max\n",
				surface_p - surfaces, surf_max - surfaces,
				r_maxsurfsseen);
	}
	if (r_numedges.value) {
		s32 edgecount = edge_p - r_edges;

		if (edgecount > r_maxedgesseen)
			r_maxedgesseen = edgecount;

		Con_Printf("Used %d of %d edges; %d max\n", edgecount,
				NUMSTACKEDGES, r_maxedgesseen);
	}
	r_refdef.ambientlight = r_ambient.value;
	if (r_refdef.ambientlight < 0)
		r_refdef.ambientlight = 0;
	if (!sv.active)
		r_draworder.value = 0;	// don't let cheaters look behind walls
	R_CheckVariables();
	R_AnimateLight();
	r_framecount++;
	// build the transformation matrix for the given view angles
	VectorCopy(r_refdef.vieworg, modelorg);
	VectorCopy(r_refdef.vieworg, r_origin);
	AngleVectors(r_refdef.viewangles, vpn, vright, vup);
	r_oldviewleaf = r_viewleaf; // current viewleaf
	r_viewleaf = Mod_PointInLeaf(r_origin, cl.worldmodel);
	r_dowarpold = r_dowarp;
	r_dowarp = r_waterwarp.value && (r_viewleaf->contents<=CONTENTS_WATER);
	vrect_t vrect;
	if ((r_dowarp != r_dowarpold) || r_viewchanged) {
		if (r_dowarp) {
			if ((vid.width <= WARP_WIDTH) &&
					(vid.height <= WARP_HEIGHT)) {
				vrect.x = 0;
				vrect.y = 0;
				vrect.width = vid.width;
				vrect.height = vid.height;
				R_ViewChanged(&vrect, sb_lines, vid.aspect);
			} else {
				f32 w = vid.width;
				f32 h = vid.height;
				if (w > WARP_WIDTH) {
					h *= (f32)WARP_WIDTH / w;
					w = WARP_WIDTH;
				}
				if (h > WARP_HEIGHT) {
					h = WARP_HEIGHT;
					w *= (f32)WARP_HEIGHT / h;
				}
				vrect.x = 0;
				vrect.y = 0;
				vrect.width = (s32)w;
				vrect.height = (s32)h;
				R_ViewChanged(&vrect, (s32)((f32)sb_lines*
					(h/(f32)vid.height)),vid.aspect*(h/w)*
					((f32)vid.width/(f32)vid.height));
			}
		} else {
			vrect.x = 0;
			vrect.y = 0;
			vrect.width = vid.width;
			vrect.height = vid.height;
			R_ViewChanged(&vrect, sb_lines, vid.aspect);
		}
		r_viewchanged = 0;
	}
	// start off with just the four screen edge clip planes
	R_TransformFrustum();
	VectorCopy(vpn, base_vpn); // save base values
	VectorCopy(vright, base_vright);
	VectorCopy(vup, base_vup);
	VectorCopy(modelorg, base_modelorg);
	R_SetSkyFrame();
	R_SetUpFrustumIndexes();
	if(r_cache_thrash) D_AllocCaches(); // allocate a bigger one
	r_cache_thrash = 0;
	c_faceclip = 0; // clear frame counts
	r_polycount = 0;
	r_drawnpolycount = 0;
	r_amodels_drawn = 0;
	r_outofsurfaces = 0;
	r_outofedges = 0;
	D_SetupFrame();
}
