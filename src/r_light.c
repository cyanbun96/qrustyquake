// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

void R_AnimateLight() // light animations
{ // 'm' is normal light, 'a' is no light, 'z' is f64 bright
	f64 t = cl.time * 10;
	s32 i = (s32)t;
	s32 f = r_lerplightstyles.value ? (t - i) * 256 : 0; // 0.8 fraction
	for(s32 j = 0; j < MAX_LIGHTSTYLES; j++){
		lightstyle_t *ls = &cl_lightstyle[j];
		if(!ls->length){
			d_lightstylevalue[j] = 256;
			continue;
		}
		s32 k = ls->map[i % ls->length] - 'a';
		s32 n = ls->map[(i + 1) % ls->length] - 'a';
		// 1 = only smooth gradual changes, keep flickers sharp
		if(r_lerplightstyles.value < 2 && abs(n - k) >= ('m'-'a')/2)
			n = k;
		d_lightstylevalue[j] = k * 22 + (((n - k) * 22 * f) >> 8);
	}
}

static s32 R_LightReachesSurf(dlight_t *light, msurface_t *surf)
{ // could R_AddDynamicLights brighten any of its lightmap samples? The same
  // distances, with room for its truncation to s32 and for rounding. A
  // surface it adds nothing to keeps its cached texels.
	f32 dist = DotProduct(light->origin, surf->plane->normal)
		- surf->plane->dist;
	f32 reach = light->radius - fabs(dist) - light->minlight + 2;
	if(reach <= 0)
		return 0;
	mtexinfo_t *tex = surf->texinfo;
	vec3_t impact;
	for(s32 i = 0; i < 3; i++)
		impact[i] = light->origin[i] - surf->plane->normal[i] * dist;
	for(s32 i = 0; i < 2; i++){
		f32 l = DotProduct(impact, tex->vecs[i]) + tex->vecs[i][3]
			- surf->texturemins[i];
		if(l <= -reach || l >= surf->extents[i] + reach)
			return 0; // samples run from 0 to extents
	}
	return 1;
}

void R_MarkLights(dlight_t *light, u64 bit, mnode_t *node)
{
	if(node->contents < 0)
		return;
	mplane_t *splitplane = node->plane;
	f32 dist = DotProduct(light->origin, splitplane->normal)
		- splitplane->dist;
	if(dist > light->radius){
		R_MarkLights(light, bit, node->children[0]);
		return;
	}
	if(dist < -light->radius){
		R_MarkLights(light, bit, node->children[1]);
		return;
	}
	// mark the polygons
	msurface_t *surf = cl.worldmodel->surfaces + node->firstsurface;
	for(u32 i = 0; i < node->numsurfaces; i++, surf++){
		if(!R_LightReachesSurf(light, surf))
			continue;
		if(surf->dlightframe != r_dlightframecount){
			surf->dlightbits = 0;
			surf->dlightframe = r_dlightframecount;
		}
		surf->dlightbits |= bit;
	}
	R_MarkLights(light, bit, node->children[0]);
	R_MarkLights(light, bit, node->children[1]);
}

void R_PushDlights()
{
	r_dlightframecount = r_framecount + 1; // because the count hasn't
					       // advanced yet for this frame
	dlight_t *l = cl_dlights;
	for(s32 i = 0; i < MAX_DLIGHTS; i++, l++){
		if(l->die < cl.time || !l->radius)
			continue;
		R_MarkLights(l, 1ULL << i, cl.worldmodel->nodes);
	}
}

s32 RecursiveLightPoint(vec3_t color, mnode_t *node, vec3_t rayorg,
			vec3_t start, vec3_t end, f32 *maxdist)
{
loc0:
	if(node->contents < 0)
		return 0; // didn't hit anything
	f32 front, back;
	if(node->plane->type < 3){ // calculate mid point
		front = start[node->plane->type] - node->plane->dist;
		back = end[node->plane->type] - node->plane->dist;
	}else{
		front = DotProduct(start,node->plane->normal)-node->plane->dist;
		back = DotProduct(end, node->plane->normal) - node->plane->dist;
	}
	if((back < 0) == (front < 0)){ // LordHavoc: optimized recursion
		node = node->children[front < 0];
		goto loc0;
	}
	f32 frac = front / (front-back);
	vec3_t mid;
	mid[0] = start[0] + (end[0] - start[0])*frac;
	mid[1] = start[1] + (end[1] - start[1])*frac;
	mid[2] = start[2] + (end[2] - start[2])*frac;
	// go down front side
	if(RecursiveLightPoint(color, node->children[front < 0],
				rayorg, start, mid, maxdist))
		return 1; // hit something
	else{
		s32 ds, dt;
		msurface_t *surf;
		static vec3_t lightspot;
		VectorCopy(mid, lightspot); // check for impact on this node
		surf = cl.worldmodel->surfaces + node->firstsurface;
		for(u32 i = 0;i < node->numsurfaces;i++, surf++)
		{
			f32 sfront, sback, dist;
			vec3_t raydelta;
			if(surf->flags & SURF_DRAWTILED)
				continue;	// no lightmaps
			//ericw -- added f64 casts to force 64-bit precision.
			//Without them the zombie at the start of jam3_ericw.bsp
			//was incorrectly being lit up in SSE builds.
			ds=(s32)((f64)DotProductF64(mid,surf->texinfo->vecs[0])
					+ surf->texinfo->vecs[0][3]);
			dt=(s32)((f64)DotProductF64(mid,surf->texinfo->vecs[1])
					+ surf->texinfo->vecs[1][3]);
			if(ds<surf->texturemins[0] || dt < surf->texturemins[1])
				continue;
			ds -= surf->texturemins[0];
			dt -= surf->texturemins[1];
			if(ds > surf->extents[0] || dt > surf->extents[1])
				continue;
			if(surf->plane->type < 3){
				sfront = rayorg[surf->plane->type]
						- surf->plane->dist;
				sback = end[surf->plane->type]
						- surf->plane->dist;
			}else{
				sfront = DotProduct(rayorg, surf->plane->normal)
					- surf->plane->dist;
				sback = DotProduct(end, surf->plane->normal)
					- surf->plane->dist;
			}
			VectorSubtract(end, rayorg, raydelta);
			dist = sfront / (sfront-sback) * VectorLength(raydelta);
			if(!surf->samples){
// We hit a surface that is flagged as lightmapped, but doesn't have actual
// lightmap info. Instead of just returning black, we'll keep looking for nearby
// surfaces that do have valid samples. This fixes occasional pitch-black models
// in otherwise well-lit areas in DOTM(e.g. mge1m1, mge4m1) caused by
// overlapping surfaces with mixed lighting data.
				const f32 nearby = 8.f;
				dist += nearby;
				*maxdist = q_min(*maxdist, dist);
				continue;
			}
			if(dist < *maxdist){
				// LordHavoc: enhanced to interpolate lighting
				// WinQuake takes the nearest sample
				s32 lit = lit_loaded && r_rgblighting.value;
				s32 dsfrac = lit ? ds&15 : 0, c[4][3] = {{0}};
				s32 dtfrac = lit ? dt&15 : 0, line3, maps;
				u8 *lightmap;
				line3 = ((surf->extents[0]>>4)+1)*3;
				// LordHavoc: *3 for color
				lightmap = surf->samples + ((dt>>4) *
					((surf->extents[0]>>4)+1) + (ds>>4))*3;
				for(maps = 0;maps < MAXLIGHTMAPS &&
						surf->styles[maps]!=255;maps++){
					s32 scale = d_lightstylevalue[
						surf->styles[maps]];
					for(s32 j = 0; j < 3; j++){
						u8 *r = lightmap + j;
						c[0][j] += r[0] * scale;
						c[1][j] += r[3] * scale;
						c[2][j] += r[line3] * scale;
						c[3][j] += r[line3+3] * scale;
					}
					// LordHavoc: *3 for colored lighting
					lightmap += ((surf->extents[0]>>4)+1) *
						((surf->extents[1]>>4)+1)*3;
				}
				for(s32 j = 0; j < 3; j++){
					s32 r00 = c[0][j]>>8, r01 = c[1][j]>>8;
					s32 r10 = c[2][j]>>8, r11 = c[3][j]>>8;
					s32 top = (((r01-r00)*dsfrac)>>4)+r00;
					s32 bot = (((r11-r10)*dsfrac)>>4)+r10;
					color[j] += (f32)
						((((bot-top)*dtfrac)>>4)+top);
				}
			}
			return 1; // success
		}
		return RecursiveLightPoint(color, node->children[front >= 0],
				rayorg, mid, end, maxdist); // go down back side
	}
}

s32 R_LightPoint(vec3_t p, f32 ofs)
{
	if(!cl.worldmodel->lightdata){
		lightcolor[0] = lightcolor[1] = lightcolor[2] = 255;
		colored_aliaslight = 0;
		return 255;
	}
	f32 maxdist = 8192; //johnfitz -- was 2048
	vec3_t start, end;
	start[0] = p[0];
	start[1] = p[1];
	start[2] = p[2] + ofs;
	end[0] = start[0];
	end[1] = start[1];
	end[2] = start[2] - maxdist;
	lightcolor[0] = lightcolor[1] = lightcolor[2] = 0;
	RecursiveLightPoint(lightcolor, cl.worldmodel->nodes, start, start, end,
								&maxdist);
	if(r_rgblighting.value && !lit_lut_initialized &&
		!(lightcolor[0]==lightcolor[1] && lightcolor[0]==lightcolor[2]))
		R_BuildLitLUT();
	colored_aliaslight = !(lightcolor[0] == lightcolor[1]
			&& lightcolor[0] == lightcolor[2]);
	return((lightcolor[0] + lightcolor[1] + lightcolor[2]) * (1.0f/3.0f));
}
