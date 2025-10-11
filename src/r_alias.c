// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// r_alias.c: routines for setting up to draw alias models
#include "quakedef.h"

mtriangle_t *ptriangles;
affinetridesc_t r_affinetridesc;
void *acolormap;
mdl_t *pmdl;
aliashdr_t *paliashdr;
finalvert_t *pfinalverts;
auxvert_t *pauxverts;
s32 r_amodels_drawn;
s32 a_skinwidth;

static trivertx_t *r_apverts;
static vec3_t r_plightvec;
static s32 r_ambientlight;
static f32 r_shadelight;
static s32 r_anumverts;
static f32 aliastransform[3][4];
static f32 ziscale;
static model_t *pmodel;
static vec3_t alias_forward, alias_right, alias_up;
static maliasskindesc_t *pskindesc;
static aedge_t aedges[12] = {
	{0,1},{1,2},{2,3},{3,0},
	{4,5},{5,6},{6,7},{7,4},
	{0,5},{1,4},{2,7},{3,6}};

void R_AliasTransformAndProjectFinalVerts(finalvert_t *fv, stvert_t *pstverts);
void R_AliasSetUpTransform(s32 trivial_accept);
void R_AliasTransformVector(vec3_t in, vec3_t out);
void R_AliasTransformFinalVert(finalvert_t *fv, auxvert_t *av,
				trivertx_t *pverts, stvert_t *pstverts);
void R_AliasProjectFinalVert(finalvert_t *fv, auxvert_t *av);

bool R_AliasCheckBBox()
{
	// expand, rotate, and translate points into worldspace
	currententity->trivial_accept = 0;
	pmodel = currententity->model;
	aliashdr_t *pahdr = Mod_Extradata(pmodel);
	pmdl = (mdl_t *) ((u8 *) pahdr + pahdr->model);
	R_AliasSetUpTransform(0);
	// construct the base bounding box for this frame
	s32 frame = currententity->frame;
	// TODO: don't repeat this check when drawing?
	if ((frame >= pmdl->numframes) || (frame < 0)) {
		Con_DPrintf("No such frame %d %s\n", frame, pmodel->name);
		frame = 0;
	}
	maliasframedesc_t *pframedesc = &pahdr->frames[frame];
	f32 bp[8][3]; // basepts
	// x worldspace coordinates
	bp[0][0]=bp[1][0]=bp[2][0]=bp[3][0]=(f32)pframedesc->bboxmin.v[0];
	bp[4][0]=bp[5][0]=bp[6][0]=bp[7][0]=(f32)pframedesc->bboxmax.v[0];
	// y worldspace coordinates
	bp[0][1]=bp[3][1]=bp[5][1]=bp[6][1]=(f32)pframedesc->bboxmin.v[1];
	bp[1][1]=bp[2][1]=bp[4][1]=bp[7][1]=(f32)pframedesc->bboxmax.v[1];
	// z worldspace coordinates
	bp[0][2]=bp[1][2]=bp[4][2]=bp[5][2]=(f32)pframedesc->bboxmin.v[2];
	bp[2][2]=bp[3][2]=bp[6][2]=bp[7][2]=(f32)pframedesc->bboxmax.v[2];
	bool zclipped = 0;
	bool zfullyclipped = 1;
	s32 minz = 9999;
	finalvert_t *pv0, *pv1, viewpts[16];
	auxvert_t *pa0, *pa1, viewaux[16];
	for (s32 i = 0; i < 8; i++) {
		R_AliasTransformVector(&bp[i][0], &viewaux[i].fv[0]);
		if (viewaux[i].fv[2] < ALIAS_Z_CLIP_PLANE) {
		// we must clip points that are closer than the near clip plane
			viewpts[i].flags = ALIAS_Z_CLIP;
			zclipped = 1;
		} else {
			if (viewaux[i].fv[2] < minz)
				minz = viewaux[i].fv[2];
			viewpts[i].flags = 0;
			zfullyclipped = 0;
		}
	}
	if (zfullyclipped)
		return 0; // everything was near-z-clipped
	s32 numv = 8;
	// organize points by edges, use edges to get new points
	for (s32 i = 0; zclipped && i < 12; i++) { // (possible trivial reject)
		// edge endpoints
		pv0 = &viewpts[aedges[i].index0];
		pv1 = &viewpts[aedges[i].index1];
		pa0 = &viewaux[aedges[i].index0];
		pa1 = &viewaux[aedges[i].index1];
		// if one end is clipped and the other isn't, make a new point
		if (!(pv0->flags ^ pv1->flags))
			continue;
		f32 f=(ALIAS_Z_CLIP_PLANE-pa0->fv[2])/(pa1->fv[2]-pa0->fv[2]);
		viewaux[numv].fv[0] = pa0->fv[0]+(pa1->fv[0]-pa0->fv[0])*f;
		viewaux[numv].fv[1] = pa0->fv[1]+(pa1->fv[1]-pa0->fv[1])*f;
		viewaux[numv].fv[2] = ALIAS_Z_CLIP_PLANE;
		viewpts[numv].flags = 0;
		numv++;
	}
	// project the vertices that remain after clipping
	u32 anyclip = 0;
	u32 allclip = ALIAS_XY_CLIP_MASK;
	for (s32 i = 0; i < numv; i++) {
		// we don't need to bother with vertices that were z-clipped
		if (viewpts[i].flags & ALIAS_Z_CLIP)
			continue;
		f32 zi = 1.0 / viewaux[i].fv[2];
		// FIXME: do with chop mode in ASM, or convert to f32
		f32 v0 = (viewaux[i].fv[0] * xscale * zi) + xcenter;
		f32 v1 = (viewaux[i].fv[1] * yscale * zi) + ycenter;
		s32 flags = 0;
		if (v0 < r_refdef.fvrectx)
			flags |= ALIAS_LEFT_CLIP;
		if (v1 < r_refdef.fvrecty)
			flags |= ALIAS_TOP_CLIP;
		if (v0 > r_refdef.fvrectright)
			flags |= ALIAS_RIGHT_CLIP;
		if (v1 > r_refdef.fvrectbottom)
			flags |= ALIAS_BOTTOM_CLIP;
		anyclip |= flags;
		allclip &= flags;
	}
	if (allclip)
		return 0; // trivial reject off one side
	currententity->trivial_accept = !anyclip & !zclipped;
	if (currententity->trivial_accept &&
		minz > (r_aliastransition + (pmdl->size * r_resfudge)))
			currententity->trivial_accept |= 2;
	return 1;
}

void R_AliasTransformVector(vec3_t in, vec3_t out)
{
	out[0] = DotProduct(in, aliastransform[0]) + aliastransform[0][3];
	out[1] = DotProduct(in, aliastransform[1]) + aliastransform[1][3];
	out[2] = DotProduct(in, aliastransform[2]) + aliastransform[2][3];
}

void R_AliasPreparePoints()
{ // General clipped case
	stvert_t *pstverts = (stvert_t*)((u8*)paliashdr + paliashdr->stverts);
	r_anumverts = pmdl->numverts;
	finalvert_t *fv = pfinalverts;
	auxvert_t *av = pauxverts;
	for (s32 i = 0; i < r_anumverts; i++,fv++,av++,r_apverts++,pstverts++) {
		R_AliasTransformFinalVert(fv, av, r_apverts, pstverts);
		if (av->fv[2] < ALIAS_Z_CLIP_PLANE)
			fv->flags |= ALIAS_Z_CLIP;
		else {
			R_AliasProjectFinalVert(fv, av);
			if (fv->v[0] < r_refdef.aliasvrect.x)
				fv->flags |= ALIAS_LEFT_CLIP;
			if (fv->v[1] < r_refdef.aliasvrect.y)
				fv->flags |= ALIAS_TOP_CLIP;
			if (fv->v[0] > r_refdef.aliasvrectright)
				fv->flags |= ALIAS_RIGHT_CLIP;
			if (fv->v[1] > r_refdef.aliasvrectbottom)
				fv->flags |= ALIAS_BOTTOM_CLIP;
		}
	}
	r_affinetridesc.numtriangles = 1; // clip and draw all triangles
	mtriangle_t *ptri = (mtriangle_t *) ((u8 *) paliashdr + paliashdr->triangles);
	finalvert_t *pfv[3];
	for (s32 i = 0; i < pmdl->numtris; i++, ptri++) {
		pfv[0] = &pfinalverts[ptri->vertindex[0]];
		pfv[1] = &pfinalverts[ptri->vertindex[1]];
		pfv[2] = &pfinalverts[ptri->vertindex[2]];
		if (pfv[0]->flags & pfv[1]-> flags & pfv[2]->flags &
				(ALIAS_XY_CLIP_MASK | ALIAS_Z_CLIP))
			continue; // completely clipped
		if (!((pfv[0]->flags | pfv[1]->flags | pfv[2]->flags) &
			(ALIAS_XY_CLIP_MASK | ALIAS_Z_CLIP))) { // unclipped
			r_affinetridesc.pfinalverts = pfinalverts;
			r_affinetridesc.ptriangles = ptri;
			D_PolysetDraw();
		} else // partially clipped
			R_AliasClipTriangle(ptri);
	}
}

void R_AliasSetUpTransform(s32 trivial_accept)
{
	f32 rotationmatrix[3][4], t2matrix[3][4];
	static f32 tmatrix[3][4];
	static f32 viewmatrix[3][4];
	vec3_t angles;
	// TODO: should really be stored with the entity instead of being reconstructed
	// TODO: should use a look-up table
	// TODO: could cache lazily, stored in the entity
	angles[ROLL] = currententity->angles[ROLL];
	angles[PITCH] = -currententity->angles[PITCH];
	angles[YAW] = currententity->angles[YAW];
	AngleVectors(angles, alias_forward, alias_right, alias_up);
	tmatrix[0][0] = pmdl->scale[0];
	tmatrix[1][1] = pmdl->scale[1];
	tmatrix[2][2] = pmdl->scale[2];
	tmatrix[0][3] = pmdl->scale_origin[0];
	tmatrix[1][3] = pmdl->scale_origin[1];
	tmatrix[2][3] = pmdl->scale_origin[2];
	// TODO: can do this with simple matrix rearrangement
	for (s32 i = 0; i < 3; i++) {
		t2matrix[i][0] = alias_forward[i];
		t2matrix[i][1] = -alias_right[i];
		t2matrix[i][2] = alias_up[i];
	}
	t2matrix[0][3] = -modelorg[0];
	t2matrix[1][3] = -modelorg[1];
	t2matrix[2][3] = -modelorg[2];
	// FIXME: can do more efficiently than full concatenation
	R_ConcatTransforms(t2matrix, tmatrix, rotationmatrix);
	// TODO: should be global, set when vright, etc., set
	VectorCopy(vright, viewmatrix[0]);
	VectorCopy(vup, viewmatrix[1]);
	VectorInverse(viewmatrix[1]);
	VectorCopy(vpn, viewmatrix[2]);
	R_ConcatTransforms(viewmatrix, rotationmatrix, aliastransform);
	// do the scaling up of x and y to screen coordinates as part of the transform
	// for the unclipped case (it would mess up clipping in the clipped case).
	// Also scale down z, so 1/z is scaled 31 bits for free, and scale down x and y
	// correspondingly so the projected x and y come out right
	// FIXME: make this work for clipped case too?
	if (!trivial_accept) 
		return;
	for (s32 i = 0; i < 4; i++) {
		aliastransform[0][i]*=aliasxscale*(1.0/((f32)0x8000*0x10000));
		aliastransform[1][i]*=aliasyscale*(1.0/((f32)0x8000*0x10000));
		aliastransform[2][i] *= 1.0 / ((f32)0x8000 * 0x10000);
	}
}

void R_AliasTransformFinalVert(finalvert_t *fv, auxvert_t *av,
		trivertx_t *pverts, stvert_t *pstverts)
{
	av->fv[0]=DotProductU8(pverts->v,aliastransform[0])+aliastransform[0][3];
	av->fv[1]=DotProductU8(pverts->v,aliastransform[1])+aliastransform[1][3];
	av->fv[2]=DotProductU8(pverts->v,aliastransform[2])+aliastransform[2][3];
	fv->v[2] = pstverts->s;
	fv->v[3] = pstverts->t;
	fv->flags = pstverts->onseam;
	// lighting
	f32 *plightnormal = r_avertexnormals[pverts->lightnormalindex];
	f32 lightcos = DotProduct(plightnormal, r_plightvec);
	s32 temp = r_ambientlight;
	if (lightcos < 0) {
		temp += (s32)(r_shadelight * lightcos);
		// because we limited the minimum ambient and shading light, we
		// don't have to clamp low light, just bright
		if (temp < 0)
			temp = 0;
	}
	fv->v[4] = temp;
}

void R_AliasTransformAndProjectFinalVerts(finalvert_t *fv, stvert_t *pstverts)
{
	trivertx_t *pverts = r_apverts;
	for (s32 i = 0; i < r_anumverts; i++, fv++, pverts++, pstverts++) {
		// transform and project
		f32 zi = 1.0 / (DotProductU8(pverts->v, aliastransform[2]) + 
				aliastransform[2][3]);
		// x, y, and z are scaled down by 1/2**31 in the transform, so
		// 1/z is scaled up by 1/2**31, and the scaling cancels out for
		// x and y in the projection
		fv->v[5] = zi;
		fv->v[0] = ((DotProductU8(pverts->v, aliastransform[0])
				+ aliastransform[0][3]) * zi) + aliasxcenter;
		fv->v[1] = ((DotProductU8(pverts->v, aliastransform[1])
				+ aliastransform[1][3]) * zi) + aliasycenter;
		fv->v[2] = pstverts->s;
		fv->v[3] = pstverts->t;
		fv->flags = pstverts->onseam;
		// lighting
		f32 *plightnormal=r_avertexnormals[pverts->lightnormalindex];
		f32 lightcos = DotProduct(plightnormal, r_plightvec);
		s32 temp = r_ambientlight;
		if (lightcos < 0) {
			temp += (s32)(r_shadelight * lightcos);
			// because we limited the minimum ambient and shading
			// light, we don't have to clamp low light, just bright
			if (temp < 0)
				temp = 0;
		}
		fv->v[4] = temp;
	}
}

void R_AliasProjectFinalVert(finalvert_t *fv, auxvert_t *av)
{
	f32 zi = 1.0 / av->fv[2]; // project points
	fv->v[5] = zi * ziscale;
	fv->v[0] = (av->fv[0] * aliasxscale * zi) + aliasxcenter;
	fv->v[1] = (av->fv[1] * aliasyscale * zi) + aliasycenter;
}

void R_AliasPrepareUnclippedPoints()
{
	stvert_t *pstverts = (stvert_t*)((u8*)paliashdr + paliashdr->stverts);
	r_anumverts = pmdl->numverts;
	// FIXME: just use pfinalverts directly?
	finalvert_t *fv = pfinalverts;
	R_AliasTransformAndProjectFinalVerts(fv, pstverts);
	if (r_affinetridesc.drawtype)
		D_PolysetDrawFinalVerts(fv, r_anumverts);
	r_affinetridesc.pfinalverts = pfinalverts;
	r_affinetridesc.ptriangles = (mtriangle_t *)
		((u8 *) paliashdr + paliashdr->triangles);
	r_affinetridesc.numtriangles = pmdl->numtris;
	D_PolysetDraw();
}

void R_AliasSetupSkin()
{
	s32 skinnum = currententity->skinnum;
	if ((skinnum >= pmdl->numskins) || (skinnum < 0)) {
		Con_DPrintf("R_AliasSetupSkin: no such skin # %d\n", skinnum);
		skinnum = 0;
	}
	pskindesc = ((maliasskindesc_t *)((u8 *) paliashdr +
				paliashdr->skindesc)) + skinnum;
	a_skinwidth = pmdl->skinwidth;
	if (pskindesc->type == ALIAS_SKIN_GROUP) {
		maliasskingroup_t *paliasskingroup = (maliasskingroup_t *)
				((u8 *) paliashdr + pskindesc->skin);
		f32 *pskinintervals = (f32 *)((u8 *) paliashdr
				+ paliasskingroup->intervals);
		s32 numskins = paliasskingroup->numskins;
		f32 fullskininterval = pskinintervals[numskins - 1];
		f32 skintime = cl.time + currententity->syncbase;
		// when loading in Mod_LoadAliasSkinGroup, we guaranteed all interval
		// values are positive, so we don't have to worry about division by 0
		f32 skintargettime = skintime - ((s32)(skintime /
					fullskininterval)) * fullskininterval;
		s32 i = 0;
		for (; i < (numskins - 1); i++) {
			if (pskinintervals[i] > skintargettime)
				break;
		}
		pskindesc = &paliasskingroup->skindescs[i];
	}
	r_affinetridesc.pskindesc = pskindesc;
	r_affinetridesc.pskin = (void *)((u8 *) paliashdr + pskindesc->skin);
	r_affinetridesc.skinwidth = a_skinwidth;
	r_affinetridesc.seamfixupX16 = (a_skinwidth >> 1) << 16;
	r_affinetridesc.skinheight = pmdl->skinheight;
}

void R_AliasSetupLighting(alight_t *plighting)
{
	// guarantee that no vertex will ever be lit below LIGHT_MIN, so we
	// don't have to clamp off the bottom
	r_ambientlight = plighting->ambientlight;
	if (r_ambientlight < LIGHT_MIN)
		r_ambientlight = LIGHT_MIN;
	r_ambientlight = (255 - r_ambientlight) << VID_CBITS;
	if (r_ambientlight < LIGHT_MIN)
		r_ambientlight = LIGHT_MIN;
	r_shadelight = plighting->shadelight;
	if (r_shadelight < 0)
		r_shadelight = 0;
	r_shadelight *= VID_GRADES;
	// rotate the lighting vector into the model's frame of reference
	r_plightvec[0] = DotProduct(plighting->plightvec, alias_forward);
	r_plightvec[1] = -DotProduct(plighting->plightvec, alias_right);
	r_plightvec[2] = DotProduct(plighting->plightvec, alias_up);
}

void R_AliasSetupFrame()
{ // set r_apverts
	s32 frame = currententity->frame;
	if ((frame >= pmdl->numframes) || (frame < 0)) {
		Con_DPrintf("R_AliasSetupFrame: no such frame %d\n", frame);
		frame = 0;
	}
	if (paliashdr->frames[frame].type == ALIAS_SINGLE) {
		r_apverts = (trivertx_t *)
			((u8 *) paliashdr + paliashdr->frames[frame].frame);
		return;
	}
	maliasgroup_t *paliasgroup = (maliasgroup_t *) ((u8 *) paliashdr +
			paliashdr->frames[frame].frame);
	f32 *pintervals = (f32 *)((u8 *)paliashdr+paliasgroup->intervals);
	s32 numframes = paliasgroup->numframes;
	f32 fullinterval = pintervals[numframes - 1];
	f32 time = cl.time + currententity->syncbase;
	// when loading in Mod_LoadAliasGroup, we guaranteed all interval values
	// are positive, so we don't have to worry about division by 0
	f32 targettime = time - ((s32)(time / fullinterval)) * fullinterval;
	s32 i = 0;
	for (; i < (numframes - 1); i++)
		if (pintervals[i] > targettime)
			break;
	r_apverts=(trivertx_t*)((u8*)paliashdr+paliasgroup->frames[i].frame);
}

void R_AliasDrawModel(alight_t *plighting)
{
	finalvert_t finalverts[MAXALIASVERTS];
	auxvert_t auxverts[MAXALIASVERTS];
	r_amodels_drawn++;
	// cache align
	pfinalverts = (finalvert_t *) (((uintptr_t) & finalverts[0]
				+ CACHE_SIZE - 1) & ~(CACHE_SIZE - 1));
	pauxverts = &auxverts[0];
	paliashdr = (aliashdr_t *) Mod_Extradata(currententity->model);
	pmdl = (mdl_t *) ((u8 *) paliashdr + paliashdr->model);
	R_AliasSetupSkin();
	R_AliasSetUpTransform(currententity->trivial_accept);
	R_AliasSetupLighting(plighting);
	R_AliasSetupFrame();
	if (!currententity->colormap)
		Sys_Error("R_AliasDrawModel: !currententity->colormap");
	r_affinetridesc.drawtype = currententity->trivial_accept == 3;
	if (r_affinetridesc.drawtype)
		D_PolysetUpdateTables(); // FIXME: precalc...
	acolormap = currententity->colormap;
	ziscale = currententity != &cl.viewent ? (f32)0x8000 * (f32)0x10000
		: (f32)0x8000 * (f32)0x10000 * 3.0;
	if (currententity->trivial_accept)
		R_AliasPrepareUnclippedPoints();
	else
		R_AliasPreparePoints();
}
