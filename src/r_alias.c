// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// r_alias.c: routines for setting up to draw alias models
#include "quakedef.h"

typedef struct { s32 index0; s32 index1; } aedge_t;
typedef struct aliaslerp_s {
        f32 v[3];
        u8 lastlightnormal;
        u8 currlightnormal;
        f32 blend;
} aliaslerp_t;

static aliaslerp_t r_aliaslerpverts[MAXALIASVERTS];
static aliaslerp_t *lerpverts;
static mdl_t *pmdl;
static aliashdr_t *paliashdr;
static vec3_t r_plightvec;
static s32 r_ambientlight;
static f32 r_shadelight;
static s32 r_anumverts;
static f32 aliastransform[3][4];
static f32 ziscale;
static vec3_t alias_forward, alias_right, alias_up;

void R_AliasSetUpTransform(s32 trivial_accept);
void R_AliasTransformVector(vec3_t in, vec3_t out);
void R_AliasProjectFinalVert(finalvert_t *fv, auxvert_t *av);


void R_AliasValidateFrame() // Clamp an out-of-range frame (it comes straight
{ // off the network) once, before anything indexes the frame table with it
	entity_t *e = currententity;
	if(e->frame >= e->model->numframes || e->frame < 0){
		Con_DPrintf("No such frame %d %s\n", e->frame, e->model->name);
		e->frame = 0;
	}
}

bool R_AliasCheckBBox()
{
	// expand, rotate, and translate points into worldspace
	currententity->trivial_accept = 0;
	static model_t *pmodel;
	pmodel = currententity->model;
	aliashdr_t *pahdr = Mod_Extradata(pmodel);
	pmdl = (mdl_t *) ((u8 *) pahdr + pahdr->model);
	R_AliasSetUpTransform(0);
	// construct the base bounding box for this frame
	s32 frame = currententity->frame; // validated by R_AliasValidateFrame
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
	for(s32 i = 0; i < 8; i++){
		R_AliasTransformVector(&bp[i][0], &viewaux[i].fv[0]);
		if(viewaux[i].fv[2] < ALIAS_Z_CLIP_PLANE){
		// we must clip points that are closer than the near clip plane
			viewpts[i].flags = ALIAS_Z_CLIP;
			zclipped = 1;
		}else{
			if(viewaux[i].fv[2] < minz)
				minz = viewaux[i].fv[2];
			viewpts[i].flags = 0;
			zfullyclipped = 0;
		}
	}
	if(zfullyclipped)
		return 0; // everything was near-z-clipped
	s32 numv = 8;
	// organize points by edges, use edges to get new points
	for(s32 i = 0; zclipped && i < 12; i++){ // (possible trivial reject)
		static aedge_t aedges[12] = {
			{0,1},{1,2},{2,3},{3,0},
			{4,5},{5,6},{6,7},{7,4},
			{0,5},{1,4},{2,7},{3,6}};
		// edge endpoints
		pv0 = &viewpts[aedges[i].index0];
		pv1 = &viewpts[aedges[i].index1];
		pa0 = &viewaux[aedges[i].index0];
		pa1 = &viewaux[aedges[i].index1];
		// if one end is clipped and the other isn't, make a new point
		if(!(pv0->flags ^ pv1->flags))
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
	for(s32 i = 0; i < numv; i++){
		// we don't need to bother with vertices that were z-clipped
		if(viewpts[i].flags & ALIAS_Z_CLIP)
			continue;
		f32 zi = 1.0 / viewaux[i].fv[2];
		// FIXME: do with chop mode in ASM, or convert to f32
		f32 v0 = (viewaux[i].fv[0] * xscale * zi) + xcenter;
		f32 v1 = (viewaux[i].fv[1] * yscale * zi) + ycenter;
		s32 flags = 0;
		if(v0 < r_refdef.fvrectx)
			flags |= ALIAS_LEFT_CLIP;
		if(v1 < r_refdef.fvrecty)
			flags |= ALIAS_TOP_CLIP;
		if(v0 > r_refdef.fvrectright)
			flags |= ALIAS_RIGHT_CLIP;
		if(v1 > r_refdef.fvrectbottom)
			flags |= ALIAS_BOTTOM_CLIP;
		anyclip |= flags;
		allclip &= flags;
	}
	if(allclip)
		return 0; // trivial reject off one side
	currententity->trivial_accept = !anyclip & !zclipped;
	if(currententity->trivial_accept &&
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

s32 R_AliasLightVert(aliaslerp_t *pverts)
{
        f32 temp = r_ambientlight;
        f32 *plightnormal = r_avertexnormals[pverts->currlightnormal];
        f32 lightcos = DotProduct(plightnormal, r_plightvec) * pverts->blend;
        if(lightcos < 0){
                temp += (s32)(r_shadelight * lightcos);
                // clamp; because we limited the minimum ambient and shading
                // light, wedon't have to clamp low light, just bright
                if(temp < 0) temp = 0;
        }
        plightnormal = r_avertexnormals[pverts->lastlightnormal];
        lightcos = DotProduct(plightnormal, r_plightvec) * (1.0f-pverts->blend);
        if(lightcos < 0){
                temp += (s32)(r_shadelight * lightcos);
                // clamp; because we limited the minimum ambient and shading
                // light, we don't have to clamp low light, just bright
                if(temp < 0)temp = 0;
        }
        return temp;
}

void R_AliasTransformFinalVert(finalvert_t *fv, auxvert_t *av,
				aliaslerp_t *pverts, stvert_t *pstverts){
        av->fv[0] =DotProduct(pverts->v,aliastransform[0])+aliastransform[0][3];
        av->fv[1] =DotProduct(pverts->v,aliastransform[1])+aliastransform[1][3];
        av->fv[2] =DotProduct(pverts->v,aliastransform[2])+aliastransform[2][3];
        fv->v[2] = pstverts->s;
        fv->v[3] = pstverts->t;
        fv->flags = pstverts->onseam;
        fv->v[4] = R_AliasLightVert(pverts);
}

void R_AliasPreparePoints()
{
        stvert_t *pstverts = (stvert_t *)((u8 *)paliashdr + paliashdr->stverts);
        r_anumverts = pmdl->numverts;
        finalvert_t *fv = pfinalverts;
        auxvert_t *av = pauxverts;
        lerpverts = r_aliaslerpverts;
        for(s32 i = 0; i < r_anumverts; i++,fv++,av++,lerpverts++,pstverts++){
                R_AliasTransformFinalVert(fv, av, lerpverts, pstverts);
                if(av->fv[2] < ALIAS_Z_CLIP_PLANE)
                        fv->flags |= ALIAS_Z_CLIP;
                else{
                        R_AliasProjectFinalVert(fv, av);
                        if(fv->v[0] < r_refdef.aliasvrect.x)
                                fv->flags |= ALIAS_LEFT_CLIP;
                        if(fv->v[1] < r_refdef.aliasvrect.y)
                                fv->flags |= ALIAS_TOP_CLIP;
                        if(fv->v[0] > r_refdef.aliasvrectright)
                                fv->flags |= ALIAS_RIGHT_CLIP;
                        if(fv->v[1] > r_refdef.aliasvrectbottom)
                                fv->flags |= ALIAS_BOTTOM_CLIP;
                }
        }
        r_affinetridesc.numtriangles = 1; // clip and draw all triangles
        mtriangle_t *ptri = (mtriangle_t*)((u8*)paliashdr+paliashdr->triangles);
        finalvert_t *pfv[3];
        for(s32 i = 0; i < pmdl->numtris; i++, ptri++){
                pfv[0] = &pfinalverts[ptri->vertindex[0]];
                pfv[1] = &pfinalverts[ptri->vertindex[1]];
                pfv[2] = &pfinalverts[ptri->vertindex[2]];
                if(pfv[0]->flags & pfv[1]->flags & pfv[2]->flags
				& (ALIAS_XY_CLIP_MASK | ALIAS_Z_CLIP))
                        continue; // completely clipped
                if(!((pfv[0]->flags|pfv[1]->flags|pfv[2]->flags)
				& (ALIAS_XY_CLIP_MASK | ALIAS_Z_CLIP))){
                        r_affinetridesc.pfinalverts = pfinalverts;
                        r_affinetridesc.ptriangles = ptri;
                        D_PolysetDraw(); // totally unclipped
                }else{ // partially clipped
                        R_AliasClipTriangle(ptri);
                }
        }
}

void R_AliasSetUpTransform(s32 trivial_accept)
{
	f32 rotationmatrix[3][4], t2matrix[3][4];
	static f32 tmatrix[3][4];
	static f32 viewmatrix[3][4];
	vec3_t angles;
	angles[ROLL] = currententity->angles[ROLL];
	angles[PITCH] = -currententity->angles[PITCH];
	angles[YAW] = currententity->angles[YAW];
	AngleVectors(angles, alias_forward, alias_right, alias_up);
	f32 gunfovscale = 1.0;
	if(currententity == &cl.viewent && 
	    scr_fov.value > 90.0 && cl_gun_fovscale.value > 0.0){
		f32 fov = scr_fov.value;
		if(r_fovmode.value){
			f32 asp = r_fovmode.value==1 ?
				((f32)vid.width / (f32)vid.height * 0.75) :
				yaspectscale.value;
			fov = 2.0f * atan(tanf(scr_fov.value *
					(M_PI/360.0f)) / asp) * (180.0f/M_PI);
		}
		gunfovscale = tanf(fov * (0.5 * M_PI / 180.0));
		gunfovscale = 1.0 + (gunfovscale - 1.0) * cl_gun_fovscale.value;
	}
	tmatrix[0][0] = pmdl->scale[0];
	tmatrix[1][1] = pmdl->scale[1] * gunfovscale;
	tmatrix[2][2] = pmdl->scale[2] * gunfovscale;
	tmatrix[0][3] = pmdl->scale_origin[0];
	tmatrix[1][3] = pmdl->scale_origin[1] * gunfovscale;
	tmatrix[2][3] = pmdl->scale_origin[2] * gunfovscale;
	for(s32 i = 0; i < 3; i++){
		t2matrix[i][0] = alias_forward[i];
		t2matrix[i][1] = -alias_right[i];
		t2matrix[i][2] = alias_up[i];
	}
	t2matrix[0][3] = -modelorg[0];
	t2matrix[1][3] = -modelorg[1];
	t2matrix[2][3] = -modelorg[2];
	R_ConcatTransforms(t2matrix, tmatrix, rotationmatrix);
	VectorCopy(vright, viewmatrix[0]);
	VectorCopy(vup, viewmatrix[1]);
	VectorInverse(viewmatrix[1]);
	VectorCopy(vpn, viewmatrix[2]);
	R_ConcatTransforms(viewmatrix, rotationmatrix, aliastransform);
// do the scaling up of x and y to screen coordinates as part of the transform
// for the unclipped case(it would mess up clipping in the clipped case).
// Also scale down z, so 1/z is scaled 31 bits for free, and scale down x and y
// correspondingly so the projected x and y come out right
	if(!trivial_accept) 
		return;
	for(s32 i = 0; i < 4; i++){
		aliastransform[0][i]*=aliasxscale*(1.0/((f32)0x8000*0x10000));
		aliastransform[1][i]*=aliasyscale*(1.0/((f32)0x8000*0x10000));
		aliastransform[2][i] *= 1.0 / ((f32)0x8000 * 0x10000);
	}
}

void R_AliasProjectFinalVert(finalvert_t *fv, auxvert_t *av)
{
	f32 zi = 1.0 / av->fv[2]; // project points
	fv->v[5] = zi * ziscale;
	fv->v[0] = (av->fv[0] * aliasxscale * zi) + aliasxcenter;
	fv->v[1] = (av->fv[1] * aliasyscale * zi) + aliasycenter;
}

void R_AliasTransformAndProjectFinalVerts(finalvert_t *fv, stvert_t *pstverts)
{
        aliaslerp_t *pverts = r_aliaslerpverts;
        for(s32 i = 0; i < r_anumverts; i++, fv++, pverts++, pstverts++){
                // transform and project
                f32 zi = 1.0 / (DotProduct(pverts->v, aliastransform[2])
				+ aliastransform[2][3]);
                // x, y, and z are scaled down by 1/2**31 in the transform,
		// so 1/z is scaled up by 1/2**31, and the scaling cancels out
		// for x and y in the projection
                fv->v[5] = zi;
                fv->v[0] = ((DotProduct(pverts->v, aliastransform[0])
				+ aliastransform[0][3]) * zi) + aliasxcenter;
                fv->v[1] = ((DotProduct(pverts->v, aliastransform[1])
				+ aliastransform[1][3]) * zi) + aliasycenter;
                fv->v[2] = pstverts->s;
                fv->v[3] = pstverts->t;
                fv->flags = pstverts->onseam;
                fv->v[4] = R_AliasLightVert(pverts);
        }
}

void R_AliasPrepareUnclippedPoints()
{
        stvert_t *pstverts = (stvert_t *)((u8 *)paliashdr + paliashdr->stverts);
        r_anumverts = pmdl->numverts;
        R_AliasTransformAndProjectFinalVerts(pfinalverts, pstverts);
        if(r_affinetridesc.drawtype)
                D_PolysetDrawFinalVerts(pfinalverts, r_anumverts);
        r_affinetridesc.pfinalverts = pfinalverts;
        r_affinetridesc.ptriangles = (mtriangle_t *)
                ((u8 *)paliashdr + paliashdr->triangles);
        r_affinetridesc.numtriangles = pmdl->numtris;
        D_PolysetDraw();
}

void R_AliasSetupSkin()
{
	s32 skinnum = currententity->skinnum;
	if((skinnum >= pmdl->numskins) || (skinnum < 0)){
		Con_DPrintf("R_AliasSetupSkin: no such skin # %d\n", skinnum);
		skinnum = 0;
	}
	static maliasskindesc_t *pskindesc;
	pskindesc = ((maliasskindesc_t *)((u8 *) paliashdr +
				paliashdr->skindesc)) + skinnum;
	static s32 a_skinwidth;
	a_skinwidth = pmdl->skinwidth;
	if(pskindesc->type == ALIAS_SKIN_GROUP){
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
		for(; i < (numskins - 1); i++){
			if(pskinintervals[i] > skintargettime)
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
	if(r_ambientlight < LIGHT_MIN)
		r_ambientlight = LIGHT_MIN;
	r_ambientlight = (255 - r_ambientlight) << VID_CBITS;
	if(r_ambientlight < LIGHT_MIN)
		r_ambientlight = LIGHT_MIN;
	r_shadelight = plighting->shadelight;
	if(r_shadelight < 0)
		r_shadelight = 0;
	r_shadelight *= VID_GRADES;
	// rotate the lighting vector into the model's frame of reference
	r_plightvec[0] = DotProduct(plighting->plightvec, alias_forward);
	r_plightvec[1] = -DotProduct(plighting->plightvec, alias_right);
	r_plightvec[2] = DotProduct(plighting->plightvec, alias_up);
}

void R_SetupAliasFrame(aliashdr_t *paliashdr, s32 frame, lerpdata_t *lerpdata)
{ // johnfitz -- rewritten to support lerping
	entity_t *e = currententity;
	if((frame >= paliashdr->numframes) || (frame < 0)){
		Con_DPrintf("R_AliasSetupFrame: no such frame %d", frame);
		frame = 0;
	}
	s32 posenum = paliashdr->frames[frame].firstpose;
	s32 numposes = paliashdr->frames[frame].numposes;
	if(numposes > 1){
		e->lerptime = paliashdr->frames[frame].interval;
		posenum += (s32)(cl.time / e->lerptime) % numposes;
	}else e->lerptime = 0.1;
	if(e->lerpflags & LERP_RESETANIM){ //kill any lerp in progress
		e->lerpstart = 0;
		e->previouspose = posenum;
		e->currentpose = posenum;
		e->lerpflags -= LERP_RESETANIM;
	}else if(e->currentpose != posenum){ // pose changed, start new lerp
		if(e->lerpflags & LERP_RESETANIM2){//defer lerping one more time
			e->lerpstart = 0;
			e->previouspose = posenum;
			e->currentpose = posenum;
			e->lerpflags -= LERP_RESETANIM2;
		}else{
			e->lerpstart = cl.time;
			e->previouspose = e->currentpose;
			e->currentpose = posenum;
		}
	}
	if(r_lerpmodels.value && !(e->model->flags & MOD_NOLERP
				&& r_lerpmodels.value != 2)){ //set up values
		lerpdata->blend = CLAMP(0, (cl.time - e->lerpstart) / 
			((e->lerpflags&LERP_FINISH && numposes==1) ?
			 (e->lerpfinish - e->lerpstart) : e->lerptime), 1);
		lerpdata->pose1 = e->previouspose;
		lerpdata->pose2 = e->currentpose;
	}else{ //don't lerp
		lerpdata->blend = 1;
		lerpdata->pose1 = posenum;
		lerpdata->pose2 = posenum;
	}
}

void R_BoundPose(aliashdr_t *hdr, lerpdata_t *lerpdata)
{ // the previous pose can come from another frame or even another model
	s32 last = hdr->numposes - 1;
	if(lerpdata->pose2 < 0) lerpdata->pose2 = 0;
	if(lerpdata->pose2 > last) lerpdata->pose2 = last;
	if(lerpdata->pose1 < 0) lerpdata->pose1 = 0;
	if(lerpdata->pose1 > last) lerpdata->pose1 = last;
}

void R_AliasSetupFrameMH(entity_t *ent)
{
	s32 frame = ent->frame; // validated by R_AliasValidateFrame
	lerpdata_t lerpdata;
	R_SetupAliasFrame(paliashdr, frame, &lerpdata);
	f32 blend = lerpdata.pose1 != lerpdata.pose2 ? lerpdata.blend : 0;
	R_BoundPose(paliashdr, &lerpdata);
	s32 *poseofs = (s32 *)((u8 *)paliashdr + paliashdr->posedata);
	trivertx_t *currv = (trivertx_t *)((u8 *)paliashdr +
			poseofs[lerpdata.pose2]);
	trivertx_t *lastv = (trivertx_t *)((u8 *)paliashdr +
			poseofs[lerpdata.pose1]);
	lerpverts = r_aliaslerpverts;
	for(s32 i = 0; i < pmdl->numverts; i++, currv++, lastv++, lerpverts++){
		lerpverts->v[0] = currv->v[0]*blend + lastv->v[0]*(1.0f-blend);
		lerpverts->v[1] = currv->v[1]*blend + lastv->v[1]*(1.0f-blend);
		lerpverts->v[2] = currv->v[2]*blend + lastv->v[2]*(1.0f-blend);
		lerpverts->currlightnormal = currv->lightnormalindex;
		lerpverts->lastlightnormal = lastv->lightnormalindex;
		lerpverts->blend = blend;
	}
}

void R_AliasDrawModel(alight_t *plighting)
{
	finalvert_t finalverts[MAXALIASVERTS];
	auxvert_t auxverts[MAXALIASVERTS];
	r_amodels_drawn++;
	pfinalverts = (finalvert_t *)(((uintptr_t) & finalverts[0] //cache align
				+ CACHE_SIZE - 1) & ~(CACHE_SIZE - 1));
	pauxverts = &auxverts[0];
	paliashdr = (aliashdr_t *)Mod_Extradata(currententity->model);
	pmdl = (mdl_t *)((u8 *) paliashdr + paliashdr->model);
	R_AliasSetupSkin();
	R_AliasSetUpTransform(currententity->trivial_accept);
	R_AliasSetupLighting(plighting);
	R_AliasSetupFrameMH(currententity);
	if(!currententity->colormap)
		Sys_Error("R_AliasDrawModel: !currententity->colormap");
	r_affinetridesc.drawtype = currententity->trivial_accept == 3;
	if(r_affinetridesc.drawtype)
		D_PolysetUpdateTables(); // FIXME: precalc...
	acolormap = currententity->colormap;
	ziscale = currententity != &cl.viewent ? (f32)0x8000 * (f32)0x10000
		: (f32)0x8000 * (f32)0x10000 * 3.0;
	if(currententity->trivial_accept)
		R_AliasPrepareUnclippedPoints();
	else
		R_AliasPreparePoints();
}

// Entity occlusion against the z-buffer (r_entocclusion)
// If every pixel of an entity's screen footprint is already covered by 
// something nearer than the nearest point the entity can have, then the alias
// rasterizer's own per-pixel test would reject every pixel of it, and we can
// skip the lighting, vertex transform and triangle setup altogether.

// Safety margins. The rasterizer works in fixed point and may touch a pixel
// or two outside the projected triangle, so the nearest depth is inflated and
// the rectangle padded. Both only ever make the test refuse to cull.
#define ENTOCC_ZMARGIN 1.0625f // nearest 1/z of the bounds, +6.25%
#define ENTOCC_ZSLACK 2 // plus two z-buffer units
#define ENTOCC_PAD 2 // pixels

// Is every z-buffer value in this (inclusive) rectangle strictly greater than
// zbound? The rasterizer draws a pixel when (zi >> 16) >= *pz, so if all
// stored values are > zbound >= any pixel it could produce, nothing is drawn.
static bool R_ZBufferRectBehind(s32 x0, s32 y0, s32 x1, s32 y1, s32 zbound)
{
	for(s32 y = y0; y <= y1; y++){
		const s16 *row = d_pzbuffer + (size_t)y * d_zwidth;
		for(s32 x = x0; x <= x1;){
			s32 n = x1 - x + 1;
			if(n > 64) n = 64;
			s16 m = row[x];
			for(s32 i = 1; i < n; i++)
				m = row[x + i] < m ? row[x + i] : m;
			if(m <= zbound)
				return 0;
			x += n;
		}
	}
	return 1;
}

// Screen rectangle (inclusive, clipped to the view) the entity can touch, and
// the largest z-buffer value any of its pixels can produce. False when we
// can't bound it (it reaches the near plane, etc.), meaning "don't cull".
static bool R_AliasScreenBounds(s32 *rx0, s32 *ry0, s32 *rx1, s32 *ry1,
		s32 *rzbound)
{
	entity_t *e = currententity;
	model_t *m = e->model;
	const f32 *mn, *mx; // bounds over every vertex of every pose, so they
			    // also hold for interpolated and grouped frames
	if(e->angles[0] || e->angles[2]){ // pitch or roll
		mn = m->rmins;
		mx = m->rmaxs;
	}else if(e->angles[1]){ // yaw
		mn = m->ymins;
		mx = m->ymaxs;
	}else{
		mn = m->mins;
		mx = m->maxs;
	}
	f32 zmin = FLT_MAX;
	f32 minx = FLT_MAX, maxx = -FLT_MAX, miny = FLT_MAX, maxy = -FLT_MAX;
	f32 tv[8][3];
	for(s32 i = 0; i < 8; i++){
		vec3_t corner, local;
		for(s32 k = 0; k < 3; k++)
			corner[k] = e->origin[k] + ((i>>k) & 1 ? mx[k] : mn[k]);
		VectorSubtract(corner, r_origin, local);
		TransformVector(local, tv[i]);
		if(tv[i][2] <= ALIAS_Z_CLIP_PLANE)
			return 0; // reaches the near plane, don't try to clip
		if(tv[i][2] < zmin)
			zmin = tv[i][2];
	}
	for(s32 i = 0; i < 8; i++){
		f32 zi = 1.0f / tv[i][2];
		f32 x = xcenter + xscale * zi * tv[i][0];
		f32 y = ycenter - yscale * zi * tv[i][1];
		if(x < minx) minx = x;
		if(x > maxx) maxx = x;
		if(y < miny) miny = y;
		if(y > maxy) maxy = y;
	}
	if(!(minx > -1e6f && maxx < 1e6f && miny > -1e6f && maxy < 1e6f))
		return 0; // absurd projection, let the normal path deal with it
	// the footprint, padded and clipped to the view
	s32 vx0 = r_refdef.vrect.x, vy0 = r_refdef.vrect.y;
	s32 vx1 = vx0 + r_refdef.vrect.width - 1;
	s32 vy1 = vy0 + r_refdef.vrect.height - 1;
	s32 x0 = (s32)floorf(minx)-ENTOCC_PAD, x1 = (s32)ceilf(maxx)+ENTOCC_PAD;
	s32 y0 = (s32)floorf(miny)-ENTOCC_PAD, y1 = (s32)ceilf(maxy)+ENTOCC_PAD;
	if(x0 < vx0) x0 = vx0;
	if(y0 < vy0) y0 = vy0;
	if(x1 > vx1) x1 = vx1;
	if(y1 > vy1) y1 = vy1;
	if(x0 > x1 || y0 > y1)
		return 0;
	// the stored value is (s16)(zi * 0x8000), same scale as the rasterizer
	f32 zbound = (1.0f / zmin) * 32768.0f * ENTOCC_ZMARGIN + ENTOCC_ZSLACK;
	if(zbound >= 32767.0f)
		return 0; // too close for the s16 buffer to prove anything
	*rx0 = x0; *ry0 = y0; *rx1 = x1; *ry1 = y1;
	*rzbound = (s32)zbound;
	return 1;
}

bool R_AliasOccluded()
{
	s32 x0, y0, x1, y1, zbound;
	if(!d_pzbuffer || !R_AliasScreenBounds(&x0, &y0, &x1, &y1, &zbound))
		return 0;
	return R_ZBufferRectBehind(x0, y0, x1, y1, zbound);
}

// Drawing an alias entity advances its animation lerp state as a side effect
// (previouspose/currentpose/lerpstart). An entity we skip must still advance,
// or it would blend from a stale pose when it comes back into view.
void R_AliasSyncLerpState()
{
	lerpdata_t lerpdata;
	R_SetupAliasFrame(Mod_Extradata(currententity->model),
			currententity->frame, &lerpdata);
}
