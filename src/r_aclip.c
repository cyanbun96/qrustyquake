// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// r_aclip.c: clip routines for drawing Alias models directly to the screen
#include "quakedef.h"

static finalvert_t fv[2][8];
static auxvert_t av[8];

void R_AliasProjectFinalVert(finalvert_t *fv, auxvert_t *av);
void R_Alias_clip_top(finalvert_t *pfv0, finalvert_t *pfv1,finalvert_t *out);
void R_Alias_clip_bottom(finalvert_t *pfv0, finalvert_t *pfv1,finalvert_t *out);
void R_Alias_clip_left(finalvert_t *pfv0, finalvert_t *pfv1, finalvert_t * out);
void R_Alias_clip_right(finalvert_t *pfv0, finalvert_t *pfv1, finalvert_t *out);

void R_Alias_clip_z(finalvert_t *pfv0, finalvert_t *pfv1, finalvert_t *out)
{ // pfv0 is the unclipped vertex, pfv1 is the z-clipped vertex
	f32 scale;
	auxvert_t *pav0, *pav1, avout;

	pav0 = &av[pfv0 - &fv[0][0]];
	pav1 = &av[pfv1 - &fv[0][0]];

	if (pfv0->v[1] >= pfv1->v[1]) {
		scale = (ALIAS_Z_CLIP_PLANE - pav0->fv[2]) /
		    (pav1->fv[2] - pav0->fv[2]);

		avout.fv[0] = pav0->fv[0] + (pav1->fv[0] - pav0->fv[0]) * scale;
		avout.fv[1] = pav0->fv[1] + (pav1->fv[1] - pav0->fv[1]) * scale;
		avout.fv[2] = ALIAS_Z_CLIP_PLANE;

		out->v[2] = pfv0->v[2] + (pfv1->v[2] - pfv0->v[2]) * scale;
		out->v[3] = pfv0->v[3] + (pfv1->v[3] - pfv0->v[3]) * scale;
		out->v[4] = pfv0->v[4] + (pfv1->v[4] - pfv0->v[4]) * scale;
	} else {
		scale = (ALIAS_Z_CLIP_PLANE - pav1->fv[2]) /
		    (pav0->fv[2] - pav1->fv[2]);

		avout.fv[0] = pav1->fv[0] + (pav0->fv[0] - pav1->fv[0]) * scale;
		avout.fv[1] = pav1->fv[1] + (pav0->fv[1] - pav1->fv[1]) * scale;
		avout.fv[2] = ALIAS_Z_CLIP_PLANE;

		out->v[2] = pfv1->v[2] + (pfv0->v[2] - pfv1->v[2]) * scale;
		out->v[3] = pfv1->v[3] + (pfv0->v[3] - pfv1->v[3]) * scale;
		out->v[4] = pfv1->v[4] + (pfv0->v[4] - pfv1->v[4]) * scale;
	}

	R_AliasProjectFinalVert(out, &avout);

	if (out->v[0] < r_refdef.aliasvrect.x)
		out->flags |= ALIAS_LEFT_CLIP;
	if (out->v[1] < r_refdef.aliasvrect.y)
		out->flags |= ALIAS_TOP_CLIP;
	if (out->v[0] > r_refdef.aliasvrectright)
		out->flags |= ALIAS_RIGHT_CLIP;
	if (out->v[1] > r_refdef.aliasvrectbottom)
		out->flags |= ALIAS_BOTTOM_CLIP;
}

void R_Alias_clip_left(finalvert_t *pfv0, finalvert_t *pfv1, finalvert_t *out)
{
	f32 scale;
	s32 i;

	if (pfv0->v[1] >= pfv1->v[1]) {
		scale = (f32)(r_refdef.aliasvrect.x - pfv0->v[0]) /
		    (pfv1->v[0] - pfv0->v[0]);
		for (i = 0; i < 6; i++)
			out->v[i] =
			    pfv0->v[i] + (pfv1->v[i] - pfv0->v[i]) * scale +
			    0.5;
	} else {
		scale = (f32)(r_refdef.aliasvrect.x - pfv1->v[0]) /
		    (pfv0->v[0] - pfv1->v[0]);
		for (i = 0; i < 6; i++)
			out->v[i] =
			    pfv1->v[i] + (pfv0->v[i] - pfv1->v[i]) * scale +
			    0.5;
	}
}

void R_Alias_clip_right(finalvert_t *pfv0, finalvert_t *pfv1, finalvert_t *out)
{
	f32 scale;
	s32 i;

	if (pfv0->v[1] >= pfv1->v[1]) {
		scale = (f32)(r_refdef.aliasvrectright - pfv0->v[0]) /
		    (pfv1->v[0] - pfv0->v[0]);
		for (i = 0; i < 6; i++)
			out->v[i] =
			    pfv0->v[i] + (pfv1->v[i] - pfv0->v[i]) * scale +
			    0.5;
	} else {
		scale = (f32)(r_refdef.aliasvrectright - pfv1->v[0]) /
		    (pfv0->v[0] - pfv1->v[0]);
		for (i = 0; i < 6; i++)
			out->v[i] =
			    pfv1->v[i] + (pfv0->v[i] - pfv1->v[i]) * scale +
			    0.5;
	}
}

void R_Alias_clip_top(finalvert_t *pfv0, finalvert_t *pfv1, finalvert_t *out)
{
	f32 scale;
	s32 i;

	if (pfv0->v[1] >= pfv1->v[1]) {
		scale = (f32)(r_refdef.aliasvrect.y - pfv0->v[1]) /
		    (pfv1->v[1] - pfv0->v[1]);
		for (i = 0; i < 6; i++)
			out->v[i] =
			    pfv0->v[i] + (pfv1->v[i] - pfv0->v[i]) * scale +
			    0.5;
	} else {
		scale = (f32)(r_refdef.aliasvrect.y - pfv1->v[1]) /
		    (pfv0->v[1] - pfv1->v[1]);
		for (i = 0; i < 6; i++)
			out->v[i] =
			    pfv1->v[i] + (pfv0->v[i] - pfv1->v[i]) * scale +
			    0.5;
	}
}

void R_Alias_clip_bottom(finalvert_t *pfv0, finalvert_t *pfv1, finalvert_t *out)
{
	f32 scale;
	s32 i;

	if (pfv0->v[1] >= pfv1->v[1]) {
		scale = (f32)(r_refdef.aliasvrectbottom - pfv0->v[1]) /
		    (pfv1->v[1] - pfv0->v[1]);

		for (i = 0; i < 6; i++)
			out->v[i] =
			    pfv0->v[i] + (pfv1->v[i] - pfv0->v[i]) * scale +
			    0.5;
	} else {
		scale = (f32)(r_refdef.aliasvrectbottom - pfv1->v[1]) /
		    (pfv0->v[1] - pfv1->v[1]);

		for (i = 0; i < 6; i++)
			out->v[i] =
			    pfv1->v[i] + (pfv0->v[i] - pfv1->v[i]) * scale +
			    0.5;
	}
}

s32 R_AliasClip(finalvert_t *in, finalvert_t *out, s32 flag, s32 count,
	void (*clip)(finalvert_t *pfv0, finalvert_t *pfv1, finalvert_t *out))
{
	s32 j = count - 1;
	s32 k = 0;
	for (s32 i = 0; i < count; j = i, i++) {
		s32 oldflags = in[j].flags & flag;
		s32 flags = in[i].flags & flag;
		if (flags && oldflags)
			continue;
		if (oldflags ^ flags) {
			clip(&in[j], &in[i], &out[k]);
			out[k].flags = 0;
			if (out[k].v[0] < r_refdef.aliasvrect.x)
				out[k].flags |= ALIAS_LEFT_CLIP;
			if (out[k].v[1] < r_refdef.aliasvrect.y)
				out[k].flags |= ALIAS_TOP_CLIP;
			if (out[k].v[0] > r_refdef.aliasvrectright)
				out[k].flags |= ALIAS_RIGHT_CLIP;
			if (out[k].v[1] > r_refdef.aliasvrectbottom)
				out[k].flags |= ALIAS_BOTTOM_CLIP;
			k++;
		}
		if (!flags) {
			out[k] = in[i];
			k++;
		}
	}
	return k;
}

void R_AliasClipTriangle(mtriangle_t *ptri)
{
	if (ptri->facesfront) { // copy vertexes and fix seam texture coords
		fv[0][0] = pfinalverts[ptri->vertindex[0]];
		fv[0][1] = pfinalverts[ptri->vertindex[1]];
		fv[0][2] = pfinalverts[ptri->vertindex[2]];
	} else {
		for (s32 i = 0; i < 3; i++) {
			fv[0][i] = pfinalverts[ptri->vertindex[i]];
			if (!ptri->facesfront
			    && (fv[0][i].flags & ALIAS_ONSEAM))
				fv[0][i].v[2] += r_affinetridesc.seamfixupX16;
		}
	}
	u32 clipflags = 
		fv[0][0].flags | fv[0][1].flags | fv[0][2].flags; // clip
	s32 k, pingpong;
	if (clipflags & ALIAS_Z_CLIP) {
		for (s32 i = 0; i < 3; i++)
			av[i] = pauxverts[ptri->vertindex[i]];
		k = R_AliasClip(fv[0], fv[1], ALIAS_Z_CLIP, 3, R_Alias_clip_z);
		pingpong = 1;
		clipflags = fv[1][0].flags | fv[1][1].flags | fv[1][2].flags;
	} else {
		pingpong = 0;
		k = 3;
	}
	if (clipflags & ALIAS_LEFT_CLIP) {
		k = R_AliasClip(fv[pingpong], fv[pingpong ^ 1],
				ALIAS_LEFT_CLIP, k, R_Alias_clip_left);
		if (k == 0) return;
		pingpong ^= 1;
	}
	if (clipflags & ALIAS_RIGHT_CLIP) {
		k = R_AliasClip(fv[pingpong], fv[pingpong ^ 1],
				ALIAS_RIGHT_CLIP, k, R_Alias_clip_right);
		if (k == 0) return;
		pingpong ^= 1;
	}
	if (clipflags & ALIAS_BOTTOM_CLIP) {
		k = R_AliasClip(fv[pingpong], fv[pingpong ^ 1],
				ALIAS_BOTTOM_CLIP, k, R_Alias_clip_bottom);
		if (k == 0) return;
		pingpong ^= 1;
	}
	if (clipflags & ALIAS_TOP_CLIP) {
		k = R_AliasClip(fv[pingpong], fv[pingpong ^ 1],
				ALIAS_TOP_CLIP, k, R_Alias_clip_top);
		if (k == 0) return;
		pingpong ^= 1;
	}
	for (s32 i = 0; i < k; i++) {
		if (fv[pingpong][i].v[0] < r_refdef.aliasvrect.x)
			fv[pingpong][i].v[0] = r_refdef.aliasvrect.x;
		else if (fv[pingpong][i].v[0] > r_refdef.aliasvrectright)
			fv[pingpong][i].v[0] = r_refdef.aliasvrectright;
		if (fv[pingpong][i].v[1] < r_refdef.aliasvrect.y)
			fv[pingpong][i].v[1] = r_refdef.aliasvrect.y;
		else if (fv[pingpong][i].v[1] > r_refdef.aliasvrectbottom)
			fv[pingpong][i].v[1] = r_refdef.aliasvrectbottom;
		fv[pingpong][i].flags = 0;
	}
	mtriangle_t mtri;
	mtri.facesfront = ptri->facesfront; // draw triangles
	r_affinetridesc.ptriangles = &mtri;
	r_affinetridesc.pfinalverts = fv[pingpong];
	mtri.vertindex[0] = 0; // FIXME: do all at once as trifan?
	for (s32 i = 1; i < k - 1; i++) {
		mtri.vertindex[1] = i;
		mtri.vertindex[2] = i + 1;
		D_PolysetDraw();
	}
}
