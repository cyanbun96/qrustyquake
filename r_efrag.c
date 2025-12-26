// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

static efrag_t **lastlink;
static entity_t *r_addent;

void R_RemoveEfrags(entity_t *ent)
{// Call when removing an object from the world or moving it to another position
	efrag_t *ef = ent->efrag;
	while (ef) {
		efrag_t **prev = &ef->leaf->efrags;
		while (1) {
			efrag_t *walk = *prev;
			if (!walk)
				break;
			if (walk == ef) { // remove this fragment
				*prev = ef->leafnext;
				break;
			} else
				prev = &walk->leafnext;
		}
		efrag_t *old = ef;
		ef = ef->entnext;
		old->entnext = cl.free_efrags; // put it on the free list
		cl.free_efrags = old;
	}
	ent->efrag = NULL;
}

static efrag_t *R_GetEfrag () // based on RMQEngine
{
        // we could just Hunk_Alloc a single efrag_t and return it, but since
        // the struct is so small (2 pointers) allocate groups of them
        // to avoid wasting too much space on the hunk allocation headers.
        if (cl.free_efrags) {
                efrag_t *ef = cl.free_efrags;
                cl.free_efrags = ef->leafnext;
                ef->leafnext = NULL;
                cl.num_efrags++;
                return ef;
        }
        else {
                cl.free_efrags = (efrag_t *) Hunk_AllocName (EXTRA_EFRAGS * sizeof (efrag_t), "efrags");
		s32 i = 0;
                for (; i < EXTRA_EFRAGS - 1; i++)
                        cl.free_efrags[i].leafnext = &cl.free_efrags[i + 1];
                cl.free_efrags[i].leafnext = NULL;
                // call recursively to get a newly allocated free efrag
                return R_GetEfrag ();
        }
}

void R_SplitEntityOnNode(mnode_t *node)
{
	if (node->contents == CONTENTS_SOLID)
		return;
	if (node->contents < 0) { // add an efrag if the node is a leaf
		if (!r_pefragtopnode)
			r_pefragtopnode = node;
		mleaf_t *leaf = (mleaf_t *)node;
		efrag_t *ef = R_GetEfrag(); // grab an efrag off the free list
		ef->entity = r_addent;
		ef->leafnext = leaf->efrags; // set the leaf links
		leaf->efrags = ef;
		return;
	}
	mplane_t *splitplane = node->plane; // NODE_MIXED
	s32 sides = BoxOnPlaneSide(r_emins, r_emaxs, splitplane);
	if (sides == 3) {
		// split on this plane
		// if this is the first splitter of this bmodel, remember it
		if (!r_pefragtopnode)
			r_pefragtopnode = node;
	}
	if (sides & 1) // recurse down the contacted sides
		R_SplitEntityOnNode(node->children[0]);
	if (sides & 2)
		R_SplitEntityOnNode(node->children[1]);
}

void R_SplitEntityOnNode2(mnode_t *node)
{
	if (node->visframe != r_visframecount)
		return;
	if (node->contents < 0) {
		if (node->contents != CONTENTS_SOLID)
			r_pefragtopnode = node;	// we've reached a non-solid
		return; // leaf, so it's  visible and not BSP clipped
	}
	mplane_t *splitplane = node->plane;
	s32 sides = BoxOnPlaneSide(r_emins, r_emaxs, splitplane);
	if (sides == 3) {
		r_pefragtopnode = node; // remember first splitter
		return;
	}
	if (sides & 1) // not split yet; recurse down the contacted side
		R_SplitEntityOnNode2(node->children[0]);
	else
		R_SplitEntityOnNode2(node->children[1]);
}

void R_AddEfrags(entity_t *ent)
{
	if (!ent->model || ent == cl_entities)
		return; // never add the world
	r_addent = ent;
	lastlink = &ent->efrag;
	r_pefragtopnode = NULL;
	model_t *entmodel = ent->model;
	for (s32 i = 0; i < 3; i++) {
		r_emins[i] = ent->origin[i] + entmodel->mins[i];
		r_emaxs[i] = ent->origin[i] + entmodel->maxs[i];
	}
	R_SplitEntityOnNode(cl.worldmodel->nodes);
	ent->topnode = r_pefragtopnode;
}

void R_StoreEfrags(efrag_t **ppefrag)
{ // FIXME: a lot of this goes away with edge-based
	efrag_t *pefrag;
	while ((pefrag = *ppefrag) != NULL) {
		entity_t *pent = pefrag->entity;
		model_t *clmodel = pent->model;
		switch (clmodel->type) {
		case mod_alias:
		case mod_brush:
		case mod_sprite:
			pent = pefrag->entity;
			if (((u32)pent->visframe != r_framecount) &&
			    (cl_numvisedicts < MAX_VISEDICTS)) {
				cl_visedicts[cl_numvisedicts++] = pent;
			// mark that we've recorded this entity for this frame
				pent->visframe = r_framecount;
			}
			ppefrag = &pefrag->leafnext;
			break;
		default:
			Sys_Error("R_StoreEfrags: Bad entity type %d\n",
				  clmodel->type);
		}
	}
}
