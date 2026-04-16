// Copyright(C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
// this is the only file outside the refresh that touches the vid buffer
#include "quakedef.h"

static s32 cliprectx0 = -1;
static s32 cliprectx1 = -1;
static s32 cliprecty0 = -1;
static s32 cliprecty1 = -1;
static rectdesc_t r_rectdesc;
static cachepic_t menu_cachepics[MAX_CACHED_PICS];
static s32 menu_numcachepics;
static u8 *draw_chars; // 8*8 graphic characters
static qpic_t *draw_backtile;
static f32 basemip[NUM_MIPS - 1] = { 1.0, 0.5 * 0.8, 0.25 * 0.8 };
static u8 dos_brightness_lut[256*3];
static s32 dos_brightness_lut_init = 0;
static u8 brightness_lut[256*3];
static s32 brightness_lut_init = 0;

s32 D_Dither(u8 *pos, f32 opacity)
{
	if(opacity >= 1.0f) return 1;
	if(opacity <= 0.01f) return 0;
	s32 dither_pat = opacity * 7.2f;
	u64 d = pos - vid.buffer;
	u64 x = d % vid.width;
	u64 y = d / vid.width;
	switch(dither_pat){
		case 0: return !(d % 6); // 1/6
		case 1: return (y&1) && ((y&3) == 3 ? (x&1) : !(x&1)); // 1/4
		case 2: return !(d % 3); // 1/3
		case 3: return (x + y) & 1; // 1/2
		case 4: return d % 3; // 2/3
		case 5: return !((y&1) && ((y&3) == 3 ? (x&1) : !(x&1))); // 3/4
		default: case 6: return d % 6; // 5/6
	}
}

void D_Init()
{
	Cvar_RegisterVariable(&d_mipcap);
	Cvar_RegisterVariable(&d_mipscale);
}

void D_SetupFrame()
{
	d_viewbuffer = r_dowarp ? r_warpbuffer : (void *)vid.buffer;
	screenwidth = vid.width;
	d_roverwrapped = 0;
	d_initial_rover = sc_rover;
	d_minmip = d_mipcap.value;
	d_minmip = CLAMP(0, d_minmip, 3);
	for(s32 i = 0; i < (NUM_MIPS - 1); i++)
		d_scalemip[i] = basemip[i] * d_mipscale.value;
}

void D_ViewChanged()
{
	s32 rowbytes = vid.width;
	scale_for_mip = xscale;
	if(yscale > xscale)
		scale_for_mip = yscale;
	d_zwidth = vid.width;
	d_pix_min = r_refdef.vrect.width / 320;
	if(d_pix_min < 1)
		d_pix_min = 1;
	d_pix_max = (s32)((f32)r_refdef.vrect.width / (320.0 / 4.0) + 0.5);
	d_pix_shift = 8 - (s32)((f32)r_refdef.vrect.width / 320.0 + 0.5);
	if(d_pix_max < 1)
		d_pix_max = 1;
	if(pixelAspect > 1.4)
		d_y_aspect_shift = 1;
	else
		d_y_aspect_shift = 0;
	d_vrectx = r_refdef.vrect.x;
	d_vrecty = r_refdef.vrect.y;
	for(u32 i = 0; i < vid.height; i++){
		d_scantable[i] = i * rowbytes;
		zspantable[i] = d_pzbuffer + i * d_zwidth;
	}
}

qpic_t *Draw_PicFromWad(s8 *name) { return W_GetLumpName(name); }

qpic_t *Draw_CachePic(s8 *path)
{
	cachepic_t *pic = menu_cachepics;
	s32 i = 0;
	for(; i < menu_numcachepics; pic++, i++)
		if(!strcmp(path, pic->name))
			break;
	if(i == menu_numcachepics){
		if(menu_numcachepics == MAX_CACHED_PICS)
			Sys_Error("menu_numcachepics == MAX_CACHED_PICS");
		menu_numcachepics++;
		strcpy(pic->name, path);
	}
	qpic_t *dat = Cache_Check(&pic->cache);
	if(dat)
		return dat;
	COM_LoadCacheFile(path, &pic->cache, NULL); // load the pic from disk
	dat = (qpic_t *) pic->cache.data;
	if(!dat)
		Sys_Error("Draw_CachePic: failed to load %s", path);
	SwapPic(dat);
	return dat;
}

qpic_t *Draw_TryCachePic(s8 *path)
{
	qpic_t *dat = (qpic_t*)COM_LoadMallocFile(path, NULL);
	if(!dat){
		Con_DPrintf("Draw_TryCachePic: failed to load %s\n", path);
		return NULL;
	}
	SwapPic(dat);
	return dat;
}

void Draw_Init()
{
	draw_chars = W_GetLumpName("conchars");
	draw_disc = W_GetLumpName("disc");
	draw_backtile = W_GetLumpName("backtile");
	r_rectdesc.width = draw_backtile->width;
	r_rectdesc.height = draw_backtile->height;
	r_rectdesc.ptexbytes = draw_backtile->data;
	r_rectdesc.rowbytes = draw_backtile->width;
}

void Draw_Character_Ex(f32 *pos, f32 *sz, s32 num, f32 *color, f32 alpha)
{ // CSQC version with alpha and RGB colors
	if(alpha <= 0) return;
	s32 al = (1-alpha) * FOG_LUT_LEVELS;
	if(al > FOG_LUT_LEVELS - 1) al = FOG_LUT_LEVELS - 1;
	if(al < 0) return;
	if(!fog_lut_built) R_BuildColorMixLUT(0);
	bool defcolor = 0;
	s32 c;
	if(color[0] == 1.0 && color[1] == 1.0 && color[2] == 1.0) defcolor = 1;
	else {
		u8(*convfunc)(u8,u8,u8)=r_labmixpal.value==1?rgbtoi_lab:rgbtoi;
		c = convfunc(color[0]*255, color[1]*255, color[2]*255);
	}
	s32 draw_w = (s32)sz[0];
	s32 draw_h = (s32)sz[1];
	num &= 255;
	s32 row = num >> 4;
	s32 col = num & 15;
	u8 *char_base = draw_chars + (row << 10) + (col << 3);
	s32 base_x = (s32)pos[0];
	s32 base_y = (s32)pos[1];
	s32 clipx0 = cliprectx0 == -1 ? 0 : cliprectx0;
	s32 clipx1 = cliprectx1 == -1 ? (s32)vid.width : cliprectx1;
	s32 clipy0 = cliprecty0 == -1 ? 0 : cliprecty0;
	s32 clipy1 = cliprecty1 == -1 ? (s32)vid.height : cliprecty1;
	s32 startx = base_x < clipx0 ? clipx0 : base_x;
        s32 starty = base_y < clipy0 ? clipy0 : base_y;
        s32 endx = (base_x + draw_w) > clipx1 ? clipx1 : (base_x + draw_w);
        s32 endy = (base_y + draw_h) > clipy1 ? clipy1 : (base_y + draw_h);
	if(startx >= endx || starty >= endy) return;
	u8 *fb = (u8*)scrbuffs[drawlayer]->pixels;
	u8 *fb0 = (u8*)scrbuffs[0]->pixels;
	for(s32 dy = starty; dy < endy; dy++){
		s32 desty = dy - base_y;
		s32 src_y = (s32)(8.0f * ((f32)desty / draw_h));
		if(src_y < 0) src_y = 0;
		if(src_y > 7) src_y = 7;
		u8 *dest = fb + dy * vid.width + startx;
		for(s32 dx = startx; dx < endx; dx++){
			s32 destx = dx - base_x;
			s32 src_x = (s32)(8.0f * ((f32)destx / draw_w));
			if(src_x < 0) src_x = 0;
			if(src_x > 7) src_x = 7;
			u8 *source = char_base + src_y * 128 + src_x;
			if(*source){
				s32 charcolor = defcolor ? *source :
				    color_mix_lut[*source][c][FOG_LUT_LEVELS/2];
				s32 c2;
				if(*dest != TRANSPARENT_COLOR)
					c2 = *dest;
				else
					c2 = fb0[dest - fb];
				*dest = color_mix_lut[charcolor][c2][al];
			}
			dest++;
		}
	}
}

void Draw_CharacterScaled(s32 x, s32 y, s32 num, s32 scale)
{ // Draws one 8*8 graphics character with 0 being transparent.
  // It can be clipped to the top of the screen to allow the console
  // to be smoothly scrolled off.
	num &= 255;
	if(y <= -8 * scale)
		return; // totally off screen
	if(y > (s32)(vid.height - 8 * scale))
		return; // don't draw past the bottom of the screen
	s32 row = num >> 4;
	s32 col = num & 15;
	u8 *source = draw_chars + (row << 10) + (col << 3);
	s32 drawline;
	s32 row_remainder = 0;
	if(y < 0){
		s32 clipped_pixels = -y;
		s32 clipped_rows = clipped_pixels / scale;
		row_remainder = clipped_pixels % scale;
		if(clipped_rows >= 8)
			return;
		source += 128 * clipped_rows;
		drawline = 8 - clipped_rows;
		y = 0;
	} else
		drawline = 8;
	u8 *dest = (u8*)scrbuffs[drawlayer]->pixels + y * vid.width + x;
	dest -= row_remainder * vid.width; // avoid jitter
	while(drawline--){
		for(s32 k = 0; k < scale; ++k){
			if(dest >= (u8*)scrbuffs[drawlayer]->pixels)
				for(s32 j = 0; j < scale; ++j)
				for(s32 i = 0; i < 8; ++i)
					if(source[i])
						dest[i * scale + j] = source[i];
			dest += vid.width;
		}
		source += 128;
	}
}

void Draw_StringScaled(s32 x, s32 y, s8 *str, s32 scale)
{
	while(*str){
		Draw_CharacterScaled(x, y, *str, scale);
		str++;
		x += 8 * scale;
	}
}

void Draw_PicScaled(s32 x, s32 y, qpic_t *pic, s32 scale)
{
	u8 *dest = (u8*)scrbuffs[drawlayer]->pixels + y * vid.width + x;
	s32 maxwidth = pic->width;
	if(x + pic->width * scale > (s32)vid.width)
		maxwidth -= (x + pic->width * scale - vid.width) / scale;
	s32 rowinc = vid.width;
	if(y + pic->height * scale > (s32)vid.height)
		return;
	for(u32 v = 0; v < (u32)pic->height; v++){
		if(v * scale + y >= vid.height)
			break;
		u8 *source = pic->data + v * pic->width;
		for(s32 k = 0; k < scale; k++){
			u8 *dest_row = dest + k * rowinc;
			for(u32 i = 0; i < (u32)maxwidth; i++){
				u8 pixel = source[i];
				for(s32 j = 0; j < scale; j++)
					dest_row[(i*scale) + j] = pixel;
			}
		}
		dest += rowinc * scale;
	}
}

void Draw_PicScaledPartial(s32 x,s32 y,s32 l,s32 t,s32 w,s32 h,qpic_t *p,s32 s)
{
        u8 *source = p->data;
		u8 *dest = (u8*)scrbuffs[drawlayer]->pixels + y * vid.width + x;
        for(u32 v = 0; v < (u32)p->height; v++){
                if(v * s + y >= vid.height || v > (u32)h)
                        return;
                if(v < (u32)t){
			source += p->width;
                        continue;
		}
                for(s32 k = 0; k < s; k++){
                        for(u32 i = 0; i < (u32)p->width; i++){
                                if(i < (u32)l || i >= (u32)w)
                                        continue;
                                for(s32 j = 0; j < s; j++)
                                        if(i * s + j + x < vid.width)
                                                dest[i * s + j] = source[i];
                        }
                        dest += vid.width;
                }
                source += p->width;
        }	
}

void Draw_TransPicScaled(s32 x, s32 y, qpic_t *pic, s32 scale)
{
	u8 *source = pic->data;
	u8 *dest = (u8*)scrbuffs[drawlayer]->pixels + y * vid.width + x;
	for(u32 v = 0; v < (u32)pic->height; v++){
		if(v * scale + y >= vid.height)
			return;
		for(s32 k = 0; k < scale; k++){
			for(u32 u = 0; u < (u32)pic->width; u++){
				u8 tbyte = source[u];
				if(tbyte == TRANSPARENT_COLOR)
					continue;
				for(s32 i = 0; i < scale; i++)
					if(u * scale + i + x < vid.width)
						dest[u * scale + i] = tbyte;
			}
			dest += vid.width;
		}
		source += pic->width;
	}
}

void Draw_TransPicTranslateScaled(s32 x, s32 y, qpic_t *p, u8 *tl, s32 scale)
{
	u8 *source = p->data;
	u8 *dest = (u8*)scrbuffs[drawlayer]->pixels + y * vid.width + x;
	for(s32 v = 0; v < p->height; v++){
		if(v * scale + y >= (s32)vid.height)
			return;
		for(s32 k = 0; k < scale; k++){
			for(s32 u = 0; u < p->width; u++){
				u8 tbyte = source[u];
				if(tbyte == TRANSPARENT_COLOR)
					continue;
				for(s32 i = 0; i < scale; i++)
					if(u * scale + i + x < (s32)vid.width)
						dest[u * scale + i] = tl[tbyte];
			}
			dest += vid.width;
		}
		source += p->width;
	}
}

void Draw_Pic_Ex(f32 *pos, f32 *sz, qpic_t *pic, f32 *srcpos, f32 *srcsz,
		f32 *color, f32 alpha)
{ // CSQC version with scaling, alpha, and RGB colors
	if(alpha <= 0) return;
	s32 al = (1-alpha) * FOG_LUT_LEVELS;
	if(al > FOG_LUT_LEVELS - 1) al = FOG_LUT_LEVELS - 1;
	if(al < 0) return;
	if(!fog_lut_built) R_BuildColorMixLUT(0);
	bool defcolor = 0;
	s32 c;
	if(color[0] == 1.0 && color[1] == 1.0 && color[2] == 1.0) defcolor = 1;
	else {
		u8(*convfunc)(u8,u8,u8)=r_labmixpal.value==1?rgbtoi_lab:rgbtoi;
		c = convfunc(color[0]*255, color[1]*255, color[2]*255);
	}
	s32 draw_w = (int)sz[0];
	s32 draw_h = (int)sz[1];
	s32 pic_w = pic->width;
	s32 pic_h = pic->height;
	s32 base_x = (s32)pos[0];
	s32 base_y = (s32)pos[1];
	s32 clipx0 = cliprectx0 == -1 ? 0 : cliprectx0;
	s32 clipx1 = cliprectx1 == -1 ? (s32)vid.width : cliprectx1;
	s32 clipy0 = cliprecty0 == -1 ? 0 : cliprecty0;
	s32 clipy1 = cliprecty1 == -1 ? (s32)vid.height : cliprecty1;
	s32 startx = base_x < clipx0 ? clipx0 : base_x;
	s32 starty = base_y < clipy0 ? clipy0 : base_y;
	s32 endx = (base_x + draw_w) > clipx1 ? clipx1 : (base_x + draw_w);
	s32 endy = (base_y + draw_h) > clipy1 ? clipy1 : (base_y + draw_h);
	if(startx >= endx || starty >= endy) return;
	u8 *fb = (u8*)scrbuffs[drawlayer]->pixels;
	u8 *fb0 = (u8*)scrbuffs[0]->pixels;
	for(s32 dy = starty; dy < endy; dy++){
		s32 desty = dy - base_y;
		s32 y = srcpos[1] * pic_h +
			(srcsz[1] * pic_h) * ((f32)desty / draw_h);
		if(y < 0) y = 0;
		if(y >= pic_h) y = pic_h - 1;
		u8 *dest = fb + dy * vid.width + startx;
		for(s32 dx = startx; dx < endx; dx++){
			s32 destx = dx - base_x;
			s32 x = srcpos[0] * pic_w +
				(srcsz[0] * pic_w) * ((f32)destx / draw_w);
			if(x < 0) x = 0;
			if(x >= pic_w) x = pic_w - 1;
			u8 *source = pic->data + y * pic_w + x;
			if(*source != TRANSPARENT_COLOR){
				s32 charcolor = defcolor ? *source :
				    color_mix_lut[*source][c][FOG_LUT_LEVELS/2];
				s32 c2;
				if(*dest != TRANSPARENT_COLOR)
					c2 = *dest;
				else
					c2 = fb0[dest - fb];
				*dest = color_mix_lut[charcolor][c2][al];
			}
			dest++;
		}
	}
}

void Draw_CharToConback(s32 num, u8 *dest)
{
	s32 row = num >> 4;
	s32 col = num & 15;
	u8 *source = draw_chars + (row << 10) + (col << 3);
	s32 drawline = 8;
	while(drawline--){
		for(s32 x = 0; x < 8; x++)
			if(source[x])
				dest[x] = 0x60 + source[x];
		source += 128;
		dest += 320;
	}
}

void Draw_CharToConbackScaled(s32 num, u8 *dest, s32 scale, s32 width)
{
	s32 row = num >> 4;
	s32 col = num & 15;
	u8 *source = draw_chars + (row << 10) + (col << 3);
	s32 drawline = 8;
	while(drawline--){
		for(s32 x = 0; x < 8*scale; x++)
			if(source[x/scale])
				for(s32 s = 0; s < scale; s++)
					dest[x+s*width] = 0x60 + source[x/scale];
		source += 128;
		dest += width*scale;
	}
}

void Draw_ConsoleBackground(s32 lines)
{
	s8 ver[100];
	qpic_t *conback = Draw_CachePic("gfx/conback.lmp");
	// hack the version number directly into the pic
	sprintf(ver, "(QrustyQuake) %4.2f", (f32)VERSION);
	s32 scale = conback->width / 320;
	u8 *dest = conback->data + conback->width*(conback->height-14*scale)
		+ conback->width - 11*scale - 8*scale * strlen(ver);
	for(u64 x = 0; x < strlen(ver); x++)
		Draw_CharToConbackScaled(ver[x], dest + x * 8 * scale, scale, conback->width);
	dest = vid.buffer; // draw the pic
	for(s32 y = 0; y < lines; y++, dest += vid.width){
		s32 v = (vid.height-lines+y)*conback->height/vid.height;
		u8 *src = conback->data + v * conback->width;
		if((s32)vid.width == conback->width)
			memcpy(dest, src, vid.width);
		else {
			s32 f = 0;
			s32 fstep = conback->width * 0x10000 / vid.width;
			for(u32 x = 0; x < vid.width; x += 4){
				dest[x] = src[f >> 16];
				f += fstep;
				dest[x + 1] = src[f >> 16];
				f += fstep;
				dest[x + 2] = src[f >> 16];
				f += fstep;
				dest[x + 3] = src[f >> 16];
				f += fstep;
			}
		}
	}
}

void R_DrawRect(vrect_t *prect, s32 rowbytes, u8 *psrc, s32 transparent)
{
	u8 *pdest = (u8*)scrbuffs[drawlayer]->pixels + (prect->y * vid.width) + prect->x;
	u64 maxdest = (u64)(unsigned long long)scrbuffs[drawlayer]->pixels+vid.width*vid.height;
	if(transparent){
		s32 srcdelta = rowbytes - prect->width;
		s32 destdelta = vid.width - prect->width;
		for(s32 i = 0; i < prect->height; i++){
			for(s32 j = 0; j < prect->width; j++){
				u8 t = *psrc;
				if(t != TRANSPARENT_COLOR)
					*pdest = t;
				psrc++;
				pdest++;
			}
			psrc += srcdelta;
			pdest += destdelta;
		}
	} else {
		for(s32 i = 0; i < prect->height; i++){
			if((u64)(unsigned long long)pdest+prect->width >= maxdest) break;
			memcpy(pdest, psrc, prect->width);
			psrc += rowbytes;
			pdest += vid.width;
		}
	}
}

void Draw_TileClear(s32 x, s32 y, s32 w, s32 h) // This repeats a 64*64
{ // tile graphic to fill the screen around a sized down refresh window
	r_rectdesc.rect.x = x;
	r_rectdesc.rect.y = y;
	r_rectdesc.rect.width = w;
	r_rectdesc.rect.height = h;
	vrect_t vr;
	vr.y = r_rectdesc.rect.y;
	s32 height = r_rectdesc.rect.height;
	s32 tileoffsety = vr.y % r_rectdesc.height;
	while(height > 0){
		vr.x = r_rectdesc.rect.x;
		s32 width = r_rectdesc.rect.width;
		vr.height = r_rectdesc.height - tileoffsety;
		if(vr.height > height)
			vr.height = height;
		s32 tileoffsetx = vr.x % r_rectdesc.width;
		while(width > 0){
			vr.width = r_rectdesc.width - tileoffsetx;
			if(vr.width > width)
				vr.width = width;
			u8 *psrc = r_rectdesc.ptexbytes +
			    (tileoffsety * r_rectdesc.rowbytes) + tileoffsetx;
			R_DrawRect(&vr, r_rectdesc.rowbytes, psrc, 0);
			vr.x += vr.width;
			width -= vr.width;
			tileoffsetx = 0; // only left tile can be left-clipped
		}
		vr.y += vr.height;
		height -= vr.height;
		tileoffsety = 0; // only top tile can be top-clipped
	}
}

void Draw_Fill(s32 x, s32 y, s32 w, s32 h, s32 c)
{ // Fills a box of pixels with a single color
	u8 *dest = (u8*)scrbuffs[drawlayer]->pixels + y * vid.width + x;
	for(s32 v = 0; v < h; v++, dest += vid.width)
		memset(dest, c, w); // Fast horizontal fill
}

void Draw_FillEx(s32 x, s32 y, s32 w, s32 h, f32 *rgb, f32 alpha)
{ // CSQC version with alpha and RGB colors
	if(alpha <= 0) return;
	u8(*convfunc)(u8,u8,u8) = r_labmixpal.value == 1 ? rgbtoi_lab : rgbtoi;
	u8 c = convfunc(rgb[0]*255, rgb[1]*255, rgb[2]*255);
	s32 al = (1-alpha) * FOG_LUT_LEVELS;
	if(al > FOG_LUT_LEVELS - 1) al = FOG_LUT_LEVELS - 1;
	if(al < 0) return;
	s32 clipx0 = cliprectx0 == -1 ? 0 : cliprectx0;
	s32 clipx1 = cliprectx1 == -1 ? (s32)vid.width : cliprectx1;
	s32 clipy0 = cliprecty0 == -1 ? 0 : cliprecty0;
	s32 clipy1 = cliprecty1 == -1 ? (s32)vid.height : cliprecty1;
	s32 startx = x < clipx0 ? clipx0 : x;
        s32 starty = y < clipy0 ? clipy0 : y;
        s32 endx = (x + w) > clipx1 ? clipx1 : (x + w);
        s32 endy = (y + h) > clipy1 ? clipy1 : (y + h);
	if(startx >= endx || starty >= endy) return;
	if(al == FOG_LUT_LEVELS - 1){
		Draw_Fill(startx, starty, endx - startx, endy - starty, c);
		return;
	}
	if(!fog_lut_built) R_BuildColorMixLUT(0);
	u8 *fb  = (u8*)scrbuffs[drawlayer]->pixels;
	u8 *fb0 = (u8*)scrbuffs[0]->pixels;
	for(s32 dy = starty; dy < endy; dy++){
		u8 *dest = fb + dy * vid.width + startx;
		for(s32 dx = startx; dx < endx; dx++){
			s32 c2;
			if(*dest != TRANSPARENT_COLOR)
				c2 = *dest;
			else
				c2 = fb0[dest - fb];
			*dest = color_mix_lut[c][c2][al];
			dest++;
		}
	}
}

void Draw_SetClipRect(s32 x0, s32 x1, s32 y0, s32 y1)
{ cliprectx0 = x0; cliprectx1 = x1; cliprecty0 = y0; cliprecty1 = y1; }
void Draw_ResetClipping()
{ cliprectx0 = -1; cliprectx1 = -1; cliprecty0 = -1; cliprecty1 = -1; }

void Draw_InitBrightnessDOSLUT()
{ // Used only for the DOS screen fade effect, thus the range [0x10-0x1F]
	dos_brightness_lut_init++;
	for(s32 i = 0; i < 256; ++i){
		dos_brightness_lut[i] = ((host_basepal[i*3] + host_basepal[i*3+1] + 
			host_basepal[i*3+2]) / 3) / 16 + 0x10;
	}
}

void Draw_InitBrightnessLUT()
{ // Mix 50-50 with black
	brightness_lut_init++;
	for(s32 i = 0; i < 256; ++i){
		brightness_lut[i] = rgbtoi(host_basepal[i*3] / 2, 
					host_basepal[i*3+1] / 2,
					host_basepal[i*3+2] / 2);
	}
}

void Draw_FadeScreen()
{
	fadescreen = 0;
	if(scr_menubgstyle.value == 3) return;
	if(scr_menubgstyle.value == 1){ // DOS
		s32 area = vid.width * vid.height;
		u8 *pdest = (u8*)scrbuffs[0]->pixels;
		if(!dos_brightness_lut_init) Draw_InitBrightnessDOSLUT();
		for(s32 i = 0; i < area; ++i)
			pdest[i] = dos_brightness_lut[pdest[i]];
		return;
	}
	if(scr_menubgstyle.value == 2){ // darken
		s32 area = vid.width * vid.height;
		u8 *pdest = (u8*)scrbuffs[0]->pixels;
		if(!brightness_lut_init) Draw_InitBrightnessLUT();
		for(s32 i = 0; i < area; ++i)
			pdest[i] = brightness_lut[pdest[i]];
		return;
	}
	for(u32 y = 0; y < vid.height / uiscale; y++) // winquake
		for(u32 i = 0; i < uiscale; i++){
			u8 *pbuf = (u8*)scrbuffs[0]->pixels + vid.width * y
				* uiscale + vid.width * i;
			u32 t = (y & 1) << 1;
			for(u32 x = 0; x < vid.width / uiscale; x++)
				if((x & 3) != t)
					for(u32 j = 0; j < uiscale; j++)
						pbuf[x * uiscale + j] = 0;
		}
}
