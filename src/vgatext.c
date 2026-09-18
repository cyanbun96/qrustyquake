// Copyright (c) 2024-2025 erysdren (it/its)
// GPLv3 See LICENSE for details.
// Originally MIT License
#include "quakedef.h"
#include "vgafont.h"

#define BLINK_HZ ((1000 / 70) * 16)

static const u8 palette[16][3] = { // ega 16-color palette
	{0x00,0x00,0x00}, {0x00,0x00,0xab}, {0x00,0xab,0x00}, {0x00,0xab,0xab},
	{0xab,0x00,0x00}, {0xab,0x00,0xab}, {0xab,0x57,0x00}, {0xab,0xab,0xab},
	{0x57,0x57,0x57}, {0x57,0x57,0xff}, {0x57,0xff,0x57}, {0x57,0xff,0xff},
	{0xff,0x57,0x57}, {0xff,0x57,0xff}, {0xff,0xff,0x57}, {0xff,0xff,0xff}
};

static void render_cell(u8 *image, s32 pitch, Uint16 cell, bool noblink)
{
	u8 code = (u8)(cell & 0xFF); // get components
	u8 attr = (u8)(cell >> 8);
	u8 blink = (attr >> 7) & 0x01; // break out attributes
	u8 bgcolor = (attr >> 4) & 0x07;
	u8 fgcolor = attr & 0x0F;
	for (s32 y = 0; y < 16; y++) {
		for (s32 x = 0; x < 8; x++) {
			image[y * pitch + x] = bgcolor; // write bgcolor
			if (blink && noblink) continue;
			u8 *bitmap = &VGA_FONT_CP437[code * 16]; // fgcolor
			if (bitmap[y] & 1 << SDL_abs(x - 7))
				image[y * pitch + x] = fgcolor;
		}
	}
}

s32 vgatext_main(SDL_Window *window, u16 *screen)
{
	u8 image0[400][640];
	u8 image1[400][640];
	if(!window || !screen)return -1;
	SDL_Renderer *renderer = SDL_GetRenderer(window); // setup renderer
	if(!renderer)return -1;
	SDL_SetRenderLogicalPresentation(
			renderer, 640, 400, SDL_LOGICAL_PRESENTATION_LETTERBOX);
	SDL_SetWindowMinimumSize(window, 640, 400);
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);
	SDL_RenderPresent(renderer);
	// setup render surface
	SDL_Surface *surface8=SDL_CreateSurface(640,400,SDL_PIXELFORMAT_INDEX8);
	if(!surface8)return -1;
	SDL_Palette *surface8_palette = SDL_CreateSurfacePalette(surface8);
	for(s32 i = 0; i < 16; i++){ // setup palette
		surface8_palette->colors[i].r = palette[i][0];
		surface8_palette->colors[i].g = palette[i][1];
		surface8_palette->colors[i].b = palette[i][2];
	}
	SDL_FillSurfaceRect(surface8, NULL, 0);
	// setup display surfaces
	SDL_PixelFormat format = SDL_GetWindowPixelFormat(window);
	SDL_Surface *windowsurface0 = SDL_CreateSurface(640, 400, format);
	SDL_Surface *windowsurface1 = SDL_CreateSurface(640, 400, format);
	if (!windowsurface0 || !windowsurface1) return -1;
	// setup render texture
	SDL_Texture *texture = SDL_CreateTexture(
			renderer, format, SDL_TEXTUREACCESS_STREAMING,640,400);
	if (!texture) return -1;
	SDL_SetTextureScaleMode(texture, SDL_SCALEMODE_NEAREST);
	SDL_Rect rect = {0, 0, 640, 400}; // setup blit rect
	for (s32 y = 0; y < 25; y++) { // render vgatext
		for (s32 x = 0; x < 80; x++) {
			s16 cell = screen[y * 80 + x];
			u8 *imgpos0 = &image0[y * 16][x * 8];
			u8 *imgpos1 = &image1[y * 16][x * 8];
			render_cell(imgpos0, 640, cell, false);
			render_cell(imgpos1, 640, cell, true);
		}
	}
	u8 *imgpos1 = &image1[24 * 16][0];
	render_cell(imgpos1, 640, 0x70dc, true); // CyanBun96: blinking cursor
	for (s32 y = 0; y < 400; y++) // create rgb image0
		SDL_memcpy(&((u8 *)surface8->pixels)[y * surface8->pitch],
				&image0[y][0], 640);
	SDL_BlitSurface(surface8, &rect, windowsurface0, &rect);
	for (s32 y = 0; y < 400; y++) // create rgb image1
		SDL_memcpy(&((u8 *)surface8->pixels)[y * surface8->pitch],
				&image1[y][0], 640);
	SDL_BlitSurface(surface8, &rect, windowsurface1, &rect);
	SDL_Surface **windowsurface = &windowsurface0; // save windowsurface ptr
	u64 next = SDL_GetTicks() + BLINK_HZ; // start counting time
	while(true){ // main loop
		SDL_Event event;
		while(SDL_PollEvent(&event)){
			if (event.type == SDL_EVENT_QUIT
				|| event.type == SDL_EVENT_KEY_DOWN
				|| event.type == SDL_EVENT_MOUSE_BUTTON_DOWN
				|| event.type == SDL_EVENT_JOYSTICK_BUTTON_DOWN)
				goto done;
		}
		u64 now = SDL_GetTicks();
		if(next <= now){
			if (windowsurface == &windowsurface0)
				windowsurface = &windowsurface1;
			else
				windowsurface = &windowsurface0;
			next += BLINK_HZ;
		}
		SDL_Delay(1000/30); // CyanBun96: FPS limit
		SDL_UpdateTexture(texture, NULL,
			(*windowsurface)->pixels, (*windowsurface)->pitch);
		SDL_RenderClear(renderer);
		SDL_RenderTexture(renderer, texture, NULL, NULL);
		SDL_RenderPresent(renderer);
	}
done:
	SDL_DestroySurface(surface8);
	SDL_DestroySurface(windowsurface0);
	SDL_DestroySurface(windowsurface1);
	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
	return 0;
}
