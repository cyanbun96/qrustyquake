// Copyright (c) 2024-2025 erysdren (it/its)
// GPLv3 See LICENSE for details.
// Originally MIT License

#include "quakedef.h"
#include "vgafont.h"

#define BLINK_HZ ((1000 / 70) * 16)

static const Uint8 palette[16][3] = { // ega 16-color palette
	{0x00,0x00,0x00}, {0x00,0x00,0xab}, {0x00,0xab,0x00}, {0x00,0xab,0xab},
	{0xab,0x00,0x00}, {0xab,0x00,0xab}, {0xab,0x57,0x00}, {0xab,0xab,0xab},
	{0x57,0x57,0x57}, {0x57,0x57,0xff}, {0x57,0xff,0x57}, {0x57,0xff,0xff},
	{0xff,0x57,0x57}, {0xff,0x57,0xff}, {0xff,0xff,0x57}, {0xff,0xff,0xff}
};

static void render_cell(Uint8 *image, int pitch, Uint16 cell, bool noblink)
{
	Uint8 code = (Uint8)(cell & 0xFF); // get components
	Uint8 attr = (Uint8)(cell >> 8);
	Uint8 blink = (attr >> 7) & 0x01; // break out attributes
	Uint8 bgcolor = (attr >> 4) & 0x07;
	Uint8 fgcolor = attr & 0x0F;
	for (int y = 0; y < 16; y++) {
		for (int x = 0; x < 8; x++) {
			image[y * pitch + x] = bgcolor; // write bgcolor
			if (blink && noblink) continue;
			Uint8 *bitmap = &VGA_FONT_CP437[code * 16]; // fgcolor
			if (bitmap[y] & 1 << SDL_abs(x - 7))
				image[y * pitch + x] = fgcolor;
		}
	}
}

int vgatext_main(SDL_Window *window, Uint16 *screen)
{
	Uint8 image0[400][640];
	Uint8 image1[400][640];
	SDL_PixelFormat format;
	SDL_Renderer *renderer;
	SDL_Texture *texture;
	SDL_Surface *surface8;
	SDL_Palette *surface8_palette;
	SDL_Surface **windowsurface;
	SDL_Surface *windowsurface0;
	SDL_Surface *windowsurface1;
	if (!window || !screen) return -1;
	renderer = SDL_GetRenderer(window); // setup renderer
	if (!renderer) return -1;
	SDL_SetRenderLogicalPresentation(renderer, 640, 400, SDL_LOGICAL_PRESENTATION_LETTERBOX);
	SDL_SetWindowMinimumSize(window, 640, 400);
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);
	SDL_RenderPresent(renderer);
	// setup render surface
	surface8 = SDL_CreateSurface(640, 400, SDL_PIXELFORMAT_INDEX8);
	if (!surface8) return -1;
	surface8_palette = SDL_CreateSurfacePalette(surface8); // setup palette
	for (int i = 0; i < 16; i++) {
		surface8_palette->colors[i].r = palette[i][0];
		surface8_palette->colors[i].g = palette[i][1];
		surface8_palette->colors[i].b = palette[i][2];
	}
	SDL_FillSurfaceRect(surface8, NULL, 0);
	// setup display surfaces
	format = SDL_GetWindowPixelFormat(window);
	windowsurface0 = SDL_CreateSurface(640, 400, format);
	windowsurface1 = SDL_CreateSurface(640, 400, format);
	if (!windowsurface0 || !windowsurface1) return -1;
	// setup render texture
	texture = SDL_CreateTexture(renderer, format, SDL_TEXTUREACCESS_STREAMING, 640, 400);
	if (!texture) return -1;
	SDL_SetTextureScaleMode(texture, SDL_SCALEMODE_NEAREST);
	SDL_Rect rect = {0, 0, 640, 400}; // setup blit rect
	for (int y = 0; y < 25; y++) { // render vgatext
		for (int x = 0; x < 80; x++) {
			uint16_t cell = screen[y * 80 + x];
			uint8_t *imgpos0 = &image0[y * 16][x * 8];
			uint8_t *imgpos1 = &image1[y * 16][x * 8];
			render_cell(imgpos0, 640, cell, false);
			render_cell(imgpos1, 640, cell, true);
		}
	}
	uint8_t *imgpos1 = &image1[24 * 16][0];
	render_cell(imgpos1, 640, 0x70dc, true); // CyanBun96: blinking cursor
	for (int y = 0; y < 400; y++) // create rgb image0
		SDL_memcpy(&((Uint8 *)surface8->pixels)[y * surface8->pitch], &image0[y][0], 640);
	SDL_BlitSurface(surface8, &rect, windowsurface0, &rect);
	for (int y = 0; y < 400; y++) // create rgb image1
		SDL_memcpy(&((Uint8 *)surface8->pixels)[y * surface8->pitch], &image1[y][0], 640);
	SDL_BlitSurface(surface8, &rect, windowsurface1, &rect);
	windowsurface = &windowsurface0; // save windowsurface ptr
	Uint64 next = SDL_GetTicks() + BLINK_HZ; // start counting time
	while (true) { // main loop
		SDL_Event event;
		while (SDL_PollEvent(&event)) {
			if (event.type == SDL_EVENT_QUIT || event.type == SDL_EVENT_KEY_DOWN || event.type == SDL_EVENT_MOUSE_BUTTON_DOWN || event.type == SDL_EVENT_JOYSTICK_BUTTON_DOWN)
				goto done;
		}
		Uint64 now = SDL_GetTicks();
		if (next <= now) {
			if (windowsurface == &windowsurface0)
				windowsurface = &windowsurface1;
			else
				windowsurface = &windowsurface0;
			next += BLINK_HZ;
		}
		SDL_UpdateTexture(texture, NULL, (*windowsurface)->pixels, (*windowsurface)->pitch);
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
