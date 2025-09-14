// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.

#include "quakedef.h"

static s32 vid_line;
static s32 drawmousemenu = 0;
static s32 m_save_demonum;
static s32 m_main_cursor;
static s32 m_singleplayer_cursor;
static s32 load_cursor; // 0 < load_cursor < MAX_SAVEGAMES
static s8 m_filenames[MAX_SAVEGAMES][MAX_OSPATH*2 + 1];
static s32 loadable[MAX_SAVEGAMES];
static s32 m_multiplayer_cursor;
static s32 setup_cursor = 4;
static s32 setup_cursor_table[] = { 40, 56, 80, 104, 140 };
static s8 setup_hostname[16];
static s8 setup_myname[16];
static s32 setup_oldtop;
static s32 setup_oldbottom;
static s32 setup_top;
static s32 setup_bottom;
static s32 m_net_cursor;
static s32 m_net_items;
static s32 options_cursor;
static s32 keys_cursor;
static s32 bind_grab;
static s32 new_cursor;
static s32 gamepad_cursor;
static s32 display_cursor;
static s32 graphics_cursor;
static s8 customwidthstr[16];
static s8 customheightstr[16];
static s32 newwinmode;
static s32 msgNumber;
static s32 m_quit_prevstate;
static bool wasInMenus;
static s32 lanConfig_port;
static s8 lanConfig_portname[6];
static s8 lanConfig_joinname[22];
static s32 lanConfig_cursor = -1;
static s32 lanConfig_cursor_table[] = { 72, 92, 124 };
static s32 startepisode;
static s32 startlevel;
static s32 maxplayers;
static bool m_serverInfoMessage = 0;
static f64 m_serverInfoMessageTime;
static s32 gameoptions_cursor_table[] = { 40, 56, 64, 72, 80, 88, 96, 112, 120 };
static s32 gameoptions_cursor;
static u8 identityTable[256];
static u8 translationTable[256];
static s32 fpslimits[] = { 30,60,72,100,120,144,165,180,200,240,300,360,500,999999 };
static s32 moden = -1;
static s32 display_sel_moden = 0;
static SDL_DisplayMode **modes = 0;
static SDL_DisplayMode *mode = 0;

level_t levels[] = {
	{ "start", "Entrance" }, // 0
	{ "e1m1", "Slipgate Complex" }, // 1
	{ "e1m2", "Castle of the Damned" },
	{ "e1m3", "The Necropolis" },
	{ "e1m4", "The Grisly Grotto" },
	{ "e1m5", "Gloom Keep" },
	{ "e1m6", "The Door To Chthon" },
	{ "e1m7", "The House of Chthon" },
	{ "e1m8", "Ziggurat Vertigo" },
	{ "e2m1", "The Installation" }, // 9
	{ "e2m2", "Ogre Citadel" },
	{ "e2m3", "Crypt of Decay" },
	{ "e2m4", "The Ebon Fortress" },
	{ "e2m5", "The Wizard's Manse" },
	{ "e2m6", "The Dismal Oubliette" },
	{ "e2m7", "Underearth" },
	{ "e3m1", "Termination Central" }, // 16
	{ "e3m2", "The Vaults of Zin" },
	{ "e3m3", "The Tomb of Terror" },
	{ "e3m4", "Satan's Dark Delight" },
	{ "e3m5", "Wind Tunnels" },
	{ "e3m6", "Chambers of Torment" },
	{ "e3m7", "The Haunted Halls" },
	{ "e4m1", "The Sewage System" }, // 23
	{ "e4m2", "The Tower of Despair" },
	{ "e4m3", "The Elder God Shrine" },
	{ "e4m4", "The Palace of Hate" },
	{ "e4m5", "Hell's Atrium" },
	{ "e4m6", "The Pain Maze" },
	{ "e4m7", "Azure Agony" },
	{ "e4m8", "The Nameless City" },
	{ "end", "Shub-Niggurath's Pit" }, // 31
	{ "dm1", "Place of Two Deaths" }, // 32
	{ "dm2", "Claustrophobopolis" },
	{ "dm3", "The Abandoned Base" },
	{ "dm4", "The Bad Place" },
	{ "dm5", "The Cistern" },
	{ "dm6", "The Dark Zone" }
};

level_t hipnoticlevels[] = { //MED 01/06/97 added hipnotic levels
	{ "start", "Command HQ" }, // 0
	{ "hip1m1", "The Pumping Station" }, // 1
	{ "hip1m2", "Storage Facility" },
	{ "hip1m3", "The Lost Mine" },
	{ "hip1m4", "Research Facility" },
	{ "hip1m5", "Military Complex" },
	{ "hip2m1", "Ancient Realms" }, // 6
	{ "hip2m2", "The Black Cathedral" },
	{ "hip2m3", "The Catacombs" },
	{ "hip2m4", "The Crypt" },
	{ "hip2m5", "Mortum's Keep" },
	{ "hip2m6", "The Gremlin's Domain" },
	{ "hip3m1", "Tur Torment" }, // 12
	{ "hip3m2", "Pandemonium" },
	{ "hip3m3", "Limbo" },
	{ "hip3m4", "The Gauntlet" },
	{ "hipend", "Armagon's Lair" }, // 16
	{ "hipdm1", "The Edge of Oblivion" } // 17
};

level_t roguelevels[] = { //PGM 01/07/97 added rogue levels
	{ "start", "Split Decision" }, //PGM 03/02/97 added dmatch level
	{ "r1m1", "Deviant's Domain" },
	{ "r1m2", "Dread Portal" },
	{ "r1m3", "Judgement Call" },
	{ "r1m4", "Cave of Death" },
	{ "r1m5", "Towers of Wrath" },
	{ "r1m6", "Temple of Pain" },
	{ "r1m7", "Tomb of the Overlord" },
	{ "r2m1", "Tempus Fugit" },
	{ "r2m2", "Elemental Fury I" },
	{ "r2m3", "Elemental Fury II" },
	{ "r2m4", "Curse of Osiris" },
	{ "r2m5", "Wizard's Keep" },
	{ "r2m6", "Blood Sacrifice" },
	{ "r2m7", "Last Bastion" },
	{ "r2m8", "Source of Evil" },
	{ "ctf1", "Division of Change" }
};

episode_t episodes[] = {
	{ "Welcome to Quake", 0, 1 },
	{ "Doomed Dimension", 1, 8 },
	{ "Realm of Black Magic", 9, 7 },
	{ "Netherworld", 16, 7 },
	{ "The Elder World", 23, 8 },
	{ "Final Level", 31, 1 },
	{ "Deathmatch Arena", 32, 6 }
};

episode_t hipnoticepisodes[] = {
	{ "Scourge of Armagon", 0, 1 },
	{ "Fortress of the Dead", 1, 5 },
	{ "Dominion of Darkness", 6, 6 },
	{ "The Rift", 12, 4 },
	{ "Final Level", 16, 1 },
	{ "Deathmatch Arena", 17, 1 }
};

episode_t rogueepisodes[] = { //PGM 01/07/97 added rogue episodes
	{ "Introduction", 0, 1 }, //PGM 03/02/97 added dmatch episode
	{ "Hell's Fortress", 1, 7 }, //MED 01/06/97  added hipnotic episodes
	{ "Corridors of Time", 8, 8 },
	{ "Deathmatch Arena", 16, 1 }
};

s8 *quitMessage[] = {
//	 .........1.........2....
	"  Are you gonna quit    ",
	"  this game just like   ",
	"   everything else?     ",
	"                        ",

	" Milord, methinks that  ",
	"   thou art a lowly     ",
	" quitter. Is this 1? ",
	"                        ",

	" Do I need to bust your ",
	"  face open for trying  ",
	"        to quit?        ",
	"                        ",

	" Man, I oughta smack you",
	"   for trying to quit!  ",
	"     Press Y to get     ",
	"      smacked out.      ",

	" Press Y to quit like a ",
	"   big loser in life.   ",
	"  Press N to stay proud ",
	"    and successful!     ",

	"   If you press Y to    ",
	"  quit, I will summon   ",
	"  Satan all over your   ",
	"      hard drive!       ",

	"  Um, Asmodeus dislikes ",
	" his children trying to ",
	" quit. Press Y to return",
	"   to your Tinkertoys.  ",

	"  If you quit now, I'll ",
	"  throw a blanket-party ",
	"   for you next time!   ",
	"                        "
};

enum { m_none, m_main, m_singleplayer, m_load, m_save, m_multiplayer, m_setup,
    m_net, m_options, m_video, m_keys, m_new, m_gamepad, m_display, m_graphics, 
    m_help, m_quit, m_lanconfig, m_gameoptions, m_search, m_slist } m_state;

void M_Menu_Main_f();
void M_Menu_SinglePlayer_f();
void M_Menu_Load_f();
void M_Menu_Save_f();
void M_Menu_MultiPlayer_f();
void M_Menu_Setup_f();
void M_Menu_Net_f();
void M_Menu_Options_f();
void M_Menu_Keys_f();
void M_Menu_Video_f();
void M_Menu_New_f();
void M_Menu_Help_f();
void M_Menu_Quit_f();
void M_Menu_LanConfig_f();
void M_Menu_GameOptions_f();
void M_Menu_Search_f();
void M_Menu_ServerList_f();
void M_Main_Draw();
void M_SinglePlayer_Draw();
void M_Load_Draw();
void M_Save_Draw();
void M_MultiPlayer_Draw();
void M_Setup_Draw();
void M_Net_Draw();
void M_Options_Draw();
void M_Keys_Draw();
void M_New_Draw();
void M_Video_Draw();
void M_Help_Draw();
void M_Quit_Draw();
void M_LanConfig_Draw();
void M_GameOptions_Draw();
void M_Search_Draw();
void M_ServerList_Draw();
void M_Main_Key(s32 key);
void M_SinglePlayer_Key(s32 key);
void M_Load_Key(s32 key);
void M_Save_Key(s32 key);
void M_MultiPlayer_Key(s32 key);
void M_Setup_Key(s32 key);
void M_Net_Key(s32 key);
void M_Options_Key(s32 key);
void M_Keys_Key(s32 key);
void M_Video_Key(s32 key);
void M_Help_Key(s32 key);
void M_Quit_Key(s32 key);
void M_LanConfig_Key(s32 key);
void M_GameOptions_Key(s32 key);
void M_ServerList_Key(s32 key);
void M_ConfigureNetSubsystem();
s32 help_page;
bool searchComplete = 0;
f64 searchCompleteTime;
s32 slist_cursor;
bool slist_sorted;
bool m_entersound; // play after drawing a frame, so caching won't disrupt the sound
bool m_recursiveDraw;
s32 m_return_state;
bool m_return_onerror;
s8 m_return_reason[32];

void M_DrawCharacter(s32 cx, s32 line, s32 num)
{ // Draws one solid graphics character
	drawlayer = lyr_menu.value;
	Draw_CharacterScaled(cx * uiscale + ((vid.width - 320 * uiscale) >> 1),
			     line * uiscale, num, uiscale);
	drawlayer = lyr_main.value;
}

void M_DrawCursor(s32 x, s32 y)
{
	M_DrawCharacter(x, y, 12 + ((s32)(realtime * 4) & 1));
}

void M_DrawCursorLine(s32 x, s32 y)
{
	M_DrawCharacter(x, y, 10 + ((s32)(realtime * 4) & 1));
}

void M_Print(s32 cx, s32 cy, s8 *str)
{
	while (*str) {
		M_DrawCharacter(cx, cy, (*str) + 128);
		str++;
		cx += 8;
	}
}

void M_PrintWhite(s32 cx, s32 cy, s8 *str)
{
	while (*str) {
		M_DrawCharacter(cx, cy, *str);
		str++;
		cx += 8;
	}
}

void M_DrawTransPic(s32 x, s32 y, qpic_t *pic)
{
	drawlayer = lyr_menu.value;
	Draw_TransPicScaled(x * uiscale + ((vid.width - 320 * uiscale) >> 1),
			    y * uiscale, pic, uiscale);
	drawlayer = lyr_main.value;
}

void M_BuildTranslationTable(s32 top, s32 bottom)
{
	for (s32 j = 0; j < 256; j++)
		identityTable[j] = j;
	u8 *dest = translationTable;
	u8 *source = identityTable;
	memcpy(dest, source, 256);
	if (top < 128) // the artists made some backwards ranges. sigh.
		memcpy(dest + TOP_RANGE, source + top, 16);
	else
		for (s32 j = 0; j < 16; j++)
			dest[TOP_RANGE + j] = source[top + 15 - j];
	if (bottom < 128)
		memcpy(dest + BOTTOM_RANGE, source + bottom, 16);
	else
		for (s32 j = 0; j < 16; j++)
			dest[BOTTOM_RANGE + j] = source[bottom + 15 - j];
}

void M_DrawTransPicTranslate(s32 x, s32 y, qpic_t *pic)
{
	drawlayer = lyr_menu.value;
	Draw_TransPicTranslateScaled(x * uiscale +
		((vid.width - 320 * uiscale) >> 1),
		y * uiscale, pic, translationTable, uiscale);
	drawlayer = lyr_main.value;
}

void M_DrawTextBox(s32 x, s32 y, s32 width, s32 lines)
{
	s32 cx = x; // draw left side
	s32 cy = y;
	qpic_t *p = Draw_CachePic("gfx/box_tl.lmp");
	M_DrawTransPic(cx, cy, p);
	p = Draw_CachePic("gfx/box_ml.lmp");
	for (s32 n = 0; n < lines; n++) {
		cy += 8;
		M_DrawTransPic(cx, cy, p);
	}
	p = Draw_CachePic("gfx/box_bl.lmp");
	M_DrawTransPic(cx, cy + 8, p);
	cx += 8; // draw middle
	while (width > 0) {
		cy = y;
		p = Draw_CachePic("gfx/box_tm.lmp");
		M_DrawTransPic(cx, cy, p);
		p = Draw_CachePic("gfx/box_mm.lmp");
		for (s32 n = 0; n < lines; n++) {
			cy += 8;
			if (n == 1)
				p = Draw_CachePic("gfx/box_mm2.lmp");
			M_DrawTransPic(cx, cy, p);
		}
		p = Draw_CachePic("gfx/box_bm.lmp");
		M_DrawTransPic(cx, cy + 8, p);
		width -= 2;
		cx += 16;
	}
	cy = y; // draw right side
	p = Draw_CachePic("gfx/box_tr.lmp");
	M_DrawTransPic(cx, cy, p);
	p = Draw_CachePic("gfx/box_mr.lmp");
	for (s32 n = 0; n < lines; n++) {
		cy += 8;
		M_DrawTransPic(cx, cy, p);
	}
	p = Draw_CachePic("gfx/box_br.lmp");
	M_DrawTransPic(cx, cy + 8, p);
}

void M_ToggleMenu_f()
{
	m_entersound = 1;
	if (key_dest == key_menu) {
		if (m_state != m_main) {
			M_Menu_Main_f();
			return;
		}
		key_dest = key_game;
		m_state = m_none;
		return;
	}
	if (key_dest == key_console)
		Con_ToggleConsole_f();
	else
		M_Menu_Main_f();
}

void M_Menu_Main_f()
{
	if (key_dest != key_menu) {
		m_save_demonum = cls.demonum;
		cls.demonum = -1;
	}
	key_dest = key_menu;
	m_state = m_main;
	m_entersound = 1;
}

void M_Main_Draw()
{
	M_DrawTransPic(16, 4, Draw_CachePic("gfx/qplaque.lmp"));
	qpic_t *p = Draw_CachePic("gfx/ttl_main.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	M_DrawTransPic(72, 32, Draw_CachePic("gfx/mainmenu.lmp"));
	s32 f = (s32)(realtime * 10) % 6;
	M_DrawTransPic(54, 32 + m_main_cursor * 20,
	       Draw_CachePic(va("gfx/menudot%i.lmp", f + 1)));
}

void M_Main_Key(s32 key)
{
	switch (key) {
	case K_ESCAPE:
		key_dest = key_game;
		m_state = m_none;
		cls.demonum = m_save_demonum;
		if (cls.demonum != -1 && !cls.demoplayback
		    && cls.state != ca_connected)
			CL_NextDemo();
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		if (++m_main_cursor >= MAIN_ITEMS)
			m_main_cursor = 0;
		break;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		if (--m_main_cursor < 0)
			m_main_cursor = MAIN_ITEMS - 1;
		break;
	case K_ENTER:
		m_entersound = 1;
		switch (m_main_cursor) {
		case 0:
			M_Menu_SinglePlayer_f();
			break;
		case 1:
			M_Menu_MultiPlayer_f();
			break;
		case 2:
			M_Menu_Options_f();
			break;
		case 3:
			M_Menu_Help_f();
			break;
		case 4:
			M_Menu_Quit_f();
			break;
		}
	}
}

void M_Menu_SinglePlayer_f()
{
	key_dest = key_menu;
	m_state = m_singleplayer;
	m_entersound = 1;
}

void M_SinglePlayer_Draw()
{
	M_DrawTransPic(16, 4, Draw_CachePic("gfx/qplaque.lmp"));
	qpic_t *p = Draw_CachePic("gfx/ttl_sgl.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	M_DrawTransPic(72, 32, Draw_CachePic("gfx/sp_menu.lmp"));
	s32 f = (s32)(realtime * 10) % 6;
	M_DrawTransPic(54, 32 + m_singleplayer_cursor * 20,
		       Draw_CachePic(va("gfx/menudot%i.lmp", f + 1)));
}

void M_SinglePlayer_Key(s32 key)
{
	switch (key) {
	case K_ESCAPE:
		M_Menu_Main_f();
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		if (++m_singleplayer_cursor >= SINGLEPLAYER_ITEMS)
			m_singleplayer_cursor = 0;
		break;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		if (--m_singleplayer_cursor < 0)
			m_singleplayer_cursor = SINGLEPLAYER_ITEMS - 1;
		break;
	case K_ENTER:
		m_entersound = 1;
		switch (m_singleplayer_cursor) {
		case 0:
			if (sv.active)
				if (!SCR_ModalMessage
				    ("Are you sure you want to\nstart a new game?\n"))
					break;
			key_dest = key_game;
			if (sv.active)
				Cbuf_AddText("disconnect\n");
			Cbuf_AddText("maxplayers 1\n");
			Cbuf_AddText("map start\n");
			break;
		case 1:
			M_Menu_Load_f();
			break;
		case 2:
			M_Menu_Save_f();
			break;
		}
	}
}

void M_ScanSaves()
{
	for (s32 i = 0; i < MAX_SAVEGAMES; i++) {
		strcpy(m_filenames[i], "--- UNUSED SLOT ---");
		loadable[i] = 0;
		s8 name[MAX_OSPATH*2];
		sprintf(name, "%s/s%i.sav", com_gamedir, i);
		FILE *f = fopen(name, "r");
		if (!f) continue;
		s32 version;
		fscanf(f, "%i\n", &version);
		fscanf(f, "%79s\n", name);
		strncpy(m_filenames[i], name, sizeof(m_filenames[i]) - 1);
		for (s32 j = 0; j < SAVEGAME_COMMENT_LENGTH; j++)
			if (m_filenames[i][j] == '_') // change _ back to space
				m_filenames[i][j] = ' ';
		loadable[i] = 1;
		fclose(f);
	}
}

void M_Menu_Load_f()
{
	m_entersound = 1;
	m_state = m_load;
	key_dest = key_menu;
	M_ScanSaves();
}

void M_Menu_Save_f()
{
	if (!sv.active || cl.intermission || svs.maxclients != 1)
		return;
	m_entersound = 1;
	m_state = m_save;
	key_dest = key_menu;
	M_ScanSaves();
}

void M_Load_Draw()
{
	qpic_t *p = Draw_CachePic("gfx/p_load.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	for (s32 i = 0; i < MAX_SAVEGAMES; i++)
		M_Print(16, 32 + 8 * i, m_filenames[i]);
	M_DrawCursor(8, 32 + load_cursor * 8);
}

void M_Save_Draw()
{
	qpic_t *p = Draw_CachePic("gfx/p_save.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	for (s32 i = 0; i < MAX_SAVEGAMES; i++)
		M_Print(16, 32 + 8 * i, m_filenames[i]);
	M_DrawCursor(8, 32 + load_cursor * 8);
}

void M_Load_Key(s32 k)
{
	switch (k) {
	case K_ESCAPE:
		M_Menu_SinglePlayer_f();
		break;
	case K_ENTER:
		S_LocalSound("misc/menu2.wav");
		if (!loadable[load_cursor])
			return;
		m_state = m_none;
		key_dest = key_game;
		// Host_Loadgame_f can't bring up the loading plaque because
		// too much stack space has been used, so do it now
		SCR_BeginLoadingPlaque();
		// issue the load command
		Cbuf_AddText(va("load s%i\n", load_cursor));
		return;
	case K_UPARROW:
	case K_LEFTARROW:
		S_LocalSound("misc/menu1.wav");
		load_cursor--;
		if (load_cursor < 0)
			load_cursor = MAX_SAVEGAMES - 1;
		break;
	case K_DOWNARROW:
	case K_RIGHTARROW:
		S_LocalSound("misc/menu1.wav");
		load_cursor++;
		if (load_cursor >= MAX_SAVEGAMES)
			load_cursor = 0;
		break;
	}
}

void M_Save_Key(s32 k)
{
	switch (k) {
	case K_ESCAPE:
		M_Menu_SinglePlayer_f();
		break;
	case K_ENTER:
		m_state = m_none;
		key_dest = key_game;
		Cbuf_AddText(va("save s%i\n", load_cursor));
		return;
	case K_UPARROW:
	case K_LEFTARROW:
		S_LocalSound("misc/menu1.wav");
		load_cursor--;
		if (load_cursor < 0)
			load_cursor = MAX_SAVEGAMES - 1;
		break;
	case K_DOWNARROW:
	case K_RIGHTARROW:
		S_LocalSound("misc/menu1.wav");
		load_cursor++;
		if (load_cursor >= MAX_SAVEGAMES)
			load_cursor = 0;
		break;
	}
}

void M_Menu_MultiPlayer_f()
{
	key_dest = key_menu;
	m_state = m_multiplayer;
	m_entersound = 1;
}

void M_MultiPlayer_Draw()
{
	M_DrawTransPic(16, 4, Draw_CachePic("gfx/qplaque.lmp"));
	qpic_t *p = Draw_CachePic("gfx/p_multi.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	M_DrawTransPic(72, 32, Draw_CachePic("gfx/mp_menu.lmp"));
	s32 f = (s32)(realtime * 10) % 6;
	M_DrawTransPic(54, 32 + m_multiplayer_cursor * 20,
		       Draw_CachePic(va("gfx/menudot%i.lmp", f + 1)));
	if (tcpipAvailable)
		return;
	M_PrintWhite((320 / 2) - ((27 * 8) / 2), 148,
		     "No Communications Available");
}

void M_MultiPlayer_Key(s32 key)
{
	switch (key) {
	case K_ESCAPE:
		M_Menu_Main_f();
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		if (++m_multiplayer_cursor >= MULTIPLAYER_ITEMS)
			m_multiplayer_cursor = 0;
		break;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		if (--m_multiplayer_cursor < 0)
			m_multiplayer_cursor = MULTIPLAYER_ITEMS - 1;
		break;
	case K_ENTER:
		m_entersound = 1;
		switch (m_multiplayer_cursor) {
		case 0:
			if (tcpipAvailable)
				M_Menu_Net_f();
			break;
		case 1:
			if (tcpipAvailable)
				M_Menu_Net_f();
			break;
		case 2:
			M_Menu_Setup_f();
			break;
		}
	}
}

void M_Menu_Setup_f()
{
	key_dest = key_menu;
	m_state = m_setup;
	m_entersound = 1;
	Q_strcpy(setup_myname, cl_name.string);
	Q_strcpy(setup_hostname, hostname.string);
	setup_top = setup_oldtop = ((s32)cl_color.value) >> 4;
	setup_bottom = setup_oldbottom = ((s32)cl_color.value) & 15;
}

void M_Setup_Draw()
{
	M_DrawTransPic(16, 4, Draw_CachePic("gfx/qplaque.lmp"));
	qpic_t *p = Draw_CachePic("gfx/p_multi.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	M_Print(64, 40, "Hostname");
	M_DrawTextBox(160, 32, 16, 1);
	M_Print(168, 40, setup_hostname);
	M_Print(64, 56, "Your name");
	M_DrawTextBox(160, 48, 16, 1);
	M_Print(168, 56, setup_myname);
	M_Print(64, 80, "Shirt color");
	M_Print(64, 104, "Pants color");
	M_DrawTextBox(64, 140 - 8, 14, 1);
	M_Print(72, 140, "Accept Changes");
	p = Draw_CachePic("gfx/bigbox.lmp");
	M_DrawTransPic(160, 64, p);
	p = Draw_CachePic("gfx/menuplyr.lmp");
	M_BuildTranslationTable(setup_top * 16, setup_bottom * 16);
	M_DrawTransPicTranslate(172, 72, p);
	M_DrawCharacter(56, setup_cursor_table[setup_cursor],
			12 + ((s32)(realtime * 4) & 1));
	if (setup_cursor == 0)
		M_DrawCharacter(168 + 8 * strlen(setup_hostname),
				setup_cursor_table[setup_cursor],
				10 + ((s32)(realtime * 4) & 1));
	if (setup_cursor == 1)
		M_DrawCharacter(168 + 8 * strlen(setup_myname),
				setup_cursor_table[setup_cursor],
				10 + ((s32)(realtime * 4) & 1));
}

void M_Setup_Key(s32 k)
{
	switch (k) {
	case K_ESCAPE:
		M_Menu_MultiPlayer_f();
		break;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		setup_cursor--;
		if (setup_cursor < 0)
			setup_cursor = NUM_SETUP_CMDS - 1;
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		setup_cursor++;
		if (setup_cursor >= NUM_SETUP_CMDS)
			setup_cursor = 0;
		break;
	case K_LEFTARROW:
		if (setup_cursor < 2)
			return;
		S_LocalSound("misc/menu3.wav");
		if (setup_cursor == 2)
			setup_top = setup_top - 1;
		if (setup_cursor == 3)
			setup_bottom = setup_bottom - 1;
		break;
	case K_RIGHTARROW:
		if (setup_cursor < 2)
			return;
forward:
		S_LocalSound("misc/menu3.wav");
		if (setup_cursor == 2)
			setup_top = setup_top + 1;
		if (setup_cursor == 3)
			setup_bottom = setup_bottom + 1;
		break;
	case K_ENTER:
		if (setup_cursor == 0 || setup_cursor == 1)
			return;
		if (setup_cursor == 2 || setup_cursor == 3)
			goto forward;
		// setup_cursor == 4 (OK)
		if (Q_strcmp(cl_name.string, setup_myname) != 0)
			Cbuf_AddText(va("name \"%s\"\n", setup_myname));
		if (Q_strcmp(hostname.string, setup_hostname) != 0)
			Cvar_Set("hostname", setup_hostname);
		if (setup_top != setup_oldtop
		    || setup_bottom != setup_oldbottom)
			Cbuf_AddText(va
				     ("color %i %i\n", setup_top,
				      setup_bottom));
		m_entersound = 1;
		M_Menu_MultiPlayer_f();
		break;
	case K_BACKSPACE:
		if (setup_cursor == 0) {
			if (strlen(setup_hostname))
				setup_hostname[strlen(setup_hostname) - 1] = 0;
		}

		if (setup_cursor == 1) {
			if (strlen(setup_myname))
				setup_myname[strlen(setup_myname) - 1] = 0;
		}
		break;
	default:
		if (k < 32 || k > 127)
			break;
		if (setup_cursor == 0) {
			s32 l = strlen(setup_hostname);
			if (l < 15) {
				setup_hostname[l + 1] = 0;
				setup_hostname[l] = k;
			}
		}
		if (setup_cursor == 1) {
			s32 l = strlen(setup_myname);
			if (l < 15) {
				setup_myname[l + 1] = 0;
				setup_myname[l] = k;
			}
		}
	}
	if (setup_top > 13)
		setup_top = 0;
	if (setup_top < 0)
		setup_top = 13;
	if (setup_bottom > 13)
		setup_bottom = 0;
	if (setup_bottom < 0)
		setup_bottom = 13;
}

s8 *net_helpMessage[] = {
//	 .........1.........2....
	"                        ",
	" Two computers connected",
	"   through two modems.  ",
	"                        ",

	"                        ",
	" Two computers connected",
	" by a null-modem cable. ",
	"                        ",

	" Novell network LANs    ",
	" or Windows 95 DOS-box. ",
	"                        ",
	"(LAN=Local Area Network)",

	" Commonly used to play  ",
	" over the Internet, but ",
	" also used on a Local   ",
	" Area Network.          "
};

void M_Menu_Net_f()
{
	key_dest = key_menu;
	m_state = m_net;
	m_entersound = 1;
	m_net_items = 4;
	if (m_net_cursor >= m_net_items)
		m_net_cursor = 0;
	m_net_cursor--;
	M_Net_Key(K_DOWNARROW);
}

void M_Net_Draw()
{
	M_DrawTransPic(16, 4, Draw_CachePic("gfx/qplaque.lmp"));
	qpic_t *p = Draw_CachePic("gfx/p_multi.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	s32 f = 32;
	p = Draw_CachePic("gfx/dim_modm.lmp");
	if (p)
		M_DrawTransPic(72, f, p);
	f += 19;
	p = Draw_CachePic("gfx/dim_drct.lmp");
	if (p)
		M_DrawTransPic(72, f, p);
	f += 19;
	p = Draw_CachePic("gfx/dim_ipx.lmp");
	M_DrawTransPic(72, f, p);
	f += 19;
	if (tcpipAvailable)
		p = Draw_CachePic("gfx/netmen4.lmp");
	else
		p = Draw_CachePic("gfx/dim_tcp.lmp");
	M_DrawTransPic(72, f, p);
	if (m_net_items == 5) { // JDC, could just be removed
		f += 19;
		p = Draw_CachePic("gfx/netmen5.lmp");
		M_DrawTransPic(72, f, p);
	}
	f = (320 - 26 * 8) / 2;
	M_DrawTextBox(f, 134, 24, 4);
	f += 8;
	M_Print(f, 142, net_helpMessage[m_net_cursor * 4 + 0]);
	M_Print(f, 150, net_helpMessage[m_net_cursor * 4 + 1]);
	M_Print(f, 158, net_helpMessage[m_net_cursor * 4 + 2]);
	M_Print(f, 166, net_helpMessage[m_net_cursor * 4 + 3]);
	f = (s32)(realtime * 10) % 6;
	M_DrawTransPic(54, 32 + m_net_cursor * 20,
		       Draw_CachePic(va("gfx/menudot%i.lmp", f + 1)));
}

void M_Net_Key(s32 k)
{
again:
	switch (k) {
	case K_ESCAPE:
		M_Menu_MultiPlayer_f();
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		if (++m_net_cursor >= m_net_items)
			m_net_cursor = 0;
		break;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		if (--m_net_cursor < 0)
			m_net_cursor = m_net_items - 1;
		break;
	case K_ENTER:
		m_entersound = 1;
		switch (m_net_cursor) {
		case 0:
			break;

		case 1:
			break;

		case 2:
			M_Menu_LanConfig_f();
			break;

		case 3:
			M_Menu_LanConfig_f();
			break;

		case 4:
			break; // multiprotocol
		}
	}
	if (m_net_cursor == 0)
		goto again;
	if (m_net_cursor == 1)
		goto again;
	if (m_net_cursor == 2)
		goto again;
	if (m_net_cursor == 3 && !tcpipAvailable)
		goto again;
}

void M_Menu_Options_f()
{
	key_dest = key_menu;
	m_state = m_options;
	m_entersound = 1;
	drawmousemenu = !(SDLWindowFlags & SDL_WINDOW_FULLSCREEN);
	if (options_cursor == 13
	    && drawmousemenu && !_windowed_mouse.value && !newoptions.value) {
		options_cursor = 0;
	}
}

void M_AdjustSliders(s32 dir)
{
	S_LocalSound("misc/menu3.wav");
	switch (options_cursor) {
	case 3: // screen size
		scr_viewsize.value += dir * 10;
		if (scr_viewsize.value < 30)
			scr_viewsize.value = 30;
		if (scr_viewsize.value > 120)
			scr_viewsize.value = 120;
		Cvar_SetValue("viewsize", scr_viewsize.value);
		break;
	case 4: // gamma
		v_gamma.value -= dir * 0.05;
		if (v_gamma.value < 0.5)
			v_gamma.value = 0.5;
		if (v_gamma.value > 1)
			v_gamma.value = 1;
		Cvar_SetValue("gamma", v_gamma.value);
		break;
	case 5: // mouse speed
		sensitivity.value += dir * 0.5;
		if (sensitivity.value < 1)
			sensitivity.value = 1;
		if (sensitivity.value > 11)
			sensitivity.value = 11;
		Cvar_SetValue("sensitivity", sensitivity.value);
		break;
	case 6: // music volume
		bgmvolume.value += dir * 0.1;
		if (bgmvolume.value < 0)
			bgmvolume.value = 0;
		if (bgmvolume.value > 1)
			bgmvolume.value = 1;
		Cvar_SetValue("bgmvolume", bgmvolume.value);
		break;
	case 7: // sfx volume
		sfxvolume.value += dir * 0.1;
		if (sfxvolume.value < 0)
			sfxvolume.value = 0;
		if (sfxvolume.value > 1)
			sfxvolume.value = 1;
		Cvar_SetValue("volume", sfxvolume.value);
		break;
	case 8: // allways run
		if (cl_forwardspeed.value > 200) {
			Cvar_SetValue("cl_forwardspeed", 200);
			Cvar_SetValue("cl_backspeed", 200);
		} else {
			Cvar_SetValue("cl_forwardspeed", 400);
			Cvar_SetValue("cl_backspeed", 400);
		}
		break;
	case 9: // invert mouse
		Cvar_SetValue("m_pitch", -m_pitch.value);
		break;
	case 10: // lookspring
		Cvar_SetValue("lookspring", !lookspring.value);
		break;
	case 11: // lookstrafe
		Cvar_SetValue("lookstrafe", !lookstrafe.value);
		break;
	case 13: // _windowed_mouse
		Cvar_SetValue("_windowed_mouse", !_windowed_mouse.value);
		break;
	}
}

void M_DrawSlider(s32 x, s32 y, f32 range)
{
	if (range < 0)
		range = 0;
	if (range > 1)
		range = 1;
	M_DrawCharacter(x - 8, y, 128);
	s32 i;
	for (i = 0; i < SLIDER_RANGE; i++)
		M_DrawCharacter(x + i * 8, y, 129);
	M_DrawCharacter(x + i * 8, y, 130);
	M_DrawCharacter(x + (SLIDER_RANGE - 1) * 8 * range, y, 131);
}

void M_DrawCheckbox(s32 x, s32 y, s32 on)
{
	M_Print(x, y, on ? "on" : "off");
}

void M_Options_Draw()
{
	M_DrawTransPic(16, 4, Draw_CachePic("gfx/qplaque.lmp"));
	qpic_t *p = Draw_CachePic("gfx/p_option.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	M_Print(16, 32, "    Customize controls");
	M_Print(16, 40, "         Go to console");
	M_Print(16, 48, "     Reset to defaults");
	M_Print(16, 56, "           Screen size");
	f32 r = (scr_viewsize.value - 30) / (120 - 30);
	M_DrawSlider(220, 56, r);
	M_Print(16, 64, "            Brightness");
	r = (1.0 - v_gamma.value) / 0.5;
	M_DrawSlider(220, 64, r);
	M_Print(16, 72, "           Mouse Speed");
	r = (sensitivity.value - 1) / 10;
	M_DrawSlider(220, 72, r);
	M_Print(16, 80, "       CD Music Volume");
	r = bgmvolume.value;
	M_DrawSlider(220, 80, r);
	M_Print(16, 88, "          Sound Volume");
	r = sfxvolume.value;
	M_DrawSlider(220, 88, r);
	M_Print(16, 96, "            Always Run");
	M_DrawCheckbox(220, 96, cl_forwardspeed.value > 200);
	M_Print(16, 104, "          Invert Mouse");
	M_DrawCheckbox(220, 104, m_pitch.value < 0);
	M_Print(16, 112, "            Lookspring");
	M_DrawCheckbox(220, 112, lookspring.value);
	M_Print(16, 120, "            Lookstrafe");
	M_DrawCheckbox(220, 120, lookstrafe.value);
	M_Print(16, 128, "         Video Options");
	if (drawmousemenu) {
		M_Print(16, 136, "             Use Mouse");
		M_DrawCheckbox(220, 136, _windowed_mouse.value);
	}
	if (newoptions.value)
		M_Print(16, drawmousemenu ? 144:136, "           New Options");
	M_DrawCursor(200, 32 + options_cursor * 8);
}

void M_Options_Key(s32 k)
{
	drawmousemenu = !(SDLWindowFlags & SDL_WINDOW_FULLSCREEN);
	switch (k) {
	case K_ESCAPE:
		M_Menu_Main_f();
		break;
	case K_ENTER:
		m_entersound = 1;
		if (options_cursor == 0) {
			M_Menu_Keys_f();
		} else if (options_cursor == 1) {
			m_state = m_none;
			Con_ToggleConsole_f();
		} else if (options_cursor == 2) {
			Cbuf_AddText("exec default.cfg\n");
		} else if (options_cursor == 12) {
			M_Menu_Video_f();
		} else if (options_cursor == 13) {
			if (drawmousemenu)
				M_AdjustSliders(1);
			else
				M_Menu_New_f();
		} else if (options_cursor == 14) {
			M_Menu_New_f();
		} else {
			M_AdjustSliders(1);
		}
		return;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		options_cursor--;
		if (options_cursor < 0)
			options_cursor =
			    12 + (drawmousemenu & 1) +
			    ((s32)newoptions.value & 1);
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		options_cursor++;
		if (options_cursor >=
		    13 + (drawmousemenu & 1) + ((s32)newoptions.value & 1))
			options_cursor = 0;
		break;
	case K_LEFTARROW:
		M_AdjustSliders(-1);
		break;
	case K_RIGHTARROW:
		M_AdjustSliders(1);
		break;
	}
}

s8 *bindnames[][2] = {
	{ "+attack", "attack" },
	{ "impulse 10", "change weapon" },
	{ "+jump", "jump / swim up" },
	{ "+forward", "walk forward" },
	{ "+back", "backpedal" },
	{ "+left", "turn left" },
	{ "+right", "turn right" },
	{ "+speed", "run" },
	{ "+moveleft", "step left" },
	{ "+moveright", "step right" },
	{ "+strafe", "sidestep" },
	{ "+lookup", "look up" },
	{ "+lookdown", "look down" },
	{ "centerview", "center view" },
	{ "+mlook", "mouse look" },
	{ "+klook", "keyboard look" },
	{ "+moveup", "swim up" },
	{ "+movedown", "swim down" }
};

void M_Menu_Keys_f()
{
	key_dest = key_menu;
	m_state = m_keys;
	m_entersound = 1;
}

void M_FindKeysForCommand(s8 *command, s32 *twokeys)
{
	twokeys[0] = twokeys[1] = -1;
	s32 l = strlen(command);
	s32 count = 0;
	for (s32 j = 0; j < 256; j++) {
		s8 *b = keybindings[j];
		if (!b)
			continue;
		if (!strncmp(b, command, l)) {
			twokeys[count] = j;
			count++;
			if (count == 2)
				break;
		}
	}
}

void M_UnbindCommand(s8 *command)
{
	s32 l = strlen(command);
	for (s32 j = 0; j < 256; j++) {
		s8 *b = keybindings[j];
		if (!b)
			continue;
		if (!strncmp(b, command, l))
			Key_SetBinding(j, "");
	}
}

void M_Keys_Draw()
{
	qpic_t *p = Draw_CachePic("gfx/ttl_cstm.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	if (bind_grab)
		M_Print(12, 32, "Press a key or button for this action");
	else
		M_Print(18, 32, "Enter to change, backspace to clear");
	for (s32 i = 0; i < NUMCOMMANDS; i++) { // search for known bindings
		s32 y = 48 + 8 * i;
		M_Print(16, y, bindnames[i][1]);
		s32 keys[2];
		M_FindKeysForCommand(bindnames[i][0], keys);
		if (keys[0] == -1) {
			M_Print(140, y, "???");
		} else {
			s8 *name = Key_KeynumToString(keys[0]);
			M_Print(140, y, name);
			s32 x = strlen(name) * 8;
			if (keys[1] != -1) {
				M_Print(140 + x + 8, y, "or");
				M_Print(140 + x + 32, y,
					Key_KeynumToString(keys[1]));
			}
		}
	}
	if (bind_grab)
		M_DrawCharacter(130, 48 + keys_cursor * 8, '=');
	else
		M_DrawCharacter(130, 48 + keys_cursor * 8,
				12 + ((s32)(realtime * 4) & 1));
}

void M_Keys_Key(s32 k)
{
	if (bind_grab) { // defining a key
		S_LocalSound("misc/menu1.wav");
		if (k == K_ESCAPE) {
			bind_grab = 0;
		} else if (k != '`') {
			s8 cmd[80];
			sprintf(cmd, "bind \"%s\" \"%s\"\n",
				Key_KeynumToString(k),
				bindnames[keys_cursor][0]);
			Cbuf_InsertText(cmd);
		}
		bind_grab = 0;
		return;
	}
	s32 keys[2]; // keep here for the OpenBSD compiler
	switch (k) {
	case K_ESCAPE:
		M_Menu_Options_f();
		break;
	case K_LEFTARROW:
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		keys_cursor--;
		if (keys_cursor < 0)
			keys_cursor = NUMCOMMANDS - 1;
		break;
	case K_DOWNARROW:
	case K_RIGHTARROW:
		S_LocalSound("misc/menu1.wav");
		keys_cursor++;
		if (keys_cursor >= NUMCOMMANDS)
			keys_cursor = 0;
		break;
	case K_ENTER: // go into bind mode
		M_FindKeysForCommand(bindnames[keys_cursor][0], keys);
		S_LocalSound("misc/menu2.wav");
		if (keys[1] != -1)
			M_UnbindCommand(bindnames[keys_cursor][0]);
		bind_grab = 1;
		break;
	case K_BACKSPACE:
	case K_DEL: // delete bindings
		S_LocalSound("misc/menu2.wav");
		M_UnbindCommand(bindnames[keys_cursor][0]);
		break;
	}
}

void M_Menu_Gamepad_f()
{
	key_dest = key_menu;
	m_state = m_gamepad;
	m_entersound = 1;
}

void M_Gamepad_Draw()
{
	s8 temp[28];
	s32 xoffs = 16;
	qpic_t *p = Draw_CachePic("gfx/ttl_cstm.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	M_Print(xoffs, 32, "Device");
	s32 count = -1;
	SDL_JoystickID *joys = SDL_GetJoysticks(&count);
	if (count <= 0) { M_Print(xoffs + 72, 32, "None"); return; }
	else {
		if (joysticknum.value > count - 1 || joysticknum.value < 0)
			joysticknum.value = 0;
		q_strlcpy(temp, SDL_GetJoystickNameForID(joys[(s32)joysticknum.value]), 28);
		M_Print(xoffs + 72, 32, temp);
	}
	M_Print(xoffs, 48,  "         Deadzone");
	M_Print(xoffs, 56,  "       Move speed");
	M_Print(xoffs, 64,  "       Look speed");
	M_Print(xoffs, 72,  "Trigger threshold");
	M_Print(xoffs, 80,  "      Move X axis");
	M_Print(xoffs, 88,  "      Move Y axis");
	M_Print(xoffs, 96,  "      Look X axis");
	M_Print(xoffs, 104, "      Look Y axis");
	M_Print(xoffs, 112, "   Trigger 1 axis");
	M_Print(xoffs, 120, "   Trigger 2 axis");
	M_Print(xoffs, 128, "     Enter button");
	M_Print(xoffs, 136, "    Escape button");
	xoffs = 184; 
	M_DrawSlider(xoffs, 48, (jdeadzone.value/(1<<14)));
	M_DrawSlider(xoffs, 56, ((jmovesens.value-0.5)/4.5));
	M_DrawSlider(xoffs, 64, ((jlooksens.value-0.5)/4.5));
	M_DrawSlider(xoffs, 72, ((jtriggerthresh.value+(1<<15))/(1<<16)));
	sprintf(temp, "%d\n", (s32)jmoveaxisx.value);
	M_Print(xoffs, 80, temp);
	sprintf(temp, "%d\n", (s32)jmoveaxisy.value);
	M_Print(xoffs, 88, temp);
	sprintf(temp, "%d\n", (s32)jlookaxisx.value);
	M_Print(xoffs, 96, temp);
	sprintf(temp, "%d\n", (s32)jlookaxisy.value);
	M_Print(xoffs, 104, temp);
	sprintf(temp, "%d\n", (s32)jtrigaxis1.value);
	M_Print(xoffs, 112, temp);
	sprintf(temp, "%d\n", (s32)jtrigaxis2.value);
	M_Print(xoffs, 120, temp);
	sprintf(temp, "%d\n", (s32)jenterbutton.value);
	M_Print(xoffs, 128, temp);
	sprintf(temp, "%d\n", (s32)jescapebutton.value);
	M_Print(xoffs, 136, temp);
	if (!gamepad_cursor) M_DrawCursor(80, 32);
	else M_DrawCursor(168, 40 + gamepad_cursor*8);
	M_Print(24, 152, "Move");
	M_DrawTextBox(16, 156, 3, 4);
	M_Print(72, 152, "Look");
	M_DrawTextBox(64, 156, 3, 4);
	f32 movex = jaxis_move_x, movey = jaxis_move_y;
	f32 lookx = jaxis_look_x, looky = jaxis_look_y;
	movex = (movex / (1<<16)) * 24;
	movey = (movey / (1<<16)) * 24;
	lookx = (lookx / (1<<16)) * 24;
	looky = (looky / (1<<16)) * 24;
	if ((!movex && !movey) || movex >= 11.9 || movex <= -11.9 || 
				movey >= 11.9 || movey <= -11.9)
		M_PrintWhite(36+movex, 176+movey, "+");
	else M_Print(36+movex, 176+movey, "+");
	if ((!lookx && !looky) || lookx >= 11.9 || lookx <= -11.9 || 
				looky >= 11.9 || looky <= -11.9)
		M_PrintWhite(84+lookx, 176+looky, "+");
	else M_Print(84+lookx, 176+looky, "+");
	if (jaxis_trig_1 > jtriggerthresh.value)M_PrintWhite(120, 152, "Trig1");
	else M_Print(120, 152, "Trig1");
	if (jaxis_trig_2 > jtriggerthresh.value)M_PrintWhite(120, 160, "Trig2");
	else M_Print(120, 160, "Trig2");
	M_DrawSlider(176, 152, ((f32)(jaxis_trig_1+(1<<15))/(1<<16)));
	M_DrawSlider(176, 160, ((f32)(jaxis_trig_2+(1<<15))/(1<<16)));
	s32 total = SDL_GetNumJoystickButtons(joystick);
	q_strlcpy(temp, "Held buttons:", 14);
	for (s32 i = 0; i < total; i++) {
		if (SDL_GetJoystickButton(joystick, i)) {
			s8 btn[4];
			snprintf(btn, 4, " %d", i);
			q_strlcat(temp, btn, 28);
		}
	}
	M_Print(120, 168, temp);
}

void M_Display_Key(s32 k)
{
	switch (k) {
	case K_ESCAPE:
		M_Menu_New_f();
		break;
	case K_LEFTARROW:
		S_LocalSound("misc/menu3.wav");
		if (display_cursor == 0 && scr_uiscale.value > 1)
			Cvar_SetValue("scr_uiscale", scr_uiscale.value - 1);
		else if (display_cursor == 1) Cvar_SetValue("scr_lockuiscale",
				scr_lockuiscale.value ? 0 : 1);
		else if (display_cursor == 2) {
			if (scr_hudstyle.value == 0)
				Cvar_SetValue("hudstyle", 9);
			else
				Cvar_SetValue("hudstyle", scr_hudstyle.value-1);
		} else if (display_cursor == 3) {
			if (crosshair.value == 0) {
				Cvar_SetValue("crosshair", 1);
				Cvar_SetValue("cl_crosschar", 'o');
			} else if (cl_crosschar.value == 'o')
				Cvar_SetValue("cl_crosschar", 'x');
			else if (cl_crosschar.value == 'x')
				Cvar_SetValue("cl_crosschar", '*');
			else if (cl_crosschar.value == '*')
				Cvar_SetValue("cl_crosschar", 5);
			else if (cl_crosschar.value == 5)
				Cvar_SetValue("cl_crosschar", '+');
			else if (cl_crosschar.value == '+')
				Cvar_SetValue("crosshair", 0);
		} else if (display_cursor == 4)
			Cvar_SetValue("aspectr", aspectr.value - 0.01);
		else if (display_cursor == 5)
			Cvar_SetValue("scr_uixscale", 
					CLAMP(0, scr_uixscale.value-0.01, 1));
		else if (display_cursor == 6) {
			s32 i = 0;
			for (; i < (s32)sizeof(fpslimits) / 4 - 1; ++i)
				if (fpslimits[i] >= (s32)host_maxfps.value)
					break;
			--i;
			if (i < 0 || i > (s32)sizeof(fpslimits) / 4 - 1)
				i = sizeof(fpslimits) / 4 - 1;
			Cvar_SetValue("host_maxfps", fpslimits[i]);
		} else if (display_cursor == 7) Cvar_SetValue("scr_showfps",
				scr_showfps.value ? 0 : 1);
		else if (display_cursor == 8)
			Cvar_SetValue("fov", scr_fov.value - 1);
		else if (display_cursor == 9)
			Cvar_SetValue("yaspectscale", yaspectscale.value-0.01);
		else if (display_cursor == 10) {
			if (newwinmode == 0) newwinmode = 2;
			else newwinmode--;
		} else if (display_cursor == 11 && newwinmode == 1) {
			if (display_sel_moden == moden-1) display_sel_moden = 0;
			else display_sel_moden++;
		}
		break;
	case K_RIGHTARROW:
	case K_ENTER:
		S_LocalSound("misc/menu3.wav");
		if (display_cursor == 0 && scr_uiscale.value < (vid.width / 320))
			Cvar_SetValue("scr_uiscale", scr_uiscale.value + 1);
		else if (display_cursor == 1) Cvar_SetValue("scr_lockuiscale",
				scr_lockuiscale.value ? 0 : 1);
		else if (display_cursor == 2) {
			if (scr_hudstyle.value == 9)
				Cvar_SetValue("hudstyle", 0);
			else
				Cvar_SetValue("hudstyle", scr_hudstyle.value+1);
		} else if (display_cursor == 3) {
			if (crosshair.value == 0) {
				Cvar_SetValue("crosshair", 1);
				Cvar_SetValue("cl_crosschar", '+');
			} else if (cl_crosschar.value == '+')
				Cvar_SetValue("cl_crosschar", 5);
			else if (cl_crosschar.value == 5)
				Cvar_SetValue("cl_crosschar", '*');
			else if (cl_crosschar.value == '*')
				Cvar_SetValue("cl_crosschar", 'x');
			else if (cl_crosschar.value == 'x')
				Cvar_SetValue("cl_crosschar", 'o');
			else if (cl_crosschar.value == 'o')
				Cvar_SetValue("crosshair", 0);
		} else if (display_cursor == 4)
			Cvar_SetValue("aspectr", aspectr.value + 0.01);
		else if (display_cursor == 5)
			Cvar_SetValue("scr_uixscale", 
					CLAMP(0, scr_uixscale.value+0.01, 1));
		else if (display_cursor == 6) {
			s32 i = 0;
			for (; i < (s32)sizeof(fpslimits) / 4 - 1; ++i)
				if (fpslimits[i] >= (s32)host_maxfps.value)
					break;
			++i;
			if (i < 0 || i > (s32)sizeof(fpslimits) / 4 - 1)
				i = 0;
			Cvar_SetValue("host_maxfps", fpslimits[i]);
		} else if (display_cursor == 7) Cvar_SetValue("scr_showfps",
				scr_showfps.value ? 0 : 1);
		else if (display_cursor == 8)
			Cvar_SetValue("fov", scr_fov.value + 1);
		else if (display_cursor == 9)
			Cvar_SetValue("yaspectscale", yaspectscale.value+0.01);
		else if (display_cursor == 10) {
			if (newwinmode == 2) newwinmode = 0;
			else newwinmode++;
		} else if (display_cursor == 11 && newwinmode == 1) {
			if (display_sel_moden == 0) display_sel_moden = moden-1;
			else display_sel_moden--;
		} else if (display_cursor == 12 && newwinmode == 1) {
			VID_SetMode(0,mode->w, mode->h, newwinmode, vid_curpal);
		} else if (display_cursor == 13
			   && Q_atoi(customwidthstr) >= 320
			   && Q_atoi(customheightstr) >= 200
			   && Q_atoi(customwidthstr) <= MAXWIDTH
			   && Q_atoi(customheightstr) <= MAXHEIGHT)
			VID_SetMode(0, Q_atoi(customwidthstr),
				    Q_atoi(customheightstr), newwinmode,
				    vid_curpal);
		break;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		if (display_cursor == 0 && newwinmode == 1) display_cursor = 12;
		else if(display_cursor==0 && newwinmode!=1) display_cursor = 13;
		else display_cursor--;
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		if (display_cursor == 12 && newwinmode == 1) display_cursor = 0;
		else if(display_cursor==13 && newwinmode!=1) display_cursor = 0;
		else display_cursor++;
		break;
	case K_BACKSPACE:
		if (display_cursor == 11 && strlen(customwidthstr))
			customwidthstr[strlen(customwidthstr) - 1] = 0;
		else if (display_cursor == 12 && strlen(customheightstr))
			customheightstr[strlen(customheightstr) - 1] = 0;
		break;
	default:
		if (k < '0' || k > '9') break;
		if (display_cursor == 11) {
			s32 l = strlen(customwidthstr);
			if (l < 7) {
				customwidthstr[l + 1] = 0;
				customwidthstr[l] = k;
			}
		} else if (display_cursor == 12) {
			s32 l = strlen(customheightstr);
			if (l < 7) {
				customheightstr[l + 1] = 0;
				customheightstr[l] = k;
			}
		}
	}
}

void M_Menu_Display_f()
{
	key_dest = key_menu;
	m_state = m_display;
	m_entersound = 1;
}

void M_Display_Draw()
{
	s8 temp[32];
	s32 xoffset = 0;
	if (newwinmode != 1) {
		if (display_cursor < 11) M_DrawCursor(192, 32+display_cursor*8);
		else if (display_cursor == 11) M_DrawCursor(192, 124);
		else if (display_cursor == 12) M_DrawCursor(192, 140);
		else if (display_cursor == 13) M_DrawCursor(192, 152);
	}
	else M_DrawCursor(192, 32 + display_cursor*8);
	M_DrawTransPic(16, 4, Draw_CachePic("gfx/qplaque.lmp"));
	qpic_t *p = Draw_CachePic("gfx/p_option.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	M_Print(xoffset, 32, "              UI Scale");
	M_Print(xoffset, 40, "         Lock UI Scale");
	M_Print(xoffset, 48, "             HUD Style");
	M_Print(xoffset, 56, "             Crosshair");
	M_Print(xoffset, 64, "          Aspect Ratio");
	M_Print(xoffset, 72, "        UI Width Ratio");
	M_Print(xoffset, 80, "             FPS Limit");
	M_Print(xoffset, 88, "           FPS Counter");
	M_Print(xoffset, 96, "         Field Of View");
	M_Print(xoffset, 104, "          Vertical FOV");
	M_Print(xoffset, 112, "           Window Mode");
	sprintf(temp, "x%d\n", (s32)scr_uiscale.value);
	M_Print(xoffset + 204, 32, temp);
	M_Print(xoffset + 204, 40, (s32)scr_lockuiscale.value==0 ? "no":"yes");
	if     (scr_hudstyle.value==0)M_Print(xoffset+204,48,"Classic");
	else if(scr_hudstyle.value==1)M_Print(xoffset+204,48,"Modern1");
	else if(scr_hudstyle.value==2)M_Print(xoffset+204,48,"Modern2");
	else if(scr_hudstyle.value==3)M_Print(xoffset+204,48,"QW");
	else if(scr_hudstyle.value==4)M_Print(xoffset+204,48,"Arcade");
	else if(scr_hudstyle.value==5)M_Print(xoffset+204,48,"Minimalist1");
	else if(scr_hudstyle.value==6)M_Print(xoffset+204,48,"Minimalist2");
	else if(scr_hudstyle.value==7)M_Print(xoffset+204,48,"Minimalist3");
	else if(scr_hudstyle.value==8)M_Print(xoffset+204,48,"Classic no BG");
	else                          M_Print(xoffset+204,48,"EZQuake");
	if (crosshair.value) { temp[0] = cl_crosschar.value+128; temp[1] = 0; }
	else sprintf(temp, "Off");
	M_Print(xoffset + 204, 56, temp);
	sprintf(temp, "%0.2f\n", aspectr.value);
	M_Print(xoffset + 204, 64, temp);
	sprintf(temp, "%0.2f\n", scr_uixscale.value);
	M_Print(xoffset + 204, 72, temp);
	sprintf(temp, "%d\n", (s32)host_maxfps.value);
	M_Print(xoffset + 204, 80, temp);
	M_Print(xoffset + 204, 88, (s32)scr_showfps.value==0 ? "off":"on");
	sprintf(temp, "%0.2f\n", aspectr.value);
	sprintf(temp, "%d\n", (s32)scr_fov.value);
	M_Print(xoffset + 204, 96, temp);
	sprintf(temp, "%0.2f\n", yaspectscale.value);
	M_Print(xoffset + 204, 104, temp);
	if (newwinmode == 0)      M_Print(xoffset + 204, 112, "Windowed");
	else if (newwinmode == 1) M_Print(xoffset + 204, 112, "Fullscreen");
	else                      M_Print(xoffset + 204, 112, "Borderless");
	if (display_cursor == 4) {
		M_DrawTextBox(52, 166, 30, 1);
		M_Print(64, 174, "          in console for auto");
		M_PrintWhite(64, 174, "aspectr 0");
	} else if (display_cursor == 6 && (s32)host_maxfps.value > 72) {
		M_DrawTextBox(52, 166, 30, 1);
		M_Print(64, 174, "Vanilla max is    expect bugs");
		M_PrintWhite(64, 174, "               72");
	} else if (display_cursor == 8 || display_cursor == 9) {
		M_DrawTextBox(100, 166, 16, 1);
		M_Print(112, 174, "Modern FOV:");
		s32 modernfov = atan(vid.height*yaspectscale.value / (vid.width
			/ tan(scr_fov.value / 360 * M_PI))) * 360 / M_PI;
		q_snprintf(temp, 32, "            %d", modernfov);
		M_PrintWhite(112, 174, temp);
	}
	if (newwinmode == 1) {
		if (!modes) { 
			modes = SDL_GetFullscreenDisplayModes(1, &moden);
			mode = modes[0];
		}
		if (moden <= 0) return;
		mode = modes[display_sel_moden];
		M_Print(xoffset, 120, "                  Mode");
		q_snprintf(temp, 16, "%dx%d@%d", mode->w, mode->h,
				(s32)(mode->refresh_rate+0.5));
		M_Print(xoffset + 204, 120, temp);
	} else {
		M_Print(xoffset, 124, "          Custom Width");
		M_DrawTextBox(xoffset + 196, 116, 8, 1);
		M_Print(xoffset + 204, 124, customwidthstr);
		if (display_cursor == 11) {
			M_DrawCursorLine(xoffset+204 + 8 * strlen(customwidthstr), 124);
			sprintf(temp, "%d", MAXWIDTH);
			M_DrawTextBox(xoffset + 68, 166, 16 + strlen(temp), 1);
			M_Print(xoffset + 80, 174, "320 <= Width <=");
			M_Print(xoffset + 208, 174, temp);
		}
		M_Print(xoffset, 140, "         Custom Height");
		M_DrawTextBox(xoffset + 196, 132, 8, 1);
		M_Print(xoffset + 204, 140, customheightstr);
		if (display_cursor == 12) {
			M_DrawCursorLine(xoffset+204+8 * strlen(customheightstr), 140);
			sprintf(temp, "%d", MAXHEIGHT);
			M_DrawTextBox(xoffset + 68, 166, 17 + strlen(temp), 1);
			M_Print(xoffset + 80, 174, "200 <= Height <=");
			M_Print(xoffset + 216, 174, temp);
		}
	}
	M_Print(xoffset + 204, newwinmode==1?128:152, "Set Mode");
}

void M_Graphics_Key(s32 k)
{
	switch (k) {
	case K_ESCAPE:
		if (graphics_cursor >= 100) { graphics_cursor /= 100; break; }
		M_Menu_New_f();
		break;
	case K_LEFTARROW:
		S_LocalSound("misc/menu3.wav");
		switch (graphics_cursor) {
		case 200: Cvar_SetValue("r_nofog",
			!r_nofog.value); break;
		case 201: Cvar_SetValue("r_fognoise",
			r_fognoise.value - 0.1); break;
		case 202: Cvar_SetValue("r_fogfactor",
			r_fogfactor.value - 0.1); break;
		case 203: Cvar_SetValue("r_fogscale",
			r_fogscale.value - 0.1); break;
		case 204: Cvar_SetValue("r_fogbrightness",
			r_fogbrightness.value - 0.1); break;
		case 205: Cvar_SetValue("r_fogstyle",
			CLAMP(0, r_fogstyle.value - 1, 3)); break;
		case 206: Cvar_SetValue("r_fogdepthcorrection",
			!r_fogdepthcorrection.value); break;
		case 207: Cvar_SetValue("r_lockfog",
			!r_lockfog.value); break;
		case 208: Cvar_SetValue("r_lockfogd",
			CLAMP(0, r_lockfogd.value - 0.05, 1)); break;
		case 209: Cvar_SetValue("r_lockfogr",
			CLAMP(0, r_lockfogr.value - 0.05, 1)); break;
		case 210: Cvar_SetValue("r_lockfogg",
			CLAMP(0, r_lockfogg.value - 0.05, 1)); break;
		case 211: Cvar_SetValue("r_lockfogb",
			CLAMP(0, r_lockfogb.value - 0.05, 1)); break;
		case 300: Cvar_SetValue("r_enableskybox",
			!r_enableskybox.value); break;
		case 301: Cvar_SetValue("r_skyfog",
			r_skyfog.value - 0.1); break;
		case 400: Cvar_SetValue("r_rgblighting",
			!r_rgblighting.value); break;
		case 401: Cvar_SetValue("r_litwater",
			!r_litwater.value); break;
		case 402: Cvar_SetValue("r_labmixpal",
			!r_labmixpal.value); break;
		case 500: Cvar_SetValue("r_twopass",
			r_twopass.value <= 0 ? 3 : r_twopass.value - 1);
			  Cvar_SetValue("r_entalpha",
			((s32)r_twopass.value)&1); break;
		case 501: Cvar_SetValue("r_alphastyle",
			!r_alphastyle.value); break;
		case 502: Cvar_SetValue("r_wateralpha",
			CLAMP(0, r_wateralpha.value - 0.1, 1)); break;
		case 503: Cvar_SetValue("r_slimealpha",
			CLAMP(0, r_slimealpha.value - 0.1, 1)); break;
		case 504: Cvar_SetValue("r_lavaalpha",
			CLAMP(0, r_lavaalpha.value - 0.1, 1)); break;
		case 505: Cvar_SetValue("r_telealpha",
			CLAMP(0, r_telealpha.value - 0.1, 1)); break;
		case 601: Cvar_SetValue("lyr_sbar",
			CLAMP(0, lyr_sbar.value - 1, 3)); break;
		case 602: Cvar_SetValue("lyr_menu",
			CLAMP(0, lyr_menu.value - 1, 3)); break;
		case 603: Cvar_SetValue("lyr_centerprint",
			CLAMP(0, lyr_centerprint.value - 1, 3)); break;
		case 604: Cvar_SetValue("lyr_console",
			CLAMP(0, lyr_console.value - 1, 3)); break;
		case 605: Cvar_SetValue("lyr_notify",
			CLAMP(0, lyr_notify.value - 1, 3)); break;
		case 606: Cvar_SetValue("lyr_crosshair",
			CLAMP(0, lyr_crosshair.value - 1, 3)); break;
		case 700: Cvar_SetValue("r_particlescale",
			CLAMP(0, r_particlescale.value - 1, 9)); break;
		case 701: Cvar_SetValue("r_mipscale",
			CLAMP(0, r_mipscale.value - 0.1, 9.9)); break;
		case 702: Cvar_SetValue("scr_menubgstyle",
			CLAMP(0, scr_menubgstyle.value - 1, 3)); break;
		case 703: Cvar_SetValue("r_rebuildmips",
			!r_rebuildmips.value); break;
		}
		break;
	case K_RIGHTARROW:
	case K_ENTER:
		if (graphics_cursor == 0) {
			Cvar_SetValue("r_twopass", 3);
			Cvar_SetValue("r_nofog", 0);
			Cvar_SetValue("r_entalpha", 1);
			Cvar_SetValue("r_litwater", 1);
			Cvar_SetValue("r_rgblighting", 1);
			Cvar_SetValue("r_enableskybox", 1);
			Cvar_SetValue("lyr_main", 0);
			Cvar_SetValue("lyr_sbar", 3);
			Cvar_SetValue("lyr_menu", 2);
			Cvar_SetValue("lyr_centerprint", 1);
			Cvar_SetValue("lyr_console", 2);
			Cvar_SetValue("lyr_notify", 2);
			Cvar_SetValue("lyr_crosshair", 3);
			Cvar_SetValue("r_rebuildmips", 1);
			Cvar_SetValue("r_mipscale", 3);
		} else if (graphics_cursor == 1) {
			Cvar_SetValue("r_twopass", 2);
			Cvar_SetValue("r_nofog", 1);
			Cvar_SetValue("r_entalpha", 0);
			Cvar_SetValue("r_litwater", 0);
			Cvar_SetValue("r_rgblighting", 0);
			Cvar_SetValue("r_enableskybox", 0);
			Cvar_SetValue("lyr_main", 0);
			Cvar_SetValue("lyr_sbar", 0);
			Cvar_SetValue("lyr_menu", 2);
			Cvar_SetValue("lyr_centerprint", 0);
			Cvar_SetValue("lyr_console", 0);
			Cvar_SetValue("lyr_notify", 0);
			Cvar_SetValue("lyr_crosshair", 0);
			Cvar_SetValue("r_rebuildmips", 0);
			Cvar_SetValue("r_mipscale", 1);
		} else if (graphics_cursor < 100) graphics_cursor *= 100;
		else switch (graphics_cursor) {
		case 200: Cvar_SetValue("r_nofog",
			!r_nofog.value); break;
		case 201: Cvar_SetValue("r_fognoise",
			r_fognoise.value + 0.1); break;
		case 202: Cvar_SetValue("r_fogfactor",
			r_fogfactor.value + 0.1); break;
		case 203: Cvar_SetValue("r_fogscale",
			r_fogscale.value + 0.1); break;
		case 204: Cvar_SetValue("r_fogbrightness",
			r_fogbrightness.value + 0.1); break;
		case 205: Cvar_SetValue("r_fogstyle",
			CLAMP(0, r_fogstyle.value + 1, 3)); break;
		case 206: Cvar_SetValue("r_fogdepthcorrection",
			!r_fogdepthcorrection.value); break;
		case 207: Cvar_SetValue("r_lockfog",
			!r_lockfog.value); break;
		case 208: Cvar_SetValue("r_lockfogd",
			CLAMP(0, r_lockfogd.value + 0.05, 1)); break;
		case 209: Cvar_SetValue("r_lockfogr",
			CLAMP(0, r_lockfogr.value + 0.05, 1)); break;
		case 210: Cvar_SetValue("r_lockfogg",
			CLAMP(0, r_lockfogg.value + 0.05, 1)); break;
		case 211: Cvar_SetValue("r_lockfogb",
			CLAMP(0, r_lockfogb.value + 0.05, 1)); break;
		case 300: Cvar_SetValue("r_enableskybox",
			!r_enableskybox.value); break;
		case 301: Cvar_SetValue("r_skyfog",
			r_skyfog.value + 0.1); break;
		case 400: Cvar_SetValue("r_rgblighting",
			!r_rgblighting.value); break;
		case 401: Cvar_SetValue("r_litwater",
			!r_litwater.value); break;
		case 402: Cvar_SetValue("r_labmixpal",
			!r_labmixpal.value); break;
		case 500: Cvar_SetValue("r_twopass",
			r_twopass.value > 2 ? 0 : r_twopass.value + 1);
			  Cvar_SetValue("r_entalpha",
			((s32)r_twopass.value)&1); break;
		case 501: Cvar_SetValue("r_alphastyle",
			!r_alphastyle.value); break;
		case 502: Cvar_SetValue("r_wateralpha",
			CLAMP(0, r_wateralpha.value + 0.1, 1)); break;
		case 503: Cvar_SetValue("r_slimealpha",
			CLAMP(0, r_slimealpha.value + 0.1, 1)); break;
		case 504: Cvar_SetValue("r_lavaalpha",
			CLAMP(0, r_lavaalpha.value + 0.1, 1)); break;
		case 505: Cvar_SetValue("r_telealpha",
			CLAMP(0, r_telealpha.value + 0.1, 1)); break;
		case 601: Cvar_SetValue("lyr_sbar",
			CLAMP(0, lyr_sbar.value + 1, 3)); break;
		case 602: Cvar_SetValue("lyr_menu",
			CLAMP(0, lyr_menu.value + 1, 3)); break;
		case 603: Cvar_SetValue("lyr_centerprint",
			CLAMP(0, lyr_centerprint.value + 1, 3)); break;
		case 604: Cvar_SetValue("lyr_console",
			CLAMP(0, lyr_console.value + 1, 3)); break;
		case 605: Cvar_SetValue("lyr_notify",
			CLAMP(0, lyr_notify.value + 1, 3)); break;
		case 606: Cvar_SetValue("lyr_crosshair",
			CLAMP(0, lyr_crosshair.value + 1, 3)); break;
		case 607:
			Cvar_SetValue("lyr_main", 0);
			Cvar_SetValue("lyr_sbar", 3);
			Cvar_SetValue("lyr_menu", 2);
			Cvar_SetValue("lyr_centerprint", 1);
			Cvar_SetValue("lyr_console", 2);
			Cvar_SetValue("lyr_notify", 2);
			Cvar_SetValue("lyr_crosshair", 2);
			break;
		case 700: Cvar_SetValue("r_particlescale",
			CLAMP(0, r_particlescale.value + 1, 9)); break;
		case 701: Cvar_SetValue("r_mipscale",
			CLAMP(0, r_mipscale.value + 0.1, 9.9)); break;
		case 702: Cvar_SetValue("scr_menubgstyle",
			CLAMP(0, scr_menubgstyle.value + 1, 3)); break;
		case 703: Cvar_SetValue("r_rebuildmips",
			!r_rebuildmips.value); break;
		}
		S_LocalSound("misc/menu3.wav");
		break;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		if (graphics_cursor == 0) graphics_cursor = 7;
		else if (graphics_cursor == 200) {
			if (r_lockfog.value) graphics_cursor = 211;
			else graphics_cursor = 207;
		}else if(graphics_cursor == 300) graphics_cursor = 301;
		else if (graphics_cursor == 400) graphics_cursor = 402;
		else if (graphics_cursor == 500) graphics_cursor = 505;
		else if (graphics_cursor == 600) graphics_cursor = 607;
		else if (graphics_cursor == 700) graphics_cursor = 703;
		else graphics_cursor--;
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		if (graphics_cursor < 100) {
			if (graphics_cursor == 7) graphics_cursor = 0;
			else graphics_cursor++;
		} else if (graphics_cursor < 300) {
			if (r_lockfog.value) {
				if (graphics_cursor==211) graphics_cursor = 200;
				else graphics_cursor++;
			} else {
				if (graphics_cursor==207) graphics_cursor = 200;
				else graphics_cursor++;
			}
		} else if (graphics_cursor < 400) {
			if (graphics_cursor == 301) graphics_cursor = 300;
			else graphics_cursor++;
		} else if (graphics_cursor < 500) {
			if (graphics_cursor == 402) graphics_cursor = 400;
			else graphics_cursor++;
		} else if (graphics_cursor < 600) {
			if (graphics_cursor == 505) graphics_cursor = 500;
			else graphics_cursor++;
		} else if (graphics_cursor < 700) {
			if (graphics_cursor == 607) graphics_cursor = 600;
			else graphics_cursor++;
		} else if (graphics_cursor < 800) {
			if (graphics_cursor == 703) graphics_cursor = 700;
			else graphics_cursor++;
		}
		break;
	default: break;
	}
	if (graphics_cursor >= 207 && graphics_cursor <= 211)
		Fog_Update(0, 0, 0, 0);
}

void M_Menu_Graphics_f()
{
	key_dest = key_menu;
	m_state = m_graphics;
	m_entersound = 1;
}

void M_Graphics_Draw()
{
	s8 temp[32];
	s32 xoffset = 0;
	switch(graphics_cursor/100){
		case 0: M_DrawCursor(6, 32 + graphics_cursor*8); break;
		default: M_DrawCursor(6+144, 32 + (graphics_cursor%100)*8); break;
	}
	qpic_t *p = Draw_CachePic("gfx/p_option.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	M_Print(xoffset, 32, "  Preset: Modern");
	M_Print(xoffset, 40, "  Preset: Classic");
	M_Print(xoffset, 48, "  Fog...");
	M_Print(xoffset, 56, "  Sky...");
	M_Print(xoffset, 64, "  Lighting...");
	M_Print(xoffset, 72, "  Translusency...");
	M_Print(xoffset, 80, "  Layers...");
	M_Print(xoffset, 88, "  Misc...");
	xoffset = 160;
	s32 x2 = 104;
	if (graphics_cursor == 0) {
		M_Print(xoffset, 32, "Enables fog,");
		M_Print(xoffset, 40, "custom skyboxes,");
		M_Print(xoffset, 48, "colored lighting,");
		M_Print(xoffset, 56, "translusency,");
		M_Print(xoffset, 64, "layered rendering,");
		M_Print(xoffset, 72, "far+rebuilt mipmaps");
	} else if (graphics_cursor == 1) {
		M_Print(xoffset, 32, "Disables fog,");
		M_Print(xoffset, 40, "custom skyboxes,");
		M_Print(xoffset, 48, "colored lighting,");
		M_Print(xoffset, 56, "translusency,");
		M_Print(xoffset, 64, "layered rendering,");
		M_Print(xoffset, 72, "far+rebuilt mipmaps");
	} else if (graphics_cursor == 2 || graphics_cursor/100 == 2) {
		M_Print(xoffset, 32, "Enabled:");
		M_Print(xoffset + x2, 32, r_nofog.value == 0 ? "On" : "Off");
		M_Print(xoffset, 40, "Noise:");
		snprintf(temp, sizeof(temp), "%0.1f", r_fognoise.value);
		M_Print(xoffset + x2, 40, temp);
		M_Print(xoffset, 48, "Thickness:");
		snprintf(temp, sizeof(temp), "%0.1f", r_fogfactor.value);
		M_Print(xoffset + x2, 48, temp);
		M_Print(xoffset, 56, "Distance:");
		snprintf(temp, sizeof(temp), "%0.1f", r_fogscale.value);
		M_Print(xoffset + x2, 56, temp);
		M_Print(xoffset, 64, "Brightness:");
		snprintf(temp, sizeof(temp), "%0.1f", r_fogbrightness.value);
		M_Print(xoffset + x2, 64, temp);
		M_Print(xoffset, 72, "Style:");
		switch ((s32)r_fogstyle.value) {
		case 0: M_Print(xoffset + x2, 72, "Noisy"); break;
		case 1: M_Print(xoffset + x2, 72, "Dither"); break;
		case 2: M_Print(xoffset + x2, 72, "Mixed"); break;
		default:M_Print(xoffset + x2, 72, "Blend"); break;
		}
		M_Print(xoffset, 80, "Depth:");
		M_Print(xoffset + x2 - 32, 80,
			r_fogdepthcorrection.value==1 ? "Corrected" : "Simple");
		M_Print(xoffset, 88, "Set By Map:");
		M_Print(xoffset + x2, 88,
			r_lockfog.value==1 ? "Off" : "On");
		if (r_lockfog.value) {
			M_Print(xoffset, 96, "Density:");
			snprintf(temp, sizeof(temp), "%0.2f", r_lockfogd.value);
			M_Print(xoffset + x2, 96, temp);
			M_Print(xoffset, 104, "Red:");
			snprintf(temp, sizeof(temp), "%0.2f", r_lockfogr.value);
			M_Print(xoffset + x2, 104, temp);
			M_Print(xoffset, 112, "Green:");
			snprintf(temp, sizeof(temp), "%0.2f", r_lockfogg.value);
			M_Print(xoffset + x2, 112, temp);
			M_Print(xoffset, 120, "Blue:");
			snprintf(temp, sizeof(temp), "%0.2f", r_lockfogb.value);
			M_Print(xoffset + x2, 120, temp);
		}
	} else if (graphics_cursor == 3 || graphics_cursor/100 == 3) {
		M_Print(xoffset, 32, "Enabled:");
		M_Print(xoffset + x2, 32, r_enableskybox.value==0 ? "On":"Off");
		M_Print(xoffset, 40, "Sky Fog:");
		snprintf(temp, sizeof(temp), "%0.1f\n", r_skyfog.value);
		M_Print(xoffset + x2, 40, temp);
	} else if (graphics_cursor == 4 || graphics_cursor/100 == 4) {
		M_Print(xoffset, 32, "Color:");
		M_Print(xoffset+x2, 32, r_rgblighting.value==1 ? "On" : "Off");
		M_Print(xoffset, 40, "Lit Water:");
		M_Print(xoffset + x2, 40, r_litwater.value == 1 ? "On" : "Off");
		M_Print(xoffset, 48, "Color Space:");
		M_Print(xoffset+x2, 48, r_labmixpal.value==1 ? "LAB" : "RGB");
	} else if (graphics_cursor == 5 || graphics_cursor/100 == 5) {
		M_Print(xoffset, 32, "Enabled:");
		if (r_twopass.value == 0)M_Print(xoffset+72, 32, "Off (Auto)");
		else if(r_twopass.value==1)M_Print(xoffset+72, 32, "On (Auto)");
		else if(r_twopass.value==2)M_Print(xoffset+72, 32, "Off");
		else M_Print(xoffset + 72, 32, "On");
		M_Print(xoffset, 40, "Style:");
		M_Print(xoffset+x2, 40, r_alphastyle.value==1?"Dither":"Blend");
		M_Print(xoffset, 48, "Water Alpha:");
		snprintf(temp, sizeof(temp), "%0.1f\n", r_wateralpha.value);
		M_Print(xoffset + x2, 48, temp);
		M_Print(xoffset, 56, "Slime Alpha:");
		snprintf(temp, sizeof(temp), "%0.1f\n", r_slimealpha.value);
		M_Print(xoffset + x2, 56, temp);
		M_Print(xoffset, 64, "Lava Alpha:");
		snprintf(temp, sizeof(temp), "%0.1f\n", r_lavaalpha.value);
		M_Print(xoffset + x2, 64, temp);
		M_Print(xoffset, 72, "Tele Alpha:");
		snprintf(temp, sizeof(temp), "%0.1f\n", r_telealpha.value);
		M_Print(xoffset + x2, 72, temp);
	} else if (graphics_cursor == 6 || graphics_cursor/100 == 6) {
		x2 += 16;
		M_Print(xoffset, 32, "World:");
		snprintf(temp, sizeof(temp), "%d", (s32)lyr_main.value);
		M_Print(xoffset + x2, 32, temp);
		M_Print(xoffset, 40, "Status Bar:");
		snprintf(temp, sizeof(temp), "%d", (s32)lyr_sbar.value);
		M_Print(xoffset + x2, 40, temp);
		M_Print(xoffset, 48, "Menu:");
		snprintf(temp, sizeof(temp), "%d", (s32)lyr_menu.value);
		M_Print(xoffset + x2, 48, temp);
		M_Print(xoffset, 56, "Centerprints:");
		snprintf(temp, sizeof(temp), "%d", (s32)lyr_centerprint.value);
		M_Print(xoffset + x2, 56, temp);
		M_Print(xoffset, 64, "Console:");
		snprintf(temp, sizeof(temp), "%d", (s32)lyr_console.value);
		M_Print(xoffset + x2, 64, temp);
		M_Print(xoffset, 72, "Notifications:");
		snprintf(temp, sizeof(temp), "%d", (s32)lyr_notify.value);
		M_Print(xoffset + x2, 72, temp);
		M_Print(xoffset, 80, "Crosshair:");
		snprintf(temp, sizeof(temp), "%d", (s32)lyr_crosshair.value);
		M_Print(xoffset + x2, 80, temp);
		M_Print(xoffset, 88, "Reset To Default");
		M_DrawTextBox(4, 150, 36, 4);
		M_Print(16, 158, "0: custom world palette, not scaled");
		M_Print(16, 166, "1: default palette, not scaled");
		M_Print(16, 174, "2: custom ui palette, scaled");
		M_Print(16, 182, "3: custom ui palette, scaled");
	} else if (graphics_cursor == 7 || graphics_cursor/100 == 7) {
		x2 += 24;
		M_Print(xoffset, 32, "Particle Scale:");
		snprintf(temp, sizeof(temp), "x%d", (s32)r_particlescale.value);
		M_Print(xoffset + x2, 32, r_particlescale.value ? temp:"Auto");
		M_Print(xoffset, 40, "Mipmap Distance:");
		snprintf(temp, sizeof(temp), " %0.1f", r_mipscale.value);
		M_Print(xoffset + x2, 40, temp);
		M_Print(xoffset, 48, "Menu BG Style:");
		M_Print(xoffset, 56, "Rebuild Mipmaps:");
		M_Print(xoffset + x2, 56, !r_rebuildmips.value?" Off":" On");
		if (graphics_cursor == 703) {
			M_DrawTextBox(12, 150, 33, 3);
			M_Print(16, 158, "  Enable if you see pink around");
			M_Print(16, 166, " cutouts or fullbright reds that");
			M_Print(16, 174, " shouldn't be there at a distance");
		}
		switch ((s32)scr_menubgstyle.value) {
		default: case 0: M_Print(xoffset + x2, 48, "Win"); break;
		case 1: M_Print(xoffset + x2, 48, "DOS"); break;
		case 2: M_Print(xoffset + x2, 48, "Dark"); break;
		case 3: M_Print(xoffset + x2, 48, "None"); break;
		}
	}
}

void M_Menu_New_f()
{
	key_dest = key_menu;
	m_state = m_new;
	m_entersound = 1;
	sprintf(customwidthstr, "%d", vid.width);
	sprintf(customheightstr, "%d", vid.height);
}

void M_New_Draw()
{
	s8 temp[32];
	s32 xoffset = 0;
	M_DrawTransPic(16, 4, Draw_CachePic("gfx/qplaque.lmp"));
	qpic_t *p = Draw_CachePic("gfx/p_option.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	M_Print(xoffset, 32, "             This Menu");
	M_DrawCheckbox(xoffset + 204, 32, newoptions.value);
	if (new_cursor == 0) {
		M_DrawTextBox(52, 158, 30, 2);
		M_Print(64, 166, "  This menu can be restored");
		M_Print(64, 174, "with the              command");
		M_PrintWhite(64, 174, "         newoptions 1");
	}
	M_Print(xoffset, 40, "         Y Mouse Speed");
	sprintf(temp, "x%0.1f\n", sensitivityyscale.value);
	M_Print(xoffset + 204, 40, temp);
	M_Print(xoffset + 204, 48, "Display...");
	M_Print(xoffset + 204, 56, "Graphics...");
	M_Print(xoffset + 204, 64, "Gamepad...");
	M_DrawCursor(xoffset + 192, 32 + new_cursor * 8);
}

void M_Gamepad_Key(s32 k)
{
	s32 axis = SDL_GetNumJoystickAxes(joystick) - 1;
	s32 buts = SDL_GetNumJoystickButtons(joystick) - 1;
	switch (k) {
	case K_ESCAPE:
		M_Menu_New_f();
		break;
	case K_LEFTARROW:
		S_LocalSound("misc/menu3.wav");
		if (gamepad_cursor == 0)
			Cvar_SetValue("joysticknum", joysticknum.value - 1);
		else if (gamepad_cursor == 1) Cvar_SetValue("jdeadzone",
			CLAMP(8, jdeadzone.value-(1<<14)/10, (1<<14)));
		else if (gamepad_cursor == 2) Cvar_SetValue("jmovesens",
			CLAMP(0.5, jmovesens.value-0.5, 5));
		else if (gamepad_cursor == 3) Cvar_SetValue("jlooksens",
			CLAMP(0.5, jlooksens.value-0.5, 5));
		else if (gamepad_cursor == 4) Cvar_SetValue("jtriggerthresh",
		    CLAMP(8-(1<<15),jtriggerthresh.value-(1<<16)/10,(1<<15)-8));
		else if (gamepad_cursor == 5 && jmoveaxisx.value >= 0)
			Cvar_SetValue("jmoveaxisx", jmoveaxisx.value - 1);
		else if (gamepad_cursor == 6 && jmoveaxisx.value >= 0)
			Cvar_SetValue("jmoveaxisy", jmoveaxisy.value - 1);
		else if (gamepad_cursor == 7 && jlookaxisx.value >= 0)
			Cvar_SetValue("jlookaxisx", jlookaxisx.value - 1);
		else if (gamepad_cursor == 8 && jlookaxisy.value >= 0)
			Cvar_SetValue("jlookaxisy", jlookaxisy.value - 1);
		else if (gamepad_cursor == 9 && jtrigaxis1.value >= 0)
			Cvar_SetValue("trigaxis1", jtrigaxis1.value - 1);
		else if (gamepad_cursor == 10 && jtrigaxis2.value >= 0)
			Cvar_SetValue("trigaxis2", jtrigaxis2.value - 1);
		else if (gamepad_cursor == 11 && jenterbutton.value >= 0)
			Cvar_SetValue("jenterbutton", jenterbutton.value - 1);
		else if (gamepad_cursor == 12 && jescapebutton.value >= 0)
			Cvar_SetValue("jescapebutton", jescapebutton.value - 1);
		break;
	case K_RIGHTARROW:
	case K_ENTER:
		S_LocalSound("misc/menu3.wav");
		if (gamepad_cursor == 0)
			Cvar_SetValue("joysticknum", joysticknum.value + 1);
		else if (gamepad_cursor == 1) Cvar_SetValue("jdeadzone",
			CLAMP(8, jdeadzone.value+(1<<14)/10, (1<<14)));
		else if (gamepad_cursor == 2) Cvar_SetValue("jmovesens",
			CLAMP(0.5, jmovesens.value+0.5, 5));
		else if (gamepad_cursor == 3) Cvar_SetValue("jlooksens",
			CLAMP(0.5, jlooksens.value+0.5, 5));
		else if (gamepad_cursor == 4) Cvar_SetValue("jtriggerthresh",
		    CLAMP(8-(1<<15),jtriggerthresh.value+(1<<16)/10,(1<<15)-8));
		else if (gamepad_cursor == 5 && jmoveaxisx.value < axis)
			Cvar_SetValue("jmoveaxisx", jmoveaxisx.value + 1);
		else if (gamepad_cursor == 6 && jmoveaxisx.value < axis)
			Cvar_SetValue("jmoveaxisy", jmoveaxisy.value + 1);
		else if (gamepad_cursor == 7 && jlookaxisx.value < axis)
			Cvar_SetValue("jlookaxisx", jlookaxisx.value + 1);
		else if (gamepad_cursor == 8 && jlookaxisy.value < axis)
			Cvar_SetValue("jlookaxisy", jlookaxisy.value + 1);
		else if (gamepad_cursor == 9 && jtrigaxis1.value < axis)
			Cvar_SetValue("trigaxis1", jtrigaxis1.value + 1);
		else if (gamepad_cursor == 10 && jtrigaxis2.value < axis)
			Cvar_SetValue("trigaxis2", jtrigaxis2.value + 1);
		else if (gamepad_cursor == 11 && jenterbutton.value < buts)
			Cvar_SetValue("jenterbutton", jenterbutton.value + 1);
		else if (gamepad_cursor == 12 && jescapebutton.value < buts)
			Cvar_SetValue("jescapebutton", jescapebutton.value + 1);
		break;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		if (gamepad_cursor == 0)
			gamepad_cursor = 12;
		else
			gamepad_cursor--;
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		if (gamepad_cursor == 12)
			gamepad_cursor = 0;
		else
			gamepad_cursor++;
		break;
	}
}

void M_New_Key(s32 k)
{
	switch (k) {
	case K_ESCAPE:
		if (!newoptions.value) options_cursor = 0;
		M_Menu_Options_f();
		break;
	case K_LEFTARROW:
		S_LocalSound("misc/menu3.wav");
		if (new_cursor == 0)
			Cvar_SetValue("newoptions", !newoptions.value);
		else if (new_cursor == 1 && sensitivityyscale.value >= 0.1)
			Cvar_SetValue("sensitivityyscale",
				      sensitivityyscale.value - 0.1);
		break;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		if (new_cursor == 0) new_cursor = 4;
		else new_cursor--;
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		if (new_cursor == 4) new_cursor = 0;
		else new_cursor++;
		break;
	case K_RIGHTARROW:
	case K_ENTER:
		S_LocalSound("misc/menu3.wav");
		if (new_cursor == 0)
			Cvar_SetValue("newoptions", !newoptions.value);
		else if (new_cursor == 1 && sensitivityyscale.value >= 0.1)
			Cvar_SetValue("sensitivityyscale",
				      sensitivityyscale.value + 0.1);
		else if (new_cursor == 2) M_Menu_Display_f();
		else if (new_cursor == 3) M_Menu_Graphics_f();
		else if (new_cursor == 4) M_Menu_Gamepad_f();
		break;
	}
}

void M_Menu_Video_f()
{
	key_dest = key_menu;
	m_state = m_video;
	m_entersound = 1;
}

void M_Video_Draw()
{
	// CyanBun96: This whole menu isn't about real resolutions anyway since
	// all of them get scaled to the window size in the end, so these modes
	// here are just some nice-looking classics from the original
	// taken from WINQUAKE.EXE ran through wine
	s8 temp[64];
	qpic_t *p = Draw_CachePic("gfx/vidmodes.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	for (s32 i = 0; i < NUM_OLDMODES; ++i)
		sprintf(modelist[i], "%dx%d", oldmodes[i * 2],
			oldmodes[i * 2 + 1]);
	M_Print(13 * 8, 36, "Windowed Modes");
	s32 column = 16;
	s32 row = 36 + 2 * 8;
	for (s32 i = 0; i < 3; i++) {
		if (i == vid_modenum)
			M_PrintWhite(column, row, modelist[i]);
		else
			M_Print(column, row, modelist[i]);
		column += 13 * 8;
	}
	M_Print(12 * 8, 36 + 4 * 8, "Fullscreen Modes");
	column = 16;
	row = 36 + 6 * 8;
	for (s32 i = 3; i < NUM_OLDMODES; i++) {
		if (i == vid_modenum)
			M_PrintWhite(column, row, modelist[i]);
		else
			M_Print(column, row, modelist[i]);
		column += 13 * 8;
		if (((i - 3) % VID_ROW_SIZE) == (VID_ROW_SIZE - 1)) {
			column = 16;
			row += 8;
		}
	}
	if (vid_testingmode) {
		sprintf(temp, "TESTING %s", modelist[vid_line]);
		M_Print(13 * 8, 36 + MODE_AREA_HEIGHT * 8 + 8 * 4, temp);
		M_Print(9 * 8, 36 + MODE_AREA_HEIGHT * 8 + 8 * 6,
			"Please wait 5 seconds...");
	} else {
		M_Print(9 * 8, 36 + MODE_AREA_HEIGHT * 8 + 8,
			"Press Enter to set mode");
		M_Print(6 * 8, 36 + MODE_AREA_HEIGHT * 8 + 8 * 3,
			"T to test mode for 5 seconds");
		s8 *ptr = VID_GetModeDescription(vid_modenum);
		if (vid_modenum >= 0 && vid_modenum < NUM_OLDMODES) {
			sprintf(temp, "D to set default: %s", ptr);
			M_Print(2 * 8, 36 + MODE_AREA_HEIGHT * 8 + 8 * 5, temp);
		}
		ptr = VID_GetModeDescription((s32)_vid_default_mode_win.value);
		if (ptr) {
			sprintf(temp, "Current default: %s", ptr);
			M_Print(3 * 8, 36 + MODE_AREA_HEIGHT * 8 + 8 * 6, temp);
		}
		M_Print(15 * 8, 36 + MODE_AREA_HEIGHT * 8 + 8 * 8,
			"Esc to exit");
		row = 36 + 2 * 8 + (vid_line / VID_ROW_SIZE) * 8;
		column = 8 + (vid_line % VID_ROW_SIZE) * 13 * 8;
		if (vid_line >= 3)
			row += 3 * 8;
		M_DrawCursor(column, row);
	}
}

void M_Video_Key(s32 key)
{
	if (vid_testingmode)
		return;
	switch (key) {
	case K_ESCAPE:
		S_LocalSound("misc/menu1.wav");
		M_Menu_Options_f();
		break;
	case K_LEFTARROW:
		S_LocalSound("misc/menu1.wav");
		vid_line = ((vid_line / VID_ROW_SIZE) * VID_ROW_SIZE) +
		    ((vid_line + 2) % VID_ROW_SIZE);
		if (vid_line >= NUM_OLDMODES)
			vid_line = NUM_OLDMODES - 1;
		break;
	case K_RIGHTARROW:
		S_LocalSound("misc/menu1.wav");
		vid_line = ((vid_line / VID_ROW_SIZE) * VID_ROW_SIZE) +
		    ((vid_line + 4) % VID_ROW_SIZE);
		if (vid_line >= NUM_OLDMODES)
			vid_line = (vid_line / VID_ROW_SIZE) * VID_ROW_SIZE;
		break;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		vid_line -= VID_ROW_SIZE;
		if (vid_line < 0) {
			vid_line += ((NUM_OLDMODES + (VID_ROW_SIZE - 1)) /
				     VID_ROW_SIZE) * VID_ROW_SIZE;
			while (vid_line >= NUM_OLDMODES)
				vid_line -= VID_ROW_SIZE;
		}
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		vid_line += VID_ROW_SIZE;
		if (vid_line >= NUM_OLDMODES) {
			vid_line -= ((NUM_OLDMODES + (VID_ROW_SIZE - 1)) /
				     VID_ROW_SIZE) * VID_ROW_SIZE;
			while (vid_line < 0)
				vid_line += VID_ROW_SIZE;
		}
		break;
	case K_ENTER:
		S_LocalSound("misc/menu1.wav");
		VID_SetMode(vid_line, 0, 0, 0, vid_curpal);
		Cvar_SetValue("_vid_default_mode_win", vid_modenum);
		Cvar_SetValue("vid_cheight", -1);
		Cvar_SetValue("vid_cwidth", -1);
		Cvar_SetValue("vid_cwmode", -1);
		break;
	case 'T':
	case 't':
		S_LocalSound("misc/menu1.wav");
		vid_realmode = vid_modenum;
		vid_testingmode = 1;
		vid_testendtime = realtime + 5.0;
		VID_SetMode(vid_line, 0, 0, 0, vid_curpal);
		printf("VID_LINE: %d\n", vid_line);
		break;
	default:
		break;
	}
}

void M_Menu_Help_f()
{
	key_dest = key_menu;
	m_state = m_help;
	m_entersound = 1;
	help_page = 0;
}

void M_Help_Draw()
{
	M_DrawTransPic(0, 0, Draw_CachePic(va("gfx/help%i.lmp", help_page)));
}

void M_Help_Key(s32 key)
{
	switch (key) {
	case K_ESCAPE:
		M_Menu_Main_f();
		break;

	case K_UPARROW:
	case K_RIGHTARROW:
		m_entersound = 1;
		if (++help_page >= NUM_HELP_PAGES)
			help_page = 0;
		break;

	case K_DOWNARROW:
	case K_LEFTARROW:
		m_entersound = 1;
		if (--help_page < 0)
			help_page = NUM_HELP_PAGES - 1;
		break;
	}

}

void M_Menu_Quit_f()
{
	if (m_state == m_quit)
		return;
	wasInMenus = (key_dest == key_menu);
	key_dest = key_menu;
	m_quit_prevstate = m_state;
	m_state = m_quit;
	m_entersound = 1;
	msgNumber = rand() & 7;
}

void M_Quit_Key(s32 key)
{
	switch (key) {
	case K_ESCAPE:
	case 'n':
	case 'N':
		if (wasInMenus) {
			m_state = m_quit_prevstate;
			m_entersound = 1;
		} else {
			key_dest = key_game;
			m_state = m_none;
		}
		break;

	case 'Y':
	case 'y':
	case K_ENTER:
		key_dest = key_console;
		Host_Quit_f();
		break;

	default:
		break;
	}

}

void M_Quit_Draw()
{
	if (wasInMenus) {
		m_state = m_quit_prevstate;
		m_recursiveDraw = 1;
		M_Draw();
		m_state = m_quit;
	}
	M_DrawTextBox(56, 76, 24, 4);
	M_Print(64, 84, quitMessage[msgNumber * 4 + 0]);
	M_Print(64, 92, quitMessage[msgNumber * 4 + 1]);
	M_Print(64, 100, quitMessage[msgNumber * 4 + 2]);
	M_Print(64, 108, quitMessage[msgNumber * 4 + 3]);
}

void M_Menu_LanConfig_f()
{
	key_dest = key_menu;
	m_state = m_lanconfig;
	m_entersound = 1;
	if (lanConfig_cursor == -1) {
		if (JoiningGame && TCPIPConfig)
			lanConfig_cursor = 2;
		else
			lanConfig_cursor = 1;
	}
	if (StartingGame && lanConfig_cursor == 2)
		lanConfig_cursor = 1;
	lanConfig_port = DEFAULTnet_hostport;
	sprintf(lanConfig_portname, "%u", lanConfig_port);
	m_return_onerror = 0;
	m_return_reason[0] = 0;
}

void M_LanConfig_Draw()
{
	M_DrawTransPic(16, 4, Draw_CachePic("gfx/qplaque.lmp"));
	qpic_t *p = Draw_CachePic("gfx/p_multi.lmp");
	s32 basex = (320 - p->width) / 2;
	M_DrawTransPic(basex, 4, p);
	s8 *startJoin;
	if (StartingGame)
		startJoin = "New Game";
	else
		startJoin = "Join Game";
	s8 *protocol;
	protocol = "TCP/IP";
	M_Print(basex, 32, va("%s - %s", startJoin, protocol));
	basex += 8;
	M_Print(basex, 52, "Address:");
	M_Print(basex + 9 * 8, 52, my_tcpip_address);
	M_Print(basex, lanConfig_cursor_table[0], "Port");
	M_DrawTextBox(basex + 8 * 8, lanConfig_cursor_table[0] - 8, 6, 1);
	M_Print(basex + 9 * 8, lanConfig_cursor_table[0], lanConfig_portname);
	if (JoiningGame) {
		M_Print(basex, lanConfig_cursor_table[1],
			"Search for local games...");
		M_Print(basex, 108, "Join game at:");
		M_DrawTextBox(basex + 8, lanConfig_cursor_table[2] - 8, 22, 1);
		M_Print(basex + 16, lanConfig_cursor_table[2],
			lanConfig_joinname);
	} else {
		M_DrawTextBox(basex, lanConfig_cursor_table[1] - 8, 2, 1);
		M_Print(basex + 8, lanConfig_cursor_table[1], "OK");
	}
	M_DrawCursor(basex - 8, lanConfig_cursor_table[lanConfig_cursor]);
	if (lanConfig_cursor == 0)
		M_DrawCursor(basex + 9 * 8 + 8 * strlen(lanConfig_portname),
				lanConfig_cursor_table[0]);
	if (lanConfig_cursor == 2)
		M_DrawCursor(basex + 16 + 8 * strlen(lanConfig_joinname),
				lanConfig_cursor_table[2]);
	if (*m_return_reason)
		M_PrintWhite(basex, 148, m_return_reason);
}

void M_LanConfig_Key(s32 key)
{
	switch (key) {
	case K_ESCAPE:
		M_Menu_Net_f();
		break;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		lanConfig_cursor--;
		if (lanConfig_cursor < 0)
			lanConfig_cursor = NUM_LANCONFIG_CMDS - 1;
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		lanConfig_cursor++;
		if (lanConfig_cursor >= NUM_LANCONFIG_CMDS)
			lanConfig_cursor = 0;
		break;
	case K_ENTER:
		if (lanConfig_cursor == 0)
			break;
		m_entersound = 1;
		M_ConfigureNetSubsystem();
		if (lanConfig_cursor == 1) {
			if (StartingGame) {
				M_Menu_GameOptions_f();
				break;
			}
			M_Menu_Search_f();
			break;
		}
		else if (lanConfig_cursor == 2) {
			m_return_state = m_state;
			m_return_onerror = 1;
			key_dest = key_game;
			m_state = m_none;
			Cbuf_AddText(va("connect \"%s\"\n",lanConfig_joinname));
			break;
		}
		break;
	case K_BACKSPACE:
		if (lanConfig_cursor == 0 && strlen(lanConfig_portname)) {
			lanConfig_portname[strlen(lanConfig_portname) - 1] = 0;
		}
		else if (lanConfig_cursor == 2 && strlen(lanConfig_joinname)) {
			lanConfig_joinname[strlen(lanConfig_joinname) - 1] = 0;
		}
		break;
	default:
		if (key < 32 || key > 127)
			break;
		if (lanConfig_cursor == 2) {
			s32 l = strlen(lanConfig_joinname);
			if (l < 21) {
				lanConfig_joinname[l + 1] = 0;
				lanConfig_joinname[l] = key;
			}
		}
		if (key < '0' || key > '9')
			break;
		if (lanConfig_cursor == 0) {
			s32 l = strlen(lanConfig_portname);
			if (l < 5) {
				lanConfig_portname[l + 1] = 0;
				lanConfig_portname[l] = key;
			}
		}
	}
	if (StartingGame && lanConfig_cursor == 2)
		lanConfig_cursor = (key == K_UPARROW);
	s32 l = Q_atoi(lanConfig_portname);
	if (l > 65535)
		l = lanConfig_port;
	else
		lanConfig_port = l;
	sprintf(lanConfig_portname, "%u", lanConfig_port);
}

void M_Menu_GameOptions_f()
{
	key_dest = key_menu;
	m_state = m_gameoptions;
	m_entersound = 1;
	if (maxplayers == 0)
		maxplayers = svs.maxclients;
	if (maxplayers < 2)
		maxplayers = svs.maxclientslimit;
}

void M_GameOptions_Draw()
{
	M_DrawTransPic(16, 4, Draw_CachePic("gfx/qplaque.lmp"));
	qpic_t *p = Draw_CachePic("gfx/p_multi.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	M_DrawTextBox(152, 32, 10, 1);
	M_Print(160, 40, "begin game");
	M_Print(0, 56, "      Max players");
	M_Print(160, 56, va("%i", maxplayers));
	M_Print(0, 64, "        Game Type");
	if (coop.value)
		M_Print(160, 64, "Cooperative");
	else
		M_Print(160, 64, "Deathmatch");
	M_Print(0, 72, "        Teamplay");
	s8 *msg;
	if (rogue) {
		switch ((s32)teamplay.value) {
			case 1: msg = "No Friendly Fire"; break;
			case 2: msg = "Friendly Fire"; break;
			case 3: msg = "Tag"; break;
			case 4: msg = "Capture the Flag"; break;
			case 5: msg = "One Flag CTF"; break;
			case 6: msg = "Three Team CTF"; break;
			default: msg = "Off"; break;
		}
		M_Print(160, 72, msg);
	} else {
		switch ((s32)teamplay.value) {
			case 1: msg = "No Friendly Fire"; break;
			case 2: msg = "Friendly Fire"; break;
			default: msg = "Off"; break;
		}
		M_Print(160, 72, msg);
	}
	M_Print(0, 80, "            Skill");
	if (skill.value == 0)
		M_Print(160, 80, "Easy difficulty");
	else if (skill.value == 1)
		M_Print(160, 80, "Normal difficulty");
	else if (skill.value == 2)
		M_Print(160, 80, "Hard difficulty");
	else
		M_Print(160, 80, "Nightmare difficulty");
	M_Print(0, 88, "       Frag Limit");
	if (fraglimit.value == 0)
		M_Print(160, 88, "none");
	else
		M_Print(160, 88, va("%i frags", (s32)fraglimit.value));
	M_Print(0, 96, "       Time Limit");
	if (timelimit.value == 0)
		M_Print(160, 96, "none");
	else
		M_Print(160, 96, va("%i minutes", (s32)timelimit.value));
	M_Print(0, 112, "         Episode");
	if (hipnotic) //MED 01/06/97 added hipnotic episodes
		M_Print(160, 112, hipnoticepisodes[startepisode].description);
	else if (rogue) //PGM 01/07/97 added rogue episodes
		M_Print(160, 112, rogueepisodes[startepisode].description);
	else
		M_Print(160, 112, episodes[startepisode].description);
	M_Print(0, 120, "           Level");
	if (hipnotic) { //MED 01/06/97 added hipnotic episodes
		M_Print(160, 120,
			hipnoticlevels[hipnoticepisodes[startepisode].
				       firstLevel + startlevel].description);
		M_Print(160, 128,
			hipnoticlevels[hipnoticepisodes[startepisode].
				       firstLevel + startlevel].name);
	}
	else if (rogue) { //PGM 01/07/97 added rogue episodes
		M_Print(160, 120,
			roguelevels[rogueepisodes[startepisode].firstLevel +
				    startlevel].description);
		M_Print(160, 128,
			roguelevels[rogueepisodes[startepisode].firstLevel +
				    startlevel].name);
	} else {
		M_Print(160, 120,
			levels[episodes[startepisode].firstLevel +
			       startlevel].description);
		M_Print(160, 128,
			levels[episodes[startepisode].firstLevel +
			       startlevel].name);
	}
	M_DrawCursor(144, gameoptions_cursor_table[gameoptions_cursor]);
	if (m_serverInfoMessage) {
		if ((realtime - m_serverInfoMessageTime) < 5.0) {
			s32 x = (320 - 26 * 8) / 2;
			M_DrawTextBox(x, 138, 24, 4);
			x += 8;
			M_Print(x, 146, "  More than 4 players   ");
			M_Print(x, 154, " requires using command ");
			M_Print(x, 162, "line parameters; please ");
			M_Print(x, 170, "   see techinfo.txt.    ");
		} else
			m_serverInfoMessage = 0;
	}
}

void M_NetStart_Change(s32 dir)
{
	s32 count;
	switch (gameoptions_cursor) {
	case 1:
		maxplayers += dir;
		if (maxplayers > svs.maxclientslimit) {
			maxplayers = svs.maxclientslimit;
			m_serverInfoMessage = 1;
			m_serverInfoMessageTime = realtime;
		}
		if (maxplayers < 2)
			maxplayers = 2;
		break;
	case 2:
		Cvar_SetValue("coop", coop.value ? 0 : 1);
		break;
	case 3:
		if (rogue)
			count = 6;
		else
			count = 2;
		Cvar_SetValue("teamplay", teamplay.value + dir);
		if (teamplay.value > count)
			Cvar_SetValue("teamplay", 0);
		else if (teamplay.value < 0)
			Cvar_SetValue("teamplay", count);
		break;
	case 4:
		Cvar_SetValue("skill", skill.value + dir);
		if (skill.value > 3)
			Cvar_SetValue("skill", 0);
		if (skill.value < 0)
			Cvar_SetValue("skill", 3);
		break;
	case 5:
		Cvar_SetValue("fraglimit", fraglimit.value + dir * 10);
		if (fraglimit.value > 100)
			Cvar_SetValue("fraglimit", 0);
		if (fraglimit.value < 0)
			Cvar_SetValue("fraglimit", 100);
		break;
	case 6:
		Cvar_SetValue("timelimit", timelimit.value + dir * 5);
		if (timelimit.value > 60)
			Cvar_SetValue("timelimit", 0);
		if (timelimit.value < 0)
			Cvar_SetValue("timelimit", 60);
		break;
	case 7:
		startepisode += dir;
		if (hipnotic) //MED 01/06/97 added hipnotic count
			count = 6;
		else if (rogue) //PGM 01/07/97 added rogue count
			count = 4; //PGM 03/02/97 added 1 for dmatch episode
		else if (registered.value)
			count = 7;
		else
			count = 2;
		if (startepisode < 0)
			startepisode = count - 1;
		if (startepisode >= count)
			startepisode = 0;
		startlevel = 0;
		break;
	case 8:
		startlevel += dir;
		if (hipnotic) //MED 01/06/97 added hipnotic episodes
			count = hipnoticepisodes[startepisode].levels;
		else if (rogue) //PGM 01/06/97 added hipnotic episodes
			count = rogueepisodes[startepisode].levels;
		else
			count = episodes[startepisode].levels;
		if (startlevel < 0)
			startlevel = count - 1;
		if (startlevel >= count)
			startlevel = 0;
		break;
	}
}

void M_GameOptions_Key(s32 key)
{
	switch (key) {
	case K_ESCAPE:
		M_Menu_Net_f();
		break;
	case K_UPARROW:
		S_LocalSound("misc/menu1.wav");
		gameoptions_cursor--;
		if (gameoptions_cursor < 0)
			gameoptions_cursor = NUM_GAMEOPTIONS - 1;
		break;
	case K_DOWNARROW:
		S_LocalSound("misc/menu1.wav");
		gameoptions_cursor++;
		if (gameoptions_cursor >= NUM_GAMEOPTIONS)
			gameoptions_cursor = 0;
		break;
	case K_LEFTARROW:
		if (gameoptions_cursor == 0)
			break;
		S_LocalSound("misc/menu3.wav");
		M_NetStart_Change(-1);
		break;
	case K_RIGHTARROW:
		if (gameoptions_cursor == 0)
			break;
		S_LocalSound("misc/menu3.wav");
		M_NetStart_Change(1);
		break;
	case K_ENTER:
		S_LocalSound("misc/menu2.wav");
		if (gameoptions_cursor == 0) {
			if (sv.active)
				Cbuf_AddText("disconnect\n");
			Cbuf_AddText("listen 0\n"); // so host_netport will be re-examined
			Cbuf_AddText(va("maxplayers %u\n", maxplayers));
			SCR_BeginLoadingPlaque();
			if (hipnotic)
				Cbuf_AddText(va
					     ("map %s\n",
					      hipnoticlevels[hipnoticepisodes
							     [startepisode].
							     firstLevel +
							     startlevel].name));
			else if (rogue)
				Cbuf_AddText(va
					     ("map %s\n",
					      roguelevels[rogueepisodes
							  [startepisode].
							  firstLevel +
							  startlevel].name));
			else
				Cbuf_AddText(va
					     ("map %s\n",
					      levels[episodes[startepisode].
						     firstLevel +
						     startlevel].name));

			return;
		}
		M_NetStart_Change(1);
		break;
	}
}

void M_Menu_Search_f()
{
	key_dest = key_menu;
	m_state = m_search;
	m_entersound = 0;
	slistSilent = 1;
	slistLocal = 0;
	searchComplete = 0;
	NET_Slist_f();
}

void M_Search_Draw()
{
	qpic_t *p = Draw_CachePic("gfx/p_multi.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	s32 x = (320 / 2) - ((12 * 8) / 2) + 4;
	M_DrawTextBox(x - 8, 32, 12, 1);
	M_Print(x, 40, "Searching...");
	if (slistInProgress) {
		NET_Poll();
		return;
	}
	if (!searchComplete) {
		searchComplete = 1;
		searchCompleteTime = realtime;
	}
	if (hostCacheCount) {
		M_Menu_ServerList_f();
		return;
	}
	M_PrintWhite((320 / 2) - ((22 * 8) / 2), 64, "No Quake servers found");
	if ((realtime - searchCompleteTime) < 3.0)
		return;
	M_Menu_LanConfig_f();
}

void M_Menu_ServerList_f()
{
	key_dest = key_menu;
	m_state = m_slist;
	m_entersound = 1;
	slist_cursor = 0;
	m_return_onerror = 0;
	m_return_reason[0] = 0;
	slist_sorted = 0;
}

void M_ServerList_Draw()
{
	if (!slist_sorted) {
		if (hostCacheCount > 1) {
			s32 i, j;
			hostcache_t temp;
			for (i = 0; i < hostCacheCount; i++)
				for (j = i + 1; j < hostCacheCount; j++)
					if (strcmp
					    (hostcache[j].name,
					     hostcache[i].name) < 0) {
						Q_memcpy(&temp, &hostcache[j],
							 sizeof(hostcache_t));
						Q_memcpy(&hostcache[j],
							 &hostcache[i],
							 sizeof(hostcache_t));
						Q_memcpy(&hostcache[i], &temp,
							 sizeof(hostcache_t));
					}
		}
		slist_sorted = 1;
	}
	qpic_t *p = Draw_CachePic("gfx/p_multi.lmp");
	M_DrawTransPic((320 - p->width) / 2, 4, p);
	for (s32 n = 0; n < hostCacheCount; n++) {
		s8 string[64];
		if (hostcache[n].maxusers)
			sprintf(string, "%-15.15s %-15.15s %2u/%2u\n",
				hostcache[n].name, hostcache[n].map,
				hostcache[n].users, hostcache[n].maxusers);
		else
			sprintf(string, "%-15.15s %-15.15s\n",
				hostcache[n].name, hostcache[n].map);
		M_Print(16, 32 + 8 * n, string);
	}
	M_DrawCharacter(0, 32 + slist_cursor * 8,
			12 + ((s32)(realtime * 4) & 1));
	if (*m_return_reason)
		M_PrintWhite(16, 148, m_return_reason);
}

void M_ServerList_Key(s32 k)
{
	switch (k) {
	case K_ESCAPE:
		M_Menu_LanConfig_f();
		break;
	case K_SPACE:
		M_Menu_Search_f();
		break;
	case K_UPARROW:
	case K_LEFTARROW:
		S_LocalSound("misc/menu1.wav");
		slist_cursor--;
		if (slist_cursor < 0)
			slist_cursor = hostCacheCount - 1;
		break;
	case K_DOWNARROW:
	case K_RIGHTARROW:
		S_LocalSound("misc/menu1.wav");
		slist_cursor++;
		if (slist_cursor >= hostCacheCount)
			slist_cursor = 0;
		break;
	case K_ENTER:
		S_LocalSound("misc/menu2.wav");
		m_return_state = m_state;
		m_return_onerror = 1;
		slist_sorted = 0;
		key_dest = key_game;
		m_state = m_none;
		Cbuf_AddText(va
			     ("connect \"%s\"\n",
			      hostcache[slist_cursor].cname));
		break;
	default:
		break;
	}

}

void M_Init()
{
	Cmd_AddCommand("togglemenu", M_ToggleMenu_f);
	Cmd_AddCommand("menu_main", M_Menu_Main_f);
	Cmd_AddCommand("menu_singleplayer", M_Menu_SinglePlayer_f);
	Cmd_AddCommand("menu_load", M_Menu_Load_f);
	Cmd_AddCommand("menu_save", M_Menu_Save_f);
	Cmd_AddCommand("menu_multiplayer", M_Menu_MultiPlayer_f);
	Cmd_AddCommand("menu_setup", M_Menu_Setup_f);
	Cmd_AddCommand("menu_options", M_Menu_Options_f);
	Cmd_AddCommand("menu_keys", M_Menu_Keys_f);
	Cmd_AddCommand("menu_video", M_Menu_Video_f);
	Cmd_AddCommand("menu_new", M_Menu_New_f);
	Cmd_AddCommand("menu_gamepad", M_Menu_Gamepad_f);
	Cmd_AddCommand("menu_display", M_Menu_Display_f);
	Cmd_AddCommand("menu_graphics", M_Menu_Graphics_f);
	Cmd_AddCommand("help", M_Menu_Help_f);
	Cmd_AddCommand("menu_quit", M_Menu_Quit_f);
}

void M_Draw()
{
	if (m_state == m_none || key_dest != key_menu)
		return;
	if (!m_recursiveDraw) {
		if (scr_con_current) {
			Draw_ConsoleBackground(vid.height);
		} else
			fadescreen = 1;
		scr_fullupdate = 0;
	} else {
		m_recursiveDraw = 0;
	}
	switch (m_state) {
		case m_none: break;
		case m_main: M_Main_Draw(); break;
		case m_singleplayer: M_SinglePlayer_Draw(); break;
		case m_load: M_Load_Draw(); break;
		case m_save: M_Save_Draw(); break;
		case m_multiplayer: M_MultiPlayer_Draw(); break;
		case m_setup: M_Setup_Draw(); break;
		case m_net: M_Net_Draw(); break;
		case m_options: M_Options_Draw(); break;
		case m_keys: M_Keys_Draw(); break;
		case m_new: M_New_Draw(); break;
		case m_gamepad: M_Gamepad_Draw(); break;
		case m_display: M_Display_Draw(); break;
		case m_graphics: M_Graphics_Draw(); break;
		case m_video: M_Video_Draw(); break;
		case m_help: M_Help_Draw(); break;
		case m_quit: M_Quit_Draw(); break;
		case m_lanconfig: M_LanConfig_Draw(); break;
		case m_gameoptions: M_GameOptions_Draw(); break;
		case m_search: M_Search_Draw(); break;
		case m_slist: M_ServerList_Draw(); break;
	}
	if (m_entersound) {
		S_LocalSound("misc/menu2.wav");
		m_entersound = 0;
	}
}

void M_Keydown(s32 key)
{
	switch (m_state) {
		case m_none: return;
		case m_main: M_Main_Key(key); return;
		case m_singleplayer: M_SinglePlayer_Key(key); return;
		case m_load: M_Load_Key(key); return;
		case m_save: M_Save_Key(key); return;
		case m_multiplayer: M_MultiPlayer_Key(key); return;
		case m_setup: M_Setup_Key(key); return;
		case m_net: M_Net_Key(key); return;
		case m_options: M_Options_Key(key); return;
		case m_keys: M_Keys_Key(key); return;
		case m_video: M_Video_Key(key); return;
		case m_new: M_New_Key(key); return;
		case m_gamepad: M_Gamepad_Key(key); return;
		case m_display: M_Display_Key(key); return;
		case m_graphics: M_Graphics_Key(key); return;
		case m_help: M_Help_Key(key); return;
		case m_quit: M_Quit_Key(key); return;
		case m_lanconfig: M_LanConfig_Key(key); return;
		case m_gameoptions: M_GameOptions_Key(key); return;
		case m_search: break;
		case m_slist: M_ServerList_Key(key); return;
	}
}

void M_ConfigureNetSubsystem()
{
	Cbuf_AddText("stopdemo\n");
	net_hostport = lanConfig_port;
}
