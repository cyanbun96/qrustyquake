// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"

#define SCL uiscale // readability macros
#define WW vid.width
#define HH vid.height
#define TT 0xff // transparent // pixel art macros
#define SS 0x45 // shade
#define BB 0x49 // base
#define LL 0x4d // lit
static qpic_t shield_r = { 8, 8, {
        TT,SS,SS,SS,SS,SS, 0,TT,
        SS,BB,BB,LL,BB,BB,SS, 0,
        SS,BB,LL,BB,BB,BB,SS, 0,
        SS,BB,LL,BB,BB,BB,SS, 0,
        SS,BB,BB,BB,BB,BB,SS, 0,
        TT,SS,BB,BB,BB,SS, 0,TT,
        TT,TT,SS,BB,SS, 0,TT,TT,
        TT,TT,TT,SS, 0,TT,TT,TT,
}};
#undef SS
#undef BB
#undef LL
#define SS 0xca // shade
#define BB 0xc6 // base
#define LL 0xc2 // lit
static qpic_t shield_y = { 8, 8, {
        TT,SS,SS,SS,SS,SS, 0,TT,
        SS,BB,BB,LL,BB,BB,SS, 0,
        SS,BB,LL,BB,BB,BB,SS, 0,
        SS,BB,LL,BB,BB,BB,SS, 0,
        SS,BB,BB,BB,BB,BB,SS, 0,
        TT,SS,BB,BB,BB,SS, 0,TT,
        TT,TT,SS,BB,SS, 0,TT,TT,
        TT,TT,TT,SS, 0,TT,TT,TT,
}};
#undef SS
#undef BB
#undef LL
#define SS 0xba // shade
#define BB 0xb6 // base
#define LL 0xb2 // lit
static qpic_t shield_g = { 8, 8, {
        TT,SS,SS,SS,SS,SS, 0,TT,
        SS,BB,BB,LL,BB,BB,SS, 0,
        SS,BB,LL,BB,BB,BB,SS, 0,
        SS,BB,LL,BB,BB,BB,SS, 0,
        SS,BB,BB,BB,BB,BB,SS, 0,
        TT,SS,BB,BB,BB,SS, 0,TT,
        TT,TT,SS,BB,SS, 0,TT,TT,
        TT,TT,TT,SS, 0,TT,TT,TT,
}};
#undef SS
#undef BB
#undef LL
static qpic_t pent8x8 = { 8, 8, {
        0x4d,0x4a,0xff,0xff,0xff,0xff,0x4a,0x4d,
        0x4a,0x4f,0x4f,0x46,0x46,0x4f,0x4f,0x4a,
        0x46,0x4a,0x4f,0xf8,0xf8,0x4f,0x4a,0x46,
        0xff,0x4d,0xf8,0x48,0x48,0xf8,0x4d,0xff,
        0x4d,0x4f,0xf8,0x4b,0x4b,0xf8,0x4d,0x4b,
        0x46,0x46,0x4a,0x4f,0x4f,0x4a,0x46,0x46,
        0xff,0xff,0x46,0x4f,0x4d,0x46,0xff,0xff,
        0xff,0xff,0xff,0x4a,0x48,0xff,0xff,0xff,
}};
static qpic_t shell8x8 = { 8, 8, {
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xff,0x5e,0x6e,0x5e,0xff,0xff,0xff,
	0xff,0xff,0x63,0x68,0x63,0xff,0xff,0xff,
	0xff,0xff,0x63,0x68,0x63,0xff,0xff,0xff,
	0xff,0xff,0x4a,0x4f,0x4a,0xff,0xff,0xff,
	0xff,0xff,0x4a,0x4f,0x4a,0xff,0xff,0xff,
	0xff,0xff,0x49,0x4d,0x49,0xff,0xff,0xff,
}};
static qpic_t nail8x8 = { 8, 8, {
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xac,0xa8,0xa4,0xa8,0xac,0xff,0xff,
	0xff,0xff,0xad,0xab,0xad,0xff,0xff,0xff,
	0xff,0xff,0xff,0xa8,0xff,0xff,0xff,0xff,
	0xff,0xff,0xff,0xa7,0xff,0xff,0xff,0xff,
	0xff,0xff,0xff,0xa6,0xff,0xff,0xff,0xff,
	0xff,0xff,0xff,0xa9,0xff,0xff,0xff,0xff,
}};
static qpic_t lavanail8x8 = { 8, 8, {
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xf8,0xfa,0xfb,0xfa,0xf8,0xff,0xff,
	0xff,0xff,0xf7,0xf8,0xf7,0xff,0xff,0xff,
	0xff,0xff,0xf8,0xfa,0xf8,0xff,0xff,0xff,
	0xff,0xff,0xf9,0xfb,0xf9,0xff,0xff,0xff,
	0xff,0xff,0xf7,0xfa,0xf7,0xff,0xff,0xff,
	0xff,0xff,0xff,0xf8,0xff,0xff,0xff,0xff,
}};
static qpic_t rocket8x8 = { 8, 8, {
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xff,0xac,0xaa,0xac,0xff,0xff,0xff,
	0xff,0xff,0xab,0xa7,0xab,0xff,0xff,0xff,
	0xff,0xff,0x62,0x64,0x62,0xff,0xff,0xff,
	0xff,0xff,0xab,0xa7,0xab,0xff,0xff,0xff,
	0xff,0xad,0xab,0xa8,0xab,0xad,0xff,0xff,
	0xff,0xaa,0xac,0xa4,0xac,0xaa,0xff,0xff,
	0xff,0xab,0xad,0xa7,0xad,0xab,0xff,0xff,
}};
static qpic_t multirocket8x8 = { 8, 8, {
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xff,0xa9,0x00,0xa7,0x00,0xaa,0xff,
	0xff,0xff,0xa7,0x00,0xa4,0x00,0xa6,0xff,
	0xff,0xff,0x4e,0x00,0x4f,0x00,0x4d,0xff,
	0xff,0xff,0x4e,0x00,0x4f,0x00,0x4d,0xff,
	0xff,0xad,0x4d,0x00,0x4e,0x00,0x4c,0xad,
	0xff,0xa5,0xa7,0x00,0xa3,0x00,0xa8,0xa6,
	0xff,0xa8,0xaa,0x00,0xa5,0x00,0xaa,0xa9,
}};
/*static qpic_t cell8x8 = { 8, 8, {
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xad,0xa9,0xad,0xff,0xad,0xa9,0xad,
	0xff,0xaa,0x17,0xaa,0xaa,0x17,0xaa,0xaa,
	0xff,0x15,0x14,0x13,0x13,0x13,0x14,0x15,
	0xff,0x17,0x15,0x36,0x39,0x3c,0x15,0x18,
	0xff,0x18,0x15,0x14,0x14,0x14,0x15,0x18,
	0xff,0x17,0x18,0x18,0x17,0x18,0x18,0x17,
}};*/
static qpic_t lightning8x8 = { 8, 8, {
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xff,0x2a,0x2c,0x2c,0xff,0xff,0xff,
	0xff,0x27,0x2d,0x2e,0x27,0xff,0xff,0xff,
	0xff,0x2b,0x2f,0x2f,0x2d,0xff,0xff,0xff,
	0xff,0xff,0x2a,0x2e,0x27,0xff,0xff,0xff,
	0xff,0xff,0x2c,0x26,0xff,0xff,0xff,0xff,
	0xff,0x27,0xff,0xff,0xff,0xff,0xff,0xff,
}};
static qpic_t plasma8x8 = { 8, 8, {
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
	0xff,0x23,0x26,0x29,0x26,0x23,0xff,0xff,
	0xff,0x26,0x2b,0x2e,0x2b,0x26,0xff,0xff,
	0xff,0x29,0x2e,0x2f,0x2e,0x29,0xff,0xff,
	0xff,0x26,0x2b,0x2e,0x2b,0x26,0xff,0xff,
	0xff,0x23,0x26,0x29,0x26,0x23,0xff,0xff,
	0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
}};
static qpic_t *sb_nums[2][11];
static qpic_t *sb_colon, *sb_slash;
static qpic_t *sb_ibar;
static qpic_t *sb_sbar;
static qpic_t *sb_scorebar;
static qpic_t *sb_weapons[7][8]; // 0 is active, 1 is owned, 2-5 are flashes
static qpic_t *sb_ammo[4];
static qpic_t *sb_sigil[4];
static qpic_t *sb_armor[3];
static qpic_t *sb_items[32];
static qpic_t *sb_faces[7][2]; // 0 gibbed 1 dead 2-6 alive 0 static 1 temporary
static qpic_t *sb_face_invis;
static qpic_t *sb_face_quad;
static qpic_t *sb_face_invuln;
static qpic_t *sb_face_invis_invuln;
static qpic_t *rsb_invbar[2];
static qpic_t *rsb_weapons[5];
static qpic_t *rsb_items[2];
static qpic_t *rsb_ammo[3];
static qpic_t *rsb_teambord; // PGM 01/19/97 - team color border
static qpic_t *hsb_weapons[7][5]; // 0 is active, 1 is owned, 2-5 are flashes
static bool sb_showscores;
static qpic_t *hsb_items[2]; //MED 01/04/97 added hipnotic items array
static s32 fragsort[MAX_SCOREBOARD];
static s32 scoreboardlines;
static s32 npos[4][2]; // ammo count num
static s32 wpos[9][2]; // weapons
static s32 kpos[2][2]; // keys
static s32 iposx[12]; // classic items
static s32 iposy;
static s32 hipweapons[4] = //MED 01/04/97 added array to simplify weapon parsing
{ HIT_LASER_CANNON_BIT, HIT_MJOLNIR_BIT, 4, HIT_PROXIMITY_GUN_BIT };

inline static s32 Sbar_ColorForMap(s32 m)
{ return m < 128 ? m + 8 : m + 8; }
inline static void Sbar_ShowScores() // Tab key down
{ if (sb_showscores) return; sb_showscores = 1; sb_updates = 0; }
inline static void Sbar_DontShowScores() // Tab key up
{ sb_showscores = 0; sb_updates = 0; }
void Sbar_Changed() // update next frame
{ sb_updates = 0; }

void Sbar_Init()
{
	for (s32 i = 0; i < 10; i++) {
		sb_nums[0][i] = Draw_PicFromWad(va("num_%i", i));
		sb_nums[1][i] = Draw_PicFromWad(va("anum_%i", i));
	}
	sb_nums[0][10] = Draw_PicFromWad("num_minus");
	sb_nums[1][10] = Draw_PicFromWad("anum_minus");
	sb_colon = Draw_PicFromWad("num_colon");
	sb_slash = Draw_PicFromWad("num_slash");
	sb_weapons[0][0] = Draw_PicFromWad("inv_shotgun");
	sb_weapons[0][1] = Draw_PicFromWad("inv_sshotgun");
	sb_weapons[0][2] = Draw_PicFromWad("inv_nailgun");
	sb_weapons[0][3] = Draw_PicFromWad("inv_snailgun");
	sb_weapons[0][4] = Draw_PicFromWad("inv_rlaunch");
	sb_weapons[0][5] = Draw_PicFromWad("inv_srlaunch");
	sb_weapons[0][6] = Draw_PicFromWad("inv_lightng");
	sb_weapons[1][0] = Draw_PicFromWad("inv2_shotgun");
	sb_weapons[1][1] = Draw_PicFromWad("inv2_sshotgun");
	sb_weapons[1][2] = Draw_PicFromWad("inv2_nailgun");
	sb_weapons[1][3] = Draw_PicFromWad("inv2_snailgun");
	sb_weapons[1][4] = Draw_PicFromWad("inv2_rlaunch");
	sb_weapons[1][5] = Draw_PicFromWad("inv2_srlaunch");
	sb_weapons[1][6] = Draw_PicFromWad("inv2_lightng");
	for (s32 i = 0; i < 5; i++) {
		sb_weapons[2+i][0] = Draw_PicFromWad(va("inva%i_shotgun",i+1));
		sb_weapons[2+i][1] = Draw_PicFromWad(va("inva%i_sshotgun",i+1));
		sb_weapons[2+i][2] = Draw_PicFromWad(va("inva%i_nailgun",i+1));
		sb_weapons[2+i][3] = Draw_PicFromWad(va("inva%i_snailgun",i+1));
		sb_weapons[2+i][4] = Draw_PicFromWad(va("inva%i_rlaunch",i+1));
		sb_weapons[2+i][5] = Draw_PicFromWad(va("inva%i_srlaunch",i+1));
		sb_weapons[2+i][6] = Draw_PicFromWad(va("inva%i_lightng",i+1));
	}
	sb_ammo[0] = Draw_PicFromWad("sb_shells");
	sb_ammo[1] = Draw_PicFromWad("sb_nails");
	sb_ammo[2] = Draw_PicFromWad("sb_rocket");
	sb_ammo[3] = Draw_PicFromWad("sb_cells");
	sb_armor[0] = Draw_PicFromWad("sb_armor1");
	sb_armor[1] = Draw_PicFromWad("sb_armor2");
	sb_armor[2] = Draw_PicFromWad("sb_armor3");
	sb_items[0] = Draw_PicFromWad("sb_key1");
	sb_items[1] = Draw_PicFromWad("sb_key2");
	sb_items[2] = Draw_PicFromWad("sb_invis");
	sb_items[3] = Draw_PicFromWad("sb_invuln");
	sb_items[4] = Draw_PicFromWad("sb_suit");
	sb_items[5] = Draw_PicFromWad("sb_quad");
	sb_sigil[0] = Draw_PicFromWad("sb_sigil1");
	sb_sigil[1] = Draw_PicFromWad("sb_sigil2");
	sb_sigil[2] = Draw_PicFromWad("sb_sigil3");
	sb_sigil[3] = Draw_PicFromWad("sb_sigil4");
	sb_faces[4][0] = Draw_PicFromWad("face1");
	sb_faces[4][1] = Draw_PicFromWad("face_p1");
	sb_faces[3][0] = Draw_PicFromWad("face2");
	sb_faces[3][1] = Draw_PicFromWad("face_p2");
	sb_faces[2][0] = Draw_PicFromWad("face3");
	sb_faces[2][1] = Draw_PicFromWad("face_p3");
	sb_faces[1][0] = Draw_PicFromWad("face4");
	sb_faces[1][1] = Draw_PicFromWad("face_p4");
	sb_faces[0][0] = Draw_PicFromWad("face5");
	sb_faces[0][1] = Draw_PicFromWad("face_p5");
	sb_face_invis = Draw_PicFromWad("face_invis");
	sb_face_invuln = Draw_PicFromWad("face_invul2");
	sb_face_invis_invuln = Draw_PicFromWad("face_inv2");
	sb_face_quad = Draw_PicFromWad("face_quad");
	sb_sbar = Draw_PicFromWad("sbar");
	sb_ibar = Draw_PicFromWad("ibar");
	sb_scorebar = Draw_PicFromWad("scorebar");
	if (hipnotic) {
		hsb_weapons[0][0] = Draw_PicFromWad("inv_laser");
		hsb_weapons[0][1] = Draw_PicFromWad("inv_mjolnir");
		hsb_weapons[0][2] = Draw_PicFromWad("inv_gren_prox");
		hsb_weapons[0][3] = Draw_PicFromWad("inv_prox_gren");
		hsb_weapons[0][4] = Draw_PicFromWad("inv_prox");
		hsb_weapons[1][0] = Draw_PicFromWad("inv2_laser");
		hsb_weapons[1][1] = Draw_PicFromWad("inv2_mjolnir");
		hsb_weapons[1][2] = Draw_PicFromWad("inv2_gren_prox");
		hsb_weapons[1][3] = Draw_PicFromWad("inv2_prox_gren");
		hsb_weapons[1][4] = Draw_PicFromWad("inv2_prox");
		hsb_items[0] = Draw_PicFromWad("sb_wsuit");
		hsb_items[1] = Draw_PicFromWad("sb_eshld");
	}
	for (s32 i = 0; hipnotic && i < 5; i++) {
		hsb_weapons[2+i][0]=Draw_PicFromWad(va("inva%i_laser",i+1));
		hsb_weapons[2+i][1]=Draw_PicFromWad(va("inva%i_mjolnir",i+1));
		hsb_weapons[2+i][2]=Draw_PicFromWad(va("inva%i_gren_prox",i+1));
		hsb_weapons[2+i][3]=Draw_PicFromWad(va("inva%i_prox_gren",i+1));
		hsb_weapons[2+i][4]=Draw_PicFromWad(va("inva%i_prox",i+1));
	}
	if (rogue) {
		rsb_invbar[0] = Draw_PicFromWad("r_invbar1");
		rsb_invbar[1] = Draw_PicFromWad("r_invbar2");
		rsb_weapons[0] = Draw_PicFromWad("r_lava");
		rsb_weapons[1] = Draw_PicFromWad("r_superlava");
		rsb_weapons[2] = Draw_PicFromWad("r_gren");
		rsb_weapons[3] = Draw_PicFromWad("r_multirock");
		rsb_weapons[4] = Draw_PicFromWad("r_plasma");
		rsb_items[0] = Draw_PicFromWad("r_shield1");
		rsb_items[1] = Draw_PicFromWad("r_agrav1");
		rsb_teambord = Draw_PicFromWad("r_teambord");//team color border
		rsb_ammo[0] = Draw_PicFromWad("r_ammolava");
		rsb_ammo[1] = Draw_PicFromWad("r_ammomulti");
		rsb_ammo[2] = Draw_PicFromWad("r_ammoplasma");
	}
	if (host_initialized) return;
	Cmd_AddCommand("+showscores", Sbar_ShowScores);
	Cmd_AddCommand("-showscores", Sbar_DontShowScores);
	Cvar_RegisterVariable(&scr_sidescore);
}

s32 Sbar_itoa(s32 num, s8 *buf)
{
	s8 *str = buf;
	if (num < 0) {
		*str++ = '-';
		num = -num;
	}
	s32 pow10 = 10;
	for (; num >= pow10; pow10 *= 10);
	do {
		pow10 /= 10;
		s32 dig = num / pow10;
		*str++ = '0' + dig;
		num -= dig * pow10;
	} while (pow10 != 1);
	*str = 0;
	return str - buf;
}

void Sbar_DrawNum(s32 x, s32 y, s32 num, s32 digits, s32 color)
{
	s8 str[12];
	s32 l = Sbar_itoa(num, str);
	s8 *ptr = str;
	if (l > digits)
		ptr += (l - digits);
	if (l < digits)
		x += (digits - l) * 24*SCL;
	while (*ptr) {
		s32 frame = *ptr == '-' ? 10 : *ptr - '0'; // 10 = '-' s8
		Draw_TransPicScaled(x, y, sb_nums[color][frame], SCL);
		x += 24*SCL;
		ptr++;
	}
}

void Sbar_DrawNumSmall(s32 x, s32 y, s32 num, s32 color)
{
	s8 buf[6];
	sprintf(buf, "%3i", num);
	if      (color == 0) color = 48;
	else if (color == 1) color = 18;
	else                 color = 176;
	if (buf[0] != ' ')
		Draw_CharacterScaled(x       , y, color+buf[0]-'0', SCL);
	if (buf[1] != ' ')
		Draw_CharacterScaled(x+8*SCL , y, color+buf[1]-'0', SCL);
	if (buf[2] != ' ')
		Draw_CharacterScaled(x+16*SCL, y, color+buf[2]-'0', SCL);
}

void Sbar_SortFrags()
{
	scoreboardlines = 0; // sort by frags
	for (s32 i = 0; i < cl.maxclients; i++) {
		if (cl.scores[i].name[0]) {
			fragsort[scoreboardlines] = i;
			scoreboardlines++;
		}
	}
	for (s32 i = 0; i < scoreboardlines; i++)
		for (s32 j = 0; j < scoreboardlines - 1 - i; j++)
			if (cl.scores[fragsort[j]].frags <
					cl.scores[fragsort[j + 1]].frags) {
				s32 k = fragsort[j];
				fragsort[j] = fragsort[j + 1];
				fragsort[j + 1] = k;
			}
}

void Sbar_SoloScoreboard()
{
	s8 str[80];
	sprintf(str, "Monsters:%3i /%3i", cl.stats[STAT_MONSTERS],
			cl.stats[STAT_TOTALMONSTERS]);
	Draw_StringScaled(WW/2-152*SCL, HH-20*SCL, str, SCL);
	sprintf(str, "Secrets :%3i /%3i", cl.stats[STAT_SECRETS],
			cl.stats[STAT_TOTALSECRETS]);
	Draw_StringScaled(WW/2-152*SCL, HH-12*SCL, str, SCL);
	s32 minutes = cl.time / 60;
	s32 seconds = cl.time - 60 * minutes;
	s32 tens = seconds / 10;
	s32 units = seconds - 10 * tens;
	sprintf(str, "Time :%3i:%i%i", minutes, tens, units);
	Draw_StringScaled(WW/2+24*SCL, HH-20*SCL, str, SCL);
	Draw_StringScaled(WW/2+24*SCL, HH-12*SCL, cl.levelname, SCL);
}

void Sbar_CalcPos()
{
	s32 arcade_offs_x = scr_hudstyle.value == 4 ? 160*SCL : 0;
	s32 arcade_offs_y = scr_hudstyle.value == 4 ? 24*SCL : 0;
	for (s32 i = 0; i < 8; i++) // items classic
		iposx[i] = WW / 2 + 32*SCL + i*16*SCL + arcade_offs_x;
	for (s32 i = 8; i < 12; i++)
		iposx[i] = WW / 2 + 64*SCL + i*8*SCL + arcade_offs_x;
	iposy = HH - 40*SCL + arcade_offs_y;
	switch ((s32)(hipnotic&&scr_hudstyle.value==3?2:scr_hudstyle.value)) {
	default: // ammo counters
	case 8: // classic (no bg)
	case 0: // classic
		for (s32 i = 0; i < 4; i++) {
			npos[i][0] = WW / 2 - 150*SCL + i*48*SCL;
			npos[i][1] = HH - 48*SCL;
		}
		break;
	case 1: // bottom center, 4x1
		for (s32 i = 0; i < 4; i++) {
			npos[i][0] = WW / 2 - 88*SCL + i*48*SCL;
			npos[i][1] = HH - 8 * SCL;
		}
		break;
	case 2: // right side, 2x2
		for (s32 i = 0; i < 4; i++) {
			npos[i][0] = WW - (scr_hudstyle.value==3?88:96)
				*SCL + 48*SCL*(i&1);
			npos[i][1] = HH - 47*SCL - 9*SCL*(i>>1);
		}
		break;
	case 3: // right side, 1x4
		for (s32 i = 0; i < 4; i++) {
			npos[i][0] = WW - 48*SCL + 8*SCL;
			npos[i][1] = HH - 78*SCL + i * 9*SCL;
		}
		break;
	case 4: // arcade
		for (s32 i = 0; i < 4; i++) {
			npos[i][0] = WW / 2 + 10*SCL + i*48*SCL;
			npos[i][1] = HH - 24*SCL;
		}
		break;
	case 9: // EZQuake
		for (s32 i = 0; i < 4; i++) {
			npos[i][0] = WW / 2 - 68*SCL + i*32*SCL;
			npos[i][1] = HH - 32*SCL;
		}
		break;
		break;
	}
	switch ((s32)scr_hudstyle.value) { // weapons
	case 8: // classic (no bg)
	case 0: // classic
	case 4: // arcade
		for (s32 i = 0; i < 9; i++) {
			wpos[i][0] = WW / 2 - 160*SCL + i*24*SCL +arcade_offs_x;
			wpos[i][1] = HH - 40*SCL + arcade_offs_y;
			if (i == 9) {
				wpos[7][0] += 8*SCL + arcade_offs_x;
				wpos[8][0] += 8*SCL + arcade_offs_y;
			}
		}
		break;
	case 9: // EZQuake
		for (s32 i = 0; i < 9; i++) {
			wpos[i][0] = WW / 2 - 72*SCL + i*16*SCL;
			wpos[i][1] = HH - 48*SCL;
			if (i == 9) {
				wpos[7][0] += 8*SCL;
				wpos[8][0] += 8*SCL;
			}
		}
		break;
	default: // modern
		for (s32 i = 0; i < 9; i++) {
			wpos[i][0] = WW - 18*SCL;	
			wpos[i][1] = HH- 16*SCL*(6+i) + hipnotic*24*SCL;
		}
		break;
	}
	switch ((s32)scr_hudstyle.value + 0x10*hipnotic) { // keys
	default:
	case 8: // classic (no bg)
	case 0:
	case 3:
		for (s32 i = 0; i < 2; i++) { // classic, qw
			kpos[i][0] = WW / 2 + 32*SCL + i*16*SCL + arcade_offs_x;
			kpos[i][1] = HH - 40*SCL + arcade_offs_y;
		}
		break;
	case 1:
		for (s32 i = 0; i < 2; i++) { // side, vertical
			kpos[i][0] = WW - 28*SCL;
			kpos[i][1] = HH - 52*SCL - i*16*SCL;
		}
		break;
	case 2:
		for (s32 i = 0; i < 2; i++) { // side, horizontal
			kpos[i][0] = WW - 24*SCL - i*16*SCL;
			kpos[i][1] = HH - 75*SCL;
		}
		break;
	case 0x10: // hipnotic
	case 0x13:
		for (s32 i = 0; i < 2; i++) { // classic, qw
			kpos[i][0] = WW / 2 + 48*SCL;
			kpos[i][1] = HH - 21*SCL + i*8*SCL;
		}
		break;
	case 0x11:
		for (s32 i = 0; i < 2; i++) { // side, vertical
			kpos[i][0] = WW - 28*SCL;
			kpos[i][1] = HH - 44*SCL - i*8*SCL;
		}
		break;
	case 0x12:
		for (s32 i = 0; i < 2; i++) { // side, horizontal
			kpos[i][0] = WW - 40*SCL - i*16*SCL;
			kpos[i][1] = HH - 68*SCL;
		}
		break;
	case 0x14:
		for (s32 i = 0; i < 2; i++) { // arcade
			kpos[i][0] = WW / 2 - 112*SCL;
			kpos[i][1] = HH - 21*SCL + i*8*SCL;
		}
		break;
	}
}

void Sbar_ItemsClassic()
{
	for (s32 i = 2; i < 6; i++)
		if (cl.items & (1 << (17 + i))) 
			Draw_TransPicScaled(iposx[i],iposy,sb_items[i],SCL);
	for (s32 i = 0; hipnotic && i < 2; i++)
		if (cl.items & (1 << (24 + i)))
			Draw_TransPicScaled(iposx[i+6],iposy,hsb_items[i],SCL);
	for (s32 i = 0; rogue && i < 2; i++)
		if (cl.items & (1 << (29 + i)))
			Draw_TransPicScaled(iposx[i+6],iposy,rsb_items[i],SCL);
	for (s32 i = 0; !rogue && i < 4; i++) // sigils
		if (cl.items & (1 << (28 + i)))
			Draw_TransPicScaled(iposx[8+i], iposy, sb_sigil[i], SCL);
}

void Sbar_ItemsModern()
{
	s32 x = 12 * SCL;
	s32 y = HH - 48 * SCL;
	if (cl.items & IT_INVULNERABILITY || cl.stats[STAT_ARMOR] > 0)
		y -= 24*SCL; // armor is visible, start above it
	for (s32 i = 2; i < 6; i++)
		if ((cl.items & (1<<(17+i))) && (!hipnotic || (i > 1))) {
			Draw_TransPicScaled(x, y, sb_items[i], SCL);
			y -= 16*SCL;
		}
	for (s32 i = 0; hipnotic && i < 2; i++)
		if (cl.items & (1<<(24+i))) {
			Draw_TransPicScaled(x, y, hsb_items[i], SCL);
			y -= 16*SCL;
		}
	for (s32 i = 0; rogue && i < 2; i++)
		if (cl.items & (1<<(29+i))) {
			Draw_TransPicScaled(x, y, rsb_items[i], SCL);
			y -= 16*SCL;
		}
}

void Sbar_DrawInventoryBg()
{
	if (scr_hudstyle.value == 9) return; // EZQuake
	qpic_t *pic = rogue ?
		rsb_invbar[cl.stats[STAT_ACTIVEWEAPON] < RIT_LAVA_NAILGUN] :
		sb_ibar;
	s32 x, y;
	switch ((s32)(hipnotic&&scr_hudstyle.value==3?2:scr_hudstyle.value)) {
	default:
	case 0: // classic
		x = WW / 2 - 160*SCL;
		y = HH - 48*SCL;
		Draw_PicScaled(x, y, pic, SCL);
		break;
	case 1: // bottom center, 4x1
		x = WW / 2 - 96*SCL;
		y = HH - 8*SCL;
		Draw_PicScaledPartial(x, y, 0, 0, 192, 8, pic, SCL);
		break;
	case 2: // right side, 2x2
		x = WW - (scr_hudstyle.value==3 ? 96 : 104)*SCL;
		y = HH - 56*SCL;
		Draw_PicScaledPartial(x, y + 9*SCL, 0,0,96,8, pic, SCL);
		Draw_PicScaledPartial(x - 96*SCL,y,96,8,192,16,pic,SCL);
		break;
	case 3: // right side, 1x4
		x = WW - 48*SCL;
		y = HH - 78*SCL;
		for (s32 i = 0; i < 4; i++)
			Draw_PicScaledPartial(x - 48*SCL*i, y + 9*SCL*i,
					i*48, 0, (i+1)*48, 8, pic, SCL);
		break;
	case 4: // arcade
		x = WW / 2;
		y = HH - 24*SCL;
		Draw_PicScaledPartial(x, y, 0, 0, 320, 24, pic, SCL);
		break;
	}
}

static inline s32 Sbar_Flash(f32 time, s32 active)
{
	s32 flash = (s32)(time * 10);
	return flash >= 10 ? active : (flash % 5) + 2; }

void Sbar_DrawWeapons()
{
	for (s32 i = 0; i < 7; i++) {
		if ((hipnotic && i==IT_GRENADE_LAUNCHER) || !(cl.items&(1<<i)))
			continue;
		s32 active = (cl.stats[STAT_ACTIVEWEAPON] == (1 << i));
		s32 flashon = Sbar_Flash(cl.time - cl.item_gettime[i], active);
		qpic_t *pic = sb_weapons[flashon][i];
		// Check for rogue's powered up weapon condition.
		if (rogue && i >= 2 && cl.stats[STAT_ACTIVEWEAPON] ==
				(RIT_LAVA_NAILGUN << (i - 2))) {
			pic = rsb_weapons[i - 2]; // powered up weapon
			active = 1;
		}
		active = scr_hudstyle.value ? active * 6 * SCL : 0;
		if (scr_hudstyle.value == 4) active = 0; // arcade
		Draw_TransPicScaled(wpos[i][0] - active, wpos[i][1], pic, SCL);
		if (flashon > 1)
			sb_updates = 0; // force update to remove flash
	}
	s32 grenadeflashing = 0;
	for (s32 i = 0; hipnotic && i < 4; i++) {
		if (!(cl.items & (1 << hipweapons[i])))
			continue;
		s32 active = (cl.stats[STAT_ACTIVEWEAPON]==(1<<hipweapons[i]));
		s32 flashon = Sbar_Flash(cl.time-cl.item_gettime[hipweapons[i]],
			active);
		active = scr_hudstyle.value ? active * 6 * SCL : 0;
		if (scr_hudstyle.value == 4) active = 0; // arcade
		s32 x = wpos[4][0] - active;
		s32 y = wpos[4][1];
		qpic_t *pic = NULL;
		if (i == 2) { 
			if (cl.items & HIT_PROXIMITY_GUN && flashon) {
				grenadeflashing = 1;
				pic = hsb_weapons[flashon][2];
			}
		} else if (i == 3) {
			if (cl.items & (1 << 4)) {
				if (flashon && !grenadeflashing)
					pic = hsb_weapons[flashon][3];
				else if (!grenadeflashing)
					pic = hsb_weapons[0][3];
			} else
				pic = hsb_weapons[flashon][4];
		} else {
			x = wpos[7 + i][0] - active;
			y = wpos[7 + i][1];
			pic = hsb_weapons[flashon][i];
		}
		if (pic)
			Draw_TransPicScaled(x, y, pic, SCL);
		if (flashon > 1)
			sb_updates = 0;
	}
}

void Sbar_DrawInventory()
{
	Sbar_DrawInventoryBg();
	Sbar_DrawWeapons();
	if (scr_hudstyle.value == 9) {
    Sbar_DrawNumSmall(npos[0][0], npos[0][1], cl.stats[6], cl.stats[6] <= 10 ? 2 : 0);
    Sbar_DrawNumSmall(npos[1][0], npos[1][1], cl.stats[7], cl.stats[7] <= 10 ? 2 : 0);
    Sbar_DrawNumSmall(npos[2][0], npos[2][1], cl.stats[8], cl.stats[8] <= 10 ? 2 : 0);
    Sbar_DrawNumSmall(npos[3][0], npos[3][1], cl.stats[9], cl.stats[9] <= 10 ? 2 : 0);
		return;
	}
	for (s32 i = 0; i < 4; i++) // ammo counters
		Sbar_DrawNumSmall(npos[i][0], npos[i][1], cl.stats[6 + i], 1);
	for (s32 i = 0; i < 2; i++) { // keys
		if (!(cl.items & (1<<(17+i))))
			continue;
		Draw_TransPicScaled(kpos[i][0], kpos[i][1], sb_items[i], SCL);
	}
	switch ((s32)scr_hudstyle.value) {
	default: case 0: case 3: Sbar_ItemsClassic(); break;
	         case 1: case 2: Sbar_ItemsModern();  break;
	}
}

void Sbar_DrawFace() // ...and HP
{
	s32 x, y;
	switch ((s32)scr_hudstyle.value) {
	case 1: // modern
	case 2:
		x = 32*SCL;
		y = HH - 32*SCL;
		break;
	case 4: // arcade
		x = WW / 2 - 184*SCL;
		y = HH - 24*SCL;
		break;
	case 9: // EZQuake
		x = WW / 2 + 12*SCL;
		y = HH - 24*SCL;
		break;
	default: // classic, qw
		x = WW / 2 - 24*SCL; // classic, qw
		y = HH - 24*SCL;
		break;
	}
	Sbar_DrawNum(x, y, cl.stats[STAT_HEALTH], 3, cl.stats[STAT_HEALTH]<=25);
	if (scr_hudstyle.value == 9) return;
	switch ((s32)scr_hudstyle.value) {
	case 1: // modern
	case 2:
		x = 8*SCL;
		y = HH - 32*SCL;
		break;
	case 4: // arcade
		x = WW / 2 - 208*SCL;
		y = HH - 24*SCL;
		break;
	default: // classic, qw
		x = WW / 2 - 48*SCL; // classic, qw
		y = HH - 24*SCL;
		break;
	}
	if (rogue && cl.maxclients!=1 && teamplay.value>3 && teamplay.value<7
			&& (!scr_hudstyle.value || scr_hudstyle.value == 3)) {
		scoreboard_t *s = &cl.scores[cl.viewentity - 1];
		s32 top = Sbar_ColorForMap(s->colors & 0xf0);
		s32 bottom = Sbar_ColorForMap((s->colors & 0x0f) << 4);
		Draw_TransPicScaled(x, y, rsb_teambord, SCL);
		Draw_Fill(WW/2 - 47*SCL, HH - 21*SCL, 22*SCL, 9*SCL, top);
		Draw_Fill(WW/2 - 47*SCL, HH - 12*SCL, 22*SCL, 9*SCL, bottom);
		Sbar_DrawNumSmall(WW/2-49*SCL, HH-21*SCL, s->frags, top==8);
		return;
	}
	if ((cl.items & (IT_INVISIBILITY | IT_INVULNERABILITY))
			== (IT_INVISIBILITY | IT_INVULNERABILITY)) {
		Draw_TransPicScaled(x, y, sb_face_invis_invuln, SCL);
		return;
	}
	if (cl.items & IT_QUAD) {
		Draw_TransPicScaled(x, y, sb_face_quad, SCL);
		return;
	}
	if (cl.items & IT_INVISIBILITY) {
		Draw_TransPicScaled(x, y, sb_face_invis, SCL);
		return;
	}
	if (cl.items & IT_INVULNERABILITY) {
		Draw_TransPicScaled(x, y, sb_face_invuln, SCL);
		return;
	}
	s32 f = cl.stats[STAT_HEALTH] >= 100 ? 4 : cl.stats[STAT_HEALTH] / 20;
	s32 anim = 0;
	if (cl.time <= cl.faceanimtime) {
		anim = 1;
		sb_updates = 0; // make sure the anim gets drawn over
	}
	Draw_TransPicScaled(x, y, sb_faces[f][anim], SCL);
}

void Sbar_DrawArmor()
{
	s32 x, y;
	switch ((s32)scr_hudstyle.value) {
		case 1: // modern
		case 2:
			x = 8 * SCL;
			y = HH - 56 * SCL;
			break;
		case 4: // arcade
			x = WW / 2 - 320*SCL;
			y = HH - 24*SCL;
			break;
		case 9: // EZQuake
			x = WW / 2 - 96*SCL;
			y = HH - 24*SCL;
			break;
		default: // classic, qw
			x = WW / 2 - 160*SCL;
			y = HH - 24*SCL;
			break;
	}
	if (cl.items & IT_INVULNERABILITY) { // armor
		Sbar_DrawNum(24*SCL + x, y, 666, 3, 1);
		Draw_TransPicScaled(x, y, draw_disc, SCL);
		return;
	}
	if (rogue) {
		if (!scr_hudstyle.value || cl.stats[STAT_ARMOR])
			Sbar_DrawNum(24*SCL + x, y, cl.stats[STAT_ARMOR], 3,
					cl.stats[STAT_ARMOR] <= 25);
		if (cl.items & RIT_ARMOR3)
			Draw_TransPicScaled(x, y, sb_armor[2], SCL);
		else if (cl.items & RIT_ARMOR2)
			Draw_TransPicScaled(x, y, sb_armor[1], SCL);
		else if (cl.items & RIT_ARMOR1)
			Draw_TransPicScaled(x, y, sb_armor[0], SCL);
	} else {
		if (!scr_hudstyle.value || cl.stats[STAT_ARMOR])
			Sbar_DrawNum(24*SCL + x, y, cl.stats[STAT_ARMOR], 3,
					cl.stats[STAT_ARMOR] <= 25);
		if (cl.items & IT_ARMOR3)
			Draw_TransPicScaled(x, y, sb_armor[2], SCL);
		else if (cl.items & IT_ARMOR2)
			Draw_TransPicScaled(x, y, sb_armor[1], SCL);
		else if (cl.items & IT_ARMOR1)
			Draw_TransPicScaled(x, y, sb_armor[0], SCL);
	}
}

void Sbar_DrawAmmo()
{
	s32 x = WW / 2 + 64*SCL; // classic, qw
	s32 y = HH - 24*SCL;
	if (scr_hudstyle.value == 1 || scr_hudstyle.value == 2) { // modern
		x = WW - 32*SCL;
		y = HH - 32*SCL;
	} else if (scr_hudstyle.value == 4) { // arcade
		x = WW / 2 - 96*SCL;
		y = HH - 24*SCL;
	}
	if (rogue) { // ammo icon
		if (cl.items & RIT_SHELLS)
			Draw_TransPicScaled(x, y, sb_ammo[0], SCL);
		else if (cl.items & RIT_NAILS)
			Draw_TransPicScaled(x, y, sb_ammo[1], SCL);
		else if (cl.items & RIT_ROCKETS)
			Draw_TransPicScaled(x, y, sb_ammo[2], SCL);
		else if (cl.items & RIT_CELLS)
			Draw_TransPicScaled(x, y, sb_ammo[3], SCL);
		else if (cl.items & RIT_LAVA_NAILS)
			Draw_TransPicScaled(x, y, rsb_ammo[0], SCL);
		else if (cl.items & RIT_PLASMA_AMMO)
			Draw_TransPicScaled(x, y, rsb_ammo[1], SCL);
		else if (cl.items & RIT_MULTI_ROCKETS)
			Draw_TransPicScaled(x, y, rsb_ammo[2], SCL);
	} else {
		if (cl.items & IT_SHELLS)
			Draw_TransPicScaled(x, y, sb_ammo[0], SCL);
		else if (cl.items & IT_NAILS)
			Draw_TransPicScaled(x, y, sb_ammo[1], SCL);
		else if (cl.items & IT_ROCKETS)
			Draw_TransPicScaled(x, y, sb_ammo[2], SCL);
		else if (cl.items & IT_CELLS)
			Draw_TransPicScaled(x, y, sb_ammo[3], SCL);
	}
	x = WW / 2 + 88*SCL; // classic, qw
	y = HH - 24*SCL;
	if (scr_hudstyle.value == 1 || scr_hudstyle.value == 2) { // modern
		x = WW - 108*SCL;
		y = HH - 32*SCL;
	} else if (scr_hudstyle.value == 4) { // arcade
		x = WW / 2 - 72*SCL;
		y = HH - 24*SCL;
	}
	Sbar_DrawNum(x, y, cl.stats[STAT_AMMO], 3, cl.stats[STAT_AMMO] <= 10);
}

static inline void Sbar_PlayerScoreBox(s32 x, s32 y, scoreboard_t *s, s32 k)
{
	s32 top = Sbar_ColorForMap(s->colors & 0xf0);
	s32 bottom = Sbar_ColorForMap((s->colors & 0x0f) << 4);
	Draw_Fill(x, y, 30*SCL, 4*SCL, top);
	Draw_Fill(x, y+4*SCL, 30*SCL, 3*SCL, bottom);
	Sbar_DrawNumSmall(x+3*SCL, y-1*SCL, s->frags, 0);
	if (k == cl.viewentity - 1) {
		Draw_CharacterScaled(x-2*SCL,y-1*SCL,16,SCL);
		Draw_CharacterScaled(x+23*SCL,y-1*SCL,17,SCL);
	}
}

void Sbar_DeathmatchOverlay()
{
	Sbar_SortFrags();
	s32 x = !scr_hudstyle.value ? WW/2 + 168*SCL : 4*SCL;
	s32 y = !scr_hudstyle.value ? HH - 48*SCL : 40*SCL;
	if (sb_showscores) {
		qpic_t *pic = Draw_CachePic("gfx/ranking.lmp");
		Draw_TransPicScaled((WW-pic->width*SCL)/2, 8*SCL, pic, SCL);
		x = WW/2 - 80*SCL;
		y = 40*SCL;
	}
	for (s32 i = 0; i < scoreboardlines; i++) {
		s32 k = fragsort[i];
		scoreboard_t *s = &cl.scores[k];
		if (!s->name[0])
			continue;
		Sbar_PlayerScoreBox(x, y, s, k);
		Draw_StringScaled(x+40*SCL, y-1*SCL, s->name, SCL);
		y += 10*SCL;
	}
}

void Sbar_DrawFrags()
{
	Sbar_SortFrags();
	s32 x = WW/2 + 33*SCL;
	for (s32 i = 0; i < (scoreboardlines<=4?scoreboardlines:4); i++) {
		s32 k = fragsort[i];
		scoreboard_t *s = &cl.scores[k];
		if (!s->name[0])
			continue;
		Sbar_PlayerScoreBox(x, HH - 47*SCL, s, k);
		x += 32*SCL;
	}
}

void Sbar_IntermissionOverlay()
{
	scr_fullupdate = 0;
	if (cl.gametype == GAME_DEATHMATCH) {
		Sbar_DeathmatchOverlay();
		return;
	}
	drawlayer = lyr_sbar.value;
	qpic_t *pic = Draw_CachePic("gfx/complete.lmp"); // plaque is 192px wide
	Draw_TransPicScaled(WW/2 - 96*SCL, 24*SCL,pic,SCL);
	s32 p = WW/2 - 160*SCL; // padding for scaling
	pic = Draw_CachePic("gfx/inter.lmp");
	Draw_TransPicScaled(p, 56*SCL, pic, SCL);
	s32 dig = cl.completed_time / 60;
	Sbar_DrawNum(WW / 2, 64*SCL, dig, 3, 0);
	s32 num = cl.completed_time - dig * 60;
	Draw_TransPicScaled(WW-p-86*SCL, 64*SCL, sb_colon, SCL);
	Draw_TransPicScaled(WW-p-74*SCL, 64*SCL, sb_nums[0][num/10], SCL);
	Draw_TransPicScaled(WW-p-54*SCL, 64*SCL, sb_nums[0][num%10], SCL);
	Sbar_DrawNum(WW/2, 104*SCL, cl.stats[STAT_SECRETS], 3, 0);
	Draw_TransPicScaled(WW-p-88*SCL, 104*SCL, sb_slash, SCL);
	Sbar_DrawNum(WW-p-80*SCL, 104*SCL, cl.stats[STAT_TOTALSECRETS], 3, 0);
	Sbar_DrawNum(WW/2, 144*SCL, cl.stats[STAT_MONSTERS], 3, 0);
	Draw_TransPicScaled(WW-p-88*SCL, 144*SCL, sb_slash, SCL);
	Sbar_DrawNum(WW-p-80*SCL, 144*SCL, cl.stats[STAT_TOTALMONSTERS], 3, 0);
	drawlayer = lyr_main.value;
}

void Sbar_FinaleOverlay()
{
	drawlayer = lyr_sbar.value;
	qpic_t *pic = Draw_CachePic("gfx/finale.lmp");
	Draw_TransPicScaled(WW/2 - pic->width/2*SCL, 16*SCL, pic, SCL);
	drawlayer = lyr_main.value;
}

void Sbar_DrawBg()
{ // only for the armor/health/ammo half
	if (scr_hudstyle.value == 0) { // classic
		Draw_TransPicScaled(WW/2 - 160*(SCL), HH-24*SCL, sb_sbar, SCL);
	} else if (scr_hudstyle.value == 4) { // arcade
		s32 x = WW / 2 - 320*SCL;
		s32 y = HH - 24*SCL;
		Draw_PicScaledPartial(x, y, 0, 0, 320, 24, sb_sbar, SCL);
	}
}

static inline s32 int_offs(s32 i)
{
	if (i >= 100) return 0;
	if (i >= 10) return 4;
	return 8;
}

void Sbar_Min(s32 color)
{
	if (cl.gametype==GAME_DEATHMATCH&&(sb_showscores||scr_sidescore.value)){
		Sbar_DeathmatchOverlay();
	}
	if (sb_showscores || cl.stats[STAT_HEALTH] <= 0) {
		Draw_TransPicScaled(WW/2-160*SCL, HH-24*SCL, sb_scorebar, SCL);
		Sbar_SoloScoreboard();
		sb_updates = 0;
		return;
	}
	s32 x = WW/2 - 12*SCL - int_offs(cl.stats[STAT_HEALTH])*SCL;
	const s32 y = HH - 12*SCL;
	Sbar_DrawNumSmall(x, y, cl.stats[STAT_HEALTH], color);
	s32 armor = cl.items & IT_INVULNERABILITY ? 666 : cl.stats[STAT_ARMOR];
	x = WW/2 - 44*SCL - int_offs(armor)*SCL;
	Sbar_DrawNumSmall(x, y, armor, color);
	x = WW/2 - 52*SCL + int_offs(armor)*SCL;
	if (cl.items&IT_INVULNERABILITY)Draw_TransPicScaled(x,y+1,&pent8x8,SCL);
	else if (rogue) {
		if (cl.items & RIT_ARMOR3)
			Draw_TransPicScaled(x, y+1, &shield_r, SCL);
		else if (cl.items & RIT_ARMOR2)
			Draw_TransPicScaled(x, y+1, &shield_y, SCL);
		else if (cl.items & RIT_ARMOR1)
			Draw_TransPicScaled(x, y+1, &shield_g, SCL);
	} else {
		if (cl.items & IT_ARMOR3)
			Draw_TransPicScaled(x, y+1, &shield_r, SCL);
		else if (cl.items & IT_ARMOR2)
			Draw_TransPicScaled(x, y+1, &shield_y, SCL);
		else if (cl.items & IT_ARMOR1)
			Draw_TransPicScaled(x, y+1, &shield_g, SCL);
	}
	x = WW/2 - int_offs(cl.stats[STAT_AMMO])*SCL + 20*SCL;
	Sbar_DrawNumSmall(x, y, cl.stats[STAT_AMMO], color);
	x = WW/2 - int_offs(cl.stats[STAT_AMMO])*SCL + 44*SCL;
	if (rogue) { // ammo icon
		if (cl.items & RIT_SHELLS)
			Draw_TransPicScaled(x, y, &shell8x8, SCL);
		else if (cl.items & RIT_NAILS)
			Draw_TransPicScaled(x, y, &nail8x8, SCL);
		else if (cl.items & RIT_ROCKETS)
			Draw_TransPicScaled(x, y, &rocket8x8, SCL);
		else if (cl.items & RIT_CELLS)
			Draw_TransPicScaled(x, y, &lightning8x8, SCL);
		else if (cl.items & RIT_LAVA_NAILS)
			Draw_TransPicScaled(x, y, &lavanail8x8, SCL);
		else if (cl.items & RIT_PLASMA_AMMO)
			Draw_TransPicScaled(x, y, &plasma8x8, SCL);
		else if (cl.items & RIT_MULTI_ROCKETS)
			Draw_TransPicScaled(x, y, &multirocket8x8, SCL);
	} else {
		if (cl.items & IT_SHELLS)
			Draw_TransPicScaled(x, y, &shell8x8, SCL);
		else if (cl.items & IT_NAILS)
			Draw_TransPicScaled(x, y, &nail8x8, SCL);
		else if (cl.items & IT_ROCKETS)
			Draw_TransPicScaled(x, y, &rocket8x8, SCL);
		else if (cl.items & IT_CELLS)
			Draw_TransPicScaled(x, y, &lightning8x8, SCL);
	}
}

void Sbar_Draw()
{
	if (scr_con_current == HH || !(scr_hudstyle.value || sb_lines)
		|| (sb_updates >= vid.numpages && !scr_hudstyle.value)
		|| cl.intermission)
		return;
	drawlayer = lyr_sbar.value;
	if (scr_hudstyle.value >= 5 && scr_hudstyle.value <= 7) {
		Sbar_Min(scr_hudstyle.value - 5);
		drawlayer = lyr_main.value;
		return;
	}
	s32 force_vanilla = (scr_hudstyle.value == 4 && WW/SCL < 640)
		|| scr_hudstyle.value == 8;
	s32 oldhudstyle = scr_hudstyle.value;
	if (force_vanilla) scr_hudstyle.value = 0;
	if (sb_lines && WW/SCL > 320)
		Draw_TileClear(0, HH - sb_lines/SCL, WW, sb_lines/SCL);
	if (scr_hudstyle.value == 0 || scr_hudstyle.value == 4)
		Sbar_DrawBg();
	if (cl.gametype==GAME_DEATHMATCH&&(sb_showscores||scr_sidescore.value))
		Sbar_DeathmatchOverlay();
	if (sb_showscores || cl.stats[STAT_HEALTH] <= 0) {
		Draw_TransPicScaled(WW/2-160*SCL, HH-24*SCL, sb_scorebar, SCL);
		Sbar_SoloScoreboard();
		sb_updates = 0;
	} else if (scr_hudstyle.value == 9) { // EZQuake
		Sbar_DrawArmor();
		Sbar_DrawFace(); // ...no HP
	} else { // bottom three, the big numbers
		Sbar_DrawArmor();
		Sbar_DrawFace(); // ...and HP
		Sbar_DrawAmmo();
	}
	if ((sb_lines/SCL>24 || oldhudstyle) && !(oldhudstyle==4 && WW/SCL<640))
		Sbar_DrawInventory();
	if (sb_lines/SCL > 24 && cl.maxclients != 1)
		Sbar_DrawFrags();
	if (force_vanilla) scr_hudstyle.value = oldhudstyle;
	drawlayer = lyr_main.value;
}

#undef SCL
#undef WW
#undef HH
