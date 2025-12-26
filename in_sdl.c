#include "quakedef.h"

static bool mouse_avail;
static f32 mouse_x;
static f32 mouse_y;
static s32 buttonremap[] = { K_MOUSE1, K_MOUSE3, K_MOUSE2, K_MOUSE4, K_MOUSE5 };

void IN_InitJoystick(SDL_UNUSED cvar_t *cvar)
{
	s32 count = -1;
	SDL_JoystickID *joys = SDL_GetJoysticks(&count);
	if (joysticknum.value > count - 1 || joysticknum.value < 0) {
		Con_Printf("Invalid joysticknum value\n");
		if (count > 0) {
			Con_Printf("Detected joysticks:\n");
			for (s32 i = 0; i < count; ++i)
				Con_Printf("%d: %s\n", i,
					SDL_GetJoystickNameForID(joys[i]));
		}
		else Con_Printf("No joysticks found\n");
	} else if (count > 0) {
		Con_Printf("Using joystick: %s\n",
			SDL_GetJoystickNameForID(joys[(s32)joysticknum.value]));
		if (joystick) { SDL_CloseJoystick(joystick); joystick = 0; }
		joystick = SDL_OpenJoystick(joys[(s32)joysticknum.value]);
		SDL_SetJoystickEventsEnabled(1);
	} else Con_Printf("No joysticks found\n");
}

void IN_RemoveJoystick()
{
	if (joystick) { SDL_CloseJoystick(joystick); joystick = 0; }
	jaxis_move_x = 0;
	jaxis_move_y = 0;
	jaxis_look_x = 0;
	jaxis_look_y = 0;
	Con_Printf("Joystick removed\n");
	IN_InitJoystick(0);
}

void Sys_SendKeyEvents()
{
	SDL_Event event;
	s32 button, winW, winH;
	s32 sym, state, mod; // keep here for OpenBSD compiler
	while (SDL_PollEvent(&event)) {
		switch (event.type) {
		case SDL_EVENT_KEY_DOWN:
		case SDL_EVENT_KEY_UP:
			sym = event.key.key;
			state = event.key.down;
			mod = SDL_GetModState() & SDL_KMOD_NUM;
			switch (sym) {
			case SDLK_DELETE: sym = K_DEL; break;
			case SDLK_BACKSPACE: sym = K_BACKSPACE; break;
			case SDLK_F1: sym = K_F1; break;
			case SDLK_F2: sym = K_F2; break;
			case SDLK_F3: sym = K_F3; break;
			case SDLK_F4: sym = K_F4; break;
			case SDLK_F5: sym = K_F5; break;
			case SDLK_F6: sym = K_F6; break;
			case SDLK_F7: sym = K_F7; break;
			case SDLK_F8: sym = K_F8; break;
			case SDLK_F9: sym = K_F9; break;
			case SDLK_F10: sym = K_F10; break;
			case SDLK_F11: sym = K_F11; break;
			case SDLK_F12: sym = K_F12; break;
			case SDLK_PAUSE: sym = K_PAUSE; break;
			// CyanBun96: vim-like keybinds that work in menus
			case SDLK_H: sym = vimmode ? K_LEFTARROW : 'h'; break;
			case SDLK_J: sym = vimmode ? K_DOWNARROW : 'j'; break;
			case SDLK_K: sym = vimmode ? K_UPARROW : 'k'; break;
			case SDLK_L: sym = vimmode ? K_RIGHTARROW : 'l'; break;
			case SDLK_UP: sym = K_UPARROW; break;
			case SDLK_DOWN: sym = K_DOWNARROW; break;
			case SDLK_RIGHT: sym = K_RIGHTARROW; break;
			case SDLK_LEFT: sym = K_LEFTARROW; break;
			case SDLK_INSERT: sym = K_INS; break;
			case SDLK_HOME: sym = K_HOME; break;
			case SDLK_END: sym = K_END; break;
			case SDLK_PAGEUP: sym = K_PGUP; break;
			case SDLK_PAGEDOWN: sym = K_PGDN; break;
			case SDLK_RSHIFT:
			case SDLK_LSHIFT: sym = K_SHIFT; break;
			case SDLK_RCTRL:
			case SDLK_LCTRL: sym = K_CTRL; break;
			case SDLK_RALT:
			case SDLK_LALT: sym = K_ALT; break;
			case SDLK_KP_0: sym = mod ? K_INS:SDLK_0;break;
			case SDLK_KP_1: sym = mod ? K_END:SDLK_1;break;
			case SDLK_KP_2: sym = mod ? K_DOWNARROW:SDLK_2;break;
			case SDLK_KP_3: sym = mod ? K_PGDN:SDLK_3;break;
			case SDLK_KP_4: sym = mod ? K_LEFTARROW:SDLK_4;break;
			case SDLK_KP_5: sym = SDLK_5; break;
			case SDLK_KP_6: sym = mod ? K_RIGHTARROW:SDLK_6;break;
			case SDLK_KP_7: sym = mod ? K_HOME:SDLK_7;break;
			case SDLK_KP_8: sym = mod ? K_UPARROW:SDLK_8;break;
			case SDLK_KP_9: sym = mod ? K_PGUP:SDLK_9;break;
			case SDLK_KP_PERIOD: sym = mod ?K_DEL:SDLK_PERIOD;break;
			case SDLK_KP_DIVIDE: sym = SDLK_SLASH; break;
			case SDLK_KP_MULTIPLY: sym = SDLK_ASTERISK; break;
			case SDLK_KP_MINUS: sym = SDLK_MINUS; break;
			case SDLK_KP_PLUS: sym = SDLK_PLUS; break;
			case SDLK_KP_ENTER: sym = SDLK_RETURN; break;
			case SDLK_KP_EQUALS: sym = SDLK_EQUALS; break;
			}
			// If we're not directly handled and still above 255
			// just force it to 0
			if (sym > 255) sym = 0;
			Key_Event(sym, state);
			break;
		// Vanilla behavior: Use Mouse OFF disables mouse input entirely
		// ON grabs the mouse, kinda like SetRelativeMouseMode(SDL_TRUE)
		// Fullscreen grabs the mouse unconditionally
		case SDL_EVENT_MOUSE_MOTION:
			mouse_x += event.motion.xrel;
			mouse_y += event.motion.yrel;
			break;
		case SDL_EVENT_MOUSE_BUTTON_DOWN:
		case SDL_EVENT_MOUSE_BUTTON_UP:
			if (event.button.button < 1 ||
				event.button.button > Q_COUNTOF(buttonremap)) {
				Con_Printf("Ignored input event for MOUSE%d\n",
						event.button.button);
				break;
			}
			button = buttonremap[event.button.button - 1];
			Key_Event(button, event.type==SDL_EVENT_MOUSE_BUTTON_DOWN);
			break;
		case SDL_EVENT_MOUSE_WHEEL:
			if (event.wheel.y > 0) {
				Key_Event(K_MWHEELUP, 1);
				Key_Event(K_MWHEELUP, 0);
			}
			else if (event.wheel.y < 0) {
				Key_Event(K_MWHEELDOWN, 1);
				Key_Event(K_MWHEELDOWN, 0);
			}
			break;
		case SDL_EVENT_WINDOW_PIXEL_SIZE_CHANGED:
			SDLWindowFlags = SDL_GetWindowFlags(window);
			SDL_GetWindowSize(window, &winW, &winH);
			Cvar_SetValue("realwidth", winW);
			Cvar_SetValue("realheight", winH);
			VID_CalcScreenDimensions(0);
			break;
		case SDL_EVENT_QUIT:
			CL_Disconnect();
			Host_ShutdownServer(0);
			Sys_Quit();
			break;
		case SDL_EVENT_JOYSTICK_BUTTON_UP:
		case SDL_EVENT_JOYSTICK_BUTTON_DOWN:
			if(event.jbutton.button == jenterbutton.value)
				Key_Event(K_ENTER, event.jbutton.down);
			else if(event.jbutton.button == jescapebutton.value)
				Key_Event(K_ESCAPE, event.jbutton.down);
			else Key_Event(K_JOY1+event.jbutton.button,
					event.jbutton.down);
			break;
		case SDL_EVENT_JOYSTICK_HAT_MOTION:
			Key_Event(K_UPARROW, event.jhat.value & 1);
			Key_Event(K_RIGHTARROW, event.jhat.value & 2);
			Key_Event(K_DOWNARROW, event.jhat.value & 4);
			Key_Event(K_LEFTARROW, event.jhat.value & 8);
			break;
		case SDL_EVENT_JOYSTICK_AXIS_MOTION:
			if      (event.jaxis.axis == jmoveaxisx.value)
				jaxis_move_x = event.jaxis.value;
			else if (event.jaxis.axis == jmoveaxisy.value)
				jaxis_move_y = event.jaxis.value;
			else if (event.jaxis.axis == jlookaxisx.value)
				jaxis_look_x = event.jaxis.value;
			else if (event.jaxis.axis == jlookaxisy.value)
				jaxis_look_y = event.jaxis.value;
			else if (event.jaxis.axis == jtrigaxis1.value) {
				Key_Event(K_AUX31, event.jaxis.value >
						jtriggerthresh.value);
				jaxis_trig_1 = event.jaxis.value; }
			else if (event.jaxis.axis == jtrigaxis2.value) {
				Key_Event(K_AUX32, event.jaxis.value >
						jtriggerthresh.value);
				jaxis_trig_2 = event.jaxis.value; }
			break;
		case SDL_EVENT_JOYSTICK_ADDED: IN_InitJoystick(0); break;
		case SDL_EVENT_JOYSTICK_REMOVED: IN_RemoveJoystick(); break;
		default:
			break;
		}
	}
}

void IN_Init()
{
	if (COM_CheckParm("-nomouse")) return;
	mouse_x = mouse_y = 0.0;
	mouse_avail = 1;
	SDL_Init(SDL_INIT_JOYSTICK);
	IN_InitJoystick(0);
}

void IN_Shutdown() { mouse_avail = 0; }

void IN_Move(usercmd_t *cmd)
{
	if (!mouse_avail) return;
	s32 jlook_x = jaxis_look_x;
	s32 jlook_y = jaxis_look_y;
	if      (jlook_x > 0 && jlook_x <  jdeadzone.value) jlook_x = 0;
	else if (jlook_x < 0 && jlook_x > -jdeadzone.value) jlook_x = 0;
	if      (jlook_y > 0 && jlook_y <  jdeadzone.value) jlook_y = 0;
	else if (jlook_y < 0 && jlook_y > -jdeadzone.value) jlook_y = 0;
	if      (jlook_x > 0) jlook_x -= jdeadzone.value;
	else if (jlook_x < 0) jlook_x += jdeadzone.value;
	if      (jlook_y > 0) jlook_y -= jdeadzone.value;
	else if (jlook_y < 0) jlook_y += jdeadzone.value;
	mouse_x += jlook_x*0.0005*jlooksens.value;
	mouse_y += jlook_y*0.0005*jlooksens.value;
	s32 jmove_x = jaxis_move_x;
	s32 jmove_y = jaxis_move_y;
	if      (jmove_x > 0 && jmove_x <  jdeadzone.value) jmove_x = 0;
	else if (jmove_x < 0 && jmove_x > -jdeadzone.value) jmove_x = 0;
	if      (jmove_y > 0 && jmove_y <  jdeadzone.value) jmove_y = 0;
	else if (jmove_y < 0 && jmove_y > -jdeadzone.value) jmove_y = 0;
	if      (jmove_x > 0) jmove_x -= jdeadzone.value;
	else if (jmove_x < 0) jmove_x += jdeadzone.value;
	if      (jmove_y > 0) jmove_y -= jdeadzone.value;
	else if (jmove_y < 0) jmove_y += jdeadzone.value;
	cmd->sidemove += jmove_x*0.005*jmovesens.value;
	cmd->forwardmove -= jmove_y*0.005*jmovesens.value;
	if ((!(SDLWindowFlags & SDL_WINDOW_FULLSCREEN
			|| SDL_GetWindowFullscreenMode(window))
			&& !_windowed_mouse.value)
				|| key_dest != key_game) {
		SDL_SetWindowRelativeMouseMode(window, 0);
		SDL_ShowCursor();
		return;
	}
	SDL_HideCursor();
	SDL_SetWindowRelativeMouseMode(window, 1);
	mouse_x *= sensitivity.value;
	mouse_y *= sensitivity.value * sensitivityyscale.value;
	if ((in_strafe.state & 1) || (lookstrafe.value && (in_mlook.state & 1)))
		cmd->sidemove += m_side.value * mouse_x;
	else
		cl.viewangles[YAW] -= m_yaw.value * mouse_x;
	if (in_mlook.state & 1)
		V_StopPitchDrift();
	if ((in_mlook.state & 1) && !(in_strafe.state & 1)) {
		cl.viewangles[PITCH] += m_pitch.value * mouse_y;
		if (cl.viewangles[PITCH] > 80)
			cl.viewangles[PITCH] = 80;
		if (cl.viewangles[PITCH] < -70)
			cl.viewangles[PITCH] = -70;
	} else {
		if ((in_strafe.state & 1) && noclip_anglehack)
			cmd->upmove -= m_forward.value * mouse_y;
		else
			cmd->forwardmove -= m_forward.value * mouse_y;
	}
	mouse_x = mouse_y = 0.0;
}
