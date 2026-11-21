#include "quakedef.h"

static bool mouse_avail;
static f32 mouse_x;
static f32 mouse_y;
static SDL_Joystick *joystick;
static SDL_Gamepad *gamepad; // opened instead of joystick when SDL maps it
static u8 jhat; // held hat or D-pad directions, SDL_HAT_* bits
static f64 jhattime[4]; // last press or repeat of each direction
static const s32 jhatkeys[4] = {
	K_UPARROW, K_RIGHTARROW, K_DOWNARROW, K_LEFTARROW };
// Gamepads get the numbers the Linux xpad driver gives an Xbox pad, which the
// j* cvar defaults and JOY/AUX binds assume. Buttons past these follow SDL's
// order from MISC1 on, the D-pad acts as the hat.
static const SDL_GamepadButton gpbuttons[] = {
	SDL_GAMEPAD_BUTTON_SOUTH, SDL_GAMEPAD_BUTTON_EAST,
	SDL_GAMEPAD_BUTTON_WEST, SDL_GAMEPAD_BUTTON_NORTH,
	SDL_GAMEPAD_BUTTON_LEFT_SHOULDER, SDL_GAMEPAD_BUTTON_RIGHT_SHOULDER,
	SDL_GAMEPAD_BUTTON_BACK, SDL_GAMEPAD_BUTTON_START,
	SDL_GAMEPAD_BUTTON_GUIDE, SDL_GAMEPAD_BUTTON_LEFT_STICK,
	SDL_GAMEPAD_BUTTON_RIGHT_STICK };
#define GP_CORE (s32)Q_COUNTOF(gpbuttons)
#define GP_BUTTONS (GP_CORE+SDL_GAMEPAD_BUTTON_COUNT-SDL_GAMEPAD_BUTTON_MISC1)
static const s32 gpaxes[SDL_GAMEPAD_AXIS_COUNT] = { 0, 1, 3, 4, 2, 5 };

static void IN_CloseJoystick()
{
	if(gamepad){ SDL_CloseGamepad(gamepad); gamepad = 0; }
	if(joystick){ SDL_CloseJoystick(joystick); joystick = 0; }
}

void IN_InitJoystick(SDL_UNUSED cvar_t *cvar)
{
	if(joysticknum.value <= -1){
		Con_Printf("Joystick disabled\n");
		SDL_SetJoystickEventsEnabled(0);
		SDL_SetGamepadEventsEnabled(0);
		return;
	}
	s32 count = -1;
	SDL_JoystickID *joys = SDL_GetJoysticks(&count);
	if(joysticknum.value > count - 1 || joysticknum.value < 0){
		Con_Printf("Invalid joysticknum value\n");
		if(count > 0){
			Con_Printf("Detected joysticks:\n");
			for(s32 i = 0; i < count; ++i)
				Con_Printf("%d: %s\n", i,
					SDL_GetJoystickNameForID(joys[i]));
		}
		else Con_Printf("No joysticks found\n");
	}else if(count > 0){
		SDL_JoystickID id = joys[(s32)joysticknum.value];
		IN_CloseJoystick();
		if(SDL_IsGamepad(id)) gamepad = SDL_OpenGamepad(id);
		if(!gamepad) joystick = SDL_OpenJoystick(id);
		Con_Printf("Using %s: %s\n", gamepad ? "gamepad" : "joystick",
			SDL_GetJoystickNameForID(id));
		SDL_SetJoystickEventsEnabled(1);
		SDL_SetGamepadEventsEnabled(1);
	}else Con_Printf("No joysticks found\n");
	SDL_free(joys);
}

static void IN_JoyHat(u8 hat) // arrows are held while the hat is
{
	for(s32 i = 0; i < 4; i++){
		if(!((hat ^ jhat) & 1 << i)) continue;
		Key_Event(jhatkeys[i], hat & 1 << i);
		jhattime[i] = Sys_DoubleTime();
	}
	jhat = hat;
}

void IN_RemoveJoystick()
{
	IN_JoyHat(0);
	IN_CloseJoystick();
	jaxis_move_x = 0;
	jaxis_move_y = 0;
	jaxis_look_x = 0;
	jaxis_look_y = 0;
	Con_Printf("Joystick removed\n");
	IN_InitJoystick(0);
}

void IN_HatRepeat()
{
	if(key_dest == key_game) return;
	if(joysticknum.value <= -1){
		jhat = 0; // jet go of the keys
		return;
	}
	f64 t = Sys_DoubleTime();
	for(s32 i = 0; i < 4; i++){
		if(!(jhat & 1 << i) || jhattime[i] + 0.5 >= t) continue;
		Key_Event(jhatkeys[i], 0); // up first so the key stays held
		Key_Event(jhatkeys[i], 1);
		jhattime[i] += 0.05;
	}
}

static void IN_JoyButton(s32 b, bool down)
{
	if(b == jenterbutton.value) Key_Event(K_ENTER, down);
	else if(b == jescapebutton.value) Key_Event(K_ESCAPE, down);
	else if(b < K_AUX31 - K_JOY1) // AUX31 and 32 are the triggers
		Key_Event(K_JOY1 + b, down);
}

static void IN_PadButton(s32 b, bool down)
{
	static const u8 dpad[4] = {
		SDL_HAT_UP, SDL_HAT_DOWN, SDL_HAT_LEFT, SDL_HAT_RIGHT };
	if(b >= SDL_GAMEPAD_BUTTON_DPAD_UP
			&& b <= SDL_GAMEPAD_BUTTON_DPAD_RIGHT){
		u8 bit = dpad[b - SDL_GAMEPAD_BUTTON_DPAD_UP];
		IN_JoyHat(down ? jhat | bit : jhat & ~bit);
		return;
	}
	s32 i = 0;
	while(i < GP_CORE && gpbuttons[i] != b) i++;
	if(i == GP_CORE) i += b - SDL_GAMEPAD_BUTTON_MISC1;
	IN_JoyButton(i, down);
}

static void IN_JoyAxis(s32 axis, s32 value)
{
	if     (axis == jmoveaxisx.value) jaxis_move_x = value;
	else if(axis == jmoveaxisy.value) jaxis_move_y = value;
	else if(axis == jlookaxisx.value) jaxis_look_x = value;
	else if(axis == jlookaxisy.value) jaxis_look_y = value;
	else if(axis == jtrigaxis1.value){
		Key_Event(K_AUX31, value > jtriggerthresh.value);
		jaxis_trig_1 = value; }
	else if(axis == jtrigaxis2.value){
		Key_Event(K_AUX32, value > jtriggerthresh.value);
		jaxis_trig_2 = value; }
}

s32 IN_JoyButtons() // in the numbering binds and cvars use
{
	if(gamepad) return GP_BUTTONS;
	return SDL_GetNumJoystickButtons(joystick);
}

s32 IN_JoyAxes()
{
	if(gamepad) return SDL_GAMEPAD_AXIS_COUNT;
	return SDL_GetNumJoystickAxes(joystick);
}

bool IN_JoyButtonHeld(s32 b)
{
	if(!gamepad) return SDL_GetJoystickButton(joystick, b);
	return SDL_GetGamepadButton(gamepad, b < GP_CORE ? gpbuttons[b] :
		(SDL_GamepadButton)(SDL_GAMEPAD_BUTTON_MISC1 + b - GP_CORE));
}

void Sys_SendKeyEvents()
{
	IN_HatRepeat();
	SDL_Event event;
	s32 button, winW, winH;
	s32 sym, state, mod; // keep here for OpenBSD compiler
	bool usemouse = SDLWindowFlags & SDL_WINDOW_FULLSCREEN
		|| SDL_GetWindowFullscreenMode(window)
		|| _windowed_mouse.value;
	if(!usemouse || (key_dest == key_menu && ui_mouse.value))
		SDL_ShowCursor();
	else SDL_HideCursor();
	SDL_SetWindowRelativeMouseMode(window, usemouse && key_dest==key_game);
	while(SDL_PollEvent(&event)){
		static s32 buttonremap[] = {
			K_MOUSE1, K_MOUSE3, K_MOUSE2, K_MOUSE4, K_MOUSE5 };
		switch(event.type){
		case SDL_EVENT_KEY_DOWN:
		case SDL_EVENT_KEY_UP:
			sym = event.key.key;
			state = event.key.down;
			mod = SDL_GetModState() & SDL_KMOD_NUM;
			switch(sym){
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
#ifdef SDL_PLATFORM_ANDROID
			case SDLK_AC_BACK: sym = K_ESCAPE; break;
#endif
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
			// key reserved for the console, independent of location
			// on keyboard Located in the top left corner(on both
			// ANSI and ISO keyboards).
			if(event.key.scancode == SDL_SCANCODE_GRAVE){
				Key_Event('`', state);
				break;
			}
			// If we're not directly handled and still above 255
			// just force it to 0
			if(sym > 255) sym = 0;
			if(event.key.repeat){
				Key_Event(sym, !state);
				Key_Event(sym, state);
			}
			else Key_Event(sym, state);
			break;
		// Vanilla behavior: Use Mouse OFF disables mouse input entirely
		// ON grabs the mouse, kinda like SetRelativeMouseMode(SDL_TRUE)
		// Fullscreen grabs the mouse unconditionally
		case SDL_EVENT_MOUSE_MOTION:
			if(key_dest == key_game){
				mouse_x += event.motion.xrel;
				mouse_y += event.motion.yrel;
			}else if(key_dest == key_menu && ui_mouse.value){
				// in points, destRect in pixels
				f32 d = SDL_GetWindowPixelDensity(window);
				s32 mx = event.motion.x * d;
				s32 my = event.motion.y * d;
				s32 rx = mx - destRect.x;
				s32 ry = my - destRect.y;
				if(rx < 0) rx = 0;
				if(ry < 0) ry = 0;
				s32 qx = rx * vid.width / destRect.w;
				s32 qy = ry * vid.height / destRect.h;
				if(qx >= (s32)vid.width-1)
					qx = vid.width-1;
				if(qy >= (s32)vid.height-1)
					qy = vid.height-1;
				M_MouseCursor(qx, qy);
			}
			break;
		case SDL_EVENT_MOUSE_BUTTON_DOWN:
		case SDL_EVENT_MOUSE_BUTTON_UP:
			if(event.button.button < 1 ||
				event.button.button > Q_COUNTOF(buttonremap)){
				Con_Printf("Ignored input event for MOUSE%d\n",
						event.button.button);
				break;
			}
			button = buttonremap[event.button.button - 1];
		      Key_Event(button,event.type==SDL_EVENT_MOUSE_BUTTON_DOWN);
			break;
		case SDL_EVENT_MOUSE_WHEEL:
			if(event.wheel.y > 0){
				Key_Event(K_MWHEELUP, 1);
				Key_Event(K_MWHEELUP, 0);
			}
			else if(event.wheel.y < 0){
				Key_Event(K_MWHEELDOWN, 1);
				Key_Event(K_MWHEELDOWN, 0);
			}
			break;
		case SDL_EVENT_WINDOW_PIXEL_SIZE_CHANGED:
			SDLWindowFlags = SDL_GetWindowFlags(window);
			SDL_GetWindowSizeInPixels(window, &winW, &winH);
			Cvar_SetValue("realwidth", winW);
			Cvar_SetValue("realheight", winH);
			VID_CalcScreenDimensions(0);
			break;
		case SDL_EVENT_WINDOW_FOCUS_LOST: Key_ClearStates(); break;
		case SDL_EVENT_QUIT:
			CL_Disconnect();
			Host_ShutdownServer(0);
			Sys_Quit();
			break;
		case SDL_EVENT_JOYSTICK_BUTTON_UP:
		case SDL_EVENT_JOYSTICK_BUTTON_DOWN:
			if(!gamepad) IN_JoyButton(event.jbutton.button,
					event.jbutton.down);
			break;
		case SDL_EVENT_JOYSTICK_HAT_MOTION:
			if(!gamepad && !event.jhat.hat)
				IN_JoyHat(event.jhat.value);
			break;
		case SDL_EVENT_JOYSTICK_AXIS_MOTION:
			if(!gamepad) IN_JoyAxis(event.jaxis.axis,
					event.jaxis.value);
			break;
		case SDL_EVENT_GAMEPAD_BUTTON_UP:
		case SDL_EVENT_GAMEPAD_BUTTON_DOWN:
			IN_PadButton(event.gbutton.button, event.gbutton.down);
			break;
		case SDL_EVENT_GAMEPAD_AXIS_MOTION:{
			s32 v = event.gaxis.value; // triggers go 0 to 32767
			if(event.gaxis.axis >= SDL_GAMEPAD_AXIS_LEFT_TRIGGER)
				v = v * 2 - 32768;
			IN_JoyAxis(gpaxes[event.gaxis.axis], v);
			break; }
		case SDL_EVENT_JOYSTICK_ADDED: IN_InitJoystick(0); break;
		case SDL_EVENT_JOYSTICK_REMOVED: IN_RemoveJoystick(); break;
		default:
			break;
		}
	}
}

void IN_Init()
{
	if(COM_CheckParm("-nomouse")) return;
	mouse_x = mouse_y = 0.0;
	mouse_avail = 1;
	SDL_Init(SDL_INIT_GAMEPAD);
	IN_InitJoystick(0);
}

void IN_Shutdown(){ mouse_avail = 0; }

void IN_Move(usercmd_t *cmd)
{
	if(!mouse_avail) return;
	s32 jlook_x = jaxis_look_x;
	s32 jlook_y = jaxis_look_y;
	if     (jlook_x > 0 && jlook_x <  jdeadzone.value) jlook_x = 0;
	else if(jlook_x < 0 && jlook_x > -jdeadzone.value) jlook_x = 0;
	if     (jlook_y > 0 && jlook_y <  jdeadzone.value) jlook_y = 0;
	else if(jlook_y < 0 && jlook_y > -jdeadzone.value) jlook_y = 0;
	if     (jlook_x > 0) jlook_x -= jdeadzone.value;
	else if(jlook_x < 0) jlook_x += jdeadzone.value;
	if     (jlook_y > 0) jlook_y -= jdeadzone.value;
	else if(jlook_y < 0) jlook_y += jdeadzone.value;
	mouse_x += jlook_x*0.0005*jlooksens.value*host_frametime*72.0;
	mouse_y += jlook_y*0.0005*jlooksens.value*host_frametime*72.0;
	s32 jmove_x = jaxis_move_x;
	s32 jmove_y = jaxis_move_y;
	if     (jmove_x > 0 && jmove_x <  jdeadzone.value) jmove_x = 0;
	else if(jmove_x < 0 && jmove_x > -jdeadzone.value) jmove_x = 0;
	if     (jmove_y > 0 && jmove_y <  jdeadzone.value) jmove_y = 0;
	else if(jmove_y < 0 && jmove_y > -jdeadzone.value) jmove_y = 0;
	if     (jmove_x > 0) jmove_x -= jdeadzone.value;
	else if(jmove_x < 0) jmove_x += jdeadzone.value;
	if     (jmove_y > 0) jmove_y -= jdeadzone.value;
	else if(jmove_y < 0) jmove_y += jdeadzone.value;
	cmd->sidemove += jmove_x*0.005*jmovesens.value;
	cmd->forwardmove -= jmove_y*0.005*jmovesens.value;
	if(CL_InCutscene()){ // drop the look, see CL_InCutscene
		mouse_x = mouse_y = 0.0;
		return;
	}
	if((!(SDLWindowFlags & SDL_WINDOW_FULLSCREEN
			|| SDL_GetWindowFullscreenMode(window))
			&& !_windowed_mouse.value) || key_dest != key_game)
				return;
	mouse_x *= sensitivity.value;
	mouse_y *= sensitivity.value * sensitivityyscale.value;
	if((in_strafe.state & 1) || (lookstrafe.value && (in_mlook.state & 1)))
		cmd->sidemove += m_side.value * mouse_x;
	else
		cl.viewangles[YAW] -= m_yaw.value * mouse_x;
	if(in_mlook.state & 1)
		V_StopPitchDrift();
	if((in_mlook.state & 1) && !(in_strafe.state & 1)){
		cl.viewangles[PITCH] += m_pitch.value * mouse_y;
		if(cl.viewangles[PITCH] > cl_maxpitch.value)
			cl.viewangles[PITCH] = cl_maxpitch.value;
		if(cl.viewangles[PITCH] < cl_minpitch.value)
			cl.viewangles[PITCH] = cl_minpitch.value;
	}else{
		if((in_strafe.state & 1) && noclip_anglehack)
			cmd->upmove -= m_forward.value * mouse_y;
		else
			cmd->forwardmove -= m_forward.value * mouse_y;
	}
	mouse_x = mouse_y = 0.0;
}
