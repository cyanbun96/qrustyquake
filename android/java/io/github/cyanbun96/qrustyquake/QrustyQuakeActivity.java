package io.github.cyanbun96.qrustyquake;

import android.view.InputDevice;
import android.view.KeyEvent;
import org.libsdl.app.SDLActivity;

public class QrustyQuakeActivity extends SDLActivity {
	@Override
	public boolean dispatchKeyEvent(KeyEvent event)
	{ // SDL takes a keyboard that also has a joystick part (the AULA F75)
	  // for a gamepad, and its arrows, D-pad keycodes, for gamepad buttons
		InputDevice dev = event.getDevice();
		int k = event.getKeyCode();
		if(mBrokenLibraries || dev == null || dev.getKeyboardType()
			!= InputDevice.KEYBOARD_TYPE_ALPHABETIC
		|| k < KeyEvent.KEYCODE_DPAD_UP || k > KeyEvent.KEYCODE_DPAD_RIGHT)
			return super.dispatchKeyEvent(event);
		if(event.getAction() == KeyEvent.ACTION_DOWN) onNativeKeyDown(k);
		else if(event.getAction() == KeyEvent.ACTION_UP) onNativeKeyUp(k);
		return true;
	}

	@Override
	protected String[] getLibraries()
	{ return new String[]{ "SDL3", "SDL3_mixer", "qrustyquake" }; }

	@Override
	protected String[] getArguments()
	{ // the command line from LauncherActivity
		String[] a = getIntent().getStringArrayExtra("argv");
		return a != null ? a : new String[]{
			"-basedir", LauncherActivity.defdir.getPath() };
	}
}
