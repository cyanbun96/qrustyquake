package io.github.cyanbun96.qrustyquake;

import android.Manifest;
import android.content.ActivityNotFoundException;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.net.Uri;
import android.os.Build;
import android.os.Bundle;
import android.os.Environment;
import android.provider.Settings;
import android.view.InputDevice;
import android.view.KeyEvent;
import java.io.File;
import org.libsdl.app.SDLActivity;

public class QrustyQuakeActivity extends SDLActivity {
	// in shared storage, where file managers and USB transfers reach it
	private static final File basedir = new File(
		Environment.getExternalStorageDirectory(), "qrustyquake");

	@Override
	protected void onCreate(Bundle state)
	{ // SDL starts the game once this window has focus, after the grant
		super.onCreate(state);
		if(Build.VERSION.SDK_INT >= 30){
			if(Environment.isExternalStorageManager()) return;
			try{ startActivity(new Intent(Settings
			.ACTION_MANAGE_APP_ALL_FILES_ACCESS_PERMISSION,
				Uri.parse("package:" + getPackageName())));
			}catch(ActivityNotFoundException e){}
		}else if(Build.VERSION.SDK_INT >= 23 && checkSelfPermission(
				Manifest.permission.WRITE_EXTERNAL_STORAGE)
				!= PackageManager.PERMISSION_GRANTED)
			requestPermissions(new String[]{
				Manifest.permission.WRITE_EXTERNAL_STORAGE}, 0);
	}

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
	{
		new File(basedir, "id1").mkdirs(); // shows where the paks go
		return new String[]{ "-basedir", basedir.getPath() };
	}
}
