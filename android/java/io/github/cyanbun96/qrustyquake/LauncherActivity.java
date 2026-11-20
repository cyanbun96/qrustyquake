package io.github.cyanbun96.qrustyquake;

import android.Manifest;
import android.app.Activity;
import android.content.ActivityNotFoundException;
import android.content.Intent;
import android.content.SharedPreferences;
import android.content.pm.PackageManager;
import android.content.res.ColorStateList;
import android.net.Uri;
import android.os.Build;
import android.os.Bundle;
import android.os.Environment;
import android.provider.DocumentsContract;
import android.provider.Settings;
import android.text.Editable;
import android.text.TextWatcher;
import android.view.KeyEvent;
import android.view.View;
import android.view.inputmethod.EditorInfo;
import android.view.inputmethod.InputMethodManager;
import android.widget.EditText;
import android.widget.TextView;
import android.widget.Toast;
import java.io.File;
import java.util.ArrayList;

// The screen before the game: the game folder, what the game will find in
// it, and the command line. Both fields are kept for the next start, Enter
// in the options field starts the game.
public class LauncherActivity extends Activity {
	// in shared storage, where file managers and USB transfers reach it
	static final File defdir = new File(
		Environment.getExternalStorageDirectory(), "qrustyquake");
	private static final String storage =
		"com.android.externalstorage.documents";
	private SharedPreferences prefs;
	private EditText dir, args;
	private TextView status;
	private ColorStateList statuscolor;
	private View grant;

	@Override
	protected void onCreate(Bundle state)
	{
		super.onCreate(state);
		if(!isTaskRoot()){ finish(); return; } // the game is under it
		setContentView(R.layout.launcher);
		prefs = getPreferences(MODE_PRIVATE);
		dir = findViewById(R.id.basedir);
		args = findViewById(R.id.args);
		status = findViewById(R.id.status);
		statuscolor = status.getTextColors();
		grant = findViewById(R.id.grant);
		if(state == null){
			String d = defdir.getPath();
			dir.setText(prefs.getString("basedir", d));
			args.setText(prefs.getString("args", ""));
			args.setSelection(args.length());
		}
		dir.addTextChangedListener(new TextWatcher(){
			public void beforeTextChanged(CharSequence s,
				int start, int count, int after){}
			public void onTextChanged(CharSequence s,
				int start, int before, int count){}
			public void afterTextChanged(Editable s){ check(); }
		});
		grant.setOnClickListener(v -> askAccess());
		findViewById(R.id.browse).setOnClickListener(v -> browse());
		findViewById(R.id.play).setOnClickListener(v -> play());
		args.setOnEditorActionListener((v, id, ev) -> {
			// Go on a soft keyboard; Enter on a real one comes
			// on press and release, the game starts after both
			if(id == EditorInfo.IME_ACTION_GO || ev != null &&
			ev.getAction() == KeyEvent.ACTION_UP) play();
			return true;
		});
	}

	@Override
	protected void onResume()
	{ // back from the access settings or the permission dialog
		super.onResume();
		check();
	}

	private boolean access()
	{ // the game reads shared storage with plain fopen
		if(Build.VERSION.SDK_INT >= 30)
			return Environment.isExternalStorageManager();
		return Build.VERSION.SDK_INT < 23 || checkSelfPermission(
			Manifest.permission.WRITE_EXTERNAL_STORAGE)
			== PackageManager.PERMISSION_GRANTED;
	}

	private void askAccess()
	{
		if(Build.VERSION.SDK_INT >= 30){
			try{ startActivity(new Intent(Settings
			.ACTION_MANAGE_APP_ALL_FILES_ACCESS_PERMISSION,
				Uri.parse("package:" + getPackageName())));
			}catch(ActivityNotFoundException e){
				toast("Allow QrustyQuake \"All files access\""
					+ " in the system settings");
			}
		}else if(Build.VERSION.SDK_INT >= 23)
			requestPermissions(new String[]{ Manifest
				.permission.WRITE_EXTERNAL_STORAGE }, 0);
	}

	private void check()
	{ // what the game would only report as an error
		boolean ok = access();
		grant.setVisibility(ok ? View.GONE : View.VISIBLE);
		if(!ok) say(false, "No file access, the game can't read "
			+ "its files");
		else if(new File(folder(), "id1/pak0.pak").isFile())
			say(true, "Found id1/pak0.pak");
		else say(false, "No id1/pak0.pak in this folder. Copy the "
			+ "pak files from Quake's id1 to id1 here");
	}

	private void say(boolean ok, String s)
	{
		status.setText(s);
		if(ok) status.setTextColor(statuscolor);
		else status.setTextColor(0xffffb74d); // amber
	}

	private String folder()
	{
		String d = dir.getText().toString().trim();
		return d.isEmpty() ? defdir.getPath() : d;
	}

	private void browse()
	{ // the system folder picker, opened at the current folder
		Intent i = new Intent(Intent.ACTION_OPEN_DOCUMENT_TREE);
		String root = Environment.getExternalStorageDirectory()
			.getPath() + "/", d = dir.getText().toString();
		if(Build.VERSION.SDK_INT >= 26 && d.startsWith(root))
			i.putExtra(DocumentsContract.EXTRA_INITIAL_URI,
				DocumentsContract.buildDocumentUri(storage,
				"primary:" + d.substring(root.length())));
		try{ startActivityForResult(i, 0); }
		catch(ActivityNotFoundException e){
			toast("No folder picker, type the path in");
		}
	}

	@Override
	protected void onActivityResult(int req, int res, Intent data)
	{
		if(res != RESULT_OK || data == null) return;
		String d = path(data.getData());
		if(d != null) dir.setText(d);
		else toast("Pick a folder on the device's storage");
	}

	private static String path(Uri tree)
	{ // the picker names a folder "primary:dir" or "1234-ABCD:dir" (an SD
	  // card), the game wants a path for fopen
		if(tree == null || !storage.equals(tree.getAuthority()))
			return null;
		String id = DocumentsContract.getTreeDocumentId(tree);
		int c = id.indexOf(':');
		if(c < 0) return null;
		String vol = id.substring(0, c);
		File root = vol.equals("primary") ?
			Environment.getExternalStorageDirectory() :
			new File("/storage", vol);
		return new File(root, id.substring(c + 1)).getPath();
	}

	private void play()
	{
		String d = folder(), a = args.getText().toString();
		// not apply(): if a game ran in this process before, the
		// next one ends it at once (SDLActivity.onCreate) and the
		// write in the background would be lost
		prefs.edit().putString("basedir", d).putString("args", a)
			.commit();
		new File(d, "id1").mkdirs(); // shows where the paks go
		ArrayList<String> argv = new ArrayList<>();
		argv.add("-basedir");
		argv.add(d);
		split(a, argv);
		// SDL takes a soft keyboard still up when its window comes for
		// its own and brings it back over the game
		((InputMethodManager)getSystemService(INPUT_METHOD_SERVICE))
			.hideSoftInputFromWindow(args.getWindowToken(), 0);
		startActivity(new Intent(this, QrustyQuakeActivity.class)
			.putExtra("argv", argv.toArray(new String[0])));
		finish();
	}

	private static void split(String s, ArrayList<String> argv)
	{ // as a shell would: spaces between arguments, quotes keep them in
		StringBuilder b = new StringBuilder();
		boolean quoted = false, arg = false;
		for(char c : s.toCharArray()){
			if(c == '"'){ quoted = !quoted; arg = true; }
			else if(quoted || !Character.isWhitespace(c)){
				b.append(c);
				arg = true;
			}else if(arg){
				argv.add(b.toString());
				b.setLength(0);
				arg = false;
			}
		}
		if(arg) argv.add(b.toString());
	}

	private void toast(String s)
	{ Toast.makeText(this, s, Toast.LENGTH_LONG).show(); }
}
