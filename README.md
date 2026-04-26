# Qrusty Quake

Quake like you remember it.

A modernized, SDL3-based WinQuake port aimed at faithfulness to the original and easy portability.

# Features

- A "New Options" menu, can be toggled with "newoptions 0/1"

- Integer/nearest neighbor scaling

- Borderless window with -borderless parameter

- Auto-resolution fullscreen with -fullscreen_desktop

- Hardware-accelerated frame > screen rendering

- r_renderscale cvar that can be used to render at lower resolutions in exclusive fullscreen
  
  - In console, use "r_renderscale 2" and "vid_setmode 320 240 1" to render at 320x240 in 640x480 fullscreen

- High resolution support
  
  - Maximum tested is 16K, 2000 times bigger that 320x200
  
  - Defined by MAXHEIGHT and MAXWIDTH in r_shared.h
  
  - Can probably be set higher for billboard gaming

- Horizontal FOV scaling for modern widescreen resolutions (yaspectscale cvar)

- Non-square pixels for 320x200 and 640x400 modes and modern aspect ratios
  
  - aspectr cvar can be used for more granular adjustment
  
  - UI can be scaled independently with scr_uixscale and scr_uiyscale
  
  - Implemented through a customizeable "layer" system
  
  - Several independently scalable layers for different UI elements

- General feature parity with the original WinQuake
  
  - "Use Mouse" option in windowed mode (also _windowed_mouse cvar)
  
  - Video configuration menu 
    
    - Mostly for show, use -width and -heigth or the new menu

- Proper UI scaling
  
  - scr_uiscale 1 for no scaling, scr_uiscale 2 for 200% size etc.

- 10 HUD styles (hudstyle 0-9)
  
  - Classic, Modern, QW, Arcade, EZQuake, Minimalist and variations
  
  - Transparent UI elements for modern mods

- VGA text blurbs after shutdown (can be disabled with the "quickexit" cvar)

- Custom menu BG fading options
  
  - DOSQuake for the DOS-styled brown fade
  
  - WinQuake for the default dithered fade
  
  - 50% black
  
  - None at all (useful for graphics settings adjustment)

- Modern console features
  
  - Centerprint logging with con_logcenterprint 1

- BGM support
  
  - All the formats you can think of (i.e. whatever sdl3_mixer supports)
  
  - Optional, requires sdl3_mixer

- The default audio rate is 11025 for more muffled WinQuake sound

- vim-like keybinds that work in menus, enable with -vimmode flag

- Modern gamepad support

- Mouse sensitivity Y-axis scaling with sensitivityyscale cvar
  
  - Affects gamepads

- FPS counter, scr_showfps 1

- Unlocked FPS with host_maxfps cvar
  
  - Render-server separation for 72<FPS

- Expanded limits, Fitzquake protocol allowing for moden mod support
  
  - 2021 rerelease support, place QuakeEX.kpf in the base folder for working localization

- Custom palette support (put the files at gfx/custompalette.lmp and gfx/palette.lmp)
  
  - Also settable through worldspawn flags in custom maps

- CSQC HUD support

  - Can be configured through the "Custom HUD" menu

  - scr_sbaralpha and scr_sbarscale cvars from modern engines work ONLY on CSQC HUDs

  - scr_qchudscale to adjust the CSQC HUD size independently of mod logic

- Modern console features

  - Movable cursor (vanilla just erased characters when you pressed "left")

  - Hold Ctrl to erase/move whole words

  - Paste from clipboard with Ctrl+V

  - Ctrl+C to clear input

  - Hold Crtl/Shift when scrolling with PgUp/PgDn to scroll faster

- Autosaving/loading

  - sv_autosave, sv_autoload, and sv_autosave_interval cvars

- Software imitations of modern rendering features
  
  - Colored lighting, .lit file support
    
    - r_rgblighting 0,1 to toggle
    
    - Lit water is supported with r_litwater cvar
  
  - Translucent liquids on supported maps (r_{water,slime,lava,tele}alpha 0-1)
  
  - Translucent entities with the "alpha" tag
  
  - r_alphastyle cvar and respective menu entry for the translucency rendering variations
    
    - Mix (default) - emulates the hardware-accelerated translucency by mixing the surface colors
    
    - Dither - much less memory-hungry option, but with less gradual translusency levels
  
  - Cubemapped skyboxes ("sky filename_"), only .tga format for now
    
    - r_skyfog [0-1] blends between sky and fog colors
    
    - r_skynoise [0-1] for noise dithering
  
  - Cutout textures (transparency)
  
  - Fog (Quakespasm-like syntax "fog" command)
    
    - r_nofog 1 to disable
    
    - r_lockfog 1 to disable forced setting by custom maps, r_lockfog{d,r,g,b} to set custom
    
    - r_fognoise to adjust noise bias level: 0 - disable, 1 - noisy
    
    - r_fogfactor to adjust fog thickness: 0.5 - lighter, 1 - full
    
    - r_fogscale to adjust the distance to the fog: 0.5 - further, 2 - closer
    
    - r_fogstyle 0/1/2/3, with 3 being the default and most "modern-looking"
    
    - r_fogbrightness (0.5 for half brightness, 2 for double), independent of the fog color dictated by the map
    
    - r_fogdepthcorrection (0,1) to avoid the old-style "fog moving with camera"
  
  - Better mipmaps
    
    - Regenerated on map load with r_rebuildmips 1 (all textures) or 2 (only cutouts)
    
    - Solves random fullbright pixels on dark textures and bright fringe around cutout textures
  
  - Dithered texture filtering, Unreal-style
  
  - HL-style liquid warping shader
  
  - Configurable particles
    
    - r_particlescale for relative scaling
    
    - r_particlesize for fixed size
    
    - r_particlestyle 0 - square, 1 - circle
    
    - r_particlealpha [0-1] for translusency, affected by r_alphastyle
  
  - r_mipscale for LOD distance adjustment

- Palette customization menu

  - v_saturation, v_contrast, v_redlevel, v_greenlevel, v_bluelevel, v_vibrance, v_brightness, and v_hue cvars

- "Custom Maps" menu

  - Sorting by name, monsters, secrets, and date

  - Quick scroll with D/U or PgDn/PgUp keys

- "Mods" menu

- "quickexit" cvar for skipping exit messages

- "resurrect" command

# Planned

- An actual design document. Lots of documentation, really.

- Modernization of the rendering engine to support more demanding custom maps

- Probably not
  
  - Multithreading (would require a total rewrite with multithreading as the foundation for everything to get any performance benefits after threading overhead)
  
  - Advanced mapping/speedrunning/development features that are not directly related to gameplay (I'm not a mapper, speedrunner or developer)

Contributions of any kind are very welcome. If someone implements their favorite speedrunning cvar or something I'll definitely try to merge it. Issue reports are also greatly appreciated.

# Building

Linux: cd src && make

Use CMakeLists.txt for building for Windows natively:

cmake -B build -DCMAKE_BUILD_TYPE=Release -DSDL3_DIR=C:/sdl3/cmake/ -DSDL3_mixer_DIR=C:/sdl3_mixer/cmake/

cd build && cmake --build .

# Successful builds

x86_64 unless specified otherwise.

VM is VirtualBox unless specified otherwise.

- Arch Linux [HW] v0.8.0
  
  - The main platform that this port is developed on. The most likely one to work
  
  - Tested with gcc, clang, and tcc

- Debian 11 [VM] v0.5.1
  
  - The oldest tested distro

- FreeBSD [HW] v0.8.0
  
  - Seemingly perfect

- OpenBSD [HW] v0.8.0
  
  - Seemingly perfect

- Ubuntu [HW, MangoPi MQ Pro, RISC-V] v0.3
  
  - Works just fine at a playable framerate (20-30~ FPS)

- HaikuOS [HW] v0.6.1
  
  - Built with cmake against SDL3 and SDL3_mixer under "haiku/"
  
  - Default heapsize crashes on launch, works with "-heapsize 100000"

- Android [HW, Termux, AARCH64, clang] v0.8.0
  
  - Ran through X11 with touch controls. *unpleasant*

- macOS [HW] v0.6.2
  
  - Built without SDL3_mixer, with cmake and SDL3 from brew, on 14.6.1
  
  - Compiles with Makefile and SDL3 from brew on MacOS 15.7.2 (mixer not tested)

- Emscripten [Browsers with WebAssembly] v0.8.0
  
  - Dark magic performed by Erysdren, please don't ask me about it

  - Updated CMakeLists by Pup Luka (since commit 5b68511)

- Windows [VM, HW] v0.8.0

# Credits

This port started out as a fork of https://github.com/atsb/sdlwinquake

Which was a fork of another fork. It's forks all the way down...

It has since evolved into the wonderfully unhinged NakedWinQuake https://github.com/atsb/NakedWinQuake

Atsb has also contributed some crucial fixes to long-standing issues and some general improvements.

Some code, including VGA text blurbs, has been taken from https://github.com/erysdren/sdl3quake/

Big thanks to Erysdren for contributing quite a lot to the backend and makefiles to make compilation on exotic systems possible, implementing BGM, mipmap regeneration, BMP screenshots, and lots of other cool stuff.

Izhido heroically made the native windows build possible, and in the process brought the Windows version to the same feature level as the others.

The FitzQuake protocol implementation, both client and server, sound system, model loading, filesystem functions, cvars and a whole lot more has been pulled directly from QuakeSpasm.

A lot of code is adapted from Ironwail. Most of the netcode is a direct copy, other chunks are adapted with minor changes. QCVM/CSQC code is mostly taken from IW too. Also the autosave system has been practically copy-pasted from Ironwail.

TGA image loading code taken from MarkV, along with lots of other software rendering code. Most of it comes from ToChriS engine by Vic.

--CyanBun96 <3
