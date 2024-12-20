# Qrusty Quake

Quake like you remember it.

A modernized, SDL2-based WinQuake port aimed at faithfulness to the original and easy portability.

# Features

- A "New Options" menu, can be toggled with "newoptions 0/1"

- Integer scaling

- Borderless window with -borderless parameter

- Auto-resolution fullscreen with -fullscreen_desktop

- Hardware-accelerated frame > screen rendering

   - Boosts performance massively on systems with GPUs

   - Tanks performance on machines without GPUs

   - Use the new -forceoldrender flag to disable

- High resolution support

   - Maximum tested is 16K, 2000 times bigger that 320x200
   
   - Defined by MAXHEIGHT and MAXWIDTH in r_shared.h and d_ifacea.h
   
   - Can probably be set higher for billboard gaming

- Non-square pixels for 320x200 and 640x400 modes

   - Can be forced on other modes with -stretchpixels

- General feature parity with the original WinQuake

   - "Use Mouse" option in windowed mode (also _windowed_mouse cvar)

   - Video configuration menu (mostly for show, use -width and -heigth)

- Proper UI scaling

- Advanced audio configuration

   - The default audio rate is 11025 for more muffled WinQuake sound

   - New flags -sndsamples and -sndpitch (try -sndpitch 5 or 15)

- vim-like keybinds that work in menus, enable with -vimmode flag

- Mouse sensitivity Y-axis scaling with sensitivityyscale cvar

# Planned

- Overhaul, modernization and trimming of the source code - removal of dead platforms and platform-specific code in favor of portable, properly formatted and readable code.

- (maybe) Modern mod support

- (probably not) Windows network fix

- (definitely not) CD Audio

Contributions of any kind are very welcome. If someone implements CD audio or something I'll definitely try to merge it.

# Building

Linux: cd src && make

Windows: don't.

The current windows binary was built with msys2 and is kinda cursed. The network doesn't work at all, for example.

make -f Makefile.w64 if you love pain i guess.

# Credits

This port started out as a fork of https://github.com/atsb/sdlwinquake

Which was a fork of another fork. It's forks all the way down...

Features some fixes from Quakespasm and other ports. I wasn't the one to implement them, but will try and give credit as I peruse the source further.

--CyanBun96 <3
