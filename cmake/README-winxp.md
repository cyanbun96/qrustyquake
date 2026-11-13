# Windows XP build

`winxp-i686-mingw.cmake` cross-compiles a 32-bit build that runs on Windows
XP, using MinGW-w64 on Linux. Tested on Arch Linux (MinGW-w64 GCC 16.2,
mingw-w64-crt 14) with Windows XP Professional SP3 in VirtualBox.

## Requirements

- MinGW-w64 for i686 with both gcc and g++, GCC 13 or newer
  (`i686-w64-mingw32-gcc --version`). On Arch: `mingw-w64-gcc`.
- CMake 3.22 or newer.
- git and network access: SDL3 and SDL3_mixer are cloned from GitHub
  during configure.

The file expects the MinGW sysroot in `/usr/i686-w64-mingw32`. If your
distribution puts it elsewhere, change `CMAKE_FIND_ROOT_PATH` in the file.

## Build

From the repository root:

    cmake -B bxp -DCMAKE_TOOLCHAIN_FILE=cmake/winxp-i686-mingw.cmake \
        -DCMAKE_BUILD_TYPE=Release -DQRUSTYQUAKE_VENDOR_SDL=ON
    cmake --build bxp --parallel

This leaves `qrustyquake.exe`, `SDL3.dll` and `SDL3_mixer.dll` in `bxp/`.

## Notes

- The default code is plain i686 with x87 floating point, so it runs on
  any Pentium Pro or later. For faster SSE2 code (Pentium 4, Athlon 64 and
  up), set the `CFLAGS` environment variable on the first configure:
  `CFLAGS=-msse2 cmake -B bxp ...`. Don't pass `-DCMAKE_C_FLAGS` or
  `-DCMAKE_CXX_FLAGS`: they replace the toolchain's flags instead of
  adding to them, and the result no longer builds for XP.
- SDL updates can break XP support. If XP reports "The procedure entry
  point ... could not be located", the build uses a function that XP's
  system DLLs don't have.
