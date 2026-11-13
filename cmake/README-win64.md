# 64-bit Windows build

`win64-mingw.cmake` cross-compiles a 64-bit Windows build on Linux with
MinGW-w64, for Windows 10 and 11. Tested on Arch Linux (MinGW-w64 GCC 16.2,
mingw-w64-crt 14), running under Wine 11.21 but not yet on real Windows.
For Windows XP, use `winxp-i686-mingw.cmake` instead; see `README-winxp.md`.

The binaries use the toolchain's default C runtime. On Arch that is the
UCRT, which Windows 10 and 11 include. The MinGW runtime is linked
statically, so there are no MinGW DLLs to ship.

## Requirements

- MinGW-w64 for x86_64 with both gcc and g++. On Arch: `mingw-w64-gcc`.
- CMake 3.22 or newer.
- git and network access: SDL3 and SDL3_mixer are cloned from GitHub
  during configure.

The file expects the MinGW sysroot in `/usr/x86_64-w64-mingw32`. If your
distribution puts it elsewhere, change `CMAKE_FIND_ROOT_PATH` in the file.

## Build

From the repository root:

    cmake -B bwin64 -DCMAKE_TOOLCHAIN_FILE=cmake/win64-mingw.cmake \
        -DCMAKE_BUILD_TYPE=Release -DQRUSTYQUAKE_VENDOR_SDL=ON
    cmake --build bwin64 --parallel

This leaves `qrustyquake.exe`, `SDL3.dll` and `SDL3_mixer.dll` in `bwin64/`.
