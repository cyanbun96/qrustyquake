# 32-bit Windows XP cross build with MinGW-w64, GCC 13+ for -mcrtdll:
# cmake -B bxp -DCMAKE_TOOLCHAIN_FILE=cmake/winxp-i686-mingw.cmake \
#	-DCMAKE_BUILD_TYPE=Release -DQRUSTYQUAKE_VENDOR_SDL=ON
set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_SYSTEM_PROCESSOR x86)
set(CMAKE_C_COMPILER i686-w64-mingw32-gcc)
set(CMAKE_CXX_COMPILER i686-w64-mingw32-g++)
set(CMAKE_RC_COMPILER i686-w64-mingw32-windres)
set(CMAKE_FIND_ROOT_PATH /usr/i686-w64-mingw32)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)
# XP has no UCRT; msvcrt-os is the system msvcrt.dll, while plain msvcrt
# is the UCRT on toolchains that default to it
set(XP_CRT "-mcrtdll=msvcrt-os")
set(CMAKE_C_FLAGS_INIT "${XP_CRT}")
# SDL always compiles its GameInput .cpp files
set(CMAKE_CXX_FLAGS_INIT
	"${XP_CRT} -include \"${CMAKE_CURRENT_LIST_DIR}/winxp-quick-exit.h\"")
set(CMAKE_EXE_LINKER_FLAGS_INIT "${XP_CRT} -static-libgcc")
set(CMAKE_SHARED_LINKER_FLAGS_INIT "${XP_CRT} -static-libgcc")
set(HAVE_GAMEINPUT_H 0) # Windows 10+ only
