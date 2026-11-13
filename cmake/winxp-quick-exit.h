/* For C++ under cmake/winxp-i686-mingw.cmake: libstdc++ configured for the
   UCRT pulls these into std::, msvcrt.dll has neither. Nothing calls them. */
extern "C" int at_quick_exit(void (*)(void));
extern "C" void quick_exit(int) __attribute__((noreturn));
