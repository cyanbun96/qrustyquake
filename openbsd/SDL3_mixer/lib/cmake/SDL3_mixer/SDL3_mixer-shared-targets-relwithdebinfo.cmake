#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "SDL3_mixer::SDL3_mixer-shared" for configuration "RelWithDebInfo"
set_property(TARGET SDL3_mixer::SDL3_mixer-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(SDL3_mixer::SDL3_mixer-shared PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELWITHDEBINFO "SDL3::SDL3-shared"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libSDL3_mixer.so.0.1.0"
  IMPORTED_SONAME_RELWITHDEBINFO "libSDL3_mixer.so.0"
  )

list(APPEND _cmake_import_check_targets SDL3_mixer::SDL3_mixer-shared )
list(APPEND _cmake_import_check_files_for_SDL3_mixer::SDL3_mixer-shared "${_IMPORT_PREFIX}/lib/libSDL3_mixer.so.0.1.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
