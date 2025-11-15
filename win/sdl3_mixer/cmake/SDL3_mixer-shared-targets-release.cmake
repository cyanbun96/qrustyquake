#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "SDL3_mixer::SDL3_mixer-shared" for configuration "Release"
set_property(TARGET SDL3_mixer::SDL3_mixer-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SDL3_mixer::SDL3_mixer-shared PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/SDL3_mixer.lib"
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "SDL3::SDL3-shared"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/SDL3_mixer.dll"
  )

list(APPEND _cmake_import_check_targets SDL3_mixer::SDL3_mixer-shared )
list(APPEND _cmake_import_check_files_for_SDL3_mixer::SDL3_mixer-shared "${_IMPORT_PREFIX}/lib/SDL3_mixer.lib" "${_IMPORT_PREFIX}/bin/SDL3_mixer.dll" )

# Import target "SDL3_mixer::ogg" for configuration "Release"
set_property(TARGET SDL3_mixer::ogg APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SDL3_mixer::ogg PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/ogg-0.dll"
  )

list(APPEND _cmake_import_check_targets SDL3_mixer::ogg )
list(APPEND _cmake_import_check_files_for_SDL3_mixer::ogg "${_IMPORT_PREFIX}/bin/ogg-0.dll" )

# Import target "SDL3_mixer::opus" for configuration "Release"
set_property(TARGET SDL3_mixer::opus APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SDL3_mixer::opus PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/opus-0.dll"
  )

list(APPEND _cmake_import_check_targets SDL3_mixer::opus )
list(APPEND _cmake_import_check_files_for_SDL3_mixer::opus "${_IMPORT_PREFIX}/bin/opus-0.dll" )

# Import target "SDL3_mixer::opusfile" for configuration "Release"
set_property(TARGET SDL3_mixer::opusfile APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SDL3_mixer::opusfile PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/opusfile-0.dll"
  )

list(APPEND _cmake_import_check_targets SDL3_mixer::opusfile )
list(APPEND _cmake_import_check_files_for_SDL3_mixer::opusfile "${_IMPORT_PREFIX}/bin/opusfile-0.dll" )

# Import target "SDL3_mixer::vorbis" for configuration "Release"
set_property(TARGET SDL3_mixer::vorbis APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SDL3_mixer::vorbis PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/vorbis-0.dll"
  )

list(APPEND _cmake_import_check_targets SDL3_mixer::vorbis )
list(APPEND _cmake_import_check_files_for_SDL3_mixer::vorbis "${_IMPORT_PREFIX}/bin/vorbis-0.dll" )

# Import target "SDL3_mixer::vorbisfile" for configuration "Release"
set_property(TARGET SDL3_mixer::vorbisfile APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SDL3_mixer::vorbisfile PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/vorbisfile-3.dll"
  )

list(APPEND _cmake_import_check_targets SDL3_mixer::vorbisfile )
list(APPEND _cmake_import_check_files_for_SDL3_mixer::vorbisfile "${_IMPORT_PREFIX}/bin/vorbisfile-3.dll" )

# Import target "SDL3_mixer::FLAC" for configuration "Release"
set_property(TARGET SDL3_mixer::FLAC APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SDL3_mixer::FLAC PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "SDL3_mixer::ogg"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/FLAC-8.dll"
  )

list(APPEND _cmake_import_check_targets SDL3_mixer::FLAC )
list(APPEND _cmake_import_check_files_for_SDL3_mixer::FLAC "${_IMPORT_PREFIX}/bin/FLAC-8.dll" )

# Import target "SDL3_mixer::gme_shared" for configuration "Release"
set_property(TARGET SDL3_mixer::gme_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SDL3_mixer::gme_shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/gme.dll"
  )

list(APPEND _cmake_import_check_targets SDL3_mixer::gme_shared )
list(APPEND _cmake_import_check_files_for_SDL3_mixer::gme_shared "${_IMPORT_PREFIX}/bin/gme.dll" )

# Import target "SDL3_mixer::xmp_shared" for configuration "Release"
set_property(TARGET SDL3_mixer::xmp_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SDL3_mixer::xmp_shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/libxmp.dll"
  )

list(APPEND _cmake_import_check_targets SDL3_mixer::xmp_shared )
list(APPEND _cmake_import_check_files_for_SDL3_mixer::xmp_shared "${_IMPORT_PREFIX}/bin/libxmp.dll" )

# Import target "SDL3_mixer::libmpg123" for configuration "Release"
set_property(TARGET SDL3_mixer::libmpg123 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SDL3_mixer::libmpg123 PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/mpg123-0.dll"
  )

list(APPEND _cmake_import_check_targets SDL3_mixer::libmpg123 )
list(APPEND _cmake_import_check_files_for_SDL3_mixer::libmpg123 "${_IMPORT_PREFIX}/bin/mpg123-0.dll" )

# Import target "SDL3_mixer::WavPack" for configuration "Release"
set_property(TARGET SDL3_mixer::WavPack APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SDL3_mixer::WavPack PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/wavpackdll.dll"
  )

list(APPEND _cmake_import_check_targets SDL3_mixer::WavPack )
list(APPEND _cmake_import_check_files_for_SDL3_mixer::WavPack "${_IMPORT_PREFIX}/bin/wavpackdll.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
