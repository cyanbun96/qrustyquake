# sdl3_mixer cmake project-config input for CMakeLists.txt script


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was SDL3_mixerConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../" ABSOLUTE)

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include(FeatureSummary)
set_package_properties(SDL3_mixer PROPERTIES
    URL "https://www.libsdl.org/projects/SDL_mixer/"
    DESCRIPTION "SDL_mixer is a sample multi-channel audio mixer library"
)

set(SDL3_mixer_FOUND                ON)

set(SDLMIXER_VENDORED              ON)

set(SDLMIXER_FLAC                  TRUE)
set(SDLMIXER_FLAC_LIBFLAC          TRUE)
set(SDLMIXER_FLAC_DRFLAC           TRUE)

set(SDLMIXER_GME                   TRUE)

set(SDLMIXER_MOD                   TRUE)
set(SDLMIXER_MOD_XMP               TRUE)
set(SDLMIXER_MOD_XMP_LITE          TRUE)

set(SDLMIXER_MP3                   TRUE)
set(SDLMIXER_MP3_DRMP3             TRUE)
set(SDLMIXER_MP3_MPG123            TRUE)

set(SDLMIXER_MIDI                  )
set(SDLMIXER_MIDI_FLUIDSYNTH       FALSE)
set(SDLMIXER_MIDI_NATIVE           )
set(SDLMIXER_MIDI_TIMIDITY         TRUE)

set(SDLMIXER_OPUS                  TRUE)

set(SDLMIXER_VORBIS                TRUE)
set(SDLMIXER_VORBIS_STB            TRUE)
set(SDLMIXER_VORBIS_TREMOR         FALSE)
set(SDLMIXER_VORBIS_VORBISFILE     TRUE)

set(SDLMIXER_WAVE                  TRUE)

set(SDLMIXER_WAVPACK               TRUE)

set(SDLMIXER_SDL3_REQUIRED_VERSION 3.3.0)

set(SDL3_mixer_SDL3_mixer-shared_FOUND FALSE)
if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/SDL3_mixer-shared-targets.cmake")
    include("${CMAKE_CURRENT_LIST_DIR}/SDL3_mixer-shared-targets.cmake")
    set(SDL3_mixer_SDL3_mixer-shared_FOUND TRUE)
endif()

set(SDL3_mixer_SDL3_mixer-static FALSE)
if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/SDL3_mixer-static-targets.cmake")

    if(SDLMIXER_VENDORED)
        if(SDLMIXER_GME)
            include(CheckLanguage)
            check_language(CXX)
            if(NOT CMAKE_CXX_COMPILER)
                message(WARNING "CXX language not enabled. Linking to SDL3_mixer::SDL3_mixer-static might fail.")
            endif()
        endif()
    else()
        set(_sdl_cmake_module_path "${CMAKE_MODULE_PATH}")
        list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

        include(CMakeFindDependencyMacro)
        include(PkgConfigHelper)

        if(NOT DEFINED CMAKE_FIND_PACKAGE_PREFER_CONFIG)
            set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)
        endif()

        if(SDLMIXER_FLAC_LIBFLAC AND NOT TARGET FLAC::FLAC)
            find_dependency(FLAC)
            # Xiph's CMake script does not handle Threads::Threads correctly:
            # https://github.com/xiph/flac/issues/820
            get_property(flac_interface_link_libraries TARGET FLAC::FLAC PROPERTY INTERFACE_LINK_LIBRARIES)
            if(Threads::Threads IN_LIST flac_interface_link_libraries)
                find_package(Threads)
            endif()
        endif()

        if(SDLMIXER_MOD_XMP AND NOT TARGET libxmp::libxmp)
            find_dependency(libxmp)
        endif()

        if(SDLMIXER_MOD_XMP_LITE AND NOT TARGET libxmp-lite::libxmp-lite)
            find_dependency(libxmp-lite)
        endif()

        if(SDLMIXER_MP3_MPG123 AND NOT TARGET MPG123::mpg123)
            find_dependency(mpg123)
        endif()

        if(SDLMIXER_MIDI_FLUIDSYNTH AND NOT TARGET FluidSynth::libfluidsynth)
            find_dependency(FluidSynth "2.2.0")
        endif()

        if(SDLMIXER_VORBIS_TREMOR AND NOT TARGET tremor::tremor)
            find_dependency(tremor)
        endif()

        if(SDLMIXER_VORBIS_VORBISFILE AND NOT TARGET Vorbis::vorbisfile)
            find_dependency(Vorbis)
        endif()

        if(SDLMIXER_OPUS AND NOT TARGET OpusFile::opusfile)
            find_dependency(OpusFile)
        endif()

        if(SDLMIXER_WAVPACK AND NOT TARGET WavPack::WavPack)
            find_dependency(wavpack)
        endif()

        set(CMAKE_MODULE_PATH "${_sdl_cmake_module_path}")
        unset(_sdl_cmake_module_path)

        if(HAIKU AND SDLMIXER_MIDI_NATIVE)
            include(CheckLanguage)
            check_language(CXX)
            if(NOT CMAKE_CXX_COMPILER)
                message(WARNING "CXX language not enabled. Linking to SDL3_mixer::SDL3_mixer-static might fail.")
            endif()
        endif()
    endif()

    include("${CMAKE_CURRENT_LIST_DIR}/SDL3_mixer-static-targets.cmake")
    set(SDL3_mixer_SDL3_mixer-static TRUE)
endif()

function(_sdl_create_target_alias_compat NEW_TARGET TARGET)
    if(CMAKE_VERSION VERSION_LESS "3.18")
        # Aliasing local targets is not supported on CMake < 3.18, so make it global.
        add_library(${NEW_TARGET} INTERFACE IMPORTED)
        set_target_properties(${NEW_TARGET} PROPERTIES INTERFACE_LINK_LIBRARIES "${TARGET}")
    else()
        add_library(${NEW_TARGET} ALIAS ${TARGET})
    endif()
endfunction()

# Make sure SDL3_mixer::SDL3_mixer always exists
if(NOT TARGET SDL3_mixer::SDL3_mixer)
    if(TARGET SDL3_mixer::SDL3_mixer-shared)
        _sdl_create_target_alias_compat(SDL3_mixer::SDL3_mixer SDL3_mixer::SDL3_mixer-shared)
    elseif(TARGET SDL3_mixer::SDL3_mixer-static)
        _sdl_create_target_alias_compat(SDL3_mixer::SDL3_mixer SDL3_mixer::SDL3_mixer-static)
    endif()
endif()

if(NOT SDL3_mixer_COMPONENTS AND NOT TARGET SDL3_mixer::SDL3_mixer-shared AND NOT TARGET SDL3_mixer::SDL3_mixer-static)
    set(SDL3_mixer_FOUND FALSE)
endif()

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was SDL3_mixerConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../" ABSOLUTE)

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################
check_required_components(SDL3_mixer)
