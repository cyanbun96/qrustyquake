cmake_minimum_required(VERSION 3.22)
project(qrustyquake LANGUAGES C HOMEPAGE_URL "https://github.com/cyanbun96/qrustyquake/" VERSION 0.5.1)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Math)

option(QRUSTYQUAKE_VENDOR_SDL "Build a local vendored copy of SDL3 and SDL3_mixer" OFF)

if(QRUSTYQUAKE_VENDOR_SDL OR MSVC)
	set(QRUSTYQUAKE_SDL3_VERSION "3.4.8")
	set(QRUSTYQUAKE_SDL3_MIXER_VERSION "3.2.2")
endif()

cmake_policy(SET CMP0135 NEW)
include(FetchContent)

if(QRUSTYQUAKE_VENDOR_SDL)
	FetchContent_Declare(
		sdl3
		GIT_REPOSITORY "https://github.com/libsdl-org/SDL.git"
		GIT_TAG "release-${QRUSTYQUAKE_SDL3_VERSION}"
		EXCLUDE_FROM_ALL
	)
	FetchContent_MakeAvailable(sdl3)
	FetchContent_Declare(
		sdl3_mixer
		GIT_REPOSITORY "https://github.com/libsdl-org/SDL_mixer.git"
		GIT_TAG "release-${QRUSTYQUAKE_SDL3_MIXER_VERSION}"
		EXCLUDE_FROM_ALL
	)
	FetchContent_MakeAvailable(sdl3_mixer)
else()
	if(MSVC)
		FetchContent_Declare(
			sdl3_msvc
			URL "https://github.com/libsdl-org/SDL/releases/download/release-${QRUSTYQUAKE_SDL3_VERSION}/SDL3-devel-${QRUSTYQUAKE_SDL3_VERSION}-VC.zip"
			URL_HASH MD5=93a24a36087db120373162fc7f0ea91f
		)
		FetchContent_MakeAvailable(sdl3_msvc)
		set(SDL3_DIR "${sdl3_msvc_SOURCE_DIR}/cmake")
		FetchContent_Declare(
			sdl3_mixer_msvc
			URL "https://github.com/libsdl-org/SDL_mixer/releases/download/release-${QRUSTYQUAKE_SDL3_MIXER_VERSION}/SDL3_mixer-devel-${QRUSTYQUAKE_SDL3_MIXER_VERSION}-VC.zip"
			URL_HASH MD5=dbe6ed719adc0b09b8882a111e12373d
		)
		FetchContent_MakeAvailable(sdl3_mixer_msvc)
		set(SDL3_mixer_DIR "${sdl3_mixer_msvc_SOURCE_DIR}/cmake")
	endif()
	if(EMSCRIPTEN)
		find_package(SDL3 REQUIRED CONFIG COMPONENTS SDL3-static)
		find_package(SDL3_mixer CONFIG)
	else()
		find_package(SDL3 REQUIRED CONFIG COMPONENTS SDL3)
		find_package(SDL3_mixer REQUIRED CONFIG)
	endif()
endif()

file(GLOB QRUSTYQUAKE_HEADERS CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.h")
file(GLOB QRUSTYQUAKE_SOURCES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.c")

add_executable(qrustyquake)

if(EMSCRIPTEN)
	target_link_options(qrustyquake PRIVATE -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=1gb -sSTACK_SIZE=10mb)
	target_link_options(qrustyquake PRIVATE --preload-file id1)
	if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		target_compile_options(qrustyquake PRIVATE -gsource-map)
	endif()
endif()

if(WIN32 OR MSVC OR CMAKE_SYSTEM_NAME STREQUAL "Windows")
	list(FILTER QRUSTYQUAKE_SOURCES EXCLUDE REGEX "net_udp.c")
	target_link_libraries(qrustyquake PRIVATE ws2_32)
else()
	list(FILTER QRUSTYQUAKE_SOURCES EXCLUDE REGEX "net_wins.c")
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Haiku")
	target_link_libraries(qrustyquake PRIVATE network)
endif()

target_sources(qrustyquake PRIVATE ${QRUSTYQUAKE_SOURCES})

target_include_directories(qrustyquake PRIVATE "${PROJECT_SOURCE_DIR}/src/")
target_compile_definitions(qrustyquake PRIVATE AVAIL_SDL3MIXER)
target_link_libraries(qrustyquake PRIVATE SDL3::SDL3 SDL3_mixer::SDL3_mixer $<TARGET_NAME_IF_EXISTS:Math::Math>)

if(EMSCRIPTEN)
	file(COPY ${PROJECT_SOURCE_DIR}/emscripten/index.html DESTINATION ${PROJECT_BINARY_DIR})
endif()

if(WIN32 OR MSVC OR CMAKE_SYSTEM_NAME STREQUAL "Windows" OR QRUSTYQUAKE_VENDOR_SDL)
	add_custom_command(TARGET qrustyquake POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			$<TARGET_FILE:SDL3::SDL3>
			$<TARGET_FILE:SDL3_mixer::SDL3_mixer>
			$<TARGET_FILE_DIR:qrustyquake>
	)
endif()
