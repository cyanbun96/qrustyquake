cmake_minimum_required(VERSION 3.22)
project(qrustyquake LANGUAGES C HOMEPAGE_URL "https://github.com/cyanbun96/qrustyquake/" VERSION 0.5.1)

if(EMSCRIPTEN)
	find_package(SDL3 REQUIRED CONFIG REQUIRED COMPONENTS SDL3-static)
	find_package(SDL3_mixer CONFIG COMPONENTS SDL3_mixer-static)
else()
	find_package(SDL3 REQUIRED CONFIG REQUIRED COMPONENTS SDL3-shared)
	find_package(SDL3_mixer CONFIG COMPONENTS SDL3_mixer-shared)
endif()

file(GLOB QRUSTYQUAKE_HEADERS "${PROJECT_SOURCE_DIR}/src/*.h")
file(GLOB QRUSTYQUAKE_SOURCES "${PROJECT_SOURCE_DIR}/src/*.c")

add_executable(qrustyquake)

if(EMSCRIPTEN)
	target_compile_options(qrustyquake PRIVATE -sUSE_SDL=3 -sUSE_SDL_MIXER=3)
	target_link_options(qrustyquake PRIVATE -sUSE_SDL=3 -sUSE_SDL_MIXER=3)
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
target_link_libraries(qrustyquake PRIVATE SDL3::SDL3)
if(SDL3_mixer_FOUND)
	target_compile_definitions(qrustyquake PRIVATE AVAIL_SDL3MIXER)
	target_link_libraries(qrustyquake PRIVATE SDL3_mixer::SDL3_mixer)
endif()

find_library(libmath m)
if(libmath)
	target_link_libraries(qrustyquake PRIVATE ${libmath})
endif()

if(EMSCRIPTEN)
	file(COPY ${PROJECT_SOURCE_DIR}/emscripten/index.html DESTINATION ${PROJECT_BINARY_DIR})
endif()
