#!/bin/bash

set -xe

SDL_VERSION=3.4.20
SDL_MIXER_VERSION=3.2.4

MAKE_OPTS="-j$(nproc) -O VERBOSE=1 V=1"

### Install cmake

CMAKE_VERSION=$1

# system arch: x86_64, aarch64, armv7l
ARCH=$(uname -m)

# https://wiki.debian.org/SupportedArchitectures
# distro packaging arch: amd64, arm64, i386, armhf
DPKG_ARCH=$(dpkg --print-architecture)

case $DPKG_ARCH in
    amd64|arm64)

        # install pre-built cmake binary hosted on Kitware
        
        CMAKE_INSTALL_SCRIPT=cmake-$CMAKE_VERSION-linux-$ARCH.sh
        wget https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/$CMAKE_INSTALL_SCRIPT
        sh ./$CMAKE_INSTALL_SCRIPT --skip-license --prefix=/usr
        rm $CMAKE_INSTALL_SCRIPT
        ;;

    i386|armhf)

        # install local pre-built binary
        CMAKE_INSTALL_ARCHIVE=cmake-$CMAKE_VERSION-$DPKG_ARCH.tar.gz
        tar zxvf $CMAKE_INSTALL_ARCHIVE -C /
        ;;

    *)
        echo "Unsupported architecure"
        exit 1
        ;;
esac

rm cmake-$CMAKE_VERSION-*.tar.gz

### Compile SDL3

wget https://github.com/libsdl-org/SDL/releases/download/release-$SDL_VERSION/SDL3-$SDL_VERSION.tar.gz
tar zxvf SDL3-$SDL_VERSION.tar.gz
rm SDL3-$SDL_VERSION.tar.gz
pushd SDL3-$SDL_VERSION

# Debian Buster ships wayland-client 1.16 while SDL3 expects >= 1.18
sed -i 's/wayland-client>=1.18/wayland-client>=1.16/' cmake/sdlchecks.cmake

mkdir build
cd build

cmake .. \
      -DCMAKE_INSTALL_PREFIX=/usr \
      -DCMAKE_BUILD_TYPE=Release \
      -DSDL_STATIC=OFF \
      -DSDL_TEST_LIBRARY=OFF\
      -DSDL_INSTALL_DOCS=OFF
make $MAKE_OPTS
make install
popd
rm -rf SDL3-$SDL_VERSION

### Compile SDL Mixer

wget https://github.com/libsdl-org/SDL_mixer/archive/refs/tags/release-$SDL_MIXER_VERSION.tar.gz
tar zxvf release-$SDL_MIXER_VERSION.tar.gz
rm release-$SDL_MIXER_VERSION.tar.gz
pushd SDL_mixer-release-$SDL_MIXER_VERSION
mkdir build
cd build
# enable only mp3 and Ogg Vorbis
cmake .. \
      -DCMAKE_INSTALL_PREFIX=/usr \
      -DCMAKE_BUILD_TYPE=Release \
      -DSDLMIXER_INSTALL_MAN=OFF \
      -DSDLMIXER_EXAMPLES_INSTALL=OFF \
      -DSDLMIXER_TESTS_INSTALL=OFF \
      -DSDLMIXER_VORBIS_VORBISFILE=ON \
      -DSDLMIXER_MP3_MPG123=ON \
      -DSDLMIXER_MIDI_TIMIDITY=OFF \
      -DSDLMIXER_AIFF=OFF \
      -DSDLMIXER_AU=OFF \
      -DSDLMIXER_VOC=OFF \
      -DSDLMIXER_WAVE=OFF \
      -DSDLMIXER_VORBIS_STB=OFF \
      -DSDLMIXER_VORBIS_TREMOR=OFF \
      -DSDLMIXER_FLAC_LIBFLAC=OFF \
      -DSDLMIXER_FLAC_DRFLAC=OFF \
      -DSDLMIXER_MP3_DRMP3=OFF
make $MAKE_OPTS
make install
popd
rm -rf SDL_mixer-release-$SDL_MIXER_VERSION
      
