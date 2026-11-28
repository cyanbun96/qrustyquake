FROM debian:buster

RUN sed -i 's|deb.debian.org|archive.debian.org|g' /etc/apt/sources.list; \
    apt-get update; \
    apt-get install -y wget; \
    mkdir /build                 

ARG CMAKE_VERSION=4.5.1

COPY cmake-$CMAKE_VERSION-*.tar.gz /build

RUN apt-get install -y \
            build-essential \
            pkg-config \
            libx11-dev \
            libxrandr-dev \
            libxcursor-dev \
            libxfixes-dev \
            libxi-dev \
            libxss-dev \
            libxtst-dev \
            libxext-dev \
            libxkbcommon-dev \
            libdrm-dev \
            libgbm-dev \
            libwayland-dev \
            libwayland-client0 \
            libxkbcommon-x11-dev \
            libgl1-mesa-dev \
            libgles2-mesa-dev \
            libegl1-mesa-dev \
            libdbus-1-dev \
            libibus-1.0-dev \
            libudev-dev \
            libfribidi-dev \
            libthai-dev \
            libasound2-dev \
            libpulse-dev \
            libaudio-dev \
            libjack-dev \
            libusb-1.0-0-dev \
            libmpg123-dev \
            libvorbis-dev
    
COPY build-SDL3.sh /build

RUN  cd /build && ./build-SDL3.sh $CMAKE_VERSION

WORKDIR /build
CMD AppImage/run-in-docker.sh 
