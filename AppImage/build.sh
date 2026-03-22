#!/bin/sh
#shellcheck disable=SC2086,SC2046

# multi-arch image (linux/amd64 and linux/arm64) hosted on dockerhub.com
IMAGE_NAME=bubblesoftapps/qrustyquake-build

# if you prefer a local image, generate it locally once with:
# docker build . -t qrustyquake-build
# and uncomment following line
#IMAGE_NAME=qrustyquake-build

if [ "$1" = "clean" ]; then
    rm -rf AppDir qrustyquake*.AppImage
    exit 0
fi

set -xe

# supported platforms: linux/amd64 and linux/arm64
# can be passed as first parameter
PLATFORM=${1:-linux/amd64}

docker run \
       --platform=$PLATFORM \
       --rm -it \
       -v ..:/build \
       -e HOST_UID=$(id -u) \
       -e HOST_GID=$(id -g) \
       $IMAGE_NAME
