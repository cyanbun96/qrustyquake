#!/bin/sh
#shellcheck disable=SC2086,SC2046

# supported platforms: linux/amd64, linux/arm64, linux/i386, linux/arm/v7
# can be passed as first parameter
PLATFORM=${1:-linux/amd64}

# multi-arch image (linux/amd64, linux/arm64, linux/i386, linux/arm/v7) hosted on Docker Hub
# it can be built and uploaded with:
# docker buildx build --progress=plain --platform=linux/amd64,linux/arm64,linux/i386,linux/arm/v7 -t <docker hub username>/qrustyquake-build-multi  --push .
IMAGE_NAME=bubblesoftapps/qrustyquake-build-multi

# if you prefer a local image, generate it locally once with:
# docker buildx build --progress=plain --platform=linux/amd64,linux/arm64,linux/i386,linux/arm/v7 -t qrustyquake-build --load .
# and uncomment following line
# this needs docker with multi-image support
# IMAGE_NAME=qrustyquake-build

if [ "$1" = "clean" ]; then
    rm -rf AppDir qrustyquake*.AppImage
    exit 0
fi

set -xe

docker run \
	--platform="$PLATFORM" \
	--rm \
	-v "$(cd .. && pwd):/build" \
	-e HOST_UID="$(id -u)" \
	-e HOST_GID="$(id -g)" \
	"$IMAGE_NAME"
