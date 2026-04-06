#!/bin/sh
#shellcheck disable=SC2086

set -xe

GAME_BIN=qrustyquake
GAME_VERSION=$(grep "#define VERSION" src/defines.h  | cut -d ' ' -f 3)
ARCH=$(uname -m)
APPDIR=AppDir
DEPLOYBIN=linuxdeploy-${ARCH}.AppImage
DEPLOYBIN_PLUGIN=linuxdeploy-plugin-appimage-${ARCH}.AppImage
# for running linuxdeploy AppImage without needing FUSE
export APPIMAGE_EXTRACT_AND_RUN=1

patch_appimage_elf_header() {
    # replace bytes 8-10 of the ELF header with zeros
    # https://refspecs.linuxfoundation.org/elf/gabi4+/ch4.eheader.html#elfid
    dd if=/dev/zero of=$1 bs=1 count=3 seek=8 conv=notrunc
}

# if there is a newer version of SDL3 from OBS (Open Build Service), install it
zypper ref -r SDL3
zypper --no-refresh in -y SDL3-devel SDL3_mixer-devel

# build qrustyquake
cd src
make clean
# -fPIC and -pie are needed on aarch64 for linuxdeploy to work (otherwise ldd on compiled binary crashes)
CFLAGS="-fPIC" LDFLAGS="-pie" $(rpm -E '%make_build') RELEASE_FLAGS="$(rpm -E '%optflags')"

cd ../AppImage

trap 'chown -R $HOST_UID:$HOST_GID ..' EXIT

if [ ! -x $DEPLOYBIN ]; then
    wget https://github.com/linuxdeploy/linuxdeploy/releases/download/1-alpha-20250213-2/$DEPLOYBIN \
         https://github.com/linuxdeploy/linuxdeploy-plugin-appimage/releases/download/1-alpha-20250213-1/$DEPLOYBIN_PLUGIN
    chmod +x $DEPLOYBIN $DEPLOYBIN_PLUGIN
fi

if [ "$ARCH" = "aarch64" ]; then
    # need to patch ELF header of linuxdeploy binaries for them to run on qemu on aarch64
    # otherwise it fails with "cannot execute binary file: Exec format error"
    # https://github.com/AppImage/AppImageKit/issues/828
    patch_appimage_elf_header $DEPLOYBIN
    patch_appimage_elf_header $DEPLOYBIN_PLUGIN
fi

rm -rf "$APPDIR"

./$DEPLOYBIN --appdir "$APPDIR" --executable ../src/$GAME_BIN --desktop-file qrustyquake.desktop --icon-file qqicon.png --output appimage

# normalize AppImage filename, including version
mv QrustyQuake-$ARCH.AppImage qrustyquake-$GAME_VERSION-linux-$ARCH.AppImage


