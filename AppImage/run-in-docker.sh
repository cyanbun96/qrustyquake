#!/bin/bash
#shellcheck disable=SC2086

set -xe

# System arch. When running in a linux/i386 Docker container it will return the host arch thus possibly x86_64 instead of i386
# Must use the packaging' sytem arch to get the actual proper arch of the distro
ARCH=$(uname -m)

# arch in the generated AppImage filename
GAME_ARCH=$ARCH
DEPLOYBIN_ARCH=$ARCH
DEPLOYBIN_PLUGIN_ARCH=$ARCH

case $ARCH in

    x86_64)
        # determine packaging' system arch if deb or rpm, or UNAME_ARCH as fallback
        # Debian: i386 or x86_64 
        ARCH=$(dpkg --print-architecture)
        
        if [ "$ARCH" = "i386" -o "$ARCH" = "i586" ]; then
            GAME_ARCH=x86
            DEPLOYBIN_ARCH=i386
            # you'd think it is the same than DEPLOYBIN_ARCH but nope
            DEPLOYBIN_PLUGIN_ARCH=i686
        fi
        ;;
    armv7l)
        GAME_ARCH=armhf
        DEPLOYBIN_ARCH=armhf
        DEPLOYBIN_PLUGIN_ARCH=armhf
        ;;
esac
        

GAME_BIN=qrustyquake
GAME_VERSION=$(grep "#define VERSION" src/defines.h  | cut -d ' ' -f 3)
APPDIR=AppDir
DEPLOYBIN=linuxdeploy-${DEPLOYBIN_ARCH}.AppImage
DEPLOYBIN_PLUGIN=linuxdeploy-plugin-appimage-${DEPLOYBIN_PLUGIN_ARCH}.AppImage
# for running linuxdeploy AppImage without needing FUSE
export APPIMAGE_EXTRACT_AND_RUN=1

patch_appimage_elf_header() {
    # replace bytes 8-10 of the ELF header with zeros
    # https://refspecs.linuxfoundation.org/elf/gabi4+/ch4.eheader.html#elfid
    dd if=/dev/zero of=$1 bs=1 count=3 seek=8 conv=notrunc
}

# build qrustyquake
cd src
make clean

# -fPIC and -pie are needed on aarch64 for linuxdeploy to work (otherwise ldd on compiled binary crashes)
CFLAGS="-fPIC" LDFLAGS="-pie" make -O V=1 VERBOSE=1 -j$(nproc) 

cd ../AppImage

trap 'rm -rf $APPDIR; chown -R $HOST_UID:$HOST_GID ..' EXIT

if [ ! -x $DEPLOYBIN ]; then
    wget https://github.com/linuxdeploy/linuxdeploy/releases/download/1-alpha-20250213-2/$DEPLOYBIN \
         https://github.com/linuxdeploy/linuxdeploy-plugin-appimage/releases/download/1-alpha-20250213-1/$DEPLOYBIN_PLUGIN
    chmod +x $DEPLOYBIN $DEPLOYBIN_PLUGIN
fi

if [ "$ARCH" = "aarch64" -o "$ARCH" = "armv7l" ]; then
    # need to patch ELF header of linuxdeploy binaries for them to run on qemu on aarch64
    # otherwise it fails with "cannot execute binary file: Exec format error"
    # https://github.com/AppImage/AppImageKit/issues/828
    patch_appimage_elf_header $DEPLOYBIN
    patch_appimage_elf_header $DEPLOYBIN_PLUGIN
fi

rm -rf "$APPDIR"

./$DEPLOYBIN --appdir "$APPDIR" --executable ../src/$GAME_BIN --desktop-file qrustyquake.desktop --icon-file qqicon.png --output appimage

# normalize AppImage filename, including version and user-facing arch
mv QrustyQuake-$DEPLOYBIN_PLUGIN_ARCH.AppImage qrustyquake-$GAME_VERSION-linux-$GAME_ARCH.AppImage


