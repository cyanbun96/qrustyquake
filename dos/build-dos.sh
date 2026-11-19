#!/bin/sh
# Reproducible MS-DOS build of Qrusty Quake
#
# usage: sh build-dos.sh [commit, default HEAD] [work directory, default
#        build next to this script]
#
# Builds a commit of the git checkout this script is in; commits from
# before the DOS support stop with an error. Fetches the rest at pinned
# versions and checks them, and leaves the game in <work>/dist, where it
# replaces only its own files, so an ID1 put there stays. Rerunning reuses
# the downloads in <work>/downloads, and two runs for the same commit give
# a byte-identical QRUSTYQ.EXE.
#
# Host: Linux x86_64 (the DJGPP below is a Linux build) with git, curl,
# tar, bzip2, unzip, cmake 3.22 or newer, make and sha256sum.
#
# The same steps by hand:
#
# 1. DJGPP cross compiler, GCC 12.2.0 from build-djgpp v3.4:
#    https://github.com/andrewwutw/build-djgpp/releases/download/v3.4/
#    djgpp-linux64-gcc1220.tar.bz2, sha256
#    8464f17017d6ab1b2bb2df4ed82357b5bf692e6e2b7fee37e315638f3d505f00
#    Unpack it and put djgpp/bin on PATH.
# 2. CWSDPMI r7, the DPMI host DJGPP programs load at startup:
#    https://www.delorie.com/pub/djgpp/current/v2misc/csdpmi7b.zip, sha256
#    deacda0488e1cdd7c4a9f32fab45662b34c0ed6b2d7d4d13bc07041b62004a8c
#    Its terms (CWSDPMI.DOC) ask whoever passes it on without its source
#    to tell users they have the right to the source and where it is,
#    which CWSDPMI.TXT does:
#    https://www.delorie.com/pub/djgpp/current/v2misc/csdpmi7s.zip
# 3. The Qrusty Quake source at that commit (git archive). Its DOS code is
#    all behind SDL_PLATFORM_DOS: it turns off the render threads (SDL's
#    DOS threads are cooperative) and networking (DJGPP has no sockets, so
#    single player only), stands in for fenv.h, fmin and fmax with the x87
#    control word and inline functions, makes the first hunk segment 24 MB
#    (DOSBox crashed when CWSDPMI swapped the usual 384 MB request),
#    includes SDL_main.h, whose DOS main locks memory for SDL's interrupt
#    handlers and enables the near pointers its video driver needs, and
#    prints a fatal error a second time, after SDL's switch back to text
#    mode has cleared the screen.
# 4. SDL at commit fff6448f523f01e705d771a9f0ae782982b71482 of main (no
#    release had the DOS port yet), and the SDL_mixer release the commit's
#    CMakeLists.txt names in QRUSTYQUAKE_SDL3_MIXER_VERSION, tag
#    release-<version>, without its submodules.
# 5. Configure and build from the source:
#    CFLAGS=-DSDL_DISABLE_SSE cmake -B build \
#        -DCMAKE_TOOLCHAIN_FILE=<SDL>/build-scripts/i586-pc-msdosdjgpp.cmake \
#        -DCMAKE_BUILD_TYPE=Release -DQRUSTYQUAKE_VENDOR_SDL=ON \
#        -DFETCHCONTENT_FULLY_DISCONNECTED=ON \
#        -DFETCHCONTENT_SOURCE_DIR_SDL3=<SDL> \
#        -DFETCHCONTENT_SOURCE_DIR_SDL3_MIXER=<SDL_mixer> \
#        -DBUILD_SHARED_LIBS=OFF -DSDL_SHARED=OFF -DSDL_STATIC=ON \
#        -DCMAKE_POSITION_INDEPENDENT_CODE=OFF
#    cmake --build build
#    Without SDL_DISABLE_SSE, SDL_mixer refuses to start on CPUs without
#    SSE, which is every DOS-era one and DOSBox's. SDL_mixer also turns on
#    position independent code for every target by default; DJGPP warns
#    that it ignores -fPIC, yet SDL_mixer.c still compiles differently, so
#    the script has it off. For identical builds the script also sets
#    SOURCE_DATE_EPOCH to the commit time (the exe prints its build time)
#    and maps the work directory out of file names with -ffile-prefix-map.
# 6. Rename build/qrustyquake.exe to QRUSTYQ.EXE, since plain DOS has 8.3
#    file names, and put CWSDPMI.EXE next to it. The script also adds
#    CWSDPMI.DOC, CWSDPMI.TXT and LICENSE.TXT.
#
# Tested with 1.0.0 and these DOS changes on an i5-12600KF, with the
# shareware ID1, timedemo demo1, against the original QUAKE.EXE 1.06. In
# DOSBox 0.74-3, Quake ran at 112 fps and this build at 95. In FreeDOS 1.4
# under QEMU with KVM (JEMMEX, 128 MB, Sound Blaster 16), Quake ran at
# 680 fps and this build at 920. Sound and music worked in both.
#
# Known problems in DOSBox, which will not be fixed:
# - Sound and music get choppy below about 22 fps. SDL's DOS Sound Blaster
#   driver queues only about 45 ms of audio and refills it when the game
#   yields to it, once a frame, so a real machine that slow has it too.
# - Modes 1024 and 1280 pixels wide give a black screen while the game
#   runs on; 640 and 800 wide work. The cause is unknown.
set -eu

SDL_URL=https://github.com/libsdl-org/SDL.git
SDL_COMMIT=fff6448f523f01e705d771a9f0ae782982b71482
MIX_URL=https://github.com/libsdl-org/SDL_mixer.git
DJGPP_TAR=djgpp-linux64-gcc1220.tar.bz2
DJGPP_URL=https://github.com/andrewwutw/build-djgpp/releases/download/v3.4
DJGPP_SHA=8464f17017d6ab1b2bb2df4ed82357b5bf692e6e2b7fee37e315638f3d505f00
CWS_URL=https://www.delorie.com/pub/djgpp/current/v2misc/csdpmi7b.zip
CWS_SHA=deacda0488e1cdd7c4a9f32fab45662b34c0ed6b2d7d4d13bc07041b62004a8c
CWS_SRC=https://www.delorie.com/pub/djgpp/current/v2misc/csdpmi7s.zip

die() { echo "build-dos.sh: $*" >&2; exit 1; }
for c in git curl tar bzip2 unzip cmake make sha256sum; do
	command -v $c >/dev/null || die "needs $c"
done
[ "$(uname -s) $(uname -m)" = "Linux x86_64" ] ||
	die "needs a Linux x86_64 host for the DJGPP toolchain"
here=$(cd "$(dirname "$0")" && pwd)
repo=$(git -C "$here" rev-parse --show-toplevel 2>/dev/null) ||
	die "$here is not in a git checkout"
ref=${1:-HEAD}
commit=$(git -C "$repo" rev-parse -q --verify "$ref^{commit}") ||
	die "no commit $ref in $repo"
git -C "$repo" show "$commit:src/sys_sdl.c" | grep -q SDL_PLATFORM_DOS ||
	die "$ref has no DOS support"
mixver=$(git -C "$repo" show "$commit:CMakeLists.txt" | sed -n \
	's/.*QRUSTYQUAKE_SDL3_MIXER_VERSION "\([0-9.]*\)".*/\1/p' | head -1)
[ -n "$mixver" ] || die "no SDL_mixer version in $ref's CMakeLists.txt"
work=${2:-$here/build}
mkdir -p "$work/downloads"
work=$(cd "$work" && pwd)
case $work in *[[:space:]]*) die "the work directory path has spaces";; esac
dl=$work/downloads

download() { # url sha256 file
	if [ ! -f "$3" ]; then
		echo "downloading $1"
		curl -fsSL -o "$3.part" "$1"
		mv "$3.part" "$3"
	fi
	echo "$2  $3" | sha256sum -c --quiet - || die "$3: wrong checksum"
}
checkout() { # url ref dir [commit it must be]
	[ -d "$3/.git" ] || git init -q "$3"
	echo "fetching $1 $2"
	git -C "$3" fetch -q --depth 1 "$1" "$2"
	git -C "$3" checkout -q --force FETCH_HEAD
	git -C "$3" clean -q -fdx
	[ -z "${4-}" ] || [ "$(git -C "$3" rev-parse HEAD)" = "$4" ] ||
		die "$1 $2 is not $4"
}

download "$DJGPP_URL/$DJGPP_TAR" $DJGPP_SHA "$dl/$DJGPP_TAR"
download "$CWS_URL" $CWS_SHA "$dl/csdpmi7b.zip"
rm -rf "$work/djgpp" "$work/cwsdpmi" "$work/qrustyquake"
tar -xjf "$dl/$DJGPP_TAR" -C "$work"
unzip -q "$dl/csdpmi7b.zip" bin/CWSDPMI.EXE bin/cwsdpmi.doc \
	-d "$work/cwsdpmi"
echo "taking $commit from $repo"
mkdir "$work/qrustyquake"
git -C "$repo" archive "$commit" | tar -x -C "$work/qrustyquake"
checkout $SDL_URL $SDL_COMMIT "$work/SDL" $SDL_COMMIT
checkout $MIX_URL refs/tags/release-$mixver "$work/SDL_mixer"

PATH=$work/djgpp/bin:$PATH
SOURCE_DATE_EPOCH=$(git -C "$repo" log -1 --format=%ct "$commit")
CFLAGS="-DSDL_DISABLE_SSE -ffile-prefix-map=$work=/qq-dos"
export PATH SOURCE_DATE_EPOCH CFLAGS
tc=$work/SDL/build-scripts/i586-pc-msdosdjgpp.cmake
rm -rf "$work/build"
cmake -S "$work/qrustyquake" -B "$work/build" -G "Unix Makefiles" \
	-DCMAKE_TOOLCHAIN_FILE="$tc" \
	-DCMAKE_BUILD_TYPE=Release -DQRUSTYQUAKE_VENDOR_SDL=ON \
	-DFETCHCONTENT_FULLY_DISCONNECTED=ON \
	-DFETCHCONTENT_SOURCE_DIR_SDL3="$work/SDL" \
	-DFETCHCONTENT_SOURCE_DIR_SDL3_MIXER="$work/SDL_mixer" \
	-DBUILD_SHARED_LIBS=OFF -DSDL_SHARED=OFF -DSDL_STATIC=ON \
	-DCMAKE_POSITION_INDEPENDENT_CODE=OFF >/dev/null
cmake --build "$work/build" --parallel "$(getconf _NPROCESSORS_ONLN)"

d=$work/dist
mkdir -p "$d"
cp "$work/build/qrustyquake.exe" "$d/QRUSTYQ.EXE"
cp "$work/cwsdpmi/bin/CWSDPMI.EXE" "$d/CWSDPMI.EXE"
sed 's/\r*$/\r/' "$work/cwsdpmi/bin/cwsdpmi.doc" > "$d/CWSDPMI.DOC"
sed 's/$/\r/' > "$d/CWSDPMI.TXT" <<EOF
CWSDPMI.EXE is CWSDPMI r7, unmodified, Copyright (C) 1995-2010 Charles W
Sandmann; see CWSDPMI.DOC. You have the right to receive the source code
and binary updates for CWSDPMI. Its source:
$CWS_SRC
EOF
sed 's/\r*$/\r/' "$work/qrustyquake/LICENSE" > "$d/LICENSE.TXT"
echo
echo "done: $d"
echo "commit $commit"
(cd "$d" && sha256sum QRUSTYQ.EXE)
