#!/bin/bash
# QrustyQuake - SteamOS Build Process for Arch Linux systems
# Do not MOVE script without changing variables as needed!
# Last Tested June 3rd, 2026

VERSION="0.7.3"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ----------------------------------------------------------------
printf "=== QrustyQuake SteamOS Build Manager ===\n"
printf "v$VERSION - Written by Pup Luka\n\n"

if [[ "$1" == "-setup" ]]; then
	# Podman Dependency verification.
	DEPS=(podman binutils glibc)
	INSTALL_QUEUE=()
	echo "Setting up system dependencies..."
	for dep in "${DEPS[@]}"; do
		if ! pacman -Qi "$dep" > /dev/null 2>&1; then
			echo " - $dep is missing!"
			INSTALL_QUEUE+=("$dep")
		fi
	done

	if [ ${#INSTALL_QUEUE[@]} -ne 0 ]; then
		echo "Installing missing package: ${INSTALL_QUEUE[*]}"
		sudo pacman -S --needed --noconfirm "${INSTALL_QUEUE[@]}"
	else
		echo "All dependencies installed and verified!"
	fi
	echo "Setup completed! Please re-run \"$0\" with no flags!"
	exit 0
fi
# ----------------------------------------------------------------
echo 'Spinning up boot container...'
# Use an Ubuntu 22.04 Podman container to isolate the build environment.
# Compiling against an older LTS release forces the engine to link with an older GLIBC,
# preventing "version not found" errors and ensuring maximum compatibility on SteamOS.
podman run --rm -v "$PROJECT_DIR:/workspace" -w /workspace ubuntu:22.04 /bin/bash -c "
	echo '1. Updating container packages & lists...'
	apt-get update -qq

	echo '2. Installing core build dependencies...'
	apt-get install -y build-essential patchelf git zip cmake pkg-config libasound2-dev \
	libpulse-dev libwayland-dev libx11-dev libxext-dev libxrandr-dev libxcursor-dev \
	libxi-dev libxss-dev libxtst-dev libxkbcommon-dev libdecor-0-dev libdbus-1-dev \
	libibus-1.0-dev libgl1-mesa-dev libegl1-mesa-dev libdrm-dev libgbm-dev wget -qq

	echo '3. Cloning and building SDL3...'
	cd /tmp
	git clone --depth 1 https://github.com/libsdl-org/SDL.git
	cd SDL
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
	cd build
	make -j$(nproc)
	make install

	echo '4. Cloning and building SDL3_mixer...'
	cd /tmp
	git clone --depth 1 https://github.com/libsdl-org/SDL_mixer.git
	cd SDL_mixer
	git submodule update --init --recursive
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
	cd build
	make -j$(nproc)
	make install

	echo '5. Compiling QrustyQuake engine...'
	cd /workspace
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
	cd build
	make -j$(nproc)

	cd /workspace
	echo '6. Staging release directory...'
	mkdir -p dist dist/lib dist/id1
	cp /usr/local/lib/libSDL3.so.0 ./dist || true
	cp /usr/local/lib/libSDL3_mixer.so.0 ./dist || true
	cp /workspace/build/qrustyquake ./dist || true

	echo '7. Patching executable to find contained libraries...'
	patchelf --set-rpath '\$ORIGIN' ./dist/qrustyquake
	readelf -d ./dist/qrustyquake | grep -E '(RUNPATH|RPATH)'

	echo '8. Bundling Quake Shareware Data...'
	wget -O ./dist/id1/pak0.pak https://github.com/puppluka/qrustyquake-puppy/raw/refs/heads/devkit-test/id1/pak0.pak
	wget -O ./dist/id1/LICENSE  https://raw.githubusercontent.com/puppluka/qrustyquake-puppy/refs/heads/devkit-test/id1/LICENSE

	echo '9. Creating distribution ZIP...'
	zip -r -9 -q qrustyquake-steamos.zip dist/**

	echo 'Build process complete! Exiting container...'
"
# ----------------------------------------------------------------
echo '================================================'
echo 'Success! Your SteamOS executable and libraries are ready.'
