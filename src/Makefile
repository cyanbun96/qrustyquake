PKGCONFIG ?= pkg-config
CC = gcc
CFLAGS = -c -Wall -Wextra -D_DEFAULT_SOURCE $(shell $(PKGCONFIG) --cflags sdl3)
LDFLAGS = -lm $(shell $(PKGCONFIG) --libs sdl3)
SDL3_MIXER_CFLAGS := $(shell $(PKGCONFIG) --cflags sdl3-mixer 2>/dev/null)
SDL3_MIXER_LIBS := $(shell $(PKGCONFIG) --libs sdl3-mixer 2>/dev/null)
TARGET = qrustyquake
SRCS = cdaudio.c chase.c cl_demo.c cl_input.c cl_main.c cl_parse.c \
	cl_tent.c cmd.c common.c console.c crc.c \
	cvar.c d_edge.c \
	d_part.c d_polyse.c d_scan.c d_sky.c d_sprite.c \
	d_surf.c d_fog.c d_common.c host.c \
	host_cmd.c in_sdl.c keys.c menu.c \
	model.c net_bsd.c net_dgrm.c net_loop.c net_main.c \
	cvarlist.c net_udp.c pr_cmds.c pr_edict.c \
	pr_exec.c r_aclip.c r_alias.c r_bsp.c r_draw.c \
	r_edge.c r_efrag.c r_light.c r_main.c r_misc.c \
	globals.c r_part.c r_sky.c r_sprite.c r_surf.c \
	sbar.c screen.c snd_dma.c snd_mem.c snd_mix.c \
	snd_sdl.c sv_main.c sv_move.c sv_phys.c sv_user.c \
	sys_sdl.c vid_sdl.c view.c wad.c world.c \
	zone.c rgbtoi.c miniz.c vgatext.c image.c

OBJS = $(SRCS:.c=.o)

BUILD_TYPE ?= release

ifeq ($(strip $(SDL3_MIXER_LIBS)),)
$(warning SDL3_Mixer not found, disabling music support)
else
CFLAGS += $(SDL3_MIXER_CFLAGS) -DAVAIL_SDL3MIXER
LDFLAGS += $(SDL3_MIXER_LIBS)
endif

# technically for Haiku we should strip -lm since it is not needed there, but
# they have a dummy library in its place so... shrug?
ifeq ($(shell uname), Haiku)
LDFLAGS+=-lnetwork
endif

ifeq ($(CC),tcc)
CFLAGS += -DSDL_DISABLE_IMMINTRIN_H
endif

DEBUG_FLAGS = -g -O0 -DDEBUG
RELEASE_FLAGS = -O3 -flto -g

ifeq ($(BUILD_TYPE), debug)
	CFLAGS += $(DEBUG_FLAGS)
else
	CFLAGS += $(RELEASE_FLAGS)
endif

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@

debug:
	$(MAKE) BUILD_TYPE=debug

release:
	$(MAKE) BUILD_TYPE=release

clean:
	rm -f $(OBJS) $(TARGET)
