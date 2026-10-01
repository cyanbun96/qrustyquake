// pl_mpeg implementation using newer SDL3
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <SDL3/SDL.h>

#define PL_MPEG_IMPLEMENTATION
#include "pl_mpeg.h"

// 1. Context to pass into our callbacks
typedef struct {
    SDL_Texture *video_tex;
    SDL_AudioStream *audio_stream;
    SDL_AudioDeviceID audio_device; // Store the device ID here
} VideoContext;

// 2. Video Callback: Push raw planes directly to the GPU
static void OnVideoDecode(plm_t *plm, plm_frame_t *frame, void *user) {
    VideoContext *ctx = (VideoContext *)user;
    SDL_UpdateYUVTexture(
        ctx->video_tex, NULL,
        frame->y.data, frame->y.width,
        frame->cb.data, frame->cb.width,
        frame->cr.data, frame->cr.width
    );
}

// 3. Audio Callback with Lazy Initialization
static void OnAudioDecode(plm_t *plm, plm_samples_t *samples, void *user) {
    VideoContext *ctx = (VideoContext *)user;

    // Create the stream dynamically on the first decoded frame once the sample rate is known
    if (!ctx->audio_stream) {
        int samplerate = plm_get_samplerate(plm);
        if (samplerate > 0) {
            SDL_AudioSpec spec;
            spec.freq = samplerate;
            spec.channels = 2;
            spec.format = SDL_AUDIO_F32;

            // SDL_OpenAudioDeviceStream safely handles the default playback constant
            ctx->audio_stream = SDL_OpenAudioDeviceStream(ctx->audio_device, &spec, NULL, NULL);
            if (ctx->audio_stream) {
                // Fetch the *actual* ID of the newly opened logical device and resume it
                SDL_ResumeAudioDevice(SDL_GetAudioStreamDevice(ctx->audio_stream));
            }
        }
    }

    if (ctx->audio_stream) {
        SDL_PutAudioStreamData(ctx->audio_stream, samples->interleaved, samples->count * 2 * sizeof(float));
    }
}

/* Plays MPEG-1 video, blocking main thread until finished/skipped.
   Returns 1 if playback finishes successfully, 0 if error/skipped. */
int PLM_PlayVideo(const char *filename, SDL_Renderer *renderer, SDL_AudioDeviceID audio_device) {
    plm_t *plm = plm_create_with_filename(filename);
    if (!plm) {
        SDL_Log("PLM_PlayVideo: Could not open file %s\n", filename);
        return 0;
    }

    int width = plm_get_width(plm);
    int height = plm_get_height(plm);

    SDL_SetRenderLogicalPresentation(renderer, width, height, SDL_LOGICAL_PRESENTATION_LETTERBOX);

    SDL_Texture *video_tex = SDL_CreateTexture(
        renderer,
        SDL_PIXELFORMAT_IYUV,
        SDL_TEXTUREACCESS_STREAMING,
        width, height
    );

    if (!video_tex) {
        SDL_Log("PLM_PlayVideo: Could not create SDL3 video texture.");
        plm_destroy(plm);
        return 0;
    }

    // Set up the context and assign callbacks
    VideoContext ctx = { video_tex, NULL, audio_device };
    plm_set_video_decode_callback(plm, OnVideoDecode, &ctx);
    plm_set_audio_decode_callback(plm, OnAudioDecode, &ctx);

    // Audio is enabled by default in pl_mpeg, but we must set the lead time
    plm_set_audio_lead_time(plm, 0.25);

    int playing = 1;
    int hard_quit = 0;
    Uint64 last_time = SDL_GetTicksNS();
    SDL_Event event;

    while (SDL_PollEvent(&event)) { /* DISCARD LINGERING EVENTS */ }

    while (playing) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_EVENT_QUIT) {
                playing = 0;
                hard_quit = 1;
            } else if (event.type == SDL_EVENT_KEY_DOWN) {
                if (event.key.key == SDLK_ESCAPE || event.key.key == SDLK_SPACE || event.key.key == SDLK_RETURN) {
                    playing = 0;
                }
            } else if (event.type == SDL_EVENT_MOUSE_BUTTON_DOWN) {
                playing = 0;
            }
        }

        Uint64 current_time = SDL_GetTicksNS();
        double dt = (current_time - last_time) / 1000000000.0;
        last_time = current_time;

        if (dt > 0.1) dt = 0.1;

        plm_decode(plm, dt);

        SDL_RenderClear(renderer);
        SDL_RenderTexture(renderer, video_tex, NULL, NULL);
        SDL_RenderPresent(renderer);

        if (plm_has_ended(plm)) {
            playing = 0;
        }

        SDL_DelayNS(1000000);
    }

    if (video_tex) SDL_DestroyTexture(video_tex);

    // Clean up the dynamically created audio stream safely
    if (ctx.audio_stream) {
        // DO NOT UNBIND: A stream created via OpenAudioDeviceStream manages its own device.
        SDL_DestroyAudioStream(ctx.audio_stream);
    }

    plm_destroy(plm);

    return !hard_quit;
}

/* Plays MPEG-1 video directly from an embedded memory array. */
int PLM_PlayVideoMemory(const unsigned char *data, size_t length, SDL_Renderer *renderer, SDL_AudioDeviceID audio_device) {
    plm_t *plm = plm_create_with_memory((uint8_t *)data, length, 0);
    if (!plm) {
        SDL_Log("PLM_PlayVideoMemory: Could not open video from memory.\n");
        return 0;
    }

    int width = plm_get_width(plm);
    int height = plm_get_height(plm);
    SDL_SetRenderLogicalPresentation(renderer, width, height, SDL_LOGICAL_PRESENTATION_LETTERBOX);

    SDL_Texture *video_tex = SDL_CreateTexture(
        renderer,
        SDL_PIXELFORMAT_IYUV,
        SDL_TEXTUREACCESS_STREAMING,
        width, height
    );

    if (!video_tex) {
        SDL_Log("PLM_PlayVideoMemory: Could not create SDL3 video texture.");
        plm_destroy(plm);
        return 0;
    }

    VideoContext ctx = { video_tex, NULL, audio_device };
    plm_set_video_decode_callback(plm, OnVideoDecode, &ctx);
    plm_set_audio_decode_callback(plm, OnAudioDecode, &ctx);
    plm_set_audio_lead_time(plm, 0.25);

    int playing = 1;
    int skipped = 0;
    Uint64 last_time = SDL_GetTicksNS();
    SDL_Event event;

    while (SDL_PollEvent(&event)) { /* DISCARD */ }

    while (playing) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_EVENT_QUIT) {
                playing = 0;
                skipped = 1;
            } else if (event.type == SDL_EVENT_KEY_DOWN) {
                if (event.key.key == SDLK_ESCAPE || event.key.key == SDLK_SPACE || event.key.key == SDLK_RETURN) {
                    playing = 0;
                    skipped = 1;
                }
            } else if (event.type == SDL_EVENT_MOUSE_BUTTON_DOWN) {
                playing = 0;
                skipped = 1;
            }
        }

        Uint64 current_time = SDL_GetTicksNS();
        double dt = (current_time - last_time) / 1000000000.0;
        last_time = current_time;

        if (dt > 0.1) dt = 0.1;

        plm_decode(plm, dt);

        SDL_RenderClear(renderer);
        SDL_RenderTexture(renderer, video_tex, NULL, NULL);
        SDL_RenderPresent(renderer);

        if (plm_has_ended(plm)) {
            playing = 0;
        }

        SDL_DelayNS(1000000);
    }

    if (video_tex) SDL_DestroyTexture(video_tex);

    if (ctx.audio_stream) {
        // DO NOT UNBIND: A stream created via OpenAudioDeviceStream manages its own device.
        SDL_DestroyAudioStream(ctx.audio_stream);
    }

    plm_destroy(plm);

    return !skipped;
}
