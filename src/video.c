#include "video.h"
#include "palette.h"

static SDL_Window* window = NULL;
static SDL_Renderer* renderer = NULL;
static SDL_Texture* texture = NULL;
static SDL_Surface* argb_buffer = NULL;
static SDL_Surface* screen_buffer = NULL;

static const SDL_PixelFormatEnum pixel_format = SDL_PIXELFORMAT_ARGB8888;


static SDL_bool VID_AllocScreenBuffer(void) {
    const Uint32 flags = 0;
    const int w = SCREEN_WIDTH;
    const int h = SCREEN_HEIGHT;
    const int depth = 8;
    const Uint32 r = 0;
    const Uint32 g = 0;
    const Uint32 b = 0;
    const Uint32 a = 0;
    screen_buffer = SDL_CreateRGBSurface(flags, w, h, depth, r, g, b, a);
    if (!screen_buffer) {
        return SDL_FALSE;
    }
    PAL_SetPalette(screen_buffer);
    SDL_FillRect(screen_buffer, NULL, PAL_COLOR_BLACK);
    return SDL_TRUE;
}

static SDL_bool VID_AllocRgbaBuffer(void) {
    const int w = SCREEN_WIDTH;
    const int h = SCREEN_HEIGHT;
    const int depth = 0;
    const int pitch = 0;
    argb_buffer = SDL_CreateRGBSurfaceWithFormatFrom(NULL, w, h, depth, pitch,
                                                     pixel_format);
    if (!argb_buffer) {
        return SDL_FALSE;
    }
    SDL_FillRect(argb_buffer, NULL, 0);
    return SDL_TRUE;
}

static SDL_bool VID_AllocTexture(void) {
    const int access = SDL_TEXTUREACCESS_STREAMING;
    const int w = SCREEN_WIDTH;
    const int h = SCREEN_HEIGHT;
    texture = SDL_CreateTexture(renderer, pixel_format, access, w, h);
    return texture != NULL;
}

static SDL_bool VID_CreateRenderer(void) {
    const int w = SCREEN_WIDTH;
    const int h = SCREEN_HEIGHT;
    const int index = -1;
    const Uint32 flags = SDL_RENDERER_TARGETTEXTURE;
    renderer = SDL_CreateRenderer(window, index, flags);
    if (!renderer) {
        return SDL_FALSE;
    }
    SDL_RenderSetLogicalSize(renderer, w, h);
    // Blank out the full screen area in case there is any junk in
    // the borders that won't otherwise be overwritten.
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);
    SDL_RenderPresent(renderer);
    return SDL_TRUE;
}

static SDL_bool VID_CreateWindow(void) {
    const char* title = "3D clipping demo";
    const int x = SDL_WINDOWPOS_CENTERED;
    const int y = SDL_WINDOWPOS_CENTERED;
    const int w = SCREEN_WIDTH;
    const int h = SCREEN_HEIGHT;
    const Uint32 flags = SDL_WINDOW_FULLSCREEN_DESKTOP;
    window = SDL_CreateWindow(title, x, y, w, h, flags);
    if (!window) {
        return SDL_FALSE;
    }
    SDL_SetRelativeMouseMode(SDL_TRUE);
    return SDL_TRUE;
}

SDL_bool VID_Init(void) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("Could not initialize video! SDL_Error: %s\n", SDL_GetError());
        return SDL_FALSE;
    }
    if (!VID_CreateWindow()) {
        const char* error = SDL_GetError();
        printf("Error creating window for video startup: %s", error);
        return SDL_FALSE;
    }
    if (!VID_CreateRenderer()) {
        const char* error = SDL_GetError();
        printf("Error creating window renderer: %s", error);
        return SDL_FALSE;
    }
    if (!VID_AllocTexture()) {
        const char* error = SDL_GetError();
        printf("Error creating renderer texture: %s", error);
        return SDL_FALSE;
    }
    if (!VID_AllocRgbaBuffer()) {
        const char* error = SDL_GetError();
        printf("Error creating RGB surface: %s", error);
        return SDL_FALSE;
    }
    if (!VID_AllocScreenBuffer()) {
        const char* error = SDL_GetError();
        printf("Error creating screen buffer: %s", error);
        return SDL_FALSE;
    }
    return SDL_TRUE;
}

void VID_Shutdown(void) {
    if (screen_buffer) {
        SDL_FreeSurface(screen_buffer);
        screen_buffer = NULL;
    }
    if (argb_buffer) {
        SDL_FreeSurface(argb_buffer);
        argb_buffer = NULL;
    }
    if (texture) {
        SDL_DestroyTexture(texture);
        texture = NULL;
    }
    if (renderer) {
        SDL_DestroyRenderer(renderer);
        renderer = NULL;
    }
    if (window) {
        SDL_DestroyWindow(window);
        window = NULL;
    }
    SDL_QuitSubSystem(SDL_INIT_VIDEO);
}

void VID_UpdateScreen(void) {
    // Update texture.
    SDL_Rect rect = {
        .x = 0,
        .y = 0,
        .w = SCREEN_WIDTH,
        .h = SCREEN_HEIGHT,
    };
    SDL_LockTexture(texture, &rect, &argb_buffer->pixels, &argb_buffer->pitch);
    SDL_LowerBlit(screen_buffer, &rect, argb_buffer, &rect);
    SDL_UnlockTexture(texture);

    // Clear the renderer's backbuffer to remove any previous contents.
    SDL_RenderClear(renderer);

    // Copy the updated texture to the backbuffer for rendering.
    SDL_RenderCopy(renderer, texture, NULL, NULL);

    // Present the backbuffer content to the screen.
    SDL_RenderPresent(renderer);
}

Uint8* VID_GetScreenBuffer(void) {
    return screen_buffer->pixels;
}
