#include "palette.h"

#define COLOR(red, green, blue)                                                \
    {.r = red, .g = green, .b = blue, .a = SDL_ALPHA_OPAQUE}

static const SDL_Color ega_colors[PAL_COLOR_NUM] = {
    [PAL_COLOR_BLACK] = COLOR(0x00, 0x00, 0x00),
    [PAL_COLOR_BLUE] = COLOR(0x00, 0x00, 0xa8),
    [PAL_COLOR_GREEN] = COLOR(0x00, 0xa8, 0x00),
    [PAL_COLOR_CYAN] = COLOR(0x00, 0xa8, 0xa8),
    [PAL_COLOR_RED] = COLOR(0xa8, 0x00, 0x00),
    [PAL_COLOR_MAGENTA] = COLOR(0xa8, 0x00, 0xa8),
    [PAL_COLOR_BROWN] = COLOR(0xa8, 0x54, 0x00),
    [PAL_COLOR_LIGHTGRAY] = COLOR(0xa8, 0xa8, 0xa8),
    [PAL_COLOR_DARKGRAY] = COLOR(0x54, 0x54, 0x54),
    [PAL_COLOR_LIGHTBLUE] = COLOR(0x54, 0x54, 0xfe),
    [PAL_COLOR_LIGHTGREEN] = COLOR(0x54, 0xfe, 0x54),
    [PAL_COLOR_LIGHTCYAN] = COLOR(0x54, 0xfe, 0xfe),
    [PAL_COLOR_LIGHTRED] = COLOR(0xfe, 0x54, 0x54),
    [PAL_COLOR_LIGHTMAGENTA] = COLOR(0xfe, 0x54, 0xfe),
    [PAL_COLOR_YELLOW] = COLOR(0xfe, 0xfe, 0x54),
    [PAL_COLOR_WHITE] = COLOR(0xfe, 0xfe, 0xfe),
};

void PAL_SetPalette(const SDL_Surface* buffer) {
    SDL_Palette* palette = buffer->format->palette;
    const int first_color = PAL_COLOR_BLACK;
    const int num_colors = PAL_COLOR_NUM;
    SDL_SetPaletteColors(palette, ega_colors, first_color, num_colors);
}
