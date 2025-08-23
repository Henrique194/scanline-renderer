#ifndef VIDEO_H
#define VIDEO_H

#include <SDL.h>

#define SCREEN_WIDTH  320
#define SCREEN_HEIGHT 240

SDL_bool VID_Init(void);
void VID_Shutdown(void);
void VID_UpdateScreen(void);
Uint8* VID_GetScreenBuffer(void);

#endif
