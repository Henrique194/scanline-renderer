#ifndef INPUT_H
#define INPUT_H

#include <SDL.h>

typedef struct {
    double v[3];
} point_t;

extern double roll;
extern double pitch;
extern double yaw;

extern double currentspeed;
extern point_t currentpos;
extern double fieldofview;

extern double xscreenscale;
extern double yscreenscale;
extern double maxscale;
extern double speedscale;

extern SDL_bool quit;

void HandleEvents(void);

#endif
