#include "input.h"
#include "video.h"

#define MOVEMENT_SPEED     3.0
#define VMOVEMENT_SPEED    3.0
#define MAX_MOVEMENT_SPEED 30.0
#define PI                 3.141592
#define ROLL_SPEED         (PI / 20.0)
#define PITCH_SPEED        (PI / 20.0)
#define YAW_SPEED          (PI / 20.0)

double roll;
double pitch;
double yaw;

double currentspeed;
point_t currentpos;
double fieldofview;

double xscreenscale;
double yscreenscale;
double maxscale;
double speedscale = 1.0;

SDL_bool quit = SDL_FALSE;


static void HandleKeyUp(const SDL_KeyboardEvent* e) {
    switch (e->keysym.scancode) {
        case SDL_SCANCODE_MINUS:
            fieldofview *= 0.9;
            xscreenscale = SCREEN_WIDTH / fieldofview;
            yscreenscale = SCREEN_HEIGHT / fieldofview;
            maxscale = SDL_max(xscreenscale, yscreenscale);
            break;
        case SDL_SCANCODE_EQUALS:
            fieldofview *= 1.1;
            xscreenscale = SCREEN_WIDTH / fieldofview;
            yscreenscale = SCREEN_HEIGHT / fieldofview;
            maxscale = SDL_max(xscreenscale, yscreenscale);
            break;
        case SDL_SCANCODE_F:
            speedscale *= 1.1;
            break;
        case SDL_SCANCODE_S:
            speedscale *= 0.9;
            break;
        default:
            break;
    }
}

static void HandleKeyDown(const SDL_KeyboardEvent* e) {
    switch (e->keysym.scancode) {
        case SDL_SCANCODE_DOWN:
            currentspeed -= MOVEMENT_SPEED * speedscale;
            currentspeed =
                SDL_max(currentspeed, -(MAX_MOVEMENT_SPEED * speedscale));
            break;
        case SDL_SCANCODE_UP:
            currentspeed += MOVEMENT_SPEED * speedscale;
            currentspeed =
                SDL_min(currentspeed, MAX_MOVEMENT_SPEED * speedscale);
            break;
        case SDL_SCANCODE_LEFT:
            yaw -= YAW_SPEED * speedscale;
            if (yaw < 0) {
                yaw += PI * 2;
            }
            break;

        case SDL_SCANCODE_RIGHT:
            yaw += YAW_SPEED * speedscale;
            if (yaw >= (PI * 2)) {
                yaw -= PI * 2;
            }
            break;
        case SDL_SCANCODE_N:
            roll += ROLL_SPEED * speedscale;
            if (roll >= (PI * 2)) {
                roll -= PI * 2;
            }
            break;
        case SDL_SCANCODE_M:
            roll -= ROLL_SPEED * speedscale;
            if (roll < 0) {
                roll += PI * 2;
            }
            break;
        case SDL_SCANCODE_A:
            pitch -= PITCH_SPEED * speedscale;
            if (pitch < 0) {
                pitch += PI * 2;
            }
            break;
        case SDL_SCANCODE_Z:
            pitch += PITCH_SPEED * speedscale;
            if (pitch >= (PI * 2)) {
                pitch -= PI * 2;
            }
            break;
        case SDL_SCANCODE_D:
            currentpos.v[1] += VMOVEMENT_SPEED;
            break;
        case SDL_SCANCODE_C:
            currentpos.v[1] -= VMOVEMENT_SPEED;
            break;
        default:
            break;
    }
}

void HandleEvents(void) {
    SDL_Event e;
    while (SDL_PollEvent(&e)) {
        switch (e.type) {
            case SDL_QUIT:
                quit = SDL_TRUE;
                break;
            case SDL_KEYDOWN:
                HandleKeyDown(&e.key);
                break;
            case SDL_KEYUP:
                HandleKeyUp(&e.key);
                break;
            default:
                break;
        }
    }
}
