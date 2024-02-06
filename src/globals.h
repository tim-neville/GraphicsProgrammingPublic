#pragma once

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <SDL.h>
#include <SDL_thread.h>

extern int num_threads;
extern SDL_Thread** threadIDs;

#define MOD(a,b) ((((a)%(b))+(b))%(b))

typedef uint32_t color_t;
