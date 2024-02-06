#pragma once

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <SDL.h>
#include <SDL_thread.h>
#include <time.h>
#include "globals.h"

extern const float deg_to_rad;
extern clock_t timer_start;

void DebugAddress(void* pointer);

void DebugText(const char* text);

void DebugInt(const char* text, int val);

void DebugFloat(char* text, float value);

void DebugError(const char* text);

void StartTimer();

void StopTimer(const char* text);

float Remap_Float(float value, float low1, float high1, float low2, float high2);

uint32_t depth_to_greyscale(float depth);
