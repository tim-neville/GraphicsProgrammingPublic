#pragma once 

#include <stdint.h>
#include "vector.h"

extern enum lighting_method
{
	NO_LIGHTING = 0,
	FLAT_SHADING = 1,
	GOURAUD_SHADING = 2,
	PHONG_SHADING = 3
} current_lighting_method;
/////

typedef struct {
	vec3_t direction;
	uint32_t color;
} light_t;

extern light_t directional_light;

uint32_t light_apply_intensity(uint32_t original_color, float percentage_factor);
