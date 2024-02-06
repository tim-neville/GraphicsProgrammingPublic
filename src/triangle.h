#pragma once

#include <stdint.h>
#include "vector.h"
#include "texture.h"
#include "display.h"
#include "swap.h"
#include "tims.h"


// TEST WHETHER IT'S FASTER JUST TO STORE an array for each vec3 vertices, uint32_t vertex_colors, vec3 normals, vec2 uvs
// Rather than keep pointering to the mesh Struct of Arrays vs Array of Structs?
typedef struct
{
	int a;
	int b;
	int c;
	int uv_a;
	int uv_b;
	int uv_c;
	int color_a;
	int color_b;
	int color_c;
	int normal_a;
	int normal_b;
	int normal_c;
} face_t;

typedef struct
{
	vec4_t projected_points[3];
	tex2_t texcoords[3];
	vec3_t vertex_normals[3];
	uint32_t vertex_colors[3];
} triangle_t;

void int_swap(int* a, int* b);
void fill_flat_bottom_triangle_flat(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color);
void fill_flat_top_triangle_flat(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color);
void fill_flat_bottom_triangle_gouraud(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t light_intensity_per_vertex[3]);
void fill_flat_top_triangle_gouraud(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t light_intensity_per_vertex[3]);
void fill_flat_bottom_triangle_phong(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t vertex_color, vec3_t normals[3]);
void fill_flat_top_triangle_phong(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t vertex_color, vec3_t normals[3]);
void draw_triangle(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color);
void draw_filled_triangle_flat(int x0, int y0, float w0, int x1, int y1, float w1, int x2, int y2, float w2, uint32_t color);
void draw_filled_triangle_gouraud(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t light_intensity_per_vertex[3]);
void draw_filled_triangle_phong(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t vertex_color, vec3_t normals[3]);

vec3_t barycentric_weights(vec2_t a, vec2_t b, vec2_t c, vec2_t p);

void draw_triangle_pixel(
	int x, int y, uint32_t color,
	vec4_t point_a, vec4_t point_b, vec4_t point_c
);

void draw_texel(
	int x, int y, uint32_t* texture,
	vec3_t point_a, vec3_t point_b, vec3_t point_c,
	tex2_t uv_a, tex2_t uv_b, tex2_t uv_c
);

void draw_textured_triangle(
	int x0, int y0, float w0, float u0, float v0,
	int x1, int y1, float w1, float u1, float v1,
	int x2, int y2, float w2, float u2, float v2,
	uint32_t* texture
);

float lerp_float(float a, float b, float f);
vec3_t lerp_vec3(vec3_t a, vec3_t b, float f);
float inverse_lerp(float a, float b, float value);
uint32_t lerp_color(uint32_t a, uint32_t b, float f);
uint32_t lerp_color_sRGB(uint32_t color1, uint32_t color2, float fraction);
uint32_t multiply_colors(uint32_t color1, uint32_t color2);
uint32_t add_colors(uint32_t color1, uint32_t color2);
uint32_t subtract_colors(uint32_t color1, uint32_t color2);
