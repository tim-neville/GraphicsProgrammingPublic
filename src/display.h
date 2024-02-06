#pragma once

#include "globals.h"
#include "light.h"
#include "tims.h"
#include "triangle.h"
#include "vector.h"
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <SDL.h>
#include <SDL_thread.h>

#define FPS 60
#define FRAME_TARGET_TIME (1000 / FPS)

///// Not currently using where as Gustavo is
extern enum cull_method
{
	CULL_NONE = 0,
	CULL_BACKFACE = 1
} current_cull_method;

extern enum render_method
{
	RENDER_WIRE = 0,
	RENDER_WIRE_VERTEX = 1,
	RENDER_FILL_TRIANGLE = 2,
	RENDER_FILL_TRIANGLE_WIRE = 3,
	RENDER_TEXTURED = 4,
	RENDER_TEXTURED_WIRE = 5
} current_render_method;

extern SDL_Window* window;
extern SDL_Renderer* renderer;

static uint32_t* color_buffer;
static float* z_buffer;
extern SDL_Texture* color_buffer_texture;

extern uint16_t window_width;
extern uint16_t window_height;

extern bool show_depth;

bool initialize_window(void);

void set_window_fullscreen(void);

void set_window_borderless(void);
void set_window_bordered(void);

void draw_pixel(int x, int y, uint32_t color);
void draw_pixel_zbuffered(int x, int y, vec3_t point_a, vec3_t point_b, vec3_t point_c, uint32_t color);

void draw_pixel_wrapped(int x, int y, uint32_t color);
void draw_pixel_in_grid(int x, int y, uint32_t color);
int threaded_grid(void* data);
void draw_grid_threaded(int grid_size);

void draw_grid_threaded(int gridSize);
void draw_grid(int gridSize);
void draw_grid_world(int grid_size, vec3_t points[], int count);
void draw_dot_grid(int gridSize);
void draw_rect(int posX, int posY, int width, int height, uint32_t color);

void draw_line3D(int x0, int y0, float w0, int x1, int y1, float w1, uint32_t color);
void draw_line(int x0, int y0, int x1, int y1, uint32_t color);
void draw_line_gouraud(int x0, int y0, int x1, int y1, uint32_t color0, uint32_t color1);
void draw_line_phong(int x0, int y0, int x1, int y1, uint32_t vertex_color, vec3_t normal0, vec3_t normal1);
void dda_line(int x0, int y0, int x1, int y1, uint32_t color);
void dda_line_gouraud(int x0, int y0, int x1, int y1, uint32_t light_intensity0, uint32_t light_intensity1);
void dda_line_phong(int x0, int y0, int x1, int y1, uint32_t vertex_color, vec3_t normal0, vec3_t normal1);
void bresenham_line(int x0, int y0, int x1, int y1, uint32_t color);
void draw_circle(int cx, int cy, int radius, uint32_t color);

void plot8points(int cx, int cy, int x, int y, uint32_t color);

void plot4points(int cx, int cy, int x, int y, uint32_t color);
void draw_rect_in_grid(int posX, int posY, int width, int height, uint32_t color);

void render_color_buffer();

void clear_color_buffer(uint32_t color);
void clear_z_buffer();

float get_zbuffer_at(int x, int y);
void update_zbuffer_at(int x, int y, float value);
float get_zbuffer_at(int x, int y);
void update_zbuffer_at(int x, int y, float value);

void display_depth();

void destroy_window(void);
