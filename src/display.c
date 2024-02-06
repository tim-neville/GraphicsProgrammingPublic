#include "display.h"

SDL_DisplayMode display_mode;

SDL_Window* window = NULL;
SDL_Renderer* renderer = NULL;

static uint32_t* color_buffer = NULL;
static float* z_buffer = NULL;

SDL_Texture* color_buffer_texture = NULL;
uint16_t window_width = 801;//2240;
uint16_t window_height = 601;//1400;

enum cull_method current_cull_method = 0;
enum render_method current_render_method = 0;

typedef struct {
	int grid_size;
	uint32_t color;
	int height_start;
	int height_end;
	int width_start;
	int width_end;
} ThreadData;

bool initialize_window(void)
{
	if (SDL_Init(SDL_INIT_EVERYTHING) != 0)
	{
		fprintf(stderr, "Error initializing SDL.\n");
		return false;
	}

	SDL_GetCurrentDisplayMode(0, &display_mode);
	// Create an SDL Window
	window = SDL_CreateWindow(
		NULL,
		SDL_WINDOWPOS_CENTERED,
		SDL_WINDOWPOS_CENTERED,
		window_width,
		window_height,
		SDL_WINDOW_BORDERLESS
	);

	if (!window)
	{
		fprintf(stderr, "Error creating SDL window.\n");
		return false;
	}

	// Create an SDL renderer
	renderer = SDL_CreateRenderer(window, -1, 0);
	if (!renderer)
	{
		fprintf(stderr, "Error creating SDL renderer.\n");
	}

	// Allocate the required memory in bytes to hold the color buffer and the z-buffer
	color_buffer = (uint32_t*)malloc(sizeof(uint32_t) * window_width * window_height);
	z_buffer = (float*)malloc(sizeof(float) * window_width * window_height);

	// Create an SDL texture that is used to display the color buffer
	color_buffer_texture = SDL_CreateTexture(
		renderer,
		SDL_PIXELFORMAT_ARGB8888,
		SDL_TEXTUREACCESS_STREAMING,
		window_width,
		window_height
	);

	return true;
}

void set_window_fullscreen(void)
{
	SDL_SetWindowPosition(window, 0, 0);
	SDL_SetWindowSize(window, display_mode.w, display_mode.h);
}

void set_window_borderless(void)
{
	SDL_SetWindowSize(window, window_width, window_height);
	SDL_SetWindowPosition(window, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED);
	SDL_SetWindowBordered(window, false);
	SDL_SetWindowFullscreen(window, false);
}

void set_window_bordered(void)
{
	SDL_SetWindowSize(window, window_width, window_height);
	SDL_SetWindowBordered(window, true);
	SDL_SetWindowFullscreen(window, false);
}

void draw_pixel(int x, int y, uint32_t color)
{
	if (x < 0 || x >= window_width || y < 0 || y >= window_height)
		return;

	color_buffer[(window_width * y) + x] = color;

}

void draw_pixel_zbuffered(int x, int y, vec3_t point_a, vec3_t point_b, vec3_t point_c, uint32_t color)
{
	int index = (window_width * y) + x;

	if (index < 0 || index > window_width * window_height)
		return;


	vec2_t point_p = { x,y };
	vec2_t a = vec2_from_vec3(point_a);
	vec2_t b = vec2_from_vec3(point_b);
	vec2_t c = vec2_from_vec3(point_c);

	vec3_t weights = barycentric_weights(a, b, c, point_p);

	float alpha = weights.x;
	float beta = weights.y;
	float gamma = weights.z;

	float interpolated_reciprocal_w;
	// Interpolate the value of 1/w for the current pixel (z of each contains the w i.e. the original Z before projection)
	interpolated_reciprocal_w = (1 / point_a.z) * alpha + (1 / point_b.z) * beta + (1 / point_c.z) * gamma;

	//int index = (window_width * y) + x;
	// Only draw the pixel if the depth value is less than the one previously stored in the z-buffer
	if ((interpolated_reciprocal_w) > z_buffer[index])
	{
		// look up point in texture and draw the pixel
		if (x > 0 && x < window_width && y > 0 && y < window_height)
		{
			// Update the z-buffer value with the 1/w of this current pixel
			color_buffer[index] = color;
			z_buffer[index] = interpolated_reciprocal_w;
		}
	}
}

void draw_pixel_wrapped(int x, int y, uint32_t color)
{
	color_buffer[(window_width * MOD(y, window_height)) + MOD(x, window_width)] = color;
}

void draw_pixel_in_grid(int x, int y, uint32_t color)
{
	color_buffer[(window_width * MOD(y, window_height - 1)) + MOD(x, window_width - 1)] = color;
}


int threaded_grid(void* data)
{
	ThreadData* tdata = data;

	for (int y = tdata->height_start; y < tdata->height_end; y++)
	{
		for (int x = tdata->width_start; x < tdata->width_end; x++)
		{
			if (x % tdata->grid_size == 0 || y % tdata->grid_size == 0)
			{
				draw_pixel_wrapped(x, y, tdata->color);
			}
		}
	}

	return 0;
}

void draw_grid_threaded(int grid_size)
{
	int startTicks = SDL_GetTicks();

	int linesPerThread = (window_height / num_threads);
	int remainingLines = window_height - linesPerThread * num_threads;

	for (int i = 0; i < num_threads; i++)
	{
		ThreadData data =
		{
			.grid_size = 20,
			.color = 0xFF333333,
			.height_start = i * linesPerThread,
			.height_end = data.height_start + linesPerThread,
			.width_start = i * window_width * linesPerThread,
			.width_end = data.width_start + window_width
		};

		if (i == num_threads - 1)
		{
			data.height_end += remainingLines;
		}

		threadIDs[i] = SDL_CreateThread(threaded_grid, "test", &data);
		//SDL_WaitThread(threadIDs[i], NULL);
	}

	for (int i = 0; i < num_threads; i++)
	{
		SDL_WaitThread(threadIDs[i], NULL);
	}

	int endTicks = SDL_GetTicks() - startTicks;
	DebugInt("Time Taken: ", endTicks);
}

void draw_grid(int grid_size)
{
	//int startTicks = SDL_GetTicks();

	uint32_t color = 0xA1A1A1FF;//0xFF333333;

	for (int y = 0; y < window_height; y++)
	{
		for (int x = 0; x < window_width; x++)
		{
			if (x % grid_size == 0 || y % grid_size == 0)
			{
				draw_pixel_wrapped(x, y, color);
				//draw_pixel(x, y, color);
				//color = 0xFF333333;
			}
		}
	}

	//int endTicks = SDL_GetTicks() - startTicks;
	//DebugText("Time Taken: ");
	//DebugInt(endTicks);
}

void draw_grid_world(int grid_size, vec3_t points[], int count)
{
	//draw_line(points[0].x, points[0].y, points[1].x, points[1].y, 0xFFFF00FF);

	for (int i = 0; i < count; i++)
	{
		draw_line(points[i].x, points[i].y, points[(i + 1) % (count - 1)].x, points[(i + 1) % (count - 1)].y, 0xFFFF00FF);
	}
}

void draw_dot_grid(int grid_size)
{
	uint32_t color = 0x404040FF;

	for (int y = 0; y < window_height; y += grid_size)
	{
		for (int x = 0; x < window_width; x += grid_size)
		{
			draw_pixel_wrapped(x, y, color);
		}
	}
}

void draw_rect(int posX, int posY, int width, int height, uint32_t color)
{
	int xmin = posX - width;
	int xmax = posX + width;
	int ymin = posY - height;
	int ymax = posY + height;

	for (int x = xmin; x < xmax; x++)
	{
		for (int y = ymin; y < ymax; y++)
		{
			draw_pixel(x, y, color);
		}
	}
}

void draw_line3D(int x0, int y0, float w0, int x1, int y1, float w1, uint32_t color)
{
	//Torbjorn Haugen from the comments on Exercise: Z-Buffer for filled triangles
	int delta_x = x1 - x0;
	int delta_y = y1 - y0;
	int delta_reciprocal_w = 1.f / w1 - 1.f / w0;
	if (abs(delta_x) == 0 && abs(delta_y) == 0) return;

	int side_len = abs(delta_x) > abs(delta_y) ? abs(delta_x) : abs(delta_y);
	float x_inc = delta_x / (float)side_len;
	float y_inc = delta_y / (float)side_len;
	float w_inc = delta_reciprocal_w / (float)side_len;
	float cx = x0;
	float cy = y0;
	float cw = 1.f / w0;
	for (int i = 0; i <= side_len; i++) {
		int x = roundf(cx);
		int y = roundf(cy);
		//float one_over_w = lerp(1.f/z0, 1.f/z1, i/(float)side_len );
		float one_over_w = cw;
		float interpolated_z = 1.0f - one_over_w;
		if (interpolated_z < z_buffer[y * window_height + x])
		{
			draw_pixel(x, y, color);
			z_buffer[y * window_height + x] = interpolated_z;
		}
		cx += x_inc;
		cy += y_inc;
		cw += w_inc;
	}
}

void draw_line(int x0, int y0, int x1, int y1, uint32_t color)
{
	dda_line(x0, y0, x1, y1, color);
	//bresenham_line(x0, y0, x1, y1, color);
}

void draw_line_gouraud(int x0, int y0, int x1, int y1, uint32_t light_intensity0, uint32_t light_intensity1)
{
	dda_line_gouraud(x0, y0, x1, y1, light_intensity0, light_intensity1);
	//bresenham_line(x0, y0, x1, y1, color);
}

void draw_line_phong(int x0, int y0, int x1, int y1, uint32_t vertex_color, vec3_t normal0, vec3_t normal1)
{
	dda_line_phong(x0, y0, x1, y1, vertex_color, normal0, normal1);
	//bresenham_line(x0, y0, x1, y1, color);
}

void dda_line(int x0, int y0, int x1, int y1, uint32_t color)
{
	int delta_x = (x1 - x0);
	int delta_y = (y1 - y0);

	int longest_side_length = abs(delta_x) >= abs(delta_y) ? abs(delta_x) : abs(delta_y);

	// Find how much we should increment in both x and y, each step
	float x_inc = delta_x / (float)longest_side_length;
	float y_inc = delta_y / (float)longest_side_length;

	float current_x = x0;
	float current_y = y0;

	for (int i = 0; i <= longest_side_length; i++)
	{
		draw_pixel(round(current_x), round(current_y), color);
		current_x += x_inc;
		current_y += y_inc;
	}
}

void dda_line_gouraud(int x0, int y0, int x1, int y1, uint32_t light_intensity0, uint32_t light_intensity1)
{
	int delta_x = (x1 - x0);
	int delta_y = (y1 - y0);

	int longest_side_length = abs(delta_x) >= abs(delta_y) ? abs(delta_x) : abs(delta_y);

	// Find how much we should increment in both x and y, each step
	float x_inc = delta_x / (float)longest_side_length;
	float y_inc = delta_y / (float)longest_side_length;

	float current_x = x0;
	float current_y = y0;

	for (int i = 0; i <= longest_side_length; i++)
	{
		float lerp = inverse_lerp(x0, x1, current_x);
		uint32_t color = lerp_color(light_intensity0, light_intensity1, lerp);
		if (x0 == x1 && y0 == y1)
		{
			color = light_intensity0;
		}

		draw_pixel(round(current_x), round(current_y), color);
		current_x += x_inc;
		current_y += y_inc;
	}
}

void dda_line_phong(int x0, int y0, int x1, int y1, uint32_t vertex_color, vec3_t normal0, vec3_t normal1)
{
	int delta_x = (x1 - x0);
	int delta_y = (y1 - y0);

	int longest_side_length = abs(delta_x) >= abs(delta_y) ? abs(delta_x) : abs(delta_y);

	// Find how much we should increment in both x and y, each step
	float x_inc = delta_x / (float)longest_side_length;
	float y_inc = delta_y / (float)longest_side_length;

	float current_x = x0;
	float current_y = y0;

	for (int i = 0; i <= longest_side_length; i++)
	{
		float lerp = inverse_lerp(x0, x1, current_x);
		vec3_t normal = lerp_vec3(normal0, normal1, lerp);
		if (x0 == x1 && y0 == y1)
		{
			normal = normal0;
		}

		float light_intensity_factor = -vec3_dot(normal, directional_light.direction);

		uint32_t color = light_apply_intensity(vertex_color, light_intensity_factor);

		draw_pixel(round(current_x), round(current_y), color);
		current_x += x_inc;
		current_y += y_inc;
	}
}

//https://gist.github.com/bert/1085538
void bresenham_line(int x0, int y0, int x1, int y1, uint32_t color)
{
	int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
	int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
	int err = dx + dy, e2; /* error value e_xy */

	for (;;) {  /* loop */
		draw_pixel(x0, y0, color);
		if (x0 == x1 && y0 == y1) break;
		e2 = 2 * err;
		if (e2 >= dy) { err += dy; x0 += sx; } /* e_xy+e_x > 0 */
		if (e2 <= dx) { err += dx; y0 += sy; } /* e_xy+e_y < 0 */
	}
}

void draw_circle(int cx, int cy, int radius, uint32_t color)
{
	int error = -radius;
	int x = radius;
	int y = 0;

	// The following while loop may altered to 'while (x > y)' for a
	// performance benefit, as long as a call to 'plot4points' follows
	// the body of the loop. This allows for the elimination of the
	// '(x != y') test in 'plot8points', providing a further benefit.
	//
	// For the sake of clarity, this is not shown here.
	while (x >= y)
	{
		plot8points(cx, cy, x, y, color);
		error += y;
		++y;
		error += y;
		// The following test may be implemented in assembly language in
		// most machines by testing the carry flag after adding 'y' to
		// the value of 'error' in the previous step, since 'error'
		// nominally has a negative value.
		if (error >= 0)
		{
			--x;
			error -= x;
			error -= x;
		}
	}
}

void plot8points(int cx, int cy, int x, int y, uint32_t color)
{
	plot4points(cx, cy, x, y, color);
	if (x != y) plot4points(cx, cy, y, x, color);
}

// The '(x != 0 && y != 0)' test in the last line of this function
// may be omitted for a performance benefit if the radius of the
// circle is known to be non-zero.
void plot4points(int cx, int cy, int x, int y, uint32_t color)
{
	draw_pixel(cx + x, cy + y, color);
	if (x != 0) draw_pixel(cx - x, cy + y, color);
	if (y != 0) draw_pixel(cx + x, cy - y, color);
	if (x != 0 && y != 0) draw_pixel(cx - x, cy - y, color);
}

void draw_rect_in_grid(int posX, int posY, int width, int height, uint32_t color)
{
	int xmin = posX + 1;
	int xmax = posX + width;
	int ymin = posY + 1;
	int ymax = posY + height;

	for (int x = xmin; x < xmax; x++)
	{
		for (int y = ymin; y < ymax; y++)
		{
			draw_pixel_in_grid(x, y, color);
		}
	}
}

void render_color_buffer()
{
	SDL_UpdateTexture(
		color_buffer_texture,
		NULL,
		color_buffer,
		(int)(window_width * sizeof(uint32_t))
	);
	SDL_RenderCopy(renderer, color_buffer_texture, NULL, NULL);

	SDL_RenderPresent(renderer);
}

void clear_color_buffer(uint32_t color)
{
	for (int i = 0; i < window_width * window_height; i++)
	{
		color_buffer[i] = color;
	}
}

void clear_z_buffer()
{
	for (int i = 0; i < window_width * window_height; i++)
	{
		z_buffer[i] = 0.0;
	}
}

float get_zbuffer_at(int x, int y)
{
	if (x < 0 || x >= window_width || y < 0 || y >= window_height)
		return 1.0;
	//index = abs((window_width * y) + x) % (window_width * window_height);
	return z_buffer[(window_width * y) + x];
}

void update_zbuffer_at(int x, int y, float value)
{
	if (x < 0 || x >= window_width || y < 0 || y >= window_height)
		return;

	z_buffer[(window_width * y) + x] = value;
}

void display_depth()
{
	for (int i = 0; i < window_width * window_height; i++)
	{
		float depth = 1.0 - z_buffer[i];
		float zNear = 0.1;
		float zFar = 2000.0;
		float value = (2.0 * zNear) / (zFar + zNear - (1.0/depth) * (zFar - zNear)); // Linearize Depth if Depth Value is reciprocal of w
		//float value = (2.0 * zNear) / (zFar + zNear - (1.0 / depth) * (zFar - zNear)); // Linearize Dpeth if Not reciprocal of W
		uint8_t val = (uint8_t)(value * 255.0);
		uint32_t argb = ((uint32_t)val << 24) | ((uint32_t)val << 16) | ((uint32_t)val << 8) | (uint32_t)val;
		color_buffer[i] = argb;
	}
}

void destroy_window(void)
{
	free(color_buffer);
	free(z_buffer);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
}
