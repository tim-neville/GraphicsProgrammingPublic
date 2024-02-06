#include "triangle.h"


void fill_flat_bottom_triangle_flat(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color)
{
	// Find the two slopes from the two triangle legs
	float inverse_slope_1 = (float)(x1 - x0) / (y1 - y0);
	float inverse_slope_2 = (float)(x2 - x0) / (y2 - y0);

	// Start x_start and x_end from the top vertex (x0, y0)
	float x_start = x0;
	float x_end = x0;

	// Loop the scanlines from top to bottom
	for (int y = y0; y <= y2; y++)
	{
		draw_line(x_start, y, x_end, y, color);

		x_start += inverse_slope_1;
		x_end += inverse_slope_2;
	}
}

void fill_flat_top_triangle_flat(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color)
{
	// Find the two slopes from the two triangle legs
	float inverse_slope_1 = (float)(x2 - x0) / (y2 - y0);
	float inverse_slope_2 = (float)(x2 - x1) / (y2 - y1);

	// Start x_start and x_end from the top vertex (x2, y2)
	float x_start = x2;
	float x_end = x2;

	// Loop the scanlines from top to bottom
	for (int y = y2; y >= y0; y--)
	{
		draw_line(x_start, y, x_end, y, color);

		x_start -= inverse_slope_1;
		x_end -= inverse_slope_2;
	}
}

void fill_flat_bottom_triangle_gouraud(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t light_intensity_per_vertex[3])
{
	// Find the two slopes from the two triangle legs
	float inverse_slope_1 = (float)(x1 - x0) / (y1 - y0);
	float inverse_slope_2 = (float)(x2 - x0) / (y2 - y0);

	// Start x_start and x_end from the top vertex (x0, y0)
	float x_start = x0;
	float x_end = x0;

	// Loop the scanlines from top to bottom
	for (int y = y0; y <= y2; y++)
	{
		float lerpStart = inverse_lerp(y0, y1, y);
		float lerpEnd = inverse_lerp(y0, y2, y);
		uint32_t colorStart = lerp_color(light_intensity_per_vertex[0], light_intensity_per_vertex[1], lerpStart);
		uint32_t colorEnd = lerp_color(light_intensity_per_vertex[0], light_intensity_per_vertex[2], lerpEnd);
		draw_line_gouraud(x_start, y, x_end, y, colorStart, colorEnd);

		x_start += inverse_slope_1;
		x_end += inverse_slope_2;
	}
}

void fill_flat_top_triangle_gouraud(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t light_intensity_per_vertex[3])
{
	// Find the two slopes from the two triangle legs
	float inverse_slope_1 = (float)(x2 - x0) / (y2 - y0);
	float inverse_slope_2 = (float)(x2 - x1) / (y2 - y1);

	// Start x_start and x_end from the top vertex (x2, y2)
	float x_start = x2;
	float x_end = x2;

	// Loop the scanlines from top to bottom
	for (int y = y2; y >= y0; y--)
	{
		float lerpStart = inverse_lerp(y2, y0, y);
		float lerpEnd = inverse_lerp(y2, y1, y);
		uint32_t colorStart = lerp_color(light_intensity_per_vertex[2], light_intensity_per_vertex[0], lerpStart);
		uint32_t colorEnd = lerp_color(light_intensity_per_vertex[2], light_intensity_per_vertex[1], lerpEnd);

		draw_line_gouraud(x_start, y, x_end, y, colorStart, colorEnd);

		x_start -= inverse_slope_1;
		x_end -= inverse_slope_2;
	}
}

void fill_flat_bottom_triangle_phong(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t vertex_color, vec3_t normals[3])
{
	// Find the two slopes from the two triangle legs
	float inverse_slope_1 = (float)(x1 - x0) / (y1 - y0);
	float inverse_slope_2 = (float)(x2 - x0) / (y2 - y0);

	// Start x_start and x_end from the top vertex (x0, y0)
	float x_start = x0;
	float x_end = x0;

	// Loop the scanlines from top to bottom
	for (int y = y0; y <= y2; y++)
	{
		float lerpStart = inverse_lerp(y0, y1, y);
		float lerpEnd = inverse_lerp(y0, y2, y);
		vec3_t normal_start = lerp_vec3(normals[0], normals[1], lerpStart);
		vec3_t normal_end = lerp_vec3(normals[0], normals[2], lerpEnd);

		draw_line_phong(x_start, y, x_end, y, vertex_color, normal_start, normal_end);

		x_start += inverse_slope_1;
		x_end += inverse_slope_2;
	}
}

void fill_flat_top_triangle_phong(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t vertex_color, vec3_t normals[3])
{
	// Find the two slopes from the two triangle legs
	float inverse_slope_1 = (float)(x2 - x0) / (y2 - y0);
	float inverse_slope_2 = (float)(x2 - x1) / (y2 - y1);

	// Start x_start and x_end from the top vertex (x2, y2)
	float x_start = x2;
	float x_end = x2;

	// Loop the scanlines from top to bottom
	for (int y = y2; y >= y0; y--)
	{
		float lerpStart = inverse_lerp(y2, y0, y);
		float lerpEnd = inverse_lerp(y2, y1, y);
		vec3_t normal_start = lerp_vec3(normals[2], normals[0], lerpStart);
		vec3_t normal_end = lerp_vec3(normals[2], normals[1], lerpEnd);

		draw_line_phong(x_start, y, x_end, y, vertex_color, normal_start, normal_end);

		x_start -= inverse_slope_1;
		x_end -= inverse_slope_2;
	}
}

void draw_triangle(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color)
{
	draw_line(x0, y0, x1, y1, color);
	draw_line(x1, y1, x2, y2, color);
	draw_line(x2, y2, x0, y0, color);
}

void draw_filled_triangle_flat(int x0, int y0, float w0, int x1, int y1, float w1, int x2, int y2, float w2, uint32_t color)
{
	// Sort the vertices by y coordinate y0 < y1 < y2
	if (y0 > y1)
	{
		int_swap(&y0, &y1);
		int_swap(&x0, &x1);
		float_swap(&w0, &w1);
	}
	if (y1 > y2)
	{
		int_swap(&y1, &y2);
		int_swap(&x1, &x2);
		float_swap(&w1, &w2);
	}
	if (y0 > y1)
	{
		int_swap(&y0, &y1);
		int_swap(&x0, &x1);
		float_swap(&w0, &w1);
	}

	vec3_t point_a = { x0, y0, w0 };
	vec3_t point_b = { x1, y1, w1 };
	vec3_t point_c = { x2, y2, w2 };

	// Render the upper part of the triangle (flat bottom)
	float inv_slope_1 = 0;
	float inv_slope_2 = 0;

	if (y1 - y0 != 0) inv_slope_1 = (float)(x1 - x0) / abs(y1 - y0); //Left leg
	if (y2 - y0 != 0) inv_slope_2 = (float)(x2 - x0) / abs(y2 - y0); //Right leg

	if (y1 - y0 != 0)
	{
		for (int y = y0; y <= y1; y++)
		{
			int x_start = x1 + (y - y1) * inv_slope_1;
			int x_end = x0 + (y - y0) * inv_slope_2;

			// Ensure scanlines are drawn left to right
			if (x_end < x_start)
			{
				int_swap(&x_start, &x_end);
			}

			for (int x = x_start; x < x_end; x++)
			{
				// Draw our pixel with the color per triangle
				draw_pixel_zbuffered(x, y, point_a, point_b, point_c, color);
			}
		}
	}

	// Render the bottom part of the triangle (flat top)
	inv_slope_1 = 0;
	inv_slope_2 = 0;

	if (y2 - y1 != 0) inv_slope_1 = (float)(x2 - x1) / abs(y2 - y1); //Left leg
	if (y2 - y0 != 0) inv_slope_2 = (float)(x2 - x0) / abs(y2 - y0); //Right leg

	if (y2 - y1 != 0)
	{
		for (int y = y1; y <= y2; y++)
		{
			int x_start = x1 + (y - y1) * inv_slope_1;
			int x_end = x0 + (y - y0) * inv_slope_2;

			// Ensure scanlines are drawn left to right
			if (x_end < x_start)
			{
				int_swap(&x_start, &x_end);
			}

			for (int x = x_start; x < x_end; x++)
			{
				// Draw our pixel with the color per triangle
				draw_pixel_zbuffered(x, y, point_a, point_b, point_c, color);
			}
		}
	}
}

vec3_t barycentric_weights(vec2_t a, vec2_t b, vec2_t c, vec2_t p)
{
	// Find the vectors between the vertices ABC and point P
	vec2_t ac = vec2_sub(c, a);
	vec2_t ab = vec2_sub(b, a);
	vec2_t ap = vec2_sub(p, a);
	vec2_t pc = vec2_sub(c, p);
	vec2_t pb = vec2_sub(b, p);

	// Compute the area of the full parallelogram/triangle ABC using 2D cross product
	float area_parallelogram_abc = (ac.x * ab.y - ac.y * ab.x); // ||AC x AB||

	// Alpha is the area of the small parallelogram/triangle PBC divided by the area of the full parallelogram/triangle ABC
	float alpha = (pc.x * pb.y - pc.y * pb.x) / area_parallelogram_abc;

	// Beta is the area of the small parallelogram/triangle PBC divided by the area of the full parallelogram/triangle ABC
	float beta = (ac.x * ap.y - ac.y * ap.x) / area_parallelogram_abc;

	// Weight gamm is easily found since barycentric coordinates always add up to 1.0
	float gamma = 1 - alpha - beta;

	vec3_t weights = { alpha, beta, gamma };
	return weights;
}

// Draw the color at pixel x and y using depth interpolation
void draw_triangle_pixel(
	int x, int y, uint32_t color,
	vec4_t point_a, vec4_t point_b, vec4_t point_c
) {
	// Create three vec2 to find the interpolation
	vec2_t p = { x, y };
	vec2_t a = vec2_from_vec4(point_a);
	vec2_t b = vec2_from_vec4(point_b);
	vec2_t c = vec2_from_vec4(point_c);

	// Calculate the barycentric coordinates of our point 'p' inside the triangle
	vec3_t weights = barycentric_weights(a, b, c, p);

	float alpha = weights.x;
	float beta = weights.y;
	float gamma = weights.z;

	// Interpolate the value of 1/w for the current pixel
	float interpolated_reciprocal_w = (1 / point_a.w) * alpha + (1 / point_b.w) * beta + (1 / point_c.w) * gamma;

	// Adjust 1/w so the pixels that are closer to the camera have smaller values
	interpolated_reciprocal_w = 1.0 - interpolated_reciprocal_w;

	// Only draw the pixel if the depth value is less than the one previously stored in the z-buffer
	if (interpolated_reciprocal_w < get_zbuffer_at(x, y)) {
		// Draw a pixel at position (x,y) with a solid color
		draw_pixel(x, y, color);

		// Update the z-buffer value with the 1/w of this current pixel
		update_zbuffer_at(x, y, interpolated_reciprocal_w);
	}
}

// Draw the texture pixel at position x and y using depth interpolation
void draw_texel(int x, int y, uint32_t* texture, vec3_t point_a, vec3_t point_b, vec3_t point_c, tex2_t uv_a, tex2_t uv_b, tex2_t uv_c)
{
	vec2_t point_p = { x,y };
	vec2_t a = vec2_from_vec3(point_a);
	vec2_t b = vec2_from_vec3(point_b);
	vec2_t c = vec2_from_vec3(point_c);

	vec3_t weights = barycentric_weights(a, b, c, point_p);

	float alpha = weights.x;
	float beta = weights.y;
	float gamma = weights.z;

	float interpolated_reciprocal_w;

	// Interpolate the value of 1/w for the current pixel
	interpolated_reciprocal_w = (1 / point_a.z) * alpha + (1 / point_b.z) * beta + (1 / point_c.z) * gamma;

	// Only draw the pixel if the depth value is less than the one previously stored in the z-buffer
	if (interpolated_reciprocal_w > get_zbuffer_at(x, y))
	{
		float interpolated_u;
		float interpolated_v;

		// Perform the interpolation of all U/w and V/w values using barycentric weights and a factor of 1/w for perspectice correct texture coord
		interpolated_u = (uv_a.u / point_a.z) * alpha + (uv_b.u / point_b.z) * beta + (uv_c.u / point_c.z) * gamma;
		interpolated_v = (uv_a.v / point_a.z) * alpha + (uv_b.v / point_b.z) * beta + (uv_c.v / point_c.z) * gamma;

		// Now we cam divide back both interpolated values by 1/w
		interpolated_u /= interpolated_reciprocal_w;
		interpolated_v /= interpolated_reciprocal_w;

		// Scale/map the UV coordinates to the full width and height of the texture.
		int tex_x = abs((int)(interpolated_u * texture_width)) % texture_width;
		int tex_y = abs((int)(interpolated_v * texture_height)) % texture_height;

		// look up point in texture and draw the pixel
		draw_pixel(x, y, texture[((texture_width * tex_y) + tex_x)]);

		// Update the z-buffer value with the 1/w of this current pixel
		update_zbuffer_at(x, y, interpolated_reciprocal_w);
	}
}

void draw_textured_triangle(
	int x0, int y0, float w0, float u0, float v0,
	int x1, int y1, float w1, float u1, float v1,
	int x2, int y2, float w2, float u2, float v2,
	uint32_t* texture)
{
	// Sort the vertices and uvs by y-coordinate y0 < y1 < y2
	if (y0 > y1)
	{
		int_swap(&y0, &y1);
		int_swap(&x0, &x1);
		float_swap(&w0, &w1);
		float_swap(&u0, &u1);
		float_swap(&v0, &v1);
	}
	if (y1 > y2)
	{
		int_swap(&y1, &y2);
		int_swap(&x1, &x2);
		float_swap(&w1, &w2);
		float_swap(&u1, &u2);
		float_swap(&v1, &v2);
	}
	if (y0 > y1)
	{
		int_swap(&y0, &y1);
		int_swap(&x0, &x1);
		float_swap(&w0, &w1);
		float_swap(&u0, &u1);
		float_swap(&v0, &v1);
	}

	vec3_t point_a = { x0, y0, w0 };
	vec3_t point_b = { x1, y1, w1 };
	vec3_t point_c = { x2, y2, w2 };
	tex2_t uv_a = { u0, v0 };
	tex2_t uv_b = { u1, v1 };
	tex2_t uv_c = { u2, v2 };

	// Render the upper part of the triangle (flat bottom)
	float inv_slope_1 = 0;
	float inv_slope_2 = 0;

	if (y1 - y0 != 0) inv_slope_1 = (float)(x1 - x0) / abs(y1 - y0); //Left leg
	if (y2 - y0 != 0) inv_slope_2 = (float)(x2 - x0) / abs(y2 - y0); //Right leg

	if (y1 - y0 != 0)
	{
		for (int y = y0; y <= y1; y++)
		{
			int x_start = x1 + (y - y1) * inv_slope_1;
			int x_end = x0 + (y - y0) * inv_slope_2;

			// Ensure scanlines are drawn left to right
			if (x_end < x_start)
			{
				int_swap(&x_start, &x_end);
			}

			for (int x = x_start; x < x_end; x++)
			{
				// Draw our pixel with the color that comes from the texture
				draw_texel(x, y, texture, point_a, point_b, point_c, uv_a, uv_b, uv_c);
			}
		}
	}

	// Render the bottom part of the triangle (flat top)
	inv_slope_1 = 0;
	inv_slope_2 = 0;

	if (y2 - y1 != 0) inv_slope_1 = (float)(x2 - x1) / abs(y2 - y1); //Left leg
	if (y2 - y0 != 0) inv_slope_2 = (float)(x2 - x0) / abs(y2 - y0); //Right leg

	if (y2 - y1 != 0)
	{
		for (int y = y1; y <= y2; y++)
		{
			int x_start = x1 + (y - y1) * inv_slope_1;
			int x_end = x0 + (y - y0) * inv_slope_2;

			// Ensure scanlines are drawn left to right
			if (x_end < x_start)
			{
				int_swap(&x_start, &x_end);
			}

			for (int x = x_start; x < x_end; x++)
			{
				// Draw our pixel with the color that comes from the texture
				draw_texel(x, y, texture, point_a, point_b, point_c, uv_a, uv_b, uv_c);
			}
		}
	}
}

void draw_filled_triangle_gouraud(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t light_intensity_per_vertex[3])
{
	// Sort the vertices by y coordinate y0 < y1 < y2
	// Sort the light intensities aswell
	if (y0 > y1)
	{
		int_swap(&y0, &y1);
		int_swap(&x0, &x1);
		color_swap(&light_intensity_per_vertex[0], &light_intensity_per_vertex[1]);
	}
	if (y1 > y2)
	{
		int_swap(&y1, &y2);
		int_swap(&x1, &x2);
		color_swap(&light_intensity_per_vertex[1], &light_intensity_per_vertex[2]);
	}
	if (y0 > y1)
	{
		int_swap(&y0, &y1);
		int_swap(&x0, &x1);
		color_swap(&light_intensity_per_vertex[0], &light_intensity_per_vertex[1]);
	}

	if (y1 == y2)
	{
		// Only draw the flat bottom triangle
		fill_flat_bottom_triangle_gouraud(x0, y0, x1, y1, x2, y2, light_intensity_per_vertex);
	}
	else if (y0 == y1)
	{
		// Only draw the flat top triangle
		fill_flat_top_triangle_gouraud(x0, y0, x1, y1, x2, y2, light_intensity_per_vertex);
	}
	else
	{
		// Calculate the new vertex (Mx, My) using triangle similarity
		int My = y1;
		int Mx = ((float)((x2 - x0) * (y1 - y0)) / (float)(y2 - y0)) + x0;

		// Lerp mid point intensity
		float lerp = inverse_lerp(y0, y2, My);
		// Use that float to Lerp the Intensity between x0,y0 and x2,y2 to get the light intensity at mx,my
		uint32_t light_intensity_My = lerp_color(light_intensity_per_vertex[0], light_intensity_per_vertex[2], lerp);

		uint32_t light_intensity_xy2 = light_intensity_per_vertex[2];
		light_intensity_per_vertex[2] = light_intensity_My;

		fill_flat_bottom_triangle_gouraud(x0, y0, x1, y1, Mx, My, light_intensity_per_vertex);

		light_intensity_per_vertex[0] = light_intensity_per_vertex[1];
		light_intensity_per_vertex[1] = light_intensity_My;
		light_intensity_per_vertex[2] = light_intensity_xy2;

		fill_flat_top_triangle_gouraud(x1, y1, Mx, My, x2, y2, light_intensity_per_vertex);
	}
}

void draw_filled_triangle_phong(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t vertex_color, vec3_t normals[3])
{
	// Sort the vertices by y coordinate y0 < y1 < y2
	if (y0 > y1)
	{
		int_swap(&y0, &y1);
		int_swap(&x0, &x1);
		vec3_swap(&normals[0], &normals[1]);
	}
	if (y1 > y2)
	{
		int_swap(&y1, &y2);
		int_swap(&x1, &x2);
		vec3_swap(&normals[1], &normals[2]);
	}
	if (y0 > y1)
	{
		int_swap(&y0, &y1);
		int_swap(&x0, &x1);
		vec3_swap(&normals[0], &normals[1]);
	}

	if (y1 == y2)
	{
		// Only draw the flat bottom triangle
		fill_flat_bottom_triangle_phong(x0, y0, x1, y1, x2, y2, vertex_color, normals);
	}
	else if (y0 == y1)
	{
		// Only draw the flat top triangle
		fill_flat_top_triangle_phong(x0, y0, x1, y1, x2, y2, vertex_color, normals);
	}
	else
	{
		// Calculate the new vertex (Mx, My) using triangle similarity
		int My = y1;
		int Mx = ((float)((x2 - x0) * (y1 - y0)) / (float)(y2 - y0)) + x0;

		// Lerp mid point intensity
		float lerp = inverse_lerp(y0, y2, My);
		// Use that float to Lerp the Intensity between x0,y0 and x2,y2 to get the light intensity at mx,my
		vec3_t normal_My = lerp_vec3(normals[0], normals[2], lerp);

		vec3_t normal_xy2 = normals[2];
		normals[2] = normal_My;

		fill_flat_bottom_triangle_phong(x0, y0, x1, y1, Mx, My, vertex_color, normals);

		normals[0] = normals[1];
		normals[1] = normal_My;
		normals[2] = normal_xy2;

		fill_flat_top_triangle_phong(x1, y1, Mx, My, x2, y2, vertex_color, normals);
	}
}

float lerp_float(float a, float b, float f)
{
	return (a * (1.0 - f)) + (b * f);
}

vec3_t lerp_vec3(vec3_t a, vec3_t b, float f)
{
	vec3_t result;
	result.x = lerp_float(a.x, b.x, f);
	result.y = lerp_float(a.y, b.y, f);
	result.z = lerp_float(a.z, b.z, f);
	return result;
}

float inverse_lerp(float a, float b, float value)
{
	return (value - a) / (b - a);
}

uint32_t lerp_color(uint32_t a, uint32_t b, float f)
{
	// Ensure f is in the range [0, 1]
	f = f < 0.0f ? 0.0f : (f > 1.0f ? 1.0f : f);

	// Extract the individual color channels
	uint8_t a_a = (a >> 24) & 0xFF;
	uint8_t a_r = (a >> 16) & 0xFF;
	uint8_t a_g = (a >> 8) & 0xFF;
	uint8_t a_b = a & 0xFF;

	uint8_t b_a = (b >> 24) & 0xFF;
	uint8_t b_r = (b >> 16) & 0xFF;
	uint8_t b_g = (b >> 8) & 0xFF;
	uint8_t b_b = b & 0xFF;

	// Perform linear interpolation for each channel
	uint8_t result_a = (uint8_t)(a_a * (1.0f - f) + b_a * f);
	uint8_t result_r = (uint8_t)(a_r * (1.0f - f) + b_r * f);
	uint8_t result_g = (uint8_t)(a_g * (1.0f - f) + b_g * f);
	uint8_t result_b = (uint8_t)(a_b * (1.0f - f) + b_b * f);

	// Combine the interpolated channels into the final color
	return ((uint32_t)result_a << 24) | ((uint32_t)result_r << 16) | ((uint32_t)result_g << 8) | result_b;
}

uint32_t lerp_color_sRGB(uint32_t color1, uint32_t color2, float fraction)
{
	uint8_t r1 = (color1 >> 16) & 0xFF;
	uint8_t g1 = (color1 >> 8) & 0xFF;
	uint8_t b1 = color1 & 0xFF;

	uint8_t r2 = (color2 >> 16) & 0xFF;
	uint8_t g2 = (color2 >> 8) & 0xFF;
	uint8_t b2 = color2 & 0xFF;

	// Convert to linear sRGB space
	float R1 = powf(r1 / 255.0f, 2.2f);
	float G1 = powf(g1 / 255.0f, 2.2f);
	float B1 = powf(b1 / 255.0f, 2.2f);

	float R2 = powf(r2 / 255.0f, 2.2f);
	float G2 = powf(g2 / 255.0f, 2.2f);
	float B2 = powf(b2 / 255.0f, 2.2f);

	// Interpolate in linear sRGB space
	float R = (1 - fraction) * R1 + fraction * R2;
	float G = (1 - fraction) * G1 + fraction * G2;
	float B = (1 - fraction) * B1 + fraction * B2;

	// Convert back to sRGB space
	uint8_t r = (uint8_t)(powf(R, 1 / 2.2f) * 255.0f + 0.5f);
	uint8_t g = (uint8_t)(powf(G, 1 / 2.2f) * 255.0f + 0.5f);
	uint8_t b = (uint8_t)(powf(B, 1 / 2.2f) * 255.0f + 0.5f);

	return (r << 16) | (g << 8) | b;
}

uint32_t multiply_colors(uint32_t color1, uint32_t color2) {
	uint32_t a1 = (color1 >> 24) & 0xff;
	uint32_t r1 = (color1 >> 16) & 0xff;
	uint32_t g1 = (color1 >> 8) & 0xff;
	uint32_t b1 = color1 & 0xff;

	uint32_t a2 = (color2 >> 24) & 0xff;
	uint32_t r2 = (color2 >> 16) & 0xff;
	uint32_t g2 = (color2 >> 8) & 0xff;
	uint32_t b2 = color2 & 0xff;

	uint32_t a = (a1 * a2) / 255;
	uint32_t r = (r1 * r2) / 255;
	uint32_t g = (g1 * g2) / 255;
	uint32_t b = (b1 * b2) / 255;

	return (a << 24) | (r << 16) | (g << 8) | b;
}

uint32_t add_colors(uint32_t color1, uint32_t color2) {
	uint32_t a1 = (color1 >> 24) & 0xff;
	uint32_t r1 = (color1 >> 16) & 0xff;
	uint32_t g1 = (color1 >> 8) & 0xff;
	uint32_t b1 = color1 & 0xff;

	uint32_t a2 = (color2 >> 24) & 0xff;
	uint32_t r2 = (color2 >> 16) & 0xff;
	uint32_t g2 = (color2 >> 8) & 0xff;
	uint32_t b2 = color2 & 0xff;

	uint32_t a = a1 + a2;
	uint32_t r = r1 + r2;
	uint32_t g = g1 + g2;
	uint32_t b = b1 + b2;

	if (a > 255) a = 255;
	if (r > 255) r = 255;
	if (g > 255) g = 255;
	if (b > 255) b = 255;

	return (a << 24) | (r << 16) | (g << 8) | b;
}

uint32_t subtract_colors(uint32_t color1, uint32_t color2) {
	uint32_t a1 = (color1 >> 24) & 0xff;
	uint32_t r1 = (color1 >> 16) & 0xff;
	uint32_t g1 = (color1 >> 8) & 0xff;
	uint32_t b1 = color1 & 0xff;

	uint32_t a2 = (color2 >> 24) & 0xff;
	uint32_t r2 = (color2 >> 16) & 0xff;
	uint32_t g2 = (color2 >> 8) & 0xff;
	uint32_t b2 = color2 & 0xff;

	uint32_t a = (a1 > a2) ? a1 - a2 : 0;
	uint32_t r = (r1 > r2) ? r1 - r2 : 0;
	uint32_t g = (g1 > g2) ? g1 - g2 : 0;
	uint32_t b = (b1 > b2) ? b1 - b2 : 0;

	return (a << 24) | (r << 16) | (g << 8) | b;
}
