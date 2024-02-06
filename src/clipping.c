#include <math.h>
#include <stdio.h>
#include "clipping.h"

#define NUM_PLANES 6
plane_t frustum_planes[NUM_PLANES];

void init_frustum_planes(float fov_x, float fov_y, float z_near, float z_far)
{
	float cos_half_fov_x = cos(fov_x * 0.5);
	float sin_half_fov_x = sin(fov_x * 0.5);
	float cos_half_fov_y = cos(fov_y * 0.5);
	float sin_half_fov_y = sin(fov_y * 0.5);

	vec3_t origin = vec3_new(0, 0, 0);

	frustum_planes[LEFT_FRUSTUM_PLANE].point = origin;
	frustum_planes[LEFT_FRUSTUM_PLANE].normal.x = cos_half_fov_x;
	frustum_planes[LEFT_FRUSTUM_PLANE].normal.y = 0;
	frustum_planes[LEFT_FRUSTUM_PLANE].normal.z = sin_half_fov_x;

	frustum_planes[RIGHT_FRUSTUM_PLANE].point = origin;
	frustum_planes[RIGHT_FRUSTUM_PLANE].normal.x = -cos_half_fov_x;
	frustum_planes[RIGHT_FRUSTUM_PLANE].normal.y = 0;
	frustum_planes[RIGHT_FRUSTUM_PLANE].normal.z = sin_half_fov_x;

	frustum_planes[TOP_FRUSTUM_PLANE].point = origin;
	frustum_planes[TOP_FRUSTUM_PLANE].normal.x = 0;
	frustum_planes[TOP_FRUSTUM_PLANE].normal.y = -cos_half_fov_y;
	frustum_planes[TOP_FRUSTUM_PLANE].normal.z = sin_half_fov_y;

	frustum_planes[BOTTOM_FRUSTUM_PLANE].point = origin;
	frustum_planes[BOTTOM_FRUSTUM_PLANE].normal.x = 0;
	frustum_planes[BOTTOM_FRUSTUM_PLANE].normal.y = cos_half_fov_y;
	frustum_planes[BOTTOM_FRUSTUM_PLANE].normal.z = sin_half_fov_y;

	frustum_planes[NEAR_FRUSTUM_PLANE].point = vec3_new(0, 0, z_near);
	frustum_planes[NEAR_FRUSTUM_PLANE].normal.x = 0;
	frustum_planes[NEAR_FRUSTUM_PLANE].normal.y = 0;
	frustum_planes[NEAR_FRUSTUM_PLANE].normal.z = 1;

	frustum_planes[FAR_FRUSTUM_PLANE].point = vec3_new(0, 0, z_far);
	frustum_planes[FAR_FRUSTUM_PLANE].normal.x = 0;
	frustum_planes[FAR_FRUSTUM_PLANE].normal.y = 0;
	frustum_planes[FAR_FRUSTUM_PLANE].normal.z = -1;
}

polygon_t create_polygon_from_triangle(vec3_t v0, vec3_t v1, vec3_t v2, tex2_t t0, tex2_t t1, tex2_t t2)
{
	polygon_t polygon = {
		.vertices = { v0, v1, v2},
		.texcoords = {t0, t1, t2},
		.num_vertices = 3
	};
	return polygon;
}

//In a fan shape
void triangles_from_polygon(polygon_t* polygon, triangle_t triangles[], int* num_triangles_after_clipping)
{
	for (int i = 0; i < polygon->num_vertices - 2; i++)
	{
		int index0 = 0;
		int index1 = i + 1;
		int index2 = i + 2;

		triangles[i].projected_points[0] = vec4_from_vec3(polygon->vertices[index0]);
		triangles[i].projected_points[1] = vec4_from_vec3(polygon->vertices[index1]);
		triangles[i].projected_points[2] = vec4_from_vec3(polygon->vertices[index2]);

		triangles[i].texcoords[0] = polygon->texcoords[index0];
		triangles[i].texcoords[1] = polygon->texcoords[index1];
		triangles[i].texcoords[2] = polygon->texcoords[index2];
	}
	*num_triangles_after_clipping = polygon->num_vertices - 2;
}

float float_lerp01(float a, float b, float t)
{
	return a + t * (b - a);
}

void clip_polygon_against_plane(polygon_t* polygon, int plane)
{
	vec3_t plane_point = frustum_planes[plane].point;
	vec3_t plane_normal = frustum_planes[plane].normal;

	// The array of inside vertices that will be part of the final polygon returned via parameter
	vec3_t inside_vertices[MAX_NUM_POLY_VERTICES];
	tex2_t inside_texcoords[MAX_NUM_POLY_VERTICES];
	int num_inside_vertices = 0;

	// Start current vertex with the first polygon vertex and texture coordinate
	vec3_t* current_vertex = &polygon->vertices[0];
	tex2_t* current_texcoord = &polygon->texcoords[0];

	// Start previous vertex with the last polygon vertex and texture coordinate
	vec3_t* previous_vertex = &polygon->vertices[polygon->num_vertices - 1];
	tex2_t* previous_texcoord = &polygon->texcoords[polygon->num_vertices - 1];

	// Start the current and previous dot product
	float current_dot = 0;
	float previous_dot = vec3_dot(vec3_sub(*previous_vertex, plane_point), plane_normal);

	// Loop while the current vertex is different than the last vertex
	while (current_vertex != &polygon->vertices[polygon->num_vertices])
	{
		// Is current vertex inside or outside of the plane
		current_dot = vec3_dot(vec3_sub(*current_vertex, plane_point), plane_normal);

		if (current_dot * previous_dot < 0)
		{
			// Calculate the interpolation factor, t = dotQ1 / (dotQ1 - dotQ2)
			float interpolation_factor = previous_dot / (previous_dot - current_dot);
			// Calculate the intersection point by lerping the vertex position, I = Q1 + t(Q2 - Q1)
			vec3_t intersection_point = {
				.x = float_lerp01(previous_vertex->x, current_vertex->x, interpolation_factor),
				.y = float_lerp01(previous_vertex->y, current_vertex->y, interpolation_factor),
				.z = float_lerp01(previous_vertex->z, current_vertex->z, interpolation_factor)
			};

			// Use the lerp formula to get the interpolated U and V texture coordinates
			tex2_t interpolated_texcoord = {
				.u = float_lerp01(previous_texcoord->u, current_texcoord->u, interpolation_factor),
				.v = float_lerp01(previous_texcoord->v, current_texcoord->v, interpolation_factor)
			};

			// Insert the new intersection point in the list of "inside vertices"
			inside_vertices[num_inside_vertices] = vec3_clone(&intersection_point);
			inside_texcoords[num_inside_vertices] = tex2_clone(&interpolated_texcoord);
			num_inside_vertices++;
		}

		// If current point is inside the plane
		if (current_dot > 0)
		{
			// Insert current vertex in the list of "inside_vertices"
			inside_vertices[num_inside_vertices] = vec3_clone(current_vertex);
			inside_texcoords[num_inside_vertices] = tex2_clone(current_texcoord);
			num_inside_vertices++;
		}

		// Move to the next vertex
		previous_dot = current_dot;
		previous_vertex = current_vertex;
		previous_texcoord = current_texcoord;
		current_vertex++; //Increment the pointer with point arithmetic (same like an array)
		current_texcoord++;
	}
	// Copy all the vertices from the inside_vertices into the destination polygon parameter
	for (int i = 0; i < num_inside_vertices; i++)
	{
		polygon->vertices[i] = vec3_clone(&inside_vertices[i]);
		polygon->texcoords[i] = tex2_clone(&inside_texcoords[i]);
	}
	polygon->num_vertices = num_inside_vertices;
}

void clip_polygon(polygon_t* polygon)
{
	clip_polygon_against_plane(polygon, LEFT_FRUSTUM_PLANE);
	clip_polygon_against_plane(polygon, RIGHT_FRUSTUM_PLANE);
	clip_polygon_against_plane(polygon, TOP_FRUSTUM_PLANE);
	clip_polygon_against_plane(polygon, BOTTOM_FRUSTUM_PLANE);
	clip_polygon_against_plane(polygon, NEAR_FRUSTUM_PLANE);
	clip_polygon_against_plane(polygon, FAR_FRUSTUM_PLANE);
}
