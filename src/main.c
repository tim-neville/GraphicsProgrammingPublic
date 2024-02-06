#include "globals.h"
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <SDL.h>
#include "upng.h"
#include "array.h"
#include "display.h"
#include "vector.h"
#include "mesh.h"
#include "camera.h"
#include "triangle.h"
#include "matrix.h"
#include "light.h"
#include "texture.h"
#include "clipping.h"
#include "tims.h"

#define USE_VERTEX_COLORS = true
const uint32_t BACKGROUND_COLOR = 0xFF79D9FC; //0xFF202020;

triangle_t* triangles_to_render = NULL;

bool is_running = false;
int previous_frame_time = 0;
float delta_time = 0;

bool show_grid = false;
bool enable_culling = true;
bool show_vertices = false;
bool show_wireframe = false;
bool show_solid = false;
bool show_normals = false;
bool show_textured = true;
bool auto_rotate = false;
bool show_depth = false;

bool is_fullscreen = false;
bool show_window_border = false;
bool is_orthographic = false;
int grid_size = 40;
vec3_t position = { 0, 0, 0 };
vec2_t rotation = { 0, 0 };
float z_near = 0.1;
float z_far = 2000;
float fov_degrees = 60;
float fov_x = 0;
float fov_y = 0;
float orthographic_size = 5;

int num_threads = 0;
SDL_Thread** threadIDs;
bool measure_time = false;

void draw_a_thing()
{
	if (show_grid)
	{
		vec3_t points[9];
		for (int y = 0; y < 3; y++)
		{
			for (int x = 0; x < 3; x++)
			{
				vec3_t p = { x * grid_size, y * grid_size, 0 };
				p = mat3x4_mul_vec3(world_matrix, p);
				vec4_t proj_p = mat4x4_mul_vec4_project(projection_matrix, p);

				// Scale into the view
				proj_p.x *= (float)(window_width >> 1);
				proj_p.y *= (float)(window_height >> 1);

				// Translate the projected points to the middle of the screen
				proj_p.x += (float)(window_width >> 1);
				proj_p.y += (float)(window_height >> 1);

				points[3 * y + x] = vec3_from_vec4(proj_p);
			}
		}

		draw_grid_world(grid_size, points, 9);
	}
	else
	{
		// Start in Model Space
		vec3_t point_a = { 0,0,0 };
		vec3_t point_b = { 1,1,0 };


		point_a = mat3x4_mul_vec3(world_matrix, point_a);
		point_b = mat3x4_mul_vec3(world_matrix, point_b);

		// Now we're in World Space

		point_a = mat3x4_mul_vec3(view_matrix, point_a);
		point_b = mat3x4_mul_vec3(view_matrix, point_b);

		// Now we're in View Space

		vec4_t point_a_proj = mat4x4_mul_vec4_project(projection_matrix, point_a);
		vec4_t point_b_proj = mat4x4_mul_vec4_project(projection_matrix, point_b);

		// Now we're in Clip Space

		// Perform perspective divide with original z-value that is now stored in w
		if (point_a_proj.w != 0.0)
		{
			point_a_proj.x /= point_a_proj.w;
			point_a_proj.y /= point_a_proj.w;
			point_a_proj.z /= point_a_proj.w;
		}

		//Flip Y because our world is Y down for some reason
		point_a_proj.y *= -1;
		point_b_proj.y *= -1;
		// Now in Image Space or Normalized Device Space (NDC)

		// After projection, scale and translate the points based on the screen width and height
		point_a_proj.x *= (float)(window_width >> 1);
		point_a_proj.y *= (float)(window_height >> 1);

		point_a_proj.x += (float)(window_width >> 1);
		point_a_proj.y += (float)(window_height >> 1);

		point_b_proj.x *= (float)(window_width >> 1);
		point_b_proj.y *= (float)(window_height >> 1);

		point_b_proj.x += (float)(window_width >> 1);
		point_b_proj.y += (float)(window_height >> 1);

		// Now in screen space

		draw_line(point_a_proj.x, point_a_proj.y, point_b_proj.x, point_b_proj.y, 0xFFFF00FF);
	}
}

void init_projection_matrix()
{
	float aspect_x = (float)window_width / (float)window_height;
	float aspect_y = (float)window_height / (float)window_width;
	fov_y = fov_degrees * deg_to_rad; //60 degrees in radians
	fov_x = atan(tan(fov_y / 2) * aspect_x) * 2.0;
	projection_matrix = mat4x4_make_perspective(fov_y, aspect_y, z_near, z_far);
}

void init_orthographic_matrix()
{
	projection_matrix = mat4x4_make_orthographic(window_width, window_height, orthographic_size, z_near, z_far);
}

void setup(void)
{
	// Get Native Screen Info
	num_threads = SDL_GetCPUCount() - 1;
	threadIDs = (SDL_Thread**)malloc(num_threads * sizeof(SDL_Thread*));

	if (threadIDs == NULL)
	{
		DebugError("Failed to allocate memory for threads");
	}

	position.z = mesh.translation.z;

	//Position shark
	//fov_degrees = 20;
	//position.z = 7;
	//position.y = 0.225;

	// Initialize the perspective projection matrix
	init_projection_matrix();

	// Initialize frustum planes with a point and a normal
	init_frustum_planes(fov_x, fov_y, z_near, z_far);

	// Manually load hardcoded texture data from the static array
	//mesh_texture = (uint32_t*)REDBRICK_TEXTURE;
	//texture_width = 64;
	//texture_height = 64;

	//load_cube_mesh_data();
	load_obj_file_data("../../assets/cube.obj");
	load_png_texture_data("../../assets/cube.png");
	// 
	//load_obj_file_data("../../assets/cube_gustavo.obj");
	//load_png_texture_data("../../assets/cube_textured.png");

	//load_obj_file_data("../../assets/f22.obj");
	//load_png_texture_data("../../assets/f22.png");
	//load_obj_file_data("../../assets/depth.obj");
	//load_obj_file_data("../../assets/sponza.obj");
	//load_png_texture_data("../../assets/pikuma.png");
	//load_obj_file_data("../../assets/f117.obj");
	//load_png_texture_data("../../assets/f117.png");
	//load_obj_file_data("../../assets/sphere.obj");
	//load_png_texture_data("../../assets/pikuma.png");
	//load_obj_file_data("../../assets/crab.obj");
	//load_png_texture_data("../../assets/crab.png");
	//load_obj_file_data("../../assets/drone.obj");
	//load_png_texture_data("../../assets/drone.png");
	//load_obj_file_data("../../assets/efa.obj");
	//load_png_texture_data("../../assets/efa.png");

	//load_obj_file_data("../../assets/test.obj");
	//load_obj_file_data("../../assets/shark.obj");
	//load_obj_file_data("../../assets/shark_smooth.obj");
	//load_obj_file_data("../../assets/shark_smooth_hp.obj");
	//load_obj_file_data("../../assets/icosphere.obj");
	//load_obj_file_data("../../assets/icosphere2.obj");
	//load_obj_file_data("../../assets/triangle.obj");
	//load_obj_file_data("../../assets/monkey.obj");
	//load_obj_file_data("../../assets/cube_textured.obj");
}

void process_input(void)
{
	SDL_Event event;

	while (SDL_PollEvent(&event))
	{
		switch (event.type)
		{
		case SDL_QUIT:
			is_running = false;
			break;
		case SDL_KEYDOWN:

			if (event.key.keysym.sym == SDLK_ESCAPE)
			{
				is_running = false;
				break;
			}
			if (event.key.keysym.sym == SDLK_f)
			{
				is_fullscreen = !is_fullscreen;

				if (is_fullscreen)
				{
					set_window_fullscreen();
					init_projection_matrix();
					init_frustum_planes(fov_x, fov_y, z_near, z_far);
					break;
				}
				else
				{
					set_window_borderless();
					init_projection_matrix();
					init_frustum_planes(fov_x, fov_y, z_near, z_far);
					break;
				}
			}

			if (event.key.keysym.sym == SDLK_RIGHT)// && grid_size * 2 <= 40)
			{
				if (is_orthographic)
				{
					orthographic_size += 0.1;
					init_orthographic_matrix();
					init_frustum_planes(fov_x, fov_y, z_near, z_far);
					break;
				}
				else
				{
					fov_degrees += 1;
					init_projection_matrix();
					init_frustum_planes(fov_x, fov_y, z_near, z_far);
					break;
				}

				//grid_size *= 2;
				//position.x = ((int)position.x * 2) % grid_size;
				//position.y = ((int)position.y * 2) % grid_size;
			}
			if (event.key.keysym.sym == SDLK_LEFT)// && grid_size / 2 >= 5)
			{
				if (is_orthographic)
				{
					orthographic_size -= 0.1;
					init_orthographic_matrix();
					init_frustum_planes(fov_x, fov_y, z_near, z_far);
					break;
				}
				else
				{
					fov_degrees -= 1;
					init_projection_matrix();
					init_frustum_planes(fov_x, fov_y, z_near, z_far);
					break;
				}

				//grid_size /= 2;
				//position.x = ((int)position.x / 2) % grid_size;
				//position.y = ((int)position.y / 2) % grid_size;
			}

			if (event.key.keysym.sym == SDLK_UP)// && grid_size / 2 >= 5)
			{
				camera.position.y += 3.0 * delta_time;
				break;
			}
			if (event.key.keysym.sym == SDLK_DOWN)// && grid_size / 2 >= 5)
			{
				camera.position.y -= 3.0 * delta_time;
				break;
			}

			if (event.key.keysym.sym == SDLK_w)
			{
				camera.forward_velocity = vec3_mul(camera.direction, 5.0 * delta_time);
				camera.position = vec3_add(camera.position, camera.forward_velocity);
				break;
			}
			if (event.key.keysym.sym == SDLK_s)
			{
				camera.forward_velocity = vec3_mul(camera.direction, 5.0 * delta_time);
				camera.position = vec3_sub(camera.position, camera.forward_velocity);
				break;
			}
			if (event.key.keysym.sym == SDLK_a)
			{
				camera.forward_velocity = vec3_mul(vec3_cross(camera.direction, VEC3_UP), 5.0 * delta_time);
				camera.position = vec3_sub(camera.position, camera.forward_velocity);
				break;
			}
			if (event.key.keysym.sym == SDLK_d)
			{
				camera.forward_velocity = vec3_mul(vec3_cross(camera.direction, VEC3_UP), 5.0 * delta_time);
				camera.position = vec3_add(camera.position, camera.forward_velocity);
				break;
			}

			if (event.key.keysym.sym == SDLK_g)
			{
				show_grid = !show_grid;
				break;
			}
			if (event.key.keysym.sym == SDLK_c)
			{
				enable_culling = !enable_culling;
				break;
			}
			if (event.key.keysym.sym == SDLK_v)
			{
				is_orthographic = !is_orthographic;

				if (is_orthographic)
				{
					init_orthographic_matrix();
					init_frustum_planes(fov_x, fov_y, z_near, z_far);
					break;
				}
				else
				{
					init_projection_matrix();
					init_frustum_planes(fov_x, fov_y, z_near, z_far);
					break;
				}
			}
			if (event.key.keysym.sym == SDLK_1)
			{
				show_solid = !show_solid;
				break;
			}
			if (event.key.keysym.sym == SDLK_2)
			{
				show_wireframe = !show_wireframe;
				break;
			}
			if (event.key.keysym.sym == SDLK_3)
			{
				show_vertices = !show_vertices;
				break;
			}
			if (event.key.keysym.sym == SDLK_4)
			{
				show_textured = !show_textured;
				break;
			}
			if (event.key.keysym.sym == SDLK_q)
			{
				camera.position.x += 3.0 * delta_time; // Is moving in world space not local camera rotation space
				break;
			}
			if (event.key.keysym.sym == SDLK_e)
			{
				camera.position.x -= 3.0 * delta_time; // Is moving in world space not local camera rotation space
				break;
			}
			if (event.key.keysym.sym == SDLK_b)
			{
				show_window_border = !show_window_border;

				if (show_window_border)
				{
					set_window_bordered();
					break;
				}
				else
				{
					set_window_borderless();
					break;
				}
			}
			if (event.key.keysym.sym == SDLK_l)
			{
				current_lighting_method = (enum lighting_method)((current_lighting_method + 1) % 4);
				break;
			}
			if (event.key.keysym.sym == SDLK_t)
			{
				measure_time = true;
				break;
			}
			if (event.key.keysym.sym == SDLK_n)
			{
				show_normals = !show_normals;
				break;
			}
			if (event.key.keysym.sym == SDLK_r)
			{
				auto_rotate = !auto_rotate;
				break;
			}
			if (event.key.keysym.sym == SDLK_0)
			{
				show_depth = !show_depth;
				break;
			}

			break;


			//case SDL_MOUSEWHEEL:
			//{
			//	if (event.wheel.y > 0)
			//	{
			//		position.z -= increment * 2;
			// break;
			//	}
			//	else if (event.wheel.y < 0)
			//	{
			//		position.z += increment * 2;
			// break;
			//	}
			//	break;
			//}
		case SDL_MOUSEMOTION:
		{
			if (SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(SDL_BUTTON_LEFT))
			{
				//vec3_t screenspace_input = { event.motion.yrel * delta_time * 0.1, -event.motion.xrel * delta_time * 0.1 , 0 };

				vec3_t cam_dir = mat3x4_mul_vec3(world_view_matrix, camera.direction);

				vec3_normalize(&cam_dir);

				vec3_t camera_right = vec3_cross(cam_dir, VEC3_UP);
				vec3_t camera_up = vec3_cross(cam_dir, camera_right);
				camera_right = vec3_mul(camera_right, -event.motion.yrel * 0.1 * delta_time);
				camera_up = vec3_mul(camera_up, event.motion.xrel * 0.1 * delta_time);
				vec3_t cam_mov_pos = vec3_add(camera_right, camera_up);
				camera.rotation = vec3_add(camera.rotation, cam_mov_pos);
				break;
			}
			if (SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(SDL_BUTTON_RIGHT))
			{
				//vec3_t cam_dir = mat3x4_mul_vec3(world_matrix, camera.direction);

				vec3_t mouse_motion = { event.motion.xrel, event.motion.yrel, 0 };

				vec3_t camera_right = vec3_cross(camera.direction, VEC3_UP);
				vec3_t camera_up = vec3_cross(camera.direction, camera_right);
				camera_right = vec3_mul(camera_right, event.motion.xrel);
				camera_up = vec3_mul(camera_up, event.motion.yrel);
				vec3_t cam_mov_pos = vec3_add(camera_right, camera_up);
				//vec3_normalize(&cam_mov_pos);

				camera.forward_velocity = vec3_mul(cam_mov_pos, 0.56 * delta_time);
				camera.position = vec3_add(camera.position, camera.forward_velocity);
				//camera.position.y += -event.motion.yrel * delta_time * 0.56;
				//camera.position.x += -event.motion.xrel * delta_time * 0.56;
				//position.x -= event.motion.xrel * increment * 0.15;
				//position.y -= event.motion.yrel * increment * 0.15;
				break;
			}
			break;
		}
		}
	}
}

void update(void)
{
	// Wait some time until we reach the target frame time in milliseconds
	int time_to_wait = FRAME_TARGET_TIME - (SDL_GetTicks() - previous_frame_time);
	// Only delay execution if we are running too fast
	if (time_to_wait > 0 && time_to_wait <= FRAME_TARGET_TIME)
	{
		SDL_Delay(time_to_wait);
	}

	// Get a delta time factor converted to seconds to be used to update out game objects
	delta_time = (SDL_GetTicks() - previous_frame_time) / 1000.0;

	previous_frame_time = SDL_GetTicks(); //how many milliseconds since SDL init at start of game

	if (measure_time)
	{
		StartTimer();
	}

	// Initialize the array of triangles to render
	/*
	// INSTEAD OF INITALIZE THE ARRAY PER FRAME WE JUST KEEP IT LIVE AND RESET THE COUNTER (ARRAY.H / ARRAY.C).
	// On Model import we should set the size of the array instead of dynamically resizing it several times on first load
	// Possibly requires two loops thru; One to get the number of faces, verts etc, and another to actually get the data
	*/

	// Apply rotation/scale/translation
	if (auto_rotate)
	{
		mesh.rotation.x += 0.1 * delta_time;
		mesh.rotation.z += 0.1 * delta_time;
	}
	else
	{
		mesh.rotation.x = rotation.y;
		mesh.rotation.y = rotation.x;
	}

	//camera.position = position;
	//camera.direction = rotation;
	mesh.translation.x = position.x;
	mesh.translation.y = position.y;
	mesh.translation.z = position.z;


	// Create scale, rotation and translation matrices
	// Matrices must be applied in order SCALE -> ROTATION -> TRANSLATION
	mat3x4_t scale_matrix = mat3x4_make_scale(mesh.scale.x, mesh.scale.y, mesh.scale.z);
	mat3x4_t rotation_matrix_x = mat3x4_make_rotation_x(mesh.rotation.x, true);
	mat3x4_t rotation_matrix_y = mat3x4_make_rotation_y(mesh.rotation.y, true);
	mat3x4_t rotation_matrix_z = mat3x4_make_rotation_z(mesh.rotation.z, true);
	mat3x4_t translation_matrix = mat3x4_make_translation(mesh.translation.x, mesh.translation.y, mesh.translation.z);

	// Combine Scale Rotation and Translation matrices into a World Matrix. Must be in that order Scale,Rotation,Translation.
	// The order looks like it's in reverse but TRS means [T]*[R]*[S]*v which is right to left.
	world_matrix = mat3x4_identity();

	world_matrix = mat3x4_mul_mat3x4(scale_matrix, world_matrix);
	world_matrix = mat3x4_mul_mat3x4(rotation_matrix_z, world_matrix);
	world_matrix = mat3x4_mul_mat3x4(rotation_matrix_y, world_matrix);
	world_matrix = mat3x4_mul_mat3x4(rotation_matrix_x, world_matrix);
	world_matrix = mat3x4_mul_mat3x4(translation_matrix, world_matrix);

	// Calculate target look at point
	// Initialize the target looking at the positive z-axis
	vec3_t target = { 0,0,1 };
	mat3x4_t camera_rotation_matrix = mat3x4_mul_mat3x4(mat3x4_make_rotation_x(camera.rotation.x, true), mat3x4_make_rotation_y(camera.rotation.y, true));

	// Rotate camera to be pointing at target
	camera.direction = mat3x4_mul_vec3(camera_rotation_matrix, target);
	target = vec3_add(camera.position, camera.direction);

	// Create the view/camera matrix
	view_matrix = mat3x4_make_look_at(camera.position, target, VEC3_UP);

	// Multiply the vector by the view matrix to transform the scene to camera space
	world_view_matrix = mat3x4_mul_mat3x4(view_matrix, world_matrix);

	int num_faces = array_length(mesh.faces);
	for (int i = 0; i < num_faces; i++)
	{
		face_t mesh_face = mesh.faces[i];

		vec3_t face_vertices[3];
		face_vertices[0] = mesh.vertices[mesh_face.a];
		face_vertices[1] = mesh.vertices[mesh_face.b];
		face_vertices[2] = mesh.vertices[mesh_face.c];

		vec3_t vertex_normals[3];
		vertex_normals[0] = mesh.normals[mesh_face.normal_a];
		vertex_normals[1] = mesh.normals[mesh_face.normal_b];
		vertex_normals[2] = mesh.normals[mesh_face.normal_c];

		uint32_t vertex_colors[3];
		vertex_colors[0] = mesh.vertex_colors[mesh_face.color_a];
		vertex_colors[1] = mesh.vertex_colors[mesh_face.color_b];
		vertex_colors[2] = mesh.vertex_colors[mesh_face.color_c];

		vec3_t transformed_vertices[3];

		// APPLY TRANSFORMATION TO VERTICES
		for (int j = 0; j < 3; j++)
		{
			vec3_t transformed_vertex = face_vertices[j];

			// Multiply the world matrix by the original vector
			transformed_vertex = mat3x4_mul_vec3(world_view_matrix, transformed_vertex);

			// Save the transformed vertex
			transformed_vertices[j] = transformed_vertex;
		}

		// Manually Calculate Normals from Triangles For Culling
		vec3_t vector_a = transformed_vertices[0];
		vec3_t vector_b = transformed_vertices[1];
		vec3_t vector_c = transformed_vertices[2];

		vec3_t vector_ab = vec3_sub(vector_b, vector_a);
		vec3_t vector_ac = vec3_sub(vector_c, vector_a);

		// CALCULATE NORMAL PER FACE
		vec3_t face_normal = vec3_cross(vector_ab, vector_ac);
		vec3_normalize(&face_normal);

		// BACKFACE CULLING
		if (enable_culling)
		{
			vec3_t camera_ray;
			vec3_t origin = { 0,0,0 };

			if (!is_orthographic)
			{
				camera_ray = vec3_sub(origin, vector_a);
			}
			else
			{
				vec3_t forward = { 0, 0, 1 };
				camera_ray = vec3_sub(origin, vec3_add(camera.position, forward));
			}

			float dot = vec3_dot(face_normal, camera_ray);

			if (dot < 0)
			{
				continue;
			}
		}

		// Calculate Inverse Transpose Matrix for Normals Transformation
		mat3x3_t inverse_transpose_matrix = mat3x3_make_transpose_inverse(world_view_matrix);

		vec3_t transformed_normals[3];

		// APPLY TRANSFORMATION TO NORMALS
		for (int j = 0; j < 3; j++)
		{
			vec3_t transformed_normal = vertex_normals[j];

			// Tranform normal
			transformed_normal = mat3x3_mul_vec3(inverse_transpose_matrix, transformed_normal);

			// Save the transformed normal
			transformed_normals[j] = transformed_normal;
		}

		// Clipping - Create a polygon from the original transformed vertices
		polygon_t polygon = create_polygon_from_triangle(
			transformed_vertices[0],
			transformed_vertices[1],
			transformed_vertices[2],
			mesh.uvs[mesh_face.uv_a],
			mesh.uvs[mesh_face.uv_b],
			mesh.uvs[mesh_face.uv_c]
		);

		// Clip the polygon and returns a new polygon with potentially new vertices
		clip_polygon(&polygon);

		// After clipping we need to break the polygon into individual triangles
		triangle_t triangles_after_clipping[MAX_NUM_POLYGON_TRIANGLES];
		int num_triangles_after_clipping = 0;

		triangles_from_polygon(&polygon, triangles_after_clipping, &num_triangles_after_clipping);

		// Loop all the assembled triangles after clipping
		for (int t = 0; t < num_triangles_after_clipping; t++)
		{
			triangle_t triangle_after_clipping = triangles_after_clipping[t];

			vec4_t projected_points[3];
			// PERFORM PROJECTION
			// Projection Matrix is responsible for:
			// - Aspect radio - adjust x and y values based on screen width and height. (a = h/w)
			// - Field of view - adjust x and y values based on the FOV angle. (fov_factor = 1 / (tan(fov_angle / 2))
			// - Normalization - adjust x and y values to sit between -1 and 1. (zfar / (zfar - znear)) - ((zfar / zfar - znear) * znear)
			// Above is conversion to screen space apparently
			// * Projection to Screen Space -> Perspective Divide -> "Image Space" or NDC Normalized Device Coordinates
			// ** Vertex normals are vectors not points, so we don't need to project them

			for (int j = 0; j < 3; j++)
			{
				// Project the current vertex taking it into CLIP SPACE
				projected_points[j] = mat4x4_mul_vec4_project(projection_matrix, vec3_from_vec4(triangle_after_clipping.projected_points[j]));

				// Perform perspective divide with original z-value that is now stored in w
				// Point is now in Image Space or Normalized Device Space (NDC)
				if (projected_points[j].w != 0.0)
				{
					projected_points[j].x /= projected_points[j].w;
					projected_points[j].y /= projected_points[j].w;
					projected_points[j].z /= projected_points[j].w;
				}

				// Invert the y values to account for flipped screen y coordinate (color buffer grows from top to bottom)
				projected_points[j].y *= -1;
				// Invert the x values to fix difference in local mapping in blender? (Models appear flipped on the X)
				projected_points[j].x *= -1;

#pragma region PROJECTED_POINTS_INVERSION

#pragma endregion

				// Scale into the view
				projected_points[j].x *= (float)(window_width >> 1);
				projected_points[j].y *= (float)(window_height >> 1);

				// Translate the projected points to the middle of the screen
				projected_points[j].x += (float)(window_width >> 1);
				projected_points[j].y += (float)(window_height >> 1);
			}

			triangle_t triangle_to_render =
			{
				.projected_points =
				{
					{ projected_points[0].x, projected_points[0].y, projected_points[0].z, projected_points[0].w },
					{ projected_points[1].x, projected_points[1].y, projected_points[1].z, projected_points[1].w },
					{ projected_points[2].x, projected_points[2].y, projected_points[2].z, projected_points[2].w },
				},
				.texcoords =
				{
					triangle_after_clipping.texcoords[0],
					triangle_after_clipping.texcoords[1],
					triangle_after_clipping.texcoords[2]
				},
				.vertex_normals =
				{
					transformed_normals[0],
					transformed_normals[1],
					transformed_normals[2]
				},
				.vertex_colors =
				{
					vertex_colors[0],
					vertex_colors[1],
					vertex_colors[2]
				}
			};

			if (current_lighting_method == FLAT_SHADING)
			{
				triangle_to_render.vertex_normals[0] = face_normal;
				triangle_to_render.vertex_normals[1] = face_normal;
				triangle_to_render.vertex_normals[2] = face_normal;
			}

			// Save the projected triangle in the array of triangles to render
			array_push(triangles_to_render, triangle_to_render);
		}
	}

	if (measure_time)
	{
		StopTimer("Update Time");
	}


	//if (measure_time)
	//{
	//	StartTimer();

	//	vec2_t result;
	//	for (int i = 0; i < 10000000; i++)
	//	{
	//		vec4_t in = { i, i, i, i };
	//		result = vec2_from_vec4(in);
	//	}

	//	StopTimer("Return no alloc");

	//	StartTimer();

	//	for (int i = 0; i < 10000000; i++)
	//	{
	//		vec4_t in = { i, i, i, i };
	//		result = vec2_from_vec4_alloc(in);
	//	}

	//	StopTimer("Return with alloc");
	//}
}

void render(void)
{
	if (measure_time)
	{
		StartTimer();
	}

	clear_color_buffer(BACKGROUND_COLOR);
	clear_z_buffer();

	//draw_grid_threaded(grid_size);

	draw_a_thing();

	// Loop all projected triangles and render them
	int num_triangles = array_length(triangles_to_render);

	for (int i = 0; i < num_triangles; i++)
	{
		triangle_t triangle = triangles_to_render[i];

		if (show_solid)
		{
			// LIGHTING

			// Calculate the triangle color base on the light angle ala intensity factor
			uint32_t triangle_color;

			switch (current_lighting_method)
			{
			case NO_LIGHTING:
			{
				if (show_textured)
				{
					draw_textured_triangle
					(
						triangle.projected_points[0].x, triangle.projected_points[0].y, triangle.projected_points[0].w, triangle.texcoords[0].u, triangle.texcoords[0].v,
						triangle.projected_points[1].x, triangle.projected_points[1].y, triangle.projected_points[1].w, triangle.texcoords[1].u, triangle.texcoords[1].v,
						triangle.projected_points[2].x, triangle.projected_points[2].y, triangle.projected_points[2].w, triangle.texcoords[2].u, triangle.texcoords[2].v,
						mesh_texture
					);
				}
				else
				{
					triangle_color = triangle.vertex_colors[0];

					draw_filled_triangle_flat
					(
						triangle.projected_points[0].x, triangle.projected_points[0].y, triangle.projected_points[0].w,
						triangle.projected_points[1].x, triangle.projected_points[1].y, triangle.projected_points[1].w,
						triangle.projected_points[2].x, triangle.projected_points[2].y, triangle.projected_points[2].w,
						triangle_color
					);
				}

				break;
			}
			case FLAT_SHADING:
			{
				if (show_textured)
				{
					// GET TEXTURES WORKING WITH FLAT SHADING
				}
				else
				{
					// Calculate Light Intensity Per Face
					//vec3_t avg_vertex_normal = vec3_add(vec3_add(triangle.vertex_normals[0], triangle.vertex_normals[1]), triangle.vertex_normals[2]);
					//avg_vertex_normal = vec3_div(avg_vertex_normal, 3);
					float light_intensity_factor = -vec3_dot(triangle.vertex_normals[0], directional_light.direction);

					triangle_color = light_apply_intensity(triangle.vertex_colors[0], light_intensity_factor);

					draw_filled_triangle_flat
					(
						triangle.projected_points[0].x, triangle.projected_points[0].y, triangle.projected_points[0].w,
						triangle.projected_points[1].x, triangle.projected_points[1].y, triangle.projected_points[1].w,
						triangle.projected_points[2].x, triangle.projected_points[2].y, triangle.projected_points[2].w,
						triangle_color
					);
				}
				break;
			}
			case GOURAUD_SHADING:
			{
				if (show_textured)
				{
					// GET TEXTURES WORKING WITH GOURAUD SHADING
				}
				else
				{
					// Calculate Light Intensity Per Vertex
					uint32_t light_intensity_per_vertex[3];

					for (int i = 0; i < 3; i++)
					{
						float light_intensity_factor = -vec3_dot(triangle.vertex_normals[i], directional_light.direction);

						light_intensity_per_vertex[i] = light_apply_intensity(triangle.vertex_colors[i], light_intensity_factor);
					}

					draw_filled_triangle_gouraud
					(
						triangle.projected_points[0].x, triangle.projected_points[0].y,
						triangle.projected_points[1].x, triangle.projected_points[1].y,
						triangle.projected_points[2].x, triangle.projected_points[2].y,
						light_intensity_per_vertex
					);
				}
				break;
			}
			case PHONG_SHADING:
			{
				if (show_textured)
				{
					// GET TEXTURES WORKING WITH PHONG SHADING
				}
				else
				{
					// Calculate Light Intensity Per Pixel (Lerp Normals and Then Calc Light Intensity Per Pixel)

					// TODO: ADD PERSPECTIVE CORRECTION TO THE NORMALS SO THEY'RE LERPED PROPERLY, Should fix out strange looking normal lines

					draw_filled_triangle_phong
					(
						triangle.projected_points[0].x, triangle.projected_points[0].y,
						triangle.projected_points[1].x, triangle.projected_points[1].y,
						triangle.projected_points[2].x, triangle.projected_points[2].y,
						triangle.vertex_colors[0],
						triangle.vertex_normals
					);
				}
				break;
			}
			default:
			{
				triangle_color = 0xFFFF00FF;
				break;
			}
			}
		}

		if (show_wireframe)
		{
			draw_triangle
			(
				triangle.projected_points[0].x, triangle.projected_points[0].y,
				triangle.projected_points[1].x, triangle.projected_points[1].y,
				triangle.projected_points[2].x, triangle.projected_points[2].y,
				0xFFF1F1F1
			);
		}

		if (show_vertices)
		{
			for (int i = 0; i < 3; i++)
			{
				draw_rect(triangle.projected_points[i].x, triangle.projected_points[i].y, 1, 1, 0xFFFF0000);
				//draw_circle(triangle.points[i].x, triangle.points[i].y, 2, 0xFFDD2222);
				//draw_pixel(triangle.points[i].x, triangle.points[i].y, 0xFFEE8800);
			}
		}

		if (show_normals)
		{
			vec2_t face_position = { 0,0 };
			for (int i = 0; i < 3; i++)
			{
				draw_line(triangle.projected_points[i].x, triangle.projected_points[i].y, triangle.projected_points[i].x - triangle.vertex_normals[0].x * 40, triangle.projected_points[i].y - triangle.vertex_normals[0].y * 40, 0xFF0000FF);
				face_position = vec2_add(face_position, vec2_from_vec4(triangle.projected_points[i]));
			}
			face_position = vec2_div(face_position, 3);

			draw_line(face_position.x, face_position.y, face_position.x - triangle.vertex_normals[0].x * 40, face_position.y - triangle.vertex_normals[0].y * 40, 0xFFFF00FF);
		}
	}

	if (show_depth)
	{
		display_depth();
	}

	// Clear the array of triangles
	array_reset(triangles_to_render);

	render_color_buffer();

	if (measure_time)
	{
		StopTimer("Render Time");
		measure_time = false;
	}
}

// FREE THE MEMORY THAT WAS DYNAMICALLY ALLOCATED BY THE PROGRAM
void free_resources(void)
{
	free(threadIDs);

	// Clear our arrays
	array_free(mesh.faces);
	array_free(mesh.vertices);
	upng_free(png_texture);
}

int main(int argc, char* argv[])
{
	StartTimer();

	is_running = initialize_window();
	setup();

	StopTimer("Setup Time");

	while (is_running)
	{
		process_input();
		update();
		render();
	}

	destroy_window();
	free_resources();

	return 0;
}
