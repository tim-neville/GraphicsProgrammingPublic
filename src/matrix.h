#pragma once

#include "vector.h"
#include <stdbool.h>
#include <stdio.h>
#include "display.h"

typedef struct
{
	float m[3][3];
} mat3x3_t;

typedef struct 
{
	float m[4][4];
} mat4x4_t;

typedef struct 
{
	float m[3][4];
} mat3x4_t;

extern mat3x4_t world_matrix;
extern mat4x4_t projection_matrix;
extern mat3x4_t view_matrix;
extern mat3x4_t world_view_matrix;
//mat4x4_t orthographic_matrix;

mat3x3_t mat3x3_identity(void);

mat4x4_t mat4x4_identity(void);
mat4x4_t mat4x4_make_scale(float sx, float sy, float sz);
mat4x4_t mat4x4_make_translation(float tx, float ty, float tz);
mat4x4_t mat4x4_make_rotation_x(float angle, bool clockwise);
mat4x4_t mat4x4_make_rotation_y(float angle, bool clockwise);
mat4x4_t mat4x4_make_rotation_z(float angle, bool clockwise);
mat4x4_t mat4x4_make_look_at(vec3_t eye, vec3_t target, vec3_t up);

mat3x4_t mat3x4_identity(void);
mat3x4_t mat3x4_make_scale(float sx, float sy, float sz);
mat3x4_t mat3x4_make_translation(float tx, float ty, float tz);
mat3x4_t mat3x4_make_rotation_x(float angle, bool clockwise);
mat3x4_t mat3x4_make_rotation_y(float angle, bool clockwise);
mat3x4_t mat3x4_make_rotation_z(float angle, bool clockwise);
mat3x4_t mat3x4_make_look_at(vec3_t eye, vec3_t target, vec3_t up);

mat4x4_t mat4x4_make_perspective(float fov, float aspect, float znear, float zfar);
mat4x4_t mat4x4_make_orthographic(float window_width, float window_height, float orthographic_size, float znear, float zfar);

vec4_t mat4x4_mul_vec4(mat4x4_t m, vec4_t v);
vec3_t mat3x4_mul_vec3(mat3x4_t m, vec3_t v);
vec4_t mat4x4_mul_vec4_project(mat4x4_t projection_matrix, vec3_t v);

mat4x4_t mat4x4_mul_mat4x4(mat4x4_t a, mat4x4_t b);
mat3x4_t mat3x4_mul_mat3x4(mat3x4_t a, mat3x4_t b);

void print_mat4x4(mat4x4_t matrix);
void print_mat3x4(mat3x4_t matrix);

float determinant3x3(mat3x3_t mat);
mat3x3_t adjugate3x3(mat3x3_t mat);
mat4x4_t adjugate4x4(mat4x4_t mat);
mat3x3_t inverse3x3(mat3x3_t mat);
float cofactor(mat4x4_t mat, int row, int col);
mat4x4_t inverse4x4(mat4x4_t mat);
mat3x3_t transpose3x3(mat3x3_t mat);
vec3_t mat3x3_mul_vec3(mat3x3_t matrix, vec3_t vec);
mat3x3_t mat3x3_make_transpose_inverse(mat3x4_t world_matrix);

mat3x3_t mat3x4_to_mat3x3(mat3x4_t mat3x4);

vec3_t local_to_world(vec3_t v);
vec3_t world_to_local(vec3_t v);
vec3_t screen_to_world(vec3_t v);
vec3_t screen_to_view(vec3_t v);
