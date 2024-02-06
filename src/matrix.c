#include "matrix.h"

mat3x4_t world_matrix;
mat3x4_t world_view_matrix;
mat4x4_t projection_matrix;
mat3x4_t view_matrix;
//mat4x4_t orthographic_matrix;

mat3x3_t mat3x3_identity(void)
{
	mat3x3_t m =
	{ {
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1}
	} };

	return m;
}

mat4x4_t mat4x4_identity(void)
{
	mat4x4_t m =
	{ {
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1}
	} };

	return m;
}

mat3x4_t mat3x4_identity(void)
{
	mat3x4_t m =
	{ {
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, 1, 0}
	} };

	return m;
}

mat4x4_t mat4x4_make_scale(float sx, float sy, float sz)
{
	mat4x4_t m =
	{ {
		{sx, 0, 0, 0},
		{0, sy, 0, 0},
		{0, 0, sz, 0},
		{0, 0, 0, 1}
	} };

	return m;
}

mat3x4_t mat3x4_make_scale(float sx, float sy, float sz)
{
	mat3x4_t m =
	{ {
		{sx, 0, 0, 0},
		{0, sy, 0, 0},
		{0, 0, sz, 0}
	} };

	return m;
}

mat4x4_t mat4x4_make_translation(float tx, float ty, float tz)
{
	mat4x4_t m =
	{ {
		{1, 0, 0, tx},
		{0, 1, 0, ty},
		{0, 0, 1, tz},
		{0, 0, 0, 1}
	} };

	return m;
}

mat3x4_t mat3x4_make_translation(float tx, float ty, float tz)
{
	mat3x4_t m =
	{ {
		{1, 0, 0, tx},
		{0, 1, 0, ty},
		{0, 0, 1, tz}
	} };

	return m;
}

// ROTATIONS BASED ON LEFT HANDED COORD SYSTEM (+Z into the screen)
mat4x4_t mat4x4_make_rotation_x(float angle, bool clockwise)
{
	int cw = clockwise ? -1 : 1;
	float c = cos(angle);
	float s = sin(angle) * cw;

	mat4x4_t m =
	{ {
		{1, 0, 0, 0},
		{0, c, s, 0},
		{0, -s, c, 0},
		{0, 0, 0, 1}
	} };

	return m;
}

mat3x4_t mat3x4_make_rotation_x(float angle, bool clockwise)
{
	int cw = clockwise ? -1 : 1;
	float c = cos(angle);
	float s = sin(angle) * cw;

	mat3x4_t m =
	{ {
		{1, 0, 0, 0},
		{0, c, s, 0},
		{0, -s, c, 0}
	} };

	return m;
}

mat4x4_t mat4x4_make_rotation_y(float angle, bool clockwise)
{
	int cw = clockwise ? -1 : 1;
	float c = cos(angle);
	float s = sin(angle) * cw;

	mat4x4_t m =
	{ {
		{c, 0, -s, 0},
		{0, 1, 0, 0},
		{s, 0, c, 0},
		{0, 0, 0, 1}
	} };

	return m;
}

mat3x4_t mat3x4_make_rotation_y(float angle, bool clockwise)
{
	int cw = clockwise ? -1 : 1;
	float c = cos(angle);
	float s = sin(angle) * cw;

	mat3x4_t m =
	{ {
		{c, 0, -s, 0},
		{0, 1, 0, 0},
		{s, 0, c, 0}
	} };

	return m;
}

mat4x4_t mat4x4_make_rotation_z(float angle, bool clockwise)
{
	int cw = clockwise ? -1 : 1;
	float c = cos(angle);
	float s = sin(angle) * cw;

	mat4x4_t m =
	{ {
		{c, s, 0, 0},
		{-s, c, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1}
	} };

	return m;
}

mat3x4_t mat3x4_make_rotation_z(float angle, bool clockwise)
{
	int cw = clockwise ? -1 : 1;
	float c = cos(angle);
	float s = sin(angle) * cw;

	mat3x4_t m =
	{ {
		{c, s, 0, 0},
		{-s, c, 0, 0},
		{0, 0, 1, 0}
	} };

	return m;
}

mat4x4_t mat4x4_make_look_at(vec3_t eye, vec3_t target, vec3_t up)
{
	vec3_t z = vec3_sub(target, eye);
	vec3_normalize(&z);
	vec3_t x = vec3_cross(up, z);
	vec3_normalize(&x);
	vec3_t y = vec3_cross(z, x);

	mat4x4_t view_matrix =
	{ {
		{x.x, x.y, x.z, -vec3_dot(x,eye)},
		{y.x, y.y, y.z, -vec3_dot(y, eye)},
		{z.x, z.y, z.z, -vec3_dot(z, eye)},
		{0, 0, 0, 1}
	} };
	return view_matrix;
}

mat3x4_t mat3x4_make_look_at(vec3_t eye, vec3_t target, vec3_t up)
{
	vec3_t z = vec3_sub(target, eye);
	vec3_normalize(&z);
	vec3_t x = vec3_cross(up, z);
	vec3_normalize(&x);
	vec3_t y = vec3_cross(z, x);

	mat3x4_t view_matrix =
	{ {
		{x.x, x.y, x.z, -vec3_dot(x,eye)},
		{y.x, y.y, y.z, -vec3_dot(y, eye)},
		{z.x, z.y, z.z, -vec3_dot(z, eye)}
	} };
	return view_matrix;
}

mat4x4_t mat4x4_make_perspective(float fov, float aspect, float znear, float zfar)
{
	float fovCalc = 1 / tan(fov * 0.5);
	float aspectCalc = aspect * fovCalc;
	float zFarCalc = zfar / (zfar - znear);
	float zNearCalc = zFarCalc * znear;

	mat4x4_t m =
	{ {
		{aspectCalc, 0, 0, 0},
		{0, fovCalc, 0, 0},
		{0, 0, zFarCalc, zNearCalc},
		{0, 0, 1, 0}
	} };

	return m;
}

mat4x4_t mat4x4_make_orthographic(float window_width, float window_height, float orthographic_size, float znear, float zfar)
{
	float windowRatio = window_width / (float)window_height;

	float half = orthographic_size * 0.5;

	float right = half * windowRatio;
	float left = -right;
	float top = half;
	float bottom = -top;

	float x = 2 / (right - left);
	float y = 2 / (top - bottom);
	float z = -2 / (zfar - znear);
	float xw = -(right + left) / (right - left);
	float yw = -(top + bottom) / (top - bottom);
	float zw = -(zfar + znear) / (zfar - znear);

	mat4x4_t m =
	{ {
		{x, 0, 0, xw},
		{0, y, 0, yw},
		{0, 0, z, zw},
		{0, 0, 0, 1}
	} };

	return m;
}

vec4_t mat4x4_mul_vec4(mat4x4_t m, vec4_t v)
{
	vec4_t result;

	result.x = m.m[0][0] * v.x + m.m[0][1] * v.y + m.m[0][2] * v.z + m.m[0][3] * v.w;
	result.y = m.m[1][0] * v.x + m.m[1][1] * v.y + m.m[1][2] * v.z + m.m[1][3] * v.w;
	result.z = m.m[2][0] * v.x + m.m[2][1] * v.y + m.m[2][2] * v.z + m.m[2][3] * v.w;
	result.w = m.m[3][0] * v.x + m.m[3][1] * v.y + m.m[3][2] * v.z + m.m[3][3] * v.w;

	return result;
}


vec3_t mat4x4_mul_vec3(mat4x4_t m, vec3_t v)
{
	vec3_t result;

	result.x = m.m[0][0] * v.x + m.m[0][1] * v.y + m.m[0][2] * v.z;
	result.y = m.m[1][0] * v.x + m.m[1][1] * v.y + m.m[1][2] * v.z;
	result.z = m.m[2][0] * v.x + m.m[2][1] * v.y + m.m[2][2] * v.z;

	return result;
}

vec3_t mat3x4_mul_vec3(mat3x4_t m, vec3_t v)
{
	vec3_t result;

	result.x = m.m[0][0] * v.x + m.m[0][1] * v.y + m.m[0][2] * v.z + m.m[0][3];
	result.y = m.m[1][0] * v.x + m.m[1][1] * v.y + m.m[1][2] * v.z + m.m[1][3];
	result.z = m.m[2][0] * v.x + m.m[2][1] * v.y + m.m[2][2] * v.z + m.m[2][3];

	return result;
}

vec4_t mat4x4_mul_vec4_project(mat4x4_t projection_matrix, vec3_t v)
{
	// Slot vec3 into a vec4
	vec4_t vec = { v.x, v.y, v.z, 1 };

	// multiply the projection matrix by our original vector
	vec4_t result = mat4x4_mul_vec4(projection_matrix, vec);

	return result;
}

mat4x4_t mat4x4_mul_mat4x4(mat4x4_t a, mat4x4_t b)
{
	mat4x4_t result;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			// There is an accumulation method that might be more cache efficient?
			result.m[i][j] =
				a.m[i][0] * b.m[0][j] +
				a.m[i][1] * b.m[1][j] +
				a.m[i][2] * b.m[2][j] +
				a.m[i][3] * b.m[3][j];
		}
	}

	return result;
}

mat3x4_t mat3x4_mul_mat3x4(mat3x4_t a, mat3x4_t b)
{
	mat3x4_t result;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			result.m[i][j] = 0.0;
			for (int k = 0; k < 3; k++) {
				result.m[i][j] += a.m[i][k] * b.m[k][j]; //Accumulation method?
			}
		}
		result.m[i][3] = a.m[i][0] * b.m[0][3] + a.m[i][1] * b.m[1][3] + a.m[i][2] * b.m[2][3] + a.m[i][3];
	}
	// Add translation component
	//for (int i = 0; i < 3; i++) {
	// result.m[i][3] = a.m[i][0] * b.m[0][3] + a.m[i][1] * b.m[1][3] + a.m[i][2] * b.m[2][3] + a.m[i][3];
	//}

	return result;
}

void print_mat4x4(mat4x4_t matrix)
{
	printf("Resultant Matrix:\n");
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			printf("%f\t", matrix.m[i][j]);
		}
		printf("\n");
	}
}

void print_mat3x4(mat3x4_t matrix)
{
	printf("Resultant Matrix:\n");
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			printf("%f\t", matrix.m[i][j]);
		}
		printf("\n");
	}
}

// Function to calculate the determinant of a 3x3 matrix
float determinant3x3(mat3x3_t mat)
{
	return mat.m[0][0] * (mat.m[1][1] * mat.m[2][2] - mat.m[1][2] * mat.m[2][1]) -
		mat.m[0][1] * (mat.m[1][0] * mat.m[2][2] - mat.m[1][2] * mat.m[2][0]) +
		mat.m[0][2] * (mat.m[1][0] * mat.m[2][1] - mat.m[1][1] * mat.m[2][0]);
}

// Function to calculate the adjugate (adjoint) of a 3x3 matrix
mat3x3_t adjugate3x3(mat3x3_t mat)
{
	mat3x3_t adj =
	{
		.m[0][0] = mat.m[1][1] * mat.m[2][2] - mat.m[1][2] * mat.m[2][1],
		.m[0][1] = -mat.m[0][1] * mat.m[2][2] + mat.m[0][2] * mat.m[2][1],
		.m[0][2] = mat.m[0][1] * mat.m[1][2] - mat.m[0][2] * mat.m[1][1],
		.m[1][0] = -mat.m[1][0] * mat.m[2][2] + mat.m[1][2] * mat.m[2][0],
		.m[1][1] = mat.m[0][0] * mat.m[2][2] - mat.m[0][2] * mat.m[2][0],
		.m[1][2] = -mat.m[0][0] * mat.m[1][2] + mat.m[0][2] * mat.m[1][0],
		.m[2][0] = mat.m[1][0] * mat.m[2][1] - mat.m[1][1] * mat.m[2][0],
		.m[2][1] = -mat.m[0][0] * mat.m[2][1] + mat.m[0][1] * mat.m[2][0],
		.m[2][2] = mat.m[0][0] * mat.m[1][1] - mat.m[0][1] * mat.m[1][0]
	};

	return adj;
}

mat4x4_t adjugate4x4(mat4x4_t mat)
{
	mat4x4_t result;

	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			// Calculate the cofactor and transpose it to get the adjugate
			result.m[j][i] = cofactor(mat, i, j);
		}
	}

	return result;
}

// Function to calculate the inverse of a 3x3 matrix
mat3x3_t inverse3x3(mat3x3_t mat)
{
	float det = determinant3x3(mat);
	mat3x3_t inverse = mat3x3_identity();

	if (det != 0) {
		mat3x3_t adj = adjugate3x3(mat);

		float invDet = 1.0f / det;

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				inverse.m[i][j] = adj.m[i][j] * invDet;
			}
		}

		return inverse;
	}
	else {
		printf("Error: Matrix is singular, inverse does not exist.\n");
		return mat3x3_identity();
	}
}

float cofactor(mat4x4_t mat, int row, int col)
{
	int sign = ((row + col) & 2) == 0 ? 1 : -1;

	mat3x3_t mat3x3 =
	{ {
		{mat.m[(row + 1) % 4][(col + 1) % 4], mat.m[(row + 1) % 4][(col + 2) % 4], mat.m[(row + 1) % 4][(col + 3) % 4]},
		{mat.m[(row + 2) % 4][(col + 1) % 4], mat.m[(row + 2) % 4][(col + 2) % 4], mat.m[(row + 2) % 4][(col + 3) % 4]},
		{mat.m[(row + 3) % 4][(col + 1) % 4], mat.m[(row + 3) % 4][(col + 2) % 4], mat.m[(row + 3) % 4][(col + 3) % 4]}
	} };

	float determinant = determinant3x3(mat3x3);

	return sign * determinant;
}

// Function to calculate the inverse of a 4x4 matrix
mat4x4_t inverse4x4(mat4x4_t mat)
{
	mat4x4_t result = mat4x4_identity();

	// Calculate the determinant of the input matrix
	float det = mat.m[0][0] * cofactor(mat, 0, 0) +
		mat.m[0][1] * cofactor(mat, 0, 1) +
		mat.m[0][2] * cofactor(mat, 0, 2) +
		mat.m[0][3] * cofactor(mat, 0, 3);

	// Check if the determinant is zero (matrix is singular)
	if (det == 0.0f) {
		// Handle the error (e.g., print an error message or return an identity matrix)
		printf("Error: Singular matrix, inverse does not exist.\n");
		// Returning an identity matrix for simplicity; you might want to handle errors differently
		return result;
	}

	// Calculate the inverse by dividing the adjugate by the determinant
	float invDet = 1.0f / det;

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			result.m[i][j] = adjugate4x4(mat).m[i][j] * invDet;
		}
	}

	return result;
}

// Function to transpose a 3x3 matrix
mat3x3_t transpose3x3(mat3x3_t mat)
{
	mat3x3_t transpose = mat3x3_identity();

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			transpose.m[i][j] = mat.m[j][i];
		}
	}

	return transpose;
}

// Function to multiply a 3x3 matrix by a 3D vector
vec3_t mat3x3_mul_vec3(mat3x3_t mat, vec3_t vec)
{
	float result[3];
	float in_vec[3];
	in_vec[0] = vec.x;
	in_vec[1] = vec.y;
	in_vec[2] = vec.z;

	for (int i = 0; i < 3; i++) {
		result[i] = 0.0f;
		for (int j = 0; j < 3; j++) {
			result[i] += mat.m[i][j] * in_vec[j];
		}
	}

	vec3_t out_vec = { .x = result[0], .y = result[1], .z = result[2] };

	return out_vec;
}

mat3x3_t mat3x3_make_transpose_inverse(mat3x4_t world_matrix)
{
	return transpose3x3(inverse3x3(mat3x4_to_mat3x3(world_matrix)));
}

mat3x3_t mat3x4_to_mat3x3(mat3x4_t mat3x4)
{
	mat3x3_t mat3x3;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			mat3x3.m[i][j] = mat3x4.m[i][j];
		}
	}

	return mat3x3;
}

vec3_t local_to_world(vec3_t v)
{
	return mat3x4_mul_vec3(world_matrix, v);
}
vec3_t world_to_local(vec3_t v)
{
	return mat3x3_mul_vec3(mat3x4_to_mat3x3(world_matrix), v);
}
vec3_t screen_to_world(vec3_t v)
{
	vec3_t screen_pos_to_clip_pos = { .x = (v.x / window_width) * 2.0 - 1.0, 1.0 - (v.y / window_height) * 2.0, v.z };
	mat4x4_t inverse_projection_matrix = inverse4x4(projection_matrix);

	return mat4x4_mul_vec3(inverse_projection_matrix, screen_pos_to_clip_pos);
}
vec3_t screen_to_view(vec3_t v)
{
	vec3_t screen_pos_to_clip_pos = { .x = (v.x / window_width) * 2.0 - 1.0, 1.0 - (v.y / window_height) * 2.0, v.z };
	mat4x4_t inverse_projection_matrix = inverse4x4(projection_matrix);
	screen_pos_to_clip_pos = mat4x4_mul_vec3(inverse_projection_matrix, screen_pos_to_clip_pos);
	mat3x3_t inverse_view_matrix = inverse3x3(mat3x4_to_mat3x3(view_matrix));

	return mat3x3_mul_vec3(inverse_view_matrix, screen_pos_to_clip_pos);
}
