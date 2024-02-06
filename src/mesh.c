#include "mesh.h"

mesh_t mesh =
{
	.vertices = NULL,
	.normals = NULL,
	.uvs = NULL,
	.faces = NULL,
	.rotation = {0, 0, 0},
	.scale = { 1, 1, 1 },
	.translation = { 0, 0, 5},
	.vertex_colors = NULL
};

vec3_t cube_vertices[N_CUBE_VERTICES] =
{
	{-1,1,-1}, //1
	{-1,-1,-1},	//2
	{-1,1,1},	//3
	{-1,-1,1},	//4
	{1,1,-1},	//5
	{1,-1,-1},	//6
	{1,1,1},	//7
	{1,-1,1}	//8
};

tex2_t cube_uvs[N_CUBE_VERTICES] =
{
	{0,0},
	{0,1},
	{1,0},
	{1,1},
	{0,0},
	{0,0},
	{0,0},
	{0,0}
};

uint32_t cube_vertex_colors[N_CUBE_VERTICES] =
{
	//Right Side
	0xFFFF00FF,
	0xFFFF00FF,
	0xFFFF00FF,
	0xFFFF00FF,
	0xFFFFFF00,
	0xFFFFFF00,
	0xFFFFFF00,
	0xFFFFFF00
};

//.OBJ Flat shaded normals output 1 normal per face, Smooth shaded outputs 1 normal per vertex
vec3_t cube_normals[N_CUBE_VERTICES] =
{
	{-1,0,0},
	{0,0,1},
	{1,0,0},
	{0,0,-1},
	{0,1,0},
	{0,-1,0},
	{0,0,0},
	{0,0,0}
};

// Face data is the indexes into the 3 arrays for vertices, uvs and normals.
//In that order so 2/4/1 3/8/1 1/1/1 means v[2]/uv[4]/n[1] v[3]/uv[8]/n[1] v[1]/uv[1]/n[1]
face_t cube_faces[N_CUBE_FACES] =
{
	// Right
	{.a = 2, .uv_a = 2, .normal_a = 1, .color_a = 2, .b = 3, .uv_b = 3, .normal_b = 1, .color_b = 3, .c = 1, .uv_c = 1, .normal_c = 1, .color_c = 1},
	{.a = 2, .uv_a = 2, .normal_a = 1, .color_a = 2, .b = 4, .uv_b = 4, .normal_b = 1, .color_b = 4, .c = 3, .uv_c = 3, .normal_c = 1, .color_c = 3},

	// Back
	{.a = 4, .uv_a = 2, .normal_a = 2, .color_a = 4, .b = 7, .uv_b = 3, .normal_b = 2, .color_b = 7, .c = 3, .uv_c = 1, .normal_c = 2, .color_c = 3},
	{.a = 4, .uv_a = 2, .normal_a = 2, .color_a = 4, .b = 8, .uv_b = 4, .normal_b = 2, .color_b = 8, .c = 7, .uv_c = 3, .normal_c = 2, .color_c = 7},

	// Left
	{.a = 8, .uv_a = 2, .normal_a = 3, .color_a = 8, .b = 5, .uv_b = 3, .normal_b = 3, .color_b = 5, .c = 7, .uv_c = 1, .normal_c = 3, .color_c = 7},
	{.a = 8, .uv_a = 2, .normal_a = 3, .color_a = 8, .b = 6, .uv_b = 4, .normal_b = 3, .color_b = 6, .c = 5, .uv_c = 3, .normal_c = 3, .color_c = 5},

	// Front
	{.a = 6, .uv_a = 2, .normal_a = 4, .color_a = 6, .b = 1, .uv_b = 3, .normal_b = 4, .color_b = 1, .c = 5, .uv_c = 1, .normal_c = 4, .color_c = 5},
	{.a = 6, .uv_a = 2, .normal_a = 4, .color_a = 6, .b = 2, .uv_b = 4, .normal_b = 4, .color_b = 2, .c = 1, .uv_c = 3, .normal_c = 4, .color_c = 1},

	// Top
	{.a = 7, .uv_a = 2, .normal_a = 5, .color_a = 7, .b = 1, .uv_b = 3, .normal_b = 5, .color_b = 1, .c = 3, .uv_c = 1, .normal_c = 5, .color_c = 3},
	{.a = 7, .uv_a = 2, .normal_a = 5, .color_a = 7, .b = 5, .uv_b = 4, .normal_b = 5, .color_b = 5, .c = 1, .uv_c = 3, .normal_c = 5, .color_c = 1},

	// Bottom
	{.a = 4, .uv_a = 2, .normal_a = 6, .color_a = 4, .b = 6, .uv_b = 3, .normal_b = 6, .color_b = 6, .c = 8, .uv_c = 1, .normal_c = 6, .color_c = 8},
	{.a = 4, .uv_a = 2, .normal_a = 6, .color_a = 4, .b = 2, .uv_b = 4, .normal_b = 6, .color_b = 2, .c = 6, .uv_c = 3, .normal_c = 6, .color_c = 6},
};

void load_cube_mesh_data(void)
{
	for (int i = 0; i < N_CUBE_VERTICES; i++)
	{
		vec3_t cube_vertex = cube_vertices[i];
		array_push(mesh.vertices, cube_vertex);

		tex2_t cube_uv = cube_uvs[i];
		array_push(mesh.uvs, cube_uv);

		vec3_t cube_normal = cube_normals[i];
		array_push(mesh.normals, cube_normal);

		uint32_t vertex_color = cube_vertex_colors[i];
		array_push(mesh.vertex_colors, vertex_color);
	}
	for (int i = 0; i < N_CUBE_FACES; i++)
	{
		face_t cube_face = cube_faces[i];
		//decrement the indexs due to obj being 1 indexed
		cube_face.a -= 1;
		cube_face.b -= 1;
		cube_face.c -= 1;
		cube_face.color_a -= 1;
		cube_face.color_b -= 1;
		cube_face.color_c -= 1;
		cube_face.normal_a -= 1;
		cube_face.normal_b -= 1;
		cube_face.normal_c -= 1;
		cube_face.uv_a -= 1;
		cube_face.uv_b -= 1;
		cube_face.uv_c -= 1;
		array_push(mesh.faces, cube_face);
	}
}

void load_obj_file_data(char* filename)
{
	FILE* file;
	file = fopen(filename, "r");

	char line[1024];

	while (fgets(line, 1024, file))
	{
		// Vertex information
		if (strncmp(line, "v ", 2) == 0)
		{
			vec3_t vertex;
			float red, green, blue;

			sscanf(line, "v %f %f %f %f %f %f", &vertex.x, &vertex.y, &vertex.z, &red, &green, &blue);

			array_push(mesh.vertices, vertex);

			uint32_t redByte = (red * 255);
			uint32_t greenByte = (green * 255);
			uint32_t blueByte = (blue * 255);

			uint32_t alpha = 0xFF;
			uint32_t vertex_color = alpha << 24 | redByte << 16 | (greenByte << 8) | blueByte;
			array_push(mesh.vertex_colors, vertex_color);
		}
		// Normal information starts with vn and provides xyz normals as signed floats
		if (strncmp(line, "vn ", 3) == 0)
		{
			vec3_t normal;

			sscanf(line, "vn %f %f %f", &normal.x, &normal.y, &normal.z);

			array_push(mesh.normals, normal);
		}
		// UV information starts with vt and provides uv coordinates in 0-1 space
		if (strncmp(line, "vt ", 3) == 0)
		{
			tex2_t texcoord;

			sscanf(line, "vt %f %f", &texcoord.u, &texcoord.v);

			// Flip/Remap the V component to account for invertex V coordinates from blender (V grows downwards), 
			// Flip/Remap the U component because we've flipped the projected X due to blender meshes coming in flipped horizontally.
			//texcoord.u = 1.0 - texcoord.u;
			texcoord.v = 1.0 - texcoord.v;

			array_push(mesh.uvs, texcoord);
		}
		// Face information starts with f. We assume they're triangles so we expect 3 sets of vertex data.
		// Vertex data is indexes into the above 3 arrays. Order is Vertex/UV/Normal
		if (strncmp(line, "f ", 2) == 0)
		{
			int vertex_indices[3];
			int texture_indices[3];
			int normal_indices[3];

			sscanf(line, "f %d/%d/%d %d/%d/%d %d/%d/%d",
				&vertex_indices[0], &texture_indices[0], &normal_indices[0],
				&vertex_indices[1], &texture_indices[1], &normal_indices[1],
				&vertex_indices[2], &texture_indices[2], &normal_indices[2]
			);
			//Construct each face
			face_t face =
			{
				//-1 because .obj indexing starts at 1
				.a = vertex_indices[0] - 1, .b = vertex_indices[1] - 1, .c = vertex_indices[2] - 1,
				.normal_a = normal_indices[0] - 1, .normal_b = normal_indices[1] - 1, .normal_c = normal_indices[2] - 1,
				.uv_a = texture_indices[0] - 1, .uv_b = texture_indices[1] - 1, .uv_c = texture_indices[2] - 1,
				.color_a = vertex_indices[0] - 1, .color_b = vertex_indices[1] - 1, .color_c = vertex_indices[2] - 1
			};

			array_push(mesh.faces, face);
		}

	}

	//Helper function to write our construction out to an .OBJ to ensure we've pulled in the original .obj correctly 
	bool write_output = false;
	if (write_output)
	{
		// WRITE OUTPUT FILE TO COMPARE WITH ORIGINAL .OBJ
		FILE* fileOut = fopen("/xOutputs/output.txt", "w");

		if (fileOut == NULL) {
			perror("Failed to open the file");
		}

		for (int i = 0; i < array_length(mesh.vertices); i++) {

			uint32_t vertex_color = mesh.vertex_colors[i];
			//uint32_t alpha = (vertex_color >> 24) & 0xFF;
			uint32_t redByte = (vertex_color >> 16) & 0xFF;
			uint32_t greenByte = (vertex_color >> 8) & 0xFF;
			uint32_t blueByte = vertex_color & 0xFF;

			float red = (float)redByte / 255.0;
			float green = (float)greenByte / 255.0;
			float blue = (float)blueByte / 255.0;

			fprintf(fileOut, "v %.6f %.6f %.6f %.4f %.4f %.4f\n", mesh.vertices[i].x, mesh.vertices[i].y, mesh.vertices[i].z, red, green, blue);
		}
		for (int i = 0; i < array_length(mesh.normals); i++) {
			fprintf(fileOut, "vn %.4f %.4f %.4f\n", mesh.normals[i].x, mesh.normals[i].y, mesh.normals[i].z);
		}
		for (int i = 0; i < array_length(mesh.uvs); i++) {
			fprintf(fileOut, "vt %.6f %.6f\n", mesh.uvs[i].u, mesh.uvs[i].v);
		}

		fprintf(fileOut, "s 1\n");

		for (int i = 0; i < array_length(mesh.faces); i++) {

			fprintf(fileOut, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", mesh.faces[i].a + 1, mesh.faces[i].uv_a + 1, mesh.faces[i].normal_a + 1, mesh.faces[i].b + 1, mesh.faces[i].uv_b + 1, mesh.faces[i].normal_b + 1, mesh.faces[i].c + 1, mesh.faces[i].uv_c + 1, mesh.faces[i].normal_c + 1);
		}

		fclose(fileOut);
	}

	fclose(file);
}
