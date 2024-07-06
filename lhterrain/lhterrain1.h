/*
Copyright (c) 2024 Ashley Rose Hale (LadyHavoc)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#pragma once

#ifndef LHTERRAIN1_H
#define LHTERRAIN1_H

#include <stdint.h>
#include <math.h>

typedef int64_t LHTerrain1_FixedPoint;

typedef struct LHTerrain1_TableInfo
{
	size_t capacity;
	size_t used;
	size_t first_unused;
	size_t after_last_used;
	size_t hash_basis;
	size_t hash_multiply;
}
LHTerrain1_TableInfo;

typedef struct LHTerrain1_Strings
{
	LHTerrain1_TableInfo info;
	uint8_t* buffer_utf8;
}
LHTerrain1_Strings;

typedef enum LHTerrain1_Flags
{
	LHTerrain1_Flag_Air    = 0x00000001,
	LHTerrain1_Flag_Opaque = 0x00000002,
	LHTerrain1_Flag_Fluid  = 0x00000004,
}
LHTerrain1_Flags;

typedef struct LHTerrain1_Materials
{
	// Name (offset in Strings buffer)
	size_t* names;
	// Flags for each material
	uint64_t* flags;
	// Pointer to some struct provided by the caller, for each material
	void* generics_pointer;
	// Number with some meaning to the caller, for each material
	int64_t* generics_number;
}
LHTerrain1_Materials;

typedef struct LHTerrain1_Points_LOD
{
	LHTerrain1_FixedPoint* x;
	LHTerrain1_FixedPoint* y;
	LHTerrain1_FixedPoint* z;
}
LHTerrain1_Points_LOD;

typedef struct LHTerrain1_Points
{
	// Metadata about the points table
	LHTerrain1_TableInfo info;
	// Bitfield representing which points are in use
	uint64_t* allocated;
	// Points data for each LOD
	LHTerrain1_Points_LOD* lod;
}
LHTerrain1_Points;

typedef struct LHTerrain1_Faces_LOD
{
	// Bitfield representing which faces exist in this LOD
	uint64_t* active;
}
LHTerrain1_Faces_LOD;

typedef struct LHTerrain1_Faces
{
	// Metadata about the faces table
	LHTerrain1_TableInfo info;
	// Bitfield representing which faces are in use
	uint64_t* allocated;
	// All 3 points used by each face uint64_t[][3]
	uint64_t* point_indices;
	// Material that is on the positive-facing side of each face
	uint64_t* material_positive;
	// Material that is on the negative-facing side of each face
	uint64_t* material_negative;
	// Faces data for each LOD
	LHTerrain1_Faces_LOD* lod;
}
LHTerrain1_Faces;

typedef struct LHTerrain1_Model
{
	/// Mesh data is stored at a fixed point integer scale as int64, this is the
	/// scaling factor from fixed point to renderable geometry (World scale)
	float fixed_point_to_world_scale;
	/// Convert from world scale to fixed point
	float fixed_point_from_world_scale;
	/// Mesh is built for this many levels of detail
	/// 
	/// How it works - each point has multiple positions (LODs), at higher LOD
	/// values some points are welded together to shared positions, causing some
	/// triangles to cease to exist (no area, due to two or more of their points
	/// being identical)
	size_t levels_of_detail;

	// Storage for strings used by the terrain information (e.g. materials)
	LHTerrain1_Strings strings;

	// Array of arrays of point positions
	LHTerrain1_Points points;

	// Array of LODs, containing bitfields representing which faces are active
	LHTerrain1_Faces faces;
}
LHTerrain1_Model;

static inline float LHTerrain1_Units_ToVertex(LHTerrain1_Model* model, LHTerrain1_FixedPoint x)
{
	return x * model->fixed_point_to_world_scale;
}
static inline LHTerrain1_FixedPoint LHTerrain1_Units_FromVertex(LHTerrain1_Model *model, float x)
{
	return (int64_t)rint(x * model->fixed_point_from_world_scale);
}

#endif
