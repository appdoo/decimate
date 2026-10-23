// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file tangents.h
/// @brief Per-corner MikkTSpace tangents - the basis a tangent-space normal map is baked and decoded in.

#include <vector>

#include <decimate/math/vec2.h>
#include <decimate/math/vec3.h>
#include <decimate/math/vec4.h>
#include <decimate/types/mesh_data.h>

namespace decimate {

/// @brief MikkTSpace tangents of a UV-mapped mesh (Mikkelsen 2008, the basis of Blender, Unity, Unreal, glTF).
///
/// Decode with the pixel-shader convention every engine uses: interpolate T and N across the triangle
/// WITHOUT normalising, `B = w * cross(N, T)`, `n = normalize(tn.x*T + tn.y*B + tn.z*N)`.
///
/// @param mesh          Triangle mesh; face `f` has corners `3*f .. 3*f+2`.
/// @param cornerUV      Per-corner UV, size `3 * faceCount` (the atlas `uv`).
/// @param vertexNormals Per-vertex normals, size `mesh.positions.size()`.
/// @returns One tangent per corner: xyz = unit T, w = +1 or -1 (bitangent handedness). Corners of
///          degenerate faces (zero UV or position area) get (0,0,0,1): no basis, consumers fall back.
///
/// @pre `cornerUV.size() == mesh.indices.size()` and `vertexNormals.size() == mesh.positions.size()`.
[[nodiscard]] std::vector<math::Vec4f> tangents(const MeshData                 &mesh,
                                                const std::vector<math::Vec2f> &cornerUV,
                                                const std::vector<math::Vec3f> &vertexNormals);

} // namespace decimate
