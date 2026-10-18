// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file normals.h
/// @brief Per-vertex normal computation (Spec 10 Step 1).
///
/// Ports `createNormals()` from `legacy/main.cpp:87-115` into pure C++20.
///
/// ## EqualUnit — the 2003 algorithm (default, golden-pinned)
///
/// ```
/// for each face (a, b, c):
///     nrm = cross(pos[b]-pos[a], pos[c]-pos[a])   // magnitude = 2·area
///     nrm = normalize(nrm)                          // discard area, equal weight
///     norm[a] += nrm; norm[b] += nrm; norm[c] += nrm
/// for each vertex v: normalize(norm[v])
/// ```
///
/// The `.norm()` before accumulation is the 2003 choice — every adjacent face
/// contributes equally regardless of area.  The golden is captured against this.
///
/// ## Area — the area-weighted alternative (NOT the golden default)
///
/// Identical except the per-face `normalize` is skipped; the 2·area magnitude
/// survives into the accumulation, so large faces dominate.  This is the
/// deliberate alternative exposed for the UI to toggle (B3 test).  It MUST NOT
/// be used in the golden or the behavior-preserving path.

#include <vector>

#include <decimate/math/vec3.h>
#include <decimate/types/mesh_data.h>

namespace decimate::bake {

/// @brief Face normal weighting strategy for `computeVertexNormals`.
///
/// `EqualUnit` — DEFAULT; reproduces 2003 `createNormals` exactly.
///              Normalize the face normal BEFORE accumulation.
///              The golden is captured against this.
///
/// `Area`      — Area-weighted alternative; skip the per-face normalize.
///              Tested separately (B3). MUST NOT replace EqualUnit in the
///              behavior-preserving pipeline (anti-pattern 11 / BLOCK).
enum class Weighting {
    EqualUnit,  ///< Equal-weight (2003 legacy default).
    Area,       ///< Area-weighted (alternative mode, not the golden path).
};

/// @brief Compute per-vertex normals for a triangle mesh.
///
/// @param data   Source mesh: `positions` (3 floats/vertex, packed) and
///               `indices` (3 ints/face, packed, 0-based).
/// @param w      Weighting strategy (default `EqualUnit` — the golden path).
/// @returns      One `Vec3<float>` per vertex (size == positions.size()).
///               Zero-length normals are possible for isolated/degenerate
///               vertices; consumers that require unit normals must normalize.
///
/// @pre  `data.indices.size() % 3 == 0`.
/// @pre  Every index in `data.indices` is in `[0, data.positions.size())`.
[[nodiscard]] std::vector<decimate::math::Vec3<float>> computeVertexNormals(
    const decimate::MeshData& data,
    Weighting             w = Weighting::EqualUnit);

}  // namespace decimate::bake
