// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file textured_mesh.h
/// @brief A triangle mesh with a separate UV index space (the OBJ `v`/`vt` model).
///
/// `MeshData` carries positions and indices only; per-corner UVs are a bare
/// `std::vector<Vec2f>` elsewhere.  `TexturedMesh` is the value type that ties the
/// two together with *topology*: faces index positions through `faceV` and UV
/// vertices through `faceU`.  Two faces that meet along a surface edge but reference
/// different UV vertices there form a texture seam (Liu et al. 2017, section 3).
///
/// Header-only and dependency-free (types + math only), so any module can take a
/// `TexturedMesh` without linking the seam library.

#include <cstddef>
#include <decimate/math/vec2.h>
#include <decimate/math/vec3.h>
#include <vector>

namespace decimate {

/// @brief Positions + UVs with independent index lists.
///
/// ## Invariants
///
/// - `faceV.size() == faceU.size()`, a multiple of 3 (triangles).
/// - `faceV[i]` indexes `positions`; `faceU[i]` indexes `uvs`.
/// - Corner `k` of face `f` is `faceV[3*f+k]` / `faceU[3*f+k]`; the face order is the
///   same in both lists.
/// - UVs are in texture space `[0,1]^2`, `v` growing with the `NormalMap` row index.
struct TexturedMesh {
    std::vector<math::Vec3f> positions; ///< Vertex positions, indexed by faceV.
    std::vector<math::Vec2f> uvs;       ///< UV vertices, indexed by faceU.
    std::vector<int>         faceV;     ///< 3 position indices per face.
    std::vector<int>         faceU;     ///< 3 UV indices per face, parallel to faceV.

    /// @brief Number of triangles.
    [[nodiscard]] std::size_t faceCount() const noexcept {
        return faceV.size() / 3;
    }

    /// @brief The 3F per-corner UV vector that `bake::bake()` consumes.
    ///
    /// Entry `3*f+k` is the UV of corner `k` of face `f`.
    ///
    /// @pre every `faceU` entry is in `[0, uvs.size())` (unchecked).
    [[nodiscard]] std::vector<math::Vec2f> cornerUV() const {
        std::vector<math::Vec2f> out;
        out.reserve(faceU.size());
        for (const int u : faceU) {
            out.push_back(uvs[static_cast<std::size_t>(u)]);
        }
        return out;
    }
};

} // namespace decimate
