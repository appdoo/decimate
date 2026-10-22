// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file mesh_types.h
/// @brief Minimal geometry types shared across decimate modules.
///
/// Only the types the core actually consumes live here. The heavy legacy
/// math3D.h header is deliberately NOT pulled in; mat3/mat4 and rendering
/// types belong to their consumer specs (08 chartification, etc.) and are
/// added when a real consumer needs them — not speculatively.

namespace decimate {

/// @brief Triangle face: three vertex indices (winding a→b→c CCW).
///
/// Matches the legacy `Face { int a, b, c; }` layout so loaders can hand a
/// `Face[]` straight to `HalfEdgeMesh::create()`.
struct Face {
    int a; ///< First  vertex index (0-based).
    int b; ///< Second vertex index.
    int c; ///< Third  vertex index.
};

} // namespace decimate
