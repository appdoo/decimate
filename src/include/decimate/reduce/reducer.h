// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file reducer.h
/// @brief Public entry-point for QEM mesh simplification (Beat-2a modernisation).
///
/// `reduce()` is the primary API consumed by tests (Spec 07), the future
/// viewer (Spec 10), and the bake pipeline.  It takes a `MeshData` (from the
/// loader) and returns a simplified `MeshData`.
///
/// ## Algorithm
///
/// Garland-Heckbert 1997 quadric error metric edge-collapse.
///
/// ## Boundary contract
///
/// Zero OS/GL dependencies.  No stdout/stderr.  No global state.
/// Suitable for both native and Emscripten (Wasm) builds.

#include <decimate/types/mesh_data.h>

namespace decimate {

/// @brief Simplify a mesh to at most @p targetVertexCount vertices using QEM
///        edge-collapse (Garland-Heckbert 1997).
///
/// If @p targetVertexCount ≥ the current vertex count, the mesh is returned
/// unchanged.  The algorithm may stop before reaching the target if the heap
/// empties (all remaining collapses would produce degenerate geometry).
///
/// @param input             Source mesh (positions + indices).
/// @param targetVertexCount Maximum number of vertices in the output.
/// @returns                 Simplified `MeshData` with packed positions and
///                          remapped face indices.  All indices are valid
///                          (< positions.size()).
MeshData reduce(const MeshData& input, int targetVertexCount);

}  // namespace decimate
