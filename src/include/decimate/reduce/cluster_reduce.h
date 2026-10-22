// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file cluster_reduce.h
/// @brief Octree vertex-clustering decimator — the heap-free reducer strategy.
///
/// `reduceClustering()` is the second simplification strategy alongside the QEM
/// edge-collapse `reduce()` (reducer.h).  Where `reduce()` is *greedy* — a
/// priority heap repeatedly collapsing the lowest-error edge — this strategy
/// has **no heap and no per-collapse deletion**: it partitions space with an
/// adaptive octree, collapses every vertex in a leaf cell to one representative,
/// and drops the triangles that degenerate.  ~O(n), embarrassingly parallel.
///
/// Lineage: Rossignac-Borrel vertex clustering (1993), the quadric-optimal
/// representative of DeCoro 2007 / meshoptimizer `simplifySloppy`.  It reuses
/// `QError` (qerror.h) verbatim for the per-cell quadric, so the *representative
/// placement* is error-optimal even though the *collapse order* is not.
///
/// Trade-off vs `reduce()`: far faster, lower quality (no error-optimal order).
///
/// ## Boundary contract
/// Zero OS/GL dependencies.  No stdout/stderr.  No global state.  Wasm-safe.

#include <decimate/types/mesh_data.h>

namespace decimate {

/// @brief Simplify a mesh toward @p targetVertexCount via octree vertex
///        clustering (heap-free).
///
/// The target sets the octree granularity (leaf size ≈ vertexCount / target),
/// so the output vertex count is **approximate** — clustering collapses by
/// spatial cell, not to an exact count.  If @p targetVertexCount ≥ the current
/// vertex count, the mesh is returned unchanged.
///
/// Each output vertex is the quadric-optimal point of its cell (centroid
/// fallback on a singular solve), clamped to the cell's bounding box so a
/// near-singular solve cannot fling a representative into a spike.  All output
/// indices are valid (< output positions.size()); no degenerate triangles.
///
/// @param input             Source mesh (positions + indices).
/// @param targetVertexCount Approximate desired vertex count.
/// @returns                 Simplified `MeshData`.
MeshData reduceClustering(const MeshData &input, int targetVertexCount);

} // namespace decimate
