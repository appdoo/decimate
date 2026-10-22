// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file halfedge_reduce.h
/// @brief Half-edge edge-collapse decimator — the data-oriented fast path.
///
/// Modernization of the author's 2003 `CHalfedgeCollapse`
/// (`legacy2003/reduce/chedgecollapse.cpp`). Connectivity and the guarded
/// collapse operator live in the shared `corner_mesh` module (the compact 8-byte
/// half-edge); this reducer is the **error-sort scheduler** over it: order edges
/// by initial QEM cost (lowest first) and collapse in active-set passes,
/// recomputing each optimum at visit. `reduce()` is the sibling that drives the
/// same `corner_mesh` collapse from a current-min heap instead. A collapse is
/// O(valence) re-linking of `sym`/`end`; faces are implicit (tombstoned,
/// reconstructed once at the end) — no per-vertex containers, no swap-back.
///
/// ## Representation
/// Each half-edge stores only `{ sym, end }` (8 bytes) — `next` is arithmetic
/// from the index (`+1`/`+1`/`-2` within the contiguous 3-per-face block) and
/// `face` is `h/3`. This is the **Directed-Edges** representation (Campagna,
/// Kobbelt & Seidel, 1998) / equivalently the **Corner Table** (Rossignac et
/// al., ~2001): the `end[]` array is the corner-vertex (V) table, `sym[]` is the
/// opposite (O) table. Half the bytes of a pointer/16-byte half-edge → ~2× the
/// half-edges per cache line (measured ~+12% at scale), and all `u32` fields are
/// 4-byte-aligned → Wasm-clean. Built via O(F) hash pairing.
///
/// The collapse + boundary-aware ring-walk (in `corner_mesh`) are the author's
/// 2003 algorithm (`chedgecollapse.cpp:122-178`), ported from raw `Halfedge*`
/// pointers to integer handles (portable, Wasm-safe, no pointer UB), with a
/// link condition + geometric quality guard (`corner::collapseLegal`) added for
/// 2-manifold safety. Quadrics + optimal placement reuse `QuadricErrorSystem`.
///
/// ## Boundary contract
/// Zero OS/GL dependencies. No stdout/stderr. No global state. Wasm-safe.

#include <decimate/types/mesh_data.h>

namespace decimate {

/// @brief Simplify a mesh toward @p targetVertexCount via in-place half-edge
///        edge collapse (QEM cost, lowest-cost-first).
///
/// Each collapse is link-condition + flip/degenerate guarded (HE-navigation
/// port of EdgeCollapse's guard), so the output is 2-manifold-preserving: a
/// clean reduced triangle mesh, all indices valid. If @p targetVertexCount ≥ the
/// current vertex count the mesh is returned unchanged. The guard is
/// conservative — it may stop above an extreme target rather than make a
/// topology-breaking collapse (e.g. it plateaus well above a 1%-of-original
/// target on the dragon).
///
/// @param input             Source mesh (positions + indices).
/// @param targetVertexCount Desired vertex count.
/// @returns                 Simplified `MeshData`.
/// @throws std::invalid_argument if a face index is outside `[0, positions.size())`.
MeshData reduceHalfEdge(const MeshData &input, int targetVertexCount);

} // namespace decimate
