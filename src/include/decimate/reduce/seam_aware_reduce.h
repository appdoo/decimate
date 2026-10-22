// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file seam_aware_reduce.h
/// @brief Seam-aware edge-collapse decimation of a textured mesh (Liu et al. 2017, section 5).
///
/// `reduce()` simplifies positions only; handed a mesh with a texture atlas it would cut
/// corners of the UV charts and reintroduce exactly the seam discontinuity that
/// `seam::erase` removed.  `reduceSeamAware()` is the textured member of the reducer
/// family: it decimates a `TexturedMesh` so that the ONE erased texture can be sampled on
/// every level of detail without re-baking and without a new seam.
///
/// ## Algorithm
///
/// - **Two corner tables, one guard.**  The surface (`faceV`) and the UV mesh (`faceU`)
///   each get a `corner::build` table; half-edge `h` names the same corner in both.  Every
///   collapse must pass `corner::collapseLegal` on the surface table and on the UV table of
///   each UV half-edge it touches (the UV pass is the UV-foldover guard).  No link or flip
///   test is repeated here.
/// - **Collapse criteria (paper section 5.1).**  An interior edge joining two UV-silhouette
///   vertices is not collapsible; with exactly one silhouette endpoint the new UV is that
///   vertex's UV.  A seam edge collapses only when at least one adjacent seam edge is
///   *unifiable* on both copies (collinear within |sin| <= 1e-5 and equal edge-length
///   ratio within 1e-5 relative).  A UV-silhouette edge that is not a seam (a 3D boundary)
///   needs collinearity only.
/// - **Metric and placement.**  Garland-Heckbert 1998 5D quadrics (x, y, z, u, v) per
///   wedge, summed on collapse.  The placement is the exact minimizer of the quadric under
///   the criteria's equalities (closed form, the single inequality `0 <= t <= 1` resolved
///   by clamping) -- no general QP solver.
/// - **Schedule.**  Lowest cost first from a lazy priority queue with stale-entry skip, as
///   `reduce()`.
///
/// ## Boundary contract
///
/// Zero OS/GL dependencies, no global state, no I/O.  Does not link `decimate::seam`:
/// seams are read off the UV corner table (a surface edge with a twin whose UV half-edges
/// have none), which is the same UV-defined criterion `seam::findSeams` uses.  Eigen is an
/// implementation detail of the `.cpp`.

#include <decimate/types/textured_mesh.h>

namespace decimate {

/// @brief Counters of one `reduceSeamAware` run.
struct SeamAwareStats {
    /// Seam edges (surface edges whose two sides use different UV vertices) of the input.
    int seamEdgesBefore = 0;
    /// Seam edges of the output.
    int seamEdgesAfter = 0;
    /// Input seam edges with at least one unifiable end (`classifySeamEdges` of the input).
    int collapsibleSeamEdges = 0;
    /// Input seam edges with no unifiable end: never collapsible.
    int uncollapsibleSeamEdges = 0;
    /// Seam edges actually collapsed.
    int collapsedSeamEdges = 0;
};


/// @brief Decimated mesh plus counters.
struct SeamAwareResult {
    TexturedMesh   mesh; ///< Same invariants as the input (`TexturedMesh` contract).
    SeamAwareStats stats;
};


/// @brief Static seam-edge classification (no quadrics, no legality test).
struct SeamCollapsibility {
    int collapsible = 0;   ///< Seam edges with at least one unifiable end.
    int uncollapsible = 0; ///< Seam edges with none.
};


/// @brief Decimate @p in to at most @p targetVertexCount surface vertices, seam-aware.
///
/// ## Contract
/// - `targetVertexCount` counts surface vertices (as `reduce`).  If it is >= the input
///   count the input is returned bit-for-bit with every counter 0 except
///   `seamEdgesBefore` / `seamEdgesAfter`.
/// - Like the reducer family it may stop above the target when no legal collapse remains.
/// - The output keeps the `TexturedMesh` invariants; every surface edge keeps 1 or 2
///   faces, every vertex link stays one fan, no UV triangle flips.
/// - Seams only ever shorten along the straight line they already follow: a UV position on
///   a silhouette stays on its input silhouette, and the two copies of a collapsed seam
///   keep a common parametrization, so the erased texture still matches across the seam.
/// - Deterministic: equal inputs give bit-equal outputs.
///
/// @param in                 Welded atlas (`seam::weld()` output).  Each UV vertex must
///                           belong to exactly one surface vertex.
/// @param targetVertexCount  Desired surface vertex count.
/// @throws std::invalid_argument if a UV vertex is shared by two surface vertices, or an
///         index is out of range.
[[nodiscard]] SeamAwareResult reduceSeamAware(const TexturedMesh &in, int targetVertexCount);


/// @brief Count the seam edges of @p tm that have at least one unifiable end vs none.
///
/// The predicate `reduceSeamAware` applies before any quadric: the adjacent seam edge at
/// an end must exist (the vertex has exactly two seam edges), be collinear with the edge
/// on both copies and have the same edge-length ratio on both.
///
/// @throws std::invalid_argument as `reduceSeamAware`.
[[nodiscard]] SeamCollapsibility classifySeamEdges(const TexturedMesh &tm);

} // namespace decimate
