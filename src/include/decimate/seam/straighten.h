// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file straighten.h
/// @brief Seam straightening (Liu et al. 2017, section 6) as an opt-in UV post-process.
///
/// Seam and boundary chains of a packed atlas are coarsened (Ramer-Douglas-Peucker), the
/// vertices of every seam are re-spaced so both copies agree on their edge-length
/// ratios, and the chart interiors follow by one harmonic displacement step.  Only UVs
/// move: the faces, the surface indices and the cutter's chart assignment are untouched.
///
/// Deliberate simplifications versus the paper: a single harmonic step instead of the
/// composite mapping, and no constrained-Delaunay ambient mesh.  A chart whose warp
/// folds over is reverted (its vertices become protected) instead.

#include <array>
#include <cstdint>
#include <decimate/types/textured_mesh.h>
#include <vector>

namespace decimate::seam {

/// @brief Tuning of `straighten`.
struct StraightenOptions {
    /// Coarsening tolerance in atlas UV units (`[0,1]^2`); the paper uses 0.01.  Seam
    /// chains measure it in the 4D concatenation of both copies, boundary chains in 2D.
    double tau = 0.01;
};


/// @brief Counters describing one `straighten` run (final pass, after any reverts).
struct StraightenStats {
    /// UV vertices protected by paper rule (1)..(4).  A vertex may satisfy several rules
    /// and is then counted under each; rule (4) counts only vertices that were protected
    /// by propagation alone (copies of a protected surface vertex, reverted charts).
    std::array<int, 4> protectedByRule {};
    int                protectedVertices = 0; ///< Distinct protected UV vertices (the union).
    int                earFaces = 0;          ///< Faces with exactly two seam/boundary half-edges (rule 3).
    /// Ears whose two half-edges are seams facing one common neighbour chart (the shape
    /// the cutter's `straighten` pass targets).
    int sameNeighbourEars = 0;
    int seamChains = 0;     ///< Seam chains (both copies count once) with an interior vertex.
    int boundaryChains = 0; ///< Boundary chains with an interior vertex.
    int revertedCharts = 0; ///< Charts rolled back because the warp folded them over.
};


/// @brief Output of `straighten`.
struct StraightenResult {
    TexturedMesh    mesh; ///< Same faces and UV indexing as the input, moved UVs.
    StraightenStats stats;
    /// One entry per UV vertex: non-zero when the vertex was protected and therefore
    /// bit-identical to the input.
    std::vector<std::uint8_t> protectedUV;
};


/// @brief Straighten the seam and boundary chains of @p tm (paper section 6).
///
/// ## Contract
/// - **Writes UVs only.** `faceV` and `faceU` of the result equal the input's; chart
///   assignment is the cutter's business.
/// - Protected vertices (stats / `protectedUV`) are bit-identical to the input.
/// - No UV triangle changes the sign of its signed area or becomes degenerate: a chart
///   violating that is reverted to its input UVs and counted in `revertedCharts`.
/// - On every chain run, adjacent edges whose 4D (both copies) directions are parallel
///   have equal edge-length ratios on the two copies.
/// - Seams and boundaries are classified by `findSeams` only.
///
/// @param tm   Welded atlas (`weld()` output, post-pack UVs).  Each UV vertex must belong
///             to exactly one surface vertex, which `weld()` guarantees.
/// @param opt  Tolerance; see `StraightenOptions`.
/// @throws std::invalid_argument if a UV vertex is shared by two surface vertices.
[[nodiscard]] StraightenResult straighten(const TexturedMesh      &tm,
                                          const StraightenOptions &opt = {});

} // namespace decimate::seam
