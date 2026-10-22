// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file straighten_chains.h
/// @brief Internal: edge-chain components and their Ramer-Douglas-Peucker coarsening
///        (Liu et al. 2017, section 6, "Edge component straightening").

#include "straighten_atlas.h"

#include <array>
#include <cstdint>
#include <vector>

namespace decimate::seam::detail {

/// A chain point: 2D UV (components 2..3 zero) or the 4D concatenation of a seam's two copies.
using Pt4 = std::array<double, 4>;


/// Target UV positions of the vertices the chain pass moves.
struct ChainTargets {
    std::vector<std::array<double, 2>> pos; ///< Per UV vertex; valid where `moved != 0`.
    std::vector<std::uint8_t>          moved;
    int                                seamChains = 0;     ///< Seam chains (a chain and its copy count once) with interior vertices.
    int                                boundaryChains = 0; ///< Boundary chains with interior vertices.
};


/// @brief Coarsen @p pts with Ramer-Douglas-Peucker; returns one keep flag per point.
///
/// Distances are Euclidean in all four components (2D chains pass zeros in 2..3).  The end
/// points are always kept; a @p closed loop (`pts.front() == pts.back()`) is split at its
/// farthest point and never drops below three kept vertices.
[[nodiscard]] std::vector<std::uint8_t> coarsen(const std::vector<Pt4> &pts, bool closed,
                                                double tau);


/// @brief Re-space the dropped points onto the coarsened polyline by normalized arc length.
///
/// Kept points keep their coordinates.  Between two kept points `a < b`, point `i` goes to
/// `lerp(pts[a], pts[b], (S_i - S_a) / (S_b - S_a))` with `S` the cumulative original
/// arc length (uniform fallback when the span has zero length).
[[nodiscard]] std::vector<Pt4> reparametrize(const std::vector<Pt4>          &pts,
                                             const std::vector<std::uint8_t> &keep);


/// @brief Straighten every seam and boundary chain between protected vertices.
/// @param isProtected  Per UV vertex, non-zero = fixed.
[[nodiscard]] ChainTargets
straightenChains(const Atlas &atlas, const std::vector<std::uint8_t> &isProtected, double tau);

} // namespace decimate::seam::detail
