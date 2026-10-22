// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file straighten_warp.h
/// @brief Internal: the harmonic interior warp and the foldover check (section 6,
///        "Warping the interior").  Eigen lives in the .cpp only.

#include "straighten_chains.h"

#include <decimate/math/vec2.h>
#include <vector>

namespace decimate::seam::detail {

/// Result of `warp`: the displaced UVs and the charts that must be reverted.
struct WarpOutput {
    std::vector<math::Vec2f> uvs;          ///< Final UVs (unmoved vertices keep their input bits).
    std::vector<int>         failedCharts; ///< Charts whose solve failed or whose triangles fold over.
};


/// @brief Apply the chain @p targets and displace every chart interior harmonically.
///
/// Per chart the interior displacement solves `L d = 0` with cotangent weights of the
/// ORIGINAL UVs and the silhouette displacements as Dirichlet data (Eigen
/// `SimplicialLDLT`, `SparseLU` when LDLT fails).  A chart fails when the solve fails, a UV
/// is not finite, or a triangle's signed area changes sign / becomes zero.
[[nodiscard]] WarpOutput warp(const Atlas &atlas, const ChainTargets &targets);

} // namespace decimate::seam::detail
