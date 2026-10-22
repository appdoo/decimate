// This file is part of the decimate library; see decimate.h for version/license details.
/// @file straighten.cpp
/// @brief seam::straighten -- protect, coarsen, warp, revert until stable (section 6).

#include "straighten_atlas.h"
#include "straighten_warp.h"

#include <algorithm>
#include <cstddef>
#include <decimate/seam/straighten.h>
#include <utility>

namespace decimate::seam {

namespace {

using detail::Atlas;

} // namespace


/// Counters of the final pass: per-rule protection counts, ears, chains.
static StraightenStats makeStats(const Atlas &atlas, const std::vector<std::uint8_t> &mask,
                                 const detail::ChainTargets &targets) {
    StraightenStats stats;
    for (const std::uint8_t m : mask) {
        for (std::size_t r = 0; r < 4; ++r) {
            stats.protectedByRule[r] += (m >> r) & 1;
        }
        stats.protectedVertices += m != 0 ? 1 : 0;
    }
    stats.earFaces = static_cast<int>(atlas.earFaces.size());
    stats.sameNeighbourEars = atlas.sameNeighbourEars;
    stats.seamChains = targets.seamChains;
    stats.boundaryChains = targets.boundaryChains;
    return stats;
}


/// Protects every UV vertex of chart @p chart; true if that protected anything new.
static bool forceChart(const Atlas &atlas, int chart, std::vector<std::uint8_t> &forced) {
    const TexturedMesh &tm = *atlas.tm;
    bool                fresh = false;
    for (const int f : atlas.chartFaces[static_cast<std::size_t>(chart)]) {
        for (std::size_t k = 0; k < 3; ++k) {
            auto &flag =
                forced[static_cast<std::size_t>(tm.faceU[(3 * static_cast<std::size_t>(f)) + k])];
            fresh = fresh || flag == 0;
            flag = 1;
        }
    }
    return fresh;
}


StraightenResult straighten(const TexturedMesh &tm, const StraightenOptions &opt) {
    const Atlas               atlas = detail::analyse(tm);
    std::vector<std::uint8_t> forced(tm.uvs.size(), 0);
    int                       reverted = 0;

    // Each round either finishes or protects at least one more vertex, so the loop terminates.
    for (;;) {
        const std::vector<std::uint8_t> mask = detail::protectionMasks(atlas, forced);
        std::vector<std::uint8_t>       isProtected(mask.size(), 0);
        std::ranges::transform(mask, isProtected.begin(),
                               [](std::uint8_t m) { return m != 0 ? 1 : 0; });

        // Straighten the chains, then let the chart interiors follow harmonically.
        const detail::ChainTargets targets = detail::straightenChains(atlas, isProtected, opt.tau);
        detail::WarpOutput         out = detail::warp(atlas, targets);

        // A chart whose warp folded over is reverted by protecting all of its vertices.
        int newlyReverted = 0;
        for (const int c : out.failedCharts) {
            newlyReverted += forceChart(atlas, c, forced) ? 1 : 0;
        }

        // Nothing left to protect: a failing chart that is already fully fixed cannot move.
        if (newlyReverted == 0) {
            StraightenResult result;
            result.mesh = tm;
            result.mesh.uvs = std::move(out.uvs);
            result.stats = makeStats(atlas, mask, targets);
            result.stats.revertedCharts = reverted;
            result.protectedUV = std::move(isProtected);
            return result;
        }
        reverted += newlyReverted;
    }
}

} // namespace decimate::seam
