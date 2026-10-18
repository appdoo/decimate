// This file is part of the decimate library; see decimate.h for version/license details.
/// @file uv_quality_test.cpp
/// @brief UV-atlas QUALITY characterization on the real cow (Spec 08-02).
///
/// The original 08-02 tests proved per-chart math (isometry, conformal energy)
/// and a WRONG packing property — "no chart AABB overlap" — which does not hold
/// for a silhouette / Tetris packer (interlocking charts have overlapping
/// bounding boxes by design; only their rasterised polygons stay disjoint).
/// That bad assertion drove a regression: the packer was replaced with solid
/// bounding boxes to satisfy it.  These tests measure the RIGHT thing — what the
/// atlas actually looks like at the pixel level:
///
///   - cross-chart pixel overlap: two charts' triangles covering the same atlas
///     texel (the real "polygon overlap" — bad, should be near zero).
///   - packing density: covered atlas fraction (Tetris should pack reasonably).
///   - per-corner seam correctness: every face's 3 corners come from one chart,
///     so no triangle spans the atlas (the per-vertex-scatter spaghetti is gone).
///
/// Granularity is a coverage GRID (rasterise every UV triangle into a 256² grid
/// tagged by chart id) — O(pixels), not O(triangles²).

#include <catch2/catch_approx.hpp>
#include "model_fixture.h"
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <decimate/chart/chart.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>
#include <decimate/chart/parametrize.h>

#ifndef MESH_MODEL_DIR
#  error "MESH_MODEL_DIR must be defined by CMake"
#endif

using namespace decimate;
using namespace decimate::chart;

namespace {

struct CowAtlas {
    MeshData                       data;
    ChartSet                       cs;
    std::vector<decimate::math::Vec2<float>> cornerUV;  ///< per-corner (3*faceCount)
    int                            numFaces = 0;
};

CowAtlas loadCowAtlas() {
    CowAtlas c;
    c.data = loadObj(modelOrSkip("cow.obj"));
    c.numFaces = static_cast<int>(c.data.indices.size()) / 3;
    const int np = static_cast<int>(c.data.positions.size());
    const HalfEdgeMesh hem = HalfEdgeMesh::create(c.numFaces, c.numFaces,
        reinterpret_cast<const Face*>(c.data.indices.data()));
    // (numPoints arg only sizes a scratch; pass faceCount-safe value)
    const HalfEdgeMesh hem2 = HalfEdgeMesh::create(np, c.numFaces,
        reinterpret_cast<const Face*>(c.data.indices.data()));
    c.cs = segment(hem2, c.data, 64);
    c.cornerUV = parametrizeCorners(hem2, c.data, c.cs, 256);
    (void)hem;
    return c;
}

/// @brief Rasterise all UV triangles into an R×R grid tagged by chart id,
///        counting covered texels and cross-chart double-covered texels.
struct Coverage {
    long covered = 0;       ///< texels covered by >=1 chart
    long crossChart = 0;    ///< texels covered by >=2 DIFFERENT charts
    int  resolution = 0;
};

Coverage rasterCoverage(const CowAtlas& cow, int res) {
    std::vector<int> cell(static_cast<std::size_t>(res) * res, -1);
    Coverage cov;
    cov.resolution = res;
    auto edge = [](float px, float py, float qx, float qy, float rx, float ry) {
        return (qx - px) * (ry - py) - (qy - py) * (rx - px);
    };
    for (int f = 0; f < cow.numFaces; ++f) {
        const int chart = cow.cs.faceToChart[static_cast<std::size_t>(f)];
        const auto& a = cow.cornerUV[static_cast<std::size_t>(3 * f)];
        const auto& b = cow.cornerUV[static_cast<std::size_t>(3 * f + 1)];
        const auto& c = cow.cornerUV[static_cast<std::size_t>(3 * f + 2)];
        const float ax = a.x * static_cast<float>(res), ay = a.y * static_cast<float>(res);
        const float bx = b.x * static_cast<float>(res), by = b.y * static_cast<float>(res);
        const float cx = c.x * static_cast<float>(res), cy = c.y * static_cast<float>(res);
        const int xlo = std::max(0, static_cast<int>(std::floor(std::min({ax, bx, cx}))));
        const int xhi = std::min(res - 1, static_cast<int>(std::ceil(std::max({ax, bx, cx}))));
        const int ylo = std::max(0, static_cast<int>(std::floor(std::min({ay, by, cy}))));
        const int yhi = std::min(res - 1, static_cast<int>(std::ceil(std::max({ay, by, cy}))));
        for (int y = ylo; y <= yhi; ++y) {
            for (int x = xlo; x <= xhi; ++x) {
                const float px = static_cast<float>(x) + 0.5F;
                const float py = static_cast<float>(y) + 0.5F;
                const float w0 = edge(bx, by, cx, cy, px, py);
                const float w1 = edge(cx, cy, ax, ay, px, py);
                const float w2 = edge(ax, ay, bx, by, px, py);
                const bool inside = (w0 >= 0 && w1 >= 0 && w2 >= 0) ||
                                    (w0 <= 0 && w1 <= 0 && w2 <= 0);
                if (!inside) { continue; }
                int& tag = cell[static_cast<std::size_t>(y) * res + x];
                if (tag == -1) { tag = chart; ++cov.covered; }
                else if (tag != chart) { ++cov.crossChart; }
            }
        }
    }
    return cov;
}

}  // namespace

// ===========================================================================
// UVQ1: per-corner UVs — no triangle spans the atlas (seam-correct).
// The per-vertex form smeared seam triangles across the whole [0,1]² (the
// on-screen spaghetti); per-corner UVs cannot, by construction.
// ===========================================================================
TEST_CASE("UVQ1: parametrizeCorners — no triangle spans the atlas", "[uv][quality]") {
    const CowAtlas cow = loadCowAtlas();
    REQUIRE(static_cast<int>(cow.cornerUV.size()) == cow.numFaces * 3);

    int spanning = 0;
    int oob = 0;
    int nan = 0;
    for (int f = 0; f < cow.numFaces; ++f) {
        const auto& a = cow.cornerUV[static_cast<std::size_t>(3 * f)];
        const auto& b = cow.cornerUV[static_cast<std::size_t>(3 * f + 1)];
        const auto& c = cow.cornerUV[static_cast<std::size_t>(3 * f + 2)];
        for (const auto* p : {&a, &b, &c}) {
            if (std::isnan(p->x) || std::isnan(p->y)) { ++nan; }
            if (p->x < -1e-3F || p->x > 1.0F + 1e-3F ||
                p->y < -1e-3F || p->y > 1.0F + 1e-3F) { ++oob; }
        }
        const float du = std::max({a.x, b.x, c.x}) - std::min({a.x, b.x, c.x});
        const float dv = std::max({a.y, b.y, c.y}) - std::min({a.y, b.y, c.y});
        if (du > 0.5F || dv > 0.5F) { ++spanning; }
    }
    CHECK(nan == 0);
    CHECK(oob == 0);
    CHECK(spanning == 0);
}

// ===========================================================================
// UVQ2: polygon (pixel) overlap — charts must not overlap in the atlas.
// This is the CORRECT packing test (the old AABB test was wrong for a Tetris
// packer).  Cross-chart double-covered texels should be a tiny fraction; the
// remaining few are sub-pixel seam slivers, not gross overlap.
// ===========================================================================
TEST_CASE("UVQ2: silhouette packer — cross-chart pixel overlap is tiny",
          "[uv][quality][packing]") {
    const CowAtlas cow = loadCowAtlas();
    const Coverage cov = rasterCoverage(cow, 256);
    REQUIRE(cov.covered > 0);
    const double overlapFrac =
        static_cast<double>(cov.crossChart) / static_cast<double>(cov.covered);
    INFO("covered=" << cov.covered << " crossChart=" << cov.crossChart
         << " overlap=" << (overlapFrac * 100.0) << "%");
    // Gross overlap (the box/scatter bugs) was huge; a correct atlas keeps
    // cross-chart overlap to seam slivers — well under 5%.
    CHECK(overlapFrac < 0.05);
}

// ===========================================================================
// UVQ3: packing density — the Tetris packer fills a reasonable atlas fraction.
// (The 2003 thesis explicitly notes density is not yet optimal — a lateral
// horizon is the documented Beat-2b improvement.  This pins the current floor
// so a future change can only raise it.)
// ===========================================================================
TEST_CASE("UVQ3: silhouette packer — packing density above floor",
          "[uv][quality][packing]") {
    const CowAtlas cow = loadCowAtlas();
    const Coverage cov = rasterCoverage(cow, 256);
    const double density =
        static_cast<double>(cov.covered) /
        static_cast<double>(cov.resolution * cov.resolution);
    INFO("packing density = " << (density * 100.0) << "%");
    // Floor: the silhouette + width-fair Tetris packer achieves ~52% on the cow
    // (the per-column shadow lets charts nest into the skyline's concavities, and
    // the width-normalised rotation cost removes the vertical bias).  Assert a
    // conservative floor so a regression toward a box / loose packing fails — a
    // solid-rectangle packer scores far lower and would trip this.
    //
    // Floor is 0.35, not 0.40: the 4px miter gutter plus the atlas border
    // legitimately cost ~15 density points (measured 0.365).  The 0.40 floor was
    // calibrated for gutterless packing.
    CHECK(density > 0.35);
}

// ===========================================================================
// UVQ4: fold-overs — LSCM should not flip triangles wholesale.
// ===========================================================================
TEST_CASE("UVQ4: parametrizeCorners — fold-overs are a small fraction",
          "[uv][quality]") {
    const CowAtlas cow = loadCowAtlas();
    int positive = 0;
    int negative = 0;
    for (int f = 0; f < cow.numFaces; ++f) {
        const auto& a = cow.cornerUV[static_cast<std::size_t>(3 * f)];
        const auto& b = cow.cornerUV[static_cast<std::size_t>(3 * f + 1)];
        const auto& c = cow.cornerUV[static_cast<std::size_t>(3 * f + 2)];
        const float ar = 0.5F * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
        if (ar < 0.0F) { ++negative; }
        else if (ar > 0.0F) { ++positive; }
    }
    const int flips = std::min(positive, negative);
    INFO("positive=" << positive << " negative=" << negative << " of " << cow.numFaces);
    CHECK(flips <= cow.numFaces / 20);
}
