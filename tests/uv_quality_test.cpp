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
#include <cstdio>

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
#error "MESH_MODEL_DIR must be defined by CMake"
#endif

using namespace decimate;
using namespace decimate::chart;

namespace {

/// The cow, chartified and flattened: the shared input of the quality tests.
struct CowAtlas {
    MeshData                           data;
    ChartSet                           cs;
    std::vector<decimate::math::Vec2f> cornerUV; ///< per-corner (3*faceCount)
    int                                numFaces = 0;
};


/// Chart cap of the cow atlas: 25, the chart count the UVQ floors were calibrated on.  Before
/// Spec 19-01 a cap of 64 produced exactly 25 charts only through the seed-history bug (the cap
/// never bound); under the TrueCentroid default cap 25 yields 25 charts.
constexpr int kCowAtlasCap = 25;

/// Loads cow.obj (SKIPs if absent), segments it into at most @p cap charts and flattens with @p flattener.
CowAtlas loadCowAtlas(int cap = kCowAtlasCap, Flattener flattener = Flattener::Bijective) {
    CowAtlas c;
    c.data = loadObj(modelOrSkip("cow.obj"));
    c.numFaces = static_cast<int>(c.data.indices.size()) / 3;
    const int          np = static_cast<int>(c.data.positions.size());
    const HalfEdgeMesh hem = HalfEdgeMesh::create(np, c.numFaces,
                                                  reinterpret_cast<const Face *>(c.data.indices.data()));

    c.cs = segment(hem, c.data, cap);
    c.cornerUV = parametrizeCorners(hem, c.data, c.cs, 256, { .flattener = flattener });
    return c;
}


/// Result of rasterising all UV triangles into an R×R grid tagged by chart id: covered texels
/// and cross-chart double-covered texels.
struct Coverage {
    long covered = 0;    ///< texels covered by >=1 chart
    long crossChart = 0; ///< texels covered by >=2 DIFFERENT charts
    int  resolution = 0;
};


/// Rasterises every face of @p cow into a @p res × @p res coverage grid (centre-sampled, either winding).
Coverage rasterCoverage(const CowAtlas &cow, int res) {
    // cell[] holds the first chart that covered a texel; a second, different chart counts as overlap.
    std::vector<int> cell(static_cast<std::size_t>(res) * res, -1);
    Coverage         cov;
    cov.resolution = res;
    auto edge = [](float px, float py, float qx, float qy, float rx, float ry) {
        return (qx - px) * (ry - py) - (qy - py) * (rx - px);
    };
    for (int f = 0; f < cow.numFaces; ++f) {
        const int   chart = cow.cs.faceToChart[static_cast<std::size_t>(f)];
        const auto &a = cow.cornerUV[static_cast<std::size_t>(3 * f)];
        const auto &b = cow.cornerUV[static_cast<std::size_t>(3 * f + 1)];
        const auto &c = cow.cornerUV[static_cast<std::size_t>(3 * f + 2)];
        const float ax = a.x * static_cast<float>(res), ay = a.y * static_cast<float>(res);
        const float bx = b.x * static_cast<float>(res), by = b.y * static_cast<float>(res);
        const float cx = c.x * static_cast<float>(res), cy = c.y * static_cast<float>(res);
        const int   xlo = std::max(0, static_cast<int>(std::floor(std::min({ ax, bx, cx }))));
        const int   xhi = std::min(res - 1, static_cast<int>(std::ceil(std::max({ ax, bx, cx }))));
        const int   ylo = std::max(0, static_cast<int>(std::floor(std::min({ ay, by, cy }))));
        const int   yhi = std::min(res - 1, static_cast<int>(std::ceil(std::max({ ay, by, cy }))));
        for (int y = ylo; y <= yhi; ++y) {
            for (int x = xlo; x <= xhi; ++x) {
                const float px = static_cast<float>(x) + 0.5f;
                const float py = static_cast<float>(y) + 0.5f;
                const float w0 = edge(bx, by, cx, cy, px, py);
                const float w1 = edge(cx, cy, ax, ay, px, py);
                const float w2 = edge(ax, ay, bx, by, px, py);
                const bool  inside = (w0 >= 0 && w1 >= 0 && w2 >= 0) ||
                                     (w0 <= 0 && w1 <= 0 && w2 <= 0);
                if (!inside) {
                    continue;
                }
                int &tag = cell[static_cast<std::size_t>(y) * res + x];
                if (tag == -1) {
                    tag = chart;
                    ++cov.covered;
                } else if (tag != chart) {
                    ++cov.crossChart;
                }
            }
        }
    }
    return cov;
}

} // namespace


/// UVQ1: per-corner UVs — no triangle spans the atlas (seam-correct).
/// The per-vertex form smeared seam triangles across the whole [0,1]² (the
/// on-screen spaghetti); per-corner UVs cannot, by construction.
TEST_CASE("UVQ1: parametrizeCorners — no triangle spans the atlas", "[uv][quality]") {
    const CowAtlas cow = loadCowAtlas();
    REQUIRE(static_cast<int>(cow.cornerUV.size()) == cow.numFaces * 3);

    int spanning = 0;
    int oob = 0;
    int nan = 0;
    for (int f = 0; f < cow.numFaces; ++f) {
        const auto &a = cow.cornerUV[static_cast<std::size_t>(3 * f)];
        const auto &b = cow.cornerUV[static_cast<std::size_t>(3 * f + 1)];
        const auto &c = cow.cornerUV[static_cast<std::size_t>(3 * f + 2)];
        for (const auto *p : { &a, &b, &c }) {
            if (std::isnan(p->x) || std::isnan(p->y)) {
                ++nan;
            }
            if (p->x < -1e-3f || p->x > 1.0f + 1e-3f ||
                p->y < -1e-3f || p->y > 1.0f + 1e-3f) {
                ++oob;
            }
        }
        const float du = std::max({ a.x, b.x, c.x }) - std::min({ a.x, b.x, c.x });
        const float dv = std::max({ a.y, b.y, c.y }) - std::min({ a.y, b.y, c.y });
        if (du > 0.5f || dv > 0.5f) {
            ++spanning;
        }
    }

    CHECK(nan == 0);
    CHECK(oob == 0);
    CHECK(spanning == 0);
}


/// UVQ2: polygon (pixel) overlap — charts must not overlap in the atlas.
/// This is the CORRECT packing test (the old AABB test was wrong for a Tetris
/// packer).  Cross-chart double-covered texels should be a tiny fraction; the
/// remaining few are sub-pixel seam slivers, not gross overlap.
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


/// UVQ3: packing density — the Tetris packer fills a reasonable atlas fraction.
/// (The 2003 thesis explicitly notes density is not yet optimal — a lateral
/// horizon is the documented Beat-2b improvement.  This pins the current floor
/// so a future change can only raise it.)
TEST_CASE("UVQ3: silhouette packer — packing density above floor",
          "[uv][quality][packing]") {
    const CowAtlas cow = loadCowAtlas();
    const Coverage cov = rasterCoverage(cow, 256);
    const double   density =
        static_cast<double>(cov.covered) /
        static_cast<double>(cov.resolution * cov.resolution);

    INFO("charts = " << cow.cs.charts.size() << ", packing density = " << (density * 100.0) << "%");
    // Floor: the silhouette + width-fair Tetris packer achieves ~52% on the cow
    // (the per-column shadow lets charts nest into the skyline's concavities, and
    // the width-normalised rotation cost removes the vertical bias).  Assert a
    // conservative floor so a regression toward a box / loose packing fails — a
    // solid-rectangle packer scores far lower and would trip this.
    //
    // Floor is 0.30, not 0.40: the 4px miter gutter plus the atlas border cost ~15 density
    // points (the 0.40 floor was calibrated for gutterless packing).  The previous 0.35 came
    // from ONE partition (measured 0.365), produced while the seed-history bug pinned the cow
    // at 25 charts.  Density is noisy in the partition, so the floor is the minimum over caps
    // 20..30 minus a 0.02 margin, evaluated at the fixed cap 25 (kCowAtlasCap).
    //
    // Measured at HEAD of the Spec-19 fix round, build-rel (Release), cow.obj, atlas resolution 256,
    // TrueCentroid; reproduce with `mesh_tests "[.uvq3-table]"` (prints the table below).
    // Density = covered texels / 256^2, per chart cap:
    //   default (Flattener::Bijective):  cap 20..30 minimum 0.3227 at cap 27, cap 25 = 0.3369
    //   Flattener::Lscm (LSCM + repair): cap 20..30 minimum 0.3413 at cap 26, cap 25 = 0.3626
    // default per cap 20..30: .3267 .3733 .3420 .3611 .3589 .3369 .3367 .3227 .3360 .3284 .3496
    // lscm    per cap 20..30: .4133 .3838 .3784 .3732 .3592 .3626 .3413 .3757 .3471 .3445 .3507
    // The default flattener costs ~7% density at cap 25 (0.363 -> 0.337): the bijective free
    // boundaries are less packable than LSCM's.  The tested cap 25 (0.3369) sits above the
    // minimum over caps (0.3227), and the floor 0.30 is that minimum minus a 0.02 margin.  A box /
    // loose packer still scores far lower, so the regression guard is kept.
    CHECK(density > 0.30);
}


/// UVQ4: fold-overs — LSCM should not flip triangles wholesale.
TEST_CASE("UVQ4: parametrizeCorners — fold-overs are a small fraction",
          "[uv][quality]") {
    const CowAtlas cow = loadCowAtlas();
    int            positive = 0;
    int            negative = 0;
    for (int f = 0; f < cow.numFaces; ++f) {
        const auto &a = cow.cornerUV[static_cast<std::size_t>(3 * f)];
        const auto &b = cow.cornerUV[static_cast<std::size_t>(3 * f + 1)];
        const auto &c = cow.cornerUV[static_cast<std::size_t>(3 * f + 2)];
        const float ar = 0.5f * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
        if (ar < 0.0f) {
            ++negative;
        } else if (ar > 0.0f) {
            ++positive;
        }
    }

    const int flips = std::min(positive, negative);
    INFO("positive=" << positive << " negative=" << negative << " of " << cow.numFaces);
    CHECK(flips <= cow.numFaces / 20);
}


/// Hidden measurement (not a check): the UVQ3 density table over the caps 20..30, both flattener modes.
TEST_CASE("UVQ3 table: density per cap, default vs Lscm", "[.uvq3-table]") {
    for (const Flattener fl : { Flattener::Bijective, Flattener::Lscm }) {
        double minD = 1.0;
        int    minCap = 0;
        for (int cap = 20; cap <= 30; ++cap) {
            const CowAtlas cow = loadCowAtlas(cap, fl);
            const Coverage cov = rasterCoverage(cow, 256);
            const double   d = static_cast<double>(cov.covered) / static_cast<double>(cov.resolution * cov.resolution);
            std::printf("%s cap %d charts %zu density %.4f\n", fl == Flattener::Bijective ? "default" : "lscm   ", cap,
                        cow.cs.charts.size(), d);
            if (d < minD) {
                minD = d;
                minCap = cap;
            }
        }
        std::printf("MIN %.4f at cap %d\n", minD, minCap);
    }
}
