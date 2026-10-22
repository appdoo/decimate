// This file is part of the decimate library; see decimate.h for version/license details.
/// @file atlas_invariant_test.cpp
/// @brief Spec 19-02 (Amendment 2): the packed atlas must not put two charts on the same texel.
///
/// For cow / spot / bunny / horse at caps 16 / 32 / 64 / 128 the full pipeline runs
/// (`segment` -> `flattenCharts` with the injectivity guard -> `packAtlas`), and every placed UV
/// triangle is rasterised into a 1024 x 1024 texel grid (texel centres, UV space):
///
///  - `intInInt`: texels covered by the UV triangles of two different charts.  MUST be 0 - this is
///    chart-on-chart overlap, the atlas bug that made bunny's baked texture wrong.  Before the
///    self-seam cut (Spec 19-02) bunny measured 58403 / 7161 / 2121 / 82 here, caused by non-disk
///    charts whose packer silhouette (first boundary loop only) missed a hole.
///  - `intInGut`: texels of one chart's triangles that fall inside ANOTHER chart's gutter outline.
///    Bounded above by the values measured at 648e6d6 (cow/128 3, spot/128 1, all others 0); it
///    may only go down.
///  - `deepest`: how far (in 256-grid pixels) a gutter vertex of one chart reaches into another
///    chart's gutter loop.  Reported, not asserted.
///
/// The table is printed on every run.  The run is heavy in the unoptimised build (labelled
/// `slow`); it is a ctest in both build directories.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>
#include <decimate/math/vec2.h>

#include "model_fixture.h"

#ifndef MESH_MODEL_DIR
#error "MESH_MODEL_DIR must be defined by CMake"
#endif

using namespace decimate;
using namespace decimate::chart;
using Vec2f = decimate::math::Vec2f;

namespace {

constexpr int kGrid = 1024; ///< Texel grid resolution in UV space.
/// A gutter loop that touches itself within this distance is simple (float rounding of the placed loop, 0.0026 of a 255-grid pixel).
constexpr float kTouchTolUv = 1e-5f;

/// Even-odd point-in-polygon test (used to find gutter vertices inside another chart's gutter).
bool insideLoop(const Vec2f &p, const std::vector<Vec2f> &poly) {
    bool in = false;
    for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        const bool straddles = (poly[i].y > p.y) != (poly[j].y > p.y);
        if (straddles && p.x < (poly[j].x - poly[i].x) * (p.y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x) {
            in = !in;
        }
    }
    return in;
}


/// Distance from @p p to the polygon outline, in 256-grid pixels.
float distToLoopPx(const Vec2f &p, const std::vector<Vec2f> &poly) {
    float best = 1e9f;
    for (std::size_t t = 0; t < poly.size(); ++t) {
        const Vec2f &a = poly[t];
        const Vec2f &b = poly[(t + 1) % poly.size()];
        const float  ex = b.x - a.x;
        const float  ey = b.y - a.y;
        const float  u = std::clamp(((p.x - a.x) * ex + (p.y - a.y) * ey) / (ex * ex + ey * ey), 0.0f, 1.0f);
        best = std::min(best, std::hypot(p.x - a.x - u * ex, p.y - a.y - u * ey));
    }
    return best * 255.0f;
}


/// Calls f(col, row) for every texel centre inside @p poly (even-odd scanline fill).
template <class F>
// NOLINTNEXTLINE(readability-function-size): scanline kernel, one pass over rows and crossings
void fillLoop(const std::vector<Vec2f> &poly, F f) {
    float y0 = 1e9f;
    float y1 = -1e9f;
    for (const Vec2f &p : poly) {
        y0 = std::min(y0, p.y);
        y1 = std::max(y1, p.y);
    }
    const int          r0 = std::max(0, static_cast<int>(std::floor(y0 * kGrid - 0.5f)));
    const int          r1 = std::min(kGrid - 1, static_cast<int>(std::ceil(y1 * kGrid - 0.5f)));
    std::vector<float> xs;
    for (int r = r0; r <= r1; ++r) {
        const float y = (static_cast<float>(r) + 0.5f) / kGrid;
        xs.clear();
        for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
            if ((poly[i].y > y) != (poly[j].y > y)) {
                xs.push_back((poly[j].x - poly[i].x) * (y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x);
            }
        }
        std::sort(xs.begin(), xs.end());
        for (std::size_t k = 0; k + 1 < xs.size(); k += 2) {
            const int c0 = std::max(0, static_cast<int>(std::ceil(xs[k] * kGrid - 0.5f)));
            const int c1 = std::min(kGrid - 1, static_cast<int>(std::floor(xs[k + 1] * kGrid - 0.5f)));
            for (int c = c0; c <= c1; ++c) {
                f(c, r);
            }
        }
    }
}


/// Texel ownership of the packed UV triangles: `owner` = first chart to cover a texel,
/// `multi` = a second, different chart covers it too.
struct Coverage {
    std::vector<int>  owner = std::vector<int>(static_cast<std::size_t>(kGrid) * kGrid, -1);
    std::vector<char> multi = std::vector<char>(static_cast<std::size_t>(kGrid) * kGrid, 0);
    long              covered = 0;
};


/// Rasterise one UV triangle (texel centres) into @p cov for chart @p chart.
// NOLINTNEXTLINE(readability-function-size): barycentric raster kernel
void rasterTriangle(const std::array<Vec2f, 3> &p, int chart, Coverage &cov) {
    const float ax = p[0].x * kGrid;
    const float ay = p[0].y * kGrid;
    const float bx = p[1].x * kGrid;
    const float by = p[1].y * kGrid;
    const float cx = p[2].x * kGrid;
    const float cy = p[2].y * kGrid;
    const float det = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy);
    if (std::fabs(det) < 1e-12f) {
        return;
    }
    const int x0 = std::max(0, static_cast<int>(std::floor(std::min({ ax, bx, cx }))));
    const int x1 = std::min(kGrid - 1, static_cast<int>(std::ceil(std::max({ ax, bx, cx }))));
    const int y0 = std::max(0, static_cast<int>(std::floor(std::min({ ay, by, cy }))));
    const int y1 = std::min(kGrid - 1, static_cast<int>(std::ceil(std::max({ ay, by, cy }))));
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            const float px = static_cast<float>(x) + 0.5f;
            const float py = static_cast<float>(y) + 0.5f;
            const float l1 = ((by - cy) * (px - cx) + (cx - bx) * (py - cy)) / det;
            const float l2 = ((cy - ay) * (px - cx) + (ax - cx) * (py - cy)) / det;
            if (l1 < 0.0f || l2 < 0.0f || l1 + l2 > 1.0f) {
                continue;
            }
            int &owner = cov.owner[static_cast<std::size_t>(y) * kGrid + static_cast<std::size_t>(x)];
            if (owner < 0) {
                owner = chart;
                ++cov.covered;
            } else if (owner != chart) {
                cov.multi[static_cast<std::size_t>(y) * kGrid + static_cast<std::size_t>(x)] = 1;
            }
        }
    }
}


/// Counters measured on one packed atlas.
struct Measure {
    std::size_t charts = 0;
    long        intInInt = 0;
    long        intInGut = 0;
    long        selfCrossingGutters = 0; ///< charts whose gutterLoop has a proper self-crossing (Spec 21)
    float       deepest = 0.0f;
    double      density = 0.0;
};


// NOLINTNEXTLINE(readability-function-size)
Coverage coverageOf(const std::vector<HalfEdgeChart> &atlas) {
    Coverage cov;
    for (std::size_t c = 0; c < atlas.size(); ++c) {
        const HalfEdgeChart &ac = atlas[c];
        for (int i = 0; i < static_cast<int>(ac.faces.size()); ++i) {
            std::array<Vec2f, 3> p;
            for (int k = 0; k < 3; ++k) {
                p[static_cast<std::size_t>(k)] = ac.uv[static_cast<std::size_t>(localCorner(ac, i, k))];
            }
            rasterTriangle(p, static_cast<int>(c), cov);
        }
    }
    return cov;
}


/// Texels of chart A's triangles that lie inside chart B's gutter outline (A != B).
// NOLINTNEXTLINE(readability-function-size)
long trianglesInGutters(const std::vector<HalfEdgeChart> &atlas, const Coverage &cov) {
    std::vector<char> hit(cov.owner.size(), 0);
    for (std::size_t b = 0; b < atlas.size(); ++b) {
        // A chart without a gutter outline (no loop or a sliver) has nothing to fill.
        if (atlas[b].gutterLoop.size() < 3) {
            continue;
        }
        fillLoop(atlas[b].gutterLoop, [&](int x, int y) {
            const std::size_t t = static_cast<std::size_t>(y) * kGrid + static_cast<std::size_t>(x);
            if (cov.owner[t] >= 0 && cov.owner[t] != static_cast<int>(b)) {
                hit[t] = 1;
            }
        });
    }
    return std::count(hit.begin(), hit.end(), static_cast<char>(1));
}


/// Deepest penetration of one chart's gutter vertex into another chart's gutter loop.
// NOLINTNEXTLINE(readability-function-size)
float deepestGutterPenetration(const std::vector<HalfEdgeChart> &atlas) {
    struct Box {
        float x0 = 1e30f, y0 = 1e30f, x1 = -1e30f, y1 = -1e30f;
    };
    std::vector<Box> box(atlas.size());

    // Per-chart bounding boxes of the gutter loops, to skip chart pairs that cannot touch.
    for (std::size_t c = 0; c < atlas.size(); ++c) {
        for (const Vec2f &g : atlas[c].gutterLoop) {
            box[c] = { std::min(box[c].x0, g.x), std::min(box[c].y0, g.y), std::max(box[c].x1, g.x),
                       std::max(box[c].y1, g.y) };
        }
    }

    float deepest = 0.0f;
    for (std::size_t i = 0; i < atlas.size(); ++i) {
        for (std::size_t j = 0; j < atlas.size(); ++j) {
            const bool apart = box[i].x1 < box[j].x0 || box[j].x1 < box[i].x0 || box[i].y1 < box[j].y0 ||
                               box[j].y1 < box[i].y0;
            if (i == j || apart) {
                continue;
            }
            for (const Vec2f &v : atlas[j].gutterLoop) {
                if (insideLoop(v, atlas[i].gutterLoop)) {
                    deepest = std::max(deepest, distToLoopPx(v, atlas[i].gutterLoop));
                }
            }
        }
    }
    return deepest;
}


/// Runs the full pipeline on one model at one chart cap and measures the packed atlas.
Measure measureAtlas(const std::string &model, int cap) {
    // Same pipeline as production: segment, flatten with the injectivity guard, pack at 256.
    const MeshData     data = load(modelOrSkip(model));
    const HalfEdgeMesh hem = HalfEdgeMesh::create(
        static_cast<int>(data.positions.size()), static_cast<int>(data.indices.size() / 3),
        reinterpret_cast<const Face *>(data.indices.data()));
    ChartSet                   cs = segment(hem, data, cap);
    std::vector<HalfEdgeChart> atlas = flattenCharts(hem, data, cs);
    packAtlas(atlas, data, 256);

    // The three counters below are what the test asserts or reports.
    const Coverage cov = coverageOf(atlas);
    Measure        m;
    m.charts = atlas.size();
    m.intInInt = std::count(cov.multi.begin(), cov.multi.end(), static_cast<char>(1));
    m.intInGut = trianglesInGutters(atlas, cov);
    m.deepest = deepestGutterPenetration(atlas);

    // A self-crossing gutter loop would make the even-odd fill above lie, so count those separately.
    for (const HalfEdgeChart &ac : atlas) {
        std::vector<Tex2> g;
        for (const Vec2f &v : ac.gutterLoop) {
            g.push_back({ v.x, v.y });
        }
        m.selfCrossingGutters += hasSelfCrossing(g, kTouchTolUv) ? 1 : 0;
    }

    m.density = static_cast<double>(cov.covered) / (static_cast<double>(kGrid) * kGrid);
    return m;
}


/// One model x cap combination with its allowed `intInGut`.  Spec 21 (a simple, self-crossing-free gutter loop)
/// brought every case to 0; the former 648e6d6 bounds (cow/128: 3, spot/128: 1) are tightened accordingly.
struct Case {
    const char *model;
    int         cap;
    long        maxInGut;
};


constexpr Case kCases[] = {
    { "cow.obj", 16, 0 },
    { "cow.obj", 32, 0 },
    { "cow.obj", 64, 0 },
    { "cow.obj", 128, 0 },
    { "spot.obj", 16, 0 },
    { "spot.obj", 32, 0 },
    { "spot.obj", 64, 0 },
    { "spot.obj", 128, 0 },
    { "bunny.iv", 16, 0 },
    { "bunny.iv", 32, 0 },
    { "bunny.iv", 64, 0 },
    { "bunny.iv", 128, 0 },
    { "horse.iv", 16, 0 },
    { "horse.iv", 32, 0 },
    { "horse.iv", 64, 0 },
    { "horse.iv", 128, 0 },
};

} // namespace

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity,readability-identifier-naming)
TEST_CASE("AI1: no two charts share a texel in the packed atlas (cow/spot/bunny/horse x 16..128)",
          "[atlas][invariant][AI1][slow]") {
    std::printf("\n%-9s %4s %6s %9s %9s %9s %9s\n", "model", "cap", "charts", "deepest", "intInGut",
                "intInInt", "density");
    for (const Case &c : kCases) {
        const Measure m = measureAtlas(c.model, c.cap);

        // Flush each row so a slow case shows progress in the ctest log.
        std::printf("%-9s %4d %6zu %9.2f %9ld %9ld %9.4f\n", c.model, c.cap, m.charts,
                    static_cast<double>(m.deepest), m.intInGut, m.intInInt, m.density);
        std::fflush(stdout);
        CAPTURE(c.model, c.cap);

        // intInInt is exact zero; intInGut only ratchets down; the gutter loops must be simple.
        CHECK(m.intInInt == 0);
        CHECK(m.intInGut <= c.maxInGut);
        CHECK(m.selfCrossingGutters == 0);
    }
}
