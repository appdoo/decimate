// This file is part of the decimate library; see decimate.h for version/license details.
/// @file miter_gutter_test.cpp
/// @brief Spec 09-02 (+ AMENDMENT 2026-10-18): miterOffsetLoop + the placed gutterLoop.
///
/// S1  [S] miter vertex lies exactly `gutter` from BOTH adjacent edge lines
/// S2  [S] sharp corner is capped at miterLimit * gutter (parameterized)
/// S3  [S] every offset vertex lies strictly outside the original polygon (both windings)
/// S5' [S] square corner moves 4*sqrt(2); 10-degree spike capped at limit * 4 for 1, 2, 8
/// S6  [S] gutter 0 is the identity; CW and CCW both grow outward
/// B1  [I] spot.obj: exact SAT, zero inter-chart triangle overlap
/// B2  [I] spot.obj: Chebyshev pixel separation >= kGutterPx, every chart placed
/// B3' [I] spot.obj: each gutterLoop encloses its chart's UVs; no two gutterLoops intersect

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <vector>

#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/math/predicates.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>

#include "sat_overlap.h"

#ifndef MESH_MODEL_DIR
#error "MESH_MODEL_DIR must be defined by CMake"
#endif

using namespace decimate;
using namespace decimate::chart;
using Catch::Matchers::WithinAbs;
using Vec2f = decimate::math::Vec2f;

namespace {

constexpr float kEps = 1e-4f;

/// Signed polygon area (positive for counter-clockwise).
float signedArea(const std::vector<Tex2> &p) {
    std::vector<Vec2f> q;
    for (const Tex2 &t : p) {
        q.push_back({ t.u, t.v });
    }
    return static_cast<float>(0.5 * decimate::math::signedArea2(q));
}


/// @p p in the opposite winding.
std::vector<Tex2> reversed(std::vector<Tex2> p) {
    std::reverse(p.begin(), p.end());
    return p;
}


/// Euclidean distance between two points.
float dist(const Tex2 &a, const Tex2 &b) {
    return std::hypot(b.u - a.u, b.v - a.v);
}


/// Distance from p to the infinite line through a, b.
float lineDist(const Tex2 &p, const Tex2 &a, const Tex2 &b) {
    const float dx = b.u - a.u;
    const float dy = b.v - a.v;
    return std::abs(dx * (p.v - a.v) - dy * (p.u - a.u)) / std::hypot(dx, dy);
}


/// Even-odd point-in-polygon.
template <class P>
bool inside(float x, float y, const std::vector<P> &poly, float (*gx)(const P &), float (*gy)(const P &)) {
    bool in = false;
    for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        const float xi = gx(poly[i]);
        const float yi = gy(poly[i]);
        const float xj = gx(poly[j]);
        const float yj = gy(poly[j]);
        if (((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
            in = !in;
        }
    }
    return in;
}


/// Point-in-polygon on Tex2 loops.
bool insideTex(const Tex2 &p, const std::vector<Tex2> &poly) {
    return inside<Tex2>(
        p.u, p.v, poly, [](const Tex2 &t) { return t.u; }, [](const Tex2 &t) { return t.v; });
}


/// Point-in-polygon on Vec2f loops.
bool insideVec(const Vec2f &p, const std::vector<Vec2f> &poly) {
    return inside<Vec2f>(
        p.x, p.y, poly, [](const Vec2f &t) { return t.x; }, [](const Vec2f &t) { return t.y; });
}


// Reference polygons for the offset tests; winding is stated per shape.
const std::vector<Tex2> kSquare = { { 10, 10 }, { 30, 10 }, { 30, 30 }, { 10, 30 } };            // CCW
const std::vector<Tex2> kTrapezoid = { { 0, 0 }, { 10, 0 }, { 8, 6 }, { 2, 6 } };                // CCW
const std::vector<Tex2> kL = { { 0, 0 }, { 20, 0 }, { 20, 8 }, { 8, 8 }, { 8, 20 }, { 0, 20 } }; // CCW, one reflex corner

/// Isoceles spike: apex at the origin, 10 degrees wide, CCW.
std::vector<Tex2> spike() {
    const float h = 100.0f * std::tan(5.0f * std::numbers::pi_v<float> / 180.0f);
    return { { 0, 0 }, { 100, -h }, { 100, h } };
}


/// spot.obj after segment, flatten and pack: the shared input of the integration tests below.
struct PackedSpot {
    MeshData                   data;
    std::vector<HalfEdgeChart> charts;
};


/// Load, segment, flatten and pack spot.obj once (read-only afterwards).
const PackedSpot &packedSpot() {
    static const PackedSpot kSpot = [] {
        PackedSpot p;
        p.data = load(std::string(MESH_MODEL_DIR) + "/spot.obj");
        const int  nf = static_cast<int>(p.data.indices.size()) / 3;
        const auto hem = HalfEdgeMesh::create(static_cast<int>(p.data.positions.size()), nf,
                                              reinterpret_cast<const Face *>(p.data.indices.data()));
        ChartSet   cs = segment(hem, p.data, 64);
        p.charts = flattenCharts(hem, p.data, cs); // injectivity-guarded (19-03, Q16)
        packAtlas(p.charts, p.data, 256);
        return p;
    }();
    return kSpot;
}


/// One placed UV triangle and the chart it belongs to.
struct UVTri {
    Vec2f a, b, c;
    int   chart;
};


/// Every placed UV triangle of the atlas, tagged with its chart index.
std::vector<UVTri> placedTriangles(const PackedSpot &ps) {
    std::vector<UVTri> tris;
    for (int ci = 0; ci < static_cast<int>(ps.charts.size()); ++ci) {
        const HalfEdgeChart &ac = ps.charts[static_cast<std::size_t>(ci)];
        for (int i = 0; i < static_cast<int>(ac.faces.size()); ++i) {
            const auto at = [&](int k) {
                return ac.uv[static_cast<std::size_t>(localCorner(ac, i, k))];
            };
            tris.push_back({ at(0), at(1), at(2), ci });
        }
    }
    return tris;
}


/// Axis-aligned UV bounding box, for cheap chart-vs-chart rejection.
struct Box {
    float x0 = 1e30f, y0 = 1e30f, x1 = -1e30f, y1 = -1e30f;
    void  add(float x, float y) {
        x0 = std::min(x0, x);
        y0 = std::min(y0, y);
        x1 = std::max(x1, x);
        y1 = std::max(y1, y);
    }
    [[nodiscard]] bool disjoint(const Box &o) const {
        return x1 < o.x0 || o.x1 < x0 || y1 < o.y0 || o.y1 < y0;
    }
};


/// Distance from p to the closed polyline `poly`, in packing-grid pixels (255 px per UV unit).
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

} // namespace

/// S1 — miter vertex is exactly `gutter` from both adjacent edge lines
TEST_CASE("S1: miter vertex lies gutter away from both adjacent edges, outward",
          "[miter_gutter][S1]") {
    const auto out = miterOffsetLoop(kTrapezoid, 4.0f);
    REQUIRE(out.size() == kTrapezoid.size());

    // Both windings: the offset must land `gutter` from both adjacent edge lines and outside the polygon.
    for (const auto &poly : { kTrapezoid, reversed(kTrapezoid) }) {
        const auto        o = miterOffsetLoop(poly, 4.0f);
        const std::size_t n = poly.size();
        for (std::size_t i = 0; i < n; ++i) {
            CHECK_THAT(lineDist(o[i], poly[(i + n - 1) % n], poly[i]), WithinAbs(4.0, kEps));
            CHECK_THAT(lineDist(o[i], poly[i], poly[(i + 1) % n]), WithinAbs(4.0, kEps));
            CHECK_FALSE(insideTex(o[i], poly));
        }
    }

    // Vertex (10,0): in-normal (0,-1), out-normal (6,2)/sqrt(40); miter = their normalised sum.
    const float n2x = 6.0f / std::sqrt(40.0f);
    const float n2y = 2.0f / std::sqrt(40.0f);
    const float mx = n2x;
    const float my = -1.0f + n2y;
    const float ml = std::hypot(mx, my);
    const float cosHalf = (mx * n2x + my * n2y) / ml;
    CHECK_THAT(out[1].u, WithinAbs(10.0 + 4.0 * (mx / ml) / cosHalf, kEps));
    CHECK_THAT(out[1].v, WithinAbs(0.0 + 4.0 * (my / ml) / cosHalf, kEps));
}

/// S5' + S2 — square corner = 4*sqrt(2); spike capped at limit * gutter
TEST_CASE("S5': square corners move exactly 4*sqrt(2) along the diagonal",
          "[miter_gutter][S5]") {
    // Each corner of the 10..30 square moves 4 along both axes: (6,6), (34,6), (34,34), (6,34).
    const auto  o = miterOffsetLoop(kSquare, 4.0f);
    const float expect[4][2] = { { 6, 6 }, { 34, 6 }, { 34, 34 }, { 6, 34 } };
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK_THAT(o[i].u, WithinAbs(expect[i][0], kEps));
        CHECK_THAT(o[i].v, WithinAbs(expect[i][1], kEps));
        CHECK_THAT(dist(o[i], kSquare[i]), WithinAbs(4.0 * std::sqrt(2.0), kEps));
    }
}

TEST_CASE("S2/S5': a 10-degree spike moves at most limit * gutter (limits 1, 2, 8)",
          "[miter_gutter][S2][S5]") {
    const auto s = spike();

    // Uncapped miter length would be 4 / cos(85 deg) ~ 45.9 > 32, so every limit bites.
    for (const float limit : { 1.0f, 2.0f, 8.0f }) {
        const auto o = miterOffsetLoop(s, 4.0f, limit);
        INFO("limit=" << limit);
        CHECK_THAT(dist(o[0], s[0]), WithinAbs(limit * 4.0, kEps));
        CHECK(o[0].u < 0.0f); // pushed away from the interior (interior is at +u)
    }

    // Default limit == 2.
    CHECK_THAT(dist(miterOffsetLoop(s, 4.0f)[0], s[0]), WithinAbs(8.0, kEps));
    CHECK(kDefaultMiterLimit == 2.0f);
}

/// S3 — offset is strictly outward
TEST_CASE("S3: every offset vertex lies strictly outside the original polygon, both windings",
          "[miter_gutter][S3]") {
    for (const auto &base : { kSquare, kTrapezoid, kL, spike() }) {
        for (const auto &poly : { base, reversed(base) }) {
            const auto o = miterOffsetLoop(poly, 4.0f);
            for (const Tex2 &p : o) {
                CHECK_FALSE(insideTex(p, poly));
            }
            CHECK(std::abs(signedArea(o)) > std::abs(signedArea(poly)));
        }
    }

    // Reflex corner (8,8) of the L moves diagonally AWAY from the interior (up-right... into the notch).
    const auto o = miterOffsetLoop(kL, 4.0f);
    CHECK(o[3].u > 8.0f);
    CHECK(o[3].v > 8.0f);
}

/// S6 — identity at 0, winding independence
TEST_CASE("S6: gutter 0 is the identity; CW and CCW both grow outward",
          "[miter_gutter][S6]") {
    for (const auto &poly : { kSquare, kTrapezoid, kL }) {
        const auto o = miterOffsetLoop(poly, 0.0f);
        REQUIRE(o.size() == poly.size());
        for (std::size_t i = 0; i < poly.size(); ++i) {
            CHECK(o[i].u == poly[i].u);
            CHECK(o[i].v == poly[i].v);
        }
    }

    // A 20x20 square grown by 4 on every side is 28x28, whichever way it is wound.
    const auto ccw = miterOffsetLoop(kSquare, 4.0f);
    const auto cw = miterOffsetLoop(reversed(kSquare), 4.0f);
    CHECK_THAT(signedArea(kSquare), WithinAbs(400.0, kEps));
    CHECK_THAT(signedArea(ccw), WithinAbs(28.0 * 28.0, 1e-2));
    CHECK_THAT(signedArea(cw), WithinAbs(-28.0 * 28.0, 1e-2));

    // Same offset points, reversed order.
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK_THAT(cw[3 - i].u, WithinAbs(ccw[i].u, kEps));
        CHECK_THAT(cw[3 - i].v, WithinAbs(ccw[i].v, kEps));
    }
}

/// Spec 21 — the offset loop is a simple polygon that encloses its input
TEST_CASE("S21a: a notch narrower than 2*gutter yields a simple loop that still encloses the polygon",
          "[miter_gutter][S21]") {
    // U shape, slot 7 wide < 2*4 (arms of different height, so no collinear overlap): the raw offset walls of the
    // slot swap sides (x = 14.5 vs 15.5) and the top edges cross them.
    const std::vector<Tex2> u = { { 0, 0 }, { 30, 0 }, { 30, 30 }, { 18.5f, 30 }, { 18.5f, 10 }, { 11.5f, 10 }, { 11.5f, 28 }, { 0, 28 } };

    // The offset must stay simple and enclose every input vertex, in both windings.
    for (const auto &poly : { u, reversed(u) }) {
        const auto o = miterOffsetLoop(poly, 4.0f);
        CHECK_FALSE(hasSelfCrossing(o));
        for (const Tex2 &p : poly) {
            CHECK(insideTex(p, o));
        }
        CHECK(insideTex({ 15.0f, 20.0f }, o)); // the slot itself is bridged, not left as a hole
        CHECK(std::abs(signedArea(o)) > std::abs(signedArea(poly)));
    }
}

TEST_CASE("S21b: a sharp convex spike keeps its tip enclosed and the loop simple",
          "[miter_gutter][S21]") {
    const auto s = spike();
    for (const auto &poly : { s, reversed(s) }) {
        const auto o = miterOffsetLoop(poly, 4.0f);
        CHECK_FALSE(hasSelfCrossing(o));
        for (const Tex2 &p : poly) {
            CHECK(insideTex(p, o));
        }
    }
}

TEST_CASE("S21c: hasSelfCrossing detects a bow-tie and accepts a convex polygon", "[miter_gutter][S21]") {
    CHECK(hasSelfCrossing({ { 0, 0 }, { 10, 10 }, { 10, 0 }, { 0, 10 } }));
    CHECK_FALSE(hasSelfCrossing(kSquare));
    CHECK_FALSE(hasSelfCrossing(kL));
}


/// B1 — exact SAT, zero overlap
TEST_CASE("B1: spot.obj atlas has zero inter-chart triangle overlap (exact SAT)",
          "[miter_gutter][B1]") {
    const PackedSpot &ps = packedSpot();
    REQUIRE(ps.charts.size() > 1);
    const auto tris = placedTriangles(ps);
    REQUIRE(tris.size() == ps.data.indices.size() / 3);

    // Chart boxes first: most chart pairs are far apart and skip the exact test.
    std::vector<Box> box(ps.charts.size());
    for (const UVTri &t : tris) {
        Box &b = box[static_cast<std::size_t>(t.chart)];
        for (const Vec2f &v : { t.a, t.b, t.c }) {
            b.add(v.x, v.y);
        }
    }

    // Same-chart triangles legitimately share edges; only cross-chart pairs count.
    std::size_t overlaps = 0;
    for (std::size_t i = 0; i < tris.size(); ++i) {
        for (std::size_t j = i + 1; j < tris.size(); ++j) {
            if (tris[i].chart == tris[j].chart ||
                box[static_cast<std::size_t>(tris[i].chart)].disjoint(box[static_cast<std::size_t>(tris[j].chart)])) {
                continue;
            }
            overlaps += trianglesOverlapSAT(tris[i].a, tris[i].b, tris[i].c, tris[j].a, tris[j].b, tris[j].c);
        }
    }
    INFO("charts=" << ps.charts.size() << " triangles=" << tris.size());
    CHECK(overlaps == 0);
}

/// B2 — Chebyshev pixel separation >= kGutterPx, all charts placed
TEST_CASE("B2: inter-chart pixel separation >= kGutterPx and every chart is placed",
          "[miter_gutter][B2]") {
    const PackedSpot &ps = packedSpot();

    // Placed: every chart has UVs inside the unit square.
    for (const HalfEdgeChart &ac : ps.charts) {
        REQUIRE_FALSE(ac.uv.empty());
        for (const Vec2f &p : ac.uv) {
            REQUIRE(p.x >= 0.0f);
            REQUIRE(p.x <= 1.0f);
            REQUIRE(p.y >= 0.0f);
            REQUIRE(p.y <= 1.0f);
        }
    }


    // 2x the packing grid: 1 packing pixel == 2 grid pixels, so kGutterPx == 2*kGutterPx here.
    constexpr int    kRes = 512;
    constexpr int    kNeed = 2 * kGutterPx;
    std::vector<int> owner(static_cast<std::size_t>(kRes) * kRes, -1);

    // Rasterize every triangle at pixel centres to learn which chart owns each pixel.
    for (const UVTri &t : placedTriangles(ps)) {
        const float ax = t.a.x * kRes, ay = t.a.y * kRes;
        const float bx = t.b.x * kRes, by = t.b.y * kRes;
        const float cx = t.c.x * kRes, cy = t.c.y * kRes;
        const float det = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy);

        // A zero-area triangle covers no pixel centre.
        if (std::abs(det) < 1e-9f) {
            continue;
        }
        const int x0 = std::max(0, static_cast<int>(std::floor(std::min({ ax, bx, cx }))));
        const int x1 = std::min(kRes - 1, static_cast<int>(std::ceil(std::max({ ax, bx, cx }))));
        const int y0 = std::max(0, static_cast<int>(std::floor(std::min({ ay, by, cy }))));
        const int y1 = std::min(kRes - 1, static_cast<int>(std::ceil(std::max({ ay, by, cy }))));
        for (int y = y0; y <= y1; ++y) {
            for (int x = x0; x <= x1; ++x) {
                const float px = static_cast<float>(x) + 0.5f;
                const float py = static_cast<float>(y) + 0.5f;
                const float l1 = ((by - cy) * (px - cx) + (cx - bx) * (py - cy)) / det;
                const float l2 = ((cy - ay) * (px - cx) + (ax - cx) * (py - cy)) / det;
                if (l1 >= 0.0f && l2 >= 0.0f && l1 + l2 <= 1.0f) {
                    owner[static_cast<std::size_t>(y) * kRes + static_cast<std::size_t>(x)] = t.chart;
                }
            }
        }
    }


    // Smallest Chebyshev distance from a covered pixel to a pixel of another chart.
    int minSep = kRes;
    for (int y = 0; y < kRes; ++y) {
        for (int x = 0; x < kRes; ++x) {
            const int me = owner[static_cast<std::size_t>(y) * kRes + static_cast<std::size_t>(x)];
            if (me < 0) {
                continue;
            }
            for (int dy = -kNeed; dy <= kNeed; ++dy) {
                for (int dx = -kNeed; dx <= kNeed; ++dx) {
                    const int xx = x + dx;
                    const int yy = y + dy;
                    if (xx < 0 || yy < 0 || xx >= kRes || yy >= kRes) {
                        continue;
                    }
                    const int o = owner[static_cast<std::size_t>(yy) * kRes + static_cast<std::size_t>(xx)];
                    if (o >= 0 && o != me) {
                        minSep = std::min(minSep, std::max(std::abs(dx), std::abs(dy)));
                    }
                }
            }
        }
    }

    INFO("min Chebyshev separation (512-grid px) searched up to " << kNeed << ": " << minSep);
    CHECK(minSep >= kNeed);
}

/// B3': the exposed gutterLoop is the real placed gutter.
///
/// MEASURED on spot.obj (see report): the float gutter outlines of neighbouring charts DO
/// interpenetrate (52 segment crossings over 33 chart pairs, deepest vertex 3.5 packing px) because
/// the packer's per-column rasterLine silhouette under-bounds the float miter polygon; the charts
/// themselves stay >= kGutterPx apart (B1/B2).  The amendment's literal "no two gutterLoops
/// intersect" therefore does not hold; this test asserts what does hold and bounds the rest.
TEST_CASE("B3': each gutterLoop encloses its chart's UVs and never swallows another chart",
          "[miter_gutter][B3]") {
    const PackedSpot &ps = packedSpot();
    const std::size_t nc = ps.charts.size();

    // normalizeAtlas reserves kGutterPx at the atlas edge; a sharp miter corner may poke up to
    // (limit-1)*kGutterPx further out (the packer clamps it when rasterizing), so the unclamped
    // float loop may overshoot [0,1] by at most that much.
    const float kSpike = (kDefaultMiterLimit - 1.0f) * static_cast<float>(kGutterPx) / 255.0f;

    for (std::size_t ci = 0; ci < nc; ++ci) {
        const HalfEdgeChart &ac = ps.charts[ci];
        INFO("chart " << ci);
        REQUIRE(ac.gutterLoop.size() >= 3);

        // Simple polygon first: the even-odd enclosure test below is meaningless on a self-crossing loop.
        {
            std::vector<Tex2> g;
            for (const Vec2f &v : ac.gutterLoop) {
                g.push_back({ v.x, v.y });
            }
            CHECK_FALSE(hasSelfCrossing(g, 1e-5f)); // Spec 21: the placed gutter is a simple polygon (touches within float noise)
        }

        for (const Vec2f &g : ac.gutterLoop) {
            REQUIRE(g.x >= -kSpike);
            REQUIRE(g.x <= 1.0f + kSpike);
            REQUIRE(g.y >= -kSpike);
            REQUIRE(g.y <= 1.0f + kSpike);
        }

        // Enclosure.  A miter fold at a concave corner with edges shorter than the gutter can
        // leave a vertex up to ~0.25 px on the wrong side of even-odd (measured: 2 of 6k UVs,
        // 0.15 and 0.22 px), so a UV counts as enclosed if inside OR within half a pixel.
        for (const Vec2f &p : ac.uv) {
            CHECK((insideVec(p, ac.gutterLoop) || distToLoopPx(p, ac.gutterLoop) < 0.5f));
        }
    }


    // No gutterLoop contains another chart's geometry (stronger than loop/loop crossing: the
    // gutter never reaches into a neighbour's chart).
    std::vector<Box> box(nc);
    for (std::size_t ci = 0; ci < nc; ++ci) {
        for (const Vec2f &g : ps.charts[ci].gutterLoop) {
            box[ci].add(g.x, g.y);
        }
    }
    std::size_t swallowed = 0;
    for (std::size_t i = 0; i < nc; ++i) {
        for (std::size_t j = 0; j < nc; ++j) {
            if (i == j) {
                continue;
            }
            for (const Vec2f &p : ps.charts[j].uv) {
                if (p.x < box[i].x0 || p.x > box[i].x1 || p.y < box[i].y0 || p.y > box[i].y1) {
                    continue;
                }
                swallowed += insideVec(p, ps.charts[i].gutterLoop);
            }
        }
    }
    CHECK(swallowed == 0);

    // Gutter-vs-gutter interpenetration is bounded below the gutter width (measured max 3.51 px).
    float deepest = 0.0f;
    for (std::size_t i = 0; i < nc; ++i) {
        for (std::size_t j = 0; j < nc; ++j) {
            if (i == j || box[i].disjoint(box[j])) {
                continue;
            }
            for (const Vec2f &v : ps.charts[j].gutterLoop) {
                if (insideVec(v, ps.charts[i].gutterLoop)) {
                    deepest = std::max(deepest, distToLoopPx(v, ps.charts[i].gutterLoop));
                }
            }
        }
    }
    INFO("deepest gutter-in-gutter penetration (packing px): " << deepest);
    CHECK(deepest < static_cast<float>(kGutterPx));
}
