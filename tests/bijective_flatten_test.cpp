// This file is part of the decimate library; see decimate.h for version/license details.
/// @file bijective_flatten_test.cpp
/// @brief Bijective free-boundary flattening (Smith & Schaefer 2015) as the guard's repair.
///
/// S1 [S]  Tutte maps a synthetic disk bijectively (circle boundary, every triangle positive)
/// S2 [S]  the gradient of E_D + E_B matches central differences
/// S3 [S]  maxStep is exact: no inversion / crossing below it, an inversion or crossing right above
/// S4 [S]  every iterate of flattenBijective is injective (the interior-point property)
/// S5 [S]  a closed chart (no boundary) is refused and left untouched
/// B1 [B]  the 3-pi saddle fan (paper Figure 2 in miniature): LSCM folds it, flattenBijective makes it
///         injective with near-isometric energy
/// I1 [I]  spot/64 chart 4 (412 faces, the chart Spec 19-03 had to split): repaired in place, no split
/// I2 [I]  Flattener::LscmSplit reproduces the 19-03 behaviour on spot/64
/// B2 [B]  Spec 19 A3: the default flattenCharts is bijective on every chart of spot/cow/bunny/horse x 16/128
/// S6 [S]  start = Current on an injective chart: fewer iterations than Tutte, energy not worse

#include "model_fixture.h"
#include "bijective_detail.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstddef>
#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>
#include <numbers>
#include <string>
#include <vector>

using namespace decimate;
using namespace decimate::chart;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

/// A model file plus its half-edge mesh.
struct Loaded {
    MeshData     data;
    HalfEdgeMesh hem;
};


/// Loads @p path and builds its half-edge mesh.
Loaded loadMesh(const std::string &path) {
    Loaded l;
    l.data = load(path);
    l.hem = HalfEdgeMesh::create(static_cast<int>(l.data.positions.size()),
                                 static_cast<int>(l.data.indices.size() / 3),
                                 reinterpret_cast<const Face *>(l.data.indices.data()));
    return l;
}


/// One chart over the whole of @p d, LSCM-flattened (so the local mesh exists).
struct Single {
    MeshData      data;
    HalfEdgeMesh  hem;
    HalfEdgeChart ac;
};


/// Wraps all faces of @p d in one chart and LSCM-flattens it.
Single singleChart(MeshData d) {
    Single s;
    s.data = std::move(d);
    const int nf = static_cast<int>(s.data.indices.size() / 3);
    s.hem = HalfEdgeMesh::create(static_cast<int>(s.data.positions.size()), nf,
                                 reinterpret_cast<const Face *>(s.data.indices.data()));

    // One chart over every face; the seed is irrelevant for a whole-mesh chart.
    Chart c;
    c.id = 0;
    c.seed = 0;
    for (int f = 0; f < nf; ++f) {
        c.faces.push_back(f);
    }
    s.ac = flatten(s.hem, s.data, c);
    return s;
}


/// Apex plus `rings` rings of `n` vertices, alternating heights so each apex angle is 45 degrees:
/// 12 spokes sum to 3 pi, the "wavy cone" of the paper's Figure 2.  Developable away from the apex.
MeshData saddleFan(int n, int rings) {
    MeshData d;
    d.positions.emplace_back(0.0f, 0.0f, 0.0f);
    const double h = 0.29;
    const double r = std::sqrt(1.0 - h * h);
    for (int k = 1; k <= rings; ++k) {
        for (int i = 0; i < n; ++i) {
            const double th = 2.0 * std::numbers::pi * i / n;
            const double s = k;
            d.positions.emplace_back(static_cast<float>(s * r * std::cos(th)),
                                     static_cast<float>(s * r * std::sin(th)),
                                     static_cast<float>(s * ((i % 2) != 0 ? h : -h)));
        }
    }
    const auto id = [n](int k, int i) {
        return 1 + (k - 1) * n + (i % n);
    };
    for (int i = 0; i < n; ++i) {
        d.indices.insert(d.indices.end(), { 0, id(1, i), id(1, i + 1) });
    }
    for (int k = 1; k < rings; ++k) {
        for (int i = 0; i < n; ++i) {
            const int a = id(k, i);
            const int b = id(k, i + 1);
            const int e = id(k + 1, i);
            const int f = id(k + 1, i + 1);
            d.indices.insert(d.indices.end(), { a, e, b, b, e, f });
        }
    }
    return d;
}


/// (w + 1) x (w + 1) grid with a bump, triangulated: a disk.
MeshData bumpyGrid(int w) {
    MeshData d;
    for (int i = 0; i <= w; ++i) {
        for (int j = 0; j <= w; ++j) {
            const double x = i / static_cast<double>(w);
            const double y = j / static_cast<double>(w);
            d.positions.emplace_back(static_cast<float>(x), static_cast<float>(y),
                                     static_cast<float>(0.3 * std::sin(std::numbers::pi * x) * std::sin(std::numbers::pi * y)));
        }
    }
    const auto id = [w](int i, int j) {
        return i * (w + 1) + j;
    };
    for (int i = 0; i < w; ++i) {
        for (int j = 0; j < w; ++j) {
            d.indices.insert(d.indices.end(), { id(i, j), id(i + 1, j), id(i, j + 1), id(i, j + 1), id(i + 1, j), id(i + 1, j + 1) });
        }
    }
    return d;
}


/// Closed octahedron.
MeshData octahedron() {
    MeshData d;
    d.positions = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
    d.indices = { 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5 };
    return d;
}


/// The chart's UVs as the optimiser's flat state vector (u0, v0, u1, v1, ...).
std::vector<double> stateOf(const HalfEdgeChart &ac) {
    std::vector<double> x(2 * ac.uv.size());
    for (std::size_t i = 0; i < ac.uv.size(); ++i) {
        x[2 * i] = ac.uv[i].x;
        x[2 * i + 1] = ac.uv[i].y;
    }
    return x;
}


/// Twice the signed UV area of face @p f in state @p x (positive = not inverted).
double signedArea(const detail::Problem &p, const std::vector<double> &x, std::size_t f) {
    const double *a = &x[2 * static_cast<std::size_t>(p.tri[3 * f])];
    const double *b = &x[2 * static_cast<std::size_t>(p.tri[3 * f + 1])];
    const double *c = &x[2 * static_cast<std::size_t>(p.tri[3 * f + 2])];
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}


/// Orientation determinant of (a, b, c): positive for a counter-clockwise turn.
double orient(const double *a, const double *b, const double *c) {
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}


/// Do segments ab and cd cross properly (touching endpoints do not count)?
bool segmentsCross(const double *a, const double *b, const double *c, const double *d) {
    const double d1 = orient(c, d, a);
    const double d2 = orient(c, d, b);
    const double d3 = orient(a, b, c);
    const double d4 = orient(a, b, d);
    return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
}


/// Do boundary edges i and j (not sharing a vertex) cross in x?
bool boundaryEdgesCross(const detail::Problem &p, const std::vector<double> &x, std::size_t i, std::size_t j) {
    const int a = p.bEdge[2 * i];
    const int b = p.bEdge[2 * i + 1];
    const int c = p.bEdge[2 * j];
    const int d = p.bEdge[2 * j + 1];
    if (a == c || a == d || b == c || b == d) {
        return false;
    }
    return segmentsCross(&x[2 * static_cast<std::size_t>(a)], &x[2 * static_cast<std::size_t>(b)],
                         &x[2 * static_cast<std::size_t>(c)], &x[2 * static_cast<std::size_t>(d)]);
}


/// Do any two boundary edges cross in x?  (One flat loop over the ordered pairs.)
bool anyBoundaryCrossing(const detail::Problem &p, const std::vector<double> &x) {
    const std::size_t ne = p.bEdge.size() / 2;
    for (std::size_t k = 0; k < ne * ne; ++k) {
        const std::size_t i = k / ne;
        const std::size_t j = k % ne;
        if (j > i && boundaryEdgesCross(p, x, i, j)) {
            return true;
        }
    }
    return false;
}


/// Is x a valid interior point: every triangle positive, no two boundary edges crossing?
bool valid(const detail::Problem &p, const std::vector<double> &x) {
    for (std::size_t f = 0; f < p.ref.size(); ++f) {
        if (!(signedArea(p, x, f) > 0.0)) {
            return false;
        }
    }
    return !anyBoundaryCrossing(p, x);
}


/// @p x moved by @p t along @p d.
std::vector<double> stepped(const std::vector<double> &x, const std::vector<double> &d, double t) {
    std::vector<double> y(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        y[i] = x[i] + t * d[i];
    }
    return y;
}


/// Deterministic pseudo-random direction in [-1, 1].
std::vector<double> randomDirection(std::size_t n, unsigned seed) {
    std::vector<double> d(n);
    unsigned            s = seed;
    for (double &v : d) {
        s = s * 1664525U + 1013904223U;
        v = static_cast<double>(s >> 8) / static_cast<double>(1U << 24) * 2.0 - 1.0;
    }
    return d;
}


/// One model x cap row of the pinned B2 energy table.
struct B2Case {
    const char *model;
    int         cap;
    int         charts;    ///< Charts after flattenCharts (cow/16 gains one by a split).
    double      energySum; ///< Sum of the per-chart isometricEnergy; pinned to 1e-4 relative.
};


// Measured in build-rel (Release), Spec 19 A3.
constexpr B2Case kB2[] = {
    { "spot.obj", 16, 16, 64.220824 },
    { "spot.obj", 128, 127, 508.130213 },
    { "cow.obj", 16, 17, 69.589173 },
    { "cow.obj", 128, 117, 496.451638 },
    { "bunny.iv", 16, 15, 60.088008 },
    { "bunny.iv", 128, 115, 460.100265 },
    { "horse.iv", 16, 16, 64.211302 },
    { "horse.iv", 128, 96, 384.49569 },
};


/// Bit-for-bit double equality (the cache must not change a single ulp).
bool bitEqual(double a, double b) {
    return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}


/// Cached barrier vs the fresh-hash overload at @p x: value and every gradient slot bit-equal.
/// @returns false (after recording the failure) on the first mismatch.
bool barrierMatches(const detail::Problem &p, const std::vector<double> &x, detail::BarrierCache &cache) {
    std::vector<double> gCached(x.size(), 0.0);
    std::vector<double> gFresh(x.size(), 0.0);
    const double        cached = detail::barrier(p, x, &gCached, cache);
    const double        fresh = detail::barrier(p, x, &gFresh);
    bool                ok = bitEqual(cached, fresh);
    for (std::size_t i = 0; ok && i < x.size(); ++i) {
        ok = bitEqual(gCached[i], gFresh[i]);
    }
    return ok;
}


/// bunny/16's largest chart as the optimiser sees it, plus `steps` + 1 consecutive optimiser states
/// (one L-BFGS iteration each, every state injective).
struct BunnyTrajectory {
    MeshData                         data;
    HalfEdgeChart                    ac;
    detail::Problem                  problem;
    std::vector<std::vector<double>> states;
};


// NOLINTNEXTLINE(readability-function-size): REQUIRE expands to many statements
BunnyTrajectory bunnyTrajectory(int steps) {
    BunnyTrajectory b;

    // The largest chart has the most boundary pairs, so the barrier's hash is exercised hardest.
    const auto     l = loadMesh(modelOrSkip("bunny.iv"));
    const ChartSet cs = segment(l.hem, l.data, 16);
    std::size_t    largest = 0;
    for (std::size_t c = 1; c < cs.charts.size(); ++c) {
        if (cs.charts[c].faces.size() > cs.charts[largest].faces.size()) {
            largest = c;
        }
    }

    b.data = l.data;
    b.ac = flatten(l.hem, l.data, cs.charts[largest]);
    BijectiveStats st;
    REQUIRE(flattenBijective(b.data, b.ac, { .maxIterations = 1 }, &st)); // Tutte start, 1 iteration
    b.problem = detail::buildProblem(b.data, b.ac);

    // Record consecutive optimiser iterates: realistic states for the cache to walk along.
    b.states.push_back(stateOf(b.ac));
    for (int k = 0; k < steps; ++k) {
        REQUIRE(flattenBijective(b.data, b.ac, { .start = BijectiveStart::Current, .maxIterations = 1 }));
        b.states.push_back(stateOf(b.ac));
    }
    return b;
}

} // namespace

/// S1 — Tutte: circle boundary, every triangle positive
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S1: flattenTutte maps a bumpy grid disk bijectively onto a circle", "[chart][bijective][S1]") {
    Single s = singleChart(bumpyGrid(8));
    REQUIRE(flattenTutte(s.data, s.ac));
    CHECK(s.ac.injective());
    CHECK(s.ac.flippedTriangles == 0);

    const detail::Problem     p = detail::buildProblem(s.data, s.ac);
    const std::vector<double> x = stateOf(s.ac);
    CHECK(valid(p, x));

    // boundary on one circle: every boundary vertex at the same distance from the centroid
    double cx = 0.0;
    double cy = 0.0;
    for (const int v : p.bVert) {
        cx += x[2 * static_cast<std::size_t>(v)];
        cy += x[2 * static_cast<std::size_t>(v) + 1];
    }
    cx /= static_cast<double>(p.bVert.size());
    cy /= static_cast<double>(p.bVert.size());
    const double r0 = std::hypot(x[2 * static_cast<std::size_t>(p.bVert[0])] - cx, x[2 * static_cast<std::size_t>(p.bVert[0]) + 1] - cy);
    for (const int v : p.bVert) {
        const double r = std::hypot(x[2 * static_cast<std::size_t>(v)] - cx, x[2 * static_cast<std::size_t>(v) + 1] - cy);
        CHECK_THAT(r, WithinRel(r0, 1e-4));
    }

    // scaled to the 3D area
    double uvArea = 0.0;
    for (std::size_t f = 0; f < p.ref.size(); ++f) {
        uvArea += 0.5 * signedArea(p, x, f);
    }
    CHECK_THAT(uvArea, WithinRel(p.area, 1e-4));
}

/// S2 — gradient of E_D + E_B against central differences
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S2: the analytic gradient of E_D and E_B matches central differences", "[chart][bijective][S2]") {
    Single s = singleChart(saddleFan(12, 3));
    REQUIRE(flattenTutte(s.data, s.ac));
    const detail::Problem p = detail::buildProblem(s.data, s.ac);
    REQUIRE(p.eps > 0.0);
    // Shrink the fan's boundary towards its centroid until the barrier is active somewhere.
    std::vector<double> x = stateOf(s.ac);
    BijectiveOptions    few;
    few.maxIterations = 30;
    REQUIRE(flattenBijective(s.data, s.ac, few));
    x = stateOf(s.ac);
    for (double &v : x) {
        v *= 0.3; // compresses every boundary gap below eps: E_B > 0 on many pairs
    }
    std::vector<double> gD(x.size(), 0.0);
    std::vector<double> gB(x.size(), 0.0);
    const double        eD = detail::distortion(p, x, &gD);
    const double        eB = detail::barrier(p, x, &gB);
    REQUIRE(std::isfinite(eD));
    REQUIRE(eB > 0.0);

    // Central differences at every coordinate; the worst relative deviation is what is pinned.
    const double h = 1e-6;
    double       worstD = 0.0;
    double       worstB = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        std::vector<double> xp = x;
        std::vector<double> xm = x;
        xp[i] += h;
        xm[i] -= h;
        const double fdD = (detail::distortion(p, xp, nullptr) - detail::distortion(p, xm, nullptr)) / (2 * h);
        const double fdB = (detail::barrier(p, xp, nullptr) - detail::barrier(p, xm, nullptr)) / (2 * h);
        worstD = std::max(worstD, std::fabs(fdD - gD[i]) / (1.0 + std::fabs(fdD)));
        worstB = std::max(worstB, std::fabs(fdB - gB[i]) / (1.0 + std::fabs(fdB)));
    }

    CHECK(worstD < 1e-5);
    CHECK(worstB < 1e-5);
}

/// S3 — the flip bound is exact
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S3: maxStep is the first inversion or boundary crossing along the direction", "[chart][bijective][S3]") {
    Single s = singleChart(saddleFan(12, 3));
    REQUIRE(flattenTutte(s.data, s.ac));
    const detail::Problem     p = detail::buildProblem(s.data, s.ac);
    const std::vector<double> x = stateOf(s.ac);
    REQUIRE(valid(p, x));

    // 40 fixed pseudo-random directions; each finite bound is probed just inside and just outside.
    int finite = 0;
    for (unsigned seed = 1; seed <= 40; ++seed) {
        const std::vector<double> d = randomDirection(x.size(), seed);
        const double              t = detail::maxStep(p, x, d, 1e6);
        CAPTURE(seed, t);
        REQUIRE(t > 0.0);
        if (t >= 1e6) {
            continue; // nothing degenerates along this direction
        }
        ++finite;
        CHECK(valid(p, stepped(x, d, 0.999 * t)));
        CHECK(valid(p, stepped(x, d, 0.5 * t)));
        // Just past the singularity a triangle is inverted or two boundary edges cross.
        CHECK_FALSE(valid(p, stepped(x, d, 1.001 * t)));
    }
    CHECK(finite >= 30);
}

/// S4 — every iterate is injective
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S4: flattenBijective is injective after every iteration count", "[chart][bijective][S4]") {
    const Single s = singleChart(saddleFan(12, 3));

    // Every prefix of the run is a valid output, so each iteration count is checked on its own.
    for (int it = 0; it <= 40; ++it) {
        HalfEdgeChart    ac = s.ac;
        BijectiveOptions o;
        o.maxIterations = it;
        REQUIRE(flattenBijective(s.data, ac, o));
        CAPTURE(it, ac.flippedTriangles, ac.boundarySimple);
        CHECK(ac.injective());
        CHECK(selfOverlapPairs(ac) == 0);
    }
}

/// S5 — a closed chart is refused
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S5: a closed chart is not a disk: refused, uv untouched", "[chart][bijective][S5]") {
    Single s = singleChart(octahedron());

    // Both entry points must refuse and leave the LSCM UVs as they were.
    const auto before = s.ac.uv;
    CHECK_FALSE(flattenTutte(s.data, s.ac));
    BijectiveStats st;
    CHECK_FALSE(flattenBijective(s.data, s.ac, {}, &st));
    CHECK(s.ac.uv == before);
}

/// B1 — the saddle fan: LSCM folds, the bijective flattener does not
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B1: the 3-pi saddle fan folds under LSCM and becomes injective and near-isometric",
          "[chart][bijective][B1]") {
    Single s = singleChart(saddleFan(12, 3));
    REQUIRE(s.ac.faces.size() == 60);
    CHECK_FALSE(s.ac.injective());
    CHECK(selfOverlapPairs(s.ac) >= 1);

    // The repair itself.
    BijectiveStats st;
    REQUIRE(flattenBijective(s.data, s.ac, {}, &st));
    CHECK(s.ac.injective());
    CHECK(s.ac.flippedTriangles == 0);
    CHECK(s.ac.boundarySimple);
    CHECK(selfOverlapPairs(s.ac) == 0);
    CHECK(st.energy >= 4.0); // 4 is a perfect isometry (paper Table 1)
    CHECK(st.energy < 4.5);
    CHECK(st.maxEnergy < 6.0);
    CHECK_THAT(isometricEnergy(s.data, s.ac), WithinRel(st.energy, 1e-9));
}

/// I1 — spot/64 chart 4: repaired in place, nothing split
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("I1: spot/64 chart 4 is repaired bijectively instead of split", "[chart][bijective][I1]") {
    const auto l = loadMesh(std::string(MESH_MODEL_DIR) + "/spot.obj");
    ChartSet   cs = segment(l.hem, l.data, 64);

    // Precondition: LSCM alone folds this chart.
    REQUIRE(cs.charts[4].faces.size() == 412);
    const HalfEdgeChart lscm = flatten(l.hem, l.data, cs.charts[4]);
    CHECK_FALSE(lscm.injective());

    // Default flattening repairs it in place: same chart count, same faces.
    InjectivityStats st;
    const auto       atlas = flattenCharts(l.hem, l.data, cs, {}, &st);
    CHECK(st.selfOverlappingBefore == 1);
    CHECK(st.repairedBijectively == 1);
    CHECK(st.bijectiveFailed == 0);
    CHECK(st.splits == 0);
    CHECK(st.selfOverlappingAfter == 0);
    CHECK(cs.charts.size() == 64);
    REQUIRE(atlas.size() == 64);
    CHECK(atlas[4].injective());
    CHECK(selfOverlapPairs(atlas[4]) == 0);
    CHECK(atlas[4].faces == cs.charts[4].faces);
    CHECK(isometricEnergy(l.data, atlas[4]) < 4.1);
    CHECK(st.bijectiveSeconds > 0.0);
}

/// I2 — Flattener::LscmSplit keeps the 19-03 path
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("I2: Flattener::LscmSplit still splits spot/64 chart 4 as Spec 19-03 did", "[chart][bijective][I2]") {
    const auto       l = loadMesh(std::string(MESH_MODEL_DIR) + "/spot.obj");
    ChartSet         cs = segment(l.hem, l.data, 64);
    InjectivityStats st;
    const auto       atlas = flattenCharts(l.hem, l.data, cs, { .flattener = Flattener::LscmSplit }, &st);
    CHECK(st.selfOverlappingBefore == 1);
    CHECK(st.repairedBijectively == 0);
    CHECK(st.splits == 2);
    CHECK(st.selfOverlappingAfter == 0);
    CHECK(cs.charts.size() == 66);
    CHECK(atlas.size() == 66);
}

/// B2 — Spec 19 Amendment 3: bijective is the default flattener

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B2: default flattenCharts never raises a chart's energy over its LSCM map and leaves no fold",
          "[chart][bijective][B2][slow][release-pinned]") {
    for (const B2Case &k : kB2) {
        CAPTURE(k.model, k.cap);
        const auto       l = loadMesh(modelOrSkip(k.model));
        ChartSet         cs = segment(l.hem, l.data, k.cap);
        InjectivityStats st;
        const auto       atlas = flattenCharts(l.hem, l.data, cs, {}, &st);
        REQUIRE(atlas.size() == cs.charts.size());

        // Per chart: injective, and never worse than its own LSCM map when that map was already injective.
        double sum = 0.0;
        for (std::size_t c = 0; c < atlas.size(); ++c) {
            const HalfEdgeChart lscm = flatten(l.hem, l.data, cs.charts[c]);
            const double        e = isometricEnergy(l.data, atlas[c]);
            CHECK(atlas[c].injective());
            CHECK(selfOverlapPairs(atlas[c]) == 0);
            if (lscm.injective()) {
                CHECK(e <= isometricEnergy(l.data, lscm) + 1e-9);
            }
            sum += e;
        }

        CHECK(st.selfOverlappingAfter == 0);
        CHECK(static_cast<int>(atlas.size()) == k.charts);
        CHECK_THAT(sum, WithinRel(k.energySum, 1e-4));
    }
}

/// S6 — start = Current: an injective LSCM map is a better start than Tutte
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S6: start = Current from an injective LSCM chart needs fewer iterations than Tutte",
          "[chart][bijective][S6]") {
    const auto l = loadMesh(std::string(MESH_MODEL_DIR) + "/spot.obj");
    ChartSet   cs = segment(l.hem, l.data, 64);

    // Chart 0..63 of spot/64: take the largest chart LSCM maps injectively.
    int best = -1;
    for (std::size_t c = 0; c < cs.charts.size(); ++c) {
        const HalfEdgeChart ac = flatten(l.hem, l.data, cs.charts[c]);
        if (ac.injective() && (best < 0 || cs.charts[c].faces.size() > cs.charts[static_cast<std::size_t>(best)].faces.size())) {
            best = static_cast<int>(c);
        }
    }
    REQUIRE(best >= 0);

    const HalfEdgeChart lscm = flatten(l.hem, l.data, cs.charts[static_cast<std::size_t>(best)]);
    REQUIRE(lscm.injective());

    // Same chart, two starting points for the optimiser.
    HalfEdgeChart  fromTutte = lscm;
    HalfEdgeChart  fromLscm = lscm;
    BijectiveStats stT;
    BijectiveStats stC;
    REQUIRE(flattenBijective(l.data, fromTutte, {}, &stT));
    REQUIRE(flattenBijective(l.data, fromLscm, { .start = BijectiveStart::Current }, &stC));

    CHECK(fromLscm.injective());
    CHECK(stC.iterations < stT.iterations);
    CHECK(stC.energy <= stT.energy + 1e-3);
}

/// 22-02 S1/S2 — the barrier's cached spatial hash is result-neutral

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S1: cached barrier == fresh-hash barrier, bit for bit, along 50 optimiser states x 20 trials",
          "[chart][bijective][S1][22-02]") {
    const BunnyTrajectory b = bunnyTrajectory(50);
    REQUIRE(b.problem.eps > 0.0);

    detail::BarrierCache cache(b.problem, b.states.front()); // ONE cache for the whole walk, as in optimise
    int                  mismatches = 0;
    int                  nonzero = 0;
    for (std::size_t k = 0; k + 1 < b.states.size(); ++k) {
        const std::vector<double> &x = b.states[k];
        std::vector<double>        d(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            d[i] = 8.0 * (b.states[k + 1][i] - x[i]); // 8 steps long: the trials reach beyond one cell
        }

        // Halving trial steps along d mimic the line search.
        double t = 1.0;
        for (int trial = 0; trial < 20; ++trial, t *= 0.5) {
            const std::vector<double> xt = stepped(x, d, t);
            mismatches += barrierMatches(b.problem, xt, cache) ? 0 : 1;
        }

        std::vector<double> g(x.size(), 0.0);
        nonzero += detail::barrier(b.problem, x, &g) > 0.0 ? 1 : 0;
    }

    CHECK(mismatches == 0);
    CHECK(nonzero > 0); // not vacuous: the barrier is active on this walk
    CHECK(cache.builds() >= 1);
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S2: the hash is rebuilt exactly when a boundary vertex moved more than one cell",
          "[chart][bijective][S2][22-02]") {
    const BunnyTrajectory      b = bunnyTrajectory(10);
    const std::vector<double> &x0 = b.states[10];
    detail::BarrierCache       cache(b.problem, x0);
    const double               cell = cache.cell();
    REQUIRE(cell > 0.0);

    // A boundary vertex in the middle of the loop is the probe: its movement decides the rebuild.
    const std::size_t i = 2 * static_cast<std::size_t>(b.problem.bVert[b.problem.bVert.size() / 2]);

    std::vector<double> x = x0;
    CHECK(barrierMatches(b.problem, x, cache));
    CHECK(cache.builds() == 1);

    x[i] += 0.99 * cell; // inside the cell: padded query, no rebuild
    CHECK(barrierMatches(b.problem, x, cache));
    CHECK(cache.builds() == 1);

    x = x0;
    x[i + 1] -= 1.01 * cell; // past the cell: rebuilt at x, results still equal
    CHECK(barrierMatches(b.problem, x, cache));
    CHECK(cache.builds() == 2);

    x[i + 1] -= 0.5 * cell; // 0.5 cell from the NEW grid: no rebuild
    CHECK(barrierMatches(b.problem, x, cache));
    CHECK(cache.builds() == 2);
}

// Timing only (never in ctest): `mesh_tests "[.bench]"`.
// NOLINTNEXTLINE(readability-function-size)
TEST_CASE("22-02 bench: flattenCharts and parametrizeCorners seconds", "[.bench][22-02-bench]") {
    // Wall-clock seconds since @p t0.
    const auto since = [](std::chrono::steady_clock::time_point t0) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };
    for (const auto &[model, cap] : { std::pair { "horse.iv", 16 }, std::pair { "bunny.iv", 16 } }) {
        const auto l = loadMesh(modelOrSkip(model));
        ChartSet   cs = segment(l.hem, l.data, cap);
        const auto t0 = std::chrono::steady_clock::now();
        const auto atlas = flattenCharts(l.hem, l.data, cs);
        std::printf("[bench] flattenCharts %s/%d: %.2f s, %zu charts\n", model, cap, since(t0), atlas.size());
    }

    // Default (bijective) against the raw LSCM baseline on the same atlas.
    const auto            l = loadMesh(modelOrSkip("horse.iv"));
    std::array<double, 2> secs {};
    int                   k = 0;
    for (const Flattener f : { Flattener::Bijective, Flattener::LscmRaw }) {
        ChartSet           cs = segment(l.hem, l.data, 64);
        InjectivityOptions o;
        o.flattener = f;
        const auto t0 = std::chrono::steady_clock::now();
        const auto uv = parametrizeCorners(l.hem, l.data, cs, 256, o);
        secs[k++] = since(t0);
        CHECK(!uv.empty());
    }
    std::printf("[bench] parametrizeCorners horse/64: default %.2f s, LscmRaw %.2f s, ratio %.2f\n", secs[0], secs[1],
                secs[0] / secs[1]);
}
