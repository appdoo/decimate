// This file is part of the decimate library; see decimate.h for version/license details.
/// @file chart_injectivity_test.cpp
/// @brief Spec 19-03: exact detection of LSCM self-overlap per chart and repair by splitChart.
///
/// S1/S2 [S]  selfOverlapPairs on hand-built UV charts (fold / touching / zero-area)
/// S3     [S]  splitChart on spot/64 chart 4
/// S4     [S]  maxRounds bounds the split path (Flattener::LscmSplit); leftovers are reported, never thrown
/// S5     [S]  CentroidMode is stored in the ChartSet
/// S6     [S]  determinism on horse/64
/// I1     [I]  spot/64 -> parametrizeCorners: chart 4 (412 faces) folds, repaired bijectively, 0 folded after
/// I2     [I]  4 meshes x caps 16/32/64/128: 0 folded charts after the guard, splits only for a non-disk chart
/// B1     [B]  guard disabled: ChartSet untouched, fold only reported

#include "model_fixture.h"

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

using namespace decimate;
using namespace decimate::chart;

/// Charts of the 1e-5 cow/16 that LSCM leaves non-injective (measured, build-rel).
constexpr int kF5RawNonInjective = 12;

namespace {

/// A mesh and its half-edge structure, loaded together.
struct Loaded {
    MeshData     data;
    HalfEdgeMesh hem;
};


/// Loads a mesh file and builds its half-edge connectivity.
Loaded loadMesh(const std::string &path) {
    Loaded l;
    l.data = load(path);
    l.hem = HalfEdgeMesh::create(static_cast<int>(l.data.positions.size()),
                                 static_cast<int>(l.data.indices.size() / 3),
                                 reinterpret_cast<const Face *>(l.data.indices.data()));
    return l;
}


/// A HalfEdgeChart built by hand: triangles over `uvs.size()` local vertices with given UVs.
HalfEdgeChart handChart(const std::vector<int>                     &tri,
                        const std::vector<std::pair<float, float>> &uvs) {
    HalfEdgeChart ac;
    const int     nf = static_cast<int>(tri.size() / 3);
    ac.mesh = HalfEdgeMesh::create(static_cast<int>(uvs.size()), nf,
                                   reinterpret_cast<const Face *>(tri.data()));
    for (const auto &[u, v] : uvs) {
        ac.uv.emplace_back(u, v);
    }
    for (int f = 0; f < nf; ++f) {
        ac.faces.push_back(f);
    }
    return ac;
}


/// faceToChart and every chart's faces list describe the same partition.
bool partitionConsistent(const ChartSet &cs, int numFaces) {
    if (static_cast<int>(cs.faceToChart.size()) != numFaces) {
        return false;
    }
    std::vector<int> seen(static_cast<std::size_t>(numFaces), 0);
    for (std::size_t c = 0; c < cs.charts.size(); ++c) {
        for (const int f : cs.charts[c].faces) {
            if (cs.faceToChart[static_cast<std::size_t>(f)] != static_cast<int>(c)) {
                return false;
            }
            ++seen[static_cast<std::size_t>(f)];
        }
    }
    return std::all_of(seen.begin(), seen.end(), [](int n) { return n == 1; });
}


/// Number of chart faces reachable from the first one across sym edges inside the chart.
int connectedInChart(const HalfEdgeMesh &hem, const ChartSet &cs, int chart) {
    const auto &faces = cs.charts[static_cast<std::size_t>(chart)].faces;
    if (faces.empty()) {
        return 0;
    }
    std::vector<char> visited(cs.faceToChart.size(), 0);
    std::queue<int>   q;
    q.push(faces.front());
    visited[static_cast<std::size_t>(faces.front())] = 1;
    int n = 0;
    while (!q.empty()) {
        const int f = q.front();
        q.pop();
        ++n;
        for (int e = 0; e < 3; ++e) {
            const HeIndex s = hem.sym(static_cast<HeIndex>(f * 3 + e));
            if (s == kNoHalfedge) {
                continue;
            }
            const int nb = static_cast<int>(hem.halfedge(s).face);
            if (cs.faceToChart[static_cast<std::size_t>(nb)] != chart ||
                visited[static_cast<std::size_t>(nb)] != 0) {
                continue;
            }
            visited[static_cast<std::size_t>(nb)] = 1;
            q.push(nb);
        }
    }
    return n;
}


/// Bitwise-equal partition: same face map, same per-chart faces and seeds.
bool sameChartSet(const ChartSet &a, const ChartSet &b) {
    if (a.faceToChart != b.faceToChart || a.charts.size() != b.charts.size()) {
        return false;
    }
    for (std::size_t c = 0; c < a.charts.size(); ++c) {
        if (a.charts[c].faces != b.charts[c].faces || a.charts[c].seed != b.charts[c].seed) {
            return false;
        }
    }
    return true;
}


/// Every corner UV is finite and inside the unit square.
bool cornersInUnitSquare(const std::vector<decimate::math::Vec2f> &uv) {
    return std::all_of(uv.begin(), uv.end(), [](const auto &p) {
        return std::isfinite(p.x) && std::isfinite(p.y) && p.x >= 0.0f && p.x <= 1.0f &&
               p.y >= 0.0f && p.y <= 1.0f;
    });
}


/// Pinned injectivity-repair cases (I2).
struct InjCase {
    const char *model;
    int         cap;
    int         before;
    int         repaired;
    int         splits;
};


/// Pinned I2 expectations, measured with the implementation (TrueCentroid default).  `before` is the
/// LSCM fold count (unchanged since the Spec 19-02 re-pin); with the bijective repair (Smith & Schaefer
/// 2015) every folded disk chart is re-flattened in place, so `splits` is 0 except cow/16, whose chart 0
/// (1531 faces) has two boundary loops (it holds the cow's hole) and is still split once.
/// Under Flattener::LscmSplit the 19-02 split counts were: cow 14/6/4/3, spot 0/0/2/1, bunny 3/2/2/3,
/// horse 10/2/3/3.
constexpr InjCase kCases[] = {
    { "cow.obj", 16, 5, 4, 1 },
    { "cow.obj", 32, 3, 3, 0 },
    { "cow.obj", 64, 4, 4, 0 },
    { "cow.obj", 128, 3, 3, 0 },
    { "spot.obj", 16, 0, 0, 0 },
    { "spot.obj", 32, 0, 0, 0 },
    { "spot.obj", 64, 1, 1, 0 },
    { "spot.obj", 128, 1, 1, 0 },
    { "bunny.iv", 16, 2, 2, 0 },
    { "bunny.iv", 32, 1, 1, 0 },
    { "bunny.iv", 64, 1, 1, 0 },
    { "bunny.iv", 128, 1, 1, 0 },
    { "horse.iv", 16, 4, 4, 0 },
    { "horse.iv", 32, 2, 2, 0 },
    { "horse.iv", 64, 2, 2, 0 },
    { "horse.iv", 128, 2, 2, 0 },
};


/// Number of charts in the atlas whose UV map is not injective (the F5 statistic).
int countNonInjective(const std::vector<HalfEdgeChart> &atlas) {
    return static_cast<int>(std::count_if(atlas.begin(), atlas.end(),
                                          [](const HalfEdgeChart &a) { return !a.injective(); }));
}


/// cow/16 with every coordinate scaled by @p s (the absolute `|dT| < 1e-12f` skip in LSCM then bites).
Loaded scaledCow(float s) {
    Loaded l = loadMesh(modelOrSkip("cow.obj"));
    for (auto &p : l.data.positions) {
        p = p * s;
    }
    return l;
}

} // namespace


/// S1 — synthetic fold (written first, red before the detector existed)
TEST_CASE("S1: a 4-triangle strip folded over itself has overlapping pairs; unfolded has none",
          "[chart][injectivity][S1]") {
    const std::vector<int>                     tri = { 0, 1, 2, 1, 3, 2, 2, 3, 4, 3, 5, 4 };
    const std::vector<std::pair<float, float>> flat = { { 0, 1 }, { 0, 0 }, { 1, 1 }, { 1, 0 }, { 2, 1 }, { 2, 0 } };
    CHECK(selfOverlapPairs(handChart(tri, flat)) == 0);

    // Triangle 3 (vertices 3,5,4) is folded back onto the left part of the strip.
    auto folded = flat;
    folded[4] = { 0.2f, 0.6f };
    folded[5] = { 0.4f, 0.1f };
    CHECK(selfOverlapPairs(handChart(tri, folded)) >= 1);
}


/// S2 — touching is not overlapping; a flip is; a zero-area triangle is skipped
TEST_CASE("S2: shared edge / shared vertex are not overlaps; a flip across the edge is",
          "[chart][injectivity][S2]") {
    const std::vector<int> quad = { 0, 1, 2, 1, 3, 2 };
    CHECK(selfOverlapPairs(handChart(quad, { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 1, 1 } })) == 0);

    const std::vector<int> bowtie = { 0, 1, 2, 0, 3, 4 };
    CHECK(selfOverlapPairs(handChart(bowtie, { { 0, 0 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } })) == 0);

    // Vertex 3 flipped to the same side of edge (1,2) as vertex 0: triangle (1,3,2) lies inside
    // (0,1,2).
    CHECK(selfOverlapPairs(handChart(quad, { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 0.2f, 0.2f } })) >= 1);

    // Zero-area triangle (collinear UVs) lying inside another one: skipped.
    const std::vector<int> degenerate = { 0, 1, 2, 3, 4, 5 };
    CHECK(selfOverlapPairs(
              handChart(degenerate, { { 0, 0 }, { 4, 0 }, { 0, 4 }, { 1, 1 }, { 2, 1 }, { 3, 1 } })) == 0);
}


/// S5 — the cutter's centroid mode travels with the ChartSet
TEST_CASE("S5: segment() stores its CentroidMode in the ChartSet", "[chart][injectivity][S5]") {
    MeshData d;
    d.positions = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 } };
    d.indices = { 0, 1, 2, 1, 3, 2 };
    const HalfEdgeMesh hem =
        HalfEdgeMesh::create(4, 2, reinterpret_cast<const Face *>(d.indices.data()));
    CHECK(segment(hem, d, 2).centroidMode == CentroidMode::TrueCentroid);
    CHECK(segment(hem, d, 2, kStraightenThreshold, CentroidMode::Legacy2003Middle).centroidMode ==
          CentroidMode::Legacy2003Middle);
}


/// I1 — spot/64: chart 4 folds, the bijective repair fixes it without a split
TEST_CASE("I1: spot/64 parametrizeCorners repairs the folded chart 4", "[chart][injectivity][I1]") {
    const auto l = loadMesh(std::string(MESH_MODEL_DIR) + "/spot.obj");
    const int  nf = static_cast<int>(l.data.indices.size() / 3);
    ChartSet   cs = segment(l.hem, l.data, 64);
    REQUIRE(cs.charts.size() == 64);
    REQUIRE(cs.charts[4].faces.size() == 412);

    InjectivityStats st;
    const auto       uv = parametrizeCorners(l.hem, l.data, cs, 256, {}, &st);
    CHECK(st.selfOverlappingBefore == 1);
    CHECK(st.selfOverlappingAfter == 0);
    CHECK(st.repairedBijectively == 1); // the bijective repair keeps the chart whole (was 2 splits)
    CHECK(st.splits == 0);
    CHECK(cs.charts.size() == 64 + static_cast<std::size_t>(st.splits));
    CHECK(partitionConsistent(cs, nf));
    REQUIRE(uv.size() == static_cast<std::size_t>(3 * nf));
    CHECK(cornersInUnitSquare(uv));
}


/// B1 — guard disabled: report only
TEST_CASE("B1: LscmRaw leaves the ChartSet untouched and reports the fold",
          "[chart][injectivity][B1]") {
    const auto     l = loadMesh(std::string(MESH_MODEL_DIR) + "/spot.obj");
    ChartSet       cs = segment(l.hem, l.data, 64);
    const ChartSet before = cs;

    InjectivityStats st;
    const auto       atlas = flattenCharts(l.hem, l.data, cs, { .flattener = Flattener::LscmRaw }, &st);

    // Raw LSCM must neither repair nor split: the set stays identical, only the fold is counted.
    CHECK(sameChartSet(cs, before));
    CHECK(atlas.size() == before.charts.size());
    CHECK(st.selfOverlappingBefore == 1);
    CHECK(st.splits == 0);
    CHECK(st.selfOverlappingAfter == 1);
}


/// S3 — splitChart is a partition-preserving restricted regrow
TEST_CASE("S3: splitChart(spot/64, chart 4) splits only chart 4 into two connected pieces",
          "[chart][injectivity][S3]") {
    const auto     l = loadMesh(std::string(MESH_MODEL_DIR) + "/spot.obj");
    const int      nf = static_cast<int>(l.data.indices.size() / 3);
    ChartSet       cs = segment(l.hem, l.data, 64);
    const ChartSet before = cs;
    const int      seed = cs.charts[4].seed;
    const int      fold = (cs.charts[4].faces.back() != seed) ? cs.charts[4].faces.back()
                                                              : cs.charts[4].faces.front();

    splitChart(l.hem, l.data, cs, 4, fold, cs.centroidMode);

    // Chart 4 gains exactly one sibling, seeded at the fold face.
    REQUIRE(cs.charts.size() == before.charts.size() + 1);
    const int newId = static_cast<int>(before.charts.size());
    CHECK(cs.charts[static_cast<std::size_t>(newId)].id == newId);
    CHECK(cs.charts[static_cast<std::size_t>(newId)].seed == fold);
    CHECK(cs.charts[4].seed == seed);

    // The two halves together are exactly the old chart 4.
    std::vector<int> merged = cs.charts[4].faces;
    merged.insert(merged.end(), cs.charts[static_cast<std::size_t>(newId)].faces.begin(),
                  cs.charts[static_cast<std::size_t>(newId)].faces.end());
    std::sort(merged.begin(), merged.end());
    CHECK(merged == before.charts[4].faces);
    CHECK(connectedInChart(l.hem, cs, 4) == static_cast<int>(cs.charts[4].faces.size()));
    CHECK(connectedInChart(l.hem, cs, newId) ==
          static_cast<int>(cs.charts[static_cast<std::size_t>(newId)].faces.size()));

    // Every other chart and its face map is untouched.
    for (std::size_t c = 0; c < before.charts.size(); ++c) {
        if (c != 4) {
            CHECK(cs.charts[c].faces == before.charts[c].faces);
        }
    }
    for (int f = 0; f < nf; ++f) {
        if (before.faceToChart[static_cast<std::size_t>(f)] != 4) {
            CHECK(cs.faceToChart[static_cast<std::size_t>(f)] ==
                  before.faceToChart[static_cast<std::size_t>(f)]);
        }
    }
    CHECK(partitionConsistent(cs, nf));

    // A fold face equal to the seed, or one outside the chart, is rejected.
    CHECK_THROWS_AS(splitChart(l.hem, l.data, cs, 4, seed, cs.centroidMode), std::invalid_argument);
    CHECK_THROWS_AS(
        splitChart(l.hem, l.data, cs, 4, before.charts[5].faces.front(), cs.centroidMode),
        std::invalid_argument);
}


/// I2 — 16 cases: after == 0, before / splits pinned
TEST_CASE("I2: every mesh x cap ends with zero self-overlapping charts",
          "[chart][injectivity][I2][release-pinned]") {
    for (const InjCase &c : kCases) {
        CAPTURE(c.model, c.cap);
        const auto        l = loadMesh(modelOrSkip(c.model));
        ChartSet          cs = segment(l.hem, l.data, c.cap);
        const std::size_t initial = cs.charts.size();

        InjectivityStats st;
        const auto       atlas = flattenCharts(l.hem, l.data, cs, {}, &st);
        CHECK(st.selfOverlappingBefore == c.before);
        CHECK(st.repairedBijectively == c.repaired);
        CHECK(st.splits == c.splits);
        CHECK(st.selfOverlappingAfter == 0);
        CHECK(cs.charts.size() == initial + static_cast<std::size_t>(st.splits));
        CHECK(atlas.size() == cs.charts.size());
        CHECK(partitionConsistent(cs, static_cast<int>(l.data.indices.size() / 3)));
        for (const HalfEdgeChart &ac : atlas) {
            CHECK(selfOverlapPairs(ac) == 0);
        }
    }
}


/// S4 — rounds are bounded; leftovers are reported (the split path, Flattener::LscmSplit)
TEST_CASE("S4: maxRounds bounds the guard and leftovers are reported, not thrown",
          "[chart][injectivity][S4]") {
    const auto l = loadMesh(modelOrSkip("cow.obj"));
    {
        ChartSet         cs = segment(l.hem, l.data, 16);
        InjectivityStats st;
        REQUIRE_NOTHROW(flattenCharts(l.hem, l.data, cs, { .flattener = Flattener::LscmSplit, .maxRounds = 1 }, &st));
        CHECK(st.rounds == 1);
        CHECK(st.selfOverlappingAfter >= 1);
    }
    {
        ChartSet         cs = segment(l.hem, l.data, 16);
        const ChartSet   before = cs;
        InjectivityStats st;
        const auto       atlas = flattenCharts(l.hem, l.data, cs, { .flattener = Flattener::LscmSplit, .maxRounds = 0 }, &st);
        CHECK(st.splits == 0);
        CHECK(st.selfOverlappingAfter == st.selfOverlappingBefore);
        CHECK(sameChartSet(cs, before));
        CHECK(atlas.size() == before.charts.size());
    }
}


/// S6 — determinism
TEST_CASE("S6: two guarded runs on horse/64 give byte-equal ChartSets and UVs",
          "[chart][injectivity][S6]") {
    const auto l = loadMesh(modelOrSkip("horse.iv"));
    ChartSet   a = segment(l.hem, l.data, 64);
    ChartSet   b = segment(l.hem, l.data, 64);
    const auto uvA = parametrizeCorners(l.hem, l.data, a, 256);
    const auto uvB = parametrizeCorners(l.hem, l.data, b, 256);
    CHECK(sameChartSet(a, b));
    REQUIRE(uvA.size() == uvB.size());
    // Exact float compare: determinism means bit-identical, not close.
    bool equal = true;
    for (std::size_t i = 0; i < uvA.size(); ++i) {
        equal = equal && uvA[i].x == uvB[i].x && uvA[i].y == uvB[i].y;
    }
    CHECK(equal);
}


/// F5 (Spec 19 Amendment 4) — nonInjectiveAfter counts every chart with !injective()
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("F5: LscmRaw on a 1e-5 cow reports every non-injective chart, not just SAT overlaps",
          "[chart][injectivity][F5]") {
    // Scale the cow down so the absolute |dT| threshold in LSCM bites and charts go non-injective silently.
    const Loaded     l = scaledCow(1e-5f);
    ChartSet         cs = segment(l.hem, l.data, 16);
    InjectivityStats st;
    const auto       atlas = flattenCharts(l.hem, l.data, cs, { .flattener = Flattener::LscmRaw }, &st);
    const int        measured = countNonInjective(atlas);
    CHECK(measured == kF5RawNonInjective);
    CHECK(st.nonInjectiveAfter == measured);
    // The SAT-based statistic misses them: this is why nonInjectiveAfter exists.
    CHECK(st.selfOverlappingAfter < st.nonInjectiveAfter);
}


// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("F5: the default guard on the unscaled cow leaves no non-injective chart",
          "[chart][injectivity][F5]") {
    const Loaded     l = loadMesh(modelOrSkip("cow.obj"));
    ChartSet         cs = segment(l.hem, l.data, 16);
    InjectivityStats st;
    const auto       atlas = flattenCharts(l.hem, l.data, cs, {}, &st);
    CHECK(st.nonInjectiveAfter == countNonInjective(atlas));
    CHECK(st.nonInjectiveAfter == 0);
}


/// F6 (Spec 19 Amendment 4) — a singular normal matrix gets a ridged retry, never silent zeros
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("F6: a chart whose K is singular (zero-area triangle) still gets finite, non-trivial UVs",
          "[chart][injectivity][F6]") {
    // Quad 0,1,2,3 in z = 0 plus vertex 4 on the edge between 1 and 2: face {1,4,2} has zero area,
    // is skipped by the LSCM assembly, and leaves vertex 4 without any equation (zero rows).
    MeshData data;
    data.positions = { { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, { 1.0f, 0.5f, 0.0f } };
    data.indices = { 0, 1, 2, 0, 2, 3, 1, 4, 2 };
    const HalfEdgeMesh hem = HalfEdgeMesh::create(5, 3, reinterpret_cast<const Face *>(data.indices.data()));

    Chart ch;
    ch.id = 0;
    ch.seed = 0;
    ch.faces = { 0, 1, 2 };
    const HalfEdgeChart ac = flatten(hem, data, ch);

    // Every vertex has a finite UV; most are placed, not zeroed.
    REQUIRE(ac.uv.size() == 5u);
    int nonOrigin = 0;
    for (const auto &uv : ac.uv) {
        CHECK(std::isfinite(uv.x));
        CHECK(std::isfinite(uv.y));
        nonOrigin += (uv.x != 0.0f || uv.y != 0.0f) ? 1 : 0;
    }
    CHECK(nonOrigin >= 3);       // pin (1,1) and the free vertices 1 and 3 are placed, not zeroed
    CHECK_FALSE(ac.injective()); // the degenerate triangle is never reported as a valid map
}
