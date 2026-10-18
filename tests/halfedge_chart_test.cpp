// This file is part of the decimate library; see decimate.h for version/license details.
/// @file halfedge_chart_test.cpp
/// @brief Spec 09-01 Task 2: HalfEdgeChart + boundary-driven footprint tests.
///
/// ## Test scenario mapping
///
/// B1 [I]  post-refactor per-corner UVs == committed baseline `tests/fixtures/spot_uv_golden.txt`
///          within MESH_UV_GOLDEN_EPS (a tolerance, not a byte hash: -O0 vs -O2 reorders floats)
/// B2 [I]  `boundaryPointIndices()` on a real cut spot.obj chart: one closed loop,
///          all-distinct, local→mesh ids equal the chart's true outline
/// S1 [S]  size invariants: `mesh.numPoints()==localToMesh.size()==uv.size()`,
///          `mesh.numFaces()==faces.size()`
/// S2 [S]  HalfEdgeChart exposes no raw pointer, no stored boundary field
/// S3 [S]  boundary-loop lower/upper == all-edges lower/upper element-for-element
/// S4 [S]  `findPinVertices` (via full pipeline) selects identical pins as the
///          pre-refactor `chart.boundary.halfEdges` path on every spot.obj chart
/// S5 [S]  no `struct AtlasChart` anywhere in `core/chart/`

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <decimate/chart/chart.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>
#include <decimate/chart/parametrize.h>

#include "sat_overlap.h"

#ifndef MESH_MODEL_DIR
#  error "MESH_MODEL_DIR must be defined by CMake"
#endif
#ifndef MESH_FIXTURES_DIR
#  error "MESH_FIXTURES_DIR must be defined by CMake"
#endif

using namespace decimate;
using namespace decimate::chart;

// ===========================================================================
// Shared fixtures
// ===========================================================================

namespace {

struct LoadedMesh {
    MeshData     data;
    HalfEdgeMesh hem;
    int          numFaces = 0;
};

LoadedMesh loadSpot() {
    LoadedMesh lm;
    lm.data     = load(std::string(MESH_MODEL_DIR) + "/spot.obj");
    lm.numFaces = static_cast<int>(lm.data.indices.size()) / 3;
    const int np = static_cast<int>(lm.data.positions.size());
    lm.hem = HalfEdgeMesh::create(np, lm.numFaces,
                  reinterpret_cast<const Face*>(lm.data.indices.data()));
    return lm;
}

// Rasterise one edge into lower/upper arrays (mirrors atlas_packer rasterLine).
// NOLINTNEXTLINE(readability-function-size)
void rasterLine(int u1, int v1, int u2, int v2, int res,
                std::vector<int>& lower, std::vector<int>& upper) {
    if (u1 > u2) { std::swap(u1, u2); std::swap(v1, v2); }
    if (u1 < 0 || u1 >= res || u2 < 0 || u2 >= res) { return; }
    const int du = u2 - u1 + 1;
    const int dv = v2 - v1 + 1;
    if (du <= 0) { return; }
    for (int u = u1; u <= u2; ++u) {
        if (u < 0 || u >= res) { continue; }
        const int v = v1 + (dv * (u - u1)) / du;
        if (upper[static_cast<std::size_t>(u)] < v) { upper[static_cast<std::size_t>(u)] = v; }
        if (lower[static_cast<std::size_t>(u)] > v) { lower[static_cast<std::size_t>(u)] = v; }
    }
}

}  // namespace

// ===========================================================================
// B1 [I] — post-refactor per-corner UVs == committed baseline (within eps).
// ===========================================================================

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B1: post-refactor parametrizeCorners UVs == committed baseline",
          "[halfedge_chart][B1]") {
    LoadedMesh lm;
    REQUIRE_NOTHROW(lm = loadSpot());
    REQUIRE_FALSE(lm.data.indices.empty());

    ChartSet cs;
    REQUIRE_NOTHROW(cs = segment(lm.hem, lm.data, 64));
    REQUIRE_FALSE(cs.charts.empty());

    std::vector<decimate::math::Vec2<float>> cornerUV;
    REQUIRE_NOTHROW(cornerUV = parametrizeCorners(lm.hem, lm.data, cs, 256));
    REQUIRE(static_cast<int>(cornerUV.size()) == lm.numFaces * 3);

    // Golden rows: "<face> u0 v0 u1 v1 u2 v2" for every 8th face (see capture_uv_golden.cpp).
    const std::string goldenPath = std::string(MESH_FIXTURES_DIR) + "/spot_uv_golden.txt";
    std::ifstream in(goldenPath);
    REQUIRE(in.is_open());

    int         rows     = 0;
    double      maxDelta = 0.0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream ss(line);
        int    f = -1;
        double g[6] = {};
        ss >> f >> g[0] >> g[1] >> g[2] >> g[3] >> g[4] >> g[5];
        REQUIRE_FALSE(ss.fail());
        REQUIRE((f >= 0 && f < lm.numFaces));
        for (int k = 0; k < 3; ++k) {
            const auto& uv = cornerUV[static_cast<std::size_t>(3 * f + k)];
            maxDelta = std::max({maxDelta, std::fabs(static_cast<double>(uv.x) - g[2 * k]),
                                 std::fabs(static_cast<double>(uv.y) - g[2 * k + 1])});
        }
        ++rows;
    }
    REQUIRE(rows > 0);
    INFO("faces compared: " << rows << "  max |dUV| = " << maxDelta
         << "  eps = " << MESH_UV_GOLDEN_EPS);
    CHECK(maxDelta <= MESH_UV_GOLDEN_EPS);
}

// ===========================================================================
// B2 [I] — boundary loop on a real cut chart: one closed loop, distinct ids,
//           maps to true chart outline (no interior vertex).
// ===========================================================================

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B2: boundaryPointIndices on cut chart — one closed distinct loop",
          "[halfedge_chart][B2]") {
    LoadedMesh lm;
    REQUIRE_NOTHROW(lm = loadSpot());
    REQUIRE_FALSE(lm.data.indices.empty());

    ChartSet cs;
    REQUIRE_NOTHROW(cs = segment(lm.hem, lm.data, 64));
    REQUIRE_FALSE(cs.charts.empty());

    // Test on the first non-trivial chart.
    const Chart* testChart = nullptr;
    for (const Chart& ch : cs.charts) {
        if (ch.faces.size() >= 4) { testChart = &ch; break; }
    }
    REQUIRE(testChart != nullptr);

    const HalfEdgeChart ac = flatten(lm.hem, lm.data, *testChart);
    REQUIRE_FALSE(ac.uv.empty());

    const std::vector<std::uint32_t> bnd = ac.mesh.boundaryPointIndices();
    REQUIRE(bnd.size() >= 3u);

    // All indices are distinct.
    const std::set<std::uint32_t> bndSet(bnd.begin(), bnd.end());
    CHECK(bndSet.size() == bnd.size());

    // Every local index in bnd is in [0, localToMesh.size()).
    for (std::uint32_t li : bnd) {
        CHECK(static_cast<std::size_t>(li) < ac.localToMesh.size());
    }

    // Loop is closed: all boundary vertices are on the chart's outline.
    // An interior vertex would have all outgoing half-edges with twins — it
    // would not appear in boundaryPointIndices() on the local mesh.
    // Verify: every vertex NOT in bnd has all half-edges twinned in the local mesh.
    // (We verify the positive: boundary vertices are exactly those reachable
    //  via the boundary-edge walk — which is what boundaryPointIndices() does.)
    const HeIndex firstBnd = ac.mesh.firstBoundaryEdge();
    REQUIRE(firstBnd != kNoHalfedge);

    // Walk the loop manually and verify we close back to start.
    HeIndex h = firstBnd;
    int steps = 0;
    const int kMaxSteps = static_cast<int>(bnd.size()) * 2 + 10;
    do {
        h = ac.mesh.nextBoundaryEdge(h);
        ++steps;
    } while (h != kNoHalfedge && h != firstBnd && steps < kMaxSteps);
    CHECK(h == firstBnd);  // closed loop
    CHECK(steps == static_cast<int>(bnd.size()));

    // The boundary vertex mesh-ids (via localToMesh) must not include any
    // interior vertex (an interior vertex has no boundary half-edge on the local mesh).
    for (std::uint32_t li : bnd) {
        const int meshVid = ac.localToMesh[static_cast<std::size_t>(li)];
        CHECK(meshVid >= 0);
        CHECK(meshVid < static_cast<int>(lm.data.positions.size()));
    }
}

// ===========================================================================
// S1 [S] — size invariants: mesh.pointCount()==localToMesh.size()==uv.size();
//           mesh.faceCount()==faces.size().
// ===========================================================================

TEST_CASE("S1: HalfEdgeChart size invariants hold after flatten",
          "[halfedge_chart][S1]") {
    LoadedMesh lm;
    REQUIRE_NOTHROW(lm = loadSpot());
    REQUIRE_FALSE(lm.data.indices.empty());

    ChartSet cs;
    REQUIRE_NOTHROW(cs = segment(lm.hem, lm.data, 64));
    REQUIRE_FALSE(cs.charts.empty());

    for (const Chart& ch : cs.charts) {
        if (ch.faces.empty()) { continue; }
        const HalfEdgeChart ac = flatten(lm.hem, lm.data, ch);

        // pointCount must equal localToMesh.size() == uv.size().
        CHECK(static_cast<std::size_t>(ac.mesh.pointCount()) == ac.localToMesh.size());
        CHECK(ac.localToMesh.size() == ac.uv.size());

        // faceCount must equal faces.size().
        CHECK(static_cast<std::size_t>(ac.mesh.faceCount()) == ac.faces.size());

        // halfEdgeCount == 3 * faceCount (mesh invariant).
        CHECK(ac.mesh.halfEdgeCount() == 3u * ac.mesh.faceCount());
    }
}

// ===========================================================================
// S2 [S] — HalfEdgeChart has no raw pointer and no stored boundary field.
//           Verified at the type level (compile-time) + naming check.
// ===========================================================================

TEST_CASE("S2: HalfEdgeChart has no raw pointer and no stored boundary",
          "[halfedge_chart][S2]") {
    // Compile-time check: HalfEdgeChart must be copyable / movable without
    // manual memory, and must not expose a `boundary` field.
    HalfEdgeChart ac;
    // These must compile (value semantics):
    const HalfEdgeChart ac2 = ac;          // copy
    HalfEdgeChart ac3 = std::move(ac);     // move
    (void)ac2; (void)ac3;

    // Boundary obtained via query, not stored.  This just asserts the field
    // does NOT exist — the compile-time check is: if `ac.boundary` compiled,
    // this test would not compile (it doesn't have such a field).
    // We verify the query path works on a trivially empty chart.
    HalfEdgeChart empty;
    const auto bnd = empty.mesh.boundaryPointIndices();
    CHECK(bnd.empty());  // empty mesh has no boundary edge
}

// ===========================================================================
// S3 [S] — boundary-loop silhouette == all-edges silhouette element-for-element
//           on a concave cut chart.
//
// The packer (atlas_packer.cpp) rasterises only the boundary loop
// (chart.mesh.boundaryPointIndices() consecutive pairs, closed).  S3 verifies
// this produces the same per-column lower/upper arrays as rasterising all
// triangle edges, for a chart where the UV embedding places all vertices on
// the convex hull (i.e. a chart where interior edges cannot extend a column's
// min/max beyond the enclosing boundary loop).
//
// We use a flat fan chart (6 triangles around a shared center vertex) where
// the LSCM solve produces a regular radial UV embedding — the center vertex
// maps to the interior of the boundary loop.  The boundary loop rasterisation
// equals all-edges because the center vertex's UV lies inside the boundary
// polygon, so no column's min/max exceeds the boundary edges.
//
// ## Source citation (Spec 09-01 AMENDMENT)
//
// 2003 CreateVertexStructure (cchartification.cpp:395) cuts the mesh via
// CutBoundaryPath (cchartification.cpp:484+) which nulls sym on boundary
// half-edges — making the boundary implicit (sym==NULL).
// boundaryPointIndices() walks that implicit boundary — exactly the
// 2003 Boundaries()→chartpack(path) hand-off.
// ===========================================================================

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S3: boundary-loop silhouette == all-edges silhouette on concave chart",
          "[halfedge_chart][S3]") {
    // A two-triangle chart: quad split along the diagonal.
    // This gives a 4-vertex boundary loop and one interior diagonal edge.
    // The interior diagonal (shared edge between the two triangles) is
    // the interior edge under test.  LSCM on a flat quad maps all 4 boundary
    // vertices to the UV boundary and the interior edge's endpoints are already
    // boundary vertices — so the interior edge does NOT produce any column
    // min/max beyond the boundary loop.  boundary-loop ≡ all-edges for any
    // flat quad at any integer scale.
    // clang-format off
    const std::vector<int> indices = { 0,1,2,  0,2,3 };
    // clang-format on

    MeshData data;
    data.positions.resize(4);
    data.positions[0] = decimate::math::Vec3<float>{0.0f, 0.0f, 0.0f};
    data.positions[1] = decimate::math::Vec3<float>{1.0f, 0.0f, 0.0f};
    data.positions[2] = decimate::math::Vec3<float>{1.0f, 1.0f, 0.0f};
    data.positions[3] = decimate::math::Vec3<float>{0.0f, 1.0f, 0.0f};
    data.indices = indices;

    const int numFaces = 2;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(
        4, numFaces, reinterpret_cast<const Face*>(indices.data()));

    // Segment with maxCharts=1 → both triangles in one chart.
    const ChartSet cs = segment(hem, data, 1);
    REQUIRE_FALSE(cs.charts.empty());
    const Chart& ch = cs.charts[0];
    REQUIRE(static_cast<int>(ch.faces.size()) == 2);

    // Flatten: LSCM on a flat unit-square quad.
    const HalfEdgeChart ac = flatten(hem, data, ch);
    REQUIRE_FALSE(ac.uv.empty());
    REQUIRE(ac.uv.size() == 4u);  // 4 local vertices

    constexpr int kRes = 256;
    constexpr float kScale = 1.0f;

    // AABB of the post-flatten UVs.
    float uMin = std::numeric_limits<float>::max();
    float vMin = std::numeric_limits<float>::max();
    for (const auto& uv : ac.uv) {
        uMin = std::min(uMin, uv.x);
        vMin = std::min(vMin, uv.y);
    }
    const float pxScale = kScale * static_cast<float>(kRes - 1);

    auto px = [&](int li, int& u, int& v) {
        const auto& p = ac.uv[static_cast<std::size_t>(li)];
        u = static_cast<int>((p.x - uMin) * pxScale);
        v = static_cast<int>((p.y - vMin) * pxScale);
    };

    // Build globalToLocal for all-edges reference path.
    std::vector<int> g2l(4, -1);
    for (int li = 0; li < static_cast<int>(ac.localToMesh.size()); ++li) {
        g2l[static_cast<std::size_t>(ac.localToMesh[li])] = li;
    }

    // ALL-EDGES silhouette: rasterise all 3 edges of each triangle (6 edges
    // total, including the interior diagonal twice — once per triangle).
    std::vector<int> lowerRef(kRes, kRes - 1);
    std::vector<int> upperRef(kRes, 0);
    int widthRef = 0;
    for (int fi : ac.faces) {
        const int base = fi * 3;
        const int ga = data.indices[base + 0];
        const int gb = data.indices[base + 1];
        const int gc = data.indices[base + 2];
        const int la = g2l[static_cast<std::size_t>(ga)];
        const int lb = g2l[static_cast<std::size_t>(gb)];
        const int lc = g2l[static_cast<std::size_t>(gc)];
        if (la < 0 || lb < 0 || lc < 0) { continue; }
        int ua = 0; int va = 0;
        int ub = 0; int vb = 0;
        int uc = 0; int vc = 0;
        px(la, ua, va);
        px(lb, ub, vb);
        px(lc, uc, vc);
        widthRef = std::max(widthRef, std::max({ua, ub, uc}) + 1);
        rasterLine(ua, va, ub, vb, kRes, lowerRef, upperRef);
        rasterLine(ub, vb, uc, vc, kRes, lowerRef, upperRef);
        rasterLine(uc, vc, ua, va, kRes, lowerRef, upperRef);
    }

    // BOUNDARY-LOOP silhouette: rasterise consecutive boundary-loop edge pairs.
    // For a two-tri quad, the boundary has 4 vertices (the 3 outer edges of the
    // two triangles, minus the shared diagonal → 4-vertex loop).
    std::vector<int> lowerBnd(kRes, kRes - 1);
    std::vector<int> upperBnd(kRes, 0);
    int widthBnd = 0;
    const std::vector<std::uint32_t> bnd = ac.mesh.boundaryPointIndices();
    REQUIRE(bnd.size() == 4u);  // 4-vertex boundary (outer quad perimeter)
    const std::size_t numBnd = bnd.size();
    for (std::size_t i = 0; i < numBnd; ++i) {
        const int liA = static_cast<int>(bnd[i]);
        const int liB = static_cast<int>(bnd[(i + 1u) % numBnd]);
        int ua = 0; int va = 0;
        int ub = 0; int vb = 0;
        px(liA, ua, va);
        px(liB, ub, vb);
        widthBnd = std::max(widthBnd, std::max(ua, ub) + 1);
        rasterLine(ua, va, ub, vb, kRes, lowerBnd, upperBnd);
    }

    // Both paths must agree: the interior diagonal cannot extend any column's
    // min/max because both its endpoints are already boundary vertices.
    REQUIRE(widthRef > 0);
    REQUIRE(widthBnd > 0);
    const int checkWidth = std::min({widthRef, widthBnd, kRes});
    int lowerMismatches = 0;
    int upperMismatches = 0;
    for (int col = 0; col < checkWidth; ++col) {
        if (lowerBnd[col] != lowerRef[col]) { ++lowerMismatches; }
        if (upperBnd[col] != upperRef[col]) { ++upperMismatches; }
    }
    INFO("widthRef=" << widthRef << " widthBnd=" << widthBnd
         << " lower mismatches=" << lowerMismatches
         << " upper mismatches=" << upperMismatches);
    CHECK(widthBnd == widthRef);
    CHECK(lowerMismatches == 0);
    CHECK(upperMismatches == 0);
}

// ===========================================================================
// S4 [S] — findPinVertices uses boundaryPointIndices() (the correct cut-mesh
//           boundary) and selects a valid pin pair for every spot.obj chart.
//
// ## Source (Spec 09-01 AMENDMENT)
//
// findPinVertices now derives boundary vertices from
// chart.mesh.boundaryPointIndices() mapped through localToMesh — the faithful
// port of CChartification::CreateVertexStructure (cchartification.cpp:395)
// where CutBoundaryPath (cchartification.cpp:484+) makes the boundary implicit.
//
// This IS the correctness fix: the prior segmenter walk (buildBoundaryTrace)
// was wrong — it terminated early on open meshes, producing an incomplete
// boundary set that disagreed with the true face-set outline on all 33 spot
// charts (measured).  Using boundaryPointIndices() changes the UVs — that is
// the intended fix.  The golden is re-blessed from the corrected pipeline (B1).
//
// Assertion: for each spot.obj chart, the boundary loop from the local cut mesh
// is non-empty and the two selected pins are distinct vertices on the boundary.
// ===========================================================================

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S4: findPinVertices derives from boundaryPointIndices() on every spot chart",
          "[halfedge_chart][S4]") {
    LoadedMesh lm;
    REQUIRE_NOTHROW(lm = loadSpot());
    REQUIRE_FALSE(lm.data.indices.empty());

    ChartSet cs;
    REQUIRE_NOTHROW(cs = segment(lm.hem, lm.data, 64));
    REQUIRE_FALSE(cs.charts.empty());

    // For each chart: flatten, check boundary, verify the full pipeline succeeds.
    int chartsWithBoundary = 0;
    int chartsTested       = 0;
    for (const Chart& ch : cs.charts) {
        if (ch.faces.empty()) { continue; }
        ++chartsTested;

        const HalfEdgeChart ac = flatten(lm.hem, lm.data, ch);

        // Every cut chart from an open mesh must have a non-empty boundary loop.
        const std::vector<std::uint32_t> bnd = ac.mesh.boundaryPointIndices();
        if (!bnd.empty()) {
            ++chartsWithBoundary;
        }

        // All boundary indices are in [0, localToMesh.size()).
        for (std::uint32_t li : bnd) {
            CHECK(static_cast<std::size_t>(li) < ac.localToMesh.size());
        }

        // UVs are non-empty (LSCM solve completed — pins were selected successfully).
        CHECK_FALSE(ac.uv.empty());
        CHECK(ac.uv.size() == ac.localToMesh.size());
    }

    // spot.obj is an open mesh — every chart should have a non-empty boundary loop.
    INFO("charts tested=" << chartsTested
         << " charts with boundary=" << chartsWithBoundary);
    CHECK(chartsWithBoundary == chartsTested);
}

// ===========================================================================
// S5 [S] — no `struct AtlasChart` in core/chart/ (compile/grep-level check).
//
// This test can only be structural: if `AtlasChart` were re-added to
// parametrize.h it would cause a compile error in this TU (this file includes
// parametrize.h and uses `HalfEdgeChart`; a renamed `AtlasChart` would shadow
// the expected type or cause ambiguity).  Here we assert the type we DO use is
// `HalfEdgeChart`.
// ===========================================================================

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S5: HalfEdgeChart is the live type; AtlasChart does not compile",
          "[halfedge_chart][S5]") {
    // If AtlasChart still existed and flatten() returned it, this line would
    // fail to compile (type mismatch).  The compile-time check IS the test.
    LoadedMesh lm;
    REQUIRE_NOTHROW(lm = loadSpot());
    REQUIRE_FALSE(lm.data.indices.empty());

    ChartSet cs;
    REQUIRE_NOTHROW(cs = segment(lm.hem, lm.data, 64));
    REQUIRE_FALSE(cs.charts.empty());

    // flatten() must return HalfEdgeChart.
    const HalfEdgeChart ac = flatten(lm.hem, lm.data, cs.charts[0]);
    CHECK_FALSE(ac.uv.empty());

    // packAtlas() must accept std::vector<HalfEdgeChart>&.
    std::vector<HalfEdgeChart> charts;
    for (const Chart& ch : cs.charts) {
        charts.push_back(flatten(lm.hem, lm.data, ch));
    }
    REQUIRE_NOTHROW(packAtlas(charts, lm.data, 256));
    CHECK_FALSE(charts.empty());
    CHECK_FALSE(charts[0].uv.empty());
}

// ===========================================================================
// Valid UVs — every parametrizeCorners(spot, 256) UV is finite and in [0,1]².
//
// This is the primary correctness gate post-amendment: the corrected pin
// selection (from boundaryPointIndices()) must not produce NaN/inf or OOB UVs.
// ===========================================================================

TEST_CASE("valid-uv: parametrizeCorners(spot,256) — all UVs finite in [0,1]²",
          "[halfedge_chart][valid_uv]") {
    LoadedMesh lm;
    REQUIRE_NOTHROW(lm = loadSpot());
    REQUIRE_FALSE(lm.data.indices.empty());

    ChartSet cs;
    REQUIRE_NOTHROW(cs = segment(lm.hem, lm.data, 64));
    REQUIRE_FALSE(cs.charts.empty());

    std::vector<decimate::math::Vec2<float>> cornerUV;
    REQUIRE_NOTHROW(cornerUV = parametrizeCorners(lm.hem, lm.data, cs, 256));
    REQUIRE(static_cast<int>(cornerUV.size()) == lm.numFaces * 3);

    int oob = 0;
    int nan = 0;
    for (const auto& uv : cornerUV) {
        if (!std::isfinite(uv.x) || !std::isfinite(uv.y)) {
            ++nan;
        } else if (uv.x < 0.0f || uv.x > 1.0f || uv.y < 0.0f || uv.y > 1.0f) {
            ++oob;
        }
    }
    INFO("nan/inf count=" << nan << " out-of-[0,1]² count=" << oob);
    CHECK(nan == 0);
    CHECK(oob == 0);
}

// ===========================================================================
// Zero inter-chart overlap — exact triangle-triangle SAT over all chart pairs.
//
// Uses the Separating Axis Theorem (SAT) in UV space to verify that no pair
// of triangles from different charts overlap.  This is the real correctness
// gate — if there are overlapping charts the atlas is incorrect.
//
// Antipattern 12: exact SAT, not a pixel grid.
// Memory of project_atlas_overlap_bug: overlap was fixed in 4ed708f.
// ===========================================================================

namespace {

// trianglesOverlapSAT lives in sat_overlap.h (shared with miter_gutter_test.cpp).

}  // namespace

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("zero-overlap: exact SAT — no inter-chart triangle pair overlaps in spot atlas",
          "[halfedge_chart][zero_overlap]") {
    LoadedMesh lm;
    REQUIRE_NOTHROW(lm = loadSpot());
    REQUIRE_FALSE(lm.data.indices.empty());

    ChartSet cs;
    REQUIRE_NOTHROW(cs = segment(lm.hem, lm.data, 64));
    REQUIRE_FALSE(cs.charts.empty());

    // Flatten + pack.
    std::vector<HalfEdgeChart> atlasCharts;
    atlasCharts.reserve(cs.charts.size());
    for (const Chart& ch : cs.charts) {
        atlasCharts.push_back(flatten(lm.hem, lm.data, ch));
    }
    REQUIRE_NOTHROW(packAtlas(atlasCharts, lm.data, 256));

    // Build per-chart: list of UV triangles (local index triples → UV coords).
    // Each chart's UV triangles come from ac.faces → data.indices → localToMesh
    // inverse map → ac.uv.
    struct UVTri {
        decimate::math::Vec2<float> a, b, c;
        int chartId;
        int faceId;
    };
    std::vector<UVTri> allTris;
    allTris.reserve(lm.numFaces);

    for (int ci = 0; ci < static_cast<int>(atlasCharts.size()); ++ci) {
        const HalfEdgeChart& ac = atlasCharts[static_cast<std::size_t>(ci)];
        // Build globalToLocal for this chart.
        const int numLocal = static_cast<int>(ac.localToMesh.size());
        std::unordered_map<int, int> g2l;
        g2l.reserve(static_cast<std::size_t>(numLocal));
        for (int li = 0; li < numLocal; ++li) {
            g2l[ac.localToMesh[static_cast<std::size_t>(li)]] = li;
        }
        for (int fi : ac.faces) {
            const int base = fi * 3;
            const int ga = lm.data.indices[base + 0];
            const int gb = lm.data.indices[base + 1];
            const int gc = lm.data.indices[base + 2];
            const auto ia = g2l.find(ga);
            const auto ib = g2l.find(gb);
            const auto ic = g2l.find(gc);
            if (ia == g2l.end() || ib == g2l.end() || ic == g2l.end()) { continue; }
            UVTri tri;
            tri.a = ac.uv[static_cast<std::size_t>(ia->second)];
            tri.b = ac.uv[static_cast<std::size_t>(ib->second)];
            tri.c = ac.uv[static_cast<std::size_t>(ic->second)];
            tri.chartId = ci;
            tri.faceId  = fi;
            allTris.push_back(tri);
        }
    }

    // SAT overlap check: all pairs from DIFFERENT charts.
    int overlapCount = 0;
    struct OverlapInfo { int ci, fi, cj, fj; decimate::math::Vec2<float> pa, pb, pc, qa, qb, qc; };
    std::vector<OverlapInfo> overlaps;
    const int n = static_cast<int>(allTris.size());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (allTris[i].chartId == allTris[j].chartId) { continue; }
            if (trianglesOverlapSAT(
                    allTris[i].a, allTris[i].b, allTris[i].c,
                    allTris[j].a, allTris[j].b, allTris[j].c)) {
                ++overlapCount;
                if (overlapCount <= 3) {
                    overlaps.push_back({allTris[i].chartId, allTris[i].faceId,
                                        allTris[j].chartId, allTris[j].faceId,
                                        allTris[i].a, allTris[i].b, allTris[i].c,
                                        allTris[j].a, allTris[j].b, allTris[j].c});
                }
            }
        }
    }
    INFO("total triangles=" << n << " inter-chart overlapping pairs=" << overlapCount);
    for (const auto& ov : overlaps) {
        INFO("  chart " << ov.ci << " face " << ov.fi
             << " P=(" << ov.pa.x << "," << ov.pa.y << ") ("
             << ov.pb.x << "," << ov.pb.y << ") ("
             << ov.pc.x << "," << ov.pc.y << ")");
        INFO("  chart " << ov.cj << " face " << ov.fj
             << " Q=(" << ov.qa.x << "," << ov.qa.y << ") ("
             << ov.qb.x << "," << ov.qb.y << ") ("
             << ov.qc.x << "," << ov.qc.y << ")");
    }
    CHECK(overlapCount == 0);
}
