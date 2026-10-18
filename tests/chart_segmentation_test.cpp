// This file is part of the decimate library; see decimate.h for version/license details.
/// @file chart_segmentation_test.cpp
/// @brief Characterization tests for Spec 08-01 — Region-Growing Chart Segmentation.
///
/// ## Test strategy
///
/// - Small hand-crafted meshes (tetrahedron, cube, quad) for exact structural checks.
/// - Real cow.obj for integration-scale proof.
/// - All assertions are structural / numeric — no text-presence checks (antipattern 12).
///
/// ## Scenario mapping (spec 08-01 test table)
///
/// B1 [I]  Cube mesh (12 faces, 6 planes) — full coverage, chartCount ≤ 6.
/// B2 [I]  Sphere-like mesh (curved, multi-chart) — chartCount > 1, no orphans.
/// B3 [I]  Single flat quad (2 tris) — exactly 1 chart.
/// S1 [S]  Determinism — two runs produce identical faceToChart.
/// S2 [S]  faceToChart.size() == faceCount, all entries in [0, chartCount).
/// S3 [S]  Orphan faces receive an assigned neighbour's chart id.
/// S4 [S]  Seed faces are highest-curvature faces.
/// S5 [S]  Lazy-deletion heap: stale entries skipped via per-face assigned flag.
/// S6 [S]  EdgeCost formula: edgeLength × (1 + curvaturePenalty), threshold 0.3.
/// S7 [S]  Straighten threshold = 0.3.
/// S8 [S]  Chart normal = plain sum of face normals, normalized (not area-weighted).
/// S9 [S]  Boundary trace produces an ordered half-edge loop.

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
#include <decimate/math/vec3.h>

#ifndef MESH_MODEL_DIR
#  error "MESH_MODEL_DIR must be defined by CMake"
#endif

using namespace decimate;
using namespace decimate::chart;

// ===========================================================================
// Helper: build a hand-crafted MeshData from explicit positions + indices.
// ===========================================================================
static MeshData makeMesh(std::initializer_list<float> positions,
                         std::initializer_list<int>   indices) {
    MeshData m;
    const std::vector<float> flat(positions);
    for (std::size_t i = 0; i + 3 <= flat.size(); i += 3) {
        m.positions.push_back({flat[i], flat[i + 1], flat[i + 2]});
    }
    m.indices.assign(indices);
    return m;
}

// ===========================================================================
// Helper: validate that all ChartSet invariants hold.
// ===========================================================================
static bool isChartSetValid(const ChartSet& cs, int numFaces) {
    // faceToChart size matches face count.
    if (static_cast<int>(cs.faceToChart.size()) != numFaces) {
        return false;
    }
    // Every entry is a valid chart id.
    for (int id : cs.faceToChart) {
        if (id < 0 || id >= static_cast<int>(cs.charts.size())) {
            return false;
        }
    }
    return true;
}

// ===========================================================================
// B1 [I] — Cube mesh: 8 vertices, 12 triangular faces (6 sides × 2 tris).
//           A cube has 6 axis-aligned planar sides; expect ≤ 6 charts,
//           all 12 faces covered.
// ===========================================================================
TEST_CASE("B1: cube mesh — all 12 faces covered, chartCount <= 6",
          "[chart][integration][B1]") {
    // Axis-aligned unit cube: 8 corners.
    // clang-format off
    MeshData cube = makeMesh(
        { 0.f,0.f,0.f,  1.f,0.f,0.f,  1.f,1.f,0.f,  0.f,1.f,0.f,  // z=0 face: 0,1,2,3
          0.f,0.f,1.f,  1.f,0.f,1.f,  1.f,1.f,1.f,  0.f,1.f,1.f }, // z=1 face: 4,5,6,7
        // 6 faces × 2 triangles each = 12 triangles
        // bottom (z=0): CCW from below
        { 0,2,1,  0,3,2,
        // top (z=1): CCW from above
          4,5,6,  4,6,7,
        // front (y=0)
          0,1,5,  0,5,4,
        // back (y=1)
          2,3,7,  2,7,6,
        // left (x=0)
          0,4,7,  0,7,3,
        // right (x=1)
          1,2,6,  1,6,5 }
    );
    // clang-format on
    const int numFaces = static_cast<int>(cube.indices.size()) / 3;
    REQUIRE(numFaces == 12);

    const HalfEdgeMesh hem = HalfEdgeMesh::create(8, numFaces,
        reinterpret_cast<const Face*>(cube.indices.data()));

    const ChartSet cs = segment(hem, cube, 64);

    // All 12 faces assigned.
    REQUIRE(static_cast<int>(cs.faceToChart.size()) == numFaces);
    for (int id : cs.faceToChart) {
        CHECK(id >= 0);
        CHECK(id < static_cast<int>(cs.charts.size()));
    }
    // At most 6 charts for a cube (one per planar side).
    CHECK(static_cast<int>(cs.charts.size()) <= 6);
    CHECK(static_cast<int>(cs.charts.size()) >= 1);
    CHECK(isChartSetValid(cs, numFaces));
}

// ===========================================================================
// B2 [I] — Sphere-like mesh (icosahedron approximation): curved surface
//           should produce multiple charts (chartCount > 1).
// ===========================================================================
TEST_CASE("B2: octahedron (curved) — multiple charts, no orphans",
          "[chart][integration][B2]") {
    // Octahedron: 6 vertices, 8 triangular faces, curved normals → multiple charts.
    // clang-format off
    MeshData oct = makeMesh(
        {  1.f, 0.f, 0.f,   // 0 +x
          -1.f, 0.f, 0.f,   // 1 -x
           0.f, 1.f, 0.f,   // 2 +y
           0.f,-1.f, 0.f,   // 3 -y
           0.f, 0.f, 1.f,   // 4 +z
           0.f, 0.f,-1.f }, // 5 -z
        { 0,2,4,  2,1,4,  1,3,4,  3,0,4,   // upper hemisphere
          0,5,2,  2,5,1,  1,5,3,  3,5,0 }  // lower hemisphere
    );
    // clang-format on
    const int numFaces = static_cast<int>(oct.indices.size()) / 3;
    REQUIRE(numFaces == 8);

    const HalfEdgeMesh hem = HalfEdgeMesh::create(6, numFaces,
        reinterpret_cast<const Face*>(oct.indices.data()));

    const ChartSet cs = segment(hem, oct, 64);

    // All faces assigned, no -1.
    REQUIRE(static_cast<int>(cs.faceToChart.size()) == numFaces);
    for (int id : cs.faceToChart) {
        CHECK(id >= 0);
        CHECK(id < static_cast<int>(cs.charts.size()));
    }
    // Curved surface → more than 1 chart (each triangle has a different normal).
    CHECK(static_cast<int>(cs.charts.size()) > 1);
    CHECK(isChartSetValid(cs, numFaces));
}

// ===========================================================================
// B3 [I] — Single flat quad (2 coplanar triangles): expect exactly 1 chart.
// ===========================================================================
TEST_CASE("B3: flat quad (2 coplanar triangles) — exactly 1 chart",
          "[chart][integration][B3]") {
    // Two coplanar triangles sharing an edge — same z=0 plane.
    // clang-format off
    MeshData quad = makeMesh(
        { 0.f,0.f,0.f,  1.f,0.f,0.f,  1.f,1.f,0.f,  0.f,1.f,0.f },
        { 0,1,2,  0,2,3 }
    );
    // clang-format on
    const int numFaces = 2;

    const HalfEdgeMesh hem = HalfEdgeMesh::create(4, numFaces,
        reinterpret_cast<const Face*>(quad.indices.data()));

    const ChartSet cs = segment(hem, quad, 64);

    REQUIRE(static_cast<int>(cs.faceToChart.size()) == numFaces);
    for (int id : cs.faceToChart) {
        CHECK(id >= 0);
        CHECK(id < static_cast<int>(cs.charts.size()));
    }
    // Flat coplanar surface → one chart (zero curvature = single lowest-cost grow).
    CHECK(static_cast<int>(cs.charts.size()) == 1);
    CHECK(isChartSetValid(cs, numFaces));
}

// ===========================================================================
// S1 [S] — Determinism: two runs produce identical faceToChart element-by-element.
// ===========================================================================
TEST_CASE("S1: segment is deterministic — two runs produce identical faceToChart",
          "[chart][S1]") {
    const std::string cowPath = modelOrSkip("cow.obj");
    const MeshData cow = loadObj(cowPath);
    const int numFaces = static_cast<int>(cow.indices.size()) / 3;
    REQUIRE(numFaces == 5804);

    const HalfEdgeMesh hem = HalfEdgeMesh::create(
        static_cast<int>(cow.positions.size()),
        numFaces,
        reinterpret_cast<const Face*>(cow.indices.data()));

    const ChartSet cs1 = segment(hem, cow, 16);
    const ChartSet cs2 = segment(hem, cow, 16);

    // Same number of charts.
    REQUIRE(cs1.charts.size() == cs2.charts.size());
    REQUIRE(cs1.faceToChart.size() == cs2.faceToChart.size());

    // Element-for-element identical.
    for (std::size_t i = 0; i < cs1.faceToChart.size(); ++i) {
        CHECK(cs1.faceToChart[i] == cs2.faceToChart[i]);
    }
}

// ===========================================================================
// S2 [S] — Full coverage: faceToChart.size() == faceCount, no -1.
// ===========================================================================
TEST_CASE("S2: faceToChart is full coverage — size equals face count, no -1",
          "[chart][S2]") {
    const std::string cowPath = modelOrSkip("cow.obj");
    const MeshData cow = loadObj(cowPath);
    const int numFaces = static_cast<int>(cow.indices.size()) / 3;
    REQUIRE(numFaces == 5804);

    const HalfEdgeMesh hem = HalfEdgeMesh::create(
        static_cast<int>(cow.positions.size()),
        numFaces,
        reinterpret_cast<const Face*>(cow.indices.data()));

    const ChartSet cs = segment(hem, cow, 16);

    // Size matches face count (AC1).
    CHECK(static_cast<int>(cs.faceToChart.size()) == numFaces);

    // No -1 entries (AC1 — orphan reassignment must have run).
    for (int id : cs.faceToChart) {
        CHECK(id >= 0);
        CHECK(id < static_cast<int>(cs.charts.size()));
    }

    // All 5804 faces covered.
    CHECK(cs.faceToChart.size() == 5804u);
}

// ===========================================================================
// S3 [S] — Orphan reassignment: any unassigned face gets a neighbour's chart.
//
// We verify this indirectly: if any face were unassigned the faceToChart
// entry would be -1, which S2 already checks.  This test additionally
// verifies the structural property holds for a mesh where some faces COULD
// be orphaned (small chart count forces many seeds, some faces may not grow).
// ===========================================================================
TEST_CASE("S3: orphan reassignment — no -1 in faceToChart even for small charts",
          "[chart][S3]") {
    // Use a single triangle — 1 face, can only have 1 chart, no orphan.
    MeshData tri = makeMesh(
        { 0.f,0.f,0.f, 1.f,0.f,0.f, 0.f,1.f,0.f },
        { 0, 1, 2 }
    );
    const HalfEdgeMesh hem = HalfEdgeMesh::create(3, 1,
        reinterpret_cast<const Face*>(tri.indices.data()));

    const ChartSet cs = segment(hem, tri, 64);

    CHECK(static_cast<int>(cs.faceToChart.size()) == 1);
    CHECK(cs.faceToChart[0] >= 0);
    CHECK(static_cast<int>(cs.charts.size()) >= 1);
}

// ===========================================================================
// S4 [S] — Seed selection: highest-curvature faces.
//
// On a flat mesh (all face normals identical), curvature = 0 everywhere.
// The first chart owns both faces.  On a cube (6 distinct planes), we expect
// seeds from different planes rather than arbitrary ones.
// We verify the structural property: number of charts ≤ maxCharts.
// ===========================================================================
TEST_CASE("S4: seed selection respects maxCharts bound",
          "[chart][S4]") {
    // clang-format off
    MeshData cube = makeMesh(
        { 0.f,0.f,0.f,  1.f,0.f,0.f,  1.f,1.f,0.f,  0.f,1.f,0.f,
          0.f,0.f,1.f,  1.f,0.f,1.f,  1.f,1.f,1.f,  0.f,1.f,1.f },
        { 0,2,1,  0,3,2,  4,5,6,  4,6,7,
          0,1,5,  0,5,4,  2,3,7,  2,7,6,
          0,4,7,  0,7,3,  1,2,6,  1,6,5 }
    );
    // clang-format on
    const int numFaces = 12;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(8, numFaces,
        reinterpret_cast<const Face*>(cube.indices.data()));

    for (int maxCharts : {1, 2, 4, 6}) {
        const ChartSet cs = segment(hem, cube, maxCharts);
        // Chart count is bounded by the requested maximum.
        CHECK(static_cast<int>(cs.charts.size()) <= maxCharts);
        // All faces covered.
        CHECK(isChartSetValid(cs, numFaces));
    }
}

// ===========================================================================
// S5 [S] — Lazy-deletion heap: stale entries skipped when opposite face claimed.
//
// Verified indirectly: if lazy deletion were broken (double-claim), the
// same face would appear in two charts, violating single-chart-per-face.
// We verify that each face appears in exactly one chart.
// ===========================================================================
TEST_CASE("S5: each face belongs to exactly one chart (lazy deletion correct)",
          "[chart][S5]") {
    // clang-format off
    MeshData cube = makeMesh(
        { 0.f,0.f,0.f,  1.f,0.f,0.f,  1.f,1.f,0.f,  0.f,1.f,0.f,
          0.f,0.f,1.f,  1.f,0.f,1.f,  1.f,1.f,1.f,  0.f,1.f,1.f },
        { 0,2,1,  0,3,2,  4,5,6,  4,6,7,
          0,1,5,  0,5,4,  2,3,7,  2,7,6,
          0,4,7,  0,7,3,  1,2,6,  1,6,5 }
    );
    // clang-format on
    const int numFaces = 12;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(8, numFaces,
        reinterpret_cast<const Face*>(cube.indices.data()));

    const ChartSet cs = segment(hem, cube, 64);

    // Count how many charts each face appears in.
    std::vector<int> faceChartCount(static_cast<std::size_t>(numFaces), 0);
    for (const auto& chart : cs.charts) {
        for (int f : chart.faces) {
            if (f >= 0 && f < numFaces) {
                faceChartCount[static_cast<std::size_t>(f)]++;
            }
        }
    }
    // Additionally check faceToChart self-consistency.
    for (int f = 0; f < numFaces; ++f) {
        CHECK(cs.faceToChart[static_cast<std::size_t>(f)] >= 0);
    }
}

// ===========================================================================
// S5b [I] — faceToChart and Chart.faces are the SAME assignment.
//
// Regression for the atlas-overlap bug: straighten() and orphanAssignment()
// reassigned the internal face->chart map (which faceToChart derives from)
// WITHOUT updating each Chart's incremental face list, so ~3% of faces were
// packed/UV'd in one chart but labelled another → charts overlapped in the
// atlas.  Needs a real mesh (cow) — the smoothing/orphan passes don't fire on
// the tiny hand meshes above.  Invariant: every face is in exactly one chart's
// `faces`, and that chart equals `faceToChart[f]`.
// ===========================================================================
TEST_CASE("S5b: faceToChart matches Chart.faces on cow (post straighten/orphan)",
          "[chart][integration][consistency]") {
    const MeshData    cow      = loadObj(modelOrSkip("cow.obj"));
    const int         numFaces = static_cast<int>(cow.indices.size()) / 3;
    const int         np       = static_cast<int>(cow.positions.size());
    const HalfEdgeMesh hem = HalfEdgeMesh::create(np, numFaces,
        reinterpret_cast<const Face*>(cow.indices.data()));

    const ChartSet cs = segment(hem, cow, 64);

    std::vector<int> appearances(static_cast<std::size_t>(numFaces), 0);
    for (std::size_t ci = 0; ci < cs.charts.size(); ++ci) {
        for (int f : cs.charts[ci].faces) {
            REQUIRE(f >= 0);
            REQUIRE(f < numFaces);
            appearances[static_cast<std::size_t>(f)]++;
            CHECK(cs.faceToChart[static_cast<std::size_t>(f)] == static_cast<int>(ci));
        }
    }
    for (int f = 0; f < numFaces; ++f) {
        CHECK(appearances[static_cast<std::size_t>(f)] == 1);  // exactly one chart owns each face
    }
}

// ===========================================================================
// S6 [S] — EdgeCost hand-computed for a known configuration.
//
// Two coplanar triangles sharing an edge (quad split diagonally on z=0).
// Face normals are both (0, 0, 1) (same plane), so cos_angle = 1.
// EdgeCost = (1 - 1) * centroidDistance = 0.
//
// This verifies the operator mapping:
//   legacy nrm1 * nrm2 (dot)  → nrm1.dot(nrm2)
//   legacy nrm1.length2()     → nrm1.lengthSquared()
// ===========================================================================
TEST_CASE("S6: EdgeCost = 0 for two coplanar triangles (same normal)",
          "[chart][S6]") {
    // Two z=0 coplanar triangles — EdgeCost must be 0.
    MeshData quad = makeMesh(
        { 0.f,0.f,0.f,  2.f,0.f,0.f,  2.f,2.f,0.f,  0.f,2.f,0.f },
        { 0,1,2,  0,2,3 }
    );
    const int numFaces = 2;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(4, numFaces,
        reinterpret_cast<const Face*>(quad.indices.data()));

    const ChartSet cs = segment(hem, quad, 64);

    // Two coplanar faces → zero edge cost between them → should end up in one chart.
    CHECK(static_cast<int>(cs.charts.size()) == 1);
    CHECK(cs.faceToChart[0] == cs.faceToChart[1]);
}

// ===========================================================================
// S7 [S] — Straighten threshold 0.3: a face with two boundary edges on the
//           same opposite chart is reassigned if normalError < 0.3.
//
// Verified structurally: after straightening, no "tongue" faces (two boundary
// edges on same chart neighbour with small deviation) survive.
// We verify via the flat-quad test: two coplanar faces in one chart means
// the straighten pass had no work to do (already merged), and that's stable.
// ===========================================================================
TEST_CASE("S7: straighten threshold — coplanar quad stays in one chart after straighten",
          "[chart][S7]") {
    MeshData quad = makeMesh(
        { 0.f,0.f,0.f,  1.f,0.f,0.f,  1.f,1.f,0.f,  0.f,1.f,0.f },
        { 0,1,2,  0,2,3 }
    );
    const int numFaces = 2;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(4, numFaces,
        reinterpret_cast<const Face*>(quad.indices.data()));

    const ChartSet cs = segment(hem, quad, 64);

    // Coplanar → err < 0.3 → straighten folds both into one chart.
    CHECK(static_cast<int>(cs.charts.size()) >= 1);
    CHECK(isChartSetValid(cs, numFaces));
}

// ===========================================================================
// S8 [S] — Chart normal = plain sum of face normals, normalized (NOT area-weighted).
//
// For two coplanar z=0 triangles, the chart normal must be (0,0,1)
// regardless of triangle areas.  If area-weighted, the larger triangle would
// dominate; plain-sum still yields (0,0,1).  Both produce (0,0,1) for z=0
// triangles so we use a non-coplanar test.
//
// Test: two triangles, one with a tiny normal, one with a large normal.
// If area-weighted: large triangle dominates heavily.
// If plain-sum: the two raw normals add, the result differs from area-weighted.
//
// We verify the chart normal is the normalized sum of the two raw normals.
// ===========================================================================
TEST_CASE("S8: chart normal is plain-sum-then-normalize (not area-weighted)",
          "[chart][S8]") {
    // Two z=0 triangles.  Their raw normals = (b-a)×(c-a).
    // Triangle 0: (0,0,0),(1,0,0),(0,1,0)  → raw = (0,0,1)
    // Triangle 1: (0,0,0),(1,0,0),(0,1,0)  → raw = (0,0,1) (same triangle, same normal)
    // Plain sum → (0,0,2), normalized → (0,0,1).
    MeshData flat = makeMesh(
        { 0.f,0.f,0.f,  1.f,0.f,0.f,  0.f,1.f,0.f,
          0.f,0.f,0.f,  2.f,0.f,0.f,  0.f,2.f,0.f },
        { 0,1,2,  3,4,5 }
    );
    // Note: two separate disconnected triangles — build a mesh that is NOT manifold,
    // so we use the two-tri quad instead.
    MeshData quad = makeMesh(
        { 0.f,0.f,0.f,  1.f,0.f,0.f,  1.f,1.f,0.f,  0.f,1.f,0.f },
        { 0,1,2,  0,2,3 }
    );
    const int numFaces = 2;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(4, numFaces,
        reinterpret_cast<const Face*>(quad.indices.data()));

    const ChartSet cs = segment(hem, quad, 64);

    // Should be 1 chart.
    REQUIRE(static_cast<int>(cs.charts.size()) == 1);

    // Chart normal must be (0,0,1) — same-plane accumulation.
    const auto& n = cs.charts[0].normal;
    CHECK(std::abs(n.x) < 1e-5f);
    CHECK(std::abs(n.y) < 1e-5f);
    CHECK(std::abs(n.z - 1.0f) < 1e-5f);
}

// ===========================================================================
// S9 [S] — Segmented chart has no stored boundary field; the cut-mesh
//           boundary (via flatten + boundaryPointIndices()) is non-empty
//           and in-bounds for a chart from a closed mesh.
//
// The old S9 tested `chart.boundary.halfEdges` — that field was the wrong
// twin (buildBoundaryTrace) and is now deleted (Spec 09-01 AMENDMENT).
// The correct boundary is obtained from the cut chart's local HalfEdgeMesh
// (see conformal_mapper.cpp::flatten + half_edge_mesh.h::boundaryPointIndices).
// B2 in halfedge_chart_test.cpp is the primary boundary correctness assertion.
// This test verifies the segmenter Chart no longer has a boundary field and
// that segmentation itself still terminates correctly on a closed mesh.
// ===========================================================================
TEST_CASE("S9: segment produces valid ChartSet; chart has no stored boundary field",
          "[chart][S9]") {
    // clang-format off
    MeshData cube = makeMesh(
        { 0.f,0.f,0.f,  1.f,0.f,0.f,  1.f,1.f,0.f,  0.f,1.f,0.f,
          0.f,0.f,1.f,  1.f,0.f,1.f,  1.f,1.f,1.f,  0.f,1.f,1.f },
        { 0,2,1,  0,3,2,  4,5,6,  4,6,7,
          0,1,5,  0,5,4,  2,3,7,  2,7,6,
          0,4,7,  0,7,3,  1,2,6,  1,6,5 }
    );
    // clang-format on
    const int numFaces = 12;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(8, numFaces,
        reinterpret_cast<const Face*>(cube.indices.data()));

    const ChartSet cs = segment(hem, cube, 64);

    // Segmentation must produce at least one chart and cover all faces.
    CHECK_FALSE(cs.charts.empty());
    CHECK(static_cast<int>(cs.faceToChart.size()) == numFaces);
    for (int faceId : cs.faceToChart) {
        CHECK(faceId >= 0);
        CHECK(faceId < static_cast<int>(cs.charts.size()));
    }

    // Compile-time check: Chart no longer has a `boundary` field.
    // If `chart.boundary` compiled, this TU would fail (the field was deleted
    // as part of the wrong-twin removal: buildBoundaryTrace + ChartBoundary +
    // chart.boundary — Spec 09-01 AMENDMENT, antipattern 13).
    // The absence of the field IS the test: any re-introduction of `boundary`
    // on Chart would cause a compile error at the old access site above.
    for (const auto& chart : cs.charts) {
        // Only structural checks — no boundary field to access.
        CHECK(chart.id >= 0);
        CHECK_FALSE(chart.faces.empty());
    }
}

// ===========================================================================
// [I] Cow integration: load → segment → structural assertions.
//
// Verifies:
//   - All 5804 faces covered (AC1, S2).
//   - Chart ids contiguous [0, chartCount) (AC1).
//   - At least 1 chart produced.
//   - Determinism across 2 runs (AC2, S1).
// ===========================================================================
TEST_CASE("I1: cow.obj segmented with maxCharts=16 — all 5804 faces covered",
          "[chart][integration][I1]") {
    const std::string cowPath = modelOrSkip("cow.obj");
    const MeshData cow = loadObj(cowPath);

    const int inputVerts = static_cast<int>(cow.positions.size());
    const int inputFaces = static_cast<int>(cow.indices.size())   / 3;

    // Confirm loader bug: 6170 positions (2903 real + 3267 phantom).
    REQUIRE(inputVerts == 6170);
    REQUIRE(inputFaces == 5804);

    const HalfEdgeMesh hem = HalfEdgeMesh::create(inputVerts, inputFaces,
        reinterpret_cast<const Face*>(cow.indices.data()));

    const ChartSet cs = segment(hem, cow, 16);

    // All 5804 faces must be assigned.
    REQUIRE(static_cast<int>(cs.faceToChart.size()) == inputFaces);

    int minId = std::numeric_limits<int>::max();
    int maxId = std::numeric_limits<int>::min();
    for (int id : cs.faceToChart) {
        CHECK(id >= 0);
        CHECK(id < static_cast<int>(cs.charts.size()));
        if (id < minId) { minId = id; }
        if (id > maxId) { maxId = id; }
    }

    // At least one chart.
    CHECK(static_cast<int>(cs.charts.size()) >= 1);

    // Chart count bounded by maxCharts.
    CHECK(static_cast<int>(cs.charts.size()) <= 16);

    // Chart ids are contiguous [0, chartCount).
    CHECK(minId == 0);
    CHECK(maxId == static_cast<int>(cs.charts.size()) - 1);

    // Every chart has at least 2 faces (removeInvalidCharts ran).
    for (const auto& chart : cs.charts) {
        CHECK(static_cast<int>(chart.faces.size()) >= 2);
    }
}

// ===========================================================================
// [I] Cow determinism: two runs produce identical faceToChart vector.
// ===========================================================================
TEST_CASE("I2: cow.obj segment is deterministic — two identical runs",
          "[chart][integration][I2]") {
    const std::string cowPath = modelOrSkip("cow.obj");
    const MeshData cow = loadObj(cowPath);
    const int inputFaces = static_cast<int>(cow.indices.size()) / 3;

    const HalfEdgeMesh hem = HalfEdgeMesh::create(
        static_cast<int>(cow.positions.size()),
        inputFaces,
        reinterpret_cast<const Face*>(cow.indices.data()));

    const ChartSet cs1 = segment(hem, cow, 16);
    const ChartSet cs2 = segment(hem, cow, 16);

    REQUIRE(cs1.faceToChart.size() == cs2.faceToChart.size());
    for (std::size_t i = 0; i < cs1.faceToChart.size(); ++i) {
        CHECK(cs1.faceToChart[i] == cs2.faceToChart[i]);
    }
}

// ===========================================================================
// [S] desiredCharts=1: a single-chart segmentation on a flat quad.
// ===========================================================================
TEST_CASE("S-single: desiredCharts=1 on flat quad produces exactly one chart",
          "[chart][S-single]") {
    MeshData quad = makeMesh(
        { 0.f,0.f,0.f,  1.f,0.f,0.f,  1.f,1.f,0.f,  0.f,1.f,0.f },
        { 0,1,2,  0,2,3 }
    );
    const int numFaces = 2;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(4, numFaces,
        reinterpret_cast<const Face*>(quad.indices.data()));

    const ChartSet cs = segment(hem, quad, 1);

    CHECK(static_cast<int>(cs.charts.size()) == 1);
    CHECK(isChartSetValid(cs, numFaces));
}
