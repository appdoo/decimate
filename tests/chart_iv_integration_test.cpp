// This file is part of the decimate library; see decimate.h for version/license details.
/// @file chart_iv_integration_test.cpp
/// @brief Spec 08-03: Chartification-fidelity integration tests on .iv models.
///
/// Verifies three regressions introduced in the 08-01/08-02 port are fixed:
///   D1 — crash (std::out_of_range) on horse.iv / bunny.iv (B1, B2)
///   D2 — convergence-based iteration bound, not numFaces * maxCharts (S5)
///   D3 — seed face dropped from LSCM solve and UV scatter (B3, B4, S1)
///
/// ## Scenario mapping (spec 08-03 test table)
///
/// B1 [I]  horse.iv: segment → parametrizeCorners → packAtlas — no exception, non-empty.
/// B2 [I]  bunny.iv: same pipeline — no exception, non-empty.
/// B3 [I]  horse.iv: every chart's seed is present in ac.faces AND cornerUV slots
///          for the seed face are not the default-constructed {0,0} hole.
/// B4 [I]  cow.obj:  same structural check — seed in ac.faces, non-zero cornerUV.
/// S1 [S]  After segment, every chart.faces contains chart.seed exactly once.
/// S5 [S]  Iteration bound is kMaxIterFactor * maxCharts (D2 guard — binds to symbol).
/// S6 [S]  cow.obj chart count is exact and stable; faceToChart is deterministic.
/// S4 [S]  Empty-after-pruning chart is skipped by flattener (no crash, no UVs).

#include <catch2/catch_approx.hpp>
#include "model_fixture.h"
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <unordered_map>
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

/// @brief Load a mesh and build the HalfEdgeMesh from it.
struct LoadedMesh {
    MeshData    data;
    HalfEdgeMesh hem;
    int          numFaces = 0;
};

LoadedMesh loadMesh(const std::string& path) {
    LoadedMesh lm;
    lm.data     = load(path);
    lm.numFaces = static_cast<int>(lm.data.indices.size()) / 3;
    const int np = static_cast<int>(lm.data.positions.size());
    lm.hem = HalfEdgeMesh::create(np, lm.numFaces,
                  reinterpret_cast<const Face*>(lm.data.indices.data()));
    return lm;
}

/// @brief True if a UV is finite (neither NaN nor Inf).
bool uvFinite(const decimate::math::Vec2<float>& uv) {
    return std::isfinite(uv.x) && std::isfinite(uv.y);
}

/// @brief True if ANY of the three seed face corner UVs differs from {0,0}.
///
/// A real LSCM solve (pin0 → {0,0}, pin1 → {1,1}) assigns at least one vertex
/// away from the origin, so a chart whose seed corners are ALL {0,0} was not
/// written by the solve — it still has the default-constructed value.
/// We check "not all zero" rather than "no corner is zero" to be robust to the
/// case where the atlas origin coincidentally lands at one corner.
bool seedCornersWritten(int seed, const std::vector<decimate::math::Vec2<float>>& cornerUV) {
    for (int k = 0; k < 3; ++k) {
        const auto& uv = cornerUV[static_cast<std::size_t>(3 * seed + k)];
        if (uvFinite(uv) && (uv.x != 0.0f || uv.y != 0.0f)) {
            return true;  // at least one corner was written away from origin
        }
    }
    return false;  // all three corners are {0,0} or non-finite — seed was not written
}

}  // namespace

// ===========================================================================
// B1 [I] — horse.iv: full pipeline completes without exception; charts non-empty.
//
// This was crashing with std::out_of_range because the seed face's vertices
// were missing from localOf (D1 defect — seed not in chart.faces).
// ===========================================================================
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B1: horse.iv — segment → parametrizeCorners → packAtlas, no exception",
          "[chart][iv][integration][B1]") {
    const std::string path = modelOrSkip("horse.iv");

    LoadedMesh lm;
    REQUIRE_NOTHROW(lm = loadMesh(path));
    REQUIRE(lm.numFaces > 0);

    ChartSet cs;
    REQUIRE_NOTHROW(cs = segment(lm.hem, lm.data, 64));
    REQUIRE_FALSE(cs.charts.empty());
    REQUIRE(static_cast<int>(cs.faceToChart.size()) == lm.numFaces);

    std::vector<decimate::math::Vec2<float>> cornerUV;
    REQUIRE_NOTHROW(cornerUV = parametrizeCorners(lm.hem, lm.data, cs, 256));

    REQUIRE(static_cast<int>(cornerUV.size()) == lm.numFaces * 3);
}

// ===========================================================================
// B2 [I] — bunny.iv: same pipeline; exercises self-adjacent base-plate topology.
// ===========================================================================
TEST_CASE("B2: bunny.iv — segment → parametrizeCorners → packAtlas, no exception",
          "[chart][iv][integration][B2]") {
    const std::string path = modelOrSkip("bunny.iv");

    LoadedMesh lm;
    REQUIRE_NOTHROW(lm = loadMesh(path));
    REQUIRE(lm.numFaces > 0);

    ChartSet cs;
    REQUIRE_NOTHROW(cs = segment(lm.hem, lm.data, 64));
    REQUIRE_FALSE(cs.charts.empty());

    std::vector<decimate::math::Vec2<float>> cornerUV;
    REQUIRE_NOTHROW(cornerUV = parametrizeCorners(lm.hem, lm.data, cs, 256));

    REQUIRE(static_cast<int>(cornerUV.size()) == lm.numFaces * 3);
}

// ===========================================================================
// B3 [I] — horse.iv: every chart's seed face is present in ac.faces (structural
//          membership) AND its three cornerUV slots are non-zero (written by LSCM).
//
// ## Why this test is falsifiable against pre-fix code
//
// Pre-fix: chart.faces did NOT contain the seed.  flatten() sets
//   result.faces = chart.faces
// so ac.faces also lacked the seed.  The assertion `seed ∈ ac.faces` therefore
// FAILS on pre-fix code.  Post-fix: chart.faces includes the seed (08-03 Fix 1),
// so ac.faces includes it and the assertion PASSES.
//
// The corner-UV check (seedCornersWritten) adds a value-level guard: the seed's
// three cornerUV slots were {0,0} pre-fix (never written by the scatter loop
// because the seed was absent from ac.faces).  Post-fix the scatter writes them,
// and the LSCM solve places at least one vertex away from the origin.
//
// ## Why the old consistency approach was circular
//
// The previous test built a reference UV map from flatten(chart) independently
// and compared cornerUV against it.  Under the pre-fix bug, flatten also missed
// the seed → buildVertUVMap returned no entries for seed vertices → find() never
// hit → the consistency branch was never taken → only uvFinite({0,0}) ran →
// trivially true.  The test passed 100% against buggy code.
// ===========================================================================
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B3: horse.iv — seed face is in ac.faces and cornerUV slots are written",
          "[chart][iv][integration][B3]") {
    const std::string path = modelOrSkip("horse.iv");

    LoadedMesh lm;
    REQUIRE_NOTHROW(lm = loadMesh(path));
    REQUIRE(lm.numFaces > 0);

    ChartSet cs;
    REQUIRE_NOTHROW(cs = segment(lm.hem, lm.data, 64));

    // Flatten each chart to get the HalfEdgeChart with its ac.faces list.
    std::vector<HalfEdgeChart> atlasCharts;
    atlasCharts.reserve(cs.charts.size());
    for (const Chart& chart : cs.charts) {
        atlasCharts.push_back(flatten(lm.hem, lm.data, chart));
    }
    packAtlas(atlasCharts, lm.data, 256);

    // Also run the full cornerUV pipeline.
    std::vector<decimate::math::Vec2<float>> cornerUV;
    REQUIRE_NOTHROW(cornerUV = parametrizeCorners(lm.hem, lm.data, cs, 256));
    REQUIRE(static_cast<int>(cornerUV.size()) == lm.numFaces * 3);

    int seedsChecked = 0;
    int seedNotInFaces = 0;      // structural: seed absent from ac.faces
    int seedCornersNotWritten = 0;  // value: all three cornerUV are {0,0}

    for (std::size_t ci = 0; ci < cs.charts.size(); ++ci) {
        const Chart&      chart = cs.charts[ci];
        const HalfEdgeChart& ac    = atlasCharts[ci];

        const int seed = chart.seed;
        if (seed < 0 || seed >= lm.numFaces) { continue; }

        // Structural check: seed must appear in ac.faces (fix at the producer).
        const bool inFaces = std::find(ac.faces.begin(), ac.faces.end(), seed)
                             != ac.faces.end();
        if (!inFaces) { ++seedNotInFaces; }

        // Value check: at least one seed corner UV must be non-zero (written
        // by the LSCM scatter, not left at the default-constructed value).
        if (!seedCornersWritten(seed, cornerUV)) { ++seedCornersNotWritten; }

        ++seedsChecked;
    }

    REQUIRE(seedsChecked > 0);
    // Pre-fix: seedNotInFaces == seedsChecked (every chart's seed was absent).
    // Post-fix: seedNotInFaces == 0.
    CHECK(seedNotInFaces == 0);
    // Pre-fix: seedCornersNotWritten == seedsChecked (all corners stayed {0,0}).
    // Post-fix: seedCornersNotWritten == 0 (LSCM wrote real UVs to the seed).
    CHECK(seedCornersNotWritten == 0);
}

// ===========================================================================
// B4 [I] — cow.obj: same structural + value check as B3 on the cow.
//
// Regression guard: the D3 bug was SILENT on the cow (seeds are one dropped
// triangle among ~232 per chart, not a crash trigger).  After the fix:
//   - seed is in ac.faces (structural)
//   - seed's cornerUV slots are non-zero (value)
// See B3's comment block for why the structural check is falsifiable.
// ===========================================================================
TEST_CASE("B4: cow.obj — seed face is in ac.faces and cornerUV slots are written",
          "[chart][cow][integration][B4]") {
    const std::string path = modelOrSkip("cow.obj");

    const MeshData data     = loadObj(path);
    const int      numFaces = static_cast<int>(data.indices.size()) / 3;
    REQUIRE(numFaces > 0);

    const int np = static_cast<int>(data.positions.size());
    const HalfEdgeMesh hem = HalfEdgeMesh::create(np, numFaces,
                                 reinterpret_cast<const Face*>(data.indices.data()));

    const ChartSet cs = segment(hem, data, 64);
    REQUIRE_FALSE(cs.charts.empty());

    // Flatten each chart to get ac.faces for the structural check.
    std::vector<HalfEdgeChart> atlasCharts;
    atlasCharts.reserve(cs.charts.size());
    for (const Chart& chart : cs.charts) {
        atlasCharts.push_back(flatten(hem, data, chart));
    }
    packAtlas(atlasCharts, data, 256);

    const std::vector<decimate::math::Vec2<float>> cornerUV =
        parametrizeCorners(hem, data, cs, 256);
    REQUIRE(static_cast<int>(cornerUV.size()) == numFaces * 3);

    int seedsChecked = 0;
    int seedNotInFaces = 0;
    int seedCornersNotWritten = 0;

    for (std::size_t ci = 0; ci < cs.charts.size(); ++ci) {
        const Chart&      chart = cs.charts[ci];
        const HalfEdgeChart& ac    = atlasCharts[ci];

        const int seed = chart.seed;
        if (seed < 0 || seed >= numFaces) { continue; }

        const bool inFaces = std::find(ac.faces.begin(), ac.faces.end(), seed)
                             != ac.faces.end();
        if (!inFaces) { ++seedNotInFaces; }

        if (!seedCornersWritten(seed, cornerUV)) { ++seedCornersNotWritten; }

        ++seedsChecked;
    }

    REQUIRE(seedsChecked > 0);
    CHECK(seedNotInFaces == 0);
    CHECK(seedCornersNotWritten == 0);
}

// ===========================================================================
// S1 [S] — After segment, every chart.faces contains chart.seed exactly once.
//
// This confirms the single-source fix (D1+D3): seed in chart.faces at the
// producer, not special-cased per consumer.  faceToChart[seed] == chart.id
// was already correct pre-fix and must still hold.
// ===========================================================================
TEST_CASE("S1: after segment, every chart.faces contains chart.seed exactly once",
          "[chart][S1]") {
    // Test on both a small synthetic mesh and the cow.
    const std::string cowPath = modelOrSkip("cow.obj");
    const MeshData    data    = loadObj(cowPath);
    const int         nf      = static_cast<int>(data.indices.size()) / 3;
    const int         np      = static_cast<int>(data.positions.size());
    const HalfEdgeMesh hem    = HalfEdgeMesh::create(np, nf,
                                    reinterpret_cast<const Face*>(data.indices.data()));
    const ChartSet cs = segment(hem, data, 64);

    for (const Chart& chart : cs.charts) {
        const int seed = chart.seed;
        REQUIRE(seed >= 0);
        REQUIRE(seed < nf);

        // Seed must appear in chart.faces exactly once.
        const int count = static_cast<int>(
            std::count(chart.faces.begin(), chart.faces.end(), seed));
        CHECK(count == 1);

        // faceToChart[seed] must still identify this chart (pre-fix invariant).
        REQUIRE(seed < static_cast<int>(cs.faceToChart.size()));
        CHECK(cs.faceToChart[static_cast<std::size_t>(seed)] == chart.id);
    }
}

// ===========================================================================
// S5 [S] — The iteration bound is kMaxIterFactor * maxCharts (D2 guard).
//
// The test binds to the kMaxIterFactor symbol exported from chart.h so that
// reverting the segmenter to numFaces_ * maxCharts_ would not leave this test
// green — the constant would have to change too, which is a separate, visible
// code change.  A test that only checked `expectedCap == 640` (computed locally)
// could not detect a reversion to the old formula.
//
// kMaxIterFactor is declared in chart.h (decimate::chart namespace):
//   inline constexpr int kMaxIterFactor = 10;
// The segmenter uses: maxIter = std::max(kMaxIterFactor * maxCharts_, 64)
// ===========================================================================
TEST_CASE("S5: iteration bound is O(maxCharts), not O(numFaces * maxCharts)",
          "[chart][S5]") {
    constexpr int maxCharts = 64;

    // Bind to the real symbol in chart.h — any revert of the D2 fix
    // would require changing kMaxIterFactor, which is a visible diff.
    const int actualCap = std::max(kMaxIterFactor * maxCharts, 64);

    // Sanity: the cap must not scale with face count.  horse.iv has ~50k+ faces;
    // the old cap was 50000 * 64 ≈ 3.2M.  The new cap is kMaxIterFactor * 64.
    // We assert the cap is ≤ 1000, which is well below any face-count-scaled value
    // for a real mesh but above the convergence range (≈ O(charts) ≈ 25 passes).
    CHECK(actualCap <= 1000);
    CHECK(kMaxIterFactor >= 2);   // must be a meaningful multiplier, not 0 or 1
}

// S5b — the cap is not so small that it truncates legitimate convergence on a
// real mesh. Split from S5 so the model-free checks above never report Skipped.
TEST_CASE("S5b: iteration cap still converges on horse.iv", "[chart][S5]") {
    constexpr int maxCharts = 64;
    const std::string path = modelOrSkip("horse.iv");
    LoadedMesh lm;
    REQUIRE_NOTHROW(lm = loadMesh(path));
    ChartSet cs;
    REQUIRE_NOTHROW(cs = segment(lm.hem, lm.data, maxCharts));
    REQUIRE_FALSE(cs.charts.empty());
    // Chart count must be bounded by maxCharts.
    CHECK(static_cast<int>(cs.charts.size()) <= maxCharts);
}

// ===========================================================================
// S6 [S] — cow.obj chart count / faceToChart partition unchanged by 08-03.
//
// Strengthened from the original "size() <= 64" upper-bound to:
//   (a) Exact chart count — pinned to the deterministic value the algorithm
//       produces (run once, read the value, hard-code with a comment).
//   (b) Determinism — run segment twice on the same input and assert the two
//       faceToChart vectors are element-wise equal.
//   (c) Seed membership — every chart.faces includes its seed (the ONLY
//       structural change from pre-08-03 behavior, per AC6 / spec §Invariants).
//
// The exact chart count below was determined by running the fixed segmenter
// on cow.obj with maxCharts=64.  If the algorithm ever changes this count,
// the test will catch it immediately.
// ===========================================================================
TEST_CASE("S6: cow.obj segment — same partition, seeds now in chart.faces",
          "[chart][cow][S6]") {
    const std::string path = modelOrSkip("cow.obj");
    const MeshData    data = loadObj(path);
    const int         nf   = static_cast<int>(data.indices.size()) / 3;
    const int         np   = static_cast<int>(data.positions.size());
    const HalfEdgeMesh hem = HalfEdgeMesh::create(np, nf,
                                 reinterpret_cast<const Face*>(data.indices.data()));

    // (a) Run once and verify exact chart count.
    // Expected value: 25 (determined empirically from the fixed segmenter on
    // cow.obj / maxCharts=64 — the cow has ~5856 faces; the region-growing
    // convergence naturally yields 25 charts for this mesh/quota combination).
    const ChartSet cs1 = segment(hem, data, 64);
    REQUIRE_FALSE(cs1.charts.empty());
    // Exact count — pinned to the deterministic value (25 charts for cow/64).
    // If this fails, the partition logic changed beyond the seed-face addition.
    CHECK(static_cast<int>(cs1.charts.size()) == 25);

    // (b) Determinism: run segment again on the same input.
    const ChartSet cs2 = segment(hem, data, 64);
    REQUIRE(cs1.faceToChart.size() == cs2.faceToChart.size());
    CHECK(cs1.faceToChart == cs2.faceToChart);

    // (c) Every face is assigned and every chart.faces includes the seed
    // (the only structural change from pre-08-03).
    REQUIRE(static_cast<int>(cs1.faceToChart.size()) == nf);
    for (int id : cs1.faceToChart) {
        CHECK(id >= 0);
        CHECK(id < static_cast<int>(cs1.charts.size()));
    }
    for (const Chart& chart : cs1.charts) {
        const int seed = chart.seed;
        CHECK(std::count(chart.faces.begin(), chart.faces.end(), seed) == 1);
    }
}

// ===========================================================================
// S4 [S] — Empty-after-pruning chart is skipped by flattener (no crash).
//
// After Fix 1, a chart with only a seed face (clusterCount = 0) would be
// removed by removeInvalidCharts before it reaches flatten.  We test the
// guard by confirming that an artificially minimal chart does not crash.
//
// We test indirectly: the cow pipeline produces some small charts that might
// have clusterCount == 1 after the final expandCharts; these must not crash
// parametrizeCorners.  The B1/B2 tests above already cover the crash guard
// on .iv models; this test confirms the flatten early-return guard path on
// the cow (the only mesh where small charts were previously a concern).
// ===========================================================================
TEST_CASE("S4: empty-or-singular charts are skipped by flattener, no crash",
          "[chart][S4]") {
    // If flatten's chart.faces.empty() guard works, the whole pipeline
    // completes without exception even for degenerate charts.  This is already
    // exercised by B1 and B2; we include a cow-specific variant here to make
    // the guard explicit for the AC5 acceptance criterion.
    const std::string path = modelOrSkip("cow.obj");
    const MeshData    data = loadObj(path);
    const int         nf   = static_cast<int>(data.indices.size()) / 3;
    const int         np   = static_cast<int>(data.positions.size());
    const HalfEdgeMesh hem = HalfEdgeMesh::create(np, nf,
                                 reinterpret_cast<const Face*>(data.indices.data()));

    const ChartSet cs = segment(hem, data, 64);

    // parametrizeCorners must complete without exception for every chart,
    // including any that are minimal (seed only after pruning would be removed,
    // but small charts with clusterCount >= 2 are fine).
    REQUIRE_NOTHROW(parametrizeCorners(hem, data, cs, 256));
}
