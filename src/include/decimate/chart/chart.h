// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file chart.h
/// @brief Region-growing chart segmentation — public API (Spec 08-01).
///
/// @note The public signatures that take `HalfEdgeMesh` / expose `HalfEdgeChart` are
///       subject to change in 0.7 (flat `MeshData` overloads; half-edge becomes an
///       implementation detail).
///
/// Port of `CChartification` / `CChart` / `CChartSet` / `CChartBoundary`
/// from `legacy/chartification/` into dependency-free C++20.
///
/// ## Design contract (Beat 2a — behavior-preserving)
///
/// The algorithm is a simultaneous multi-source region grow keyed by
/// `EdgeCost = (1 - cos(angle)) * centroidDistance`, with:
///   - Highest-curvature seed selection (curvature = 1 - minDot over 1-ring)
///   - Lazy-deletion min-heap (no decrease-key — Beat-2b forbidden path)
///   - Lloyd seed relocation via deepest-face BFS; cycle detection via seedHistory
///   - Boundary straightening at threshold 0.3
///   - Chart normal = plain sum of face normals, normalized (not area-weighted)
///   - Orphan reassignment: any unassigned face takes lowest-id neighbour's chart
///
/// ## Boundary contracts
///
/// - No `#include <GL/...>`, no GLUT, no `windows.h`, no global mutable state.
/// - No `#include` reaching into `legacy/`.
/// - No file I/O — operates on in-memory `HalfEdgeMesh` / `MeshData` only.
/// - `legacy/split/` is untouched (Spec 09 territory).

#include <set>
#include <vector>

#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/math/vec3.h>

namespace decimate::chart {

// ---------------------------------------------------------------------------
// Chart — one segmented region of the mesh
// ---------------------------------------------------------------------------

/// @brief A single chart: a contiguous set of faces with a shared normal proxy.
///
/// ## Normal accumulation (Beat 2a — NOT area-weighted)
///
/// Every time a face is grown into the chart, its raw face normal is added
/// to `normalAccum` (the sum) and `normal` is the running normalised sum.
/// This matches `CChart::Grow()` / `CChart::Normal()` exactly.
struct Chart {
    int                  id           = -1;  ///< Chart index in the parent ChartSet.
    int                  seed         = -1;  ///< Seed face index for this chart.
    int                  clusterCount = 0;   ///< Number of Grow() calls (mirrors legacy CChart::clusterCount).
    std::vector<int>     faces;              ///< Complete triangle set for this chart (grown faces + seed).
                                             ///  The seed face is always present exactly once (08-03 contract).
    decimate::math::Vec3<float>    normal;             ///< Normalised sum of face normals (plain sum, Beat 2a).
    decimate::math::Vec3<float>    normalAccum;        ///< Running unnormalized normal accumulator.
    std::set<int>        seedHistory;        ///< Cycle-detection set for Lloyd relocation.
    int                  searchIndex = -1;   ///< Deepest-face BFS scratch (ComputeSeeds).

    /// @brief Grow this chart by one face: accumulate normal, record face.
    ///
    /// Mirrors `CChart::Grow(Vec3& normal)` exactly.
    void grow(const decimate::math::Vec3<float>& faceNormal);

    /// @brief Update the seed and return true if a cycle was detected.
    ///
    /// Mirrors `CChart::SetSeed(int seed, Vec3& normal)`:
    /// sets seed + resets normalAccum to faceNormal; returns true if
    /// `seed` is already in seedHistory (cycle = convergence).
    bool setSeed(int newSeed, const decimate::math::Vec3<float>& faceNormal);

    /// @brief Reset face list and accumulators, keeping the seed and history.
    ///
    /// Mirrors `CChart::Clear()` — called at the start of each `ExpandCharts`.
    void clear();
};

// ---------------------------------------------------------------------------
// ChartSet — collection of Chart objects produced by segment()
// ---------------------------------------------------------------------------

/// @brief Output of `segment()`: per-face chart assignment + chart metadata.
struct ChartSet {
    std::vector<Chart> charts;          ///< One entry per chart.
    std::vector<int>   faceToChart;     ///< Size == mesh.faceCount(); every entry in [0, charts.size()).
};

// ---------------------------------------------------------------------------
// segment — public entry point
// ---------------------------------------------------------------------------

/// @brief Safety-cap multiplier for the segment driver loop (08-03 Fix 3 / D2).
///
/// `maxIter = std::max(kMaxIterFactor * maxCharts, 64)`.  The natural exit
/// is `computeSeeds()` convergence (≈ O(charts) passes); this cap is a
/// backstop only.  Tests bind to this constant to ensure the cap cannot
/// silently revert to the O(numFaces × maxCharts) form that caused D2.
inline constexpr int kMaxIterFactor = 10;

/// @brief Partition @p mesh into at most @p maxCharts charts using the
///        region-growing algorithm from the 2003 `CChartification`.
///
/// ## Algorithm stages (behavior-preserving port — Beat 2a)
///
/// 1. **Seed select** — highest-curvature faces (curvature = 1 - minDot over
///    1-ring), capped at `maxCharts`.
/// 2. **Region grow** — simultaneous multi-source BFS keyed by
///    `EdgeCost = (1 - cos_angle) * centroidDistance`, lazy-deletion heap
///    (stale entries skipped when the opposite face is already claimed).
/// 3. **Seed relocate** — Lloyd BFS to find deepest face per chart; accept
///    when `seedHistory` detects a cycle.
/// 4. **Invalid-chart cleanup** — charts with fewer than 2 faces are removed
///    and `ExpandCharts` re-run.
/// 5. **Boundary straighten** — faces with two boundary edges on the same
///    opposite chart are reassigned if normalError < 0.3.
/// 6. **Orphan assignment** — any face still unassigned is given its
///    lowest-id assigned neighbour's chart id.
///
/// @param mesh       Half-edge connectivity (from `HalfEdgeMesh::create`).
/// @param data       Vertex positions (one Vec3f per vertex).
/// @param maxCharts  Maximum number of charts to produce (default 64).
/// @returns          `ChartSet` with `faceToChart.size() == mesh.faceCount()`,
///                   every entry in `[0, charts.size())`, no entry == -1.
///
/// @pre  `mesh.faceCount() > 0`.
/// @pre  `data.positions.size() >= mesh.pointCount()`.
[[nodiscard]] ChartSet segment(const decimate::HalfEdgeMesh& mesh,
                                const decimate::MeshData&     data,
                                int                        maxCharts = 64);

}  // namespace decimate::chart
