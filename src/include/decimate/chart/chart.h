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
///   - Seeding: one chart at the lowest face of every connected component (face 0 on a
///     single-component mesh), then one new chart per iteration at the last face grown
///     (`lastFace_`), up to `maxCharts` (no curvature term anywhere)
///   - Lazy-deletion min-heap (no decrease-key — Beat-2b forbidden path)
///   - Lloyd seed relocation via deepest-face BFS; cycle detection via seedHistory
///   - Boundary straightening at threshold 0.3 (`kStraightenThreshold`, overridable per call)
///   - Chart normal = sum of RAW face normals (cross products, so area-weighted by
///     construction, `legacy2003/math/algorithm.h:55-67`), normalized
///   - A component left without any chart becomes one chart covering the whole component
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

/// @brief A single chart: a contiguous set of faces with a shared normal proxy.
///
/// ## Normal accumulation (Beat 2a)
///
/// Every time a face is grown into the chart, its raw (unnormalised cross-product, hence
/// area-weighted) face normal is added to `normalAccum` (the sum) and `normal` is the
/// running normalised sum.
/// This matches `CChart::Grow()` / `CChart::Normal()` exactly.
struct Chart {
    int              id = -1;          ///< Chart index in the parent ChartSet.
    int              seed = -1;        ///< Seed face index for this chart.
    int              clusterCount = 0; ///< Number of Grow() calls (mirrors 2003 `CChart::clusterCount`, `legacy2003/chartification/cchart.h:29`).
    std::vector<int> faces;            ///< Complete triangle set for this chart (grown faces + seed).
                                       ///  The seed face is always present exactly once (08-03 contract).
    math::Vec3f   normal;              ///< Normalised sum of raw face normals (area-weighted).
    math::Vec3f   normalAccum;         ///< Running unnormalized normal accumulator.
    std::set<int> seedHistory;         ///< Cycle-detection set for Lloyd relocation.
    int           searchIndex = -1;    ///< Deepest-face BFS scratch (ComputeSeeds).

    /// Intra-chart seam half-edges (Spec 19-02): every half-edge `h` whose `sym` exists, whose
    /// two faces both belong to this chart, and which the cutter left flagged as boundary after
    /// straightening, i.e. the edges 2003 `CutBoundaryPath` cuts open inside a chart
    /// (`legacy2003/chartification/cchartification.cpp:488-500`).  Both halves of an edge are
    /// listed, ascending.  `flatten()` does not merge corners across them, so a chart that
    /// touches itself is cut into a disk.  Inter-chart edges are never listed.
    std::vector<int> seamHalfEdges;

    /// @brief Grow this chart by one face: accumulate normal, record face.
    ///
    /// Mirrors `CChart::Grow(Vec3& normal)` exactly.
    void grow(const math::Vec3f &faceNormal);

    /// @brief Update the seed and return true if a cycle was detected.
    ///
    /// Mirrors `CChart::SetSeed(int seed, Vec3& normal)`:
    /// sets seed + resets normalAccum to faceNormal; returns true if
    /// `seed` is already in seedHistory (cycle = convergence).
    bool setSeed(int newSeed, const math::Vec3f &faceNormal);

    /// @brief Reset face list and accumulators, keeping the seed and history.
    ///
    /// Mirrors `CChart::Clear()` — called at the start of each `ExpandCharts`.
    void clear();
};


/// @brief Which face "midpoint" the cutter uses for centroid distances (Spec 19-01, D-a).
///
/// The 2003 `middle()` returns `((a + b) * 0.5f + c) * 0.5f` = a/4 + b/4 + c/2
/// (`legacy2003/math/algorithm.h:70-73`), which is not the centroid: a 2003 bug.
enum class CentroidMode {
    TrueCentroid,    ///< `(a + b + c) / 3`.  The default for production and all tests.
    Legacy2003Middle ///< Reproduces 2003 `middle()` bit-for-bit.  Only the 2003 golden uses it.
};


/// @brief Output of `segment()`: per-face chart assignment + chart metadata.
struct ChartSet {
    std::vector<Chart> charts;      ///< One entry per chart.
    std::vector<int>   faceToChart; ///< Size == mesh.faceCount(); every entry in [0, charts.size()).
    /// Face-centre formula `segment()` used; `splitChart` reuses it so a split costs faces
    /// exactly as the original cut did (Spec 19-03).
    CentroidMode centroidMode = CentroidMode::TrueCentroid;
};


/// @brief Split chart @p chartId in two, in place (Spec 19-03).
///
/// A restricted two-seed region grow over the faces of the chart only: the seeds are the
/// chart's `seed` and @p foldFace, the cost is the segmenter's edge cost (same heap, same
/// running chart normals), no Lloyd relocation, no straightening, no other chart touched.
/// The piece grown from @p foldFace becomes a new chart (`id == charts.size()` before the
/// call, `seed == foldFace`); `faceToChart` and both `faces` lists (ascending) are rewritten.
/// Both pieces stay dual-connected.  This is the ONLY way a `ChartSet` gains charts after
/// `segment()`.
///
/// @throws std::invalid_argument if @p chartId is out of range, @p foldFace is not a face of
///         the chart, or @p foldFace equals the chart's seed.
void splitChart(const decimate::HalfEdgeMesh &mesh,
                const decimate::MeshData     &data,
                ChartSet                     &charts,
                int                           chartId,
                int                           foldFace,
                CentroidMode                  centroidMode);


/// @brief Default boundary-straighten threshold of `segment()`.
///
/// Source: `if (err < 0.3)` at `legacy2003/chartification/cchartification.cpp:542`.
inline constexpr double kStraightenThreshold = 0.3;


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
/// 1. **Seed select** — one chart at the lowest face of each connected component,
///    then one chart per iteration at the last face grown, capped at `maxCharts`
///    (the cap can be exceeded only when there are more components than `maxCharts`).
///    A newborn chart's seed is NOT in its `seedHistory`
///    (`legacy2003/chartification/cchart.cpp:23-31`), so the cap binds.
/// 2. **Region grow** — simultaneous multi-source BFS keyed by
///    `EdgeCost = (1 - cos_angle) * centroidDistance`, lazy-deletion heap
///    (stale entries skipped when the opposite face is already claimed).
/// 3. **Seed relocate** — Lloyd BFS to find deepest face per chart; accept
///    when `seedHistory` detects a cycle.
/// 4. **Invalid-chart cleanup** — charts with fewer than 2 faces are removed
///    and `ExpandCharts` re-run.
/// 5. **Boundary straighten** — faces with two boundary edges on the same
///    opposite chart are reassigned if normalError < `straightenThreshold`
///    (default `kStraightenThreshold` = 0.3).
/// 6. **Component fallback** — a component that lost every chart (all of its charts
///    had fewer than 2 grown faces) becomes one chart owning the whole component.
///
/// @param mesh       Half-edge connectivity (from `HalfEdgeMesh::create`).
/// @param data       Vertex positions (one Vec3f per vertex).
/// @param maxCharts  Maximum number of charts to produce (default 64).
/// @param straightenThreshold  Normal-error bound (`1 - cos`, in `[0, 2]`) below which the
///                   boundary-straighten pass (stage 5) moves a face to the chart across its
///                   two boundary edges.  Default `kStraightenThreshold` (0.3, the constant
///                   at `legacy2003/chartification/cchartification.cpp:542`) reproduces the
///                   pre-Spec-17 port byte-for-byte (fingerprints, `seam_straighten_test.cpp` B1);
///                   the scan order is the 2003 one (`i*3 + (j+1)%3`, `legacy2003/chartification/cchartification.cpp:519-521`),
///                   and a moved face also erases the stale edge (see `segmenter.cpp`).  `0.0` disables the
///                   pass (the test is strict `<`), `2.0` accepts every eligible face.
/// @param centroidMode  Face-centre formula feeding the edge cost; see `CentroidMode`.
///                   Default `TrueCentroid`.  `Legacy2003Middle` exists for the 2003 golden only.
/// @returns          `ChartSet` with `faceToChart.size() == mesh.faceCount()`,
///                   every entry in `[0, charts.size())`, no entry == -1.
///
/// @pre  `mesh.faceCount() > 0`.
/// @pre  `data.positions.size() >= mesh.pointCount()`.
[[nodiscard]] ChartSet segment(const decimate::HalfEdgeMesh &mesh,
                               const decimate::MeshData     &data,
                               int                           maxCharts = 64,
                               double                        straightenThreshold = kStraightenThreshold,
                               CentroidMode                  centroidMode = CentroidMode::TrueCentroid);

} // namespace decimate::chart
