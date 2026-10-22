// This file is part of the decimate library; see decimate.h for version/license details.
/// @file segmenter.cpp
/// @brief Region-growing chart segmentation — port of `CChartification`.
///
/// ## Behavior-preservation notes (Beat 2a)
///
/// ### Flagged deviations from the 2003 cutter (Spec 19-01)
///
///  - Centroid: the DEFAULT is the true centroid `(a+b+c)/3`, a deliberate fix of the
///    2003 `middle()` bug, which returns a/4 + b/4 + c/2 (`legacy2003/math/algorithm.h:70-73`,
///    used at `legacy2003/chartification/cchartification.cpp:53`).  `CentroidMode::
///    Legacy2003Middle` reproduces the 2003 formula; only the `chart_legacy_golden` test uses it.
///  - Seeding: one seed per connected component (ascending lowest face) before the loop, and a
///    per-component fallback chart after the removal loop.  2003 seeds only face 0 and skips
///    unreached faces (`legacy2003/chartification/cchartification.cpp:437`).  Identical on single-component meshes.
///  - Straighten stale edge (Spec 19-02): after a face moves, 2003 calls `SetBoundaryFlag(h->next->next)`
///    and `Erase(h)` (`legacy2003/chartification/cchartification.cpp:547-548`) but leaves `h->next` flagged.  That edge is
///    now between two faces of the same chart, so it stays behind as a spurious self-seam.  The
///    port also erases `next(h)`.  Measured effect on the 2003 golden (`chart_legacy_golden`,
///    Legacy2003Middle): the partition is unchanged; summed chart-local vertex counts are lower by
///    cow/24 -1, cow/32 -1, cow/64 -1, horse/128 -4, bunny/8 -21, bunny/32 -18 (the six cases pinned as `localDelta`).  Scan order is the 2003 one (`i*3 + (j+1)%3`,
///    `legacy2003/chartification/cchartification.cpp:519-521`).
///  - `Chart::seamHalfEdges` (Spec 19-02): the intra-chart edges still flagged after straightening
///    are exported (2003 `CutBoundaryPath` cuts them inside `CreateVertexStructure`, `:488-500`).
///  - `lastFace_` guard: no new chart when `lastFace_` is already a seed (2003 has none: `legacy2003/chartification/cchartification.cpp:92-95`).
///  - D-d: 2003 never initialises `CChart::searchIndex` (`legacy2003/chartification/cchart.h:42`,
///    ctor `legacy2003/chartification/cchart.cpp:23-31`); on a closed mesh whose single chart has no boundary left,
///    `ApplySeeds` reads garbage.  The port keeps `searchIndex = -1` (seed unchanged).  Not
///    emulated; spot.obj is therefore excluded from the golden.
///  - Edge-cost denominator is float (product and sqrt), as `legacy2003/chartification/cchartification.cpp:131`;
///    a zero-length-normal guard (returns 0) is a port addition, unreachable on the goldens.
///  - No final regrow after the Lloyd loop unless `RemoveInvalidCharts` removed a chart
///    (`legacy2003/chartification/cchartification.cpp:97-105`): the chart map is the grow from the pre-move seeds.
///
/// ### Legacy Chartify loop (legacy2003/chartification/cchartification.cpp:74-114)
///
/// ```
/// CChart *chart = charts->Create(0);   // seed at face 0
/// do {
///   ExpandCharts(charts);
///   if (charts->Size() < numDesiredCharts)
///     chart = charts->Create(lastface); // add one per iteration at lastface
/// } while (!ComputeSeeds(charts, boundary));
/// while (charts->RemoveInvalidCharts()) ExpandCharts(charts);
/// while (Straighten(boundary));
/// ```
///
/// ### EdgeCost (legacy2003/chartification/cchartification.cpp:120-135)
///
/// Legacy:
/// ```cpp
/// Vec3 nrm1 = faceNormals[h->face];
/// Vec3 nrm2 = faceNormals[h->sym->face];
/// if (chartMap[h->face]) nrm1 = chartMap[h->face]->Normal();
/// else if (h->sym && chartMap[h->sym->face]) nrm2 = chartMap[h->sym->face]->Normal();
/// double nrmFactor = (nrm1 * nrm2) / sqrt(nrm1.length2() * nrm2.length2());
/// return (1.0 - nrmFactor) * distance;
/// ```
/// Operator mapping:
///   legacy `nrm1 * nrm2` (Vec3 dot)  → `nrm1.dot(nrm2)`
///   legacy `nrm1.length2()`           → `nrm1.lengthSquared()`
///   legacy `(a - b).length()`         → `(a - b).length()`
///
/// ### Grow-loop pop/claim logic (legacy2003/chartification/cchartification.cpp:310-345)
///
/// ```
/// while (heap.pop(&h, &value)) {
///   int f    = h->face;
///   int fopp = h->sym->face;
///   if (chartMap[fopp] == NULL) {       // lazy skip if already claimed
///     chartMap[f]->Grow(faceNormals[fopp]);
///     chartMap[fopp] = chartMap[f];
///     boundary->Erase(h);
///     for i in {0,1}: push h->sym->next+i if sym && IsBoundary(loop)
///     // RepairHeap commented out — stale entries stay (lazy deletion)
///   }
/// }
/// ```
///
/// ### Iteration order preservation
///
/// 2003 `CChartSet` uses `std::set<CChart*>` — pointer-ordered (`legacy2003/chartification/cchartset.h:24`).
/// We use `std::vector<unique_ptr<Chart>>` with sequential access matching
/// chart creation order (face 0 first, lastface after each iteration).
/// Creation order is deterministic for a fixed input — output is identical
/// across runs (T-S1).
///
/// 2003 `CSet<int> seedHistory` (`legacy2003/chartification/cchart.h:28`) is sorted — mapped to `std::set<int>`.
///
/// ### ClusterCount vs faces
///
/// 2003 `CChart::clusterCount` is zeroed at `legacy2003/chartification/cchart.cpp:28` and incremented at `legacy2003/chartification/cchart.cpp:157`.
/// 2003 `RemoveInvalidCharts` removes charts with `clusterCount < 2` (`legacy2003/chartification/cchartset.cpp:78`).
/// A chart starts with just a seed (clusterCount=0); each face grown in
/// increments clusterCount by 1.  We track this as `Chart::clusterCount`.

#include <decimate/chart/chart.h>

#include "segmenter_detail.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <queue>
#include <vector>

#include <decimate/math/vec3.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>

namespace decimate::chart {

namespace {

/// Lazy-deletion min-heap element.  A `std::priority_queue` with `std::greater` is a min-heap, like 2003 `Heap<Halfedge*, double>` (pop returns the minimum cost).
struct HeapEntry {
    double cost;
    int    heIndex; // half-edge index in HalfEdgeMesh

    bool operator>(const HeapEntry &o) const noexcept {
        return cost > o.cost;
    }
};


using MinHeap = std::priority_queue<HeapEntry,
                                    std::vector<HeapEntry>,
                                    std::greater<HeapEntry>>;

/// Pre-computed per-face normals and centroids.
struct FaceGeometry {
    math::Vec3f normal;   ///< Raw (unnormalised) face normal.
    math::Vec3f centroid; ///< Face centroid.
};


/// Compute raw face normal and centroid for each face.
/// Raw normal = (b-a) × (c-a) — unnormalised (legacy `normal_raw`, `legacy2003/math/algorithm.h:55-67`).
/// Centroid: the ONLY centroid formula in the cutter.  `TrueCentroid` = (a+b+c)/3 (a deliberate
/// fix); `Legacy2003Middle` = ((a+b)*0.5f + c)*0.5f, the 2003 `middle()` bug
/// (`legacy2003/math/algorithm.h:70-73`).
std::vector<FaceGeometry> computeFaceGeometry(
    const decimate::HalfEdgeMesh &mesh,
    const decimate::MeshData     &data,
    CentroidMode                  mode) {
    const int                 numFaces = static_cast<int>(mesh.faceCount());
    std::vector<FaceGeometry> geom(static_cast<std::size_t>(numFaces));

    for (int f = 0; f < numFaces; ++f) {
        // Legacy WriteFace wiring (half_edge_mesh.cpp):
        //   he[f*3+0]: end=c, next=f*3+1
        //   he[f*3+1]: end=b, next=f*3+2
        //   he[f*3+2]: end=a, next=f*3+0
        // So: a=he[f*3+2].end, b=he[f*3+1].end, c=he[f*3+0].end.
        const auto &he0 = mesh.halfedge(static_cast<decimate::HeIndex>(f * 3 + 0));
        const auto &he1 = mesh.halfedge(static_cast<decimate::HeIndex>(f * 3 + 1));
        const auto &he2 = mesh.halfedge(static_cast<decimate::HeIndex>(f * 3 + 2));

        auto getPos = [&](std::uint32_t vi) -> math::Vec3f {
            return data.positions[static_cast<std::size_t>(vi)];
        };

        const math::Vec3f va = getPos(he2.end);
        const math::Vec3f vb = getPos(he1.end);
        const math::Vec3f vc = getPos(he0.end);

        geom[static_cast<std::size_t>(f)].normal = (vb - va).cross(vc - va);
        geom[static_cast<std::size_t>(f)].centroid =
            (mode == CentroidMode::Legacy2003Middle) ? ((va + vb) * 0.5f + vc) * 0.5f
                                                     : (va + vb + vc) * (1.0f / 3.0f);
    }
    return geom;
}


/// Per-face chart assignment (nullptr = unassigned).
using ChartMap = std::vector<Chart *>;


using detail::BoundaryFlags;


/// All state of one `segment()` run (port of 2003 `CChartification`).
class Segmenter {
public:
    Segmenter(const decimate::HalfEdgeMesh &mesh,
              const decimate::MeshData     &data,
              int                           maxCharts,
              double                        straightenThreshold,
              CentroidMode                  centroidMode)
        : mesh_(mesh),
          data_(data),
          maxCharts_(maxCharts),
          straightenThreshold_(straightenThreshold),
          centroidMode_(centroidMode),
          numFaces_(static_cast<int>(mesh.faceCount())),
          geom_(computeFaceGeometry(mesh, data, centroidMode)),
          chartMap_(static_cast<std::size_t>(numFaces_), nullptr) {
        boundary_.init(mesh);
        labelComponents();
    }

    /// The whole Chartify run: seed, grow and relocate, clean up, straighten, assemble (2003
    /// `CChartification::Chartify`, `legacy2003/chartification/cchartification.cpp:74-114`).
    ChartSet run();

    /// Spec 19-03: split chart @p chartId of @p cs at @p foldFace (see `splitChart`).
    void split(ChartSet &cs, int chartId, int foldFace);

private:
    const decimate::HalfEdgeMesh &mesh_;
    const decimate::MeshData     &data_;
    int                           maxCharts_;
    double                        straightenThreshold_;
    CentroidMode                  centroidMode_;
    int                           numFaces_;

    std::vector<FaceGeometry>           geom_;
    std::vector<std::unique_ptr<Chart>> charts_;
    ChartMap                            chartMap_;
    BoundaryFlags                       boundary_;
    int                                 lastFace_ = 0;
    std::vector<int>                    compOf_;    ///< Dual-graph component id per face.
    std::vector<int>                    compFirst_; ///< Lowest face of each component (ascending).

    /// Push the face across half-edge @p h if it exists and has no component yet.
    void claimNeighbour(decimate::HeIndex h, int comp, std::vector<int> &stack) {
        const decimate::HeIndex s = mesh_.sym(h);
        if (s == decimate::kNoHalfedge) {
            return;
        }
        const int nb = static_cast<int>(mesh_.halfedge(s).face);
        if (compOf_[static_cast<std::size_t>(nb)] >= 0) {
            return;
        }
        compOf_[static_cast<std::size_t>(nb)] = comp;
        stack.push_back(nb);
    }

    /// Depth-first flood of the component of @p first over sym edges.
    void floodComponent(int first, int comp) {
        std::vector<int> stack { first };
        compOf_[static_cast<std::size_t>(first)] = comp;
        while (!stack.empty()) {
            const int cur = stack.back();
            stack.pop_back();
            for (int e = 0; e < 3; ++e) {
                claimNeighbour(static_cast<decimate::HeIndex>(cur * 3 + e), comp, stack);
            }
        }
    }

    /// Label dual-graph components; component ids ascend with their lowest face.
    void labelComponents() {
        compOf_.assign(static_cast<std::size_t>(numFaces_), -1);
        for (int f = 0; f < numFaces_; ++f) {
            if (compOf_[static_cast<std::size_t>(f)] >= 0) {
                continue;
            }
            compFirst_.push_back(f);
            floodComponent(f, static_cast<int>(compFirst_.size()) - 1);
        }
    }

    /// Cost of crossing half-edge @p h: (1 − cosAngle) × centroidDistance (2003 `EdgeCost`,
    /// `legacy2003/chartification/cchartification.cpp:120-135`).
    ///
    /// The normal for the claimed face side uses the chart's running average
    /// (adapts dynamically as the chart grows, mirroring legacy exactly).
    [[nodiscard]] double edgeCost(decimate::HeIndex h) const {
        const decimate::HeIndex s = mesh_.sym(h);
        assert(s != decimate::kNoHalfedge);

        const int f = static_cast<int>(mesh_.halfedge(h).face);
        const int fopp = static_cast<int>(mesh_.halfedge(s).face);

        math::Vec3f nrm1 = geom_[static_cast<std::size_t>(f)].normal;
        math::Vec3f nrm2 = geom_[static_cast<std::size_t>(fopp)].normal;

        // Switch to chart running normal (legacy EdgeCost exact port).
        // Branches update different variables (nrm1 vs nrm2) — not a clone.
        // NOLINT(bugprone-branch-clone): distinct assignments; symmetric structure
        // mirrors the legacy code exactly.
        if (chartMap_[static_cast<std::size_t>(f)] != nullptr) { // NOLINT(bugprone-branch-clone)
            nrm1 = chartMap_[static_cast<std::size_t>(f)]->normal;
        } else if (chartMap_[static_cast<std::size_t>(fopp)] != nullptr) {
            nrm2 = chartMap_[static_cast<std::size_t>(fopp)]->normal;
        }

        // Cos(angle) between normals via dot product / geometric mean of lengths.
        // Denominator is a FLOAT product and float sqrt, as `legacy2003/chartification/cchartification.cpp:131`
        // (`sqrt(nrm1.length2() * nrm2.length2())`, Vec3 length2() returns float).
        const float l2a = nrm1.lengthSquared();
        const float l2b = nrm2.lengthSquared();
        double      nrmFactor = 0.0;
        if (l2a > 0.0f && l2b > 0.0f) {
            nrmFactor = static_cast<double>(nrm1.dot(nrm2)) /
                        static_cast<double>(std::sqrt(l2a * l2b));
        }

        // Centroid-to-centroid distance.
        const math::Vec3f delta =
            geom_[static_cast<std::size_t>(fopp)].centroid -
            geom_[static_cast<std::size_t>(f)].centroid;
        const double dist = static_cast<double>(delta.length());

        return (1.0 - nrmFactor) * dist;
    }

    /// Centroid-to-centroid distance across half-edge @p h (2003 `CChartification::Distance`,
    /// `legacy2003/chartification/cchartification.cpp:141-144`).
    [[nodiscard]] double centroidDistance(decimate::HeIndex h) const {
        const decimate::HeIndex s = mesh_.sym(h);
        assert(s != decimate::kNoHalfedge);
        const int         f = static_cast<int>(mesh_.halfedge(h).face);
        const int         fopp = static_cast<int>(mesh_.halfedge(s).face);
        const math::Vec3f delta =
            geom_[static_cast<std::size_t>(fopp)].centroid -
            geom_[static_cast<std::size_t>(f)].centroid;
        return static_cast<double>(delta.length());
    }

    /// Allocate a new chart seeded at @p seedFace and return it.
    ///
    /// The seed is included in `chart->faces` immediately (08-03 Fix 1 / D1+D3), as the 2003
    /// `AddInfoToChart` added the seed first; every consumer iterating `chart.faces` then sees the
    /// complete triangle set without special-casing.
    Chart *createChart(int seedFace) {
        auto chart = std::make_unique<Chart>();
        chart->id = static_cast<int>(charts_.size());
        // Newborn chart: seed + normal set directly, seedHistory stays EMPTY, as the
        // 2003 CChart ctor does (legacy2003/chartification/cchart.cpp:23-31).
        const auto &seedNormal = geom_[static_cast<std::size_t>(seedFace)].normal;
        chart->seed = seedFace;
        chart->normalAccum = seedNormal;
        chart->normal = seedNormal.normalized();
        chart->faces.push_back(seedFace); // seed is part of the complete face list
        Chart *raw = chart.get();
        charts_.push_back(std::move(chart));
        return raw;
    }

    /// Reset the chart map entries of all charts to their seed faces and push the seeds' outgoing
    /// half-edges (with twins) onto a fresh heap, each with its EdgeCost.
    MinHeap initExpandHeap() {
        MinHeap heap;
        for (auto &cptr : charts_) {
            Chart *chart = cptr.get();
            chart->clear();
            const int seed = chart->seed;
            // clear() drops the seed too; it must come back first (08-03 Fix 1), as 2003
            // `AddInfoToChart` processed the seed face first.
            chart->faces.push_back(seed);
            chartMap_[static_cast<std::size_t>(seed)] = chart;
            for (int e = 0; e < 3; ++e) {
                const decimate::HeIndex h = static_cast<decimate::HeIndex>(seed * 3 + e);
                if (mesh_.sym(h) != decimate::kNoHalfedge) {
                    heap.push({ edgeCost(h), static_cast<int>(h) });
                }
            }
        }
        return heap;
    }

    /// Full simultaneous multi-source region grow (2003 `ExpandCharts`,
    /// `legacy2003/chartification/cchartification.cpp:284-345`).  The full reset per call is
    /// preserved (Beat 2a).  RepairHeap is commented out in 2003 (`:342`), so stale entries stay
    /// (lazy deletion).
    void expandCharts() {
        // Every run starts from scratch: no chart owns a face and every edge is boundary again.
        std::fill(chartMap_.begin(), chartMap_.end(), nullptr);
        boundary_.reset();
        lastFace_ = 0;
        growFromSeeds();
    }

    /// The grow loop proper: seed every chart in `charts_` and claim faces until the heap is
    /// dry.  Faces already claimed in `chartMap_` (by `split` as outside-the-chart blockers)
    /// are skipped by the lazy-deletion test, so the grow cannot leave its region.
    void growFromSeeds() {
        MinHeap heap = initExpandHeap();

        // Lazy deletion: entries are never updated, a stale one is recognised by its claimed target face.
        while (!heap.empty()) {
            const HeapEntry top = heap.top();
            heap.pop();

            const decimate::HeIndex h = static_cast<decimate::HeIndex>(top.heIndex);
            const decimate::HeIndex s = mesh_.sym(h);
            if (s == decimate::kNoHalfedge) {
                continue; // mesh boundary: nothing to grow into
            }

            const int f = static_cast<int>(mesh_.halfedge(h).face);
            const int fopp = static_cast<int>(mesh_.halfedge(s).face);

            if (chartMap_[static_cast<std::size_t>(fopp)] != nullptr) {
                continue; // stale: the face across was claimed after this entry was pushed
            }

            Chart *chart = chartMap_[static_cast<std::size_t>(f)];
            if (chart == nullptr) {
                continue; // a pre-claimed blocker or an unseeded face has no chart to grow
            }

            chart->grow(geom_[static_cast<std::size_t>(fopp)].normal);
            chart->faces.push_back(fopp);
            chartMap_[static_cast<std::size_t>(fopp)] = chart;
            lastFace_ = fopp;

            // The shared edge is interior now; the claimed face's two other edges may be new frontier.
            boundary_.erase(h);

            // Push the two new boundary half-edges from fopp's remaining edges.
            decimate::HeIndex loop = mesh_.next(s);
            for (int i = 0; i < 2; ++i) {
                if (mesh_.sym(loop) != decimate::kNoHalfedge &&
                    boundary_.isBoundary(loop)) {
                    heap.push({ edgeCost(loop), static_cast<int>(loop) });
                }
                loop = mesh_.next(loop);
            }
            // RepairHeap commented out in 2003 (`legacy2003/chartification/cchartification.cpp:342`) — stale entries remain.
        }
    }

    /// Seed the Lloyd search heap with the intra-chart neighbours of every boundary half-edge
    /// (2003 `InitializeSeedSearch`, `legacy2003/chartification/cchartification.cpp:150-176`).
    MinHeap initializeSeedSearch(std::vector<bool> &visited) {
        const int nhe = static_cast<int>(mesh_.halfEdgeCount());
        std::fill(visited.begin(), visited.end(), false);

        MinHeap heap;
        for (int i = 0; i < nhe; ++i) {
            const decimate::HeIndex h = static_cast<decimate::HeIndex>(i);
            if (!boundary_.isBoundary(h)) {
                continue; // the search starts from the chart borders and works inwards
            }
            const int f = static_cast<int>(mesh_.halfedge(h).face);
            Chart    *chart = chartMap_[static_cast<std::size_t>(f)];
            if (chart == nullptr) {
                continue; // an unassigned face has no chart whose search could start here
            }
            visited[static_cast<std::size_t>(f)] = true;
            chart->searchIndex = f;

            // Push intra-chart neighbours (next and next->next of h).
            decimate::HeIndex loop = mesh_.next(h);
            for (int k = 0; k < 2; ++k) {
                const decimate::HeIndex ls = mesh_.sym(loop);
                if (ls == decimate::kNoHalfedge) {
                    loop = mesh_.next(loop);
                    continue;
                }
                const int loopFace = static_cast<int>(mesh_.halfedge(loop).face);
                const int loopSymFace = static_cast<int>(mesh_.halfedge(ls).face);
                if (!visited[static_cast<std::size_t>(loopSymFace)] &&
                    chartMap_[static_cast<std::size_t>(loopFace)] ==
                        chartMap_[static_cast<std::size_t>(loopSymFace)]) {
                    heap.push({ centroidDistance(loop), static_cast<int>(loop) });
                }
                loop = mesh_.next(loop);
            }
        }
        return heap;
    }

    /// Adopt the searched face of each chart as its new seed; true when every chart hit a seed cycle
    /// (2003 `ApplySeeds`, `legacy2003/chartification/cchartification.cpp:182-204`).
    bool applySeeds() {
        bool cycleDetected = true;
        for (auto &cptr : charts_) {
            Chart    *chart = cptr.get();
            const int seed = chart->searchIndex;
            if (seed >= 0 && seed < numFaces_) {
                const bool cycle =
                    chart->setSeed(seed, geom_[static_cast<std::size_t>(seed)].normal);
                cycleDetected = cycleDetected && cycle;
            } else {
                cycleDetected = false;
            }
        }
        return cycleDetected;
    }

    /// One Lloyd relocation pass: move each chart's seed to its deepest face (2003 `ComputeSeeds`,
    /// `legacy2003/chartification/cchartification.cpp:205-247`).  True when the seeds converged.
    bool computeSeeds() {
        std::vector<bool> visited(static_cast<std::size_t>(numFaces_), false);
        MinHeap           heap = initializeSeedSearch(visited);

        while (!heap.empty()) {
            const HeapEntry top = heap.top();
            heap.pop();
            const double value = top.cost;

            const decimate::HeIndex h = static_cast<decimate::HeIndex>(top.heIndex);
            const decimate::HeIndex s = mesh_.sym(h);
            if (s == decimate::kNoHalfedge) {
                continue;
            }
            const int fopp = static_cast<int>(mesh_.halfedge(s).face);
            if (visited[static_cast<std::size_t>(fopp)]) {
                continue; // a shorter path reached this face first
            }
            visited[static_cast<std::size_t>(fopp)] = true;

            Chart *chart = chartMap_[static_cast<std::size_t>(fopp)];
            if (chart == nullptr) {
                continue;
            }

            // Pops arrive in ascending distance, so the last face popped per chart is its deepest one.
            chart->searchIndex = fopp;

            decimate::HeIndex loop = mesh_.next(s);
            for (int i = 0; i < 2; ++i) {
                const decimate::HeIndex ls = mesh_.sym(loop);
                if (ls == decimate::kNoHalfedge) {
                    loop = mesh_.next(loop);
                    continue;
                }
                const int loopFace = static_cast<int>(mesh_.halfedge(loop).face);
                const int loopSymFace = static_cast<int>(mesh_.halfedge(ls).face);
                if (!visited[static_cast<std::size_t>(loopSymFace)] &&
                    chartMap_[static_cast<std::size_t>(loopFace)] ==
                        chartMap_[static_cast<std::size_t>(loopSymFace)]) {
                    const double dist =
                        static_cast<double>(
                            (geom_[static_cast<std::size_t>(loopSymFace)].centroid -
                             geom_[static_cast<std::size_t>(loopFace)].centroid)
                                .length());
                    heap.push({ value + dist, static_cast<int>(loop) });
                }
                loop = mesh_.next(loop);
            }
        }

        // A small chart has no meaningful deepest face: keep its seed (mirrors the 2003 `ClusterCount() < 4` check, `legacy2003/chartification/cchartification.cpp:244`).
        for (auto &cptr : charts_) {
            Chart *chart = cptr.get();
            if (chart->clusterCount < 4) {
                chart->searchIndex = chart->seed;
            }
        }

        return applySeeds();
    }

    /// Remove charts with clusterCount < 2 (fewer than 2 grow calls); port of 2003
    /// `CChartSet::RemoveInvalidCharts` (`legacy2003/chartification/cchartset.cpp:71-91`).
    /// Returns true if any were removed.
    bool removeInvalidCharts() {
        const std::size_t before = charts_.size();
        charts_.erase(
            std::remove_if(charts_.begin(), charts_.end(),
                           [](const std::unique_ptr<Chart> &c) {
                               return c->clusterCount < 2;
                           }),
            charts_.end());
        // `assemble` indexes charts by id, so the ids must stay contiguous.
        for (std::size_t i = 0; i < charts_.size(); ++i) {
            charts_[i]->id = static_cast<int>(i);
        }
        return charts_.size() < before;
    }

    /// Cosine error `1 - cos` between two face normals; 1 when either is degenerate.
    [[nodiscard]] double normalDeviation(int faceA, int faceB) const {
        const auto  &na = geom_[static_cast<std::size_t>(faceA)].normal;
        const auto  &nb = geom_[static_cast<std::size_t>(faceB)].normal;
        const double la = static_cast<double>(na.lengthSquared());
        const double lb = static_cast<double>(nb.lengthSquared());
        if (la <= 0.0 || lb <= 0.0) {
            return 1.0;
        }
        return 1.0 - static_cast<double>(na.dot(nb)) / std::sqrt(la * lb);
    }

    /// Try to hand face @p i to the chart across the two boundary edges starting at @p h.
    /// Returns true if the face was reassigned (boundary modified).
    [[nodiscard]] bool tryStraightenEdge(int i, decimate::HeIndex h) {
        const decimate::HeIndex hnxt = mesh_.next(h);
        if (!(boundary_.isBoundary(h) && boundary_.isBoundary(hnxt))) {
            return false; // only a face with two consecutive border edges is a notch
        }
        const decimate::HeIndex hs = mesh_.sym(h);
        const decimate::HeIndex hnxs = mesh_.sym(hnxt);
        if (hs == decimate::kNoHalfedge || hnxs == decimate::kNoHalfedge) {
            return false;
        }
        Chart *ochart1 = chartMap_[static_cast<std::size_t>(mesh_.halfedge(hs).face)];
        Chart *ochart2 = chartMap_[static_cast<std::size_t>(mesh_.halfedge(hnxs).face)];
        if (ochart1 == nullptr || ochart2 == nullptr) {
            return false;
        }
        if (ochart1 != ochart2 || chartMap_[static_cast<std::size_t>(i)] == ochart1) {
            return false;
        }
        // Default 0.3 = legacy2003/chartification/cchartification.cpp:542 (`if (err < 0.3)`).
        if (normalDeviation(i, static_cast<int>(mesh_.halfedge(hs).face)) < straightenThreshold_) {
            // Move the face over: its third edge becomes border, the two it was bordered by dissolve.
            chartMap_[static_cast<std::size_t>(i)] = ochart1;
            boundary_.setFlag(mesh_.next(mesh_.next(h)), true);
            boundary_.erase(h);
            // Divergence from 2003 (`legacy2003/chartification/cchartification.cpp:547-548` stops after Erase(h)): next(h)
            // now lies inside the chart and must not stay behind as a self-seam.
            boundary_.erase(hnxt);
            return true;
        }
        return false;
    }

    /// Boundary-smoothing pass (2003 `Straighten`, `legacy2003/chartification/cchartification.cpp:505-558`):
    /// reassign a face with two consecutive boundary edges on the same opposite chart iff
    /// normalError < straightenThreshold_.  True when any face moved.
    bool straighten() {
        bool changed = false;
        for (int i = 0; i < numFaces_; ++i) {
            const Chart *chart = chartMap_[static_cast<std::size_t>(i)];
            if (chart == nullptr || i == chart->seed) {
                continue; // the seed anchors its chart
            }

            // 2003 starts at h->next (`legacy2003/chartification/cchartification.cpp:519-521`): order 1, 2, 0.
            for (int j = 0; j < 3; ++j) {
                if (tryStraightenEdge(i, static_cast<decimate::HeIndex>(i * 3 + (j + 1) % 3))) {
                    changed = true;
                    break;
                }
            }
        }
        return changed;
    }

    /// One chart per dual-graph component at its lowest face.  For a single-component mesh this is
    /// exactly `createChart(0)`, the 2003 bootstrap
    /// (`legacy2003/chartification/cchartification.cpp:94`).  Created regardless of `maxCharts_`.
    void seedComponents() {
        for (const int first : compFirst_) {
            createChart(first);
        }
    }

    /// Give every component that lost all its charts one chart owning it whole.
    /// After the removal loop `chartMap_` is fresh; a component whose lowest face is
    /// unassigned lost every chart (all had < 2 grown faces) and gets one chart owning
    /// all its faces.  Port addition: 2003 leaves such faces unassigned (`legacy2003/chartification/cchartification.cpp:437`).
    void fallbackPerComponent() {
        std::vector<Chart *> fallback(compFirst_.size(), nullptr);
        for (std::size_t c = 0; c < compFirst_.size(); ++c) {
            if (chartMap_[static_cast<std::size_t>(compFirst_[c])] == nullptr) {
                fallback[c] = createChart(compFirst_[c]);
            }
        }
        for (int f = 0; f < numFaces_; ++f) {
            Chart *chart = fallback[static_cast<std::size_t>(compOf_[static_cast<std::size_t>(f)])];
            if (chart == nullptr) {
                continue;
            }
            chartMap_[static_cast<std::size_t>(f)] = chart;
            // The last regrow had no charts here, so every flag is still set; an unflagged
            // fallback chart is one piece, not a heap of single-face islands (Spec 19-02).
            for (int e = 0; e < 3; ++e) {
                boundary_.setFlag(static_cast<decimate::HeIndex>(f * 3 + e), false);
            }
        }
    }

    /// Build the final `ChartSet` from `chartMap_`.
    ///
    ChartSet assemble() {
        ChartSet result;
        result.charts.reserve(charts_.size());
        result.centroidMode = centroidMode_;

        for (auto &cptr : charts_) {
            result.charts.push_back(*cptr);
        }


        // Chart::faces from the grow goes stale (straighten and the fallback move faces), so
        // single source of truth: derive both Chart.faces and faceToChart from
        // chartMap_ (the authoritative face->chart assignment). The Chart.faces
        // built incrementally during grow goes stale when straighten() and
        // fallbackPerComponent() reassign chartMap_ without touching it; rebuilding
        // here keeps the two representations consistent by construction.
        // (chart->id == index in result.charts — see createChart / compaction.)
        for (Chart &ch : result.charts) {
            ch.faces.clear();
        }
        result.faceToChart.resize(static_cast<std::size_t>(numFaces_), -1);

        for (int f = 0; f < numFaces_; ++f) {
            const Chart *chart = chartMap_[static_cast<std::size_t>(f)];
            if (chart != nullptr) {
                result.faceToChart[static_cast<std::size_t>(f)] = chart->id;
                result.charts[static_cast<std::size_t>(chart->id)].faces.push_back(f);
            }
        }


        exportSeams(result);
        return result;
    }

    /// The half-edges of @p seams whose two faces both belong to chart @p id of @p cs.
    [[nodiscard]] std::vector<int> seamsInside(const std::vector<int> &seams, const ChartSet &cs,
                                               int id) const {
        std::vector<int> out;
        for (const int h : seams) {
            const auto he = static_cast<decimate::HeIndex>(h);
            const bool both = cs.faceToChart[mesh_.halfedge(he).face] == id &&
                              cs.faceToChart[mesh_.halfedge(mesh_.sym(he)).face] == id;
            if (both) {
                out.push_back(h);
            }
        }
        return out;
    }

    /// Fill `Chart::seamHalfEdges`: half-edges with a sym, both faces in the same chart, still
    /// flagged boundary.  Half-edge indices ascend, so each list is sorted.
    void exportSeams(ChartSet &cs) const {
        const int nhe = static_cast<int>(mesh_.halfEdgeCount());
        for (int i = 0; i < nhe; ++i) {
            const auto              h = static_cast<decimate::HeIndex>(i);
            const decimate::HeIndex s = mesh_.sym(h);
            if (s == decimate::kNoHalfedge || !boundary_.isBoundary(h)) {
                continue;
            }
            const int c = cs.faceToChart[mesh_.halfedge(h).face];
            if (c >= 0 && c == cs.faceToChart[mesh_.halfedge(s).face]) {
                cs.charts[static_cast<std::size_t>(c)].seamHalfEdges.push_back(i);
            }
        }
    }
};


/// Restricted two-seed regrow inside one chart (Spec 19-03): the old chart's faces are re-grown from its seed and @p foldFace, nothing else moves.
void Segmenter::split(ChartSet &cs, int chartId, int foldFace) {
    Chart &old = cs.charts[static_cast<std::size_t>(chartId)];

    // Two fresh charts seeded at the old seed and the fold face; every face outside the old
    // chart is pre-claimed by a sentinel so the grow cannot leave it.
    Chart *keep = createChart(old.seed);
    Chart *fold = createChart(foldFace);
    Chart  blocker;
    std::fill(chartMap_.begin(), chartMap_.end(), &blocker);
    for (const int f : old.faces) {
        chartMap_[static_cast<std::size_t>(f)] = nullptr;
    }
    boundary_.reset();
    growFromSeeds();

    // Rewrite both face lists ascending; a face the grow could not reach stays with the old chart.
    const int        newId = static_cast<int>(cs.charts.size());
    std::vector<int> keepFaces;
    std::vector<int> foldFaces;
    for (const int f : old.faces) {
        const bool toFold = chartMap_[static_cast<std::size_t>(f)] == fold;
        (toFold ? foldFaces : keepFaces).push_back(f);
    }
    for (const int f : foldFaces) {
        cs.faceToChart[static_cast<std::size_t>(f)] = newId;
    }

    // Each piece keeps the parent's seams that lie wholly inside it (Spec 19-03 A1); the edges
    // between the pieces are inter-chart now and never seams.
    keep->seamHalfEdges = seamsInside(old.seamHalfEdges, cs, chartId);
    fold->seamHalfEdges = seamsInside(old.seamHalfEdges, cs, newId);
    keep->id = chartId;
    keep->faces = std::move(keepFaces);
    keep->seedHistory = old.seedHistory;
    fold->id = newId;
    fold->faces = std::move(foldFaces);
    cs.charts[static_cast<std::size_t>(chartId)] = *keep;
    cs.charts.push_back(*fold);
}


ChartSet Segmenter::run() { // NOLINT(readability-function-size): direct port of the Chartify loop
    if (numFaces_ <= 0) {
        return {}; // nothing to segment
    }

    // Bootstrap: one seed per dual-graph component (== createChart(0) on a connected mesh,
    // as 2003 `legacy2003/chartification/cchartification.cpp:94`).
    seedComponents();

    // Main chartification loop.
    // Guards vs legacy: skip creating a new chart when lastFace_ is already a
    // seed (avoids duplicate-seed infinite loops on tiny meshes).
    //
    // Iteration cap (08-03 Fix 3 / D2): a SMALL CONSTANT multiple of maxCharts_,
    // not scaled by numFaces_.  The natural exit is computeSeeds() convergence
    // (the per-chart seedHistory cycle detection is the authority); this cap is a
    // backstop against a genuine infinite loop only.  Legacy (do { } while
    // (!ComputeSeeds())) had NO cap and converged in ~O(charts) passes (~25).
    // Using 10× maxCharts_ is more than enough for any realistic convergence
    // count and cannot reintroduce D2: even at maxCharts_=64 the cap is 640
    // outer iterations (vs the old 8.6M for Armadillo).
    const int maxIter = std::max(kMaxIterFactor * maxCharts_, 64);
    int       iter = 0;
    do {
        expandCharts();
        if (static_cast<int>(charts_.size()) < maxCharts_) {
            // A new chart at an existing seed would make two identical charts.
            bool alreadySeed = false;
            for (const auto &cptr : charts_) {
                if (cptr->seed == lastFace_) {
                    alreadySeed = true;
                    break;
                }
            }
            if (!alreadySeed) {
                createChart(lastFace_);
            }
        }
        ++iter;
    } while (!computeSeeds() && iter < maxIter);

    // Remove charts with clusterCount < 2; re-grow until stable (2003 `legacy2003/chartification/cchartification.cpp:101-104`).  The regrow
    // also runs with zero charts left, which clears the dangling chartMap_ entries.
    // There is deliberately NO unconditional regrow after this loop: the last computeSeeds()
    // moved the seeds, and 2003 keeps the chart map grown from the pre-move seeds
    // (`legacy2003/chartification/cchartification.cpp:97-105`), so a final expandCharts() would change the partition.
    while (removeInvalidCharts()) {
        expandCharts();
    }

    fallbackPerComponent();

    // Smooth the chart borders until no face wants to move.  2003 loops unbounded
    // (`legacy2003/chartification/cchartification.cpp:105`); a threshold near 2.0 can let two faces trade
    // places forever, so the loop gets the main loop's cap (deliberate deviation, unreachable on the goldens).
    detail::runBounded([this] { return straighten(); }, maxIter);

    return assemble();
}

} // namespace


void Chart::grow(const math::Vec3f &faceNormal) {
    normalAccum += faceNormal;
    normal = normalAccum.normalized();
    clusterCount++;
}


bool Chart::setSeed(int newSeed, const math::Vec3f &faceNormal) {
    seed = newSeed;
    normalAccum = faceNormal;
    normal = faceNormal.normalized();

    const bool cycle = seedHistory.count(newSeed) > 0;
    seedHistory.insert(newSeed);
    return cycle; // true = cycle detected = convergence
}


void Chart::clear() {
    faces.clear();
    clusterCount = 0;
    // keep seed, normal, normalAccum, seedHistory, id, searchIndex
}


ChartSet segment(const decimate::HalfEdgeMesh &mesh,
                 const decimate::MeshData     &data,
                 int                           maxCharts,
                 double                        straightenThreshold,
                 CentroidMode                  centroidMode) {
    Segmenter seg(mesh, data, maxCharts, straightenThreshold, centroidMode);
    return seg.run();
}


void splitChart(const decimate::HalfEdgeMesh &mesh,
                const decimate::MeshData     &data,
                ChartSet                     &charts,
                int                           chartId,
                int                           foldFace,
                CentroidMode                  centroidMode) {
    if (chartId < 0 || chartId >= static_cast<int>(charts.charts.size())) {
        throw std::invalid_argument("splitChart: chartId out of range");
    }

    // Both pieces need their own seed: the old seed stays with the first, the fold face starts the second.
    const Chart &c = charts.charts[static_cast<std::size_t>(chartId)];
    if (foldFace == c.seed) {
        throw std::invalid_argument("splitChart: foldFace is the chart's seed");
    }
    if (std::find(c.faces.begin(), c.faces.end(), foldFace) == c.faces.end()) {
        throw std::invalid_argument("splitChart: foldFace is not a face of the chart");
    }

    Segmenter seg(mesh, data, /*maxCharts=*/2, kStraightenThreshold, centroidMode);
    seg.split(charts, chartId, foldFace);
}

} // namespace decimate::chart
