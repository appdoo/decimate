// This file is part of the decimate library; see decimate.h for version/license details.
/// @file segmenter.cpp
/// @brief Region-growing chart segmentation — port of `CChartification`.
///
/// ## Behavior-preservation notes (Beat 2a)
///
/// ### Legacy Chartify loop (cchartification.cpp:74-114)
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
/// ### EdgeCost (legacy cchartification.cpp:120-135)
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
/// ### Grow-loop pop/claim logic (legacy cchartification.cpp:310-345)
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
/// Legacy `CChartSet` uses `std::set<CChart*>` — pointer-ordered.
/// We use `std::vector<unique_ptr<Chart>>` with sequential access matching
/// chart creation order (face 0 first, lastface after each iteration).
/// Creation order is deterministic for a fixed input — output is identical
/// across runs (T-S1).
///
/// Legacy `CSet<int>` (seedHistory) is sorted — mapped to `std::set<int>`.
///
/// ### ClusterCount vs faces
///
/// Legacy `CChart::clusterCount` = number of `Grow()` calls.
/// Legacy `RemoveInvalidCharts` removes charts with `clusterCount < 2`.
/// A chart starts with just a seed (clusterCount=0); each face grown in
/// increments clusterCount by 1.  We track this as `Chart::clusterCount`.

#include <decimate/chart/chart.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <vector>

#include <decimate/math/vec3.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>

namespace decimate::chart {

// ===========================================================================
// Chart member implementations
// ===========================================================================

void Chart::grow(const decimate::math::Vec3<float>& faceNormal) {
    normalAccum += faceNormal;
    normal       = normalAccum.normalized();
    clusterCount++;
}

bool Chart::setSeed(int newSeed, const decimate::math::Vec3<float>& faceNormal) {
    seed        = newSeed;
    normalAccum = faceNormal;
    normal      = faceNormal.normalized();

    const bool cycle = seedHistory.count(newSeed) > 0;
    seedHistory.insert(newSeed);
    return cycle;  // true = cycle detected = convergence
}

void Chart::clear() {
    faces.clear();
    clusterCount = 0;
    // keep seed, normal, normalAccum, seedHistory, id, searchIndex
}

// ===========================================================================
// Internal types and helpers
// ===========================================================================

namespace {

// ---------------------------------------------------------------------------
// HeapEntry — lazy-deletion min-heap element
//
// Legacy: Heap<Halfedge*, double> — min-heap, pop returns the minimum cost.
// Port: std::priority_queue with std::greater → min-heap.
// ---------------------------------------------------------------------------
struct HeapEntry {
    double  cost;
    int     heIndex;  // half-edge index in HalfEdgeMesh

    bool operator>(const HeapEntry& o) const noexcept { return cost > o.cost; }
};

using MinHeap = std::priority_queue<HeapEntry,
                                    std::vector<HeapEntry>,
                                    std::greater<HeapEntry>>;

// ---------------------------------------------------------------------------
// Pre-computed per-face normals and centroids
// ---------------------------------------------------------------------------
struct FaceGeometry {
    decimate::math::Vec3<float> normal;    ///< Raw (unnormalised) face normal.
    decimate::math::Vec3<float> centroid;  ///< Face centroid.
};

/// Compute raw face normal and centroid for each face.
/// Raw normal = (b-a) × (c-a) — unnormalised.
/// Matches legacy `normal_raw(points, &faces[i])` and `middle(points, &faces[i])`.
std::vector<FaceGeometry> computeFaceGeometry(
        const decimate::HalfEdgeMesh& mesh,
        const decimate::MeshData&     data) {

    const int numFaces = static_cast<int>(mesh.faceCount());
    std::vector<FaceGeometry> geom(static_cast<std::size_t>(numFaces));

    for (int f = 0; f < numFaces; ++f) {
        // Legacy WriteFace wiring (half_edge_mesh.cpp):
        //   he[f*3+0]: end=c, next=f*3+1
        //   he[f*3+1]: end=b, next=f*3+2
        //   he[f*3+2]: end=a, next=f*3+0
        // So: a=he[f*3+2].end, b=he[f*3+1].end, c=he[f*3+0].end.
        const auto& he0 = mesh.halfedge(static_cast<decimate::HeIndex>(f * 3 + 0));
        const auto& he1 = mesh.halfedge(static_cast<decimate::HeIndex>(f * 3 + 1));
        const auto& he2 = mesh.halfedge(static_cast<decimate::HeIndex>(f * 3 + 2));

        auto getPos = [&](std::uint32_t vi) -> decimate::math::Vec3<float> {
            return data.positions[static_cast<std::size_t>(vi)];
        };

        const decimate::math::Vec3<float> va = getPos(he2.end);
        const decimate::math::Vec3<float> vb = getPos(he1.end);
        const decimate::math::Vec3<float> vc = getPos(he0.end);

        geom[static_cast<std::size_t>(f)].normal   = (vb - va).cross(vc - va);
        geom[static_cast<std::size_t>(f)].centroid = (va + vb + vc) * (1.0f / 3.0f);
    }
    return geom;
}

// ---------------------------------------------------------------------------
// ChartMap — per-face chart assignment (nullptr = unassigned)
// ---------------------------------------------------------------------------
using ChartMap = std::vector<Chart*>;

// ---------------------------------------------------------------------------
// BoundaryFlags — port of CChartBoundary (boolean flag per half-edge + vertex degrees)
//
// The per-ExpandCharts full reset (O(numHalfedges)) is preserved (Beat 2a).
// Generation-counter optimization is a Beat-2b forbidden path.
// ---------------------------------------------------------------------------
struct BoundaryFlags {
    std::vector<bool> flags;     ///< One flag per half-edge: true = boundary.
    std::vector<int>  degrees;   ///< Per-vertex boundary degree count.
    const decimate::HalfEdgeMesh* mesh = nullptr;

    void init(const decimate::HalfEdgeMesh& m) {
        mesh = &m;
        const int nhe = static_cast<int>(m.halfEdgeCount());
        const int np  = static_cast<int>(m.pointCount());
        flags.assign(static_cast<std::size_t>(nhe), false);
        degrees.assign(static_cast<std::size_t>(np), 0);
    }

    /// @brief Reset to full-boundary state — O(numHalfedges).
    ///
    /// Mirrors CChartBoundary::Reset(): all flags = true; degrees computed
    /// from unique undirected edges (each edge processed once via the "lower" rule).
    void reset() {
        const int nhe = static_cast<int>(mesh->halfEdgeCount());
        const int np  = static_cast<int>(mesh->pointCount());

        std::fill(flags.begin(), flags.end(), true);
        std::fill(degrees.begin(), degrees.end(), 0);

        std::vector<bool> visited(static_cast<std::size_t>(nhe), false);

        for (int i = 0; i < nhe; ++i) {
            if (visited[static_cast<std::size_t>(i)]) {
                continue;
            }
            decimate::HeIndex lo = static_cast<decimate::HeIndex>(i);
            const decimate::HeIndex s = mesh->sym(lo);
            if (s != decimate::kNoHalfedge &&
                mesh->halfedge(s).face < mesh->halfedge(lo).face) {
                lo = s;
            }
            visited[static_cast<std::size_t>(lo)] = true;
            if (s != decimate::kNoHalfedge) {
                visited[static_cast<std::size_t>(s)] = true;
            }

            const auto& he   = mesh->halfedge(lo);
            const std::uint32_t endV =
                mesh->halfedge(mesh->next(mesh->next(lo))).end;
            degrees[static_cast<std::size_t>(he.end)]++;
            degrees[static_cast<std::size_t>(endV)]++;
        }
    }

    [[nodiscard]] bool isBoundary(decimate::HeIndex h) const noexcept {
        return flags[static_cast<std::size_t>(lower(h))];
    }

    /// @brief Erase half-edge h from boundary — mirrors CChartBoundary::Erase.
    void erase(decimate::HeIndex h) {
        const decimate::HeIndex lo = lower(h);
        if (!flags[static_cast<std::size_t>(lo)]) {
            return;
        }
        flags[static_cast<std::size_t>(lo)] = false;
        const decimate::HeIndex s = mesh->sym(lo);
        if (s != decimate::kNoHalfedge) {
            flags[static_cast<std::size_t>(s)] = false;
        }

        const auto& he           = mesh->halfedge(lo);
        const std::uint32_t endV   = he.end;
        const std::uint32_t startV = mesh->halfedge(mesh->next(mesh->next(lo))).end;

        if (degrees[static_cast<std::size_t>(endV)] > 0) {
            degrees[static_cast<std::size_t>(endV)]--;
        }
        if (degrees[static_cast<std::size_t>(startV)] > 0) {
            degrees[static_cast<std::size_t>(startV)]--;
        }

        if (degrees[static_cast<std::size_t>(endV)] < 2) {
            decimate::HeIndex inv = searchBoundary(lo);
            if (inv != decimate::kNoHalfedge) {
                erase(inv);
            }
        }
        if (degrees[static_cast<std::size_t>(startV)] < 2) {
            decimate::HeIndex nn  = mesh->next(mesh->next(lo));
            decimate::HeIndex inv = searchBoundary(nn);
            if (inv != decimate::kNoHalfedge) {
                erase(inv);
            }
        }
    }

    /// @brief Direct set — mirrors CChartBoundary::SetBoundaryFlag.
    void setFlag(decimate::HeIndex h, bool boundary) {
        if (boundary == isBoundary(h)) {
            return;
        }
        const decimate::HeIndex lo = lower(h);
        flags[static_cast<std::size_t>(lo)] = boundary;
        const decimate::HeIndex s = mesh->sym(lo);
        if (s != decimate::kNoHalfedge) {
            flags[static_cast<std::size_t>(s)] = boundary;
        }
        const auto&         he     = mesh->halfedge(lo);
        const int           delta  = boundary ? 1 : -1;
        const std::uint32_t endV   = he.end;
        const std::uint32_t startV = mesh->halfedge(mesh->next(mesh->next(lo))).end;
        degrees[static_cast<std::size_t>(endV)]   += delta;
        degrees[static_cast<std::size_t>(startV)] += delta;
    }

private:
    [[nodiscard]] decimate::HeIndex lower(decimate::HeIndex h) const noexcept {
        const decimate::HeIndex s = mesh->sym(h);
        if (s != decimate::kNoHalfedge &&
            mesh->halfedge(s).face < mesh->halfedge(h).face) {
            return s;
        }
        return h;
    }

    /// Find a boundary half-edge in the 1-ring around `start`.
    /// Mirrors CChartBoundary::SearchBoundary.
    [[nodiscard]] decimate::HeIndex searchBoundary(decimate::HeIndex start) const {
        // Direction 1: loop via next->sym
        decimate::HeIndex loop = start;
        bool done1 = false;
        do {
            if (isBoundary(loop)) {
                return loop;
            }
            const decimate::HeIndex nxt = mesh->next(loop);
            if (isBoundary(nxt)) {
                return nxt;
            }
            const decimate::HeIndex nxtSym = mesh->sym(nxt);
            if (nxtSym == decimate::kNoHalfedge) {
                done1 = true;
                break;
            }
            loop = nxtSym;
        } while (loop != start);

        if (!done1 && loop == start) {
            return decimate::kNoHalfedge;
        }

        // Direction 2: loop via sym->next->next
        loop = start;
        do {
            if (isBoundary(loop)) {
                return loop;
            }
            const decimate::HeIndex s = mesh->sym(loop);
            if (s != decimate::kNoHalfedge && isBoundary(s)) {
                return s;
            }
            if (s != decimate::kNoHalfedge) {
                loop = mesh->next(mesh->next(s));
            } else {
                break;
            }
        } while (loop != start);

        return decimate::kNoHalfedge;
    }
};

// ---------------------------------------------------------------------------
// Segmenter — encapsulates all state for one run of segment()
// ---------------------------------------------------------------------------
class Segmenter {
public:
    Segmenter(const decimate::HalfEdgeMesh& mesh,
              const decimate::MeshData&     data,
              int                        maxCharts)
        : mesh_(mesh),
          data_(data),
          maxCharts_(maxCharts),
          numFaces_(static_cast<int>(mesh.faceCount())),
          geom_(computeFaceGeometry(mesh, data)),
          chartMap_(static_cast<std::size_t>(numFaces_), nullptr) {
        boundary_.init(mesh);
    }

    ChartSet run();

private:
    const decimate::HalfEdgeMesh& mesh_;
    const decimate::MeshData&     data_;
    int                        maxCharts_;
    int                        numFaces_;

    std::vector<FaceGeometry>            geom_;
    std::vector<std::unique_ptr<Chart>>  charts_;
    ChartMap                             chartMap_;
    BoundaryFlags                        boundary_;
    int                                  lastFace_ = 0;

    // -----------------------------------------------------------------------
    // edgeCost — port of CChartification::EdgeCost (cchartification.cpp:120-135)
    // -----------------------------------------------------------------------
    /// Cost of crossing half-edge @p h: (1 − cosAngle) × centroidDistance.
    ///
    /// The normal for the claimed face side uses the chart's running average
    /// (adapts dynamically as the chart grows, mirroring legacy exactly).
    [[nodiscard]] double edgeCost(decimate::HeIndex h) const {
        const decimate::HeIndex s = mesh_.sym(h);
        assert(s != decimate::kNoHalfedge);

        const int f    = static_cast<int>(mesh_.halfedge(h).face);
        const int fopp = static_cast<int>(mesh_.halfedge(s).face);

        decimate::math::Vec3<float> nrm1 = geom_[static_cast<std::size_t>(f)].normal;
        decimate::math::Vec3<float> nrm2 = geom_[static_cast<std::size_t>(fopp)].normal;

        // Switch to chart running normal (legacy EdgeCost exact port).
        // Branches update different variables (nrm1 vs nrm2) — not a clone.
        // NOLINT(bugprone-branch-clone): distinct assignments; symmetric structure
        // mirrors the legacy code exactly.
        if (chartMap_[static_cast<std::size_t>(f)] != nullptr) {         // NOLINT(bugprone-branch-clone)
            nrm1 = chartMap_[static_cast<std::size_t>(f)]->normal;
        } else if (chartMap_[static_cast<std::size_t>(fopp)] != nullptr) {
            nrm2 = chartMap_[static_cast<std::size_t>(fopp)]->normal;
        }

        // Cos(angle) between normals via dot product / geometric mean of lengths.
        // Legacy: nrm1 * nrm2 was the Vec3 dot (operator*); nrm1.length2() was lengthSquared.
        const double l2a = static_cast<double>(nrm1.lengthSquared());
        const double l2b = static_cast<double>(nrm2.lengthSquared());
        double nrmFactor = 0.0;
        if (l2a > 0.0 && l2b > 0.0) {
            nrmFactor = static_cast<double>(nrm1.dot(nrm2)) / std::sqrt(l2a * l2b);
        }

        // Centroid-to-centroid distance.
        const decimate::math::Vec3<float> delta =
            geom_[static_cast<std::size_t>(fopp)].centroid -
            geom_[static_cast<std::size_t>(f)].centroid;
        const double dist = static_cast<double>(delta.length());

        return (1.0 - nrmFactor) * dist;
    }

    // -----------------------------------------------------------------------
    // centroidDistance — port of CChartification::Distance (:136-140)
    // -----------------------------------------------------------------------
    [[nodiscard]] double centroidDistance(decimate::HeIndex h) const {
        const decimate::HeIndex s = mesh_.sym(h);
        assert(s != decimate::kNoHalfedge);
        const int f    = static_cast<int>(mesh_.halfedge(h).face);
        const int fopp = static_cast<int>(mesh_.halfedge(s).face);
        const decimate::math::Vec3<float> delta =
            geom_[static_cast<std::size_t>(fopp)].centroid -
            geom_[static_cast<std::size_t>(f)].centroid;
        return static_cast<double>(delta.length());
    }

    // -----------------------------------------------------------------------
    // createChart — allocate a new chart with a given seed face
    //
    // The seed is included in chart->faces immediately (08-03 Fix 1 / D1+D3).
    // Legacy AddInfoToChart added the seed first; we mirror that here at the
    // producer so every consumer iterating chart.faces sees the complete
    // triangle set (seed + grown faces), without any per-consumer special-casing.
    // -----------------------------------------------------------------------
    Chart* createChart(int seedFace) {
        auto chart = std::make_unique<Chart>();
        chart->id  = static_cast<int>(charts_.size());
        chart->setSeed(seedFace, geom_[static_cast<std::size_t>(seedFace)].normal);
        chart->faces.push_back(seedFace);  // seed is part of the complete face list
        Chart* raw = chart.get();
        charts_.push_back(std::move(chart));
        return raw;
    }

    // -----------------------------------------------------------------------
    // expandCharts — port of ExpandCharts (cchartification.cpp:284-345)
    // -----------------------------------------------------------------------
    /// Full simultaneous multi-source region grow.
    ///
    /// Per-ExpandCharts full reset preserved (Beat 2a).
    /// @brief Seed the heap from all chart seed faces; reset chartMap_ and boundary_.
    ///
    /// Each chart's seed face is placed in chartMap_; all outgoing half-edges
    /// with twins are pushed onto the heap with their EdgeCost.
    MinHeap initExpandHeap() {
        MinHeap heap;
        for (auto& cptr : charts_) {
            Chart* chart = cptr.get();
            chart->clear();
            const int seed = chart->seed;
            // Re-add seed to faces after clear() (08-03 Fix 1: seed is always
            // the first member of the complete face list, mirroring legacy
            // AddInfoToChart which processed the seed face first).
            chart->faces.push_back(seed);
            chartMap_[static_cast<std::size_t>(seed)] = chart;
            for (int e = 0; e < 3; ++e) {
                const decimate::HeIndex h = static_cast<decimate::HeIndex>(seed * 3 + e);
                if (mesh_.sym(h) != decimate::kNoHalfedge) {
                    heap.push({edgeCost(h), static_cast<int>(h)});
                }
            }
        }
        return heap;
    }

    /// RepairHeap is commented out in legacy → stale entries stay (lazy deletion).
    void expandCharts() {
        // Full reset — O(numFaces) + O(numHalfedges).
        std::fill(chartMap_.begin(), chartMap_.end(), nullptr);
        boundary_.reset();
        lastFace_ = 0;

        MinHeap heap = initExpandHeap();

        // Grow phase — lazy deletion via per-face chartMap_ assignment.
        while (!heap.empty()) {
            const HeapEntry top = heap.top();
            heap.pop();

            const decimate::HeIndex h = static_cast<decimate::HeIndex>(top.heIndex);
            const decimate::HeIndex s = mesh_.sym(h);
            if (s == decimate::kNoHalfedge) {
                continue;
            }

            const int f    = static_cast<int>(mesh_.halfedge(h).face);
            const int fopp = static_cast<int>(mesh_.halfedge(s).face);

            // Lazy deletion: skip stale entries (fopp already claimed).
            if (chartMap_[static_cast<std::size_t>(fopp)] != nullptr) {
                continue;
            }

            Chart* chart = chartMap_[static_cast<std::size_t>(f)];
            if (chart == nullptr) {
                continue;
            }

            chart->grow(geom_[static_cast<std::size_t>(fopp)].normal);
            chart->faces.push_back(fopp);
            chartMap_[static_cast<std::size_t>(fopp)] = chart;
            lastFace_ = fopp;

            boundary_.erase(h);

            // Push the two new boundary half-edges from fopp's remaining edges.
            decimate::HeIndex loop = mesh_.next(s);
            for (int i = 0; i < 2; ++i) {
                if (mesh_.sym(loop) != decimate::kNoHalfedge &&
                    boundary_.isBoundary(loop)) {
                    heap.push({edgeCost(loop), static_cast<int>(loop)});
                }
                loop = mesh_.next(loop);
            }
            // RepairHeap commented out in legacy — stale entries remain.
        }
    }

    // -----------------------------------------------------------------------
    // initializeSeedSearch — port of InitializeSeedSearch (:150-176)
    // -----------------------------------------------------------------------
    MinHeap initializeSeedSearch(std::vector<bool>& visited) {
        const int nhe = static_cast<int>(mesh_.halfEdgeCount());
        std::fill(visited.begin(), visited.end(), false);

        MinHeap heap;
        for (int i = 0; i < nhe; ++i) {
            const decimate::HeIndex h = static_cast<decimate::HeIndex>(i);
            if (!boundary_.isBoundary(h)) {
                continue;
            }
            const int f      = static_cast<int>(mesh_.halfedge(h).face);
            Chart*    chart  = chartMap_[static_cast<std::size_t>(f)];
            if (chart == nullptr) {
                continue;
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
                const int loopFace    = static_cast<int>(mesh_.halfedge(loop).face);
                const int loopSymFace = static_cast<int>(mesh_.halfedge(ls).face);
                if (!visited[static_cast<std::size_t>(loopSymFace)] &&
                    chartMap_[static_cast<std::size_t>(loopFace)] ==
                    chartMap_[static_cast<std::size_t>(loopSymFace)]) {
                    heap.push({centroidDistance(loop), static_cast<int>(loop)});
                }
                loop = mesh_.next(loop);
            }
        }
        return heap;
    }

    // -----------------------------------------------------------------------
    // applySeeds — port of ApplySeeds (:185-204)
    // -----------------------------------------------------------------------
    bool applySeeds() {
        bool cycleDetected = true;
        for (auto& cptr : charts_) {
            Chart* chart = cptr.get();
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

    // -----------------------------------------------------------------------
    // computeSeeds — port of ComputeSeeds (:205-247)
    // -----------------------------------------------------------------------
    bool computeSeeds() {
        std::vector<bool> visited(static_cast<std::size_t>(numFaces_), false);
        MinHeap heap = initializeSeedSearch(visited);

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
                continue;
            }
            visited[static_cast<std::size_t>(fopp)] = true;

            Chart* chart = chartMap_[static_cast<std::size_t>(fopp)];
            if (chart == nullptr) {
                continue;
            }
            chart->searchIndex = fopp;

            decimate::HeIndex loop = mesh_.next(s);
            for (int i = 0; i < 2; ++i) {
                const decimate::HeIndex ls = mesh_.sym(loop);
                if (ls == decimate::kNoHalfedge) {
                    loop = mesh_.next(loop);
                    continue;
                }
                const int loopFace    = static_cast<int>(mesh_.halfedge(loop).face);
                const int loopSymFace = static_cast<int>(mesh_.halfedge(ls).face);
                if (!visited[static_cast<std::size_t>(loopSymFace)] &&
                    chartMap_[static_cast<std::size_t>(loopFace)] ==
                    chartMap_[static_cast<std::size_t>(loopSymFace)]) {
                    const double dist =
                        static_cast<double>(
                            (geom_[static_cast<std::size_t>(loopSymFace)].centroid -
                             geom_[static_cast<std::size_t>(loopFace)].centroid).length());
                    heap.push({value + dist, static_cast<int>(loop)});
                }
                loop = mesh_.next(loop);
            }
        }

        // Small charts: keep current seed (mirrors legacy clusterCount < 4 check).
        for (auto& cptr : charts_) {
            Chart* chart = cptr.get();
            if (chart->clusterCount < 4) {
                chart->searchIndex = chart->seed;
            }
        }

        return applySeeds();
    }

    // -----------------------------------------------------------------------
    // removeInvalidCharts — port of CChartSet::RemoveInvalidCharts (:71-91)
    // -----------------------------------------------------------------------
    /// Remove charts with clusterCount < 2 (fewer than 2 grow calls).
    /// Returns true if any were removed.
    bool removeInvalidCharts() {
        const std::size_t before = charts_.size();
        charts_.erase(
            std::remove_if(charts_.begin(), charts_.end(),
                [](const std::unique_ptr<Chart>& c) {
                    return c->clusterCount < 2;
                }),
            charts_.end());
        // Re-number ids to remain contiguous after removal.
        for (std::size_t i = 0; i < charts_.size(); ++i) {
            charts_[i]->id = static_cast<int>(i);
        }
        return charts_.size() < before;
    }

    // -----------------------------------------------------------------------
    // normalDeviation — cosine error between two face normals
    // -----------------------------------------------------------------------
    [[nodiscard]] double normalDeviation(int faceA, int faceB) const {
        const auto& na = geom_[static_cast<std::size_t>(faceA)].normal;
        const auto& nb = geom_[static_cast<std::size_t>(faceB)].normal;
        const double la = static_cast<double>(na.lengthSquared());
        const double lb = static_cast<double>(nb.lengthSquared());
        if (la <= 0.0 || lb <= 0.0) {
            return 1.0;
        }
        return 1.0 - static_cast<double>(na.dot(nb)) / std::sqrt(la * lb);
    }

    // -----------------------------------------------------------------------
    // tryStraigthenEdge — attempt to reassign face i via half-edge h
    // -----------------------------------------------------------------------
    /// Returns true if the face was reassigned (boundary modified).
    [[nodiscard]] bool tryStraightenEdge(int i, decimate::HeIndex h) {
        const decimate::HeIndex hnxt  = mesh_.next(h);
        if (!(boundary_.isBoundary(h) && boundary_.isBoundary(hnxt))) {
            return false;
        }
        const decimate::HeIndex hs   = mesh_.sym(h);
        const decimate::HeIndex hnxs = mesh_.sym(hnxt);
        if (hs == decimate::kNoHalfedge || hnxs == decimate::kNoHalfedge) {
            return false;
        }
        Chart* ochart1 = chartMap_[static_cast<std::size_t>(mesh_.halfedge(hs).face)];
        Chart* ochart2 = chartMap_[static_cast<std::size_t>(mesh_.halfedge(hnxs).face)];
        if (ochart1 == nullptr || ochart2 == nullptr) {
            return false;
        }
        if (ochart1 != ochart2 || chartMap_[static_cast<std::size_t>(i)] == ochart1) {
            return false;
        }
        // Threshold 0.3: mirrors legacy Straighten exactly.
        if (normalDeviation(i, static_cast<int>(mesh_.halfedge(hs).face)) < 0.3) {
            chartMap_[static_cast<std::size_t>(i)] = ochart1;
            boundary_.setFlag(mesh_.next(mesh_.next(h)), true);
            boundary_.erase(h);
            return true;
        }
        return false;
    }

    // -----------------------------------------------------------------------
    // straighten — port of Straighten (:505-558)
    // -----------------------------------------------------------------------
    /// Boundary-smoothing pass: reassign a face with two consecutive boundary
    /// edges on the same opposite chart iff normalError < 0.3.
    bool straighten() {
        bool changed = false;
        for (int i = 0; i < numFaces_; ++i) {
            const Chart* chart = chartMap_[static_cast<std::size_t>(i)];
            if (chart == nullptr || i == chart->seed) {
                continue;
            }
            for (int j = 0; j < 3; ++j) {
                if (tryStraightenEdge(i, static_cast<decimate::HeIndex>(i * 3 + j))) {
                    changed = true;
                    break;
                }
            }
        }
        return changed;
    }

    // -----------------------------------------------------------------------
    // bestNeighbourChart — lowest-id assigned chart among face f's neighbours
    // -----------------------------------------------------------------------
    [[nodiscard]] Chart* bestNeighbourChart(int f) const {
        Chart* best = nullptr;
        for (int e = 0; e < 3; ++e) {
            const decimate::HeIndex s = mesh_.sym(
                static_cast<decimate::HeIndex>(f * 3 + e));
            if (s == decimate::kNoHalfedge) {
                continue;
            }
            Chart* nc = chartMap_[static_cast<std::size_t>(mesh_.halfedge(s).face)];
            if (nc != nullptr && (best == nullptr || nc->id < best->id)) {
                best = nc;
            }
        }
        return best;
    }

    // -----------------------------------------------------------------------
    // orphanAssignment — ensure every face has a chart assignment
    // -----------------------------------------------------------------------
    /// Any face left unassigned (nullptr in chartMap_) gets the lowest-id
    /// chart among its assigned neighbours. Mirrors CChartSet orphan pass.
    void orphanAssignment() {
        bool progress = true;
        while (progress) {
            progress = false;
            for (int f = 0; f < numFaces_; ++f) {
                if (chartMap_[static_cast<std::size_t>(f)] != nullptr) {
                    continue;
                }
                Chart* best = bestNeighbourChart(f);
                if (best != nullptr) {
                    chartMap_[static_cast<std::size_t>(f)] = best;
                    progress = true;
                }
            }
        }
        // Fallback: isolated faces with no assigned neighbours get chart 0.
        if (!charts_.empty()) {
            Chart* fallback = charts_[0].get();
            for (int f = 0; f < numFaces_; ++f) {
                if (chartMap_[static_cast<std::size_t>(f)] == nullptr) {
                    chartMap_[static_cast<std::size_t>(f)] = fallback;
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // assemble — build the final ChartSet from chartMap_
    // -----------------------------------------------------------------------
    ChartSet assemble() {
        ChartSet result;
        result.charts.reserve(charts_.size());

        for (auto& cptr : charts_) {
            result.charts.push_back(*cptr);
        }

        // Single source of truth: derive both Chart.faces and faceToChart from
        // chartMap_ (the authoritative face->chart assignment). The Chart.faces
        // built incrementally during grow goes stale when straighten() and
        // orphanAssignment() reassign chartMap_ without touching it; rebuilding
        // here keeps the two representations consistent by construction.
        // (chart->id == index in result.charts — see createChart / compaction.)
        for (Chart& ch : result.charts) {
            ch.faces.clear();
        }
        result.faceToChart.resize(static_cast<std::size_t>(numFaces_), -1);
        for (int f = 0; f < numFaces_; ++f) {
            const Chart* chart = chartMap_[static_cast<std::size_t>(f)];
            if (chart != nullptr) {
                result.faceToChart[static_cast<std::size_t>(f)] = chart->id;
                result.charts[static_cast<std::size_t>(chart->id)].faces.push_back(f);
            }
        }

        return result;
    }
};

// ---------------------------------------------------------------------------
// Segmenter::run — mirrors legacy CChartification::Chartify (:74-114)
// ---------------------------------------------------------------------------
ChartSet Segmenter::run() {  // NOLINT(readability-function-size): mirrors legacy Chartify (:74-114) — direct port
    if (numFaces_ <= 0) {
        return {};
    }

    // Bootstrap: start with face 0 as the first chart (mirrors legacy exactly).
    createChart(0);

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
    int iter = 0;
    do {
        expandCharts();
        if (static_cast<int>(charts_.size()) < maxCharts_) {
            // Check that lastFace_ is not already a seed for any chart.
            bool alreadySeed = false;
            for (const auto& cptr : charts_) {
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

    // Remove charts with clusterCount < 2; re-grow until stable.
    while (removeInvalidCharts()) {
        if (charts_.empty()) {
            break;
        }
        expandCharts();
    }

    // Final ExpandCharts after invalid removal (mirrors legacy pattern).
    if (!charts_.empty()) {
        expandCharts();
    }

    // If all charts were removed (degenerate small mesh), create a fallback
    // chart covering all faces via orphan assignment.
    if (charts_.empty()) {
        createChart(0);
        // Mark all faces as belonging to chart 0 directly.
        for (int f = 0; f < numFaces_; ++f) {
            chartMap_[static_cast<std::size_t>(f)] = charts_[0].get();
        }
    }

    // Boundary straightening.
    while (straighten()) {}

    // Orphan assignment: fill any unassigned faces.
    orphanAssignment();

    return assemble();
}

}  // namespace (anonymous)

// ===========================================================================
// Public entry point
// ===========================================================================

ChartSet segment(const decimate::HalfEdgeMesh& mesh,
                  const decimate::MeshData&     data,
                  int                        maxCharts) {
    Segmenter seg(mesh, data, maxCharts);
    return seg.run();
}

}  // namespace decimate::chart
