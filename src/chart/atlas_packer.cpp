// This file is part of the decimate library; see decimate.h for version/license details.
/// @file atlas_packer.cpp
/// @brief Atlas packing + full parametrize() orchestrator — Spec 08-02, Beat 2a.
///
/// ## Behavior-preservation notes
///
/// ### ScaleProportional (legacy catlaspacker.cpp ~398-427)
///
/// Legacy formula per chart:
/// ```cpp
/// float c_allArea  = charts->Area(points);          // sum of all 3D face areas
/// float chartArea  = chart->Area(points);            // this chart's 3D area
/// float textureArea = TextureArea(chart);            // this chart's UV area
/// float cscale = sqrtf(chartArea / (c_allArea * textureArea));
/// vertices[i].u *= cscale;
/// vertices[i].v *= cscale;
/// ```
/// This normalises each chart so that post-scale UV area == chartArea / totalArea3D.
/// After applying it to every chart the sum of UV areas equals 1.0.
///
/// ### Sort order
///
/// Legacy `Cmp_Charts` sorts descending by `FaceCount()` using `std::list::sort`
/// (stable merge sort, preserves insertion order for equal-count charts).
/// We use `std::stable_sort` with descending face-count and ascending-id
/// tie-break.  Note: the spec says "area desc" but the legacy criterion is face
/// count.  Beat 2a preserves the legacy face-count criterion exactly.
///
/// ### Pack loop
///
/// Legacy tries a descending `ratio` (start 4.0, step -0.01) until all charts
/// fit into the grid.  The `ratio` controls the UV-to-pixel scale factor
/// (`scale = ratio / sqrt(1)`).  We preserve this outer loop exactly.
/// `atlasResolution` (default 256) replaces the legacy hardcoded `discreetSize`.
///
/// ### 16-rotation search
///
/// For each rotation step (0°, 22.5°, …, 337.5°) the chart boundary UVs are
/// discretised into the `upper` / `lower` column arrays, then `FindOptimal`
/// scans for the placement with least wasted horizon area.  The placement
/// with the globally smallest waste across all 16 rotations wins.
///
/// Triangle vertices for the silhouette footprint are derived from the local
/// `HalfEdgeMesh` in `buildSilhouetteFootprint`.  The chart is a valid 2D
/// embedding (LSCM with no foldovers — measured: 0 flipped triangles, 0 holes
/// across all spot charts), so the per-column min/max silhouette of the boundary
/// loop equals that of all triangle edges (measured: identical on every chart).
/// All-edges is rasterised here only because the boundary-loop walk becomes
/// load-bearing in Spec 09-02 (the miter gutter offsets that loop); both yield
/// the same footprint, gated by the zero-overlap SAT test.
///
/// ### parametrize() UV scatter
///
/// Per-mesh-vertex UVs are written chart by chart in order (chart 0, 1, …).
/// A seam vertex shared by multiple charts is overwritten each time, ending
/// with the LAST chart's UV — matching legacy `CConformalMapper::Parametrize`.
///
/// ### Boundary contract
/// No GL, no OS, no I/O, no global mutable state.

#include <decimate/chart/parametrize.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <vector>
#include <span>

#include <decimate/math/vec2.h>
#include <decimate/math/vec3.h>
#include <decimate/chart/chart.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>

namespace decimate::chart {

// ===========================================================================
// Internal UV type
// ===========================================================================

// ===========================================================================
// Horizon — discrete column-based atlas horizon (port of legacy CHorizon)
// ===========================================================================

/// @brief Tracks the current packed height at each column of the atlas grid.
///
/// `data_[i]` is the height at column `i` (0 = empty column).  After each
/// chart is placed, its upper silhouette is merged into the horizon via
/// `update()`.
class Horizon {
public:
    explicit Horizon(int resolution) : res_(resolution), data_(resolution, 0) {}

    /// @brief Reset all columns to zero (empty atlas).
    void reset() { std::fill(data_.begin(), data_.end(), 0); }

    /// @brief Height at column `index`.
    int operator[](int index) const noexcept { return data_[index]; }

    /// @brief Total discrete area of the atlas (res²).
    [[nodiscard]] int totalArea() const noexcept { return res_ * res_; }

    /// @brief Wasted area if a chart with lower silhouette `lower` is placed at (x, y).
    ///
    /// `lower[i]` is the per-column minimum v of the chart (its bottom silhouette).
    /// Returns the sum of pixel gaps between the horizon top and the chart bottom,
    /// or -1 if any column of the chart falls outside the atlas boundary or
    /// collides with already-packed content.
    ///
    /// Port of `CHorizon::Area(int x, int y, int size, int *chart)` where `chart`
    /// is the `lower` (bottom silhouette) array.
    [[nodiscard]] int wastedArea(int x, int y, int size,  // NOLINT(bugprone-easily-swappable-parameters): (x, y, size) are distinct grid coordinates
                                 const std::vector<int>& lower) const {
        int wasted = 0;
        for (int i = 0; i < size; ++i) {
            int col    = x + i;
            if (col >= res_) { return -1; }             // past atlas edge
            int bottom = y + lower[i];                  // chart bottom at this column
            if (data_[col] >= bottom) { return -1; }   // horizon already above chart bottom
            wasted += bottom - data_[col];              // gap between horizon and chart bottom
        }
        return wasted;
    }

    /// @brief Merge a chart's upper silhouette into the horizon.
    ///
    /// Port of `CHorizon::Update(int xoffs, int yoffs, int size, int *chart)`.
    /// Gutter is NOT applied here — it lives solely in `findOptimal` (applied to
    /// the landing y), so the placed chart already sits `gutter` above the
    /// skyline.  Adding it again here would double-count the separation (it did:
    /// it cost ~5% density).  We write the crown at its true height `upper[i] +
    /// yoffs`.
    void update(int xoffs, int yoffs, int size,  // NOLINT(bugprone-easily-swappable-parameters): column-offset, row-offset, span-width
                const std::vector<int>& upper) {
        for (int i = 0; i < size; ++i) {
            int col = xoffs + i;
            if (col < res_) {
                data_[col] = upper[i] + yoffs;  // gutter already in yoffs (findOptimal)
            }
        }
    }

    [[nodiscard]] int resolution() const noexcept { return res_; }

private:
    int              res_;
    std::vector<int> data_;
};

// ===========================================================================
// Geometry helpers
// ===========================================================================

/// @brief Half-cross-product magnitude — 3D area of a triangle.
static float triArea3D(const decimate::math::Vec3<float>& a,
                       const decimate::math::Vec3<float>& b,
                       const decimate::math::Vec3<float>& c) {
    return (b - a).cross(c - a).length() * 0.5f;
}

/// @brief 3D surface area of an HalfEdgeChart summed over its face list.
static float chartArea3D(const HalfEdgeChart& ac, const MeshData& data) {
    float total = 0.0f;
    const auto& p = data.positions;
    for (int fi : ac.faces) {
        const decimate::math::Vec3<float>& a = p[static_cast<std::size_t>(data.indices[3 * fi])];
        const decimate::math::Vec3<float>& b = p[static_cast<std::size_t>(data.indices[3 * fi + 1])];
        const decimate::math::Vec3<float>& c = p[static_cast<std::size_t>(data.indices[3 * fi + 2])];
        total += triArea3D(a, b, c);
    }
    return total;
}

/// @brief UV (shoelace) area of a chart, using a globalToLocal index map.
///
/// Port of `CAtlasPacker::TextureArea(CChart *)`.
///
/// @param ac           Atlas chart (uv array indexed by local vertex id).
/// @param data         Mesh face indices.
/// @param globalToLocal  Map[globalVertexId] → local vertex id (-1 if absent).
static float chartAreaUV(const HalfEdgeChart&       ac,
                         const MeshData&         data,
                         const std::vector<int>& globalToLocal) {
    float area = 0.0f;
    for (int fi : ac.faces) {
        int ga = data.indices[3 * fi];
        int gb = data.indices[3 * fi + 1];
        int gc = data.indices[3 * fi + 2];
        int la = globalToLocal[ga];
        int lb = globalToLocal[gb];
        int lc = globalToLocal[gc];
        if (la < 0 || lb < 0 || lc < 0) { continue; }
        auto& a = ac.uv[la];
        auto& b = ac.uv[lb];
        auto& c = ac.uv[lc];
        // Shoelace determinant (matches legacy TextureArea formula).
        float det = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
        area += std::abs(det) * 0.5f;
    }
    return area;
}

// ===========================================================================
// ScaleProportional
// ===========================================================================

/// @brief Scale each chart's UVs so that UV area ratios match 3D area ratios.
///
/// Port of `CAtlasPacker::ScaleProportional(CChartSet *)`.
/// After scaling, each chart's UV area equals `chartArea3D / totalArea3D`.
/// A chart with zero 3D or UV area is skipped (matches legacy `IsFZero` guard).
static void scaleProportional(std::vector<HalfEdgeChart>& charts, const MeshData& data) {
    // 1. Total 3D surface area across all charts.
    float totalArea3D = 0.0f;
    for (const HalfEdgeChart& ac : charts) {
        totalArea3D += chartArea3D(ac, data);
    }
    if (totalArea3D < 1e-12f) { return; }  // degenerate — nothing to scale

    // 2. Per-chart scaling.
    const int numVerts = static_cast<int>(data.positions.size());
    std::vector<int> globalToLocal(numVerts, -1);

    for (HalfEdgeChart& ac : charts) {
        // Build globalToLocal for this chart.
        const int numLocal = static_cast<int>(ac.localToMesh.size());
        for (int li = 0; li < numLocal; ++li) {
            globalToLocal[ac.localToMesh[li]] = li;
        }

        float area3D = chartArea3D(ac, data);
        float areaUV = chartAreaUV(ac, data, globalToLocal);

        // Reset globalToLocal (only the entries touched by this chart).
        for (int gv : ac.localToMesh) {
            globalToLocal[gv] = -1;
        }

        if (areaUV < 1e-12f || area3D < 1e-12f) { continue; }

        // Legacy formula: cscale = sqrt(chartArea / (totalArea * textureArea)).
        float cscale = std::sqrt(area3D / (totalArea3D * areaUV));

        for (decimate::math::Vec2<float>& uv : ac.uv) {
            uv.x *= cscale;
            uv.y *= cscale;
        }
    }
}

// ===========================================================================
// Tex2 AABB helpers
// ===========================================================================

static Tex2 minTex(const std::vector<Tex2>& pts) {
    Tex2 m = pts[0];
    for (const Tex2& p : pts) {
        m.u = std::min(m.u, p.u);
        m.v = std::min(m.v, p.v);
    }
    return m;
}

static Tex2 maxTex(const std::vector<Tex2>& pts) {
    Tex2 m = pts[0];
    for (const Tex2& p : pts) {
        m.u = std::max(m.u, p.u);
        m.v = std::max(m.v, p.v);
    }
    return m;
}

/// @brief Centre of the AABB spanned by `pts`.
///
/// Port of `CAtlasPacker::Mid(int size, Vertex *vertices)`.
/// Legacy computes `(max - min) / 2`, which is the half-range, not the centre.
/// We preserve that exact computation for rotation-centre consistency.
static Tex2 midTex(const std::vector<Tex2>& pts) {
    Tex2 mn = minTex(pts);
    Tex2 mx = maxTex(pts);
    return { (mx.u - mn.u) * 0.5f, (mx.v - mn.v) * 0.5f };
}

// ===========================================================================
// FindOptimal
// ===========================================================================

/// @brief Result of a single placement scan.
struct OptimalPlacement {
    int  x;       ///< Best column (−1 = not found).
    int  y;       ///< Best row (−1 = not found).
    int  waste;   ///< Wasted pixel area at (x, y); horizon.totalArea() if not found.
};

/// @brief Scan the atlas for the position of MINIMUM wasted shadow area.
///
/// This is the Tetris best-fit from the thesis (Greissl 2003, §Positionsbestimmung):
/// "Der optimale Einfügepunkt ist gefunden, wenn diese Fläche minimal ist" — the
/// optimal insertion point is where the discrete shadow area between the chart's
/// lower horizon and the atlas horizon is MINIMAL.  Minimising that area makes a
/// chart settle into the deepest matching notch of the current skyline, so charts
/// interlock densely.
///
/// ## TETRIS DIRECT-DROP, FUSED (experiment/tetris-findoptimal)
///
/// The `y` for a column is NOT searched — it falls out of the area sum for free.
/// Per placement column `x`, walk the chart width ONCE accumulating, over
/// chart-column `k` (skyline window `s = horizon[x+k]`, chart underside
/// `t = lower[k]`):
///   colGap = t − s            (= lower[k] − skyline[x+k], may be negative)
///   rawArea += colGap
///   minGap   = min(minGap, colGap)
///   maxUp    = max(maxUp, upper[k])     (chart crown, for the height-fit check)
///
/// The Tetris landing lifts the chart so its DEEPEST column just touches the
/// skyline: `y = −minGap` (+ gutter).  Because the shadow area is linear in y —
///   area(y) = Σ(y + lower[k] − skyline[x+k]) = rawArea + width·y
/// — the area at the landing is `rawArea − width·minGap` in the SAME pass; no
/// second walk, no collision branch (landing-to-touch can't penetrate).  So each
/// column is ONE width-pass: O(res · width), and y + area both come for free.
///
/// SoC: `findOptimal` is told nothing about how the silhouette was built — it
/// receives both arrays + the gutter and owns the drop AND the fit decision.
///
/// @param width   Chart width in grid columns.
/// @param lower   Chart lower silhouette (underside): per-column min v offset.
/// @param upper   Chart upper silhouette (crown): per-column max v offset — for
///                the height-fit check (`y + maxUpper + gutter` must fit res).
/// @param horizon Current atlas horizon.
/// @param gutter  Boundary in pixels applied to the landing y AND the fit check.
/// @returns Placement of minimum wasted area; x==−1 when no valid position fits.
static OptimalPlacement findOptimal(int width,
                                    const std::vector<int>& lower,
                                    const std::vector<int>& upper,
                                    const Horizon& horizon,
                                    int gutter) {
    const int res = horizon.resolution();
    OptimalPlacement best{-1, -1, horizon.totalArea()};

    for (int x = 0; x < res - width; ++x) {
        // ONE fused pass: area sum, deepest column (→ landing y), and crown.
        int rawArea = 0;
        int minGap  = std::numeric_limits<int>::max();
        int maxUp   = 0;
        for (int k = 0; k < width; ++k) {
            const int colGap = lower[static_cast<std::size_t>(k)] - horizon[x + k];
            rawArea += colGap;
            if (colGap < minGap) { minGap = colGap; }
            if (upper[static_cast<std::size_t>(k)] > maxUp) { maxUp = upper[static_cast<std::size_t>(k)]; }
        }

        const int y = -minGap + gutter;             // drop deepest column to touch + boundary
        if (y < 0 || y + maxUp + gutter >= res) {   // height fit (gutter inside, SoC)
            continue;
        }
        const int area = rawArea - width * minGap;  // area at the landing — free, same pass
        if (area < best.waste) {
            best = {x, y, area};  // deepest-notch (min-shadow) column so far
        }
    }
    return best;
}

// ===========================================================================
// packSingleChart — 16-rotation horizon fit
// ===========================================================================

/// @brief Chart footprint in pixel space for one rotation.
struct RotationFootprint {
    int width  = 0;  ///< AABB width in grid columns (0 = degenerate / skip).
    int height = 0;  ///< AABB height in grid rows    (0 = degenerate / skip).
};

/// @brief Rasterise one boundary edge into the per-column lower/upper silhouette.
///
/// Port of `CAtlasPacker::RasterLine` (2003).  For every column the edge crosses,
/// records the min v in `lower[u]` and max v in `upper[u]`.  After all the
/// chart's triangle edges are rasterised, `lower` is the chart's TRUE bottom
/// silhouette (the "shadow" underside) and `upper` its top — not a box.  One
/// sample per column via integer interpolation `v = v1 + (dv*(u-u1))/du`.
///
/// ## Alternative explored: Bresenham (Andor, 2025) — kept here for discussion
///
/// A per-pixel Bresenham DDA on `std::span` (add-only inner loop, single
/// unsigned bounds compare).  Originally swapped in because it made the atlas
/// pack tighter — but a fixture-driven benchmark (`tests/raster_bench.cpp`,
/// cow→1000, 43 charts) corrected the picture:
///   - SPEED: Bresenham is NOT faster here.  It steps once per PIXEL of edge
///     length; interp steps once per COLUMN.  On real (short) chart edges it
///     does ~1.3× more inner iterations; on long edges far more.  The per-column
///     divide is cheap/predictable on a modern CPU — never the bottleneck.
///   - HOT PATH: rasterLine is NOT the packer's cost at all.  `findOptimal`'s
///     O(res²) position scan dominates (packAtlas: res128 1.2s → res256 5.5s →
///     res512 27s ≈ O(res²)).  Packer speed work = skyline findOptimal (Spec 13).
///   - WHY THE ATLAS LOOKED BETTER: interp DROPS any edge with an out-of-bounds
///     endpoint (`if (u1<0||u1>=res||u2<0||u2>=res) return;`) — whole edge gone,
///     hole in the silhouette.  Bresenham clips per-pixel and keeps the in-bounds
///     span.  So the win was CORRECTNESS (OOB clipping), not the algorithm.
///   → Open question to discuss: the cheapest correct form is THIS interp loop
///     with the OOB endpoint-drop replaced by per-column CLAMP (keep the cheap
///     per-column cost, gain Bresenham's clipping).  Bresenham verbatim:
///
/// @code
/// static void rasterLine(int u1, int v1, int u2, int v2, int res,
///                        std::span<int> lower, std::span<int> upper) {
///     const int dx =  std::abs(u2 - u1);
///     const int sx =  u1 < u2 ? 1 : -1;
///     const int dy = -std::abs(v2 - v1);
///     const int sy =  v1 < v2 ? 1 : -1;
///     int err = dx + dy;
///     int u = u1, v = v1;
///     for (;;) {
///         if (static_cast<unsigned>(u) < static_cast<unsigned>(res)) {  // 0<=u<res, one compare
///             if (v > upper[u]) upper[u] = v;
///             if (v < lower[u]) lower[u] = v;
///         }
///         if (u == u2 && v == v2) { break; }
///         const int e2 = 2 * err;
///         if (e2 >= dy) { err += dy; u += sx; }
///         if (e2 <= dx) { err += dx; v += sy; }
///     }
/// }
/// @endcode
static void rasterLine(int u1, int v1, int u2, int v2, int res,
                       std::vector<int>& lower, std::vector<int>& upper) {
    if (u1 > u2) { std::swap(u1, u2); std::swap(v1, v2); }
    if (u1 < 0 || u1 >= res || u2 < 0 || u2 >= res) { return; }

    const int du = u2 - u1 + 1;
    const int dv = v2 - v1 + 1;
    if (du <= 0) { return; }

    for (int u = u1; u <= u2; ++u) {
        if (u < 0 || u >= res) {
            continue;
        }
        const int v = v1 + (dv * (u - u1)) / du;  // integer interpolation (legacy)
        if (upper[static_cast<std::size_t>(u)] < v) { upper[static_cast<std::size_t>(u)] = v; }
        if (lower[static_cast<std::size_t>(u)] > v) { lower[static_cast<std::size_t>(u)] = v; }
    }
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters): (gutter, miterLimit) signature fixed by the spec
std::vector<Tex2> miterOffsetLoop(const std::vector<Tex2>& loop, float gutter, float miterLimit) {
    const std::size_t n = loop.size();
    if (n < 3) {
        return loop;
    }

    // Winding (shoelace): interior is on the left of each edge when area > 0 (CCW).
    float area2 = 0.0f;
    for (std::size_t i = 0; i < n; ++i) {
        const Tex2& a = loop[i];
        const Tex2& b = loop[(i + 1) % n];
        area2 += a.u * b.v - b.u * a.v;
    }
    const float wind = (area2 >= 0.0f) ? 1.0f : -1.0f;
    // Outward unit normal of a directed edge (right-hand perp, winding-signed).
    auto outNormal = [wind](const Tex2& a, const Tex2& b, float& nx, float& ny) {
        const float dx  = b.u - a.u;
        const float dy  = b.v - a.v;
        const float len = std::sqrt(dx * dx + dy * dy);
        if (len < 1e-6f) {
            nx = 0.0f;
            ny = 0.0f;
            return;
        }
        nx = wind * dy / len;
        ny = -wind * dx / len;
    };

    // Move along the bisector of the two edge normals, length gutter/cos(half-angle),
    // clamped at miterLimit x.  Outward-only => the silhouette only grows.
    const float minCos = 1.0f / std::max(miterLimit, 1.0f);
    std::vector<Tex2> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        const Tex2& c = loop[i];
        float n1x = 0.0f;
        float n1y = 0.0f;
        float n2x = 0.0f;
        float n2y = 0.0f;
        outNormal(loop[(i + n - 1) % n], c, n1x, n1y);   // incoming edge
        outNormal(c, loop[(i + 1) % n], n2x, n2y);       // outgoing edge
        float mx = n1x + n2x;
        float my = n1y + n2y;
        const float ml = std::sqrt(mx * mx + my * my);
        if (ml < 1e-6f) {
            mx = n1x;
            my = n1y;
        } else {
            mx /= ml;
            my /= ml;
        }
        const float cosHalf = std::max(mx * n1x + my * n1y, minCos);
        const float d       = gutter / cosHalf;
        out[i] = {c.u + mx * d, c.v + my * d};
    }
    return out;
}

/// @brief Chart boundary loop in padded pixel space for one rotation.
///
/// Pixel = (uv - mn) * pxScale + kGutterPx: the +kGutterPx pad keeps the outward
/// miter offset in-grid on all four sides; applyBestPlacement adds the matching
/// +kGutterPx so the chart sits centred inside its padded footprint.
/// Empty when the chart has fewer than 3 boundary vertices.
static std::vector<Tex2> boundaryLoopPx(const std::vector<Tex2>& dst,
                                        const HalfEdgeMesh&      localMesh,
                                        const Tex2&              mn,
                                        float                    pxScale) {
    const std::vector<std::uint32_t> loop = localMesh.boundaryPointIndices();
    if (loop.size() < 3) {
        return {};
    }
    const auto pad = static_cast<float>(kGutterPx);
    std::vector<Tex2> px(loop.size());
    for (std::size_t i = 0; i < loop.size(); ++i) {
        const Tex2& p = dst[static_cast<std::size_t>(loop[i])];
        px[i] = {(p.u - mn.u) * pxScale + pad, (p.v - mn.v) * pxScale + pad};
    }
    return px;
}

/// @brief Build the chart's per-column silhouette footprint for one rotation.
///
/// Rasterises the MITER-OFFSET boundary loop into `lower`/`upper` via `rasterLine`.
/// The loop is the half-edges of the cut mesh with no twin — the faithful 2003
/// silhouette (`legacy2003/chartification/cchartification.cpp:395`,
/// `legacy2003/split/splitter.cpp:1600`).  The chart is a valid 2D embedding, so
/// the boundary bounds the interior.  The offset itself is `miterOffsetLoop`.
static RotationFootprint buildSilhouetteFootprint(const std::vector<Tex2>&   dst,
                                                  const HalfEdgeMesh&        localMesh,
                                                  const Tex2&                mn,
                                                  float                      scale,
                                                  int                        res,
                                                  std::vector<int>&          lower,
                                                  std::vector<int>&          upper) {
    for (int col = 0; col < res; ++col) {
        upper[static_cast<std::size_t>(col)] = 0;
        lower[static_cast<std::size_t>(col)] = res - 1;
    }

    const float pxScale = scale * static_cast<float>(res - 1);
    const std::vector<Tex2> loop = boundaryLoopPx(dst, localMesh, mn, pxScale);
    if (loop.empty()) {
        return {0, 0};
    }
    const std::vector<Tex2> off = miterOffsetLoop(loop, static_cast<float>(kGutterPx));

    const std::size_t n = off.size();
    int width = 0;
    std::vector<int> qx(n);
    std::vector<int> qy(n);
    const auto clampPx = [res](float a) {
        return std::clamp(static_cast<int>(a), 0, res - 1);
    };
    for (std::size_t i = 0; i < n; ++i) {
        qx[i] = clampPx(off[i].u);
        qy[i] = clampPx(off[i].v);
        width = std::max(width, qx[i] + 1);
    }
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = (i + 1) % n;
        rasterLine(qx[i], qy[i], qx[j], qy[j], res, lower, upper);
    }

    if (width <= 0 || width > res) { return {0, 0}; }

    // Height = ymax(upper) − ymin(lower) over the spanned columns.
    int ymin = res;
    int ymax = 0;
    for (int col = 0; col < width; ++col) {
        ymin = std::min(ymin, lower[static_cast<std::size_t>(col)]);
        ymax = std::max(ymax, upper[static_cast<std::size_t>(col)]);
    }
    const int height = ymax - ymin;
    if (height <= 0 || height >= res) { return {0, 0}; }
    return {width, height};
}

/// @brief Winning placement result for the 16-rotation search.
struct BestRotation {
    float angleDeg = 0.0f;  ///< Rotation angle that minimises wasted area.
    int   x        = -1;    ///< Grid column of the best placement (−1 = none found).
    int   y        = -1;    ///< Grid row   of the best placement.
    int   width    = 0;     ///< Chart width in grid columns at this rotation.
};

/// @brief Apply the winning rotation + translation to `ac.uv` (write pass).
///
/// Re-rotates `src` by `best.angleDeg` around `mid`, shifts so the chart's
/// pixel-space minimum lands at `(best.x + 0.5, best.y + 0.5)` in [0,1]²
/// UV, and stores the result in `ac.uv`.
///
/// @param ac    Chart whose uv array is updated in place.
/// @param src   Pre-rotation UV vertices (same layout as ac.uv before packing).
/// @param mid   Rotation centre (AABB half-range, from midTex).
/// @param best  Winning placement from the 16-rotation search.
/// @param scale UV-to-pixel scale factor.
/// @param res   Atlas grid resolution.
static void applyBestPlacement(HalfEdgeChart&              ac,
                               const std::vector<Tex2>& src,
                               const Tex2&              mid,
                               const BestRotation&      best,
                               float                    scale,
                               int                      res) {
    const int   numUV = static_cast<int>(src.size());
    const float rad   = best.angleDeg * (std::numbers::pi_v<float> / 180.0f);
    const float sina  = std::sin(rad);
    const float cosa  = std::cos(rad);

    // Re-rotate src around mid.
    std::vector<Tex2> dst(numUV);
    for (int i = 0; i < numUV; ++i) {
        const float u = src[i].u - mid.u;
        const float v = src[i].v - mid.v;
        dst[i].u = mid.u + u * cosa - v * sina;
        dst[i].v = mid.v + u * sina + v * cosa;
    }

    const Tex2  mn     = minTex(dst);
    const float invRes = 1.0f / static_cast<float>(res - 1);
    // Half-pixel offset: matches legacy `(best_x + 0.5f) / (discreetSize - 1)`.
    // +kGutterPx: the footprint was padded by kGutterPx (buildSilhouetteFootprint),
    // so the chart's own min sits kGutterPx inside its placed footprint — the gutter.
    const float offX   = (static_cast<float>(best.x + kGutterPx) + 0.5f) * invRes;
    const float offY   = (static_cast<float>(best.y + kGutterPx) + 0.5f) * invRes;

    for (int i = 0; i < numUV; ++i) {
        ac.uv[i].x = (dst[i].u - mn.u) * scale + offX;
        ac.uv[i].y = (dst[i].v - mn.v) * scale + offY;
    }

    // Placed gutter outline: the miter-offset loop is in padded pixel space
    // q = (dst-mn)*pxScale + kGutterPx, so uv = (q - kGutterPx)*invRes + offX
    //                                         = (q + best.x + 0.5)*invRes.
    const float pxScale = scale * static_cast<float>(res - 1);
    const std::vector<Tex2> off = miterOffsetLoop(
        boundaryLoopPx(dst, ac.mesh, mn, pxScale), static_cast<float>(kGutterPx));
    const auto bx = static_cast<float>(best.x) + 0.5f;
    const auto by = static_cast<float>(best.y) + 0.5f;
    ac.gutterLoop.resize(off.size());
    for (std::size_t i = 0; i < off.size(); ++i) {
        ac.gutterLoop[i] = {(off[i].u + bx) * invRes, (off[i].v + by) * invRes};
    }
}

/// @brief Run the 16-rotation search and return the best placement found.
///
/// Tries each of the 16 evenly-spaced rotation angles (0°, 22.5°, …, 337.5°),
/// computes the chart's pixel-space AABB footprint at each rotation, and
/// searches the horizon for the placement with minimum wasted area.
///
/// @param src     Chart UV vertices in pre-pack UV coordinates.
/// @param mid     Rotation centre (AABB half-range, from midTex).
/// @param scale   UV-to-pixel scale factor.
/// @param res     Atlas grid resolution.
/// @param horizon Current atlas horizon (read-only during the search).
/// @param lower   Scratch column buffer (written by buildRectFootprint).
/// @brief Scratch column buffers shared across all rotation steps.
///
/// `lower` and `upper` are reused per rotation (reset in `buildRectFootprint`).
/// `best` captures the winning rotation's upper silhouette for the horizon update.
struct ScratchBuffers {
    std::vector<int> lower;   ///< Per-column minimum v (reset each rotation).
    std::vector<int> upper;   ///< Per-column maximum v (reset each rotation).
    std::vector<int> best;    ///< Upper silhouette of the best rotation found.
};

/// @param[out] scratch  Reusable column buffers; `scratch.best` holds result.
/// @returns Best rotation + placement found; BestRotation::x == -1 when none.
static BestRotation findBestRotation(const std::vector<Tex2>& src,
                                     const HalfEdgeMesh&      localMesh,
                                     const Tex2&              mid,
                                     float                    scale,
                                     int                      res,
                                     Horizon&                 horizon,
                                     ScratchBuffers&          scratch) {
    constexpr int   kSteps    = 16;
    constexpr float kAngleInc = 360.0f / static_cast<float>(kSteps);

    // Compare rotations by shadow area PER COLUMN (waste / width), not raw
    // shadow.  A tall-narrow orientation spans fewer columns and so accrues a
    // smaller TOTAL shadow almost by construction — biasing the search toward
    // vertical charts regardless of how well they actually nest.  Normalising
    // by width measures the average gap a column leaves under the chart, which
    // is the orientation-fair "how snugly does this underside fit the skyline"
    // criterion (Greissl, original design intent).  Width is constant across
    // all x-positions of a fixed rotation, so findOptimal's per-position argmin
    // is unaffected — only the rotation choice changes.
    float        bestPerCol = std::numeric_limits<float>::max();
    BestRotation best;
    std::vector<Tex2> dst(src.size());

    for (int step = 0; step < kSteps; ++step) {
        const float angleDeg = static_cast<float>(step) * kAngleInc;
        const float rad      = angleDeg * (std::numbers::pi_v<float> / 180.0f);
        const float sina     = std::sin(rad);
        const float cosa     = std::cos(rad);

        for (int i = 0; i < static_cast<int>(src.size()); ++i) {
            const float u = src[i].u - mid.u;
            const float v = src[i].v - mid.v;
            dst[i].u = mid.u + u * cosa - v * sina;
            dst[i].v = mid.v + u * sina + v * cosa;
        }

        const Tex2             mn = minTex(dst);
        const RotationFootprint fp =
            buildSilhouetteFootprint(dst, localMesh, mn, scale, res,
                                     scratch.lower, scratch.upper);
        if (fp.width <= 0) { continue; }

        // Vertical-only gutter is now 0: the all-around miter offset
        // (buildSilhouetteFootprint, kGutterPx) provides the chart separation on
        // every side, so adding the legacy vertical +2 here would double-count it.
        const OptimalPlacement pl =
            findOptimal(fp.width, scratch.lower, scratch.upper, horizon, 0);
        if (pl.x == -1) { continue; }

        const float perCol = static_cast<float>(pl.waste) / static_cast<float>(fp.width);
        if (perCol < bestPerCol) {
            bestPerCol = perCol;
            best       = {angleDeg, pl.x, pl.y, fp.width};
            scratch.best.assign(scratch.upper.begin(),
                                scratch.upper.begin() + fp.width);
        }
    }
    return best;
}

/// @brief Try 16 rotations of `ac` and pick the best placement in `horizon`.
///
/// Port of `CAtlasPacker::Pack(CChart*, float scale, bool write)`.
/// Each rotation's chart AABB is treated as a rectangle for the horizon
/// collision check — see `buildRectFootprint` for the rationale.
///
/// @param ac      Chart — `uv` mutated in-place when `write` is true.
/// @param scale   UV-to-pixel scale (pixels per UV unit).
/// @param write   When true, apply the best rotation + translation to `ac.uv`.
/// @param horizon Atlas horizon (updated to include the placed chart).
/// @param res     Atlas grid resolution.
/// @param scratch Reusable column buffers (shared across all charts in a pass).
/// @returns true if a valid placement was found.
static bool packSingleChart(HalfEdgeChart&     ac,
                            float           scale,
                            bool            write,
                            Horizon&        horizon,
                            int             res,
                            ScratchBuffers& scratch) {
    const int numUV = static_cast<int>(ac.uv.size());
    if (numUV == 0) { return true; }  // degenerate chart

    // Copy ac.uv to src Tex2 array and compute the rotation centre.
    std::vector<Tex2> src(numUV);
    for (int i = 0; i < numUV; ++i) {
        src[i].u = ac.uv[i].x;
        src[i].v = ac.uv[i].y;
    }
    const Tex2 mid = midTex(src);

    // Run the 16-rotation search (silhouette footprint + width-fair cost).
    // The local mesh is passed; buildSilhouetteFootprint iterates its faces
    // to rasterise all triangle edges — equivalent to the old localTri path.
    const BestRotation best =
        findBestRotation(src, ac.mesh, mid, scale, res, horizon, scratch);

    // Update horizon with the winning rotation's silhouette.
    if (best.x != -1) {
        horizon.update(best.x, best.y, best.width, scratch.best);
    }

    // Write best placement UVs back to the chart (only in write pass).
    if (write && best.x != -1) {
        applyBestPlacement(ac, src, mid, best, scale, res);
    }

    return (best.x != -1);
}

// ===========================================================================
// Normalize: map the packed layout into [0,1]²
// ===========================================================================

/// @brief Translate + uniformly scale all chart UVs so the aggregate AABB is [0,1]².
static void normalizeAtlas(std::vector<HalfEdgeChart>& charts, float border) {
    if (charts.empty()) { return; }

    float uMin = std::numeric_limits<float>::max();
    float vMin = std::numeric_limits<float>::max();
    float uMax = std::numeric_limits<float>::lowest();
    float vMax = std::numeric_limits<float>::lowest();

    for (const HalfEdgeChart& ac : charts) {
        for (const auto& uv : ac.uv) {
            uMin = std::min(uMin, uv.x);
            vMin = std::min(vMin, uv.y);
            uMax = std::max(uMax, uv.x);
            vMax = std::max(vMax, uv.y);
        }
    }

    // Reserve the gutter at the atlas edges too: extend the bounds outward by
    // `border` so the outermost charts keep the same margin from UV 0/1 that
    // they keep from each other (otherwise normalisation stretches them to the
    // border and the edge charts touch it).
    uMin -= border; vMin -= border;
    uMax += border; vMax += border;

    float rangeU = uMax - uMin;
    float rangeV = vMax - vMin;
    if (rangeU < 1e-12f || rangeV < 1e-12f) { return; }

    float invU = 1.0f / rangeU;
    float invV = 1.0f / rangeV;

    for (HalfEdgeChart& ac : charts) {
        for (auto& uv : ac.uv) {
            uv.x = (uv.x - uMin) * invU;
            uv.y = (uv.y - vMin) * invV;
        }
        for (auto& g : ac.gutterLoop) {
            g.x = (g.x - uMin) * invU;
            g.y = (g.y - vMin) * invV;
        }
    }
}

// ===========================================================================
// packAtlas — public API
// ===========================================================================

void packAtlas(std::vector<HalfEdgeChart>& charts,
               const MeshData&          data,
               int                      atlasResolution) {
    if (charts.empty()) { return; }

    // 1. Proportional UV scaling to match 3D area ratios.
    scaleProportional(charts, data);

    // 2. Sort by face count descending (matches legacy Cmp_Charts / std::list::sort).
    //    std::stable_sort preserves insertion (chart-id) order for equal counts.
    std::vector<HalfEdgeChart*> order;
    order.reserve(charts.size());
    for (HalfEdgeChart& ac : charts) { order.push_back(&ac); }
    std::stable_sort(order.begin(), order.end(),
        [](const HalfEdgeChart* a, const HalfEdgeChart* b) {
            return a->faces.size() > b->faces.size();
        });

    int res = atlasResolution;
    Horizon        horizon(res);
    ScratchBuffers scratch;
    scratch.lower.assign(static_cast<std::size_t>(res * 2), res - 1);
    scratch.upper.assign(static_cast<std::size_t>(res * 2), 0);

    // 3. Descending-ratio loop: shrink until all charts fit.
    //    Matches legacy: ratio starts at 4.0, decrements by 0.01.
    float ratio   = 4.0f;
    bool  success = false;
    while (!success) {
        success = true;
        ratio  -= 0.01f;
        float scale = ratio / std::sqrt(1.0f);
        horizon.reset();

        for (HalfEdgeChart* ac : order) {
            if (!packSingleChart(*ac, scale, /*write=*/false, horizon, res, scratch)) {
                success = false;
                break;
            }
        }

        if (ratio < 0.01f) { break; }  // safety: degenerate input
    }

    // 4. Final write pass with the winning ratio.
    horizon.reset();
    for (HalfEdgeChart* ac : order) {
        packSingleChart(*ac, ratio / std::sqrt(1.0f), /*write=*/true, horizon, res, scratch);
    }

    // 5. Normalise to [0,1]², reserving the gutter as a border (kGutterPx in the
    //    res-column packing grid → kGutterPx/(res-1) in UV).
    normalizeAtlas(charts, static_cast<float>(kGutterPx) / static_cast<float>(res - 1));
}

// ===========================================================================
// parametrize — public API
// ===========================================================================

namespace {

/// Flatten every chart then scale+pack into [0,1]² — shared by both public forms.
std::vector<HalfEdgeChart> flattenAndPack(const HalfEdgeMesh& mesh, const MeshData& data,
                                       const ChartSet& charts, int atlasResolution) {
    std::vector<HalfEdgeChart> atlasCharts;
    atlasCharts.reserve(charts.charts.size());
    for (const Chart& ch : charts.charts) {
        atlasCharts.push_back(flatten(mesh, data, ch));
    }
    packAtlas(atlasCharts, data, atlasResolution);
    return atlasCharts;
}

}  // namespace

std::vector<decimate::math::Vec2<float>>
parametrize(const HalfEdgeMesh& mesh,
            const MeshData&     data,
            const ChartSet&     charts,
            int                 atlasResolution) {
    const std::vector<HalfEdgeChart> atlasCharts =
        flattenAndPack(mesh, data, charts, atlasResolution);

    // Scatter per-local UV back to per-mesh-vertex UV.
    // LOSSY: a seam vertex (in multiple charts) keeps the last chart's UV.
    // parametrizeCorners() is the seam-correct form; this stays for the legacy contract.
    const int numVerts = static_cast<int>(data.positions.size());
    std::vector<decimate::math::Vec2<float>> perVertexUV(numVerts, decimate::math::Vec2<float>(0.0f, 0.0f));
    for (const HalfEdgeChart& ac : atlasCharts) {
        const int numLocal = static_cast<int>(ac.localToMesh.size());
        for (int li = 0; li < numLocal; ++li) {
            const int gv = ac.localToMesh[li];
            if (gv >= 0 && gv < numVerts) {
                perVertexUV[gv] = ac.uv[static_cast<std::size_t>(li)];
            }
        }
    }
    return perVertexUV;
}

std::vector<decimate::math::Vec2<float>>
parametrizeCorners(const HalfEdgeMesh& mesh,
                   const MeshData&     data,
                   const ChartSet&     charts,
                   int                 atlasResolution) {
    const std::vector<HalfEdgeChart> atlasCharts =
        flattenAndPack(mesh, data, charts, atlasResolution);

    const int numFaces = static_cast<int>(data.indices.size()) / 3;
    std::vector<decimate::math::Vec2<float>> cornerUV(
        static_cast<std::size_t>(numFaces) * 3, decimate::math::Vec2<float>(0.0f, 0.0f));

    // Per chart: mesh-vertex -> local-index (inverse of HalfEdgeChart::localToMesh).
    std::vector<std::unordered_map<int, int>> meshToLocal(atlasCharts.size());
    for (std::size_t ci = 0; ci < atlasCharts.size(); ++ci) {
        const HalfEdgeChart& ac = atlasCharts[ci];
        meshToLocal[ci].reserve(ac.localToMesh.size());
        for (int li = 0; li < static_cast<int>(ac.localToMesh.size()); ++li) {
            meshToLocal[ci][ac.localToMesh[static_cast<std::size_t>(li)]] = li;
        }
    }

    // Each face's three corners take UVs from its OWNING chart only — so no
    // triangle ever straddles a chart seam.
    for (int f = 0; f < numFaces; ++f) {
        const int ci = (f < static_cast<int>(charts.faceToChart.size()))
                           ? charts.faceToChart[static_cast<std::size_t>(f)]
                           : -1;
        if (ci < 0 || ci >= static_cast<int>(atlasCharts.size())) {
            continue;  // unassigned face (should not happen post-segment)
        }
        const HalfEdgeChart& ac  = atlasCharts[static_cast<std::size_t>(ci)];
        const auto&       map = meshToLocal[static_cast<std::size_t>(ci)];
        for (int k = 0; k < 3; ++k) {
            const int v  = data.indices[static_cast<std::size_t>(3 * f + k)];
            const auto it = map.find(v);
            if (it != map.end()) {
                cornerUV[static_cast<std::size_t>(3 * f + k)] =
                    ac.uv[static_cast<std::size_t>(it->second)];
            }
        }
    }
    return cornerUV;
}

}  // namespace decimate::chart
