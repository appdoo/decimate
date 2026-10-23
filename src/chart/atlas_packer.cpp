// This file is part of the decimate library; see decimate.h for version/license details.
/// @file atlas_packer.cpp
/// @brief Atlas packing + full parametrize() orchestrator — Spec 08-02, Beat 2a.
///
/// ## Notes on the 2003 packer (`legacy2003/chartification/catlaspacker.cpp`)
///
/// This file started as a port of the 2003 packer.  It is NOT behavior-preserving any more: the
/// silhouette (`coverSegment`, Spec 19 Amendment 4) and the miter gutter differ from 2003.
///
/// ### ScaleProportional (`legacy2003/chartification/catlaspacker.cpp:398-427`)
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
/// 2003 `Cmp_Charts` sorts descending by `FaceCount()` (`legacy2003/chartification/catlaspacker.cpp:26-34`) with
/// `list.sort` (`:165`; `std::list::sort` is stable, so equal counts keep insertion order).
/// We use `std::stable_sort` with descending face-count and ascending-id
/// tie-break.  Note: the spec says "area desc" but the 2003 criterion is face
/// count (`legacy2003/chartification/catlaspacker.cpp:32`), which is kept.
///
/// ### Pack loop
///
/// 2003 tries a descending `ratio` (`legacy2003/chartification/catlaspacker.cpp:153-179`) (start 4.0, step -0.01) until all charts
/// fit into the grid.  The `ratio` controls the UV-to-pixel scale factor
/// (`scale = ratio / sqrt(1)`).  We preserve this outer loop exactly.
/// `atlasResolution` (default 256) replaces the 2003 constructor argument `discreetSize` (`legacy2003/chartification/catlaspacker.cpp:39-42`).
///
/// ### 16-rotation search
///
/// For each rotation step (0°, 22.5°, …, 337.5°) the chart boundary UVs are
/// discretised into the `upper` / `lower` column arrays, then `FindOptimal`
/// scans for the placement with least wasted horizon area.  The placement
/// with the globally smallest waste across all 16 rotations wins.
///
/// The silhouette footprint is built in `buildSilhouetteFootprint` from the chart's boundary loop
/// (`boundaryPointIndices()` of the local `HalfEdgeMesh`), offset outward by the miter gutter and
/// rasterised per column.  The chart is a valid 2D embedding (no foldovers), so the boundary bounds
/// the interior; the zero-overlap SAT test gates the result.
///
/// ### parametrize() UV scatter
///
/// Per-mesh-vertex UVs are written chart by chart in order (chart 0, 1, …).
/// A seam vertex shared by multiple charts is overwritten each time, ending
/// with the LAST chart's UV (the port's choice for this lossy per-vertex form).
///
/// ### Boundary contract
/// No GL, no OS, no I/O, no global mutable state.

#include <decimate/chart/parametrize.h>

#include "chart_quality.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <span>

#include <decimate/math/vec2.h>
#include <decimate/math/vec3.h>
#include <decimate/chart/chart.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>

namespace decimate::chart {

namespace {

/// @brief Tracks the current packed height at each column of the atlas grid.
///
/// `data_[i]` is the height at column `i` (0 = empty column).  After each
/// chart is placed, its upper silhouette is merged into the horizon via
/// `update()`.
class Horizon {
public:
    explicit Horizon(int resolution)
        : res_(resolution), data_(resolution, 0) {
    }

    /// @brief Reset all columns to zero (empty atlas).
    void reset() {
        std::fill(data_.begin(), data_.end(), 0);
    }

    /// @brief Height at column `index`.
    int operator[](int index) const noexcept {
        return data_[index];
    }

    /// @brief Total discrete area of the atlas (res²).
    [[nodiscard]] int totalArea() const noexcept {
        return res_ * res_;
    }

    /// @brief Wasted area if a chart with lower silhouette `lower` is placed at (x, y).
    ///
    /// `lower[i]` is the per-column minimum v of the chart (its bottom silhouette).
    /// Returns the sum of pixel gaps between the horizon top and the chart bottom,
    /// or -1 if any column of the chart falls outside the atlas boundary or
    /// collides with already-packed content.
    ///
    /// Port of `CHorizon::Area(int x, int y, int size, int *chart)` where `chart`
    /// is the `lower` (bottom silhouette) array.
    [[nodiscard]] int wastedArea(int x, int y, int size, // NOLINT(bugprone-easily-swappable-parameters): (x, y, size) are distinct grid coordinates
                                 const std::vector<int> &lower) const {
        int wasted = 0;
        for (int i = 0; i < size; ++i) {
            const int col = x + i;
            if (col >= res_) {
                return -1; // the chart pokes past the atlas edge
            }

            const int bottom = y + lower[i]; // chart bottom at this column
            if (data_[col] >= bottom) {
                return -1; // the horizon already reaches the chart bottom: collision
            }

            wasted += bottom - data_[col]; // gap between horizon and chart bottom
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
    void update(int xoffs, int yoffs, int size, // NOLINT(bugprone-easily-swappable-parameters): column-offset, row-offset, span-width
                const std::vector<int> &upper) {
        for (int i = 0; i < size; ++i) {
            const int col = xoffs + i;
            if (col < res_) {
                data_[col] = upper[i] + yoffs; // gutter already in yoffs (findOptimal)
            }
        }
    }

    /// @brief Atlas grid side length in columns.
    [[nodiscard]] int resolution() const noexcept {
        return res_;
    }

private:
    int              res_;
    std::vector<int> data_;
};


/// @brief Result of a single placement scan.
struct OptimalPlacement {
    int x;     ///< Best column (−1 = not found).
    int y;     ///< Best row (−1 = not found).
    int waste; ///< Wasted pixel area at (x, y); horizon.totalArea() if not found.
};


/// @brief Chart footprint in pixel space for one rotation.
struct RotationFootprint {
    int width = 0;  ///< AABB width in grid columns (0 = degenerate / skip).
    int height = 0; ///< AABB height in grid rows    (0 = degenerate / skip).
};


/// @brief Winning placement result for the 16-rotation search.
struct BestRotation {
    float angleDeg = 0.0f; ///< Rotation angle that minimises wasted area.
    int   x = -1;          ///< Grid column of the best placement (−1 = none found).
    int   y = -1;          ///< Grid row   of the best placement.
    int   width = 0;       ///< Chart width in grid columns at this rotation.
};


/// @brief Scratch column buffers shared across all rotation steps.
///
/// `lower` and `upper` are reused per rotation (reset in `buildSilhouetteFootprint`).
/// `best` captures the winning rotation's upper silhouette for the horizon update.
struct ScratchBuffers {
    std::vector<int> lower; ///< Per-column minimum v (reset each rotation).
    std::vector<int> upper; ///< Per-column maximum v (reset each rotation).
    std::vector<int> best;  ///< Upper silhouette of the best rotation found.
};

} // namespace


/// @brief Half-cross-product magnitude — 3D area of a triangle.
static float triArea3D(const math::Vec3f &a,
                       const math::Vec3f &b,
                       const math::Vec3f &c) {
    return (b - a).cross(c - a).length() * 0.5f;
}


/// @brief 3D surface area of an HalfEdgeChart summed over its face list.
static float chartArea3D(const HalfEdgeChart &ac, const MeshData &data) {
    float       total = 0.0f;
    const auto &p = data.positions;
    for (int fi : ac.faces) {
        const math::Vec3f &a = p[static_cast<std::size_t>(data.indices[3 * fi])];
        const math::Vec3f &b = p[static_cast<std::size_t>(data.indices[3 * fi + 1])];
        const math::Vec3f &c = p[static_cast<std::size_t>(data.indices[3 * fi + 2])];
        total += triArea3D(a, b, c);
    }
    return total;
}


/// @brief UV (shoelace) area of a chart, summed over its local triangles.
///
/// Port of `CAtlasPacker::TextureArea(CChart *)`.  Corners come from the chart's local
/// triangulation (`localCorner`), so a self-seam cut chart is measured per corner.
///
/// @param ac  Atlas chart (uv array indexed by local vertex id).
static float chartAreaUV(const HalfEdgeChart &ac) {
    float     area = 0.0f;
    const int numFaces = static_cast<int>(ac.faces.size());
    for (int i = 0; i < numFaces; ++i) {
        const auto &a = ac.uv[static_cast<std::size_t>(localCorner(ac, i, 0))];
        const auto &b = ac.uv[static_cast<std::size_t>(localCorner(ac, i, 1))];
        const auto &c = ac.uv[static_cast<std::size_t>(localCorner(ac, i, 2))];
        // Shoelace determinant, |det|/2 per face.  NOT the 2003 expression (`TextureArea`,
        // `legacy2003/chartification/catlaspacker.cpp:464` applies fabs to the first term only).
        const float det = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
        area += std::abs(det) * 0.5f;
    }
    return area;
}


/// @brief Scale each chart's UVs so that UV area ratios match 3D area ratios.
///
/// Port of `CAtlasPacker::ScaleProportional(CChartSet *)`.
/// After scaling, each chart's UV area equals `chartArea3D / totalArea3D`.
/// A chart with zero 3D or UV area is skipped (the 2003 `IsFZero` guard, `legacy2003/chartification/catlaspacker.cpp:411`).
static void scaleProportional(std::vector<HalfEdgeChart> &charts, const MeshData &data) {
    // Each chart's share is its area over the total, so the total comes first.
    float totalArea3D = 0.0f;
    for (const HalfEdgeChart &ac : charts) {
        totalArea3D += chartArea3D(ac, data);
    }
    if (!(totalArea3D >= 1e-12f)) {
        return; // degenerate input (zero, or NaN from a non-finite position): the share would divide by zero
    }

    for (HalfEdgeChart &ac : charts) {
        const float area3D = chartArea3D(ac, data);
        const float areaUV = chartAreaUV(ac);

        if (!(areaUV >= 1e-12f) || !(area3D >= 1e-12f)) {
            continue; // a zero-area or non-finite chart has no scale to match (NaN fails `>=`)
        }

        // 2003 formula (`legacy2003/chartification/catlaspacker.cpp:417`): cscale = sqrt(chartArea / (totalArea * textureArea)).
        const float cscale = std::sqrt(area3D / (totalArea3D * areaUV));

        for (math::Vec2f &uv : ac.uv) {
            uv.x *= cscale;
            uv.y *= cscale;
        }
    }
}


/// Component-wise minimum corner of the AABB of `pts` (non-empty).
static Tex2 minTex(const std::vector<Tex2> &pts) {
    Tex2 m = pts[0];
    for (const Tex2 &p : pts) {
        m.u = std::min(m.u, p.u);
        m.v = std::min(m.v, p.v);
    }
    return m;
}


/// Component-wise maximum corner of the AABB of `pts` (non-empty).
static Tex2 maxTex(const std::vector<Tex2> &pts) {
    Tex2 m = pts[0];
    for (const Tex2 &p : pts) {
        m.u = std::max(m.u, p.u);
        m.v = std::max(m.v, p.v);
    }
    return m;
}


/// @brief Centre of the AABB spanned by `pts`.
///
/// Port of `CAtlasPacker::Mid(int size, Vertex *vertices)`.
/// 2003 computes `(max - min) / 2` (`legacy2003/chartification/catlaspacker.cpp:118-119`), which is the half-range, not the centre.
/// We preserve that exact computation for rotation-centre consistency.
static Tex2 midTex(const std::vector<Tex2> &pts) {
    const Tex2 mn = minTex(pts);
    const Tex2 mx = maxTex(pts);
    return { (mx.u - mn.u) * 0.5f, (mx.v - mn.v) * 0.5f };
}


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
static OptimalPlacement findOptimal(int                     width,
                                    const std::vector<int> &lower,
                                    const std::vector<int> &upper,
                                    const Horizon          &horizon,
                                    int                     gutter) {
    const int        res = horizon.resolution();
    OptimalPlacement best { -1, -1, horizon.totalArea() };

    for (int x = 0; x < res - width; ++x) {
        // One fused pass over the chart width: area sum, deepest column (-> landing y) and crown.
        int rawArea = 0;
        int minGap = std::numeric_limits<int>::max();
        int maxUp = 0;
        for (int k = 0; k < width; ++k) {
            const int colGap = lower[static_cast<std::size_t>(k)] - horizon[x + k];
            rawArea += colGap;
            if (colGap < minGap) {
                minGap = colGap;
            }
            if (upper[static_cast<std::size_t>(k)] > maxUp) {
                maxUp = upper[static_cast<std::size_t>(k)];
            }
        }

        // Drop the deepest column onto the skyline, plus the gutter.
        const int y = -minGap + gutter;
        if (y < 0 || y + maxUp + gutter >= res) {
            continue; // the chart would not fit under the atlas ceiling
        }

        // The area at the landing falls out of the same pass.
        const int area = rawArea - width * minGap;
        if (area < best.waste) {
            best = { x, y, area }; // deepest-notch (min-shadow) column so far
        }
    }
    return best;
}


/// @brief Cover one boundary edge in the per-column lower/upper silhouette, conservatively.
///
/// Replaces the 2003 `CAtlasPacker::RasterLine` sampler (`legacy2003/chartification/catlaspacker.cpp:339-367`;
/// Spec 19 Amendment 4).  That sampler took ONE integer-interpolated row per column
/// (`v = v1 + (dv*(u-u1))/du`, `:362-367`), so a steep edge (many rows inside one column) and the
/// column's end points were under-covered: the footprint sat inside the exact float miter
/// polygon and neighbouring charts' gutters, even their interiors, could poke into it.  Here
/// every column the segment touches gets `lower = min floor(v)` and `upper = max floor(v)` over
/// the part of the segment inside that column (closed interval, so an end point on a column
/// boundary covers both columns).  The silhouette therefore contains every pixel the polygon's
/// boundary touches.  End points must already lie in [0, res-1].
///
/// ## Why conservative, and what the deleted Bresenham experiment still teaches
///
/// Still valid from the benchmark behind the deleted `rasterLine` comment (`tests/raster_bench.cpp`,
/// cow -> 1000 faces, 43 charts, git 648e6d6^):
///  - SPEED: a per-pixel Bresenham DDA is NOT faster than a per-column loop.  It steps once per
///    pixel of edge length, a per-column loop once per column, and the per-column work is cheap
///    and predictable.  `coverSegment` is per-column, so it keeps that property.
///  - HOT PATH: the silhouette is not the packer's cost; `findOptimal`'s O(res^2) position scan is.
///  - CORRECTNESS, not speed, is what made the atlas look better: the 2003 sampler DROPPED a
///    whole edge when an end point was out of bounds (`legacy2003/chartification/catlaspacker.cpp:347-349`, a hole in the
///    silhouette).  `buildSilhouetteFootprint` clamps the end points into the grid instead, so no
///    edge is lost.
/// Conservative here means: the covered cell set is a SUPERSET of the pixels the exact float
/// segment passes through (floor of the extreme rows per column, closed column interval), so the
/// footprint never lies inside the exact miter polygon.  It over-covers by at most one row at
/// either end of a column; that costs a little density (measured, Spec 19 report) and buys
/// gutters that cannot interpenetrate.
// NOLINTNEXTLINE(readability-function-size,bugprone-easily-swappable-parameters): segment end points + silhouette arrays; a struct would only wrap them
void coverSegment(float x1, float y1, float x2, float y2, int res,
                  std::vector<int> &lower, std::vector<int> &upper) {
    if (x1 > x2) {
        std::swap(x1, x2);
        std::swap(y1, y2);
    }
    const int   cLo = static_cast<int>(std::floor(x1));
    const int   cHi = std::min(res - 1, static_cast<int>(std::floor(x2)));
    const float slope = (x2 > x1) ? (y2 - y1) / (x2 - x1) : 0.0f;
    for (int u = cLo; u <= cHi; ++u) {
        const float xa = std::max(x1, static_cast<float>(u));
        const float xb = std::min(x2, static_cast<float>(u + 1));
        const float ya = y1 + slope * (xa - x1);
        const float yb = (x2 > x1) ? y1 + slope * (xb - x1) : y2;
        const int   lo = static_cast<int>(std::floor(std::min(ya, yb)));
        const int   hi = static_cast<int>(std::floor(std::max(ya, yb)));
        upper[static_cast<std::size_t>(u)] = std::max(upper[static_cast<std::size_t>(u)], std::min(hi, res - 1));
        lower[static_cast<std::size_t>(u)] = std::min(lower[static_cast<std::size_t>(u)], std::max(lo, 0));
    }
}


/// Twice the signed area of a closed loop (shoelace); > 0 is counter-clockwise.
static float signedArea2(const std::vector<Tex2> &p) {
    float a2 = 0.0f;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const Tex2 &a = p[i];
        const Tex2 &b = p[(i + 1) % p.size()];
        a2 += a.u * b.v - b.u * a.v;
    }
    return a2;
}


/// Crossing point of edges i and j of @p p, or false when they do not cross properly.
static bool edgesCrossAt(const std::vector<Tex2> &p, std::size_t i, std::size_t j, float tol, Tex2 &x) {
    const std::size_t n = p.size();
    const Tex2       &a = p[i];
    const Tex2       &b = p[(i + 1) % n];
    const Tex2       &c = p[j];
    const Tex2       &d = p[(j + 1) % n];
    if (!math::segmentsCross(math::Vec2f { a.u, a.v }, math::Vec2f { b.u, b.v }, math::Vec2f { c.u, c.v },
                             math::Vec2f { d.u, d.v })) {
        return false;
    }

    // The crossing parameter along a-b comes from the signed distances of a and b to line c-d.
    const double da = orientValue(c.u, c.v, d.u, d.v, a.u, a.v);
    const double db = orientValue(c.u, c.v, d.u, d.v, b.u, b.v);

    // A touch, not a crossing: an end point within `tol` of the other segment's line.
    if (tol > 0.0f) {
        const double lcd = std::hypot(static_cast<double>(d.u) - c.u, static_cast<double>(d.v) - c.v);
        const double lab = std::hypot(static_cast<double>(b.u) - a.u, static_cast<double>(b.v) - a.v);
        const double dc = orientValue(a.u, a.v, b.u, b.v, c.u, c.v);
        const double dd = orientValue(a.u, a.v, b.u, b.v, d.u, d.v);
        if (std::min(std::abs(da) / lcd, std::abs(db) / lcd) <= tol || std::min(std::abs(dc) / lab, std::abs(dd) / lab) <= tol) {
            return false;
        }
    }

    const auto t = static_cast<float>(da / (da - db));
    x = { a.u + t * (b.u - a.u), a.v + t * (b.v - a.v) };
    return true;
}


/// Can edges e and f (indices into @p p) cross at all: non-adjacent, y-intervals overlapping?
static bool mayCross(const std::vector<Tex2> &p, std::size_t e, std::size_t f) {
    const std::size_t n = p.size();
    const std::size_t lo = std::min(e, f);
    const std::size_t hi = std::max(e, f);
    if (hi - lo < 2 || (lo == 0 && hi == n - 1)) {
        return false; // same or neighbouring edges share an end point, they cannot cross properly
    }

    const auto yLo = [&](std::size_t k) {
        return std::min(p[k].v, p[(k + 1) % n].v);
    };
    const auto yHi = [&](std::size_t k) {
        return std::max(p[k].v, p[(k + 1) % n].v);
    };
    return yLo(f) <= yHi(e) && yHi(f) >= yLo(e);
}


/// A proper crossing of the closed polyline @p p: edges i < j (non-adjacent), crossing point @p x.
/// Sweep over the edges sorted by their left end (x-interval overlap, then y-interval overlap), so a
/// thousand-vertex boundary costs O(B log B + candidate pairs) instead of O(B^2).
// NOLINTNEXTLINE(readability-function-size): sweep kernel, sorted edge x outer loop x active-window inner loop
static bool firstCrossing(const std::vector<Tex2> &p, float tol, std::size_t &i, std::size_t &j, Tex2 &x) {
    const std::size_t        n = p.size();
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t { 0 });
    const auto xMin = [&](std::size_t e) {
        return std::min(p[e].u, p[(e + 1) % n].u);
    };
    const auto xMax = [&](std::size_t e) {
        return std::max(p[e].u, p[(e + 1) % n].u);
    };

    // Sorted by left end, the inner loop can stop at the first edge that starts right of this one.
    std::ranges::sort(order, [&](std::size_t l, std::size_t r) { return xMin(l) < xMin(r); });
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t e = order[k];
        for (std::size_t m = k + 1; m < n && xMin(order[m]) <= xMax(e); ++m) {
            const std::size_t f = order[m];
            if (mayCross(p, e, f) && edgesCrossAt(p, std::min(e, f), std::max(e, f), tol, x)) {
                i = std::min(e, f);
                j = std::max(e, f);
                return true;
            }
        }
    }
    return false;
}


/// Cut every fold lobe out of a closed offset polyline.  At a crossing of edges i and j the loop splits into
/// the inner part {X, p[i+1..j]} and the outer part {p[0..i], X, p[j+1..]}; the part with the smaller signed
/// area (the inverted lobe a reflex vertex folds back over the chart, whose even-odd region is winding 2)
/// is dropped.  Each cut removes at least one vertex, so the loop terminates with a simple polygon.
static void removeFoldLoops(std::vector<Tex2> &p) {
    std::size_t i = 0;
    std::size_t j = 0;
    Tex2        x;
    while (p.size() >= 4 && firstCrossing(p, 0.0f, i, j, x)) {
        // Split the loop at the crossing X into the inner part {X, p[i+1..j]} and the outer part.
        std::vector<Tex2> inner(p.begin() + static_cast<std::ptrdiff_t>(i) + 1,
                                p.begin() + static_cast<std::ptrdiff_t>(j) + 1);
        inner.push_back(x);
        std::vector<Tex2> outer(p.begin(), p.begin() + static_cast<std::ptrdiff_t>(i) + 1);
        outer.push_back(x);
        outer.insert(outer.end(), p.begin() + static_cast<std::ptrdiff_t>(j) + 1, p.end());

        // Keep the part that winds like the whole loop; the dropped one is the inverted lobe.
        const float wind = signedArea2(p) >= 0.0f ? 1.0f : -1.0f;
        p = (wind * signedArea2(inner) > wind * signedArea2(outer)) ? std::move(inner) : std::move(outer);
    }
}


bool hasSelfCrossing(const std::vector<Tex2> &loop, float touchTol) {
    std::size_t i = 0;
    std::size_t j = 0;
    Tex2        x;
    return loop.size() >= 4 && firstCrossing(loop, touchTol, i, j, x);
}


/// The raw miter offset: one moved vertex per input vertex, fold lobes still in (`miterOffsetLoop` cuts them).
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters): (gutter, miterLimit) mirror miterOffsetLoop
static std::vector<Tex2> miterOffsetRaw(const std::vector<Tex2> &loop, float gutter, float miterLimit) {
    const std::size_t n = loop.size();
    if (n < 3) {
        return loop; // no polygon to offset
    }

    // The outward side depends on the winding: interior is left of each edge when the area is > 0 (CCW).
    float area2 = 0.0f;
    for (std::size_t i = 0; i < n; ++i) {
        const Tex2 &a = loop[i];
        const Tex2 &b = loop[(i + 1) % n];
        area2 += a.u * b.v - b.u * a.v;
    }
    const float wind = (area2 >= 0.0f) ? 1.0f : -1.0f;

    // Outward unit normal of a directed edge (right-hand perp, winding-signed).
    auto outNormal = [wind](const Tex2 &a, const Tex2 &b, float &nx, float &ny) {
        const float dx = b.u - a.u;
        const float dy = b.v - a.v;
        const float len = std::sqrt(dx * dx + dy * dy);
        if (len < 1e-6f) { // a zero-length edge has no direction; it contributes no normal
            nx = 0.0f;
            ny = 0.0f;
            return;
        }
        nx = wind * dy / len;
        ny = -wind * dx / len;
    };

    // Move along the bisector of the two edge normals, length gutter/cos(half-angle),
    // clamped at miterLimit x.  Outward-only => the silhouette only grows.
    const float       minCos = 1.0f / std::max(miterLimit, 1.0f);
    std::vector<Tex2> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        const Tex2 &c = loop[i];
        float       n1x = 0.0f;
        float       n1y = 0.0f;
        float       n2x = 0.0f;
        float       n2y = 0.0f;
        outNormal(loop[(i + n - 1) % n], c, n1x, n1y); // incoming edge
        outNormal(c, loop[(i + 1) % n], n2x, n2y);     // outgoing edge

        // Straight continuation or a 180-degree spike has no bisector: fall back to the incoming normal.
        float       mx = n1x + n2x;
        float       my = n1y + n2y;
        const float ml = std::sqrt(mx * mx + my * my);
        if (ml < 1e-6f) {
            mx = n1x;
            my = n1y;
        } else {
            mx /= ml;
            my /= ml;
        }

        const float cosHalf = std::max(mx * n1x + my * n1y, minCos);
        const float d = gutter / cosHalf;
        out[i] = { c.u + mx * d, c.v + my * d };
    }
    return out;
}


std::vector<Tex2> miterOffsetLoop(const std::vector<Tex2> &loop, float gutter, float miterLimit) {
    std::vector<Tex2> out = miterOffsetRaw(loop, gutter, miterLimit);
    removeFoldLoops(out);
    return out;
}


/// @brief Chart boundary loop in padded pixel space for one rotation.
///
/// Pixel = (uv - mn) * pxScale + kGutterPx: the +kGutterPx pad keeps the plain outward
/// offset in-grid on all four sides (a clamped miter spike can still leave it by up to
/// (limit-1) * kGutterPx; `spikeShift` translates the footprint for that); applyBestPlacement
/// adds the matching +kGutterPx so the chart sits inside its padded footprint.
/// Empty when the chart has fewer than 3 boundary vertices.
static std::vector<Tex2> boundaryLoopPx(const std::vector<Tex2> &dst,
                                        const HalfEdgeMesh      &localMesh,
                                        const Tex2              &mn,
                                        float                    pxScale) {
    const std::vector<std::uint32_t> loop = localMesh.boundaryPointIndices();
    if (loop.size() < 3) {
        return {}; // no outline to offset
    }

    const auto        pad = static_cast<float>(kGutterPx);
    std::vector<Tex2> px(loop.size());
    for (std::size_t i = 0; i < loop.size(); ++i) {
        const Tex2 &p = dst[static_cast<std::size_t>(loop[i])];
        px[i] = { (p.u - mn.u) * pxScale + pad, (p.v - mn.v) * pxScale + pad };
    }
    return px;
}


/// Translation that brings a raw miter offset back into the footprint grid.  A clamped-miter spike leaves the
/// chart by up to `kDefaultMiterLimit * kGutterPx`, but the chart pad is only `kGutterPx`, so such a spike
/// ends at a negative pixel coordinate.  The footprint used to clamp it to 0, cutting the spike flat while
/// the published `gutterLoop` kept it; the cut spike then reached into the neighbouring chart (the last
/// horse/128 intInGut texel, measured).  Shifting instead keeps the whole spike in the footprint.
static Tex2 spikeShift(const std::vector<Tex2> &off) {
    Tex2 s;
    for (const Tex2 &p : off) {
        s.u = std::max(s.u, -p.u);
        s.v = std::max(s.v, -p.v);
    }
    return s;
}


/// Rasterise the offset loop @p off (shifted by `spikeShift`, clamped into the grid) into the per-column
/// silhouette.  @returns the footprint width in columns.
static int coverOffsetLoop(const std::vector<Tex2> &off, int res, std::vector<int> &lower, std::vector<int> &upper) {
    const std::size_t n = off.size();
    const Tex2        shift = spikeShift(off);
    const float       hi = static_cast<float>(res - 1);

    // Clamp into the grid so coverSegment never indexes outside the column arrays.
    std::vector<Tex2> pt(n);
    int               width = 0;
    for (std::size_t i = 0; i < n; ++i) {
        pt[i] = { std::clamp(off[i].u + shift.u, 0.0f, hi), std::clamp(off[i].v + shift.v, 0.0f, hi) };
        width = std::max(width, static_cast<int>(pt[i].u) + 1);
    }

    for (std::size_t i = 0; i < n; ++i) {
        const Tex2 &a = pt[i];
        const Tex2 &b = pt[(i + 1) % n];
        coverSegment(a.u, a.v, b.u, b.v, res, lower, upper);
    }
    return width;
}


/// @brief Build the chart's per-column silhouette footprint for one rotation.
///
/// Rasterises the MITER-OFFSET boundary loop into `lower`/`upper` via `coverSegment` (conservative, Spec 19-03).
/// The loop is the half-edges of the cut mesh with no twin — the faithful 2003
/// silhouette (`legacy2003/chartification/cchartification.cpp:395`,
/// `legacy2003/split/splitter.cpp:1600`).  The chart is a valid 2D embedding, so
/// the boundary bounds the interior.  The offset itself is `miterOffsetLoop`.
static RotationFootprint buildSilhouetteFootprint(const std::vector<Tex2> &dst,
                                                  const HalfEdgeMesh      &localMesh,
                                                  const Tex2              &mn,
                                                  float                    scale,
                                                  int                      res,
                                                  std::vector<int>        &lower,
                                                  std::vector<int>        &upper) {
    // Reset to the empty silhouette (lower above everything, upper below everything) so coverSegment can min/max into it.
    for (int col = 0; col < res; ++col) {
        upper[static_cast<std::size_t>(col)] = 0;
        lower[static_cast<std::size_t>(col)] = res - 1;
    }

    const float             pxScale = scale * static_cast<float>(res - 1);
    const std::vector<Tex2> loop = boundaryLoopPx(dst, localMesh, mn, pxScale);
    if (loop.empty()) {
        return { 0, 0 }; // no boundary loop: the chart has no footprint
    }

    // The RAW offset, fold lobes included: every vertex of the lobe-free `gutterLoop` lies on a raw segment, so
    // the raw per-column envelope contains it (the footprint bounds the published gutter), and the cut
    // (a global O(cuts * B log B) pass, thousands of times over 16 rotations) is paid once per chart in
    // applyBestPlacement instead.  Measured: cutting per rotation made horse/16 pack 0.4 s -> 18 s.
    const std::vector<Tex2> off = miterOffsetRaw(loop, static_cast<float>(kGutterPx), kDefaultMiterLimit);

    const int width = coverOffsetLoop(off, res, lower, upper);
    if (width <= 0 || width > res) {
        return { 0, 0 }; // degenerate or wider than the atlas
    }

    // The height decides whether the chart can fit at all.
    int ymin = res;
    int ymax = 0;
    for (int col = 0; col < width; ++col) {
        ymin = std::min(ymin, lower[static_cast<std::size_t>(col)]);
        ymax = std::max(ymax, upper[static_cast<std::size_t>(col)]);
    }
    const int height = ymax - ymin;
    if (height <= 0 || height >= res) {
        return { 0, 0 }; // flat, or taller than the atlas
    }

    return { width, height };
}


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
static void applyBestPlacement(HalfEdgeChart           &ac,
                               const std::vector<Tex2> &src,
                               const Tex2              &mid,
                               const BestRotation      &best,
                               float                    scale,
                               int                      res) {
    const int   numUV = static_cast<int>(src.size());
    const float rad = best.angleDeg * (std::numbers::pi_v<float> / 180.0f);
    const float sina = std::sin(rad);
    const float cosa = std::cos(rad);

    // Replay the winning rotation (the search kept only the angle).
    std::vector<Tex2> dst(numUV);
    for (int i = 0; i < numUV; ++i) {
        const float u = src[i].u - mid.u;
        const float v = src[i].v - mid.v;
        dst[i].u = mid.u + u * cosa - v * sina;
        dst[i].v = mid.v + u * sina + v * cosa;
    }


    const Tex2  mn = minTex(dst);
    const float invRes = 1.0f / static_cast<float>(res - 1);
    // Half-pixel offset, as 2003 `(best_x + 0.5f) / (discreetSize - 1)` (`legacy2003/chartification/catlaspacker.cpp:325`).
    // +kGutterPx: the footprint was padded by kGutterPx (buildSilhouetteFootprint),
    // so the chart's own min sits kGutterPx inside its placed footprint — the gutter.
    const float             pxScale = scale * static_cast<float>(res - 1);
    const std::vector<Tex2> loop = boundaryLoopPx(dst, ac.mesh, mn, pxScale);
    const Tex2              shift = spikeShift(miterOffsetRaw(loop, static_cast<float>(kGutterPx), kDefaultMiterLimit));
    const float             offX = (static_cast<float>(best.x + kGutterPx) + shift.u + 0.5f) * invRes;
    const float             offY = (static_cast<float>(best.y + kGutterPx) + shift.v + 0.5f) * invRes;


    for (int i = 0; i < numUV; ++i) {
        ac.uv[i].x = (dst[i].u - mn.u) * scale + offX;
        ac.uv[i].y = (dst[i].v - mn.v) * scale + offY;
    }

    // Publish the placed gutter outline.  The miter-offset loop is in padded pixel space
    // q = (dst-mn)*pxScale + kGutterPx, so uv = (q - kGutterPx)*invRes + offX
    //                                         = (q + best.x + 0.5)*invRes.
    const std::vector<Tex2> off = miterOffsetLoop(loop, static_cast<float>(kGutterPx));
    const auto              bx = static_cast<float>(best.x) + shift.u + 0.5f;
    const auto              by = static_cast<float>(best.y) + shift.v + 0.5f;
    ac.gutterLoop.resize(off.size());
    for (std::size_t i = 0; i < off.size(); ++i) {
        ac.gutterLoop[i] = { (off[i].u + bx) * invRes, (off[i].v + by) * invRes };
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
/// @param[out] scratch  Reusable column buffers; `scratch.best` holds result.
/// @returns Best rotation + placement found; BestRotation::x == -1 when none.
static BestRotation findBestRotation(const std::vector<Tex2> &src,
                                     const HalfEdgeMesh      &localMesh,
                                     const Tex2              &mid,
                                     float                    scale,
                                     int                      res,
                                     const Horizon           &horizon,
                                     ScratchBuffers          &scratch) {
    constexpr int   kSteps = 16;
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
    float             bestPerCol = std::numeric_limits<float>::max();
    BestRotation      best;
    std::vector<Tex2> dst(src.size());

    for (int step = 0; step < kSteps; ++step) {
        const float angleDeg = static_cast<float>(step) * kAngleInc;
        const float rad = angleDeg * (std::numbers::pi_v<float> / 180.0f);
        const float sina = std::sin(rad);
        const float cosa = std::cos(rad);

        // Rotate around the AABB centre so every step starts from the same pre-pack UVs.
        for (int i = 0; i < static_cast<int>(src.size()); ++i) {
            const float u = src[i].u - mid.u;
            const float v = src[i].v - mid.v;
            dst[i].u = mid.u + u * cosa - v * sina;
            dst[i].v = mid.v + u * sina + v * cosa;
        }

        const Tex2              mn = minTex(dst);
        const RotationFootprint fp =
            buildSilhouetteFootprint(dst, localMesh, mn, scale, res,
                                     scratch.lower, scratch.upper);
        if (fp.width <= 0) {
            continue; // degenerate footprint at this angle
        }

        // Vertical-only gutter is now 0: the all-around miter offset
        // (buildSilhouetteFootprint, kGutterPx) provides the chart separation on
        // every side, so a further vertical gutter here would double-count it.
        const OptimalPlacement pl =
            findOptimal(fp.width, scratch.lower, scratch.upper, horizon, 0);
        if (pl.x == -1) {
            continue; // does not fit anywhere at this angle
        }

        const float perCol = static_cast<float>(pl.waste) / static_cast<float>(fp.width);
        if (perCol < bestPerCol) {
            bestPerCol = perCol;
            best = { angleDeg, pl.x, pl.y, fp.width };
            scratch.best.assign(scratch.upper.begin(),
                                scratch.upper.begin() + fp.width);
        }
    }
    return best;
}


/// @brief Try 16 rotations of `ac` and pick the best placement in `horizon`.
///
/// Port of `CAtlasPacker::Pack(CChart*, float scale, bool write)`.
/// Each rotation is placed by its per-column silhouette footprint (`buildSilhouetteFootprint`).
///
/// @param ac      Chart — `uv` mutated in-place when `write` is true.
/// @param scale   UV-to-pixel scale (pixels per UV unit).
/// @param write   When true, apply the best rotation + translation to `ac.uv`.
/// @param horizon Atlas horizon (updated to include the placed chart).
/// @param res     Atlas grid resolution.
/// @param scratch Reusable column buffers (shared across all charts in a pass).
/// @returns true if a valid placement was found.
static bool packSingleChart(HalfEdgeChart  &ac,
                            float           scale,
                            bool            write,
                            Horizon        &horizon,
                            int             res,
                            ScratchBuffers &scratch) {
    const int numUV = static_cast<int>(ac.uv.size());
    if (numUV == 0) {
        return true; // nothing to place, nothing to fail
    }

    // The search works on a Tex2 copy so the chart's own UVs stay untouched in the dry pass.
    std::vector<Tex2> src(numUV);
    for (int i = 0; i < numUV; ++i) {
        src[i].u = ac.uv[i].x;
        src[i].v = ac.uv[i].y;
    }
    const Tex2 mid = midTex(src);

    // The local mesh is passed so buildSilhouetteFootprint can walk its boundary loop.
    const BestRotation best =
        findBestRotation(src, ac.mesh, mid, scale, res, horizon, scratch);

    // The next chart must sit above this one's crown, so the horizon takes its upper silhouette.
    if (best.x != -1) {
        horizon.update(best.x, best.y, best.width, scratch.best);
    }

    // The dry pass only asks "does everything fit"; UVs are written in the final pass.
    if (write && best.x != -1) {
        applyBestPlacement(ac, src, mid, best, scale, res);
    }

    return (best.x != -1);
}


/// @brief Translate + uniformly scale all chart UVs so the aggregate AABB is [0,1]².
static void normalizeAtlas(std::vector<HalfEdgeChart> &charts, float border) {
    if (charts.empty()) {
        return;
    }

    float uMin = std::numeric_limits<float>::max();
    float vMin = std::numeric_limits<float>::max();
    float uMax = std::numeric_limits<float>::lowest();
    float vMax = std::numeric_limits<float>::lowest();

    // The gutter outlines lie inside the UV bounds, so the UVs alone give the aggregate box.
    for (const HalfEdgeChart &ac : charts) {
        for (const auto &uv : ac.uv) {
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
    uMin -= border;
    vMin -= border;
    uMax += border;
    vMax += border;

    const float rangeU = uMax - uMin;
    const float rangeV = vMax - vMin;
    if (rangeU < 1e-12f || rangeV < 1e-12f) {
        return; // all UVs coincide: nothing to stretch
    }

    const float invU = 1.0f / rangeU;
    const float invV = 1.0f / rangeV;

    for (HalfEdgeChart &ac : charts) {
        for (auto &uv : ac.uv) {
            uv.x = (uv.x - uMin) * invU;
            uv.y = (uv.y - vMin) * invV;
        }
        for (auto &g : ac.gutterLoop) {
            g.x = (g.x - uMin) * invU;
            g.y = (g.y - vMin) * invV;
        }
    }
}


void packAtlas(std::vector<HalfEdgeChart> &charts,
               const MeshData             &data,
               int                         atlasResolution) {
    // A grid under two columns has no width to pack into (res - 1 divides the gutter, res * 2 sizes the buffers).
    if (atlasResolution < 2) {
        throw std::invalid_argument("packAtlas: atlasResolution must be >= 2");
    }
    if (charts.empty()) {
        return;
    }

    // A closed chart has no boundary loop, so no rotation can be placed and every ratio would fail:
    // the cutter left it without a seam, which is a segmenter defect, not a packing case.
    for (std::size_t i = 0; i < charts.size(); ++i) {
        const HalfEdgeChart &ac = charts[i];
        if (ac.mesh.faceCount() > 0 && ac.mesh.firstBoundaryEdge() == kNoHalfedge) {
            throw std::invalid_argument("packAtlas: chart " + std::to_string(i) +
                                        " has no boundary loop (closed chart, no seam)");
        }
    }

    // Charts must carry texel density proportional to their 3D area before packing.
    scaleProportional(charts, data);

    // Big charts first: they are the hardest to fit, small ones fill the notches left behind.
    // Sort by face count descending (2003 `Cmp_Charts`, `legacy2003/chartification/catlaspacker.cpp:26-34`, `:165`).
    //    std::stable_sort preserves insertion (chart-id) order for equal counts.
    std::vector<HalfEdgeChart *> order;
    order.reserve(charts.size());
    for (HalfEdgeChart &ac : charts) {
        order.push_back(&ac);
    }
    std::stable_sort(order.begin(), order.end(),
                     [](const HalfEdgeChart *a, const HalfEdgeChart *b) {
                         return a->faces.size() > b->faces.size();
                     });

    int            res = atlasResolution;
    Horizon        horizon(res);
    ScratchBuffers scratch;
    scratch.lower.assign(static_cast<std::size_t>(res * 2), res - 1);
    scratch.upper.assign(static_cast<std::size_t>(res * 2), 0);

    // Shrink the pack scale until every chart fits.
    //    As 2003: ratio starts at 4.0 and decrements by 0.01 per attempt (`legacy2003/chartification/catlaspacker.cpp:153,171`).
    float ratio = 4.0f;
    bool  success = false;
    while (!success) {
        success = true;
        ratio -= 0.01f;
        float scale = ratio / std::sqrt(1.0f);
        horizon.reset();

        for (HalfEdgeChart *ac : order) {
            if (!packSingleChart(*ac, scale, /*write=*/false, horizon, res, scratch)) {
                success = false;
                break;
            }
        }

        if (ratio < 0.01f) {
            break; // safety: degenerate input never fits
        }
    }

    // The search passes were dry; this one writes the UVs at the winning ratio.
    horizon.reset();
    for (HalfEdgeChart *ac : order) {
        packSingleChart(*ac, ratio / std::sqrt(1.0f), /*write=*/true, horizon, res, scratch);
    }

    // Normalise to [0,1]², reserving the gutter as a border (kGutterPx in the
    //    res-column packing grid → kGutterPx/(res-1) in UV).
    normalizeAtlas(charts, static_cast<float>(kGutterPx) / static_cast<float>(res - 1));
}


/// Guarded-flatten every chart (Spec 19-03) then scale+pack into [0,1]² — shared by both public forms.
static std::vector<HalfEdgeChart> flattenAndPack(const HalfEdgeMesh &mesh, const MeshData &data,
                                                 ChartSet &charts, int atlasResolution,
                                                 const InjectivityOptions &opts, InjectivityStats *stats) {
    std::vector<HalfEdgeChart> atlasCharts = flattenCharts(mesh, data, charts, opts, stats);
    packAtlas(atlasCharts, data, atlasResolution);
    return atlasCharts;
}


std::vector<math::Vec2f>
parametrize(const HalfEdgeMesh       &mesh,
            const MeshData           &data,
            ChartSet                 &charts,
            int                       atlasResolution,
            const InjectivityOptions &opts,
            InjectivityStats         *stats) {
    const std::vector<HalfEdgeChart> atlasCharts =
        flattenAndPack(mesh, data, charts, atlasResolution, opts, stats);

    // Scatter per-local UV back to per-mesh-vertex UV.
    // LOSSY: a seam vertex (in multiple charts) keeps the last chart's UV.
    // parametrizeCorners() is the seam-correct form; this stays for the legacy contract.
    const int                numVerts = static_cast<int>(data.positions.size());
    std::vector<math::Vec2f> perVertexUV(numVerts, math::Vec2f(0.0f, 0.0f));
    for (const HalfEdgeChart &ac : atlasCharts) {
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


std::vector<math::Vec2f>
parametrizeCorners(const HalfEdgeMesh       &mesh,
                   const MeshData           &data,
                   ChartSet                 &charts,
                   int                       atlasResolution,
                   const InjectivityOptions &opts,
                   InjectivityStats         *stats) {
    const std::vector<HalfEdgeChart> atlasCharts =
        flattenAndPack(mesh, data, charts, atlasResolution, opts, stats);

    return cornerUVs(atlasCharts, static_cast<int>(data.indices.size()) / 3);
}


/// Write the corner UVs of every face chart @p ac owns (local face i of the chart is ac.faces[i]).
static void scatterChart(const HalfEdgeChart &ac, std::vector<math::Vec2f> &cornerUV) {
    for (int i = 0; i < static_cast<int>(ac.faces.size()); ++i) {
        const std::size_t base = 3 * static_cast<std::size_t>(ac.faces[static_cast<std::size_t>(i)]);
        for (int k = 0; k < 3; ++k) {
            cornerUV[base + static_cast<std::size_t>(k)] = ac.uv[static_cast<std::size_t>(localCorner(ac, i, k))];
        }
    }
}


std::vector<math::Vec2f> cornerUVs(const std::vector<HalfEdgeChart> &atlas, int numFaces) {
    std::vector<math::Vec2f> cornerUV(
        static_cast<std::size_t>(numFaces) * 3, math::Vec2f(0.0f, 0.0f));

    // Each face's three corners take UVs from its OWNING chart only, so no triangle ever straddles
    // a chart seam, and a mesh vertex the self-seam cut split resolves per corner through the
    // chart's local triangulation.
    for (const HalfEdgeChart &ac : atlas) {
        scatterChart(ac, cornerUV);
    }
    return cornerUV;
}

} // namespace decimate::chart
