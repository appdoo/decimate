// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file parametrize.h
/// @brief LSCM chart flattening + atlas packing — public surface for Spec 08-02.
///
/// @note The public signatures that take `HalfEdgeMesh` / expose `HalfEdgeChart` are
///       subject to change in 0.7 (flat `MeshData` overloads; half-edge becomes an
///       implementation detail).
///
/// Consumes the `ChartSet` from 08-01 (`chart.h`), flattens each chart to 2D
/// with least-squares conformal maps (LSCM, Lévy et al. 2002) via a sparse
/// Cholesky solve, scales each chart by surface-area ratio, and packs
/// the charts into one normalized [0,1]² atlas — the UV domain Spec 10 (bake)
/// raytraces into.  Derived from the 2003 `cconformalmapper` + `catlaspacker`; the solver (Cholesky for CG), the
/// silhouette and the injectivity guard differ from 2003.
///
/// Boundary contract: no GL/OS, no I/O, no global mutable state.

#include <vector>

#include <decimate/math/vec2.h>
#include <decimate/chart/chart.h>         // decimate::chart::Chart, ChartSet (Spec 08-01)
#include <decimate/mesh/half_edge_mesh.h> // decimate::HalfEdgeMesh
#include <decimate/types/mesh_data.h>     // decimate::MeshData

namespace decimate::chart {

/// @brief 2D point (u, v) used by the packer and the gutter-loop API.
struct Tex2 {
    float u = 0.0f;
    float v = 0.0f;
};


/// Atlas gutter in packing-grid pixels: how far each chart boundary is offset
/// outward so packed charts never touch (Greissl miter boundary offset, Spec 09-02).
inline constexpr int kGutterPx = 4;


/// Default miter limit: a corner never moves further than `limit * gutter`.
inline constexpr float kDefaultMiterLimit = 2.0f;


/// @brief Offset a closed boundary loop outward by `gutter` (miter join).
///
/// Each vertex moves along the bisector of its two edge normals by
/// `gutter / cos(half-angle)`, capped at `miterLimit * gutter` (sharp spikes are
/// bevelled, never extended).  Winding-independent: the outward side is derived
/// from the loop's signed area.  Pure float in/out; no pixel clamping (callers
/// clamp/rasterize).  `gutter == 0` returns the input unchanged.
///
/// The raw offset folds back over itself where a notch is narrower than 2 x `gutter` (reflex vertices
/// whose offset edges cross); `miterOffsetLoop` cuts those fold lobes out, so the result is a simple
/// polygon (`hasSelfCrossing` is false) whenever `loop` is.
///
/// @param loop        Ordered closed loop (last vertex connects to the first).
/// @param gutter      Outward offset distance, in the loop's own units (>= 0).
/// @param miterLimit  Cap factor (values below 1 are treated as 1).
/// @returns One offset vertex per input vertex (fewer only when fold lobes were cut out, plus the
///          crossing points); loops with < 3 vertices are returned as is.
[[nodiscard]] std::vector<Tex2> miterOffsetLoop(const std::vector<Tex2> &loop,
                                                float                    gutter,
                                                float                    miterLimit = kDefaultMiterLimit);

/// @brief True when two non-adjacent edges of the closed polyline @p loop cross properly (O(B^2)).
///        The same predicate as `HalfEdgeChart::boundarySimple`; false means @p loop is a simple polygon.
///
/// @param touchTol  Crossings where an end point lies within this distance of the other edge count as
///                  touches, not crossings (default 0 = exact).  A loop that was simple in exact arithmetic
///                  can touch itself where `miterOffsetLoop` cut a fold, and float rounding in a later
///                  transform (atlas placement, normalisation) turns such a touch into a ~1e-6 crossing.
[[nodiscard]] bool hasSelfCrossing(const std::vector<Tex2> &loop, float touchTol = 0.0f);


/// @brief Cover one segment in a per-column lower/upper silhouette, conservatively (packer primitive,
///        public for its unit test, like `miterOffsetLoop`).
///
/// For every pixel column `u` the segment touches, `lower[u]` becomes `min(lower[u], floor(v_min))`
/// and `upper[u]` becomes `max(upper[u], floor(v_max))`, where `v_min`/`v_max` are the segment's
/// extreme rows inside the closed interval `[u, u+1]`.  A steep edge therefore fills all its rows
/// in one column, and an end point exactly on a column boundary covers both columns when the
/// segment extends into the left one.  Rows are clamped to `[0, res-1]`.
///
/// @param x1,y1,x2,y2  End points in pixel units, each already inside `[0, res-1]`.
/// @param res          Grid resolution; `lower` and `upper` have at least `res` entries.
/// @param lower,upper  Per-column silhouette, updated in place (caller initialises: lower = res, upper = -1).
void coverSegment(float x1, float y1, float x2, float y2, int res,
                  std::vector<int> &lower, std::vector<int> &upper);

/// @brief One chart's flattened (pre-pack) result.
///
/// Carries an owned `HalfEdgeMesh` built from the chart's local triangulation,
/// so boundary and adjacency are derivable as queries rather than cached fields.
/// `localToMesh[li]` maps local vertex index → mesh vertex id (== the old
/// `vertexOfLocal`); many-to-one across a chart's seams, so never invert it (`localCorner`).
/// `faces` holds the mesh face ids (needed for 3D-area accounting in the packer).
///
/// ## Boundary contract
/// The boundary is obtained via `mesh.boundaryPointIndices()` — there is NO
/// stored `boundary` field on this struct (Spec 09-01 Boundary Contracts).
struct HalfEdgeChart {
    HalfEdgeMesh             mesh;        ///< Owned per-chart half-edge mesh (local, index-based).
    std::vector<math::Vec2f> uv;          ///< Per local vertex UV (size == mesh.pointCount()).
    std::vector<int>         localToMesh; ///< Local vertex index → mesh vertex id.
    std::vector<int>         faces;       ///< Mesh face ids (3D-area term; parallel to mesh faces).
    /// Placed gutter outline (miter-offset boundary loop) in atlas UV [0,1]^2, same
    /// transform as `uv`.  Filled by `packAtlas`; empty if the chart has no boundary loop.
    std::vector<math::Vec2f> gutterLoop;
    /// Flattening quality, measured by `flatten` in O(n + B^2), B = boundary edges: triangles whose UV orientation is
    /// opposite to the chart's (or degenerate), and whether the boundary edges cross each other.
    /// A disk chart with no flipped triangle and a simple boundary is injective (no fold).
    int  flippedTriangles = 0;
    bool boundarySimple = true;

    /// True when the flattening is a valid embedding: no flipped triangle and a simple boundary.
    [[nodiscard]] bool injective() const noexcept {
        return flippedTriangles == 0 && boundarySimple;
    }
};


/// @brief Local vertex at corner @p k of local face @p localFace of a flattened chart.
///
/// Corner `k` = 0, 1, 2 is the face's vertex a, b, c (the order of `MeshData::indices`).  Read
/// from the chart's own local triangulation, so a mesh vertex that the self-seam cut split
/// into several local vertices resolves to the one this face's corner uses.
///
/// This is the ONLY sanctioned way to go from a face corner to its local vertex / UV
/// (`uv[localCorner(ac, i, k)]`).  `localToMesh` is many-to-one on a cut chart and must not be
/// inverted.
///
/// @pre `0 <= localFace < ac.faces.size()` and `0 <= k < 3`.
[[nodiscard]] inline int localCorner(const HalfEdgeChart &ac, int localFace, int k) noexcept {
    // writeFace stores the corners as half-edges 3f+0 (end c), 3f+1 (end b), 3f+2 (end a).
    return static_cast<int>(ac.mesh.halfedge(static_cast<HeIndex>(3 * localFace + 2 - k)).end);
}


/// @brief LSCM-flatten a single chart.
///
/// The chart's local mesh is built by corner: corners are merged only across intra-chart edges
/// that are not in `Chart::seamHalfEdges`, so a chart that touches itself (a tube, the bunny
/// base plate) is cut open into a disk along its seams, as 2003 `CreateVertexStructure` did (`legacy2003/chartification/cchartification.cpp:395-407`).
/// `localToMesh` may therefore repeat a mesh vertex id.
///
/// Pins the two farthest boundary vertices to fix the gauge freedom, assembles
/// the per-triangle Cauchy-Riemann rows, solves the normal equations by sparse Cholesky,
/// and returns per-local-vertex UVs (pre-pack, pre-scale).
///
/// The returned `HalfEdgeChart::mesh` is the owned local triangulation;
/// `boundaryPointIndices()` on it yields the chart's outline.
///
/// @param mesh   Half-edge connectivity for the whole mesh.
/// @param data   Vertex positions (for the LSCM 3D geometry).
/// @param chart  The chart to flatten (its face list).
/// @returns The flattened chart; `uv` finite for a non-degenerate chart.
[[nodiscard]] HalfEdgeChart flatten(const HalfEdgeMesh &mesh,
                                    const MeshData     &data,
                                    const Chart        &chart);

/// @brief Number of triangle pairs of a flattened chart whose UV interiors overlap.
///
/// The LSCM map of a chart is not guaranteed injective: it can fold over itself
/// (the defect the 2003 thesis names, `diploma_tex/chartification.tex:98`).  A pair
/// of local triangles counts when a separating-axis test over the six edge normals
/// finds no axis with a gap >= -eps, `eps = 1e-6 * (UV bounding-box diagonal)`;
/// triangles that merely touch (shared edge or vertex) are not overlaps, and
/// triangles with |area| <= 1e-12 * diag^2 are skipped.  Corners are read from the
/// chart's local triangulation (`ac.mesh`), local face `i` being `ac.faces[i]`.
/// Candidate pairs come from a uniform grid over the UV box, each pair counted once.
///
/// @param ac  A flattened chart (pre-pack or packed; the test is scale invariant).
/// @returns   0 for an injective chart; > 0 means the chart folds.
[[nodiscard]] int selfOverlapPairs(const HalfEdgeChart &ac);


/// @brief Where `flattenBijective` starts its descent.
enum class BijectiveStart {
    Tutte,   ///< Tutte embedding: bijective for every disk chart (the paper's start).
    Current, ///< The chart's current `uv`, rescaled to the 3D area; used only when it is already
             ///< `injective()`, else Tutte.  Any bijective start is valid for the interior-point
             ///< method; LSCM is far closer to the optimum than Tutte.
};


/// @brief Knobs of `flattenBijective` (Smith & Schaefer 2015, section 3.3).
struct BijectiveOptions {
    BijectiveStart start = BijectiveStart::Tutte; ///< Starting map; Tutte keeps direct callers unchanged.
    int            maxIterations = 3000;          ///< L-BFGS iterations (each one keeps the map bijective).
    double         tolerance = 1e-6;              ///< Stop when the last 10 iterations together lowered the energy by less than this fraction.
};


/// @brief What `flattenBijective` measured.
struct BijectiveStats {
    int    iterations = 0;  ///< L-BFGS iterations run.
    double energy = 0.0;    ///< Isometric error E_D per unit 3D area (paper Table 1 "avg"; 4 = isometric).
    double maxEnergy = 0.0; ///< Largest per-triangle isometric error (paper Table 1 "max").
    double seconds = 0.0;   ///< Wall time of the call.
};


/// @brief Tutte embedding of a disk chart: boundary on the unit circle by 3D arc length, interior
///        by mean-value weights (Floater 2003), scaled to the chart's 3D area.  Bijective by Tutte's
///        theorem — the paper's starting point (section 3.3, Figure 1 left).
///
/// Rewrites `ac.uv` and the quality fields; `ac.mesh` / `localToMesh` / `faces` are read only.
/// @returns false (and leaves `ac` untouched) when the local mesh is not a topological disk
///          (closed, a hole, or several components) — then the paper's section-4 hole filling
///          would be needed, which this port does not do.
bool flattenTutte(const MeshData &data, HalfEdgeChart &ac);


/// @brief Re-flatten a disk chart so the map is bijective and the isometric distortion low
///        (Smith & Schaefer 2015: Tutte start, L-BFGS on E_D + E_B with the flip-avoiding line
///        search, free boundary kept from self-intersecting by the barrier term).
///
/// Ignores the current `ac.uv` unless `opts.start == Current` and it is injective (an LSCM fold is
/// no valid start); overwrites it either way.  The result is
/// `injective()` whenever the function returns true.  Deterministic.
/// @returns false when the chart is not a disk (see `flattenTutte`); `ac` is then untouched.
bool flattenBijective(const MeshData &data, HalfEdgeChart &ac,
                      const BijectiveOptions &opts = {}, BijectiveStats *stats = nullptr);

/// @brief Isometric distortion of a flattened chart (paper Eq. 2 integrated over the 3D area,
///        Eq. 3): `sum_T A_T (1 + A_T^2/A_uv^2) ||J_T||_F^2 / sum_T A_T`.  4 is a perfect
///        isometry; +inf when a UV triangle is degenerate or flipped.  Not scale invariant.
/// @param maxTriangle  Optional out: the largest per-triangle value.
[[nodiscard]] double isometricEnergy(const MeshData &data, const HalfEdgeChart &ac,
                                     double *maxTriangle = nullptr);

/// @brief Sander et al. 2001 stretch of one flattened chart, as `packAtlas` will scale it.
///
/// Per triangle, Gamma >= gamma are the singular values of the uv -> 3D map (Gamma = 1 everywhere is an
/// isometry).  `packAtlas` multiplies every chart by `sqrt(area3D / areaUV)`, so the values are
/// normalised by the chart's OWN area ratio `k = sqrt(sum A_uv / sum A_3D)`:
/// `l2 = k sqrt(sum_T A_T (Gamma^2 + gamma^2)/2 / sum_T A_T)` (A_T = 3D area), `linf = k max_T Gamma`.
/// Flipped or degenerate UV triangles are skipped (the fold guard handles those), so the numbers are
/// only meaningful for an `injective()` chart.  Unit-free and scale invariant.
struct ChartStretch {
    double l2 = 0.0;       ///< Area-weighted RMS stretch; 0 for a chart with no measurable triangle.
    double linf = 0.0;     ///< Worst-triangle stretch (max Gamma).
    int    worstFace = -1; ///< Local face index of that triangle (ties: lowest), -1 if none.
};


/// @brief Measure the stretch of a flattened chart (O(faces), no solve).
/// @param perFace  Optional out: normalised Gamma per local face, -1 for a skipped triangle.
[[nodiscard]] ChartStretch chartStretch(const MeshData &data, const HalfEdgeChart &ac,
                                        std::vector<double> *perFace = nullptr);


/// @brief Strategy of `flattenCharts`: which flattener runs on a chart and what happens to a chart
///        that stays folded.  One enum, no flag matrix.
///
/// Every mode starts by LSCM-flattening and testing each chart.  `BijectiveOptions::start` is chosen
/// per mode and documented here, never silently overridden:
/// - `Bijective` always starts `flattenBijective` from `BijectiveStart::Current` (the LSCM map when it
///   is injective, otherwise `flattenBijective` itself falls back to Tutte); the caller's
///   `InjectivityOptions::bijective.start` is ignored in this mode.
/// - `Lscm` repairs only a folded chart and uses the caller's `bijective.start` (default Tutte,
///   because an LSCM fold is no valid start).
enum class Flattener {
    /// Default.  Every disk chart is re-flattened by `flattenBijective`; a folded chart it refuses
    /// (not a disk) is split at the fold.
    Bijective,
    /// LSCM; a folded chart is re-flattened by `flattenBijective`, and split when that refuses it.
    Lscm,
    /// LSCM; every folded chart is split at the fold (the Spec 19-03 path), no bijective step.
    LscmSplit,
    /// LSCM, flatten and detect only: the `ChartSet` is not touched, nothing is repaired.
    LscmRaw,
};


/// @brief Knobs of the injectivity guard in `flattenCharts`.
struct InjectivityOptions {
    Flattener        flattener = Flattener::Bijective; ///< Strategy, see `Flattener`.
    int              maxRounds = 8;                    ///< Upper bound on detect-and-split rounds; leftovers are reported, never thrown.
    BijectiveOptions bijective;                        ///< Knobs of `flattenBijective` (`start` only honoured by `Flattener::Lscm`).

    /// Stretch bound (Spec 20): a chart whose Sander L2 stretch exceeds this is cut open inside itself
    /// (see `flattenCharts`).  `0` = off (the 2003 behaviour, every golden unchanged).  See `ChartStretch`
    /// for the metric; 1 is an isometry.
    double maxStretchL2 = 0.0;
    /// Same for the chart's worst triangle (L-infinity, max Gamma).  `0` = off.
    double maxStretchLinf = 0.0;
};


/// @brief What the injectivity guard found and did.
struct InjectivityStats {
    int    selfOverlappingBefore = 0; ///< Charts folding after the first flattening.
    int    selfOverlappingAfter = 0;  ///< Charts still folding when the guard stopped.
    int    splits = 0;                ///< Charts added by `splitChart` (== growth of ChartSet::charts).
    int    rounds = 0;                ///< Detect-and-split rounds run (0 if nothing folded).
    int    repairedBijectively = 0;   ///< Folded charts `flattenBijective` made injective.
    int    bijectiveFailed = 0;       ///< Folded charts it refused (non-disk) — these went to the split.
    double bijectiveSeconds = 0.0;    ///< Wall time spent in `flattenBijective`.
    int    flattenedBijectively = 0;  ///< Disk charts `flattenBijective` flattened (Flattener::Bijective).
    /// Charts with `!HalfEdgeChart::injective()` when the guard stopped.  Unlike `selfOverlappingAfter`
    /// (SAT overlap pairs only) this also counts flipped triangles, a crossing boundary and charts a
    /// strategy never touched (`LscmRaw`), so a caller can trust "0" to mean an injective atlas.
    int nonInjectiveAfter = 0;
    int    closedOpened = 0;   ///< Closed (boundary-less) charts opened by an initial cut; with the bound on or off.
    int    cuts = 0;           ///< Cut rounds kept by the stretch bound (each adds up to four paths).
    double maxChartL2 = 0.0;   ///< Largest per-chart L2 stretch of the final charts (measured with the bound on or off).
    double maxChartLinf = 0.0; ///< Largest per-chart L-infinity stretch of the final charts.
};


/// @brief Flatten every chart of @p charts, repairing self-overlapping charts.
///
/// The single entry for building an atlas: all `parametrize*` forms and the viewer go
/// through it.  Each chart is LSCM-flattened and tested (`injective()`, then `selfOverlapPairs`).
/// What happens next is the `Flattener` strategy: `Bijective` (default) re-flattens every disk chart
/// by `flattenBijective`; `Lscm` only a folding one; `LscmSplit` and `LscmRaw` never use it.  A
/// folding chart that is still folded (`flattenBijective` refuses a non-disk; or `LscmSplit`) is
/// split by `splitChart` at the involved face farthest (dual BFS inside the chart) from the chart's
/// seed, ties to the lowest face index; only the two pieces are re-flattened and re-tested
/// (bijectively again when the strategy uses it).  At most `opts.maxRounds` split rounds run;
/// charts that still fold are only counted in `stats->selfOverlappingAfter` / `nonInjectiveAfter`.
///
/// Stretch bound (Spec 20, off by default; Gu, Gortler, Hoppe 2002, section 3): a chart whose LSCM map
/// breaks `opts.maxStretchL2` / `maxStretchLinf` (a bound of 0 is off) gets up to four vertex-disjoint
/// shortest interior paths added to `Chart::seamHalfEdges`, from its most-stretched faces (worst first, from
/// the face corner farthest from the boundary) to the boundary, and is re-flattened.  Rounds repeat until the
/// bound is met, no path is left, or a round does not lower the chart's L2 stretch (that round is undone).
/// The heal above then runs once on the final charts.  Charts are never split for stretch; an unmet bound is
/// reported through `stats->maxChartL2` / `maxChartLinf`.  With both bounds 0 nothing here runs.
///
/// Closed charts (no boundary, no seam; Spec 24 H4) are opened first, bound or not, so they flatten and
/// pack: genus 0 along the shortest path between the two mutually farthest vertices (at least two edges),
/// genus > 0 along the pruned complement of a dual spanning tree from the chart's seed face (Gu et al.
/// 2002, section 3.1); `stats->closedOpened` counts them.  Charts with a boundary are never touched by this
/// step, so every existing result is unchanged.  `Flattener::LscmRaw` is detect-only and opens nothing.
/// Deterministic.
///
/// @param charts  Updated IN PLACE when splits occur: appended charts, rewritten
///                `faceToChart` and `faces` lists.  Untouched under `Flattener::LscmRaw`.
/// @param opts    Guard knobs; default ON.
/// @param stats   Optional out-parameter.
/// @returns One flattened (unpacked) chart per entry of the updated `charts.charts`.
[[nodiscard]] std::vector<HalfEdgeChart>
flattenCharts(const HalfEdgeMesh &mesh, const MeshData &data, ChartSet &charts,
              const InjectivityOptions &opts = {}, InjectivityStats *stats = nullptr);

/// @brief Scale + 16-rotation horizon-pack a set of flattened charts into [0,1]².
///
/// Applies `sqrt(area3D / areaUV)` per chart (area preservation), then places
/// each chart (largest first) trying 16 rotations against a column horizon,
/// then normalizes the whole layout into the unit square.  Mutates the charts'
/// `uv` in place to their final atlas positions.
///
/// @param charts          Flattened charts (uv updated in place to atlas coords).
/// @param data            Vertex positions (for the 3D area term).
/// @param atlasResolution Packing grid resolution (default 256), at least 2.
/// @throws std::invalid_argument if `atlasResolution < 2`, or if a chart with faces has no
///         boundary loop (a closed chart without a seam cannot be packed; the message names it).
void packAtlas(std::vector<HalfEdgeChart> &charts,
               const MeshData             &data,
               int                         atlasResolution = 256);

/// @brief Full pipeline: flatten every chart (injectivity-guarded), scale, pack, emit per-vertex UVs.
///
/// @param mesh            Half-edge connectivity.
/// @param data            Vertex positions.
/// @param charts          The segmentation from 08-01; updated in place by the guard (see `flattenCharts`).
/// @param atlasResolution Packing resolution (default 256).
/// @param opts            Injectivity guard knobs (default ON).
/// @param stats           Optional guard report.
/// @returns Per-mesh-vertex UV in [0,1]² (size == data.positions.size()).
///
/// @warning LOSSY for a seamed atlas: a vertex shared by N charts can only hold
///          ONE UV, so the last chart written wins and triangles straddling a
///          chart seam get inconsistent corner UVs.  Use `parametrizeCorners`
///          for a correct atlas; this overload is kept for the legacy contract.
[[nodiscard]] std::vector<math::Vec2f>
parametrize(const HalfEdgeMesh &mesh, const MeshData &data,
            ChartSet &charts, int atlasResolution = 256,
            const InjectivityOptions &opts = {}, InjectivityStats *stats = nullptr);

/// @brief Full pipeline emitting per-FACE-CORNER UVs (the seam-correct form).
///
/// Returns `3 * faceCount` UVs: entry `3*f + k` is the atlas UV of corner `k`
/// of face `f`.  Because every face belongs to exactly one chart, all three
/// corners of a face come from that one chart's flattening — so no triangle
/// ever straddles a chart seam (the defect of the per-vertex form).  This is
/// the standard "wedge UV" representation a texture atlas and the bake (Spec 10)
/// consume.
///
/// @param mesh            Half-edge connectivity.
/// @param data            Vertex positions.
/// @param charts          The segmentation from 08-01 (uses `faceToChart`); updated in place by
///                        the injectivity guard (see `flattenCharts`).
/// @param atlasResolution Packing resolution (default 256).
/// @param opts            Injectivity guard knobs (default ON).
/// @param stats           Optional guard report.
/// @returns Per-corner UVs in [0,1]², size `3 * (data.indices.size()/3)`.
[[nodiscard]] std::vector<math::Vec2f>
parametrizeCorners(const HalfEdgeMesh &mesh, const MeshData &data,
                   ChartSet &charts, int atlasResolution = 256,
                   const InjectivityOptions &opts = {}, InjectivityStats *stats = nullptr);

/// @brief Scatter a packed atlas into per-FACE-CORNER UVs (the seam-correct wedge form).
///
/// Entry `3*f + k` is the UV of corner `k` of face `f`, read from the ONE chart that owns face
/// `f` (`HalfEdgeChart::faces`), so a mesh vertex the cut split resolves per corner.  This is the
/// scatter `parametrizeCorners` ends with; callers that keep the packed charts (for the gutter or
/// the outlines) use it to get the same UVs without flattening twice.  Faces owned by no chart
/// keep (0, 0).
///
/// @param atlas     Flattened + packed charts (`flattenCharts` then `packAtlas`).
/// @param numFaces  Face count of the mesh; the result has `3 * numFaces` entries.
[[nodiscard]] std::vector<math::Vec2f>
cornerUVs(const std::vector<HalfEdgeChart> &atlas, int numFaces);

} // namespace decimate::chart
