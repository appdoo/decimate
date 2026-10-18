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
/// with least-squares conformal maps (LSCM, Lévy et al. 2002) via the
/// `ConjugateGradientSolver`, scales each chart by surface-area ratio, and packs
/// the charts into one normalized [0,1]² atlas — the UV domain Spec 10 (bake)
/// raytraces into.  Behavior-preserving port of `cconformalmapper` + `catlaspacker`.
///
/// Boundary contract: no GL/OS, no I/O, no global mutable state.

#include <vector>

#include <decimate/math/vec2.h>
#include <decimate/chart/chart.h>          // decimate::chart::Chart, ChartSet (Spec 08-01)
#include <decimate/mesh/half_edge_mesh.h> // decimate::HalfEdgeMesh
#include <decimate/types/mesh_data.h>      // decimate::MeshData

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
/// @param loop        Ordered closed loop (last vertex connects to the first).
/// @param gutter      Outward offset distance, in the loop's own units (>= 0).
/// @param miterLimit  Cap factor (values below 1 are treated as 1).
/// @returns Same vertex count as `loop`; loops with < 3 vertices are returned as is.
[[nodiscard]] std::vector<Tex2> miterOffsetLoop(const std::vector<Tex2>& loop,
                                                float gutter,
                                                float miterLimit = kDefaultMiterLimit);

/// @brief One chart's flattened (pre-pack) result.
///
/// Carries an owned `HalfEdgeMesh` built from the chart's local triangulation,
/// so boundary and adjacency are derivable as queries rather than cached fields.
/// `localToMesh[li]` maps local vertex index → mesh vertex id (== the old
/// `vertexOfLocal`).  `faces` holds the mesh face ids (needed for 3D-area
/// accounting in the packer).
///
/// ## Boundary contract
/// The boundary is obtained via `mesh.boundaryPointIndices()` — there is NO
/// stored `boundary` field on this struct (Spec 09-01 Boundary Contracts).
struct HalfEdgeChart {
    HalfEdgeMesh                   mesh;         ///< Owned per-chart half-edge mesh (local, index-based).
    std::vector<decimate::math::Vec2<float>> uv;           ///< Per local vertex UV (size == mesh.pointCount()).
    std::vector<int>               localToMesh;  ///< Local vertex index → mesh vertex id.
    std::vector<int>               faces;        ///< Mesh face ids (3D-area term; parallel to mesh faces).
    /// Placed gutter outline (miter-offset boundary loop) in atlas UV [0,1]^2, same
    /// transform as `uv`.  Filled by `packAtlas`; empty if the chart has no boundary loop.
    std::vector<decimate::math::Vec2<float>> gutterLoop;
};

/// @brief LSCM-flatten a single chart.
///
/// Pins the two farthest boundary vertices to fix the gauge freedom, assembles
/// the per-triangle Cauchy-Riemann rows, solves the normal equations via CG,
/// and returns per-local-vertex UVs (pre-pack, pre-scale).
///
/// The returned `HalfEdgeChart::mesh` is the owned local triangulation;
/// `boundaryPointIndices()` on it yields the chart's outline.
///
/// @param mesh   Half-edge connectivity for the whole mesh.
/// @param data   Vertex positions (for the LSCM 3D geometry).
/// @param chart  The chart to flatten (its face list).
/// @returns The flattened chart; `uv` finite for a non-degenerate chart.
[[nodiscard]] HalfEdgeChart flatten(const HalfEdgeMesh& mesh,
                                    const MeshData&     data,
                                    const Chart&        chart);

/// @brief Scale + 16-rotation horizon-pack a set of flattened charts into [0,1]².
///
/// Applies `sqrt(area3D / areaUV)` per chart (area preservation), then places
/// each chart (largest first) trying 16 rotations against a column horizon,
/// then normalizes the whole layout into the unit square.  Mutates the charts'
/// `uv` in place to their final atlas positions.
///
/// @param charts          Flattened charts (uv updated in place to atlas coords).
/// @param data            Vertex positions (for the 3D area term).
/// @param atlasResolution Packing grid resolution (legacy hardcode 256).
void packAtlas(std::vector<HalfEdgeChart>& charts,
               const MeshData&             data,
               int                         atlasResolution = 256);

/// @brief Full pipeline: flatten every chart, scale, pack, emit per-vertex UVs.
///
/// @param mesh            Half-edge connectivity.
/// @param data            Vertex positions.
/// @param charts          The segmentation from 08-01.
/// @param atlasResolution Packing resolution (default 256).
/// @returns Per-mesh-vertex UV in [0,1]² (size == data.positions.size()).
///
/// @warning LOSSY for a seamed atlas: a vertex shared by N charts can only hold
///          ONE UV, so the last chart written wins and triangles straddling a
///          chart seam get inconsistent corner UVs.  Use `parametrizeCorners`
///          for a correct atlas; this overload is kept for the legacy contract.
[[nodiscard]] std::vector<decimate::math::Vec2<float>>
parametrize(const HalfEdgeMesh& mesh, const MeshData& data,
            const ChartSet& charts, int atlasResolution = 256);

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
/// @param charts          The segmentation from 08-01 (uses `faceToChart`).
/// @param atlasResolution Packing resolution (default 256).
/// @returns Per-corner UVs in [0,1]², size `3 * (data.indices.size()/3)`.
[[nodiscard]] std::vector<decimate::math::Vec2<float>>
parametrizeCorners(const HalfEdgeMesh& mesh, const MeshData& data,
                   const ChartSet& charts, int atlasResolution = 256);

}  // namespace decimate::chart
