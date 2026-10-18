// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file normal_bake.h
/// @brief Raytraced normal bake — low←high (Spec 10 Step 2).
///
/// `bake()` is the modernized `normalmap::mapnormals` from
/// `legacy/split/texturesample.cpp:339-383`.
///
/// ## Algorithm summary (three passes)
///
/// 1. **Clear**: fill the texel buffer with the background sentinel (128,128,128).
/// 2. **Rasterize + raytrace**: for each low-poly face, rasterize its UV
///    triangle into the atlas, interpolate the low-poly ray origin/direction,
///    cast the ray against every high-poly face, range-compress the hit normal
///    (or fall back to the low normal on a miss) into the texel.
/// 3. **Dilate**: 6-pass dilation fills background texels from their filled
///    neighbours (seam bleed / gutter fill).
///
/// ## Seam-correct gutter fill (`bakeSeamFill`)
///
/// @note `bakeSeamFill` takes a `HalfEdgeMesh`; that signature is subject to change
///       in 0.7 (flat `MeshData` overload; half-edge becomes an implementation detail).
///
/// A second entry point (`bakeSeamFill`) replaces the blur dilate at chart
/// seams with a cross-chart reference-map sample (Andor 2024, Ray et al. 2010
/// "invisible seams" indirection map — NOT the MS patent US10,839,590 polygon-
/// projection mechanism).
///
/// Pass 1: identical interior bake (bake() with dilatePasses=0).
/// Pass 2: seam-correct gutter fill — each gutter texel just outside a
///         chart-seam edge samples the frozen pass-1 atlas at the
///         corresponding point inside the neighbour chart (t→1-t param reversal,
///         perpendicular flip).  Open-mesh boundary edges fall back to classic
///         blur dilate (same as bake()).
///
/// ## Boundary contracts
///
/// - No file I/O — returns the `NormalMap` in memory.
/// - No GL / GLUT / windows.h / global state.
/// - Does NOT reach into `legacy/`.
///
/// ## getTestFaces
///
/// The legacy 6-plane prism cull (`getTestFaces`) is implemented (faithful port):
/// per low face, the high mesh is culled to the faces inside the low triangle's
/// 1.5×-scaled swept prism, and only those are raytraced.  This is the dense→
/// sparse normal-transfer optimisation, not a general RT accel.  Output is
/// identical when high is not denser than low; the speedup applies when high ≫
/// low.  (A BVH/SAH accel for the per-face ray test is a separate Spec-13 item.)

#include <vector>

#include <decimate/math/vec2.h>
#include <decimate/math/vec3.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/bake/normal_map.h>

namespace decimate::bake {

/// @brief Raytrace the high-poly normals onto the low-poly UV atlas.
///
/// @param low          Reduced (low-poly) mesh — positions + indices.
/// @param lowCornerUV  Per-corner atlas UVs in [0,1]², size `3*lowFaceCount`
///                     (entry `3*f+k` = corner `k` of face `f`). From
///                     `decimate::chart::parametrizeCorners`.
/// @param lowNormals   Per-vertex normals on the low mesh (size == nVerts_low).
///                     From `computeVertexNormals(low, EqualUnit)`.
/// @param high         Full-detail (high-poly) mesh — positions + indices.
/// @param highNormals  Per-vertex normals on the high mesh (size == nVerts_high).
///                     From `computeVertexNormals(high, EqualUnit)`.
/// @param texSize      Texture resolution (e.g. 256).  Must be > 0.
/// @param dilatePasses Opt-in gutter extension: pixels to grow each chart's
///                     baked region outward (8-neighbour fill, one px/pass) to
///                     avoid seams under bilinear filtering on the mesh.
///                     DEFAULT 0 = the original 2003 output (no dilation; no
///                     normals outside a chart). When > 0 the atlas packer must
///                     reserve matching clearance, else skirts bleed between
///                     charts (see packAtlas's gutter).
/// @returns            `NormalMap` of size `texSize × texSize` containing the
///                     range-encoded normals; background texels are (128,128,128).
///
/// @pre `lowCornerUV.size() == 3 * (low.indices.size() / 3)`.
/// @pre `lowNormals.size()  == low.positions.size()`.
/// @pre `highNormals.size() == high.positions.size()`.
[[nodiscard]] NormalMap bake(const decimate::MeshData&                   low,
                             const std::vector<decimate::math::Vec2<float>>&   lowCornerUV,
                             const std::vector<decimate::math::Vec3<float>>&   lowNormals,
                             const decimate::MeshData&                   high,
                             const std::vector<decimate::math::Vec3<float>>&   highNormals,
                             int                                      texSize,
                             int                                      dilatePasses = 0,
                             bool                                     tangentSpace = false);

/// @brief Bake + seam-correct gutter fill (2-pass, cross-chart sampling).
///
/// Identical to `bake(dilatePasses=0)` for the interior (Pass 1).  Pass 2
/// fills each gutter texel by sampling the neighbour chart's pass-1 value
/// across the seam edge, using the cross-chart reference map derived from
/// `hem.sym()` and `cornerUV`.
///
/// Open-mesh boundary edges (sym == kNoHalfedge) fall back to classic blur
/// dilate (dilatePasses controls how many blur passes after the seam fill).
///
/// ## Seam detection
///
/// A mesh edge is a CHART SEAM iff the two faces sharing it have different
/// atlas UVs at their shared vertices (the same 3D edge sits at two distinct
/// atlas locations). The refMap is built at bake time from `hem.sym()` +
/// `cornerUV`; no change to the parametrize API is needed.
///
/// ## Correctness invariants
///
/// 1. Param reversal: t along e_a maps to (1-t) along e_b (sym runs opposite).
/// 2. Perpendicular flip: depth d OUTSIDE chart A = depth d INSIDE chart B.
///    Sample at `pointOnEdge(e_b, 1-t) + inwardNormal(e_b) * d`.
///
/// @param low          Low-poly mesh.
/// @param hem          Half-edge mesh built from `low` (provides sym()).
/// @param lowCornerUV  Per-corner atlas UVs in [0,1]², size 3*faceCount.
/// @param faceToChart  Per-face chart index, size faceCount (from ChartSet).
/// @param lowNormals   Per-vertex normals on the low mesh.
/// @param high         High-poly mesh.
/// @param highNormals  Per-vertex normals on the high mesh.
/// @param texSize      Texture resolution.
/// @param dilatePasses Blur-dilate passes applied after seam fill (open edges).
/// @returns            NormalMap with seam-correct gutter.
///
/// @pre `lowCornerUV.size() == 3 * (low.indices.size() / 3)`.
/// @pre `faceToChart.size() == low.indices.size() / 3`.
[[nodiscard]] NormalMap bakeSeamFill(
    const decimate::MeshData&                   low,
    const decimate::HalfEdgeMesh&               hem,
    const std::vector<decimate::math::Vec2<float>>&   lowCornerUV,
    const std::vector<int>&                 faceToChart,
    const std::vector<decimate::math::Vec3<float>>&   lowNormals,
    const decimate::MeshData&                   high,
    const std::vector<decimate::math::Vec3<float>>&   highNormals,
    int                                      texSize,
    int                                      dilatePasses = 2,
    bool                                     tangentSpace = false);

}  // namespace decimate::bake
