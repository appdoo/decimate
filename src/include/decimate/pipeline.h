// This file is part of the decimate library; see decimate.h for version/license details.
/// @file pipeline.h
/// @brief The two-call pipeline: `chartify` (charts + packed UVs) and `bakeAtlas` (normal map + seam metric).
///
/// Core only: no GL, no printing. Stats and timings are returned; the caller logs them.
#pragma once

#include <array>
#include <vector>

#include <decimate/bake/normal_map.h>
#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/math/vec2.h>
#include <decimate/math/vec4.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/seam/straighten.h>
#include <decimate/types/mesh_data.h>

namespace decimate {

/// @brief Knobs of `chartify`. The defaults are the 2003 behaviour.
struct ChartOptions {
    int                       maxCharts = 64;        ///< Segmenter cap. A stretch bound usually wants 1 (it cuts the rest).
    chart::InjectivityOptions injectivity;           ///< Flattener and Spec 20 stretch bound; default bound off.
    int                       atlasResolution = 256; ///< Packing grid, at least 2.
    bool                      straighten = false;    ///< Seam straightening (Liu et al. 2017, section 6) after packing.
};


/// @brief Everything `chartify` produces: the charted mesh, its charts, the packed atlas.
struct Atlas {
    MeshData                 mesh;   ///< Owned copy of the input.
    HalfEdgeMesh             hem;    ///< Connectivity of `mesh`.
    chart::ChartSet          charts; ///< Segmentation after the injectivity guard (may hold more charts than asked).
    std::vector<math::Vec2f> uv;     ///< Packed atlas UV per corner (3 per face).
    /// MikkTSpace tangent per corner (xyz = T, w = bitangent sign) of `uv` and the EqualUnit vertex normals
    /// of `mesh`: the basis the tangent-space bake writes in, ready to ship as glTF `TANGENT`.
    std::vector<math::Vec4f> tangents;
    /// The packed charts `uv` was scattered from. Valid only when `ChartOptions::straighten` was false:
    /// straightening rewrites `uv` and leaves these stale.
    std::vector<chart::HalfEdgeChart> packed;

    chart::InjectivityStats stats;
    seam::StraightenStats   straightenStats;                                              ///< All zero unless `straighten` ran.
    double                  hemSeconds = 0.0, segmentSeconds = 0.0, flattenSeconds = 0.0; ///< Per stage; flatten includes pack + straighten.
};


/// @brief Segment, flatten, pack (and optionally straighten) `low` into an atlas.
/// @pre `low` has at least one face. An empty mesh yields an empty Atlas.
[[nodiscard]] Atlas chartify(const MeshData &low, const ChartOptions &opt = {});


enum class BakeSpace { Object,
                       Tangent };

/// @brief How the gutter around the charts is made seam-safe.
enum class Gutter {
    Seamless, ///< The one that hides the seams: Erased, in both spaces.
    Dilate,   ///< Plain blur dilation.
    SeamFill, ///< Real neighbour-chart texels in the gutter.
    Erased,   ///< Seam erasure on top of the dilate bake; in tangent space it is frame-aware (`atlas.tangents`).
};

struct BakeOptions {
    BakeSpace space = BakeSpace::Object;
    Gutter    gutter = Gutter::Seamless;
    int       resolution = 256;
    int       dilate = 2; ///< Blur passes of the bake.
};


/// @brief The baked map and how it came out.
struct BakeResult {
    bake::NormalMap       map;
    Gutter                gutterUsed = Gutter::Dilate; ///< Never Seamless: what was actually built (see `Gutter`).
    std::array<double, 3> dTotal {};                   ///< Seam mismatch D_total of `map` per channel (paper eq. 9).
    std::array<double, 3> dDilate {};                  ///< Same for the plain dilate map; equals `dTotal` unless Erased.
    int                   unknownTexels = 0;           ///< Erased: texels in the solve; else 0.
    double                rescale = 1.0;               ///< Erased, tangent space: factor that kept the solution in [-1, 1].
};


/// @brief Bake the normals of `high` onto `atlas.mesh` (vertex normals EqualUnit on both meshes).
/// If erasure throws, the dilate map is returned and `gutterUsed` says Dilate.
[[nodiscard]] BakeResult bakeAtlas(const Atlas &atlas, const MeshData &high, const BakeOptions &opt = {});

} // namespace decimate
