// This file is part of the decimate library; see decimate.h for version/license details.
#include <decimate/pipeline.h>

#include <chrono>
#include <exception>
#include <utility>

#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normals.h>
#include <decimate/seam/erase.h>
#include <decimate/seam/seam_energy.h>
#include <decimate/seam/seam_topology.h>

namespace decimate {

static double secondsSince(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}


Atlas chartify(const MeshData &low, const ChartOptions &opt) {
    Atlas a;
    a.mesh = low;
    const int numPoints = static_cast<int>(low.positions.size());
    const int numFaces = static_cast<int>(low.indices.size() / 3);
    if (numFaces <= 0 || numPoints <= 0) {
        return a;
    }

    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    const auto *faces = reinterpret_cast<const Face *>(low.indices.data()); // flat index triples are Face[]

    auto t = std::chrono::steady_clock::now();
    a.hem = HalfEdgeMesh::create(numPoints, numFaces, faces);
    a.hemSeconds = secondsSince(t);

    t = std::chrono::steady_clock::now();
    a.charts = chart::segment(a.hem, low, opt.maxCharts);
    a.segmentSeconds = secondsSince(t);

    // Flatten and pack ONCE: the UV scatter, the straightening and any gutter reader all use this result.
    t = std::chrono::steady_clock::now();
    a.packed = chart::flattenCharts(a.hem, low, a.charts, opt.injectivity, &a.stats);
    chart::packAtlas(a.packed, low, opt.atlasResolution);
    a.uv = chart::cornerUVs(a.packed, numFaces);

    // Straightening rewrites `uv`, so it runs after the scatter; `packed` then no longer matches.
    if (opt.straighten) {
        auto res = seam::straighten(seam::weld(low, a.uv));
        a.straightenStats = res.stats;
        a.uv = res.mesh.cornerUV();
    }
    a.flattenSeconds = secondsSince(t);
    return a;
}


/// Seamless is a decision, not a mode: object space erases, tangent space (seams not well-defined,
/// Liu et al. 2017 footnote 3) fills from the neighbour chart. Erased in tangent space falls back the same way.
static Gutter resolveGutter(Gutter g, BakeSpace space) {
    if (g == Gutter::Seamless) {
        g = Gutter::Erased;
    }
    return g == Gutter::Erased && space == BakeSpace::Tangent ? Gutter::SeamFill : g;
}


BakeResult bakeAtlas(const Atlas &atlas, const MeshData &high, const BakeOptions &opt) {
    // Low and high poly use the same vertex-normal weighting so the bake compares like with like.
    const auto lowN = bake::computeVertexNormals(atlas.mesh, bake::Weighting::EqualUnit);
    const auto highN = bake::computeVertexNormals(high, bake::Weighting::EqualUnit);
    const bool tangent = opt.space == BakeSpace::Tangent;

    // D_total (the seam metric) is defined on the welded textured mesh, so every gutter builds it.
    const auto tm = seam::weld(atlas.mesh, atlas.uv);

    BakeResult r;
    r.gutterUsed = resolveGutter(opt.gutter, opt.space);
    if (r.gutterUsed == Gutter::SeamFill) {
        r.map = bake::bakeSeamFill(atlas.mesh, atlas.hem, atlas.uv, atlas.charts.faceToChart, lowN, high, highN,
                                   opt.resolution, opt.dilate, tangent);
        r.dTotal = r.dDilate = seam::seamEnergy(tm, r.map);
        return r;
    }

    // Dilate and Erased share one bake; erasure is applied on top so the two differ only by the erasure.
    r.map = bake::bake(atlas.mesh, atlas.uv, lowN, high, highN, opt.resolution, opt.dilate, tangent);
    if (r.gutterUsed == Gutter::Erased) {
        try {
            auto e = seam::erase(tm, r.map, seam::NormalSpace::Object);
            r.map = std::move(e.map);
            r.dDilate = e.dBefore;
            r.dTotal = e.dAfter;
            r.unknownTexels = e.unknownTexels;
            return r;
        } catch (const std::exception &) {
            r.gutterUsed = Gutter::Dilate; // erasure failed: keep the dilate map, say so
        }
    }
    r.dTotal = r.dDilate = seam::seamEnergy(tm, r.map);
    return r;
}

} // namespace decimate
