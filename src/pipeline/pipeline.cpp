// This file is part of the decimate library; see decimate.h for version/license details.
#include <decimate/pipeline.h>

#include <chrono>
#include <exception>
#include <utility>

#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normals.h>
#include <decimate/bake/tangents.h>
#include <decimate/seam/erase.h>
#include <decimate/seam/seam_energy.h>
#include <decimate/seam/seam_topology.h>

namespace decimate {

static double secondsSince(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}


/// Low-poly vertex normals of the bake. Tangents and bake share this one definition: EqualUnit, like the high mesh.
static std::vector<math::Vec3f> lowNormals(const MeshData &mesh) {
    return bake::computeVertexNormals(mesh, bake::Weighting::EqualUnit);
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

    // The tangents read the same EqualUnit normals bakeAtlas bakes with (lowNormals), so frame and bake agree.
    a.tangents = tangents(low, a.uv, lowNormals(low));
    return a;
}


/// Seamless is a decision, not a mode: it is Erased, in both spaces.
static Gutter resolveGutter(Gutter g) {
    return g == Gutter::Seamless ? Gutter::Erased : g;
}


/// Erase the seams of @p map in place; the frames of the atlas make the tangent space well-defined.
static seam::EraseResult eraseSeams(const Atlas &atlas, const TexturedMesh &tm, const std::vector<math::Vec3f> &lowN,
                                    const bake::NormalMap &map, BakeSpace space) {
    if (space == BakeSpace::Object) {
        return seam::erase(tm, map, seam::NormalSpace::Object);
    }
    return seam::erase(tm, map, seam::TangentFrames { atlas.tangents, lowN });
}


BakeResult bakeAtlas(const Atlas &atlas, const MeshData &high, const BakeOptions &opt) {
    // Low and high poly use the same vertex-normal weighting so the bake compares like with like.
    const auto                     lowN = lowNormals(atlas.mesh);
    const auto                     highN = bake::computeVertexNormals(high, bake::Weighting::EqualUnit);
    const std::vector<math::Vec4f> objectSpace; // empty tangents select object space
    const auto                    &tangent = opt.space == BakeSpace::Tangent ? atlas.tangents : objectSpace;

    // D_total (the seam metric) is defined on the welded textured mesh, so every gutter builds it.
    const auto tm = seam::weld(atlas.mesh, atlas.uv);

    BakeResult r;
    r.gutterUsed = resolveGutter(opt.gutter);
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
            auto e = eraseSeams(atlas, tm, lowN, r.map, opt.space);
            r.map = std::move(e.map);
            r.dDilate = e.dBefore;
            r.dTotal = e.dAfter;
            r.unknownTexels = e.unknownTexels;
            r.rescale = e.rescale;
            return r;
        } catch (const std::exception &) {
            r.gutterUsed = Gutter::Dilate; // erasure failed: keep the dilate map, say so
        }
    }
    r.dTotal = r.dDilate = seam::seamEnergy(tm, r.map);
    return r;
}

} // namespace decimate
