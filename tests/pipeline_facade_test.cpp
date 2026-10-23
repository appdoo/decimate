// This file is part of the decimate library; see decimate.h for version/license details.
/// @file pipeline_facade_test.cpp
/// @brief Spec 25: `chartify` + `bakeAtlas` are the viewer's hand-assembled pipeline, moved.
///
/// a [I] cow, viewer settings: facade uv and map bytes == the hand-assembled sequence, 3 gutters x 2 spaces + Seamless
/// b [S] Seamless: Erased in both spaces, object D_total < 1e-8
/// c [S] Erased in tangent space is built frame-aware (Spec 28), reported in gutterUsed
/// d [S] library default = 2003: charts == segment(64), plain flatten, no cuts

#include "model_fixture.h"

#include <catch2/catch_test_macros.hpp>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normals.h>
#include <decimate/pipeline.h>
#include <decimate/seam/erase.h>
#include <decimate/seam/seam_topology.h>

using namespace decimate;

namespace {

/// The viewer's 16 charts, stretch bound off (D3).
ChartOptions viewerOptions() {
    ChartOptions o;
    o.maxCharts = 16;
    return o;
}


const MeshData &cow() {
    static const MeshData m = load(modelOrSkip("cow.obj"));
    return m;
}


const Atlas &viewerAtlas() {
    static const Atlas a = chartify(cow(), viewerOptions());
    return a;
}


/// The viewer's makeGutterMap, copied from before the move (SeamFill / dilate / erase on top of dilate).
bake::NormalMap handBake(const Atlas &a, const MeshData &high, Gutter g, const std::vector<math::Vec4f> &tangent) {
    const auto lowN = bake::computeVertexNormals(a.mesh, bake::Weighting::EqualUnit);
    const auto highN = bake::computeVertexNormals(high, bake::Weighting::EqualUnit);
    if (g == Gutter::SeamFill) {
        return bake::bakeSeamFill(a.mesh, a.hem, a.uv, a.charts.faceToChart, lowN, high, highN, 256, 2, tangent);
    }
    auto map = bake::bake(a.mesh, a.uv, lowN, high, highN, 256, 2, tangent);
    if (g == Gutter::Erased) {
        const auto tm = seam::weld(a.mesh, a.uv);
        if (tangent.empty()) {
            return seam::erase(tm, map, seam::NormalSpace::Object).map;
        }
        return seam::erase(tm, map, seam::TangentFrames { tangent, lowN }).map;
    }
    return map;
}


bool sameMap(const bake::NormalMap &x, const bake::NormalMap &y) {
    return x.width == y.width && x.height == y.height && x.pixels == y.pixels;
}

} // namespace


TEST_CASE("a: facade equals the hand-assembled viewer path on the cow", "[pipeline][integration]") {
    const Atlas &a = viewerAtlas();

    // The hand-assembled chartification (viewer's computeChartification before the move).
    HalfEdgeMesh hem = HalfEdgeMesh::create(static_cast<int>(cow().positions.size()),
                                            static_cast<int>(cow().indices.size() / 3),
                                            reinterpret_cast<const Face *>(cow().indices.data()));
    auto         charts = chart::segment(hem, cow(), 16);
    auto         packed = chart::flattenCharts(hem, cow(), charts, {});
    chart::packAtlas(packed, cow(), 256);
    CHECK(a.uv == chart::cornerUVs(packed, static_cast<int>(cow().indices.size() / 3)));
    CHECK(a.charts.faceToChart == charts.faceToChart);

    struct Row {
        Gutter    requested, built;
        BakeSpace space;
    };
    const Row rows[] = {
        { Gutter::Dilate, Gutter::Dilate, BakeSpace::Object },
        { Gutter::Dilate, Gutter::Dilate, BakeSpace::Tangent },
        { Gutter::SeamFill, Gutter::SeamFill, BakeSpace::Object },
        { Gutter::SeamFill, Gutter::SeamFill, BakeSpace::Tangent },
        { Gutter::Erased, Gutter::Erased, BakeSpace::Object },
        { Gutter::Erased, Gutter::Erased, BakeSpace::Tangent },
        { Gutter::Seamless, Gutter::Erased, BakeSpace::Object },
        { Gutter::Seamless, Gutter::Erased, BakeSpace::Tangent },
    };
    // Spec 26: chartify fills one MikkTSpace tangent per corner (empty would silently bake object space).
    REQUIRE(a.tangents.size() == a.uv.size());
    for (const Row &r : rows) {
        BakeOptions opt;
        opt.space = r.space;
        opt.gutter = r.requested;
        const BakeResult res = bakeAtlas(a, cow(), opt);
        CHECK(res.gutterUsed == r.built);
        CHECK(sameMap(res.map, handBake(a, cow(), r.built, r.space == BakeSpace::Tangent ? a.tangents : std::vector<math::Vec4f> {})));
    }
}


TEST_CASE("b: Seamless picks the gutter that hides the seams per space", "[pipeline][state]") {
    const Atlas &a = viewerAtlas();

    const BakeResult object = bakeAtlas(a, cow(), {});
    CHECK(object.gutterUsed == Gutter::Erased);
    for (const double d : object.dTotal) {
        CHECK(d < 1e-8);
    }
    CHECK(object.dDilate[0] > 1e-3); // the dilate map it replaced has a visible seam

    BakeOptions tangent;
    tangent.space = BakeSpace::Tangent;
    const BakeResult t = bakeAtlas(a, cow(), tangent);
    CHECK(t.gutterUsed == Gutter::Erased);
    CHECK(t.dTotal[0] + t.dTotal[1] + t.dTotal[2] < 1e-8);
}


TEST_CASE("c: Erased in tangent space is frame-aware and says so", "[pipeline][state]") {
    BakeOptions opt;
    opt.space = BakeSpace::Tangent;
    opt.gutter = Gutter::Erased;
    const BakeResult r = bakeAtlas(viewerAtlas(), cow(), opt);
    CHECK(r.gutterUsed == Gutter::Erased);
    CHECK(r.unknownTexels > 0);
    CHECK(r.rescale >= 0.7);
}


TEST_CASE("d: the library default is the 2003 pipeline", "[pipeline][state]") {
    const Atlas a = chartify(cow());

    HalfEdgeMesh hem = HalfEdgeMesh::create(static_cast<int>(cow().positions.size()),
                                            static_cast<int>(cow().indices.size() / 3),
                                            reinterpret_cast<const Face *>(cow().indices.data()));
    const auto   plain = chart::segment(hem, cow(), 64);
    CHECK(a.charts.charts.size() == plain.charts.size());
    CHECK(a.charts.faceToChart == plain.faceToChart);
    auto charts = plain;
    auto packed = chart::flattenCharts(hem, cow(), charts, {});
    chart::packAtlas(packed, cow(), 256);
    CHECK(a.uv == chart::cornerUVs(packed, static_cast<int>(cow().indices.size() / 3)));
    CHECK(a.stats.cuts == 0);
    CHECK(a.stats.closedOpened == 0);
    CHECK(a.stats.splits == 0);
}
