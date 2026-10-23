// This file is part of the decimate library; see decimate.h for version/license details.
/// @file tangent_erase_test.cpp
/// @brief Spec 28: seam erasure with the tangent frames in the mismatch -- tangent-space seams reach the object-space floor.
///
/// a [I] spot 256, 16 charts: seam angle of the decoded normals (mean <= 0.20, max <= 0.6 deg)
/// b [I] cow 256: mean <= 0.25, max <= 1.0 deg
/// c [I] spot 1024 (hidden): mean <= 0.20 deg
/// d [S] spot 256: framed int|r|^2 <= 1e-8, no texel outside [-1, 1], rescale >= 0.7 reported
/// e [S] object space is byte-identical to before the frames existed (pinned hash)
/// f [S] Seamless / Erased in tangent space -> gutterUsed == Erased

#include "model_fixture.h"
#include "seam_metric.h"

#include <catch2/catch_test_macros.hpp>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normals.h>
#include <decimate/pipeline.h>
#include <decimate/reduce/reducer.h>
#include <decimate/seam/erase.h>
#include <decimate/seam/seam_topology.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

using namespace decimate;

namespace {

/// The viewer's pipeline on `high`: reduce to `keep` of the vertices, 16 charts.
struct Fixture {
    MeshData                 high;
    Atlas                    atlas;
    std::vector<math::Vec3f> vn;
    seam::SeamTopology       topo;
};


Fixture makeFixture(const std::string &path, float keep) {
    Fixture f;
    f.high = load(path);
    const MeshData low = reduce(f.high, static_cast<int>(static_cast<float>(f.high.positions.size()) * keep));
    ChartOptions   co;
    co.maxCharts = 16;
    f.atlas = chartify(low, co);
    f.vn = bake::computeVertexNormals(f.atlas.mesh);
    f.topo = seam::findSeams(seam::weld(f.atlas.mesh, f.atlas.uv));
    return f;
}


const Fixture &spot() {
    static const Fixture f = makeFixture(std::string(MESH_MODEL_DIR) + "/spot.obj", 0.24f);
    return f;
}


const Fixture &cow() {
    static const Fixture f = makeFixture(modelOrSkip("cow.obj"), 0.24f);
    return f;
}


BakeResult bakeWith(const Fixture &f, BakeSpace space, Gutter gutter, int res) {
    BakeOptions o;
    o.space = space;
    o.gutter = gutter;
    o.resolution = res;
    return bakeAtlas(f.atlas, f.high, o);
}


seammetric::Stats seamAngle(const Fixture &f, BakeSpace space, Gutter gutter, int res) {
    const bool tangent = space == BakeSpace::Tangent;
    return seammetric::seamStats(f.atlas, f.vn, bakeWith(f, space, gutter, res).map, tangent, f.topo);
}


void printTable(const char *name, const Fixture &f, int res) {
    const auto e = seamAngle(f, BakeSpace::Tangent, Gutter::Erased, res);
    const auto s = seamAngle(f, BakeSpace::Tangent, Gutter::SeamFill, res);
    const auto d = seamAngle(f, BakeSpace::Tangent, Gutter::Dilate, res);
    const auto o = seamAngle(f, BakeSpace::Object, Gutter::Erased, res);
    std::printf("%s %d: tangent erased %.3f/%.2f  seam-fill %.3f/%.2f  dilate %.3f/%.2f  object erased %.3f/%.2f\n", name,
                res, e.mean, e.max, s.mean, s.max, d.mean, d.max, o.mean, o.max);
}


unsigned long long fnv(const bake::NormalMap &m) {
    unsigned long long h = 1469598103934665603ull;
    for (const bake::Rgb &p : m.pixels) {
        for (const int b : { p.r, p.g, p.b }) {
            h = (h ^ static_cast<unsigned>(b)) * 1099511628211ull;
        }
    }
    return h;
}

} // namespace


TEST_CASE("Spec 28: spot 256 tangent seams reach the object floor", "[seam][spec28][I]") {
    printTable("spot", spot(), 256);
    const auto s = seamAngle(spot(), BakeSpace::Tangent, Gutter::Erased, 256);
    CHECK(s.mean <= 0.20);
    CHECK(s.max <= 0.6);
}


TEST_CASE("Spec 28: cow 256 tangent seams reach the object floor", "[seam][spec28][I]") {
    printTable("cow", cow(), 256);
    const auto s = seamAngle(cow(), BakeSpace::Tangent, Gutter::Erased, 256);
    CHECK(s.mean <= 0.25);
    CHECK(s.max <= 1.0);
}


TEST_CASE("Spec 28: spot 1024 tangent seams", "[.][seam][spec28][I]") {
    printTable("spot", spot(), 1024);
    CHECK(seamAngle(spot(), BakeSpace::Tangent, Gutter::Erased, 1024).mean <= 0.20);
}


TEST_CASE("Spec 28: the framed residual is gone and the solution stays in range", "[seam][spec28][S]") {
    const Fixture &f = spot();
    const auto     tm = seam::weld(f.atlas.mesh, f.atlas.uv);
    const auto     dilated = bakeWith(f, BakeSpace::Tangent, Gutter::Dilate, 256).map;
    const auto     r = seam::erase(tm, dilated, seam::TangentFrames { f.atlas.tangents, f.vn });

    CHECK(r.dBefore[0] + r.dBefore[1] + r.dBefore[2] > 1e-4); // the dilate map has a framed seam to erase
    CHECK(r.dAfter[0] + r.dAfter[1] + r.dAfter[2] <= 1e-8);
    CHECK(r.rescale >= 0.7);
    CHECK(r.rescale <= 1.0);

    // Not one solved texel was clipped by the quantizer: the solution is inside [-1, 1] before rounding.
    const auto clipped = std::count_if(r.texels.begin(), r.texels.end(), [](double v) { return std::fabs(v) > 1.0; });
    CHECK(clipped == 0);

    // The facade reports the same rescale.
    CHECK(bakeWith(f, BakeSpace::Tangent, Gutter::Erased, 256).rescale == r.rescale);
}


TEST_CASE("Spec 28: object space is byte-identical to the scalar erasure", "[seam][spec28][S][release-pinned]") {
    // FNV-1a of the spot 256 object-erased map, taken before the frame code existed (build-rel).
    const auto r = bakeWith(spot(), BakeSpace::Object, Gutter::Seamless, 256);
    CHECK(r.rescale == 1.0);
    CHECK(fnv(r.map) == 0x30928efe7fe08abbull);
}


TEST_CASE("Spec 28: Seamless and Erased in tangent space are Erased", "[seam][spec28][S]") {
    for (const Gutter g : { Gutter::Seamless, Gutter::Erased }) {
        const auto r = bakeWith(spot(), BakeSpace::Tangent, g, 128);
        CHECK(r.gutterUsed == Gutter::Erased);
        CHECK(r.unknownTexels > 0);
    }
}


/// Wall time of the framed vs the object erase of one fixture at one resolution (receipt, never asserted).
static void printTiming(const char *name, const Fixture &f, int res) {
    const auto dilate = bakeWith(f, BakeSpace::Tangent, Gutter::Dilate, res).map;
    const auto objMap = bakeWith(f, BakeSpace::Object, Gutter::Dilate, res).map;
    const auto tm = seam::weld(f.atlas.mesh, f.atlas.uv);
    auto       secs = [](auto &&fn) {
        const auto t = std::chrono::steady_clock::now();
        fn();
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
    };
    int          unknownTan = 0;
    double       rescale = 0.0;
    const double tTan = secs([&] { const auto r = seam::erase(tm, dilate, seam::TangentFrames { f.atlas.tangents, f.vn });
        unknownTan = r.unknownTexels;
        rescale = r.rescale; });
    const double tObj = secs([&] { (void)seam::erase(tm, objMap, seam::NormalSpace::Object); });
    std::printf("%s %d: tangent erase %.2f s (%d unknowns, rescale %.4f), object erase %.2f s\n", name, res, tTan, unknownTan, rescale, tObj);
}


TEST_CASE("Spec 28: timing of tangent vs object erasure", "[.][seam][spec28][timing][measure]") {
    printTiming("spot", spot(), 256);
    printTiming("cow", cow(), 256);
    printTiming("spot", spot(), 1024);
}
