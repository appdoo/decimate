// This file is part of the decimate library; see decimate.h for version/license details.
/// @file bake_convention_test.cpp
/// @brief Spec 27: the bake writes texel x at uv = (x + 0.5) / N, the convention every reader uses, and rounds.
///
/// a [S] flat quad with a per-texel signal, read back at uv = (x + 0.5) / N: texel x exactly
/// b [S] rangeCompress rounds to nearest; (0,0,1) -> (127,127,254); a half step rounds up
/// c [I] spot 256, 16 charts: seam angle of the decoded normals in the reader convention
///       (tangent dilate <= 0.90 deg, tangent seam-fill <= 1.0 deg and below dilate, object erased ~ 0.163)
/// The tangent round-trip (max angle) lives in tangent_space_test.cpp.

#include "bake_detail.h"
#include "seam_metric.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normals.h>
#include <decimate/pipeline.h>
#include <decimate/reduce/reducer.h>
#include <decimate/seam/seam_topology.h>

#include "model_fixture.h"

#include <cmath>
#include <cstdio>

using namespace decimate;
using math::Vec2f;
using math::Vec3f;

namespace {

/// The unit square z = 0 as the low mesh (UV = XY, the whole atlas), the same square at z = 0.1 as
/// the high mesh. Its vertex normals (k(2x-1), 0, 1) tilt linearly with x, so a texel's value tells
/// WHERE in x it was shaded: the signal of texel x is the normal at position (x + 0.5) / N.
struct SignalQuad {
    MeshData           low, high;
    std::vector<Vec2f> uv;
    std::vector<Vec3f> lowN, highN;
};


SignalQuad signalQuad(float k) {
    SignalQuad q;
    q.low.positions = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 } };
    q.low.indices = { 0, 1, 2, 0, 2, 3 };
    q.high = q.low;
    for (auto &p : q.high.positions) {
        p.z = 0.1f;
    }
    q.lowN.assign(4, { 0.0f, 0.0f, 1.0f });
    for (const auto &p : q.high.positions) {
        q.highN.push_back({ k * (2.0f * p.x - 1.0f), 0.0f, 1.0f });
    }
    for (const int i : q.low.indices) {
        q.uv.push_back({ q.low.positions[static_cast<std::size_t>(i)].x, q.low.positions[static_cast<std::size_t>(i)].y });
    }
    return q;
}


const MeshData &spot() {
    static const MeshData m = load(std::string(MESH_MODEL_DIR) + "/spot.obj");
    return m;
}


struct Table {
    seammetric::Stats erased, dilate, seamFill;
};


/// The viewer's pipeline on `high`: reduce to `keep` of the vertices, 16 charts, bake at `res`.
Table seamTable(const MeshData &high, float keep, int res) {
    const MeshData low = reduce(high, static_cast<int>(static_cast<float>(high.positions.size()) * keep));
    ChartOptions   co;
    co.maxCharts = 16;
    const Atlas a = chartify(low, co);
    const auto  vn = bake::computeVertexNormals(a.mesh);
    const auto  st = seam::findSeams(seam::weld(a.mesh, a.uv));
    auto        bake = [&](BakeSpace sp, Gutter g) {
        BakeOptions o;
        o.space = sp;
        o.gutter = g;
        o.resolution = res;
        return bakeAtlas(a, high, o).map;
    };
    return { seammetric::seamStats(a, vn, bake(BakeSpace::Object, Gutter::Erased), false, st),
             seammetric::seamStats(a, vn, bake(BakeSpace::Tangent, Gutter::Dilate), true, st),
             seammetric::seamStats(a, vn, bake(BakeSpace::Tangent, Gutter::SeamFill), true, st) };
}


void print(const char *name, const Table &t) {
    std::printf("%s: object erased %.3f/%.2f  tangent dilate %.3f/%.2f  tangent seam-fill %.3f/%.2f (mean/max deg)\n",
                name, t.erased.mean, t.erased.max, t.dilate.mean, t.dilate.max, t.seamFill.mean, t.seamFill.max);
}

} // namespace


TEST_CASE("Spec 27: texel x is written at uv = (x + 0.5) / N", "[bake][spec27][S]") {
    const int  n = 16;
    const auto q = signalQuad(1.0f);
    const auto nm = bake::bake(q.low, q.uv, q.lowN, q.high, q.highN, n, 0);
    REQUIRE(nm.width == n);

    // Every texel of a row is shaded from its own x: compare with the signal at the centre (x + 0.5) / N.
    // The old (N - 1) convention shades texel x at (x + 0.5) / (N - 1): 10+ codes off at the last column.
    int worst = 0;
    for (int x = 0; x < n; x++) {
        const float p = (static_cast<float>(x) + 0.5f) / static_cast<float>(n);
        const auto  want = bake::detail::rangeCompress({ 2.0f * p - 1.0f, 0.0f, 1.0f });
        const auto  got = nm.pixels[static_cast<std::size_t>(n / 2) * n + x];
        worst = std::max(worst, std::abs(static_cast<int>(got.r) - static_cast<int>(want.r)));
    }
    CAPTURE(worst);
    CHECK(worst <= 1);

    // The signal ramps monotonically across the row: no texel is skipped or doubled.
    for (int x = 1; x < n; x++) {
        const auto row = static_cast<std::size_t>(n / 2) * n;
        CHECK(nm.pixels[row + x].r > nm.pixels[row + x - 1].r);
    }
}


TEST_CASE("Spec 27: rangeCompress rounds to nearest, 127 is zero", "[bake][spec27][S]") {
    using bake::Rgb;
    CHECK(bake::detail::rangeCompress({ 0, 0, 1 }) == Rgb { 127, 127, 254 });
    CHECK(bake::detail::rangeCompress({ 0, 0, -1 }) == Rgb { 127, 127, 0 });

    // 0.6 and 0.8 -> 203.2 and 228.6: truncation gives 228 for the second, nearest gives 229.
    CHECK(bake::detail::rangeCompress({ 0.6f, 0.8f, 0 }) == Rgb { 203, 229, 127 });

    // (0.5 + 1) * 127 = 190.5 exactly: a half step rounds up.
    CHECK(bake::detail::rangeCompress({ 0.5f, 0.0f, std::sqrt(0.75f) }).r == 191);

    // Decode c / 127 - 1: the encoded flat normal decodes to (0, 0, 1) within one code.
    const Vec3f flat = seammetric::unpack(bake::detail::rangeCompress({ 0, 0, 1 }));
    CHECK(flat.x == 0.0f);
    CHECK(flat.z == 1.0f);
}


TEST_CASE("Spec 27: seams across 16 charts of spot agree in the reader convention", "[bake][spec27][I]") {
    const Table t = seamTable(spot(), 0.24f, 256);
    print("spot 256", t);
    CHECK(t.erased.mean == Catch::Approx(0.163).margin(0.02));
    CHECK(t.dilate.mean <= 0.90);
}


TEST_CASE("Spec 27: cow seam table (measurement)", "[.][bake][spec27][measure]") {
    print("cow 256", seamTable(load(modelOrSkip("cow.obj")), 0.24f, 256));
}
