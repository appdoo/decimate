// This file is part of the decimate library; see decimate.h for version/license details.
/// @file bake_bvh_test.cpp
/// @brief Spec 22-01: exact BVH nearest-hit + deterministic parallel bake (I1, I3 cow, I4, B1, S1-S7, bench).
///
/// The oracle is the 2003 brute-force loop (bvhtest::brute, tests/bake_bvh_common.h). "Identical"
/// means same face id and t/u/v bit-equal. Bunny/horse exactness (I2, I3) lives in the `slow` exe
/// bake_bvh_equality_slow.

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cstdio>

#include "bake_bvh_common.h"
#include "model_fixture.h"

using namespace bvhtest;
using decimate::MeshData;
using decimate::bake::NormalMap;
using decimate::bake::Rgb;
namespace detail = decimate::bake::detail;

namespace {

/// The bake's |t| bound (a hit at exactly this distance is a miss).
constexpr float kMax = 100000.0f;

/// A triangle soup from explicit triples (3 vertices each, one face per triple).
MeshData soup(const std::vector<std::array<Vec3, 3>> &tris) {
    MeshData m;
    for (const auto &t : tris) {
        for (const Vec3 &v : t) {
            m.indices.push_back(static_cast<int>(m.positions.size()));
            m.positions.push_back(v);
        }
    }
    return m;
}


/// Axis-aligned-in-XY triangle at height z covering (0.2, 0.2) in a ray's way; `scale` sets its size
/// and `ox` its centroid offset (powers of two keep t exact).
std::array<Vec3, 3> tri(float z, float scale, float ox) {
    return { Vec3 { ox, 0.0f, z }, Vec3 { ox + scale, 0.0f, z }, Vec3 { ox, scale, z } };
}


/// Far-away filler faces so the BVH has more than one leaf.
void addFillers(std::vector<std::array<Vec3, 3>> &t, int n) {
    for (int i = 0; i < n; i++) {
        t.push_back(tri(1000.0f + static_cast<float>(i), 1.0f, 500.0f + 8.0f * static_cast<float>(i)));
    }
}


/// The cow reduced to 1500 faces, charted at cap 64 and parametrized on a 256 grid: the real-pipeline fixture.
Pipeline cowPipeline() {
    return buildPipeline(modelOrSkip("cow.obj"), 1500, 64, 256);
}

} // namespace

/// S1-S6: tie-break, sign, bound, degenerate rays, empty mesh.
TEST_CASE("S1 tie on |t|: the lowest face index wins", "[bake][bvh][S1]") {
    // Face 0 small, face 5 large; both are hit at t = 4 by the ray from (0.25, 0.25, 4) down -z. They differ in
    // centroid, so they can sit in different leaves; faces 1-4 and 6-13 are far fillers.
    for (const bool swapped : { false, true }) {
        std::vector<std::array<Vec3, 3>> t;
        const auto                       small = tri(0.0f, 1.0f, 0.0f);
        const auto                       big = std::array<Vec3, 3> { Vec3 { -16.0f, -16.0f, 0.0f }, Vec3 { 16.0f, -16.0f, 0.0f },
                                                                     Vec3 { 0.0f, 16.0f, 0.0f } };
        t.push_back(swapped ? big : small);
        addFillers(t, 4);
        t.push_back(swapped ? small : big);
        addFillers(t, 8);
        const MeshData m = soup(t);
        const Vec3     o { 0.25f, 0.25f, 4.0f };
        const Vec3     d { 0.0f, 0.0f, -1.0f };

        // Both tracers must pick the lowest face id, whichever leaf holds it.
        const Hit b = brute(m, o, d, kMax);
        REQUIRE(b.face == 0);
        const Hit h = rtrace::TriangleBvh(m).nearest(o, d, kMax);
        REQUIRE(h.face == 0);
        REQUIRE(sameHit(h, b));
        // the precondition of a real tie: face 5 also hits at the identical t
        Vec3  a = m.positions[m.indices[15]], b1 = m.positions[m.indices[16]], c = m.positions[m.indices[17]];
        Vec3  oo = o, dd = d;
        float t5 = 0, u = 0, v = 0;
        REQUIRE(rtrace::intersect_triangle(oo, dd, a, b1, c, &t5, &u, &v));
        REQUIRE(bits(t5) == bits(h.t));
    }
}

TEST_CASE("S2 tie is on |t|, sign-blind", "[bake][bvh][S2]") {
    for (const bool swapped : { false, true }) {
        std::vector<std::array<Vec3, 3>> t;
        t.push_back(tri(swapped ? 2.0f : -2.0f, 4.0f, 0.0f));
        t.push_back(tri(swapped ? -2.0f : 2.0f, 4.0f, 0.0f));
        const MeshData m = soup(t);
        const Vec3     o { 0.25f, 0.25f, 0.0f };
        const Vec3     d { 0.0f, 0.0f, 1.0f };
        const Hit      h = rtrace::TriangleBvh(m).nearest(o, d, kMax);
        REQUIRE(h.face == 0);
        REQUIRE(sameHit(h, brute(m, o, d, kMax)));
        REQUIRE(std::fabs(h.t) == 2.0F);
    }
}

TEST_CASE("S3 a negative t counts (line test)", "[bake][bvh][S3]") {
    std::vector<std::array<Vec3, 3>> t;
    t.push_back(tri(2.0f, 4.0f, 0.0f));  // face 0 at t = +2
    t.push_back(tri(-1.0f, 4.0f, 0.0f)); // face 1 at t = -1
    const MeshData m = soup(t);
    const Hit      h = rtrace::TriangleBvh(m).nearest({ 0.25f, 0.25f, 0.0f }, { 0.0f, 0.0f, 1.0f }, kMax);
    REQUIRE(h.face == 1);
    REQUIRE(h.t == -1.0F);
}

TEST_CASE("S4 |t| = 100000 misses, 99999 hits", "[bake][bvh][S4]") {
    for (const auto &[z, face] : { std::pair { 100000.0f, -1 }, std::pair { 99999.0f, 0 } }) {
        std::vector<std::array<Vec3, 3>> t;
        t.push_back(tri(z, 4.0f, 0.0f));
        const MeshData m = soup(t);
        const Vec3     o { 0.25f, 0.25f, 0.0f };
        const Vec3     d { 0.0f, 0.0f, 1.0f };
        const Hit      h = rtrace::TriangleBvh(m).nearest(o, d, kMax);
        REQUIRE(h.face == face);
        REQUIRE(sameHit(h, brute(m, o, d, kMax)));
    }
}

TEST_CASE("S5 degenerate rays miss like the brute force", "[bake][bvh][S5]") {
    std::vector<std::array<Vec3, 3>> t;
    t.push_back(tri(0.0f, 4.0f, 0.0f));
    addFillers(t, 8);
    const MeshData            m = soup(t);
    const rtrace::TriangleBvh bvh(m);
    const float               nan = std::numeric_limits<float>::quiet_NaN();
    const float               inf = std::numeric_limits<float>::infinity();
    const std::vector<Ray>    rays = {
        { { 0.25f, 0.25f, 1.0f }, { 0.0f, 0.0f, 0.0f } },
        { { nan, 0.25f, 1.0f }, { 0.0f, 0.0f, -1.0f } },
        { { 0.25f, 0.25f, 1.0f }, { nan, 0.0f, -1.0f } },
        { { 0.25f, nan, nan }, { 0.0f, nan, 0.0f } },
        { { inf, 0.25f, 1.0f }, { 0.0f, 0.0f, -1.0f } },
        { { 0.25f, 0.25f, 1.0f }, { inf, 0.0f, -1.0f } },
    };
    for (const Ray &r : rays) {
        const Hit h = bvh.nearest(r.o, r.d, kMax);
        REQUIRE(sameHit(h, brute(m, r.o, r.d, kMax)));
    }
    REQUIRE(bvh.nearest(rays[0].o, rays[0].d, kMax).face == -1);
    REQUIRE(bvh.nearest(rays[1].o, rays[1].d, kMax).face == -1);
    REQUIRE(bvh.nearest(rays[2].o, rays[2].d, kMax).face == -1);
}

TEST_CASE("S6 empty high mesh: every query misses, the bake falls back to the low normal", "[bake][bvh][S6]") {
    const MeshData            empty;
    const rtrace::TriangleBvh bvh(empty);
    REQUIRE(bvh.nearest({ 0, 0, 0 }, { 0, 0, 1 }, kMax).face == -1);

    // One low triangle facing +Z covering the whole atlas.
    MeshData low;
    low.positions = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 } };
    low.indices = { 0, 1, 2 };
    const std::vector<decimate::math::Vec2f> uv = { { 0.0f, 0.0f }, { 1.0f, 0.0f }, { 0.0f, 1.0f } };
    const std::vector<Vec3>                  n(3, Vec3 { 0, 0, 1 });
    const NormalMap                          nm = decimate::bake::bake(low, uv, n, empty, {}, 16, 0, false);
    REQUIRE(nm.pixels[2 * 16 + 2] == Rgb { 127, 127, 254 });
}

/// S7: last-writer-wins ownership (serial order), any thread count.
TEST_CASE("S7 overlapping UV triangles: face 1 owns the shared texels", "[bake][bvh][S7]") {
    MeshData low;
    low.positions = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 } };
    low.indices = { 0, 1, 2, 3, 4, 5 };
    const std::vector<decimate::math::Vec2f> uv = { { 0.1f, 0.1f }, { 0.9f, 0.1f }, { 0.1f, 0.9f }, { 0.1f, 0.1f }, { 0.9f, 0.1f }, { 0.1f, 0.9f } };
    const std::vector<Vec3>                  n = { { 0, 0, 1 }, { 0, 0, 1 }, { 0, 0, 1 }, { 1, 0, 0 }, { 1, 0, 0 }, { 1, 0, 0 } };
    const MeshData                           none;
    for (const int threads : { 1, 0 }) {
        const NormalMap nm = decimate::bake::bake(low, uv, n, none, {}, 32, 0, false, threads);
        int             covered = 0;
        for (const Rgb &p : nm.pixels) {
            if (p == Rgb { 128, 128, 128 }) {
                continue;
            }
            covered++;
            REQUIRE(p == Rgb { 254, 127, 127 }); // face 1's +X normal, never face 0's (127,127,254)
        }
        REQUIRE(covered > 100);
    }
}

/// I1 / I3 (cow) / I4 / B1: the real pipeline on the cow model.
TEST_CASE("I3 cow: BVH nearest == brute nearest on 10000 random rays", "[bake][bvh][I3]") {
    const MeshData            cow = decimate::load(modelOrSkip("cow.obj"));
    const rtrace::TriangleBvh bvh(cow);
    const std::vector<Ray>    rays = randomRays(cow, 10000);
    int                       hits = 0;
    int                       misses = 0;
    int                       mismatches = 0;
    for (const Ray &r : rays) {
        const Hit b = brute(cow, r.o, r.d, kMax);
        const Hit h = bvh.nearest(r.o, r.d, kMax);
        mismatches += sameHit(h, b) ? 0 : 1;
        (b.face >= 0 ? hits : misses)++;
    }
    REQUIRE(mismatches == 0);
    REQUIRE(hits >= 1);
    REQUIRE(misses >= 1);
}

TEST_CASE("I3 cow: BVH == brute on the real bake rays (texel origins + normals)", "[bake][bvh][I3]") {
    const Pipeline            p = cowPipeline();
    const rtrace::TriangleBvh bvh(p.high);
    int                       mismatches = 0;
    int                       count = 0;
    // Every 7th low-vertex ray along its vertex normal, both directions of the line.
    for (std::size_t v = 0; v < p.low.positions.size(); v += 7) {
        for (const float s : { 1.0f, -1.0f }) {
            const Vec3 d = p.lowN[v] * s;
            mismatches += sameHit(bvh.nearest(p.low.positions[v], d, kMax), brute(p.high, p.low.positions[v], d, kMax)) ? 0 : 1;
            count++;
        }
    }
    REQUIRE(count > 100);
    REQUIRE(mismatches == 0);
}

TEST_CASE("I1 cow: BVH bake == brute-force bake, both tangent modes", "[bake][bvh][I1]") {
    const Pipeline    p = cowPipeline();
    const BruteTracer oracle(p.high);
    for (const bool tangent : { false, true }) {
        const NormalMap want = detail::bakeWith(oracle, p.low, p.cornerUV, p.lowN, p.high, p.highN, 256, 0, tangent, 1);
        const NormalMap got = decimate::bake::bake(p.low, p.cornerUV, p.lowN, p.high, p.highN, 256, 0, tangent);
        REQUIRE(want.pixels.size() == 65536);
        REQUIRE(samePixels(got, want));
    }
}

TEST_CASE("I4 cow: bakeSeamFill is thread-count independent", "[bake][bvh][I4]") {
    const Pipeline  p = cowPipeline();
    const NormalMap a = decimate::bake::bakeSeamFill(p.low, p.hem, p.cornerUV, p.faceToChart, p.lowN, p.high, p.highN,
                                                     256, 2, false, 1);
    const NormalMap b = decimate::bake::bakeSeamFill(p.low, p.hem, p.cornerUV, p.faceToChart, p.lowN, p.high, p.highN,
                                                     256, 2, false, 0);
    REQUIRE(samePixels(a, b));
}

TEST_CASE("B1 cow: bake is byte-identical for threads 1, 2, 3, 7, 0", "[bake][bvh][B1]") {
    const Pipeline  p = cowPipeline();
    const NormalMap serial = decimate::bake::bake(p.low, p.cornerUV, p.lowN, p.high, p.highN, 256, 0, false, 1);
    for (const int threads : { 2, 3, 7, 0 }) {
        const NormalMap par = decimate::bake::bake(p.low, p.cornerUV, p.lowN, p.high, p.highN, 256, 0, false, threads);
        REQUIRE(samePixels(par, serial));
    }
}

/// Bench (hidden): bake seconds for brute / BVH 1 thread / BVH all threads.
TEST_CASE("bake timing: brute vs BVH vs BVH+threads", "[.bench]") {
    using Clock = std::chrono::steady_clock;
    const auto sec = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };
    for (const char *name : { "cow.obj", "bunny.iv", "horse.iv" }) {
        const Pipeline            p = buildPipeline(modelOrSkip(name), 1500, 64, 256);
        const auto                t0 = Clock::now();
        const rtrace::TriangleBvh bvh(p.high);
        const auto                t1 = Clock::now();
        const BruteTracer         oracle(p.high);
        const NormalMap           ref = detail::bakeWith(oracle, p.low, p.cornerUV, p.lowN, p.high, p.highN, 256, 0, false, 1);
        const auto                t2 = Clock::now();
        const NormalMap           one = decimate::bake::bake(p.low, p.cornerUV, p.lowN, p.high, p.highN, 256, 0, false, 1);
        const auto                t3 = Clock::now();
        const NormalMap           all = decimate::bake::bake(p.low, p.cornerUV, p.lowN, p.high, p.highN, 256, 0, false, 0);
        const auto                t4 = Clock::now();
        std::printf("[bench] %-9s faces=%zu  bvh-build %.1f ms | brute %.3f s | bvh x1 %.3f s | bvh xN %.3f s | equal=%d\n",
                    name, p.high.indices.size() / 3, 1000.0 * sec(t0, t1), sec(t1, t2), sec(t2, t3), sec(t3, t4),
                    static_cast<int>(samePixels(ref, one) && samePixels(ref, all)));
    }
}
