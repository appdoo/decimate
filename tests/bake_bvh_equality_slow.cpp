// This file is part of the decimate library; see decimate.h for version/license details.
/// @file bake_bvh_equality_slow.cpp
/// @brief Spec 22-01 I2 + I3 on the dense models (bunny, horse): the brute-force oracle is minutes
///        long here, hence the own executable and the `slow` ctest label.

#include <catch2/catch_test_macros.hpp>

#include "bake_bvh_common.h"
#include "model_fixture.h"

using namespace bvhtest;
namespace detail = decimate::bake::detail;

TEST_CASE("I3 dense: BVH nearest == brute nearest on 10000 random rays", "[bake][bvh][I3][slow]") {
    for (const char *name : { "bunny.iv", "horse.iv" }) {
        const decimate::MeshData  m = decimate::load(modelOrSkip(name));
        const rtrace::TriangleBvh bvh(m);
        int                       hits = 0;
        int                       misses = 0;
        int                       mismatches = 0;
        for (const Ray &r : randomRays(m, 10000)) {
            const Hit b = brute(m, r.o, r.d, 100000.0f);
            mismatches += sameHit(bvh.nearest(r.o, r.d, 100000.0f), b) ? 0 : 1;
            (b.face >= 0 ? hits : misses)++;
        }
        INFO(name);
        REQUIRE(mismatches == 0);
        REQUIRE(hits >= 1);
        REQUIRE(misses >= 1);
    }
}


TEST_CASE("I2 dense: BVH bake == brute-force bake (bunny, horse)", "[bake][bvh][I2][slow]") {
    for (const char *name : { "bunny.iv", "horse.iv" }) {
        const Pipeline    p = buildPipeline(modelOrSkip(name), 1500, 64, 256);
        const BruteTracer oracle(p.high);
        const auto        want = detail::bakeWith(oracle, p.low, p.cornerUV, p.lowN, p.high, p.highN, 256, 0, false, 0);
        const auto        got = decimate::bake::bake(p.low, p.cornerUV, p.lowN, p.high, p.highN, 256, 0, false);
        INFO(name);
        REQUIRE(samePixels(got, want));
    }
}
