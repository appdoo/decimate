// This file is part of the decimate library; see decimate.h for version/license details.
// Characterization test for core/raytracing (Spec 04, T1–T4).
//
// Two kinds of assertions:
//   1. Bit-exact anchors (T1–T3): hand-derived values independent of the
//      golden fixture — they pin the golden capture itself against arithmetic.
//   2. Golden comparison (T4): for every (function, input) row in the fixed
//      input set the core port's output must equal the committed golden value
//      within epsilon = 1e-6 (booleans exact).

#include <decimate/raytracing/raytracing.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

using V3 = decimate::math::Vec3f;

namespace {

/// Relative tolerance for the golden comparison (booleans are compared exactly).
constexpr float kEps = 1e-6f;


/// One golden fixture line: the hit flag plus whichever of t/u/v/alpha/beta the function reports.
struct GoldenRow {
    int   hit = 0;
    float t = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    float alpha = 0.0f;
    float beta = 0.0f;
};


/// Parses a golden fixture file into a key→row map.  Throws if the file cannot be opened.
/// Each line looks like:
///   <key> hit=<0|1> [t=<float>] [u=<float>] [v=<float>] [alpha=<float>] [beta=<float>]
std::unordered_map<std::string, GoldenRow> loadGolden(const std::string &path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open golden fixture: " + path);
    }
    std::unordered_map<std::string, GoldenRow> rows;
    std::string                                line;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        std::istringstream ss(line);
        std::string        key;
        ss >> key;

        GoldenRow   row {};
        std::string token;
        while (ss >> token) {
            auto eq = token.find('=');
            if (eq == std::string::npos) {
                continue;
            }
            std::string name = token.substr(0, eq);
            std::string val = token.substr(eq + 1);
            if (name == "hit") {
                row.hit = std::stoi(val);
            } else if (name == "t") {
                row.t = std::stof(val);
            } else if (name == "u") {
                row.u = std::stof(val);
            } else if (name == "v") {
                row.v = std::stof(val);
            } else if (name == "alpha") {
                row.alpha = std::stof(val);
            } else if (name == "beta") {
                row.beta = std::stof(val);
            }
        }
        rows[key] = row;
    }
    return rows;
}


/// Resolves the golden fixture path relative to the source tree.
/// CTest sets CTEST_SOURCE_DIRECTORY; fall back to common relative paths.
std::string goldenPath() {
    // Check environment variable first (set by CTest or CI)
    const char *srcDir = std::getenv("CTEST_SOURCE_DIRECTORY");
    if (srcDir != nullptr) {
        return std::string(srcDir) + "/tests/fixtures/raytracing_golden.txt";
    }

    // Walk up from __FILE__ to find the fixtures directory
    std::filesystem::path here = std::filesystem::path(__FILE__).parent_path();
    for (int i = 0; i < 4; ++i) {
        auto candidate = here / "tests" / "fixtures" / "raytracing_golden.txt";
        if (std::filesystem::exists(candidate)) {
            return candidate.string();
        }
        here = here.parent_path();
    }

    // Last resort: check relative to cwd
    return "tests/fixtures/raytracing_golden.txt";
}

} // namespace


/// Bit-exact anchors (T1, T2, T3)
TEST_CASE("evaluate_t: axis-aligned plane/ray hit is exact", "[raytracing][T1]") {
    // N=(0,0,1), d=0, O=(0,0,1), D=(0,0,-1)
    // t = -(0 + 1) / (-1) = 1.0 exactly
    V3    nrm { 0.0f, 0.0f, 1.0f };
    V3    origin { 0.0f, 0.0f, 1.0f };
    V3    dir { 0.0f, 0.0f, -1.0f };
    float t = 0.0f;

    bool hit = rtrace::evaluate_t(t, nrm, 0.0f, origin, dir);

    REQUIRE(hit == true);
    REQUIRE(t == 1.0f); // bit-exact: rational result representable in float
}


TEST_CASE("intersect_triangle: interior hit exact barycentric coordinates", "[raytracing][T2]") {
    // orig=(0.25,0.25,1), dir=(0,0,-1), T=[v0(0,0,0),v1(1,0,0),v2(0,1,0)]
    // Expected: hit=true, t=1.0, u=0.25, v=0.25  (exact in float)
    V3    orig { 0.25f, 0.25f, 1.0f };
    V3    dir { 0.0f, 0.0f, -1.0f };
    V3    v0 { 0.0f, 0.0f, 0.0f };
    V3    v1 { 1.0f, 0.0f, 0.0f };
    V3    v2 { 0.0f, 1.0f, 0.0f };
    float t = 0.0f;
    float u = 0.0f;
    float v = 0.0f;

    bool hit = rtrace::intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &v);

    REQUIRE(hit == true);
    REQUIRE(t == 1.0f);
    REQUIRE(u == 0.25f);
    REQUIRE(v == 0.25f);
}


TEST_CASE("intersect_triangle: point outside triangle returns false", "[raytracing][T3]") {
    // orig=(2,2,1), dir=(0,0,-1) — point (2,2,0) is outside the triangle
    V3    orig { 2.0f, 2.0f, 1.0f };
    V3    dir { 0.0f, 0.0f, -1.0f };
    V3    v0 { 0.0f, 0.0f, 0.0f };
    V3    v1 { 1.0f, 0.0f, 0.0f };
    V3    v2 { 0.0f, 1.0f, 0.0f };
    float t = 0.0f;
    float u = 0.0f;
    float v = 0.0f;

    bool hit = rtrace::intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &v);

    REQUIRE(hit == false);
}


/// Golden comparison tests (T4)
TEST_CASE("raytracing golden: evaluate_t interior hit", "[raytracing][T4][golden]") {
    auto  golden = loadGolden(goldenPath());
    auto &g = golden.at("evaluate_t");

    V3    nrm { 0.0f, 0.0f, 1.0f };
    V3    origin { 0.0f, 0.0f, 1.0f };
    V3    dir { 0.0f, 0.0f, -1.0f };
    float t = 0.0f;
    bool  hit = rtrace::evaluate_t(t, nrm, 0.0f, origin, dir);

    REQUIRE(static_cast<int>(hit) == g.hit);
    if (g.hit != 0) {
        REQUIRE(t == Catch::Approx(g.t).epsilon(kEps));
    }
}


TEST_CASE("raytracing golden: evaluate_t parallel ray (no hit)", "[raytracing][T4][golden]") {
    auto  golden = loadGolden(goldenPath());
    auto &g = golden.at("evaluate_t_parallel");

    V3    nrm { 0.0f, 0.0f, 1.0f };
    V3    origin { 0.0f, 0.0f, 1.0f };
    V3    dir { 1.0f, 0.0f, 0.0f };
    float t = 0.0f;
    bool  hit = rtrace::evaluate_t(t, nrm, 0.0f, origin, dir);

    REQUIRE(static_cast<int>(hit) == g.hit);
}


TEST_CASE("raytracing golden: intersect_triangle interior hit", "[raytracing][T4][golden]") {
    auto  golden = loadGolden(goldenPath());
    auto &g = golden.at("intersect_triangle_interior");

    V3    orig { 0.25f, 0.25f, 1.0f };
    V3    dir { 0.0f, 0.0f, -1.0f };
    V3    v0 { 0.0f, 0.0f, 0.0f };
    V3    v1 { 1.0f, 0.0f, 0.0f };
    V3    v2 { 0.0f, 1.0f, 0.0f };
    float t = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    bool  hit = rtrace::intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &v);

    REQUIRE(static_cast<int>(hit) == g.hit);
    REQUIRE(t == Catch::Approx(g.t).epsilon(kEps));
    REQUIRE(u == Catch::Approx(g.u).epsilon(kEps));
    REQUIRE(v == Catch::Approx(g.v).epsilon(kEps));
}


TEST_CASE("raytracing golden: intersect_triangle outside miss", "[raytracing][T4][golden]") {
    auto  golden = loadGolden(goldenPath());
    auto &g = golden.at("intersect_triangle_outside");

    V3    orig { 2.0f, 2.0f, 1.0f };
    V3    dir { 0.0f, 0.0f, -1.0f };
    V3    v0 { 0.0f, 0.0f, 0.0f };
    V3    v1 { 1.0f, 0.0f, 0.0f };
    V3    v2 { 0.0f, 1.0f, 0.0f };
    float t = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    bool  hit = rtrace::intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &v);

    REQUIRE(static_cast<int>(hit) == g.hit);
}


TEST_CASE("raytracing golden: intersect_triangle parallel ray", "[raytracing][T4][golden]") {
    auto  golden = loadGolden(goldenPath());
    auto &g = golden.at("intersect_triangle_parallel");

    V3    orig { 0.25f, 0.25f, 1.0f };
    V3    dir { 1.0f, 0.0f, 0.0f };
    V3    v0 { 0.0f, 0.0f, 0.0f };
    V3    v1 { 1.0f, 0.0f, 0.0f };
    V3    v2 { 0.0f, 1.0f, 0.0f };
    float t = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    bool  hit = rtrace::intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &v);

    REQUIRE(static_cast<int>(hit) == g.hit);
}


TEST_CASE("raytracing golden: intersect_triangle origin behind triangle",
          "[raytracing][T4][golden]") {
    auto  golden = loadGolden(goldenPath());
    auto &g = golden.at("intersect_triangle_behind");

    V3    orig { 0.25f, 0.25f, -1.0f };
    V3    dir { 0.0f, 0.0f, -1.0f };
    V3    v0 { 0.0f, 0.0f, 0.0f };
    V3    v1 { 1.0f, 0.0f, 0.0f };
    V3    v2 { 0.0f, 1.0f, 0.0f };
    float t = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    bool  hit = rtrace::intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &v);

    REQUIRE(static_cast<int>(hit) == g.hit);
    if (g.hit != 0) {
        REQUIRE(t == Catch::Approx(g.t).epsilon(kEps));
        REQUIRE(u == Catch::Approx(g.u).epsilon(kEps));
        REQUIRE(v == Catch::Approx(g.v).epsilon(kEps));
    }
}


TEST_CASE("raytracing golden: intersect_triangle_bfc backface is culled",
          "[raytracing][T4][golden]") {
    auto  golden = loadGolden(goldenPath());
    auto &g = golden.at("intersect_triangle_bfc_backface");

    V3    orig { 0.25f, 0.25f, -1.0f };
    V3    dir { 0.0f, 0.0f, 1.0f };
    V3    v0 { 0.0f, 0.0f, 0.0f };
    V3    v1 { 1.0f, 0.0f, 0.0f };
    V3    v2 { 0.0f, 1.0f, 0.0f };
    float t = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    bool  hit = rtrace::intersect_triangle_bfc(orig, dir, v0, v1, v2, &t, &u, &v);

    REQUIRE(static_cast<int>(hit) == g.hit);
}


TEST_CASE("raytracing golden: intersect_triangle_bfc front face hit", "[raytracing][T4][golden]") {
    auto  golden = loadGolden(goldenPath());
    auto &g = golden.at("intersect_triangle_bfc_front");

    V3    orig { 0.25f, 0.25f, 1.0f };
    V3    dir { 0.0f, 0.0f, -1.0f };
    V3    v0 { 0.0f, 0.0f, 0.0f };
    V3    v1 { 1.0f, 0.0f, 0.0f };
    V3    v2 { 0.0f, 1.0f, 0.0f };
    float t = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    bool  hit = rtrace::intersect_triangle_bfc(orig, dir, v0, v1, v2, &t, &u, &v);

    REQUIRE(static_cast<int>(hit) == g.hit);
    REQUIRE(t == Catch::Approx(g.t).epsilon(kEps));
    REQUIRE(u == Catch::Approx(g.u).epsilon(kEps));
    REQUIRE(v == Catch::Approx(g.v).epsilon(kEps));
}


TEST_CASE("raytracing golden: IntersectTriangle interior hit", "[raytracing][T4][golden]") {
    auto  golden = loadGolden(goldenPath());
    auto &g = golden.at("IntersectTriangle_interior");

    V3    orig { 0.25f, 0.25f, 1.0f };
    V3    dir { 0.0f, 0.0f, -1.0f };
    V3    v0 { 0.0f, 0.0f, 0.0f };
    V3    v1 { 1.0f, 0.0f, 0.0f };
    V3    v2 { 0.0f, 1.0f, 0.0f };
    float t = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    bool  hit = rtrace::IntersectTriangle(orig, dir, v0, v1, v2, &t, &u, &v);

    REQUIRE(static_cast<int>(hit) == g.hit);
    REQUIRE(t == Catch::Approx(g.t).epsilon(kEps));
    REQUIRE(u == Catch::Approx(g.u).epsilon(kEps));
    REQUIRE(v == Catch::Approx(g.v).epsilon(kEps));
}


TEST_CASE("raytracing golden: barycentricBadouel interior", "[raytracing][T4][golden]") {
    auto  golden = loadGolden(goldenPath());
    auto &g = golden.at("barycentricBadouelf_interior");

    float p[3] = { 0.25f, 0.25f, 0.0f };
    float nrm[3] = { 0.0f, 0.0f, 1.0f };
    float pp0[3] = { 0.0f, 0.0f, 0.0f };
    float pp1[3] = { 1.0f, 0.0f, 0.0f };
    float pp2[3] = { 0.0f, 1.0f, 0.0f };
    float alpha = 0.0f;
    float beta = 0.0f;
    bool  hit = rtrace::barycentricBadouel(alpha, beta, p, nrm, pp0, pp1, pp2);

    REQUIRE(static_cast<int>(hit) == g.hit);
    REQUIRE(alpha == Catch::Approx(g.alpha).epsilon(kEps));
    REQUIRE(beta == Catch::Approx(g.beta).epsilon(kEps));
}


TEST_CASE("raytracing golden: barycentricBadouel outside miss", "[raytracing][T4][golden]") {
    auto  golden = loadGolden(goldenPath());
    auto &g = golden.at("barycentricBadouelf_outside");

    float p[3] = { 2.0f, 0.0f, 0.0f };
    float nrm[3] = { 0.0f, 0.0f, 1.0f };
    float pp0[3] = { 0.0f, 0.0f, 0.0f };
    float pp1[3] = { 1.0f, 0.0f, 0.0f };
    float pp2[3] = { 0.0f, 1.0f, 0.0f };
    float alpha = 0.0f;
    float beta = 0.0f;
    bool  hit = rtrace::barycentricBadouel(alpha, beta, p, nrm, pp0, pp1, pp2);

    REQUIRE(static_cast<int>(hit) == g.hit);
}
