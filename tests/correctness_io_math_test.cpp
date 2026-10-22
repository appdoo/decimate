// This file is part of the decimate library; see decimate.h for version/license details.
/// @file correctness_io_math_test.cpp
/// @brief Spec 24-03 (hidden state and math): M1 OBJ loader thread-safety and number
///        parsing, M2 lexer number parsing, M3 sparse-matrix singular pivots, M4 quaternion
///        matrix round-trip.  Exact values throughout.

#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <clocale>
#include <cmath>
#include <decimate/io/lexer.h>
#include <decimate/math/quat.h>
#include <decimate/math/sparse_matrix.h>
#include <decimate/types/mesh_data.h>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

#ifndef MESH_MODEL_DIR
#error "MESH_MODEL_DIR must be defined by CMake"
#endif
#ifndef MESH_FIXTURES_DIR
#error "MESH_FIXTURES_DIR must be defined by CMake"
#endif

namespace {

/// Writes `text` to a per-process temp file named after `tag` and returns its path.
std::string writeTemp(const std::string &tag, const std::string &text) {
    const auto path = std::filesystem::temp_directory_path() /
                      ("decimate_24_03_" + std::to_string(::getpid()) + "_" + tag);
    std::ofstream(path, std::ios::binary) << text;
    return path.string();
}

/// True when both meshes hold identical indices and bit-equal positions.
bool sameMesh(const decimate::MeshData &a, const decimate::MeshData &b) {
    return a.indices == b.indices && a.positions.size() == b.positions.size() &&
           std::equal(a.positions.begin(), a.positions.end(), b.positions.begin(),
                      [](const auto &p, const auto &q) {
                          return p[0] == q[0] && p[1] == q[1] && p[2] == q[2];
                      });
}

/// Loads `path` `rounds` times and returns how many loads differ from `expected`.
int countMismatches(const std::string &path, const decimate::MeshData &expected, int rounds) {
    int bad = 0;
    for (int i = 0; i < rounds; ++i) {
        try {
            bad += sameMesh(decimate::loadObj(path), expected) ? 0 : 1;
        } catch (const std::exception &) {
            ++bad;
        }
    }
    return bad;
}

} // namespace

// ===========================================================================
// M1 — OBJ loader: no hidden state, no locale
// ===========================================================================

TEST_CASE("M1-I: two threads loading cow.obj and spot.obj 50x each equal the single-thread load",
          "[io][m1]") {
    const std::string cow = std::string(MESH_MODEL_DIR) + "/cow.obj";
    const std::string spot = std::string(MESH_MODEL_DIR) + "/spot.obj";
    if (!std::filesystem::exists(cow)) {
        SKIP("cow.obj not present in MESH_MODEL_DIR");
    }

    const decimate::MeshData expectA = decimate::loadObj(cow);
    const decimate::MeshData expectB = decimate::loadObj(spot);

    int         badA = -1;
    int         badB = -1;
    std::thread ta([&] { badA = countMismatches(cow, expectA, 50); });
    std::thread tb([&] { badB = countMismatches(spot, expectB, 50); });
    ta.join();
    tb.join();

    CHECK(badA == 0);
    CHECK(badB == 0);
}

TEST_CASE("M1-S: OBJ number parsing keeps atof/atoi leniency", "[io][m1]") {
    // '+' sign, leading '.', exponent, trailing garbage, missing coordinate,
    // FROZEN BUG 2 slash reference.
    const std::string        path = writeTemp("lenient.obj", "v +1.5 .5 1e2\n"
                                                             "v 1.0abc 2.5\n"
                                                             "v -3 +4 5.25\n"
                                                             "f 1//4 +2/9/9 3\n");
    const decimate::MeshData m = decimate::loadObj(path);
    std::filesystem::remove(path);

    REQUIRE(m.positions.size() == 3);
    CHECK(m.positions[0][0] == 1.5F);
    CHECK(m.positions[0][1] == 0.5F);
    CHECK(m.positions[0][2] == 100.0F);
    CHECK(m.positions[1][0] == 1.0F);
    CHECK(m.positions[1][1] == 2.5F);
    CHECK(m.positions[1][2] == 0.0F); // missing token reads 0
    CHECK(m.positions[2][0] == -3.0F);
    CHECK(m.positions[2][1] == 4.0F);
    CHECK(m.positions[2][2] == 5.25F);
    CHECK(m.indices == std::vector<int> { 0, 1, 2 });
}

TEST_CASE("M1-S: OBJ parsing ignores the process locale", "[io][m1]") {
    const std::string old = std::setlocale(LC_NUMERIC, nullptr);
    if (std::setlocale(LC_NUMERIC, "de_DE.UTF-8") == nullptr) {
        WARN("de_DE.UTF-8 locale unavailable: locale independence NOT exercised");
        return;
    }

    const std::string  path = writeTemp("locale.obj", "v 1.5 2.5 3.5\nv 0 0 0\nv 1 1 1\nf 1 2 3\n");
    decimate::MeshData m;
    try {
        m = decimate::loadObj(path);
    } catch (...) {
        std::setlocale(LC_NUMERIC, old.c_str());
        std::filesystem::remove(path);
        throw;
    }
    std::setlocale(LC_NUMERIC, old.c_str());
    std::filesystem::remove(path);

    REQUIRE(m.positions.size() == 3);
    CHECK(m.positions[0][0] == 1.5F);
    CHECK(m.positions[0][1] == 2.5F);
    CHECK(m.positions[0][2] == 3.5F);
}

// ===========================================================================
// M2 — Inventor lexer: from_chars, lone '-' is an error
// ===========================================================================

TEST_CASE("M2-S: a lone '-' is a parse error, not INT 0", "[io][m2]") {
    const std::string path = writeTemp("lone_minus.iv", "- 5");
    decimate::Lexer   lexer(path);
    REQUIRE(lexer);
    CHECK_THROWS_AS(lexer.match(static_cast<int>(decimate::Token::INT)), decimate::ParseError);
    std::filesystem::remove(path);
}

TEST_CASE("M2-S: lexer numbers convert to exact values", "[io][m2]") {
    const std::string path = writeTemp("numbers.iv", "-12 3.5 -0.25 1E2 7.");
    decimate::Lexer   lexer(path);
    REQUIRE(lexer);

    lexer.match(static_cast<int>(decimate::Token::INT));
    CHECK(lexer.getInt() == -12);
    lexer.match(static_cast<int>(decimate::Token::FLOAT));
    CHECK(lexer.getFloat() == 3.5F);
    lexer.match(static_cast<int>(decimate::Token::FLOAT));
    CHECK(lexer.getFloat() == -0.25F);
    lexer.match(static_cast<int>(decimate::Token::FLOAT));
    CHECK(lexer.getFloat() == 100.0F);
    lexer.match(static_cast<int>(decimate::Token::FLOAT));
    CHECK(lexer.getFloat() == 7.0F);
    std::filesystem::remove(path);
}

// ===========================================================================
// M3 — SMat: every pivot validated, singular -> false
// ===========================================================================

TEST_CASE("M3-S: forwardEliminate rejects a matrix whose LAST pivot is empty", "[math][m3]") {
    // [[1,0],[0,0]]: pivot 0 is fine, the last row has no entry at all.
    decimate::math::SMat a(2, 2);
    a.insert(0, 0, 1.0);

    decimate::math::SMat aug(2, 4);
    aug.copy(a);
    CHECK_FALSE(aug.forwardEliminate(aug, aug.rows.data()));
}

TEST_CASE("M3-S: inverse of a singular matrix returns without a null dereference", "[math][m3]") {
    // Rank 1 with an empty last row: backSubstitute used to dereference it.
    decimate::math::SMat a(2, 2);
    a.insert(0, 0, 1.0);
    a.insert(0, 1, 2.0);
    CHECK_NOTHROW(a.inverse(nullptr));
}

TEST_CASE("M3-S: inverse of a regular matrix is still exact", "[math][m3]") {
    // [[2,0],[0,4]]^-1 = [[0.5,0],[0,0.25]]
    decimate::math::SMat a(2, 2);
    a.insert(0, 0, 2.0);
    a.insert(1, 1, 4.0);
    const decimate::math::SMat inv = a.inverse(nullptr);
    CHECK(inv.getvalue(0, 0) == 0.5);
    CHECK(inv.getvalue(1, 1) == 0.25);
    CHECK(inv.getvalue(0, 1) == 0.0);
    CHECK(inv.getvalue(1, 0) == 0.0);
}

// ===========================================================================
// M4 — Quat <-> Mat4 round trip
// ===========================================================================

TEST_CASE("M4-S: fromMatrix(toMatrix(q)) equals q up to sign on 100 random unit quaternions",
          "[math][m4]") {
    std::mt19937                     rng(24003);
    std::normal_distribution<double> gauss(0.0, 1.0);

    for (int i = 0; i < 100; ++i) {
        const decimate::math::Quatd q =
            decimate::math::Quatd { gauss(rng), gauss(rng), gauss(rng), gauss(rng) }.normalized();
        const decimate::math::Quatd r = decimate::math::Quatd::fromMatrix(q.toMatrix());

        const double dot = q.x * r.x + q.y * r.y + q.z * r.z + q.w * r.w;
        INFO("sample " << i);
        CHECK(std::abs(dot) == Catch::Approx(1.0).margin(1e-9));
    }
}


TEST_CASE("M4-S: lookRotation still maps +z onto the forward direction", "[math][m4]") {
    // lookRotation builds its matrix for fromMatrix(); it must keep working when the
    // column-major read is fixed.
    const decimate::math::Vec3d forward { 1.0, 2.0, 3.0 };
    const decimate::math::Vec3d n = forward.normalized();
    const decimate::math::Vec3d v =
        decimate::math::Quatd::lookRotation(forward).rotate(decimate::math::Vec3d { 0.0, 0.0, 1.0 });

    CHECK(v.x == Catch::Approx(n.x).margin(1e-12));
    CHECK(v.y == Catch::Approx(n.y).margin(1e-12));
    CHECK(v.z == Catch::Approx(n.z).margin(1e-12));
}
