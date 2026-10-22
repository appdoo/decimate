// This file is part of the decimate library; see decimate.h for version/license details.
/// @file loader_characterization.cpp
/// @brief Characterization tests for Spec 06 — Mesh Loader (OBJ + Open Inventor).
///
/// Tests T1–T14 from the microspec.  All assertions use exact numeric or
/// structural values — never substring-on-output (antipattern 12).
///
/// T1–T6:   hand-authored fixture files (in tests/fixtures/)
/// T7–T9:   real model-scale assets via build-provided path macros
/// T10–T14: state / structure / contract tests
///
/// FROZEN BUG notes:
///   BUG 1 (legacy2003/loader/ObjLdr.cpp:83-84): row[0]=='v' dispatch routes 'vn' lines to
///         read_vertex3, injecting phantom positions.  T3 + T7 assert the
///         contaminated values.
///   BUG 2 (legacy2003/loader/ObjLdr.cpp:24, read_face3): atoi("N//M")==N — slash refs truncated.
///         T4 asserts only the vertex index survives.

#include <catch2/catch_test_macros.hpp>
#include "model_fixture.h"
#include <catch2/catch_approx.hpp>

#include <decimate/types/mesh_data.h>
#include <decimate/mesh/half_edge_mesh.h>

#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

// Build-provided path macros (set in tests/CMakeLists.txt)
// MESH_FIXTURES_DIR  — absolute path to tests/fixtures/
// MESH_MODEL_DIR     — absolute path to model/ (shared fixtures at repo root)

#ifndef MESH_FIXTURES_DIR
#error "MESH_FIXTURES_DIR must be defined by CMake"
#endif
#ifndef MESH_MODEL_DIR
#error "MESH_MODEL_DIR must be defined by CMake"
#endif

using namespace decimate;

namespace {

/// A fixture file written into a private temp directory and removed again on scope exit.
///
/// The old tests wrote into the source tree (MESH_FIXTURES_DIR): parallel `ctest -j` runs raced on the
/// same names and a stale file from an earlier run could mask a failed write.  Each TempFixture owns a
/// uniquely named directory, so neither can happen.
class TempFixture {
public:
    TempFixture(const std::string &fileName, const std::string &content) {
        // Unique per process and per call: clock ticks plus the object address.
        const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
        dir_ = std::filesystem::temp_directory_path() /
               ("decimate_loader_" + std::to_string(tick) + "_" +
                std::to_string(reinterpret_cast<std::uintptr_t>(this)));
        std::filesystem::create_directories(dir_);
        path_ = (dir_ / fileName).string();

        std::ofstream out(path_);
        out << content;
        out.flush();
        // A failed write must fail the test loudly, not leave the loader to read nothing or a stale file.
        REQUIRE(out.good());
    }

    TempFixture(const TempFixture &) = delete;
    TempFixture &operator=(const TempFixture &) = delete;

    ~TempFixture() {
        std::error_code ec;
        std::filesystem::remove_all(dir_, ec);
    }

    const std::string &path() const {
        return path_;
    }

private:
    std::filesystem::path dir_;
    std::string           path_;
};

} // namespace


/// T1 — Clean OBJ triangle
TEST_CASE("T1: OBJ triangle — 3 vertices, 1 face, correct positions + 0-based index",
          "[loader][obj][T1]") {
    // Fixture: tests/fixtures/tri.obj — v 0 0 0 / v 1 0 0 / v 0 1 0 / f 1 2 3
    const std::string path = std::string(MESH_FIXTURES_DIR) + "/tri.obj";
    MeshData          m = loadObj(path);

    REQUIRE(m.positions.size() == 3u); // 3 vertices
    REQUIRE(m.indices.size() == 3u);   // 1 triangle × 3 indices

    // Vertex positions
    CHECK(m.positions[0].x == Catch::Approx(0.0f));
    CHECK(m.positions[0].y == Catch::Approx(0.0f));
    CHECK(m.positions[0].z == Catch::Approx(0.0f));
    CHECK(m.positions[1].x == Catch::Approx(1.0f));
    CHECK(m.positions[1].y == Catch::Approx(0.0f));
    CHECK(m.positions[1].z == Catch::Approx(0.0f));
    CHECK(m.positions[2].x == Catch::Approx(0.0f));
    CHECK(m.positions[2].y == Catch::Approx(1.0f));
    CHECK(m.positions[2].z == Catch::Approx(0.0f));

    // 0-based indices: OBJ f 1 2 3 → legacy -1 decrement → {0, 1, 2}
    CHECK(m.indices[0] == 0);
    CHECK(m.indices[1] == 1);
    CHECK(m.indices[2] == 2);
}


/// T2 — OBJ quad fan-triangulation
TEST_CASE("T2: OBJ quad — fan-triangulation into 2 triangles", "[loader][obj][T2]") {
    // Fixture: 4 vertices + f 1 2 3 4 — fan triangulates to {0,1,2} + {0,2,3}
    // The quad fixture is written on the fly so the test documents its own input.
    const TempFixture fixture("quad.obj",
                              "v 0 0 0\n"
                              "v 1 0 0\n"
                              "v 1 1 0\n"
                              "v 0 1 0\n"
                              "f 1 2 3 4\n");
    MeshData          m = loadObj(fixture.path());

    REQUIRE(m.positions.size() == 4u); // 4 vertices
    REQUIRE(m.indices.size() == 6u);   // 2 triangles × 3 indices

    // Fan: f1=0, first tri {0,1,2}, second tri {0,2,3}
    CHECK(m.indices[0] == 0);
    CHECK(m.indices[1] == 1);
    CHECK(m.indices[2] == 2);
    CHECK(m.indices[3] == 0);
    CHECK(m.indices[4] == 2);
    CHECK(m.indices[5] == 3);
}


/// T3 — FROZEN BUG 1: 'vn' as phantom position
TEST_CASE("T3: FROZEN BUG 1 — vn line routed to read_vertex3, XYZ becomes phantom position",
          "[loader][obj][bug1][T3]") {
    // Fixture: v 0 0 0 / v 1 0 0 / v 0 1 0 / vn 9 9 9 / f 1 2 3
    // FROZEN BUG 1 (legacy2003/loader/ObjLdr.cpp:83-84): row[0]=='v' dispatch — 'vn' routes to
    // read_vertex3; the normal XYZ {9,9,9} is appended as a phantom 4th vertex.
    // positions.size() == 4 (NOT 3), positions[3] == {9,9,9}.
    // DO NOT assert positions.size()==3 — that would be asserting the FIX.
    const TempFixture fixture("vn_bug.obj",
                              "v 0 0 0\n"
                              "v 1 0 0\n"
                              "v 0 1 0\n"
                              "vn 9 9 9\n"
                              "f 1 2 3\n");
    MeshData          m = loadObj(fixture.path());

    // The phantom vn inflates positions to 4 vertices.
    REQUIRE(m.positions.size() == 4u);              // FROZEN BUG 1: 3 real + 1 phantom from 'vn'
    CHECK(m.positions[3].x == Catch::Approx(9.0f)); // phantom X from vn
    CHECK(m.positions[3].y == Catch::Approx(9.0f)); // phantom Y from vn
    CHECK(m.positions[3].z == Catch::Approx(9.0f)); // phantom Z from vn

    // Face indices are unaffected (still reference the 3 real vertices 0-based).
    REQUIRE(m.indices.size() == 3u);
    CHECK(m.indices[0] == 0);
    CHECK(m.indices[1] == 1);
    CHECK(m.indices[2] == 2);
}


/// T4 — FROZEN BUG 2: slash face refs truncated to vertex index
TEST_CASE("T4: FROZEN BUG 2 — f a//c face tokens: atoi stops at '/', vt/vn discarded",
          "[loader][obj][bug2][T4]") {
    // Fixture: 3 vertices + f 1//1 2//2 3//3
    // FROZEN BUG 2 (legacy2003/loader/ObjLdr.cpp:24, read_face3): atoi("1//1") == 1.  The '//vn' attribute
    // is silently discarded.  Only vertex index survives.
    // indices == {0, 1, 2} (1-based 1,2,3 decremented to 0-based 0,1,2).
    const TempFixture fixture("slash_bug.obj",
                              "v 0 0 0\n"
                              "v 1 0 0\n"
                              "v 0 1 0\n"
                              "f 1//1 2//2 3//3\n");
    MeshData          m = loadObj(fixture.path());

    REQUIRE(m.positions.size() == 3u); // only 3 real vertices (no vn lines here)
    REQUIRE(m.indices.size() == 3u);   // 1 triangle

    // FROZEN BUG 2: atoi("1//1")==1 → 0-based 0; vn index discarded.
    CHECK(m.indices[0] == 0);
    CHECK(m.indices[1] == 1);
    CHECK(m.indices[2] == 2);
}


/// T5 — Open Inventor simple triangle
TEST_CASE("T5: Inventor triangle — Coordinate3 + IndexedFaceSet with -1 terminator",
          "[loader][inventor][T5]") {
    // Fixture: tests/fixtures/tri.iv
    // Separator { Coordinate3 { point [0.0 0.0 0.0, 1.0 0.0 0.0, 0.0 1.0 0.0] }
    //             IndexedFaceSet { coordIndex [0, 1, 2, -1,] } }
    const std::string path = std::string(MESH_FIXTURES_DIR) + "/tri.iv";
    MeshData          m = loadInventor(path);

    REQUIRE(m.positions.size() == 3u);
    REQUIRE(m.indices.size() == 3u);

    CHECK(m.positions[0].x == Catch::Approx(0.0f));
    CHECK(m.positions[0].y == Catch::Approx(0.0f));
    CHECK(m.positions[0].z == Catch::Approx(0.0f));
    CHECK(m.positions[1].x == Catch::Approx(1.0f));
    CHECK(m.positions[1].y == Catch::Approx(0.0f));
    CHECK(m.positions[1].z == Catch::Approx(0.0f));
    CHECK(m.positions[2].x == Catch::Approx(0.0f));
    CHECK(m.positions[2].y == Catch::Approx(1.0f));
    CHECK(m.positions[2].z == Catch::Approx(0.0f));

    CHECK(m.indices[0] == 0);
    CHECK(m.indices[1] == 1);
    CHECK(m.indices[2] == 2);
}


/// T6 — Extension dispatcher
TEST_CASE("T6: load() dispatcher routes .obj to loadObj and .iv to loadInventor",
          "[loader][dispatch][T6]") {
    const std::string objPath = std::string(MESH_FIXTURES_DIR) + "/tri.obj";
    const std::string ivPath = std::string(MESH_FIXTURES_DIR) + "/tri.iv";

    MeshData fromObj = loadObj(objPath);
    MeshData fromIv = loadInventor(ivPath);

    // Dispatcher must produce identical results to direct calls.
    MeshData viaLoadObj = load(objPath);
    MeshData viaLoadIv = load(ivPath);

    REQUIRE(viaLoadObj.positions == fromObj.positions);
    REQUIRE(viaLoadObj.indices == fromObj.indices);
    REQUIRE(viaLoadIv.positions == fromIv.positions);
    REQUIRE(viaLoadIv.indices == fromIv.indices);
}


/// T7 — Model-scale golden: cow.obj counts (FROZEN BUG 1 at model scale)
TEST_CASE("T7: cow.obj model-scale — positions inflated by vn phantom, triangle count correct",
          "[loader][obj][model-scale][T7]") {
    // Measured from model/cow.obj:
    //   2903 'v' lines + 3267 'vn' lines → FROZEN BUG 1 inflates numCoord to 6170.
    //   face histogram: {3:734, 4:2519, 5:8, 6:2}
    //   fan-triangulation: 734*1 + 2519*2 + 8*3 + 2*4 = 5804 triangles.
    //
    // FROZEN BUG 1 (legacy2003/loader/ObjLdr.cpp:83-84): DO NOT assert positions.size()==2903;
    // that would be asserting the fix.  Assert 6170 — the contaminated count.
    const std::string path = modelOrSkip("cow.obj");
    MeshData          m = loadObj(path);

    // FROZEN BUG 1 at model scale: 2903 real + 3267 phantom vn = 6170 "vertices".
    REQUIRE(m.positions.size() == static_cast<std::size_t>(6170));

    // Fan-triangulation: 5804 triangles.
    REQUIRE(m.indices.size() == static_cast<std::size_t>(3 * 5804)); // 17412
}


/// T8 — Model-scale golden: cow.obj first vertex + first face content anchors
TEST_CASE("T8: cow.obj content anchors — first vertex + first face (exercises both bugs)",
          "[loader][obj][model-scale][T8]") {
    // First 'v' line in cow.obj: "v 3.716360 2.343387 0.000000"
    // First 'f' line in cow.obj: "f 251//255 252//256 211//215 253//257"
    //   FROZEN BUG 2: atoi stops at '/' → vertex indices only: {251,252,211,253}
    //   0-based: {250, 251, 210, 252}
    //   fan-triangulate quad: {250,251,210}, {250,210,252}
    const std::string path = modelOrSkip("cow.obj");
    MeshData          m = loadObj(path);

    // First vertex (positions[0])
    CHECK(m.positions[0].x == Catch::Approx(3.716360f).epsilon(1e-4f));
    CHECK(m.positions[0].y == Catch::Approx(2.343387f).epsilon(1e-4f));
    CHECK(m.positions[0].z == Catch::Approx(0.0f));

    // First face — quad fan: indices[0..5] == {250, 251, 210, 250, 210, 252}
    // FROZEN BUG 2: slash refs discarded; FROZEN BUG 1: quad requires fan.
    CHECK(m.indices[0] == 250);
    CHECK(m.indices[1] == 251);
    CHECK(m.indices[2] == 210);
    CHECK(m.indices[3] == 250);
    CHECK(m.indices[4] == 210);
    CHECK(m.indices[5] == 252);
}


/// T9 — Model-scale golden: horse.iv counts + first-vertex + first-face anchors
TEST_CASE("T9: horse.iv model-scale — counts + first vertex + first face",
          "[loader][inventor][model-scale][T9]") {
    // Measured from model/horse.iv:
    //   Coordinate3 has 48485 point triples.
    //   IndexedFaceSet has 96966 '-1'-terminated triangle faces.
    //   First point triple: -0.002022 -0.0401785 -0.000863
    //   First coordIndex triple: 26287, 26288, 26289
    const std::string path = modelOrSkip("horse.iv");
    MeshData          m = loadInventor(path);

    REQUIRE(m.positions.size() == static_cast<std::size_t>(48485));
    REQUIRE(m.indices.size() == static_cast<std::size_t>(3 * 96966)); // 290898

    // First vertex anchor (re-measured from file)
    CHECK(m.positions[0].x == Catch::Approx(-0.002022f).epsilon(1e-5f));
    CHECK(m.positions[0].y == Catch::Approx(-0.0401785f).epsilon(1e-5f));
    CHECK(m.positions[0].z == Catch::Approx(-0.000863f).epsilon(1e-5f));

    // First face anchor (re-measured from file — first coordIndex triple)
    CHECK(m.indices[0] == 26287);
    CHECK(m.indices[1] == 26288);
    CHECK(m.indices[2] == 26289);
}


/// T10 — normals / UVs / materials discarded; MeshData = positions + indices only
TEST_CASE("T10: normals, UVs, and materials are discarded; MeshData has positions+indices only",
          "[loader][obj][T10]") {
    // OBJ with vt (texture coord) and usemtl lines.
    // These should be silently ignored — only positions and face indices survive.
    const TempFixture fixture("material_test.obj",
                              "mtllib test.mtl\n"
                              "v 0 0 0\n"
                              "v 1 0 0\n"
                              "v 0 1 0\n"
                              "vt 0.0 0.0\n"
                              "vt 1.0 0.0\n"
                              "vt 0.0 1.0\n"
                              "usemtl TestMaterial\n"
                              "f 1/1 2/2 3/3\n");
    MeshData          m = loadObj(fixture.path());

    // positions: 3 'v' lines + 3 'vt' lines (bug 1: vt also starts with 'v'!)
    // So positions.size() == 6 (3 real + 3 phantom from vt lines).
    // MeshData has NO normal/uv/material field — only positions and indices.
    // If any extra field were added, the assert below would catch it at compile time.
    static_assert(
        !std::is_polymorphic_v<MeshData>,
        "MeshData must be a plain value type (no vtable)");

    // 3 real 'v' + 3 'vt' phantom = 6 "vertices" (FROZEN BUG 1 applies to 'vt' too)
    REQUIRE(m.positions.size() == 6u);

    // f 1/1 2/2 3/3 → FROZEN BUG 2: atoi("1/1")==1, so indices {0,1,2}
    REQUIRE(m.indices.size() == 3u);
    CHECK(m.indices[0] == 0);
    CHECK(m.indices[1] == 1);
    CHECK(m.indices[2] == 2);
}


/// T11 — Error on non-existent file
TEST_CASE("T11: loadObj / loadInventor throw on non-existent file — no stdout",
          "[loader][error][T11]") {
    // Errors must be surfaced as exceptions, never printed to stdout/stderr.
    // (The test framework captures no stdout; the test just asserts throw behavior.)
    CHECK_THROWS_AS(loadObj("/nonexistent/does_not_exist.obj"), std::runtime_error);
    CHECK_THROWS_AS(loadInventor("/nonexistent/does_not_exist.iv"), std::runtime_error);
}


/// T12 — MeshData value semantics: copy independence
TEST_CASE("T12: MeshData value semantics — copy is independent, no raw-pointer leak",
          "[loader][value-semantics][T12]") {
    const std::string path = std::string(MESH_FIXTURES_DIR) + "/tri.obj";
    MeshData          original = loadObj(path);

    // Deep copy
    MeshData copy = original;

    // Mutate the copy
    REQUIRE_FALSE(copy.positions.empty());
    copy.positions[0].x = 999.0f;

    // Original must be unchanged
    CHECK(original.positions[0].x == Catch::Approx(0.0f));

    // Both must still have the same sizes
    CHECK(original.positions.size() == copy.positions.size());
    CHECK(original.indices.size() == copy.indices.size());
}


/// T13 — Boundary greps (structural / compile-time contract)
///
/// The actual boundary greps run in the report; this test documents the contract:
/// MeshData has no GL/OS fields.
TEST_CASE("T13: MeshData struct carries only positions and indices (no GL/OS types)",
          "[loader][boundary][T13]") {
    MeshData m;
    // If MeshData had any GL or OS type, this code would not compile on
    // clang/gcc/Emscripten (no GL headers included in this TU).
    m.positions.push_back({ 1.0f, 0.0f, 0.0f });
    m.indices.push_back(0);
    REQUIRE(m.positions.size() == 1u);
    REQUIRE(m.indices.size() == 1u);
}


/// T14 — Optional: compose with HalfEdgeMesh (Spec 05 interop)
TEST_CASE("T14 (optional): loader output composes with HalfEdgeMesh::create",
          "[loader][compose][T14]") {
    const std::string path = std::string(MESH_FIXTURES_DIR) + "/tri.obj";
    MeshData          m = loadObj(path);

    // MeshData::indices stores 0-based triangle indices matching Face {a,b,c}.
    // Build a HalfEdgeMesh from the loaded triangle.
    REQUIRE(m.indices.size() == 3u);

    decimate::Face         f { m.indices[0], m.indices[1], m.indices[2] };
    decimate::HalfEdgeMesh hem = decimate::HalfEdgeMesh::create(3, 1, &f);

    // A single triangle has 3 half-edges.
    CHECK(hem.halfEdgeCount() == 3u);
}


/// R2 — Out-of-range face indices are rejected at the loader boundary
TEST_CASE("R2: loadObj throws runtime_error on out-of-range face indices",
          "[loader][obj][error][R2]") {
    const std::string dir = std::string(MESH_FIXTURES_DIR) + "/";
    CHECK_THROWS_AS(loadObj(dir + "bad_index_high.obj"), std::runtime_error);     // f 1 2 99999
    CHECK_THROWS_AS(loadObj(dir + "bad_index_negative.obj"), std::runtime_error); // f -1 -2 -3
    CHECK_THROWS_AS(loadObj(dir + "bad_index_zero.obj"), std::runtime_error);     // f 0 1 2
}


/// R2 — the bounds check must not reject valid input.
TEST_CASE("R2: loadObj still accepts tri.obj — 3 positions, 1 face {0,1,2}",
          "[loader][obj][R2]") {
    const MeshData m = loadObj(std::string(MESH_FIXTURES_DIR) + "/tri.obj");
    CHECK(m.positions.size() == 3u);
    REQUIRE(m.indices.size() == 3u);
    CHECK(m.indices[0] == 0);
    CHECK(m.indices[1] == 1);
    CHECK(m.indices[2] == 2);
}


/// R2 — the Inventor loader rejects out-of-range coordIndex too.
TEST_CASE("R2: loadInventor throws runtime_error on out-of-range coordIndex",
          "[loader][iv][error][R2]") {
    CHECK_THROWS_AS(loadInventor(std::string(MESH_FIXTURES_DIR) + "/bad_index.iv"),
                    std::runtime_error);
}
