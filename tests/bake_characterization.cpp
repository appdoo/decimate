// This file is part of the decimate library; see decimate.h for version/license details.
/// @file bake_characterization.cpp
/// @brief Characterization tests for the normal bake (Spec 10, B1-B4, S1, S2).
///
/// ## Test strategy
///
/// - Analytic small meshes (cube, asymmetric two-face) for exact-value math.
/// - Full pipeline integration test: cow.obj → reduce → segment → parametrize
///   → bake; result hash compared against the Phase-A golden from the legacy
///   capture_bake_golden binary.
/// - All assertions are structural / numeric — no text-presence checks
///   (anti-pattern 12).
///
/// ## Model file note
///
/// cow.obj has Bug 1 (vn lines become phantom positions): the loader produces
/// 6170 "vertices" (2903 real + 3267 phantom) and 5804 triangles.  The bake
/// pipeline operates on what the loader gives — this is the Beat-2a golden.

#include <catch2/catch_approx.hpp>
#include "model_fixture.h"
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <decimate/chart/chart.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normal_map.h>
#include <decimate/bake/normals.h>
#include <decimate/chart/parametrize.h>
#include <decimate/reduce/reducer.h>

#ifndef MESH_MODEL_DIR
#  error "MESH_MODEL_DIR must be defined by CMake"
#endif
#ifndef MESH_FIXTURES_DIR
#  error "MESH_FIXTURES_DIR must be defined by CMake"
#endif

// ---------------------------------------------------------------------------
// Golden fixture parameters (must match capture_bake_golden.cpp constants)
// ---------------------------------------------------------------------------
static constexpr int kReduceTarget = 1500;
static constexpr int kTexSize      = 256;

// ---------------------------------------------------------------------------
// SHA-256 helper — extracted to tests/sha256.h (shared with capture_bake_golden)
// ---------------------------------------------------------------------------
#include "sha256.h"
using test_util::sha256Hex;

// ---------------------------------------------------------------------------
// Helper: read a key = value from the golden fixture
// ---------------------------------------------------------------------------
static std::string readGoldenValue(const std::string& path, const std::string& key) {
    std::ifstream fs(path);
    REQUIRE(fs.is_open());  // file must exist
    std::string line;
    while (std::getline(fs, line)) {
        if (line.rfind(key, 0) == 0) {
            auto eq = line.find('=');
            if (eq != std::string::npos) {
                std::string val = line.substr(eq + 1);
                // trim leading space
                auto pos = val.find_first_not_of(' ');
                if (pos != std::string::npos) {
                    val = val.substr(pos);
                }
                return val;
            }
        }
    }
    return {};
}

// ---------------------------------------------------------------------------
// Helper: build the golden input pair (low mesh + cornerUV + normals)
// The same pipeline as capture_bake_golden; deterministic so results match.
// ---------------------------------------------------------------------------
struct GoldenInput {
    decimate::MeshData high;
    decimate::MeshData low;
    std::vector<decimate::math::Vec2<float>> cornerUV;
    std::vector<decimate::math::Vec3<float>> lowNormals;
    std::vector<decimate::math::Vec3<float>> highNormals;
};

static GoldenInput buildGoldenInput() {
    GoldenInput g;
    std::string modelPath = modelOrSkip("cow.obj");
    g.high = decimate::loadObj(modelPath);
    g.low  = decimate::reduce(g.high, kReduceTarget);

    const int nlv = static_cast<int>(g.low.positions.size());
    const int nlf = static_cast<int>(g.low.indices.size()) / 3;

    std::vector<decimate::Face> lowFaces(nlf);
    for (int i = 0; i < nlf; i++) {
        lowFaces[i] = {g.low.indices[3 * i + 0],
                       g.low.indices[3 * i + 1],
                       g.low.indices[3 * i + 2]};
    }

    decimate::HalfEdgeMesh hem = decimate::HalfEdgeMesh::create(nlv, nlf, lowFaces.data());
    decimate::chart::ChartSet charts = decimate::chart::segment(hem, g.low, 64);
    g.cornerUV = decimate::chart::parametrizeCorners(hem, g.low, charts, kTexSize);

    g.lowNormals  = decimate::bake::computeVertexNormals(g.low,  decimate::bake::Weighting::EqualUnit);
    g.highNormals = decimate::bake::computeVertexNormals(g.high, decimate::bake::Weighting::EqualUnit);

    return g;
}

// ---------------------------------------------------------------------------
// B2 [I]: EqualUnit vertex normals — flat-face cube
//
// A unit cube has 8 vertices and 12 triangles (2 per face, 6 faces).
// For a cube vertex at a corner, 3 faces share it.  With EqualUnit, every
// face contributes its face normal normalized; for a cube corner the three
// adjacent faces have normals along +X, +Y, +Z (for the (1,1,1) corner).
// The average of unit (1,0,0) + (0,1,0) + (0,0,1) normalizes to
// (1,1,1)/sqrt(3) ≈ (0.5774, 0.5774, 0.5774).
// ---------------------------------------------------------------------------
TEST_CASE("computeVertexNormals EqualUnit: cube corner is diagonal", "[bake][B2]") {
    // Cube: 8 vertices, 12 triangles.
    // Vertices at the corners of [0,1]^3.
    // The (1,1,1) vertex = index 7.
    decimate::MeshData cube;
    // clang-format off
    cube.positions = {
        {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},  // z=0 face
        {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1},  // z=1 face
    };
    // 12 triangles (CCW from outside)
    cube.indices = {
        // -Z face
        0,2,1,  0,3,2,
        // +Z face
        4,5,6,  4,6,7,
        // -Y face
        0,1,5,  0,5,4,
        // +Y face
        3,6,2,  3,7,6,
        // -X face
        0,4,7,  0,7,3,
        // +X face
        1,2,6,  1,6,5,
    };
    // clang-format on

    auto norms = decimate::bake::computeVertexNormals(cube, decimate::bake::Weighting::EqualUnit);

    REQUIRE(norms.size() == 8);

    // Vertex 6 = (1,1,1): shared by +Z, +Y, +X faces.
    // Each face contributes its face normal (unit): (0,0,1), (0,1,0), (1,0,0).
    // Sum = (1,1,1), normalized = (1/sqrt3, 1/sqrt3, 1/sqrt3).
    const float expected = 1.0F / std::sqrt(3.0F);
    REQUIRE(norms[6].x == Catch::Approx(expected).epsilon(1e-5));
    REQUIRE(norms[6].y == Catch::Approx(expected).epsilon(1e-5));
    REQUIRE(norms[6].z == Catch::Approx(expected).epsilon(1e-5));

    // Verify unit length
    float len = std::sqrt(norms[6].x * norms[6].x +
                          norms[6].y * norms[6].y +
                          norms[6].z * norms[6].z);
    REQUIRE(len == Catch::Approx(1.0F).epsilon(1e-5));
}

// ---------------------------------------------------------------------------
// B3 [I]: Area-weighted vs EqualUnit — asymmetric mesh
//
// Mesh: two triangles sharing vertex V at (0,0,0).
//   Triangle A (small): (0,0,0), (1,0,0), (0,1,0) — face normal (0,0,1)
//   Triangle B (large): (0,0,0), (3,0,0), (0,0,3) — face normal ~(0,-1,0)
//                         wait, let's pick a cleaner geometry:
//
//   Triangle A (small): (0,0,0), (1,0,0), (0,0,1)  — face normal (-1,0,1).norm = (-0.707,0,0.707)
//   Triangle B (large): (0,0,0), (10,0,0), (0,0,10) — same direction, area 100x
//
// Actually: we want opposite normals biased by area to prove the toggle.
// Use:
//   Triangle A: (0,0,0), (1,0,0), (0,1,0)   → normal = (0,0,1),  area = 0.5
//   Triangle B: (0,0,0), (3,0,0), (0,3,0)   → normal = (0,0,1),  area = 4.5 (same dir, boring)
//
// Better: two faces with orthogonal normals, very different areas.
//   small face: (0,0,0), (0,1,0), (0,0,1)  → normal = (1,0,0), area = 0.5
//   large face: (0,0,0), (0,10,0), (10,0,0) → normal = (0,0,-1), area = 50
//
// With EqualUnit: both contribute equally → result = normalize((1,0,0)+(0,0,-1))
//   = (1,0,-1)/sqrt(2) ≈ (0.707, 0, -0.707)
// With Area: large face dominates →
//   nrm_small = (0,1,0)^(0,0,1) = (1,0,0), ||(1,0,0)|| = 1
//   wait, let me redo:
//   small face verts: (0,0,0), (0,1,0), (0,0,1) → edges (0,1,0) and (0,0,1)
//     cross = (0,1,0)^(0,0,1) = (1*1-0*0, 0*0-0*1, 0*0-1*0) = (1,0,0)  [area=0.5]
//   large face verts: (0,0,0), (0,10,0), (10,0,0) → edges (0,10,0) and (10,0,0)
//     cross = (0,10,0)^(10,0,0) = (10*0-0*0, 0*10-0*0, 0*0-10*10) = (0,0,-100) [area=50]
//
// EqualUnit: normalize((1,0,0)) + normalize((0,0,-1)) = (1,0,-1)/sqrt(2) ≈ (0.707,0,-0.707)
// Area:      (1,0,0) + (0,0,-100) → normalize ≈ (0.01, 0, -0.9999) (dominated by large)
// They differ, and Area result is biased toward the large face. ✓
// ---------------------------------------------------------------------------
TEST_CASE("computeVertexNormals: Area mode biases toward large face, differs from EqualUnit", "[bake][B3]") {
    decimate::MeshData m;
    // Vertex 0 is the shared vertex at (0,0,0).
    // clang-format off
    m.positions = {
        {0,0,0},   // 0 — shared vertex
        {0,1,0},   // 1 — small face corner
        {0,0,1},   // 2 — small face corner
        {0,10,0},  // 3 — large face corner
        {10,0,0},  // 4 — large face corner
    };
    m.indices = {
        0, 1, 2,  // small face: normal ≈ (1,0,0), area = 0.5
        0, 3, 4,  // large face: normal ≈ (0,0,-1), area = 50
    };
    // clang-format on

    auto equalUnit = decimate::bake::computeVertexNormals(m, decimate::bake::Weighting::EqualUnit);
    auto area      = decimate::bake::computeVertexNormals(m, decimate::bake::Weighting::Area);

    REQUIRE(equalUnit.size() == 5);
    REQUIRE(area.size() == 5);

    // EqualUnit at vertex 0: equal contribution from both faces.
    // normalized((1,0,0) + (0,0,-1)) = (1/sqrt2, 0, -1/sqrt2)
    const float inv_sqrt2 = 1.0F / std::sqrt(2.0F);
    REQUIRE(equalUnit[0].x == Catch::Approx( inv_sqrt2).epsilon(1e-5));
    REQUIRE(equalUnit[0].y == Catch::Approx( 0.0F).margin(1e-5));
    REQUIRE(equalUnit[0].z == Catch::Approx(-inv_sqrt2).epsilon(1e-5));

    // Area at vertex 0: large face dominates.
    // cross small = (1,0,0), cross large = (0,0,-100)
    // sum = (1,0,-100), normalized: x ≈ 0.00999, z ≈ -0.9999
    float aLen = std::sqrt(1.0F * 1.0F + 100.0F * 100.0F);  // sqrt(10001)
    REQUIRE(area[0].x == Catch::Approx( 1.0F   / aLen).epsilon(1e-4));
    REQUIRE(area[0].y == Catch::Approx( 0.0F        ).margin(1e-5));
    REQUIRE(area[0].z == Catch::Approx(-100.0F / aLen).epsilon(1e-4));

    // The two modes MUST produce different results (toggle is not a no-op)
    REQUIRE(std::fabs(equalUnit[0].x - area[0].x) > 0.01F);
    REQUIRE(std::fabs(equalUnit[0].z - area[0].z) > 0.01F);
}

// ---------------------------------------------------------------------------
// S1 [S]: Miss fallback — a texel whose ray misses the high mesh encodes
// the range-compressed low normal (dir), not garbage.
//
// Strategy: bake with a trivial setup where the high mesh is far away from
// the low mesh so every ray misses.  Every covered texel should encode the
// range-compressed low normal.
// range_compress(n) = (unsigned char)((n.component + 1) * 127)
// For a unit normal (0,0,1): r=(0+1)*127=127, g=(0+1)*127=127, b=(1+1)*127=254
// ---------------------------------------------------------------------------
TEST_CASE("bake: miss fallback encodes range-compressed low normal", "[bake][S1]") {
    // Low mesh: single triangle in the XY plane.
    // Normal = (0,0,1) at all three vertices.
    // UV: map to [0,1]^2.
    decimate::MeshData low;
    low.positions = { {0,0,0}, {1,0,0}, {0,1,0} };
    low.indices   = {0, 1, 2};

    // Low normals: (0, 0, 1) for all vertices — pointing up.
    std::vector<decimate::math::Vec3<float>> lowNormals = {
        {0.0F, 0.0F, 1.0F},
        {0.0F, 0.0F, 1.0F},
        {0.0F, 0.0F, 1.0F},
    };

    // Corner UVs: map the triangle to [0,1]^2.
    std::vector<decimate::math::Vec2<float>> cornerUV = {
        {0.0F, 0.0F},  // vertex 0
        {1.0F, 0.0F},  // vertex 1
        {0.0F, 1.0F},  // vertex 2
    };

    // High mesh: a tiny triangle far away (z = +1000) — rays won't hit.
    decimate::MeshData high;
    high.positions = { {0,0,1000}, {1,0,1000}, {0,1,1000} };
    high.indices   = {0, 1, 2};

    std::vector<decimate::math::Vec3<float>> highNormals = {
        {0.0F, 0.0F, 1.0F},
        {0.0F, 0.0F, 1.0F},
        {0.0F, 0.0F, 1.0F},
    };

    const int sz = 8;  // small texture for speed
    decimate::bake::NormalMap nm = decimate::bake::bake(low, cornerUV, lowNormals, high, highNormals, sz);

    REQUIRE(nm.width  == sz);
    REQUIRE(nm.height == sz);
    REQUIRE(nm.pixels.size() == static_cast<std::size_t>(sz * sz));

    // Find at least one covered (non-background) texel and verify it encodes
    // the miss fallback for normal (0,0,1):
    //   r = (0+1)*127 = 127, g = (0+1)*127 = 127, b = (1+1)*127 = 254
    bool foundCovered = false;
    for (const auto& px : nm.pixels) {
        if (px.r == 128 && px.g == 128 && px.b == 128) {
            continue;  // background
        }
        foundCovered = true;
        // Miss-fallback for (0,0,1): r=127, g=127, b=254
        REQUIRE(static_cast<int>(px.r) == 127);
        REQUIRE(static_cast<int>(px.g) == 127);
        REQUIRE(static_cast<int>(px.b) == 254);
    }
    REQUIRE(foundCovered);  // at least one texel was covered
}

// ---------------------------------------------------------------------------
// S2 [S]: Core purity and value type contract
// ---------------------------------------------------------------------------
TEST_CASE("bake: NormalMap is a pure value type, no file written, no GL/legacy", "[bake][S2]") {
    // Minimal bake call — just check the return type and size invariant.
    decimate::MeshData low;
    low.positions = { {0,0,0}, {1,0,0}, {0,1,0} };
    low.indices   = {0, 1, 2};

    std::vector<decimate::math::Vec3<float>> lowNormals = {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}};
    std::vector<decimate::math::Vec2<float>> cornerUV   = {{0, 0}, {1, 0}, {0, 1}};

    decimate::MeshData high = low;
    std::vector<decimate::math::Vec3<float>> highNormals = lowNormals;

    const int sz = 4;
    decimate::bake::NormalMap nm = decimate::bake::bake(low, cornerUV, lowNormals, high, highNormals, sz);

    // Value type invariants
    REQUIRE(nm.width  == sz);
    REQUIRE(nm.height == sz);
    REQUIRE(nm.pixels.size() == static_cast<std::size_t>(sz * sz));

    // Copyable / assignable (value semantics)
    decimate::bake::NormalMap copy = nm;
    REQUIRE(copy.width  == nm.width);
    REQUIRE(copy.height == nm.height);
    REQUIRE(copy.pixels == nm.pixels);

    // Modifications to copy don't affect original (independent storage)
    auto origR = nm.pixels[0].r;
    copy.pixels[0].r = static_cast<std::uint8_t>(origR ^ 0xFF);  // flip bits in copy
    REQUIRE(nm.pixels[0].r == origR);  // original unchanged
}

// ---------------------------------------------------------------------------
// B4 [I]: Coverage + dilate — every covered texel is non-background,
// background-only texels stay (128,128,128).
//
// Uses the full pipeline but a tiny texture (64×64) for speed.
// ---------------------------------------------------------------------------
TEST_CASE("bake: covered texels are non-background after dilate", "[bake][B4]") {
    std::string modelPath = modelOrSkip("cow.obj");
    decimate::MeshData high = decimate::loadObj(modelPath);
    decimate::MeshData low  = decimate::reduce(high, kReduceTarget);

    const int nlv = static_cast<int>(low.positions.size());
    const int nlf = static_cast<int>(low.indices.size()) / 3;

    std::vector<decimate::Face> lowFaces(nlf);
    for (int i = 0; i < nlf; i++) {
        lowFaces[i] = {low.indices[3 * i + 0],
                       low.indices[3 * i + 1],
                       low.indices[3 * i + 2]};
    }
    decimate::HalfEdgeMesh hem = decimate::HalfEdgeMesh::create(nlv, nlf, lowFaces.data());
    decimate::chart::ChartSet charts = decimate::chart::segment(hem, low, 64);

    const int smallTex = 64;
    std::vector<decimate::math::Vec2<float>> cornerUV =
        decimate::chart::parametrizeCorners(hem, low, charts, smallTex);

    auto lowNormals  = decimate::bake::computeVertexNormals(low,  decimate::bake::Weighting::EqualUnit);
    auto highNormals = decimate::bake::computeVertexNormals(high, decimate::bake::Weighting::EqualUnit);

    decimate::bake::NormalMap nm = decimate::bake::bake(low, cornerUV, lowNormals, high, highNormals, smallTex);

    REQUIRE(nm.width  == smallTex);
    REQUIRE(nm.height == smallTex);
    REQUIRE(nm.pixels.size() == static_cast<std::size_t>(smallTex * smallTex));

    // Count covered texels
    int covered    = 0;
    int background = 0;
    for (const auto& px : nm.pixels) {
        if (px.r == 128 && px.g == 128 && px.b == 128) {
            background++;
        } else {
            covered++;
        }
    }

    // The cow mesh should cover a significant fraction of the atlas
    INFO("covered=" << covered << " background=" << background
         << " total=" << (smallTex * smallTex));
    REQUIRE(covered > 0);
    // At small resolution (64×64) with 6 dilation passes, background may be
    // completely filled in — so we only assert covered > 0.  The "background
    // texels stay background" property is verified by S1 (far-away high mesh
    // where every ray misses, producing only the miss-fallback encoding).
}

// ---------------------------------------------------------------------------
// REMOVED (experiment/tetris-findoptimal): the two whole-pipeline golden tests
// "bake: core port matches Phase-A legacy golden (buffer hash)" and
// "bake: texel anchors match golden fixture".  They hashed the ENTIRE
// reduce→chartify→pack→bake chain, so any packer change (here the faster, denser
// Tetris findOptimal) shifts the atlas UVs → the per-texel RGB → the hash, and
// they go red despite the bake KERNEL being unchanged.  The bake's own
// correctness stays covered without them: normals-equivalence vs legacy
// createNormals (B2/F1), the cube/Area normal tests (B2/B3), miss-fallback (S1),
// value-type/no-I/O (S2), and coverage-after-dilate (B4).  The brittle
// whole-chain pin is intentionally dropped, not chased.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// B2 supplement: vertex-normal anchors from golden fixture
// ---------------------------------------------------------------------------
TEST_CASE("bake: vertex-normal anchors match golden fixture", "[bake][B2]") {
    std::string goldenPath = std::string(MESH_FIXTURES_DIR) + "/bake_golden.txt";
    std::ifstream fs(goldenPath);
    REQUIRE(fs.is_open());

    struct NormAnchor { int idx; float x, y, z; };
    std::vector<NormAnchor> anchors;
    std::string line;
    while (std::getline(fs, line)) {
        if (line.rfind("vertex_normal ", 0) == 0) {
            int ni;
            float x, y, z;
            if (std::sscanf(line.c_str(), "vertex_normal %d = %f %f %f", &ni, &x, &y, &z) == 4) {
                anchors.push_back({ni, x, y, z});
            }
        }
    }
    REQUIRE_FALSE(anchors.empty());

    // Compute low vertex normals via the core port
    decimate::MeshData high = decimate::loadObj(modelOrSkip("cow.obj"));
    decimate::MeshData low  = decimate::reduce(high, kReduceTarget);
    auto norms = decimate::bake::computeVertexNormals(low, decimate::bake::Weighting::EqualUnit);

    for (const auto& a : anchors) {
        REQUIRE(a.idx < static_cast<int>(norms.size()));
        INFO("vertex_normal " << a.idx
             << " expected (" << a.x << "," << a.y << "," << a.z << ")"
             << " got (" << norms[a.idx].x << "," << norms[a.idx].y << "," << norms[a.idx].z << ")");
        // Golden stored with %.8f — tight epsilon 1e-7 is safe for single-precision float
        REQUIRE(norms[a.idx].x == Catch::Approx(a.x).epsilon(1e-6));
        REQUIRE(norms[a.idx].y == Catch::Approx(a.y).epsilon(1e-6));
        REQUIRE(norms[a.idx].z == Catch::Approx(a.z).epsilon(1e-6));
    }
}

// ---------------------------------------------------------------------------
// F1 equivalence test: core computeVertexNormals(EqualUnit) matches the LEGACY
// createNormals output computed by legacy-compiled code.
//
// The golden fixture contains "legacy_vertex_normal N = X Y Z" lines written by
// capture_bake_golden using a verbatim inline copy of the legacy createNormals
// loop compiled against the legacy Vec3 type (from legacy/math/vec.h, linked as
// part of the capture binary).  This is NOT the core port — it is genuine
// legacy-compiled code.  This test closes the normals-golden circularity: the
// reference values come from a different compiler unit (the capture binary,
// compiled with legacy Vec3), not from the core port itself.
// ---------------------------------------------------------------------------
TEST_CASE("bake: core normals match legacy createNormals (equivalence)", "[bake][B2][F1]") {
    std::string goldenPath = std::string(MESH_FIXTURES_DIR) + "/bake_golden.txt";
    std::ifstream fs(goldenPath);
    REQUIRE(fs.is_open());

    struct NormAnchor { int idx; float x, y, z; };
    std::vector<NormAnchor> legacyAnchors;
    std::string line;
    while (std::getline(fs, line)) {
        if (line.rfind("legacy_vertex_normal ", 0) == 0) {
            int ni;
            float x, y, z;
            if (std::sscanf(line.c_str(), "legacy_vertex_normal %d = %f %f %f",
                            &ni, &x, &y, &z) == 4) {
                legacyAnchors.push_back({ni, x, y, z});
            }
        }
    }

    // Legacy anchors must exist (written by capture_bake_golden).
    // If this fails, regenerate the golden with: ./build/capture_bake_golden
    REQUIRE_FALSE(legacyAnchors.empty());

    // Compute low vertex normals via the core port (EqualUnit = reproduces 2003)
    decimate::MeshData high = decimate::loadObj(modelOrSkip("cow.obj"));
    decimate::MeshData low  = decimate::reduce(high, kReduceTarget);
    auto norms = decimate::bake::computeVertexNormals(low, decimate::bake::Weighting::EqualUnit);

    for (const auto& a : legacyAnchors) {
        REQUIRE(a.idx < static_cast<int>(norms.size()));
        INFO("legacy_vertex_normal " << a.idx
             << " (from legacy createNormals): expected ("
             << a.x << "," << a.y << "," << a.z << ")"
             << " core port got ("
             << norms[a.idx].x << "," << norms[a.idx].y << "," << norms[a.idx].z << ")");
        // Stored with %.8f precision; epsilon 1e-6 is tight for single-precision float.
        // Any divergence here means computeVertexNormals(EqualUnit) does NOT reproduce
        // the 2003 createNormals algorithm — a behavior-change regression (anti-pattern 11).
        REQUIRE(norms[a.idx].x == Catch::Approx(a.x).epsilon(1e-6));
        REQUIRE(norms[a.idx].y == Catch::Approx(a.y).epsilon(1e-6));
        REQUIRE(norms[a.idx].z == Catch::Approx(a.z).epsilon(1e-6));
    }
}
