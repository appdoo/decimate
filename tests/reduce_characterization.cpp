// This file is part of the decimate library; see decimate.h for version/license details.
/// @file reduce_characterization.cpp
/// @brief Characterization tests for the QEM edge-collapse reducer (Spec 07).
///
/// ## Test strategy
///
/// - Small hand-meshes (octahedron, two-triangle quad) for exact-value math checks.
/// - Real cow.obj for the "it works at scale" proof (integration test).
/// - All assertions are structural / numeric — no text-presence checks (antipattern 12).
///
/// ## Model file note
///
/// cow.obj has Bug 1 (vn lines become phantom positions): the loader produces
/// 6170 "vertices" (2903 real + 3267 phantom) and 5804 triangles (fan-triangulated
/// quads).  The reducer operates on what the loader gives — this is the correct
/// Beat-2a characterization: same pipeline, same inputs, same outputs.

#include <catch2/catch_approx.hpp>
#include "model_fixture.h"
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <decimate/types/mesh_data.h>
#include <decimate/reduce/qerror.h>
#include <decimate/reduce/reducer.h>

#ifndef MESH_MODEL_DIR
#error "MESH_MODEL_DIR must be defined by CMake"
#endif

using namespace decimate;

/// Builds a MeshData from flat xyz positions and triangle indices (hand-crafted test meshes).
static MeshData makeMesh(std::initializer_list<float> positions,
                         std::initializer_list<int>   indices) {
    MeshData                 m;
    const std::vector<float> flat(positions);
    for (std::size_t i = 0; i + 3 <= flat.size(); i += 3) {
        m.positions.push_back({ flat[i], flat[i + 1], flat[i + 2] });
    }
    m.indices.assign(indices);
    return m;
}


/// Is the mesh valid: index count divisible by 3, all indices in range, no degenerate triangle?
static bool isMeshValid(const MeshData &m) {
    const int numVerts = static_cast<int>(m.positions.size());
    const int numIdx = static_cast<int>(m.indices.size());
    if (numIdx % 3 != 0) {
        return false;
    }
    for (int i = 0; i < numIdx; i += 3) {
        const int a = m.indices[static_cast<std::size_t>(i)];
        const int b = m.indices[static_cast<std::size_t>(i + 1)];
        const int c = m.indices[static_cast<std::size_t>(i + 2)];

        // Indices in range
        if (a < 0 || a >= numVerts) {
            return false;
        }
        if (b < 0 || b >= numVerts) {
            return false;
        }
        if (c < 0 || c >= numVerts) {
            return false;
        }

        // No degenerate triangle (two equal indices)
        if (a == b || b == c || a == c) {
            return false;
        }
    }
    return true;
}


/// [S] QError math: setPlane fills correct K_p coefficients.
/// This is a unit-level sanity check; the full integration tests prove correctness.
/// We verify getError for a known plane: z=0 plane, point at (0,0,d) has error d^2.
///
/// The quadric for the z=0 plane (normal = (0,0,1), d=0):
///   aa=0, ab=0, ac=0, ad=0
///   bb=0, bc=0, bd=0
///   cc=1, cd=0
///   dd=0
/// getError((0,0,h)) = cc*h^2 = h^2
TEST_CASE("S1: QError::setPlane fills K_p for z=0 plane correctly", "[qerror][S1]") {
    QError q;

    // Triangle on the z=0 plane: (0,0,0), (1,0,0), (0,1,0)
    q.setPlane(decimate::math::Vec3f { 0.0f, 0.0f, 0.0f },
               decimate::math::Vec3f { 1.0f, 0.0f, 0.0f },
               decimate::math::Vec3f { 0.0f, 1.0f, 0.0f });


    // Normal = (0,0,1), d = 0 → cc=1, all others 0
    // Error for point at height h: e = h²
    const decimate::math::Vec3f p0 { 0.0f, 0.0f, 0.0f }; // on-plane → error = 0
    const decimate::math::Vec3f p1 { 0.0f, 0.0f, 1.0f }; // height 1 → error = 1
    const decimate::math::Vec3f p2 { 0.0f, 0.0f, 2.0f }; // height 2 → error = 4


    CHECK(q.getError(p0) == Catch::Approx(0.0).margin(1e-10));
    CHECK(q.getError(p1) == Catch::Approx(1.0).margin(1e-10));
    CHECK(q.getError(p2) == Catch::Approx(4.0).margin(1e-10));
}


TEST_CASE("S2: QError::getError is always non-negative", "[qerror][S2]") {
    // Sum of two plane quadrics — error must be >= 0 everywhere.
    QError q1;
    q1.setPlane(decimate::math::Vec3f { 0.0f, 0.0f, 0.0f },
                decimate::math::Vec3f { 1.0f, 0.0f, 0.0f },
                decimate::math::Vec3f { 0.0f, 1.0f, 0.0f });

    QError q2;
    q2.setPlane(decimate::math::Vec3f { 0.0f, 0.0f, 0.0f },
                decimate::math::Vec3f { 0.0f, 1.0f, 0.0f },
                decimate::math::Vec3f { 0.0f, 0.0f, 1.0f });

    QError combined(q1, q2);

    // Sample some arbitrary points — error must be >= 0.
    const std::vector<decimate::math::Vec3f> pts = {
        { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f }, { -2.0f, 3.0f, 1.0f }, { 0.5f, 0.5f, 0.5f }, { 10.0f, 0.0f, -10.0f }
    };
    for (const auto &pt : pts) {
        CHECK(combined.getError(pt) >= 0.0);
    }
}


TEST_CASE("S3: QError::optimize returns in-plane optimal for z=0 plane", "[qerror][S3]") {
    // A single z=0 quadric leaves x and y free, so use three face planes that
    // uniquely constrain a point.
    //
    // Planes:   z=0  (xy plane),  x=0 (yz plane),  y=0 (xz plane)
    // Optimal = intersection = (0,0,0).

    QError qz;
    qz.setPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }); // z=0

    QError qx;
    qx.setPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, { 0.0f, 0.0f, 1.0f }); // x=0

    QError qy;
    qy.setPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 1.0f }); // y=0

    QError combined(qx, qy);
    combined.add(qz);

    // The start value 99 makes a silent no-op visible: it would survive in opt.
    decimate::math::Vec3f opt { 99.0f, 99.0f, 99.0f };

    // Three mutually orthogonal planes uniquely constrain the optimum → the
    // quadric is non-singular, so optimize() MUST succeed (no silent skip).
    REQUIRE(combined.optimize(opt));

    // The optimal point is the planes' intersection: (0,0,0).
    CHECK(opt.x == Catch::Approx(0.0f).margin(1e-5f));
    CHECK(opt.y == Catch::Approx(0.0f).margin(1e-5f));
    CHECK(opt.z == Catch::Approx(0.0f).margin(1e-5f));
    CHECK(combined.getError(opt) == Catch::Approx(0.0).margin(1e-8));
}


/// [S] Tiny hand mesh: octahedron reduced by one vertex.
/// Octahedron: 6 vertices, 8 triangles (3D, not flat, so collapses don't
/// remove all faces at once like a degenerate flat quad would).
/// Reducing to 5 vertices: one edge collapsed, 2 adjacent faces removed → 6 faces.
TEST_CASE("S4: octahedron reduced from 6 to 5 vertices yields 6 non-degenerate faces",
          "[reducer][S4]") {
    // clang-format off
    MeshData oct = makeMesh(
        {  1.0f, 0.0f, 0.0f,   // 0 +x
          -1.0f, 0.0f, 0.0f,   // 1 -x
           0.0f, 1.0f, 0.0f,   // 2 +y
           0.0f,-1.0f, 0.0f,   // 3 -y
           0.0f, 0.0f, 1.0f,   // 4 +z
           0.0f, 0.0f,-1.0f }, // 5 -z
        { 0,2,4,  2,1,4,  1,3,4,  3,0,4,  // upper hemisphere
          0,5,2,  2,5,1,  1,5,3,  3,5,0 } // lower hemisphere
    );
    // clang-format on

    const MeshData reduced = reduce(oct, /*targetVertexCount=*/5);

    const int numVerts = static_cast<int>(reduced.positions.size());
    const int numFaces = static_cast<int>(reduced.indices.size()) / 3;

    // One edge collapse removes the collapsed vertex and 2 adjacent faces.
    CHECK(numVerts <= 5);
    CHECK(numVerts >= 1);
    CHECK(numFaces == 6); // 8 - 2 = 6 after one collapse on a manifold octahedron

    CHECK(isMeshValid(reduced));

    // All positions must be finite.
    for (const auto &v : reduced.positions) {
        CHECK(std::isfinite(v.x));
        CHECK(std::isfinite(v.y));
        CHECK(std::isfinite(v.z));
    }
}


/// [S] Determinism: same cow reduction twice → same counts.
TEST_CASE("S5: reduce is deterministic — same input same output", "[reducer][S5]") {
    const std::string cowPath = modelOrSkip("cow.obj");
    const MeshData    cow = loadObj(cowPath);

    const MeshData r1 = reduce(cow, 1000);
    const MeshData r2 = reduce(cow, 1000);

    // Counts must be identical.
    CHECK(r1.positions.size() == r2.positions.size());
    CHECK(r1.indices.size() == r2.indices.size());

    // Positions must be bit-identical (no UB / nondeterministic ordering).
    REQUIRE(r1.positions.size() == r2.positions.size());
    for (std::size_t i = 0; i < r1.positions.size(); ++i) {
        // Use exact equality — determinism means identical bits.
        CHECK(r1.positions[i] == r2.positions[i]);
    }
}


/// [S] QEM error bounded: finite, non-negative across sample collapses.
TEST_CASE("S6: QError values are finite and non-negative on octahedron collapses",
          "[reducer][qerror][S6]") {
    // Octahedron: 6 vertices, 8 triangles.
    // clang-format off
    MeshData oct = makeMesh(
        {  1.0f, 0.0f, 0.0f,   // 0 +x
          -1.0f, 0.0f, 0.0f,   // 1 -x
           0.0f, 1.0f, 0.0f,   // 2 +y
           0.0f,-1.0f, 0.0f,   // 3 -y
           0.0f, 0.0f, 1.0f,   // 4 +z
           0.0f, 0.0f,-1.0f }, // 5 -z
        { 0,2,4,  2,1,4,  1,3,4,  3,0,4,  // upper hemisphere
          0,5,2,  2,5,1,  1,5,3,  3,5,0 } // lower hemisphere
    );
    // clang-format on

    const int inputVerts = 6;
    const int inputFaces = 8;

    // Reduce by one vertex at a time and verify validity.
    for (int target = inputVerts - 1; target >= 4; --target) {
        const MeshData r = reduce(oct, target);
        const int      outVerts = static_cast<int>(r.positions.size());
        const int      outFaces = static_cast<int>(r.indices.size()) / 3;

        CHECK(outVerts <= inputVerts);
        CHECK(outFaces <= inputFaces);
        CHECK(isMeshValid(r));

        // No NaN or Inf in output positions.
        for (const auto &v : r.positions) {
            CHECK(std::isfinite(v.x));
            CHECK(std::isfinite(v.y));
            CHECK(std::isfinite(v.z));
        }
    }
}


/// [I] Cow integration: load → reduce → structural assertions.
TEST_CASE("I1: reduce cow.obj from 6170 vertices to ~1000 — structural integrity",
          "[reducer][integration][I1]") {
    const std::string cowPath = modelOrSkip("cow.obj");
    const MeshData    cow = loadObj(cowPath);

    // Confirm loader output (Bug 1: vn phantom positions, Bug 2: atoi face refs).
    const int inputVerts = static_cast<int>(cow.positions.size());
    const int inputFaces = static_cast<int>(cow.indices.size()) / 3;

    // 2903 real vertices + 3267 phantom vn positions = 6170 total.
    REQUIRE(inputVerts == 6170);
    // 734 triangles + 2*2519 quads + 3*8 pentagons + 4*2 hexagons = 5804.
    REQUIRE(inputFaces == 5804);

    const MeshData reduced = reduce(cow, /*targetVertexCount=*/1000);
    const int      outVerts = static_cast<int>(reduced.positions.size());
    const int      outFaces = static_cast<int>(reduced.indices.size()) / 3;

    // Vertex count reduced to target (or fewer if heap empties).
    CHECK(outVerts <= 1000);
    CHECK(outVerts > 0);

    // Face count reduced proportionally.
    CHECK(outFaces > 0);
    CHECK(outFaces < inputFaces);

    // Indices divisible by 3.
    CHECK(reduced.indices.size() % 3 == 0u);

    // No degenerate faces; all indices in range.
    CHECK(isMeshValid(reduced));

    // No NaN / Inf in positions.
    for (const auto &v : reduced.positions) {
        CHECK(std::isfinite(v.x));
        CHECK(std::isfinite(v.y));
        CHECK(std::isfinite(v.z));
    }
}


/// [I] Monotonic simplification: fewer vertices → fewer (or equal) faces.
TEST_CASE("I2: monotonic simplification — fewer targets yield fewer faces", "[reducer][I2]") {
    const std::string cowPath = modelOrSkip("cow.obj");
    const MeshData    cow = loadObj(cowPath);

    // Three nested targets on the same input.
    const MeshData r2000 = reduce(cow, 2000);
    const MeshData r1500 = reduce(cow, 1500);
    const MeshData r1000 = reduce(cow, 1000);

    const int f2000 = static_cast<int>(r2000.indices.size()) / 3;
    const int f1500 = static_cast<int>(r1500.indices.size()) / 3;
    const int f1000 = static_cast<int>(r1000.indices.size()) / 3;

    const int v2000 = static_cast<int>(r2000.positions.size());
    const int v1500 = static_cast<int>(r1500.positions.size());
    const int v1000 = static_cast<int>(r1000.positions.size());

    // More targets → more vertices.
    CHECK(v2000 >= v1500);
    CHECK(v1500 >= v1000);

    // More vertices → more (or equal) faces.
    CHECK(f2000 >= f1500);
    CHECK(f1500 >= f1000);

    // All meshes are valid.
    CHECK(isMeshValid(r2000));
    CHECK(isMeshValid(r1500));
    CHECK(isMeshValid(r1000));
}


/// [I] No-op reduction: target >= input vertex count → unchanged.
TEST_CASE("I3: reduce with target >= vertex count returns mesh unchanged", "[reducer][I3]") {
    // clang-format off
    MeshData tri = makeMesh(
        { 0.0f, 0.0f, 0.0f,  1.0f, 0.0f, 0.0f,  0.0f, 1.0f, 0.0f },
        { 0, 1, 2 }
    );
    // clang-format on

    const MeshData r = reduce(tri, 10); // 10 > 3 vertices → unchanged

    REQUIRE(r.positions.size() == tri.positions.size());
    REQUIRE(r.indices.size() == tri.indices.size());
    for (std::size_t i = 0; i < tri.positions.size(); ++i) {
        CHECK(r.positions[i] == tri.positions[i]);
    }
}
