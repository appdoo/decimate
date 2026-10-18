// This file is part of the decimate library; see decimate.h for version/license details.
/// @file cluster_reduce_test.cpp
/// @brief Tests for the octree vertex-clustering decimator (reduceClustering).
///
/// Strategy: hand-meshes for exact structural behavior (separation, coplanar
/// collapse, clamp), spot.obj (committed, CC0) for at-scale integration. All
/// assertions are structural/numeric — no text-presence checks (antipattern 12).

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>
#include <vector>

#include <decimate/reduce/cluster_reduce.h>
#include <decimate/types/mesh_data.h>

#ifndef MESH_MODEL_DIR
#  error "MESH_MODEL_DIR must be defined by CMake"
#endif

using namespace decimate;

namespace {

MeshData makeMesh(std::vector<float> xyz, std::vector<int> tris) {
    MeshData m;
    for (std::size_t i = 0; i + 3 <= xyz.size(); i += 3) {
        m.positions.push_back({xyz[i], xyz[i + 1], xyz[i + 2]});
    }
    m.indices = std::move(tris);
    return m;
}

MeshData loadSpot() {
    return loadObj(std::string(MESH_MODEL_DIR) + "/spot.obj");
}

}  // namespace

// C1 [S] — two triangles far apart stay two distinct representatives; they must
// not merge into one point even under maximal clustering.
TEST_CASE("C1: separated geometry does not merge", "[cluster][state]") {
    // Triangle A near x=0, triangle B near x=1000.
    MeshData m = makeMesh(
        {0, 0, 0, 1, 0, 0, 0, 1, 0,           // A
         1000, 0, 0, 1001, 0, 0, 1000, 1, 0},  // B
        {0, 1, 2, 3, 4, 5});

    const MeshData out = reduceClustering(m, 2);

    // The gap (x in 1..1000) must stay empty: no cell may straddle both
    // triangles, so every representative sits in A's or B's region — clamping
    // to a cell's vertex AABB guarantees it. Exact cluster count depends on how
    // the octree splits each triangle, so we assert the separation, not a count.
    bool nearOrigin = false;
    bool nearFar    = false;
    for (const auto& p : out.positions) {
        REQUIRE((p.x <= 1.001F || p.x >= 999.999F));  // never stranded mid-gap
        nearOrigin = nearOrigin || (p.x <= 1.001F);
        nearFar    = nearFar || (p.x >= 999.999F);
    }
    REQUIRE(nearOrigin);   // the origin triangle survived
    REQUIRE(nearFar);      // the far triangle survived
    REQUIRE(out.indices.empty());  // each triangle collapsed within itself -> dropped
}

// C2 [S] — coplanar vertices in one cell collapse to a single representative
// that lies on their shared plane.
TEST_CASE("C2: coplanar cluster collapses onto its plane", "[cluster][state]") {
    // Unit quad in the z=0 plane, two triangles.
    MeshData m = makeMesh({0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0},
                          {0, 1, 2, 0, 2, 3});

    const MeshData out = reduceClustering(m, 1);

    REQUIRE(out.positions.size() == 1);
    REQUIRE(out.positions[0].z == Catch::Approx(0.0).margin(1e-5));  // on the plane
    REQUIRE(out.indices.empty());                                    // all faces degenerate
}

// C3 [I] — aggressive clustering of a real mesh yields a valid, reduced mesh:
// fewer verts/faces, every index in range, no degenerate triangle, no NaN.
TEST_CASE("C3: spot.obj clusters to a valid reduced mesh", "[cluster][integration]") {
    const MeshData in = loadSpot();
    REQUIRE(in.positions.size() > 100);

    const int target = static_cast<int>(in.positions.size()) / 4;
    const MeshData out = reduceClustering(in, target);

    const int outV = static_cast<int>(out.positions.size());
    REQUIRE(outV > 0);
    REQUIRE(outV < static_cast<int>(in.positions.size()));            // reduced
    REQUIRE(out.indices.size() < in.indices.size());                  // fewer faces

    for (const auto& p : out.positions) {
        REQUIRE(std::isfinite(p.x));
        REQUIRE(std::isfinite(p.y));
        REQUIRE(std::isfinite(p.z));
    }

    const int triCount = static_cast<int>(out.indices.size()) / 3;
    for (int t = 0; t < triCount; ++t) {
        const std::size_t base = static_cast<std::size_t>(t) * 3;
        const int a = out.indices[base];
        const int b = out.indices[base + 1];
        const int c = out.indices[base + 2];
        REQUIRE((a >= 0 && a < outV));
        REQUIRE((b >= 0 && b < outV));
        REQUIRE((c >= 0 && c < outV));
        REQUIRE((a != b && b != c && a != c));   // no degenerate survived
    }
}

// C4 [S] — clamp invariant: every representative is clamped to its cell box, so
// the output bounding box can never exceed the input bounding box.
TEST_CASE("C4: output stays within the input bounding box", "[cluster][state]") {
    const MeshData in = loadSpot();

    math::Vec3f lo = in.positions[0];
    math::Vec3f hi = in.positions[0];
    for (const auto& p : in.positions) {
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }

    const MeshData out = reduceClustering(in, static_cast<int>(in.positions.size()) / 8);

    const float eps = 1e-4F;
    for (const auto& p : out.positions) {
        REQUIRE(p.x >= lo.x - eps);
        REQUIRE(p.x <= hi.x + eps);
        REQUIRE(p.y >= lo.y - eps);
        REQUIRE(p.y <= hi.y + eps);
        REQUIRE(p.z >= lo.z - eps);
        REQUIRE(p.z <= hi.z + eps);
    }
}

// C5 [I] — determinism: identical input + target -> identical output.
TEST_CASE("C5: reduceClustering is deterministic", "[cluster][integration]") {
    const MeshData in = loadSpot();
    const int target = static_cast<int>(in.positions.size()) / 4;

    const MeshData a = reduceClustering(in, target);
    const MeshData b = reduceClustering(in, target);

    REQUIRE(a.indices == b.indices);
    REQUIRE(a.positions.size() == b.positions.size());
    for (std::size_t i = 0; i < a.positions.size(); ++i) {
        REQUIRE(a.positions[i].x == b.positions[i].x);
        REQUIRE(a.positions[i].y == b.positions[i].y);
        REQUIRE(a.positions[i].z == b.positions[i].z);
    }
}

// C6 [I] — survival: at a mild target the mesh is not annihilated; a substantial
// fraction of faces survives (the decimator produces a usable LOD, not a blob).
TEST_CASE("C6: moderate reduction keeps the mesh usable", "[cluster][integration]") {
    const MeshData in = loadSpot();
    const int inFaces = static_cast<int>(in.indices.size()) / 3;

    const MeshData out = reduceClustering(in, static_cast<int>(in.positions.size()) / 2);
    const int outFaces = static_cast<int>(out.indices.size()) / 3;

    REQUIRE(outFaces > inFaces / 4);                 // mesh largely survives
    REQUIRE(outFaces < inFaces);                     // but was reduced
}
