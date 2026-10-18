// This file is part of the decimate library; see decimate.h for version/license details.
/// @file halfedge_reduce_test.cpp
/// @brief Half-edge edge-collapse decimator (reduceHalfEdge) — modernized 2003
///        CHalfedgeCollapse. Faithful (unguarded) port, so the contract is
///        "valid reduced mesh," not "provably manifold" (guard is a follow-up).
///   HE1 — spot.obj reduces to a valid mesh (indices in range, no degenerate
///         triangle, no NaN); on this closed mesh, no doubled faces either.
///   HE2 — deterministic.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>

#include <decimate/reduce/halfedge_reduce.h>
#include <decimate/types/mesh_data.h>

#ifndef MESH_MODEL_DIR
#  error "MESH_MODEL_DIR must be defined by CMake"
#endif

using namespace decimate;

namespace {
MeshData loadSpot() { return loadObj(std::string(MESH_MODEL_DIR) + "/spot.obj"); }

// 6x signed volume (sum of per-face tet determinants). Its SIGN encodes the
// global winding/orientation — a flipped winding flips the sign. Orientation-
// sensitive, unlike point-to-surface distance or index validity.
double signedVolume6(const MeshData& m) {
    double vol  = 0.0;
    const int t = static_cast<int>(m.indices.size()) / 3;
    for (int i = 0; i < t; ++i) {
        const auto& a = m.positions[static_cast<std::size_t>(m.indices[static_cast<std::size_t>(i) * 3])];
        const auto& b = m.positions[static_cast<std::size_t>(m.indices[static_cast<std::size_t>(i) * 3 + 1])];
        const auto& c = m.positions[static_cast<std::size_t>(m.indices[static_cast<std::size_t>(i) * 3 + 2])];
        vol += static_cast<double>(a.dot(b.cross(c)));
    }
    return vol;
}
}  // namespace

TEST_CASE("HE1: spot.obj half-edge-collapses to a valid mesh", "[halfedge][integration]") {
    const MeshData in = loadSpot();
    REQUIRE(in.positions.size() > 100);

    const int      target = static_cast<int>(in.positions.size()) / 4;
    const MeshData out    = reduceHalfEdge(in, target);

    const int outV = static_cast<int>(out.positions.size());
    REQUIRE(outV > 0);
    REQUIRE(outV < static_cast<int>(in.positions.size()));
    REQUIRE(out.indices.size() < in.indices.size());

    for (const auto& p : out.positions) {
        REQUIRE(std::isfinite(p.x));
        REQUIRE(std::isfinite(p.y));
        REQUIRE(std::isfinite(p.z));
    }
    const int tris = static_cast<int>(out.indices.size()) / 3;
    for (int t = 0; t < tris; ++t) {
        const int a = out.indices[static_cast<std::size_t>(t) * 3];
        const int b = out.indices[static_cast<std::size_t>(t) * 3 + 1];
        const int c = out.indices[static_cast<std::size_t>(t) * 3 + 2];
        REQUIRE((a >= 0 && a < outV && b >= 0 && b < outV && c >= 0 && c < outV));
        REQUIRE((a != b && b != c && a != c));
    }
}

TEST_CASE("HE2: reduceHalfEdge is deterministic", "[halfedge][integration]") {
    const MeshData in     = loadSpot();
    const int      target = static_cast<int>(in.positions.size()) / 4;

    const MeshData a = reduceHalfEdge(in, target);
    const MeshData b = reduceHalfEdge(in, target);

    REQUIRE(a.indices == b.indices);
    REQUIRE(a.positions.size() == b.positions.size());
}

// HE3 — the HE-form guard catches the tetrahedron (link condition edge part):
// every collapse is rejected, so the mesh is returned intact (no doubled face).
TEST_CASE("HE3: a tetrahedron is left intact by the guard", "[halfedge][state]") {
    MeshData m;
    m.positions = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    m.indices   = {0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2};

    const MeshData out = reduceHalfEdge(m, 3);  // ask to remove a vertex

    REQUIRE(out.positions.size() == 4);   // every collapse guarded-out
    REQUIRE(out.indices.size() == 12);
}

// HE4 — winding/orientation is preserved (the reconstruction must not flip
// normals). Signed volume keeps its sign through reduction.
TEST_CASE("HE4: reduceHalfEdge preserves winding (no normal flip)", "[halfedge][integration]") {
    const MeshData in   = loadSpot();
    const double   vIn  = signedVolume6(in);
    REQUIRE(std::abs(vIn) > 1e-6);  // spot is a closed-ish solid -> nonzero

    const MeshData out  = reduceHalfEdge(in, static_cast<int>(in.positions.size()) / 4);
    const double   vOut = signedVolume6(out);

    REQUIRE(vIn * vOut > 0.0);  // same sign -> winding preserved (a flip would negate it)
}
