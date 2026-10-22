// This file is part of the decimate library; see decimate.h for version/license details.
/// @file half_edge_boundary_test.cpp
/// @brief Tests for the HalfEdgeMesh boundary-walk API (Phase A) —
///        firstBoundaryEdge / nextBoundaryEdge / boundaryPointIndices.
///
/// The boundary is implicit: a half-edge with `sym == kNoHalfedge` is a boundary
/// edge. These verify derivation on closed (no boundary) and open (single loop)
/// meshes, with exact vertex sets and loop length.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstdint>
#include <vector>

#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_types.h>

using namespace decimate;

namespace {

/// Half-edge mesh over a flat triangle index list (three entries per face).
HalfEdgeMesh makeHem(int numPoints, const std::vector<int> &idx) {
    const int numFaces = static_cast<int>(idx.size()) / 3;
    return HalfEdgeMesh::create(numPoints, numFaces,
                                reinterpret_cast<const Face *>(idx.data()));
}

} // namespace


TEST_CASE("HE-B1: closed cube has no boundary", "[halfedge][boundary]") {
    // 8-vertex cube, 12 triangles (the closed mesh used by chart S5).
    const std::vector<int> idx = {
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6,
        0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5
    };
    const HalfEdgeMesh hem = makeHem(8, idx);

    CHECK(hem.firstBoundaryEdge() == kNoHalfedge); // closed ⇒ no boundary edge
    CHECK(hem.boundaryPointIndices().empty());
}


TEST_CASE("HE-B2: single triangle boundary loop is its 3 vertices", "[halfedge][boundary]") {
    const HalfEdgeMesh hem = makeHem(3, { 0, 1, 2 });

    REQUIRE(hem.firstBoundaryEdge() != kNoHalfedge);

    // The loop is unordered in the API contract, so compare sorted.
    std::vector<std::uint32_t> b = hem.boundaryPointIndices();
    REQUIRE(b.size() == 3u);
    std::sort(b.begin(), b.end());
    CHECK(b == std::vector<std::uint32_t> { 0u, 1u, 2u });
}


TEST_CASE("HE-B3: quad (two triangles, shared diagonal) — 4-vertex boundary loop",
          "[halfedge][boundary]") {
    // tris {0,1,2} and {0,2,3} share diagonal 0–2 (the only twinned edge).
    const HalfEdgeMesh hem = makeHem(4, { 0, 1, 2, 0, 2, 3 });

    std::vector<std::uint32_t> b = hem.boundaryPointIndices();
    REQUIRE(b.size() == 4u);
    std::sort(b.begin(), b.end());
    CHECK(b == std::vector<std::uint32_t> { 0u, 1u, 2u, 3u });

    // The boundary is a single closed loop: walking returns to the start in
    // exactly 4 steps (4 boundary edges).
    const HeIndex start = hem.firstBoundaryEdge();
    REQUIRE(start != kNoHalfedge);

    // The step cap keeps a broken nextBoundaryEdge from hanging the test.
    HeIndex h = start;
    int     steps = 0;
    do {
        h = hem.nextBoundaryEdge(h);
        ++steps;
    } while (h != kNoHalfedge && h != start && steps < 100);
    CHECK(h == start); // closed loop
    CHECK(steps == 4);
}
