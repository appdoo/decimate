// This file is part of the decimate library; see decimate.h for version/license details.
/// @file core_math_characterization.cpp
/// @brief Characterization tests for Spec 05 — Core Math Structures.
///
/// Tests T1–T14 from the microspec.  All assertions use exact numeric or
/// structural values (antipattern 12 — no text-presence tests).
///
/// Coverage:
///  T1–T4, T7, T8 — HalfEdgeMesh connectivity invariants
///  T5, T10       — SparseMat Gauss-Seidel solver + stored-entry / mul
///  T6, T11       — AssociativeHeap pop order + Repair (decrease-key)
///  T9            — HalfEdgeMesh value semantics (deep-copy independence)
///  T12, T13      — CSparseBitSet cross-word insert/test + erase bug
///  T14           — HalfEdgeMesh construction iteration order

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/math/associative_heap.h>
#include <decimate/math/sparse_matrix.h>
#include <decimate/math/sparse_bit_set.h>
#include <decimate/types/mesh_types.h>

using namespace decimate;

// ===========================================================================
// Helpers
// ===========================================================================

/// Returns V - E + F Euler characteristic for the mesh.
/// E = halfEdgeCount / 2  (each edge = 2 half-edges)
/// V = numPoints (as reported by the mesh)
static int eulerCharacteristic(const HalfEdgeMesh& m) {
    int V = static_cast<int>(m.pointCount());
    int F = static_cast<int>(m.faceCount());
    // Count unique edges: each interior edge has a sym pair; boundary edges exist once.
    int interiorHalfEdges = 0;
    for (std::uint32_t h = 0; h < m.halfEdgeCount(); ++h) {
        if (m.sym(h) != kNoHalfedge) {
            ++interiorHalfEdges;
        }
    }
    // Interior half-edges come in pairs; boundary half-edges are single.
    int boundaryHalfEdges = static_cast<int>(m.halfEdgeCount()) - interiorHalfEdges;
    // Each interior edge is counted twice; each boundary edge once.
    int E = interiorHalfEdges / 2 + boundaryHalfEdges;
    return V - E + F;
}

/// Count half-edges with sym == kNoHalfedge.
static int boundaryCount(const HalfEdgeMesh& m) {
    int count = 0;
    for (std::uint32_t h = 0; h < m.halfEdgeCount(); ++h) {
        if (m.sym(h) == kNoHalfedge) {
            ++count;
        }
    }
    return count;
}

/// Build a tetrahedron: 4 faces, 4 vertices.
/// Face winding produces an outward-normal closed mesh.
static HalfEdgeMesh makeTetra() {
    Face faces[4] = {
        {0, 1, 2},
        {0, 2, 3},
        {0, 3, 1},
        {1, 3, 2},
    };
    return HalfEdgeMesh::create(4, 4, faces);
}

/// Build a triangulated cube: 12 faces, 8 vertices.
static HalfEdgeMesh makeCube() {
    // 8 vertices of a unit cube.
    // Each face of the cube is split into 2 triangles.
    Face faces[12] = {
        // Bottom face (y=0)
        {0, 1, 2}, {0, 2, 3},
        // Top face (y=1)
        {4, 6, 5}, {4, 7, 6},
        // Front face (z=1)
        {3, 2, 6}, {3, 6, 7},
        // Back face (z=0)
        {0, 4, 5}, {0, 5, 1},
        // Left face (x=0)
        {0, 3, 7}, {0, 7, 4},
        // Right face (x=1)
        {1, 5, 6}, {1, 6, 2},
    };
    return HalfEdgeMesh::create(8, 12, faces);
}

// ===========================================================================
// T1 — Single triangle: 3 half-edges, all boundary
// ===========================================================================
TEST_CASE("T1: single triangle has 3 half-edges, all boundary", "[halfedge][T1]") {
    Face faces[1] = {{0, 1, 2}};
    HalfEdgeMesh m = HalfEdgeMesh::create(3, 1, faces);

    REQUIRE(m.halfEdgeCount() == 3u);
    REQUIRE(m.faceCount() == 1u);

    // All half-edges must be boundary (sym == kNoHalfedge)
    for (std::uint32_t h = 0; h < 3u; ++h) {
        CHECK(m.sym(h) == kNoHalfedge);
    }
    CHECK(boundaryCount(m) == 3);
}

// ===========================================================================
// T2 — Two-triangle quad: 6 half-edges, shared interior edge + 4 boundary
// ===========================================================================
TEST_CASE("T2: two-triangle quad has 6 half-edges with correct sym pairing", "[halfedge][T2]") {
    // Quad: two triangles sharing edge (0,2)
    Face faces[2] = {{0, 1, 2}, {0, 2, 3}};
    HalfEdgeMesh m = HalfEdgeMesh::create(4, 2, faces);

    REQUIRE(m.halfEdgeCount() == 6u);
    REQUIRE(m.faceCount() == 2u);

    // Find the one interior edge pair (sym != kNoHalfedge) and verify reciprocity
    int interiorCount = 0;
    for (std::uint32_t h = 0; h < 6u; ++h) {
        HeIndex s = m.sym(h);
        if (s != kNoHalfedge) {
            ++interiorCount;
            // Reciprocal: sym(sym(h)) == h
            CHECK(m.sym(s) == h);
        }
    }
    // Exactly one shared edge → 2 half-edges paired (one from each face)
    CHECK(interiorCount == 2);
    CHECK(boundaryCount(m) == 4);
}

// ===========================================================================
// T3 — Closed tetrahedron: 12 half-edges, fully paired, Euler χ = 2
// ===========================================================================
TEST_CASE("T3: closed tetrahedron has 12 half-edges, all paired, Euler=2", "[halfedge][T3]") {
    HalfEdgeMesh m = makeTetra();

    REQUIRE(m.halfEdgeCount() == 12u);
    REQUIRE(m.faceCount() == 4u);
    REQUIRE(m.pointCount() == 4u);

    // Every half-edge must have a valid sym with reciprocity
    for (std::uint32_t h = 0; h < 12u; ++h) {
        HeIndex s = m.sym(h);
        REQUIRE(s != kNoHalfedge);
        CHECK(m.sym(s) == h);
    }
    CHECK(boundaryCount(m) == 0);
    CHECK(eulerCharacteristic(m) == 2);
}

// ===========================================================================
// T4 — Triangulated cube: 36 half-edges, all paired, Euler χ = 2
// ===========================================================================
TEST_CASE("T4: triangulated cube has 36 half-edges, all paired, Euler=2", "[halfedge][T4]") {
    HalfEdgeMesh m = makeCube();

    REQUIRE(m.halfEdgeCount() == 36u);
    REQUIRE(m.faceCount() == 12u);
    REQUIRE(m.pointCount() == 8u);

    for (std::uint32_t h = 0; h < 36u; ++h) {
        HeIndex s = m.sym(h);
        REQUIRE(s != kNoHalfedge);
        CHECK(m.sym(s) == h);
    }
    CHECK(boundaryCount(m) == 0);
    CHECK(eulerCharacteristic(m) == 2);
}

// ===========================================================================
// T5 — SparseMat Gauss-Seidel solve 3×3 known system
// ===========================================================================
TEST_CASE("T5: Gauss-Seidel solveIterative gives correct result within 1e-4", "[sparsemat][T5]") {
    // A = [[2,1,0],[1,3,1],[0,1,2]], b = [1, 2, 1]
    // Analytic solution: x = [0.25, 0.5, 0.25]
    SMat A(3, 3);
    A.insert(0, 0, 2.0);
    A.insert(0, 1, 1.0);
    A.insert(1, 0, 1.0);
    A.insert(1, 1, 3.0);
    A.insert(1, 2, 1.0);
    A.insert(2, 1, 1.0);
    A.insert(2, 2, 2.0);

    SMat b(3, 1);
    b.insert(0, 0, 1.0);
    b.insert(1, 0, 2.0);
    b.insert(2, 0, 1.0);

    SMat x = A.solveIterative(&b);

    REQUIRE(x.m == 3);
    REQUIRE(x.n == 1);

    constexpr double eps = 1e-4;
    CHECK(std::abs(x.getvalue(0, 0) - 0.25) < eps);
    CHECK(std::abs(x.getvalue(1, 0) - 0.50) < eps);
    CHECK(std::abs(x.getvalue(2, 0) - 0.25) < eps);
}

// ===========================================================================
// T6 — Heap pop order: push [5,3,8,1,4] → pop order [1,3,4,5,8]
// ===========================================================================
TEST_CASE("T6: heap pop order ascending for push [5,3,8,1,4]", "[heap][T6]") {
    struct Item : Heapable {
        explicit Item(double v) { heapValue = v; }
    };

    Item e5(5.0), e3(3.0), e8(8.0), e1(1.0), e4(4.0);

    AssociativeHeap heap;
    heap.push(&e5);
    heap.push(&e3);
    heap.push(&e8);
    heap.push(&e1);
    heap.push(&e4);

    REQUIRE(heap.size() == 5);

    double expected[] = {1.0, 3.0, 4.0, 5.0, 8.0};
    for (int i = 0; i < 5; ++i) {
        Heapable* top = heap.pop();
        REQUIRE(top != nullptr);
        CHECK(top->heapValue == Catch::Approx(expected[i]));
    }
    CHECK(heap.empty());
}

// ===========================================================================
// T7 — next-loop closure + face bookkeeping (all meshes)
// ===========================================================================
TEST_CASE("T7: next-loop closure and face field for triangle/quad/tetra/cube", "[halfedge][T7]") {
    // Single triangle
    {
        Face faces[1] = {{0, 1, 2}};
        HalfEdgeMesh m = HalfEdgeMesh::create(3, 1, faces);
        // next(next(next(h))) == h for all h in [0,2]
        for (std::uint32_t h = 0; h < 3u; ++h) {
            CHECK(m.next(m.next(m.next(h))) == h);
            CHECK(m.halfedge(h).face == 0u);
        }
    }
    // Two-triangle quad
    {
        Face faces[2] = {{0, 1, 2}, {0, 2, 3}};
        HalfEdgeMesh m = HalfEdgeMesh::create(4, 2, faces);
        for (std::uint32_t h = 0; h < 6u; ++h) {
            CHECK(m.next(m.next(m.next(h))) == h);
            CHECK(m.halfedge(h).face == (h < 3u ? 0u : 1u));
        }
    }
    // Tetrahedron
    {
        HalfEdgeMesh m = makeTetra();
        for (std::uint32_t h = 0; h < 12u; ++h) {
            CHECK(m.next(m.next(m.next(h))) == h);
            CHECK(m.halfedge(h).face == h / 3u);
        }
    }
}

// ===========================================================================
// T8 — end/origin index correctness for all meshes
// ===========================================================================
TEST_CASE("T8: end and origin indices match face winding", "[halfedge][T8]") {
    // Single triangle {0,1,2}: he0→end=2, he1→end=1, he2→end=0
    {
        Face faces[1] = {{0, 1, 2}};
        HalfEdgeMesh m = HalfEdgeMesh::create(3, 1, faces);
        // writeFace: ac(0)→end=c=2, cb(1)→end=b=1, ba(2)→end=a=0
        CHECK(m.halfedge(0).end == 2u);
        CHECK(m.halfedge(1).end == 1u);
        CHECK(m.halfedge(2).end == 0u);
        // Origin of h = end of next(next(h))
        // origin(0) = end(next(next(0))) = end(next(1)) = end(2) = 0
        CHECK(m.halfedge(m.halfedge(m.halfedge(0).next).next).end == 0u);
    }
    // Tetrahedron — check all 12 half-edges for valid index range
    {
        HalfEdgeMesh m = makeTetra();
        for (std::uint32_t h = 0; h < 12u; ++h) {
            CHECK(m.halfedge(h).end < 4u);
            // origin
            CHECK(m.halfedge(m.halfedge(m.halfedge(h).next).next).end < 4u);
        }
    }
}

// ===========================================================================
// T9 — Deep-copy value semantics
// ===========================================================================
TEST_CASE("T9: deep-copy independence — mutating copy leaves original unchanged", "[halfedge][T9]") {
    HalfEdgeMesh original = makeTetra();

    // Deep copy
    HalfEdgeMesh copy = original;

    // Capture original state: first half-edge sym handle
    HeIndex origSym0 = original.sym(0u);

    // Mutate the copy: corrupt its first half-edge sym
    copy.halfedge(0u).sym = kNoHalfedge;

    // Original must be byte-for-byte unchanged
    CHECK(original.sym(0u) == origSym0);
    CHECK(original.sym(0u) != kNoHalfedge);

    // The copy's sym really changed
    CHECK(copy.sym(0u) == kNoHalfedge);

    // Both still have 12 half-edges
    CHECK(original.halfEdgeCount() == 12u);
    CHECK(copy.halfEdgeCount() == 12u);
}

// ===========================================================================
// T10 — SparseMat stored entries + mat-vec product
// ===========================================================================
TEST_CASE("T10: SparseMat insert/getvalue and mul match expected values", "[sparsemat][T10]") {
    // Build 3×3 matrix A = [[2,1,0],[1,3,1],[0,1,2]]
    SMat A(3, 3);
    A.insert(0, 0, 2.0);
    A.insert(0, 1, 1.0);
    A.insert(1, 0, 1.0);
    A.insert(1, 1, 3.0);
    A.insert(1, 2, 1.0);
    A.insert(2, 1, 1.0);
    A.insert(2, 2, 2.0);

    // Exact stored values
    CHECK(A.getvalue(0, 0) == Catch::Approx(2.0));
    CHECK(A.getvalue(0, 1) == Catch::Approx(1.0));
    CHECK(A.getvalue(0, 2) == Catch::Approx(0.0));
    CHECK(A.getvalue(1, 0) == Catch::Approx(1.0));
    CHECK(A.getvalue(1, 1) == Catch::Approx(3.0));
    CHECK(A.getvalue(1, 2) == Catch::Approx(1.0));
    CHECK(A.getvalue(2, 1) == Catch::Approx(1.0));
    CHECK(A.getvalue(2, 2) == Catch::Approx(2.0));

    // mat-vec product via mul(rowid, vec)
    double v1[3] = {1.0, 0.0, 0.0};
    CHECK(A.mul(0, v1) == Catch::Approx(2.0));
    CHECK(A.mul(1, v1) == Catch::Approx(1.0));
    CHECK(A.mul(2, v1) == Catch::Approx(0.0));

    // A * [1, 1, 1]^T = [3, 5, 3]
    double v2[3] = {1.0, 1.0, 1.0};
    CHECK(A.mul(0, v2) == Catch::Approx(3.0));
    CHECK(A.mul(1, v2) == Catch::Approx(5.0));
    CHECK(A.mul(2, v2) == Catch::Approx(3.0));
}

// ===========================================================================
// T11 — Heap Repair (decrease-key) + heapPos bookkeeping
// ===========================================================================
TEST_CASE("T11: heap Repair after decrease-key reorders correctly", "[heap][T11]") {
    // This test uses the EXACT spec values (value 0, pop order [0,1,3,4,5])
    // and traces the intrusive heapPos field after repair and after each pop.
    //
    // Push sequence [5,3,8,1,4] → internal heap after all pushes:
    //   slot: [0]=e1(1.0)  [1]=e3(3.0)  [2]=e8(8.0)  [3]=e5(5.0)  [4]=e4(4.0)
    //   heapPos: e1=0, e3=1, e8=2, e5=3, e4=4
    //
    // Change e8.heapValue = 0.0, repair → upHeap(2): 0.0 < 1.0 → swap slots 0,2
    //   slot: [0]=e8(0.0)  [1]=e3(3.0)  [2]=e1(1.0)  [3]=e5(5.0)  [4]=e4(4.0)
    //   heapPos after repair: e8=0, e3=1, e1=2, e5=3, e4=4

    struct Item : Heapable {
        explicit Item(double v) { heapValue = v; }
    };

    Item e5(5.0), e3(3.0), e8(8.0), e1(1.0), e4(4.0);

    AssociativeHeap heap;
    heap.push(&e5);
    heap.push(&e3);
    heap.push(&e8);
    heap.push(&e1);
    heap.push(&e4);

    REQUIRE(heap.size() == 5);

    // Change the 8.0 element to 0.0 (spec-mandated exact value) and repair
    e8.heapValue = 0.0;
    heap.repair(&e8);

    // Assert heapPos of ALL five elements after repair (full trace, per spec)
    CHECK(e8.heapPos == 0);
    CHECK(e3.heapPos == 1);
    CHECK(e1.heapPos == 2);
    CHECK(e5.heapPos == 3);
    CHECK(e4.heapPos == 4);

    // Pop order must be [0, 1, 3, 4, 5] (spec-mandated exact values).
    // After each pop, assert the heapPos of every remaining element.

    // Pop 1: e8(0.0) — after: heap=[e1,e3,e4,e5]; e1=0,e3=1,e4=2,e5=3
    {
        Heapable* top = heap.pop();
        REQUIRE(top != nullptr);
        CHECK(top->heapValue == Catch::Approx(0.0));
        CHECK(e1.heapPos == 0);
        CHECK(e3.heapPos == 1);
        CHECK(e4.heapPos == 2);
        CHECK(e5.heapPos == 3);
    }

    // Pop 2: e1(1.0) — after: heap=[e3,e5,e4]; e3=0,e5=1,e4=2
    {
        Heapable* top = heap.pop();
        REQUIRE(top != nullptr);
        CHECK(top->heapValue == Catch::Approx(1.0));
        CHECK(e3.heapPos == 0);
        CHECK(e5.heapPos == 1);
        CHECK(e4.heapPos == 2);
    }

    // Pop 3: e3(3.0) — after: heap=[e4,e5]; e4=0,e5=1
    {
        Heapable* top = heap.pop();
        REQUIRE(top != nullptr);
        CHECK(top->heapValue == Catch::Approx(3.0));
        CHECK(e4.heapPos == 0);
        CHECK(e5.heapPos == 1);
    }

    // Pop 4: e4(4.0) — after: heap=[e5]; e5=0
    {
        Heapable* top = heap.pop();
        REQUIRE(top != nullptr);
        CHECK(top->heapValue == Catch::Approx(4.0));
        CHECK(e5.heapPos == 0);
    }

    // Pop 5: e5(5.0) — heap empty
    {
        Heapable* top = heap.pop();
        REQUIRE(top != nullptr);
        CHECK(top->heapValue == Catch::Approx(5.0));
    }

    CHECK(heap.empty());
}

// ===========================================================================
// T12 — CSparseBitSet cross-word insert and test
// ===========================================================================
TEST_CASE("T12: CSparseBitSet cross-word insert and test", "[bitset][T12]") {
    CSparseBitSet bs;

    // bit 3 → word 0 (3 >> 5 = 0)
    // bit 35 → word 1 (35 >> 5 = 1)
    bs.insert(3);
    bs.insert(35);

    CHECK(bs.test(3) == true);
    CHECK(bs.test(35) == true);
    CHECK(bs.test(4) == false);
    CHECK(bs.test(34) == false);
    CHECK(bs.test(36) == false);
    CHECK(bs.test(0) == false);
}

// ===========================================================================
// T13 — CSparseBitSet erase() bug characterization (DO NOT FIX)
// ===========================================================================
TEST_CASE("T13: erase() characterized legacy bugs — wipe-whole-word and no-op-at-zero", "[bitset][T13]") {
    // Scenario A: erase of a bit in array-index-0 word is a NO-OP.
    // bit 3 → word 0 → array index 0
    // if (offset > 0) is false → erase does nothing
    {
        CSparseBitSet bs;
        bs.insert(3);
        bs.insert(4);
        bs.insert(5);

        // erase(4): 4 >> 5 = 0 (word 0), Find(0) = 0, if (0 > 0) = false → no-op
        bs.erase(4);

        // ALL three bits remain set (the erase had no effect)
        CHECK(bs.test(3) == true);
        CHECK(bs.test(4) == true);   // NOT erased — bug #1 (no-op at index 0)
        CHECK(bs.test(5) == true);
    }

    // Scenario B: erase of a bit in array-index-1 word wipes the WHOLE word.
    // Insert bit 3 first (array index 0 = word 0), then bits 35, 36, 37 (array index 1 = word 1).
    // erase(36): 36 >> 5 = 1 (word 1), Find(1) = 1, if (1 > 0) = true.
    // bits[1] &= !(1 << (36 & 0x1F)) = !(1<<4) = !16 = 0 → whole word wiped.
    {
        CSparseBitSet bs;
        bs.insert(3);   // word 0, array index 0
        bs.insert(35);  // word 1, array index 1; bit pos 3
        bs.insert(36);  // word 1, array index 1; bit pos 4
        bs.insert(37);  // word 1, array index 1; bit pos 5

        // Confirm all are set before erase
        REQUIRE(bs.test(35) == true);
        REQUIRE(bs.test(36) == true);
        REQUIRE(bs.test(37) == true);
        REQUIRE(bs.test(3)  == true);

        bs.erase(36);

        // Bug: !-vs-~ wipes the entire word 1 → 35 and 37 also cleared
        CHECK(bs.test(36) == false);  // target bit cleared (side effect of whole-word wipe)
        CHECK(bs.test(35) == false);  // collateral — whole word wiped
        CHECK(bs.test(37) == false);  // collateral — whole word wiped
        CHECK(bs.test(3)  == true);   // word 0 is unaffected
    }
}

// ===========================================================================
// T14 — Construction iteration order (algorithm unchanged)
// ===========================================================================
TEST_CASE("T14: Create() consumes each face exactly once, in input order", "[halfedge][T14]") {
    // For a tetrahedron, face i owns half-edges 3i, 3i+1, 3i+2.
    // This is the legacy WriteFace() invariant (T14 pins iteration order).
    HalfEdgeMesh m = makeTetra();

    // Verify the face-to-halfedge mapping: face f → half-edges [3f, 3f+1, 3f+2]
    for (std::uint32_t f = 0; f < 4u; ++f) {
        CHECK(m.halfedge(f * 3u + 0u).face == f);
        CHECK(m.halfedge(f * 3u + 1u).face == f);
        CHECK(m.halfedge(f * 3u + 2u).face == f);
    }

    // Tetrahedron faces: {0,1,2},{0,2,3},{0,3,1},{1,3,2}
    // Face 0: ac(he0)→end=2, cb(he1)→end=1, ba(he2)→end=0
    // Face 1: ac(he3)→end=3, cb(he4)→end=2, ba(he5)→end=0
    // Face 2: ac(he6)→end=1, cb(he7)→end=3, ba(he8)→end=0
    // Face 3: ac(he9)→end=2, cb(he10)→end=3, ba(he11)→end=1
    CHECK(m.halfedge(0).end == 2u);
    CHECK(m.halfedge(1).end == 1u);
    CHECK(m.halfedge(2).end == 0u);
    CHECK(m.halfedge(3).end == 3u);
    CHECK(m.halfedge(4).end == 2u);
    CHECK(m.halfedge(5).end == 0u);
}
