// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file corner_mesh.h
/// @brief Compact 8-byte half-edge (Directed-Edges / Corner-Table) with a
///        guarded edge-collapse operator — the shared reducer connectivity.
///
/// This is the reusable representation underneath the edge-collapse reducers
/// (`reduceHalfEdge`, and the heap-driven QEM `reduce`). Each half-edge stores
/// only `{ sym, end }` (8 bytes): `next` is arithmetic from the index
/// (`+1`/`+1`/`-2` within the contiguous 3-per-face block) and `face` is `h/3`.
/// This is the **Directed-Edges** representation (Campagna, Kobbelt & Seidel,
/// 1998) / the **Corner Table** (Rossignac et al., ~2001): `end[]` is the
/// corner-vertex (V) table, `sym[]` the opposite (O) table. Half the bytes of a
/// 16-byte/pointer half-edge → ~2× the half-edges per cache line; all `u32`
/// fields 4-byte-aligned → Wasm-clean. Built via O(F) hash pairing.
///
/// A collapse mutates `sym`/`end` in place and tombstones faces (`deadFace`);
/// the triangle list is reconstructed once at the end. The scheduler (heap or
/// error-order) lives in the reducer; this module is the operator + guard only.
///
/// Boundary contract: zero OS/GL dependencies; no global state; Wasm-safe.

#include <decimate/math/vec3.h>
#include <decimate/types/mesh_data.h>

#include <array>
#include <cstdint>
#include <vector>

namespace decimate::corner {

using HeIndex = std::uint32_t;
inline constexpr HeIndex kNoHE = 0xFFFFFFFFu;

/// @brief 8-byte half-edge: twin + tip vertex. `next`/`face` are derived from h.
struct He8 {
    HeIndex sym = kNoHE;  ///< Opposite half-edge (twin), or kNoHE on a boundary.
    HeIndex end = 0;      ///< Tip (end) vertex of this directed edge.
};

/// @brief `next` in the face loop — pure arithmetic (never stored/mutated):
///        within the 3-block the cyclic successor is +1, +1, then -2.
inline HeIndex nextA(HeIndex h) noexcept { return (h % 3u == 2u) ? h - 2u : h + 1u; }

/// @brief Tip vertex of half-edge h.
inline HeIndex tipOf(const std::vector<He8>& he, HeIndex h) noexcept { return he[h].end; }

/// @brief Origin (start) vertex of half-edge h — the CCW predecessor's tip.
inline HeIndex originOf(const std::vector<He8>& he, HeIndex h) noexcept {
    return he[nextA(nextA(h))].end;
}

/// @brief The three vertices of face `faceBase` (= 3·f) in ORIGINAL CCW winding.
///
/// build() stores corner ends reversed (c, b, a), so this is the single place
/// that knows the convention. INVARIANT: cross(v[1]-v[0], v[2]-v[0]) equals the
/// input face's normal (b-a)×(c-a) — half-edge and face agree on orientation.
inline std::array<HeIndex, 3> faceVertices(const std::vector<He8>& he,
                                           HeIndex faceBase) noexcept {
    return {he[faceBase + 2].end, he[faceBase + 1].end, he[faceBase].end};
}

/// @brief Visit every half-edge incident to (ending at) the vertex `incoming`
///        ends at — both ring directions, boundary-aware.
template <class F>
void forEachIncident(const std::vector<He8>& he, HeIndex incoming, F visit) {
    HeIndex loop = incoming;
    do {
        visit(loop);
        loop = he[nextA(loop)].sym;
    } while (loop != kNoHE && loop != incoming);
    if (loop == incoming) {
        return;  // closed ring
    }
    loop = incoming;  // boundary fan — other direction
    do {
        visit(loop);
        const HeIndex s = he[loop].sym;
        loop            = (s != kNoHE) ? nextA(nextA(s)) : kNoHE;
    } while (loop != kNoHE && loop != incoming);
}

/// @brief Build the 8-byte corner-table from a face mesh via O(F) hash pairing.
std::vector<He8> build(const MeshData& mesh);

/// @brief Collapse half-edge h: merge its tip into p1, tombstone both incident
///        faces (sets their `deadFace` entries). Caller moves p1 to its target.
void collapse(std::vector<He8>& he, HeIndex h, HeIndex p1, std::vector<char>& deadFace);

/// @brief Guard: true iff collapsing h (tip p2 → origin p1, p1 moving to
///        @p target) preserves 2-manifoldness (link condition) and flips/
///        degenerates no incident triangle.
bool collapseLegal(const std::vector<He8>& he, const std::vector<math::Vec3f>& pos, HeIndex h,
                   HeIndex p1, HeIndex p2, const math::Vec3f& target);

/// @brief Compact live vertices and emit live faces (tombstone GC), in original
///        CCW winding (via faceVertices — normals consistent with the input).
MeshData reconstruct(const std::vector<He8>& he, const std::vector<math::Vec3f>& pos,
                     const std::vector<char>& deadFace, const std::vector<char>& deadVert,
                     int numFaces);

}  // namespace decimate::corner
