// This file is part of the decimate library; see decimate.h for version/license details.
/// @file half_edge_mesh.cpp
/// @brief HalfEdgeMesh implementation — behavior-preserving port of
///        legacy/math/halfedge.cpp (Spec 05, Beat 2a).
///
/// ## Algorithm provenance
///
/// `WriteFace()` winding is preserved byte-for-byte from legacy. The twin (sym)
/// pairing, however, replaces the legacy O(F·valence) `Symmetric`-search with
/// **O(F) hash pairing** over directed edges (key (origin,tip) → twin (tip,origin)).
/// This is the sort/hash-based construction the findings dossier deferred — now
/// realized. It is NOT a silent swap: verified output-identical to the legacy
/// scan across the full corpus (T14 + every reduce/UV/bake golden green), because
/// on a manifold mesh each edge has a unique twin, so both methods agree.

#include <decimate/mesh/half_edge_mesh.h>

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace decimate {

// ---------------------------------------------------------------------------
// writeFace — mirrors legacy WriteFace()
// ---------------------------------------------------------------------------
// Wire the three half-edges for face `faceIndex` in the CCW loop:
//
//   index+0  (ac): end=c, next=index+1
//   index+1  (cb): end=b, next=index+2
//   index+2  (ba): end=a, next=index+0
//
// The "end" vertex naming follows legacy convention: ac ends at c, cb ends
// at b, ba ends at a.  The sym pointers are left at kNoHalfedge here;
// hash pairing (step 2) fills them in afterwards.
void HalfEdgeMesh::writeFace(int faceIndex, const Face* faces) {
    const HeIndex base = static_cast<HeIndex>(faceIndex) * 3u;
    HeIndex ac = base;
    HeIndex cb = base + 1u;
    HeIndex ba = base + 2u;

    // Face membership
    halfedges_[ac].face = static_cast<std::uint32_t>(faceIndex);
    halfedges_[cb].face = static_cast<std::uint32_t>(faceIndex);
    halfedges_[ba].face = static_cast<std::uint32_t>(faceIndex);

    // End-vertex indices (matches legacy: ac->end = faces[f].c, etc.)
    halfedges_[ac].end = static_cast<std::uint32_t>(faces[faceIndex].c);
    halfedges_[cb].end = static_cast<std::uint32_t>(faces[faceIndex].b);
    halfedges_[ba].end = static_cast<std::uint32_t>(faces[faceIndex].a);

    // CCW next loop
    halfedges_[ac].next = cb;
    halfedges_[cb].next = ba;
    halfedges_[ba].next = ac;

    // sym left at kNoHalfedge (filled by hash pairing below)
}

// ---------------------------------------------------------------------------
// create — mirrors legacy HalfEdgeMesh::Create()
// ---------------------------------------------------------------------------
HalfEdgeMesh HalfEdgeMesh::create(int numPoints, int numFaces, const Face* faces) {
    HalfEdgeMesh mesh;
    if (numFaces <= 0 || numPoints <= 0) {
        return mesh;
    }

    // Allocate flat half-edge storage — all fields default to kNoHalfedge / 0
    const int numHalfEdges = numFaces * 3;
    mesh.halfedges_.resize(static_cast<std::size_t>(numHalfEdges));
    mesh.numFaces_  = static_cast<std::uint32_t>(numFaces);
    mesh.numPoints_ = static_cast<std::uint32_t>(numPoints);

    // Step 1: wire faces without sym (sets next/end; mirrors legacy WriteFace loop).
    for (int i = 0; i < numFaces; ++i) {
        mesh.writeFace(i, faces);
    }

    // Step 2: pair twins by hashing directed edges — O(F), replacing the legacy
    // O(F*valence) connect[]-scan pairing. For a half-edge (origin a ->
    // tip b), its twin is the half-edge (b -> a): look it up, else insert (a->b).
    std::unordered_map<std::uint64_t, HeIndex> edgeMap;
    edgeMap.reserve(static_cast<std::size_t>(numHalfEdges));
    for (int i = 0; i < numHalfEdges; ++i) {
        const HeIndex       h = static_cast<HeIndex>(i);
        const std::uint32_t b = mesh.halfedges_[h].end;  // tip
        const std::uint32_t a =
            mesh.halfedges_[mesh.halfedges_[mesh.halfedges_[h].next].next].end;  // origin
        const std::uint64_t keyBA = (static_cast<std::uint64_t>(b) << 32) | a;
        const auto          it    = edgeMap.find(keyBA);
        if (it != edgeMap.end()) {
            const HeIndex s        = it->second;
            mesh.halfedges_[h].sym = s;
            mesh.halfedges_[s].sym = h;
            edgeMap.erase(it);
        } else {
            const std::uint64_t keyAB = (static_cast<std::uint64_t>(a) << 32) | b;
            edgeMap.emplace(keyAB, h);
        }
    }

    return mesh;
}

}  // namespace decimate
