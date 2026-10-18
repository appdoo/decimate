// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file half_edge_mesh.h
/// @brief Index-handle half-edge mesh — the connectivity backbone of mesh-core.
///
/// ## Design
///
/// The legacy implementation stored connectivity as raw Halfedge* pointers
/// into a `new Halfedge[3F]` slab.  This modernisation replaces pointers with
/// 32-bit integer *handles* (indices into a flat `std::vector<Halfedge>`).
///
/// ### Why index handles instead of pointers?
/// - **Portability / Wasm:** indices are 4 bytes regardless of pointer width.
/// - **Value semantics:** the mesh is a plain copyable value type; pointers
///   into a heap slab invalidate on copy.
/// - **No UB:** the legacy `h->sym - halfedges` pointer-subtraction is
///   fragile undefined behaviour if the array ever moves.
/// - **Smaller footprint:** 16 bytes/halfedge vs 32 bytes on 64-bit.
///
/// ### API contract (specs 06–10 are built against this)
///
/// - `HalfEdgeMesh::create()` is the sole factory.
/// - Traversal: `firstOutgoing(vertex)` / `nextOutgoing(he)` replace the
///   pointer-walk inlines in the legacy header.
/// - `sym(h)` and `next(h)` return handles; `kNoHalfedge` is the boundary
///   sentinel (was `NULL` in legacy).
/// - The construction *algorithm* is preserved exactly (Symmetric-search
///   pairing, same iteration order) — only the representation changes.
///   Beat-2b sort-based O(F log F) rebuild is explicitly deferred.
///
/// ### Value semantics
///
/// `HalfEdgeMesh` satisfies the Rule of Zero: the compiler-generated copy
/// constructor / assignment do a deep copy of both vectors; independent
/// storage is guaranteed (T9).

#include <cstdint>
#include <vector>

#include <decimate/types/mesh_types.h>

namespace decimate {

// ---------------------------------------------------------------------------
// HeIndex — half-edge handle type
// ---------------------------------------------------------------------------

/// @brief Index handle for a half-edge.  Lightweight; pass by value.
using HeIndex = std::uint32_t;

/// @brief Sentinel value indicating "no half-edge" (boundary / unset).
///
/// Replaces the legacy `NULL` pointer.  Any `sym` or `next` field holding
/// this value denotes a boundary edge or an uninitialized entry.
inline constexpr HeIndex kNoHalfedge = 0xFFFFFFFFu;

// ---------------------------------------------------------------------------
// Halfedge — value object
// ---------------------------------------------------------------------------

/// @brief A single directed half-edge in the mesh.
///
/// Each triangular face `f` owns exactly three half-edges, wired in a CCW
/// `next` loop: `next(next(next(h))) == h`.
///
/// | Field | Semantics |
/// |-------|-----------|
/// | `sym` | Twin half-edge on the shared edge, or `kNoHalfedge` if this edge lies on a boundary. |
/// | `next`| Next half-edge in the CCW face loop. |
/// | `face`| Owning face index (same for all three half-edges of a face). |
/// | `end` | Index of the vertex at the *tip* (end) of this directed edge. |
///
/// **Origin resolution:** the vertex at the *start* of half-edge `h` equals
/// `end` of `next(next(h))` (the CCW predecessor in the triangle loop).
struct Halfedge {
    HeIndex       sym  = kNoHalfedge;  ///< Twin, or kNoHalfedge on boundary.
    HeIndex       next = kNoHalfedge;  ///< Next in CCW face loop.
    std::uint32_t face = 0;            ///< Owning face index.
    std::uint32_t end  = 0;            ///< End-vertex index.
};

// ---------------------------------------------------------------------------
// HalfEdgeMesh
// ---------------------------------------------------------------------------

/// @brief Half-edge mesh — index-handle connectivity over value storage.
///
/// Constructed via `create()`; the mesh is a value type (copyable,
/// independent storage).  Vertex *positions* are not stored here — only
/// connectivity indices.  The loader (Spec 06) holds the position array and
/// passes faces to `create()`.
///
/// ## Invariants (maintained by `create()`)
/// - `halfEdgeCount() == 3 * faceCount()`.
/// - For every half-edge `h` on a closed mesh: `sym(sym(h)) == h`.
/// - For every half-edge `h`: `next(next(next(h))) == h`.
/// - On an open mesh: boundary edges have `sym(h) == kNoHalfedge`.
class HalfEdgeMesh {
public:
    /// @brief Constructs an empty (invalid) mesh.
    ///
    /// A default-constructed mesh holds no half-edges; useful as a target for
    /// copy/move assignment.
    HalfEdgeMesh() = default;

    /// @brief Rule of Zero — compiler-generated copy/move/dtor are correct.
    ///
    /// Copy deep-copies both vectors; independent storage guaranteed (T9).
    HalfEdgeMesh(const HalfEdgeMesh&)            = default;
    HalfEdgeMesh& operator=(const HalfEdgeMesh&) = default;
    HalfEdgeMesh(HalfEdgeMesh&&) noexcept            = default;
    HalfEdgeMesh& operator=(HalfEdgeMesh&&) noexcept = default;
    ~HalfEdgeMesh()                               = default;

    // -----------------------------------------------------------------------
    // Factory
    // -----------------------------------------------------------------------

    /// @brief Build the half-edge mesh from a flat face array.
    ///
    /// Preserves the legacy `Symmetric`-search pairing algorithm exactly —
    /// only the representation changes (pointer → index handle).
    ///
    /// @param numPoints  Number of distinct vertex positions (used to size
    ///                   the internal connectivity scratchpad).
    /// @param numFaces   Number of triangular faces.
    /// @param faces      Pointer to an array of `numFaces` `Face` structs
    ///                   (vertex indices, CCW winding).
    /// @return           A fully wired `HalfEdgeMesh`.
    ///
    /// @pre  `faces` is non-null when `numFaces > 0`.
    /// @pre  Every vertex index in `faces` is in `[0, numPoints)`.
    static HalfEdgeMesh create(int numPoints, int numFaces, const Face* faces);

    // -----------------------------------------------------------------------
    // Accessors — handle-based traversal API
    // -----------------------------------------------------------------------

    /// @brief Total number of half-edges (== 3 × faceCount).
    [[nodiscard]] std::uint32_t halfEdgeCount() const noexcept {
        return static_cast<std::uint32_t>(halfedges_.size());
    }

    /// @brief Number of faces.
    [[nodiscard]] std::uint32_t faceCount() const noexcept {
        return numFaces_;
    }

    /// @brief Number of vertex positions (as given to `create()`).
    [[nodiscard]] std::uint32_t pointCount() const noexcept {
        return numPoints_;
    }

    /// @brief Return the `Halfedge` value for handle `h`.
    ///
    /// @pre  `h < halfEdgeCount()`.
    [[nodiscard]] const Halfedge& halfedge(HeIndex h) const noexcept {
        return halfedges_[h];
    }

    /// @brief Return a mutable reference to the half-edge for handle `h`.
    ///
    /// @pre  `h < halfEdgeCount()`.
    [[nodiscard]] Halfedge& halfedge(HeIndex h) noexcept {
        return halfedges_[h];
    }

    /// @brief Convenience: sym handle of `h`.
    ///
    /// Returns `kNoHalfedge` if `h` is a boundary edge.
    [[nodiscard]] HeIndex sym(HeIndex h) const noexcept {
        return halfedges_[h].sym;
    }

    /// @brief Convenience: next handle in the face loop.
    [[nodiscard]] HeIndex next(HeIndex h) const noexcept {
        return halfedges_[h].next;
    }

    // -----------------------------------------------------------------------
    // Traversal helpers (mirror the legacy free functions)
    // -----------------------------------------------------------------------

    /// @brief First outgoing half-edge from the *incoming* half-edge `h`.
    ///
    /// Equivalent to legacy `first_outgoing(h)` which returned `h->next`.
    /// The returned handle is the next CCW edge in the same face.
    ///
    /// @return `next(h)` — always valid (never sentinel) for a well-formed mesh.
    [[nodiscard]] HeIndex firstOutgoing(HeIndex h) const noexcept {
        return halfedges_[h].next;
    }

    /// @brief Next outgoing half-edge around the vertex, starting from an
    ///        already-outgoing half-edge `h`.
    ///
    /// Equivalent to legacy `next_outgoing(h)`:
    /// - If `h` has a twin (`sym != kNoHalfedge`), returns `next(sym(h))`.
    /// - Otherwise returns `kNoHalfedge` (boundary reached — no more outgoing).
    ///
    /// Use this to walk the one-ring of a vertex:
    /// ```cpp
    /// HeIndex out = mesh.firstOutgoing(incoming);
    /// do {
    ///     // process out ...
    ///     out = mesh.nextOutgoing(out);
    /// } while (out != kNoHalfedge);
    /// ```
    [[nodiscard]] HeIndex nextOutgoing(HeIndex h) const noexcept {
        const HeIndex s = halfedges_[h].sym;
        if (s == kNoHalfedge) {
            return kNoHalfedge;
        }
        return halfedges_[s].next;
    }

    // -----------------------------------------------------------------------
    // Boundary walk (index port of mudvis HalfEdgeMesh::firstBoundaryEdge /
    // nextBoundaryEdge / boundaryPointIndices). The boundary is IMPLICIT — a
    // half-edge with `sym == kNoHalfedge` is a boundary edge — so it is derived
    // by walking, not stored. Retires the need for a cached boundary set.
    // -----------------------------------------------------------------------

    /// @brief First boundary half-edge (the first with `sym == kNoHalfedge`).
    /// @return Its handle, or `kNoHalfedge` if the mesh is closed (no boundary).
    [[nodiscard]] HeIndex firstBoundaryEdge() const noexcept {
        for (HeIndex h = 0; h < static_cast<HeIndex>(halfedges_.size()); ++h) {
            if (halfedges_[h].sym == kNoHalfedge) {
                return h;
            }
        }
        return kNoHalfedge;
    }

    /// @brief Next boundary half-edge along the boundary loop from `h`.
    ///
    /// Fans around `h`'s end-vertex (`next`, then `next(sym(·))`) until it hits
    /// the next half-edge with no twin. Mirrors mudvis `nextBoundaryEdge`.
    ///
    /// @param h A boundary half-edge (`sym(h) == kNoHalfedge`).
    /// @return The next boundary half-edge, or `kNoHalfedge` if none/degenerate.
    [[nodiscard]] HeIndex nextBoundaryEdge(HeIndex h) const noexcept {
        HeIndex n = halfedges_[h].next;
        while (n != kNoHalfedge && n != h) {
            if (halfedges_[n].sym == kNoHalfedge) {
                return n;
            }
            n = halfedges_[halfedges_[n].sym].next;
        }
        return kNoHalfedge;
    }

    /// @brief Ordered vertex indices of the boundary loop containing the first
    ///        boundary edge (the `end` vertex of each boundary half-edge).
    ///
    /// Empty for a closed mesh. NOTE: returns ONE loop (the one containing
    /// `firstBoundaryEdge()`); a mesh with multiple boundary loops (holes) has
    /// the others omitted — matches mudvis `boundaryPointIndices`. For the
    /// per-chart use (build a sub-mesh from one chart's faces, then call this),
    /// that loop is the chart's outer boundary.
    [[nodiscard]] std::vector<std::uint32_t> boundaryPointIndices() const {
        std::vector<std::uint32_t> indices;
        const HeIndex start = firstBoundaryEdge();
        if (start == kNoHalfedge) {
            return indices;  // closed mesh — no boundary
        }
        HeIndex h = start;
        do {
            indices.push_back(halfedges_[h].end);
            h = nextBoundaryEdge(h);
        } while (h != kNoHalfedge && h != start);
        return indices;
    }

    /// @brief Raw read-only access to the flat half-edge array.
    ///
    /// Intended for downstream modules (reduce, parametrize) that iterate
    /// over all half-edges.  Prefer the accessor methods where possible.
    [[nodiscard]] const std::vector<Halfedge>& halfedges() const noexcept {
        return halfedges_;
    }

private:
    std::vector<Halfedge> halfedges_;    ///< Flat value storage (3 × numFaces entries).
    std::uint32_t         numFaces_  = 0;
    std::uint32_t         numPoints_ = 0;

    // -----------------------------------------------------------------------
    // Construction helpers (matching legacy static helpers, representation-only change)
    // -----------------------------------------------------------------------

    /// @brief Wire the three half-edges for face `faceIndex`.
    ///
    /// Preserves legacy `WriteFace()` exactly:
    /// - edge 0 (`ac`): end = faces[f].c, next = edge 1.
    /// - edge 1 (`cb`): end = faces[f].b, next = edge 2.
    /// - edge 2 (`ba`): end = faces[f].a, next = edge 0.
    void writeFace(int faceIndex, const Face* faces);
};

}  // namespace decimate
