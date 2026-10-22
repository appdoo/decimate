// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file seam_topology.h
/// @brief Texture-seam topology of a `TexturedMesh` (Liu et al. 2017, section 3).
///
/// `weld()` lifts a mesh plus per-corner UVs into a `TexturedMesh`; `findSeams()`
/// classifies every surface edge.  Seams are defined by the UVs alone -- never by a
/// chart label -- so the same predicate serves the atlas the cutter produced and any
/// later UV post-process.

#include <decimate/math/vec2.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/textured_mesh.h>
#include <vector>

namespace decimate::seam {

/// @brief One directed edge of a face: corner `corner` -> corner `(corner+1) % 3`.
struct HalfEdgeRef {
    int face = 0;   ///< Face index.
    int corner = 0; ///< Start corner (0..2); the edge ends at corner `(corner+1) % 3`.
};


/// @brief A manifold surface edge whose two half-edges carry different UVs.
///
/// `a` runs from the smaller to the larger surface vertex id, `b` is its opposite
/// half-edge in the neighbouring face.  Parameter gamma in [0,1] walks `a` from its
/// start; `b` is traversed end-to-start so both sides meet the same 3D point.
struct SeamPair {
    HalfEdgeRef a;
    HalfEdgeRef b;
    double      length3D = 0.0; ///< Euclidean edge length in object space.
};


/// @brief The UV silhouette: seams, 3D boundary edges and UV foldovers.
///
/// The seam energy (`seamEnergy`) weighs `seams` only; the erase solver additionally
/// treats `boundary` and `foldovers` as places where texels leave the parametrization.
struct SeamTopology {
    std::vector<SeamPair>    seams;     ///< Manifold edges with differing UVs.
    std::vector<HalfEdgeRef> boundary;  ///< Edges with a single incident face.
    std::vector<HalfEdgeRef> foldovers; ///< Matching-UV edges whose two triangles overlap.
};


/// @brief Lift @p mesh and per-corner UVs into a `TexturedMesh`.
///
/// UV vertex identity is (surface vertex id, exact UV bit pattern): two corners share
/// a UV vertex iff they sit on the same surface vertex with bit-identical UVs.  The
/// surface indices are kept verbatim (`faceV == mesh.indices`) and every corner UV is
/// preserved bit-for-bit, so `cornerUV()` of the result equals @p cornerUV.  No chart
/// labels are needed, so UVs from any post-process (straightening) weld the same way.
///
/// @param mesh      Triangle mesh.
/// @param cornerUV  3 UVs per face, entry `3*f+k` = corner `k` of face `f`.
/// @throws std::invalid_argument if `cornerUV.size() != mesh.indices.size()` or the
///         index count is not a multiple of 3.
[[nodiscard]] TexturedMesh weld(const MeshData &mesh, const std::vector<math::Vec2f> &cornerUV);


/// @brief Classify every surface edge of @p tm.
///
/// An undirected surface edge with both half-edges present is a *seam* when the UVs of
/// its two half-edges differ (by UV vertex id AND by UV value); otherwise it is
/// *matching*, and a matching edge whose two triangles lie on the same side of the
/// edge in UV space is a *foldover*.  An edge with one half-edge is *boundary*.
/// Pairs are ordered by (smaller, larger) surface vertex id.  This is the only seam
/// predicate in `src/seam/`.
///
/// @pre `tm` satisfies the `TexturedMesh` invariants.
[[nodiscard]] SeamTopology findSeams(const TexturedMesh &tm);

} // namespace decimate::seam
