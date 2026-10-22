// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file straighten_atlas.h
/// @brief Internal: the read-only analysis of an atlas that `straighten` works from.
///
/// Silhouette edges (seam half-edges and 3D boundary half-edges) as UV-vertex graph edges,
/// the charts (face components), the ear faces, and the paper's protection rules.  Not
/// installed, no Eigen.

#include <array>
#include <cstdint>
#include <decimate/seam/straighten.h>
#include <decimate/types/textured_mesh.h>
#include <vector>

namespace decimate::seam::detail {

/// One silhouette half-edge between two UV vertices.
struct SilEdge {
    int u0 = 0;
    int u1 = 0;
    int twin = -1; ///< Index of the opposite seam half-edge; -1 for a 3D boundary edge.
    int face = 0;  ///< Face the half-edge belongs to.
};

/// Read-only topology of the input atlas.
struct Atlas {
    const TexturedMesh           *tm = nullptr;
    std::vector<int>              svOf;          ///< Surface vertex of each UV vertex (-1: unused).
    std::vector<std::vector<int>> copies;        ///< UV vertices of each surface vertex.
    std::vector<SilEdge>          edges;         ///< Seam half-edges (twin-linked) and boundary edges.
    std::vector<std::vector<int>> incident;      ///< Silhouette edges at each UV vertex.
    std::vector<int>              seamCount;     ///< Seam edges at each surface vertex.
    std::vector<int>              boundaryCount; ///< 3D boundary edges at each surface vertex.
    std::vector<int>              chartOfFace;   ///< Chart id per face.
    std::vector<std::vector<int>> chartFaces;
    std::vector<int>              earFaces; ///< Faces with exactly two silhouette half-edges.
    int                           sameNeighbourEars = 0;
};


/// @brief Analyse @p tm (classification through `findSeams`).
/// @throws std::invalid_argument if a UV vertex belongs to two surface vertices.
[[nodiscard]] Atlas analyse(const TexturedMesh &tm);


/// Protection-mask layout: bit r is set when paper rule r+1 protects the UV vertex.  This is the
/// union of the bits of rules (1)-(3); bit 3 (value 8) is rule (4), set only by propagation to
/// vertices that no other rule protected.
inline constexpr std::uint8_t kRuleBits = 0x7;


/// @brief Protection masks for rules (1)-(4).
/// @param forced  Vertices of reverted charts (per UV vertex); protected via rule (4).
[[nodiscard]] std::vector<std::uint8_t> protectionMasks(const Atlas                     &atlas,
                                                        const std::vector<std::uint8_t> &forced);

} // namespace decimate::seam::detail
