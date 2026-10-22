// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file erase_domain.h
/// @brief Internal: which texels the erase solve optimizes (paper section 4.1).

#include "seam_pieces.h"

#include <decimate/seam/seam_topology.h>
#include <decimate/types/textured_mesh.h>
#include <vector>

namespace decimate::seam::detail {

/// @brief Per-texel classification, row-major (`y * width + x`).
///
/// - `unknown`: the texel takes part in the system.  It is either inside a UV triangle
///   (centre in the closed triangle, tolerance 1e-8) or in a bilinear cell crossed by
///   the UV silhouette (seam, boundary and foldover edges).
/// - `free`: silhouette-cell texels outside every UV triangle.  They carry no
///   trustworthy value (the bake never wrote them) so they are not anchored to the input.
/// Every other texel is untouched by the solve.
struct EraseDomain {
    std::vector<char> unknown;
    std::vector<char> free;
    int               freeCount = 0;
};


/// @brief Classify every texel of @p grid against the UV mesh and the silhouette of @p topo.
///
/// @pre `grid` is at least 2x2 and `tm` satisfies the `TexturedMesh` invariants.
/// @returns masks of size `grid.width * grid.height`; `freeCount` counts the `free` texels.
[[nodiscard]] EraseDomain buildEraseDomain(const TexturedMesh &tm, const SeamTopology &topo,
                                           TexelGrid grid);

} // namespace decimate::seam::detail
