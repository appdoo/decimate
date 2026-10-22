// This file is part of the decimate library; see decimate.h for version/license details.
/// @file seam_energy.cpp
/// @brief seamEnergy -- D_total of paper eq. 9 over the exact pieces of seam_pieces.h.

#include "seam_pieces.h"

#include <cstddef>
#include <decimate/seam/seam_energy.h>
#include <stdexcept>
#include <vector>

namespace decimate::seam {

/// Reject a texture the exact integration cannot address (extent below 2x2 or a wrong value count).
static void requireTexture(int width, int height, std::size_t valueCount) {
    if (width < 2 || height < 2) {
        throw std::invalid_argument("seamEnergy: texture must be at least 2x2");
    }
    if (valueCount != static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 3) {
        throw std::invalid_argument("seamEnergy: texel count must be width * height * 3");
    }
}


std::array<double, 3> seamEnergy(const TexturedMesh &tm, int width, int height,
                                 std::span<const double> texels) {
    requireTexture(width, height, texels.size());

    // No seams (or only zero-length ones): the energy is exactly zero, and the normalization by seam length would divide by 0.
    const SeamTopology topo = findSeams(tm);
    if (!detail::hasSeamLength(topo)) {
        return { 0.0, 0.0, 0.0 };
    }

    const std::vector<detail::SeamPiece> pieces =
        detail::seamPieces(tm, topo, { width, height }, detail::PieceKind::Value);
    return detail::integrateEnergy(pieces, detail::totalSeamLength(topo), texels);
}


std::array<double, 3> seamEnergy(const TexturedMesh &tm, const bake::NormalMap &map) {
    if (map.pixels.size() !=
        static_cast<std::size_t>(map.width) * static_cast<std::size_t>(map.height)) {
        throw std::invalid_argument("seamEnergy: NormalMap pixel count does not match its size");
    }

    // Bytes to [0,1], the unit D_total is reported in.
    std::vector<double> texels;
    texels.reserve(map.pixels.size() * 3);
    for (const bake::Rgb &p : map.pixels) {
        texels.push_back(static_cast<double>(p.r) / 255.0);
        texels.push_back(static_cast<double>(p.g) / 255.0);
        texels.push_back(static_cast<double>(p.b) / 255.0);
    }
    return seamEnergy(tm, map.width, map.height, texels);
}

} // namespace decimate::seam
