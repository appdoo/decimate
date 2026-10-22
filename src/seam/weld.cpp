// This file is part of the decimate library; see decimate.h for version/license details.
/// @file weld.cpp
/// @brief seam::weld -- per-corner UVs -> TexturedMesh with UV vertex identity.

#include <bit>
#include <cstddef>
#include <cstdint>
#include <decimate/seam/seam_topology.h>
#include <map>
#include <stdexcept>
#include <tuple>

namespace decimate::seam {

namespace {

/// (surface vertex, UV bit pattern): the identity of a UV vertex.
using UvKey = std::tuple<int, std::uint32_t, std::uint32_t>;

} // namespace


/// Key of the UV vertex a corner belongs to; bit patterns, so only bit-identical UVs merge.
static UvKey makeKey(int surfaceVertex, math::Vec2f uv) {
    return { surfaceVertex, std::bit_cast<std::uint32_t>(uv.x), std::bit_cast<std::uint32_t>(uv.y) };
}


TexturedMesh weld(const MeshData &mesh, const std::vector<math::Vec2f> &cornerUV) {
    if (mesh.indices.size() % 3 != 0 || cornerUV.size() != mesh.indices.size()) {
        throw std::invalid_argument("weld: cornerUV must hold 3 UVs per triangle");
    }

    // Surface indices are kept verbatim; only the UV side is re-indexed.
    TexturedMesh tm;
    tm.positions = mesh.positions;
    tm.faceV = mesh.indices;
    tm.faceU.resize(mesh.indices.size());

    // First corner with a given key creates the UV vertex; later ones reuse it.
    std::map<UvKey, int> uvIndex;
    for (std::size_t i = 0; i < mesh.indices.size(); ++i) {
        const auto [it, fresh] = uvIndex.try_emplace(makeKey(mesh.indices[i], cornerUV[i]),
                                                     static_cast<int>(tm.uvs.size()));
        if (fresh) {
            tm.uvs.push_back(cornerUV[i]);
        }
        tm.faceU[i] = it->second;
    }
    return tm;
}

} // namespace decimate::seam
