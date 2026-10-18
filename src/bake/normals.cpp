// This file is part of the decimate library; see decimate.h for version/license details.
/// @file normals.cpp
/// @brief Per-vertex normal computation — port of legacy createNormals().
///
/// Behavior-preserving port of `legacy/main.cpp:87-108`.
/// The ONLY observable difference is the container type (MeshData/Vec3<float>
/// vs raw arrays); the arithmetic is identical, float by float.

#include <decimate/bake/normals.h>

#include <cmath>

namespace decimate::bake {

std::vector<decimate::math::Vec3<float>> computeVertexNormals(const decimate::MeshData& data,
                                                     Weighting             w) {
    const int nv = static_cast<int>(data.positions.size());
    const int nf = static_cast<int>(data.indices.size()) / 3;

    std::vector<decimate::math::Vec3<float>> normals(nv);  // zero-initialized by Vec3 ctor

    const auto& p = data.positions;

    for (int fi = 0; fi < nf; fi++) {
        int va = data.indices[3 * fi + 0];
        int vb = data.indices[3 * fi + 1];
        int vc = data.indices[3 * fi + 2];

        // Edge vectors (b-a) and (c-a)
        decimate::math::Vec3<float> edgeAb = p[vb] - p[va];
        decimate::math::Vec3<float> edgeAc = p[vc] - p[va];

        // Face normal (cross product; magnitude = 2·area)
        decimate::math::Vec3<float> nrm = edgeAb.cross(edgeAc);

        // EqualUnit: normalize BEFORE accumulation — reproduces 2003 exactly.
        // Area: skip normalize; let 2·area magnitude survive into accumulation.
        if (w == Weighting::EqualUnit) {
            nrm.normalize();  // same as legacy Vec3::norm()
        }

        // Accumulate into the three corner vertices
        normals[va] += nrm;
        normals[vb] += nrm;
        normals[vc] += nrm;
    }

    // Final normalize: bring every vertex normal to unit length
    for (auto& n : normals) {
        n.normalize();
    }

    return normals;
}

}  // namespace decimate::bake
