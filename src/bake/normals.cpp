// This file is part of the decimate library; see decimate.h for version/license details.
/// @file normals.cpp
/// @brief Per-vertex normal computation - port of legacy createNormals().
///
/// Port of `createNormals` (legacy2003/main.cpp:87-108): same arithmetic, float by float, on
/// MeshData/Vec3f instead of raw arrays. The Area weighting is an addition, not 2003 behavior.

#include <decimate/bake/normals.h>

#include <cmath>

namespace decimate::bake {

std::vector<math::Vec3f> computeVertexNormals(const decimate::MeshData &data,
                                              Weighting                 w) {
    const int nv = static_cast<int>(data.positions.size());
    const int nf = static_cast<int>(data.indices.size()) / 3;

    std::vector<math::Vec3f> normals(nv); // zero-initialized by Vec3 ctor

    const auto &p = data.positions;

    // Accumulate every face normal into its three corner vertices.
    for (int fi = 0; fi < nf; fi++) {
        int va = data.indices[3 * fi + 0];
        int vb = data.indices[3 * fi + 1];
        int vc = data.indices[3 * fi + 2];

        math::Vec3f edgeAb = p[vb] - p[va];
        math::Vec3f edgeAc = p[vc] - p[va];

        // Cross product: magnitude = 2*area, which is the weight when Area skips the normalize below.
        math::Vec3f nrm = edgeAb.cross(edgeAc);

        // EqualUnit: normalize BEFORE accumulation, as 2003 does (legacy2003/main.cpp:97).
        // Area: skip it; the 2*area magnitude survives into the sum.
        if (w == Weighting::EqualUnit) {
            nrm.normalize();
        }

        normals[va] += nrm;
        normals[vb] += nrm;
        normals[vc] += nrm;
    }

    // Unit length per vertex; an isolated vertex stays zero (see the header contract).
    for (auto &n : normals) {
        n.normalize();
    }

    return normals;
}

} // namespace decimate::bake
