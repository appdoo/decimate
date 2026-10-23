// This file is part of the decimate library; see decimate.h for version/license details.
/// @file seam_metric.h
/// @brief Test helper: the angle between the decoded normals on the two sides of every chart seam.
///
/// This is what a renderer shows at a seam: each side is read the way an engine reads a texture
/// (bilinear, texel x centred at uv = (x + 0.5) / N, i.e. uv * N - 0.5 in texel units), decoded with
/// `c / 127 - 1`, and, in tangent space, turned into a world normal by that side's own frame. A bake
/// whose texel convention differs from this reader shows up as a seam angle, not as a missing texel.
#pragma once

#include <decimate/bake/normal_map.h>
#include <decimate/bake/normals.h>
#include <decimate/pipeline.h>
#include <decimate/seam/seam_topology.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace seammetric {

struct Stats {
    double mean = 0.0; ///< seam-length weighted mean angle, degrees
    double max = 0.0;  ///< worst sample, degrees
};


inline decimate::math::Vec3f unpack(decimate::bake::Rgb c) {
    return { static_cast<float>(c.r) / 127.0f - 1.0f, static_cast<float>(c.g) / 127.0f - 1.0f,
             static_cast<float>(c.b) / 127.0f - 1.0f };
}


inline decimate::math::Vec3f unit(decimate::math::Vec3f v) {
    const float l = std::sqrt(v.dot(v));
    return l > 0.0f ? v * (1.0f / l) : v;
}


/// Bilinear fetch in the reader convention, clamped to the border.
inline decimate::math::Vec3f fetch(const decimate::bake::NormalMap &m, decimate::math::Vec2f uv) {
    const float x = uv.x * static_cast<float>(m.width) - 0.5f;
    const float y = uv.y * static_cast<float>(m.height) - 0.5f;
    const int   x0 = static_cast<int>(std::floor(x));
    const int   y0 = static_cast<int>(std::floor(y));
    const float fx = x - static_cast<float>(x0);
    const float fy = y - static_cast<float>(y0);
    auto        at = [&](int i, int j) {
        i = std::clamp(i, 0, m.width - 1);
        j = std::clamp(j, 0, m.height - 1);
        return unpack(m.pixels[static_cast<std::size_t>(j) * m.width + i]);
    };
    return (at(x0, y0) * (1 - fx) + at(x0 + 1, y0) * fx) * (1 - fy) +
           (at(x0, y0 + 1) * (1 - fx) + at(x0 + 1, y0 + 1) * fx) * fy;
}


/// World normal of the texel read at barycentrics (w0, w1, w2) of face `f`.
inline decimate::math::Vec3f decode(const decimate::Atlas &a, const std::vector<decimate::math::Vec3f> &vn,
                                    const decimate::bake::NormalMap &m, bool tangent, int f, float w0, float w1,
                                    float w2) {
    using decimate::math::Vec2f;
    using decimate::math::Vec3f;
    const Vec3f s = fetch(m, a.uv[3 * f] * w0 + a.uv[3 * f + 1] * w1 + a.uv[3 * f + 2] * w2);
    if (!tangent) {
        return unit(s);
    }

    // The engine's frame: T and N interpolated, B = w * cross(N, T).
    const decimate::math::Vec4f *c = &a.tangents[3 * f];
    const int                   *ix = &a.mesh.indices[3 * f];
    const Vec3f                  t = Vec3f { c[0].x, c[0].y, c[0].z } * w0 + Vec3f { c[1].x, c[1].y, c[1].z } * w1 +
                                     Vec3f { c[2].x, c[2].y, c[2].z } * w2;
    const Vec3f                  n = vn[ix[0]] * w0 + vn[ix[1]] * w1 + vn[ix[2]] * w2;
    return unit(t * s.x + n.cross(t) * (c[0].w * s.y) + n * s.z);
}


/// 16 samples per seam edge; the two sides run in opposite directions along the shared edge.
inline Stats seamStats(const decimate::Atlas &a, const std::vector<decimate::math::Vec3f> &vn,
                       const decimate::bake::NormalMap &m, bool tangent, const decimate::seam::SeamTopology &st) {
    double sum = 0.0;
    double len = 0.0;
    double worst = 0.0;
    for (const auto &sp : st.seams) {
        for (int k = 0; k < 16; k++) {
            const float t = (static_cast<float>(k) + 0.5f) / 16.0f;
            float       wa[3] = { 0, 0, 0 };
            float       wb[3] = { 0, 0, 0 };
            wa[sp.a.corner] = 1 - t;
            wa[(sp.a.corner + 1) % 3] = t;
            wb[sp.b.corner] = t;
            wb[(sp.b.corner + 1) % 3] = 1 - t;
            const auto   na = decode(a, vn, m, tangent, sp.a.face, wa[0], wa[1], wa[2]);
            const auto   nb = decode(a, vn, m, tangent, sp.b.face, wb[0], wb[1], wb[2]);
            const double deg = std::acos(std::clamp(static_cast<double>(na.dot(nb)), -1.0, 1.0)) * 180.0 /
                               3.14159265358979323846;
            sum += deg * sp.length3D / 16;
            len += sp.length3D / 16;
            worst = std::max(worst, deg);
        }
    }
    return { sum / len, worst };
}

} // namespace seammetric
