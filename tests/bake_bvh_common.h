// This file is part of the decimate library; see decimate.h for version/license details.
/// @file bake_bvh_common.h
/// @brief Shared fixtures for the Spec 22-01 tests: the brute-force oracle (test-only), the
///        reduce -> segment -> parametrize pipeline, and the random-ray generator.
#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normals.h>
#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/raytracing/bvh.h>
#include <decimate/raytracing/raytracing.h>
#include <decimate/reduce/reducer.h>
#include <decimate/types/mesh_data.h>

#include "bake_detail.h"

namespace bvhtest {

using Hit = rtrace::TriangleBvh::Hit;
using Vec3 = decimate::math::Vec3f;

/// The 2003 brute-force nearest-hit loop, verbatim in behaviour: every face ascending,
/// intersect_triangle on mutable copies of the three vertices, keep a hit iff |t| < |mint|,
/// starting from mint = maxAbsT (legacy2003/split/texturesample.cpp:84,97).
inline Hit brute(const decimate::MeshData &m, Vec3 orig, Vec3 dir, float maxAbsT) {
    Hit       best;
    float     mint = maxAbsT;
    const int nf = static_cast<int>(m.indices.size()) / 3;
    for (int i = 0; i < nf; i++) {
        Vec3  v0 = m.positions[m.indices[3 * i + 0]];
        Vec3  v1 = m.positions[m.indices[3 * i + 1]];
        Vec3  v2 = m.positions[m.indices[3 * i + 2]];
        float t = 0.0f;
        float u = 0.0f;
        float v = 0.0f;
        if (rtrace::intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &v) && std::fabs(t) < std::fabs(mint)) {
            mint = t;
            best = { i, t, u, v };
        }
    }
    return best;
}


/// Test-only Tracer: the brute-force loop behind the production seam.
class BruteTracer final : public decimate::bake::detail::Tracer {
public:
    explicit BruteTracer(const decimate::MeshData &high)
        : high_(high) {
    }
    [[nodiscard]] Hit nearest(Vec3 orig, Vec3 dir, float maxAbsT) const override {
        return brute(high_, orig, dir, maxAbsT);
    }

private:
    const decimate::MeshData &high_;
};


/// The IEEE bit pattern of @p f, so hits compare bit-equal instead of within a tolerance.
inline std::uint32_t bits(float f) {
    return std::bit_cast<std::uint32_t>(f);
}


/// Same face, and t/u/v bit-equal (a miss compares only the face).
inline bool sameHit(const Hit &a, const Hit &b) {
    if (a.face != b.face) {
        return false;
    }
    return a.face < 0 || (bits(a.t) == bits(b.t) && bits(a.u) == bits(b.u) && bits(a.v) == bits(b.v));
}


/// Same size and byte-identical pixels (the bake must be deterministic, not just close).
inline bool samePixels(const decimate::bake::NormalMap &a, const decimate::bake::NormalMap &b) {
    return a.width == b.width && a.height == b.height && a.pixels.size() == b.pixels.size() &&
           std::memcmp(a.pixels.data(), b.pixels.data(), a.pixels.size() * sizeof(decimate::bake::Rgb)) == 0;
}


/// reduce(high, target) -> HalfEdgeMesh -> segment(cap) -> parametrizeCorners(texSize), EqualUnit normals.
struct Pipeline {
    decimate::MeshData                 high;
    decimate::MeshData                 low;
    std::vector<decimate::math::Vec2f> cornerUV;
    std::vector<Vec3>                  lowN;
    std::vector<Vec3>                  highN;
    std::vector<int>                   faceToChart;
    decimate::HalfEdgeMesh             hem;
};


/// Builds the low-poly bake input from the model at @p path; fixture for every real-pipeline bake test.
inline Pipeline buildPipeline(const std::string &path, int reduceTarget, int chartCap, int texSize) {
    Pipeline p;
    p.high = decimate::load(path);
    p.low = decimate::reduce(p.high, reduceTarget);

    // The chart segmenter works on a half-edge mesh of the reduced geometry.
    const int                   nlv = static_cast<int>(p.low.positions.size());
    const int                   nlf = static_cast<int>(p.low.indices.size()) / 3;
    std::vector<decimate::Face> faces(nlf);
    for (int i = 0; i < nlf; i++) {
        faces[i] = { p.low.indices[3 * i + 0], p.low.indices[3 * i + 1], p.low.indices[3 * i + 2] };
    }
    p.hem = decimate::HalfEdgeMesh::create(nlv, nlf, faces.data());

    // Per-corner UVs and equal-weight vertex normals: the inputs bake() consumes.
    decimate::chart::ChartSet charts = decimate::chart::segment(p.hem, p.low, chartCap);
    p.cornerUV = decimate::chart::parametrizeCorners(p.hem, p.low, charts, texSize);
    p.faceToChart.assign(charts.faceToChart.begin(), charts.faceToChart.end());
    p.lowN = decimate::bake::computeVertexNormals(p.low, decimate::bake::Weighting::EqualUnit);
    p.highN = decimate::bake::computeVertexNormals(p.high, decimate::bake::Weighting::EqualUnit);
    return p;
}


/// One query ray: origin and direction.
struct Ray {
    Vec3 o;
    Vec3 d;
};


/// A reproducible ray set (std::mt19937(22)): half the origins uniform in the 1.2x-scaled bbox,
/// half on random vertices; directions uniform on the sphere, 5% axis-aligned (two zero components).
inline std::vector<Ray> randomRays(const decimate::MeshData &m, int count) {
    std::mt19937 rng(22);

    // Bounding box of the mesh, to place origins around it.
    Vec3 lo = m.positions[0];
    Vec3 hi = m.positions[0];
    for (const auto &p : m.positions) {
        for (int k = 0; k < 3; k++) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    }
    const Vec3                            mid { (lo.x + hi.x) * 0.5f, (lo.y + hi.y) * 0.5f, (lo.z + hi.z) * 0.5f };
    std::uniform_real_distribution<float> uni(-1.0f, 1.0f);
    std::uniform_real_distribution<float> uni01(0.0f, 1.0f);
    std::normal_distribution<float>       gauss;
    std::uniform_int_distribution<int>    pickV(0, static_cast<int>(m.positions.size()) - 1);
    std::uniform_int_distribution<int>    pick3(0, 2);


    std::vector<Ray> rays;
    rays.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; i++) {
        Ray r;
        if (i % 2 == 0) {
            for (int k = 0; k < 3; k++) {
                r.o[k] = mid[k] + uni(rng) * 0.6f * (hi[k] - lo[k]); // 1.2x bbox
            }
        } else {
            r.o = m.positions[pickV(rng)];
        }

        // Axis-aligned directions stress the parallel-ray paths of the triangle test.
        if (uni01(rng) < 0.05f) {
            r.d = Vec3 { 0.0f, 0.0f, 0.0f };
            r.d[pick3(rng)] = uni01(rng) < 0.5f ? -1.0f : 1.0f;
        } else {
            r.d = Vec3 { gauss(rng), gauss(rng), gauss(rng) };
            r.d.normalize();
        }
        rays.push_back(r);
    }
    return rays;
}

} // namespace bvhtest
