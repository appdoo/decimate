// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file bake_detail.h
/// @brief NON-PUBLIC test seam of the normal bake (not installed; include via the src/bake path).
///
/// `bakeWith` is the whole bake with the high-mesh ray query injected. Production (`bake`) hands it
/// the BVH; the equality tests hand it the 2003 brute-force loop and demand byte-identical pixels.

#include <vector>

#include <decimate/bake/normal_bake.h>
#include <decimate/raytracing/bvh.h>

namespace decimate::bake::detail {

/// The single path for high-mesh ray queries in the bake.
class Tracer {
public:
    Tracer() = default;
    Tracer(const Tracer &) = delete;
    Tracer &operator=(const Tracer &) = delete;
    virtual ~Tracer() = default;

    /// Nearest face along the line orig + t*dir by smallest |t| (ties: lowest face index);
    /// `maxAbsT` is the exclusive miss bound. Must be safe to call from several threads.
    [[nodiscard]] virtual rtrace::TriangleBvh::Hit nearest(math::Vec3f orig,
                                                           math::Vec3f dir,
                                                           float       maxAbsT) const = 0;
};

/// Normalize `nrm` and encode it as 8-bit RGB: code = round((c + 1) * 127); decode is code / 127 - 1.
[[nodiscard]] Rgb rangeCompress(math::Vec3f nrm);

/// `bake()` with an injected tracer. Same contract and parameters as `bake()`; `high` supplies the
/// face indices for the normal blend and must be the mesh the tracer was built over.
[[nodiscard]] NormalMap bakeWith(const Tracer &tracer, const decimate::MeshData &low,
                                 const std::vector<math::Vec2f> &lowCornerUV,
                                 const std::vector<math::Vec3f> &lowNormals,
                                 const decimate::MeshData       &high,
                                 const std::vector<math::Vec3f> &highNormals,
                                 int texSize, int dilatePasses,
                                 const std::vector<math::Vec4f> &tangents, int threads);

} // namespace decimate::bake::detail
