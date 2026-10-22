// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file erase.h
/// @brief Seam erasure: a texture post-process that removes the residual seam (Liu et
///        al. 2017, section 4).
///
/// Gutter filling (dilate, seam-fill) only narrows the seam: along a seam edge the two
/// UV half-edges still sample different texels, so bilinear reconstruction disagrees.
/// `erase()` solves for new texel values that make the two reconstructions agree
/// exactly along every seam edge, while staying as close as possible to the input and
/// keeping the texture smooth.  It changes neither the mesh nor the UVs.
///
/// ## Contract
///
/// - **Object-space normal maps only** (paper footnote 3): a tangent-space seam is not
///   well-defined, so `NormalSpace::Tangent` throws `std::invalid_argument`.
/// - **GPU sampling convention:** texel coordinate = `uv * N - 0.5`, clamp-to-edge
///   (see seam_energy.h).  The erasure is exact in what a GL_LINEAR sampler fetches.
/// - **Units:** the solve runs in double on texels in [0,1] (byte / 255).
/// - **8-bit output:** the double solution is clamped to [0,1], scaled by 255 and
///   rounded.  The result is NOT renormalized -- the shader normalizes the fetched
///   normal.
/// - **Unknowns:** only texels that belong to the UV mesh (centre inside a UV triangle)
///   or to a bilinear cell crossed by the UV silhouette are optimized.  Every other
///   texel is copied bit-identically.
/// - A mesh without seam pairs is a no-op.
///
/// ## Energy (paper section 4.1)
///
/// `w_seam * D + w_change * E_change + w_C1 * E_C1 + w_grad * E_grad` with
/// `w_seam = 1e10`, `w_change = 1e4` (`1e2` in global mode), `w_C1 = 1e2`,
/// `w_grad = 1`.  One sparse factorization serves all three channels.

#include <array>
#include <decimate/bake/normal_map.h>
#include <decimate/types/textured_mesh.h>
#include <vector>

namespace decimate::seam {

/// @brief Which frame the normals of the map are expressed in.
///
/// `NormalMap` carries no space tag, so the caller states it.
enum class NormalSpace { Object,
                         Tangent };


/// @brief Solver options.
struct EraseOptions {
    /// Global mode: weaker anchor to the input (`w_change = 1e2`) so the correction may
    /// spread further from the seams.
    bool global = false;
};


/// @brief Outcome of `erase()`.
struct EraseResult {
    bake::NormalMap     map;                 ///< Quantized result, same size as the input.
    std::vector<double> texels;              ///< `w*h*3` solution values, pre-quantization, row-major,
                                             ///< channel-interleaved; non-unknown texels are byte/255.
    std::array<double, 3> dBefore {};        ///< D_total of the input, per channel.
    std::array<double, 3> dAfter {};         ///< D_total of `texels`, per channel.
    int                   unknownTexels = 0; ///< Texels that took part in the solve.
};


/// @brief Erase the texture seams of @p map.
///
/// @param tm       Textured mesh whose UVs address @p map.
/// @param map      Baked object-space normal map (>= 2x2).
/// @param space    Frame of @p map; must be `NormalSpace::Object`.
/// @param options  Solver options.
/// @returns the erased map plus before/after seam energy.
/// @throws std::invalid_argument if @p space is `Tangent` or @p map is malformed;
///         @p map is never modified.
/// @throws std::runtime_error if the linear system cannot be factorized.
[[nodiscard]] EraseResult erase(const TexturedMesh &tm, const bake::NormalMap &map,
                                NormalSpace space, const EraseOptions &options = {});

} // namespace decimate::seam
