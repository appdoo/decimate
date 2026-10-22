// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file seam_energy.h
/// @brief The seam mismatch metric D_total (Liu et al. 2017, eq. 3-9).
///
/// For each seam pair the two half-edges sample the texture at different places; the
/// energy is the integral of the squared difference of the two bilinear
/// reconstructions along the edge, normalized by the 3D length of all seam edges.
///
/// ## Sampling convention (the GPU's)
///
/// A texture of N texels maps uv to texel coordinate `uv * N - 0.5`; lookups clamp to
/// the edge texels (GL_LINEAR + GL_CLAMP_TO_EDGE).  Texel centres sit at integer
/// coordinates.  `NormalMap` row `y` corresponds to `v`, as GL texture row 0 does.
///
/// ## Units
///
/// Texel values are in [0,1] (byte / 255), per channel.  D_total is reported in those
/// units, so a quoted value is directly comparable with the reference implementation.
///
/// ## Exactness
///
/// Each seam edge is split at every parameter where either half-edge crosses a texel
/// grid line; on each piece both reconstructions are polynomials of degree <= 2 in the
/// edge parameter, so the squared difference (degree <= 4) is integrated exactly.

#include <array>
#include <decimate/bake/normal_map.h>
#include <decimate/types/textured_mesh.h>
#include <span>

namespace decimate::seam {

/// @brief D_total of an 8-bit normal map, per channel (R, G, B).
///
/// Bit-identical to the `texels` overload fed with `byte / 255.0`.  Returns exact
/// zeros when @p tm has no seam pair.
///
/// @pre `map.width >= 2`, `map.height >= 2`, `map.pixels.size() == width * height`.
/// @throws std::invalid_argument when a precondition fails.
[[nodiscard]] std::array<double, 3> seamEnergy(const TexturedMesh &tm, const bake::NormalMap &map);


/// @brief D_total of a double-precision texture, per channel.
///
/// @param tm      Textured mesh; its UVs address the texture.
/// @param width   Texture width in texels (>= 2).
/// @param height  Texture height in texels (>= 2).
/// @param texels  `width * height * 3` values in [0,1], row-major, channel-interleaved.
/// @throws std::invalid_argument when the size or the extent is invalid.
[[nodiscard]] std::array<double, 3> seamEnergy(const TexturedMesh &tm, int width, int height,
                                               std::span<const double> texels);

} // namespace decimate::seam
