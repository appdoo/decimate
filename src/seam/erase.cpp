// This file is part of the decimate library; see decimate.h for version/license details.
/// @file erase.cpp
/// @brief seam::erase -- seam erasure on an object-space or tangent-space normal map (paper section 4).

#include "block_ldlt.h"
#include "erase_domain.h"
#include "erase_system.h"
#include "seam_pieces.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <decimate/seam/erase.h>
#include <decimate/seam/seam_energy.h>
#include <decimate/seam/seam_topology.h>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

namespace decimate::seam {

namespace {

constexpr std::size_t kChannels = 3;

/// Ridge added to K's diagonal only when the factorization refuses K (the 1e10 seam penalty can make it
/// singular to working precision, e.g. on the straightened bijective bunny atlas).
constexpr double kRidge = 1e-9;

} // namespace


/// Reject maps the solver cannot address (smaller than 2x2 or a pixel buffer of the wrong size).
static void requireMap(const bake::NormalMap &map) {
    const auto count = static_cast<std::size_t>(std::max(map.width, 0)) *
                       static_cast<std::size_t>(std::max(map.height, 0));
    if (map.width < 2 || map.height < 2 || map.pixels.size() != count) {
        throw std::invalid_argument("erase: NormalMap must be at least 2x2 with matching pixels");
    }
}


/// Byte <-> texel value: `v = byte / scale - offset`.  Object space keeps the bytes in [0,1]; tangent
/// space uses the engine's signed decode `byte / 127 - 1`, because the frames act on the signed vector.
struct Codec {
    double offset = 0.0;
    double scale = 255.0;
};

constexpr Codec kObjectCodec { 0.0, 255.0 };
constexpr Codec kTangentCodec { 1.0, 127.0 };


/// The map as the solver's native texels, channel-interleaved.
static std::vector<double> toTexels(const bake::NormalMap &map, Codec codec) {
    std::vector<double> t;
    t.reserve(map.pixels.size() * kChannels);
    for (const bake::Rgb &p : map.pixels) {
        t.push_back(static_cast<double>(p.r) / codec.scale - codec.offset);
        t.push_back(static_cast<double>(p.g) / codec.scale - codec.offset);
        t.push_back(static_cast<double>(p.b) / codec.scale - codec.offset);
    }
    return t;
}


/// Back to 8 bits: scale, clamp to [0,255], round to nearest.
static std::uint8_t quantize(double v, Codec codec) {
    return static_cast<std::uint8_t>(std::lround(std::clamp((v + codec.offset) * codec.scale, 0.0, 255.0)));
}


/// Solve the normal equations for all three channels with one factorization (`3*u + channel`).
///
/// The seam constraint enters as the penalty w_seam * D (paper section 4.1), exactly as
/// the reference does.  The penalty solution is deliberately NOT sharpened into an exact
/// constraint: the seam operator has near-null directions, along which the exact
/// constrained minimizer moves texels by tens of 8-bit levels for a negligible change of
/// D (measured: up to 26 levels on the golden fixtures); the finite penalty regularizes
/// them.  If the factorization refuses K, one retry with a tiny ridge on the diagonal; then it throws.
template <int B>
static std::vector<double> solveSystem(const detail::EraseSystem<B> &sys) {
    detail::BlockLdlt<B> ldlt(sys.blocks);
    if (!ldlt.factorize(sys.blocks)) {
        detail::BlockMatrix<B> ridged = sys.blocks;
        for (std::size_t u = 0; u < static_cast<std::size_t>(sys.n); ++u) {
            for (int k = 0; k < B; ++k) {
                ridged.diag[B * B * u + static_cast<std::size_t>(k * (B + 1))] += kRidge;
            }
        }
        if (!ldlt.factorize(ridged)) {
            throw std::runtime_error("erase: linear system could not be factorized");
        }
    }

    std::vector<double> x(sys.rhs.size());
    ldlt.solve(sys.rhs, x);

    // A NaN would be quantized to an arbitrary byte and pass for texture.
    if (!std::ranges::all_of(x, [](double v) { return std::isfinite(v); })) {
        throw std::runtime_error("erase: linear solve produced non-finite values");
    }
    return x;
}


/// A result that leaves @p map untouched (no seams: nothing to erase).
static EraseResult passthrough(const bake::NormalMap &map, std::vector<double> texels) {
    EraseResult r;
    r.map = map;
    r.texels = std::move(texels);
    return r;
}


/// Write solution values (times @p scale) into the double texels and the quantized map.
template <int B>
static void applySolution(const detail::EraseSystem<B> &sys, const std::vector<double> &x, double scale, Codec codec,
                          EraseResult &out) {
    for (int u = 0; u < sys.n; ++u) {
        const auto                          texel = static_cast<std::size_t>(sys.texelOf[static_cast<std::size_t>(u)]);
        bake::Rgb                          &px = out.map.pixels[texel];
        std::array<std::uint8_t, kChannels> bytes {};
        for (std::size_t ch = 0; ch < kChannels; ++ch) {
            const double v = scale * x[kChannels * static_cast<std::size_t>(u) + ch];
            out.texels[kChannels * texel + ch] = v;
            bytes[ch] = quantize(v, codec);
        }
        px = { bytes[0], bytes[1], bytes[2] };
    }
}


/// The erasure of one space: B = 1 object space (@p frames null), B = 3 tangent space (@p frames set).
template <int B>
static EraseResult eraseSpace(const TexturedMesh &tm, const bake::NormalMap &map, const detail::FrameField *frames,
                              const EraseOptions &options) {
    requireMap(map);

    const Codec             codec = B == 3 ? kTangentCodec : kObjectCodec;
    const detail::TexelGrid grid { map.width, map.height };
    std::vector<double>     texels = toTexels(map, codec);
    const SeamTopology      topo = findSeams(tm);

    // No seam pairs (or only zero-length ones): the texture is already consistent, hand it back unchanged.
    if (!detail::hasSeamLength(topo)) {
        return passthrough(map, std::move(texels));
    }

    // The exact piecewise-polynomial seam mismatch and the texels the solve may touch.
    const auto                valuePieces = detail::seamPieces(tm, topo, grid, detail::PieceKind::Value, frames);
    const auto                gradPieces = detail::seamPieces(tm, topo, grid, detail::PieceKind::Gradient);
    const double              total = detail::totalSeamLength(topo);
    const detail::EraseDomain domain = detail::buildEraseDomain(tm, topo, grid);

    // Global mode trades fidelity to the input for a correction that spreads further.
    detail::EraseWeights weights;
    weights.change = options.global ? 1e2 : 1e4;

    const auto sys = detail::buildEraseSystem<B>({ domain, grid, texels, valuePieces, gradPieces, total, weights });

    EraseResult result = passthrough(map, texels);
    const auto  energy = [&](std::span<const double> t) { return detail::integrateEnergy(valuePieces, total, t); };
    result.dBefore = energy(texels);
    result.unknownTexels = sys.n - sys.pinned;

    // r = 0 fixes only the direction: the framed solve may leave [-1, 1], and a global shrink keeps the range.
    const std::vector<double> x = solveSystem(sys);
    if (B == 3) {
        const double peak = std::ranges::max(x | std::views::transform([](double v) { return std::abs(v); }));
        result.rescale = std::min(1.0, 1.0 / peak);
    }
    applySolution(sys, x, result.rescale, codec, result);
    result.dAfter = energy(result.texels);
    return result;
}


EraseResult erase(const TexturedMesh &tm, const bake::NormalMap &map, NormalSpace space,
                  const EraseOptions &options) {
    // Without frames a tangent-space seam is not well-defined, so reject before doing any work.
    if (space == NormalSpace::Tangent) {
        throw std::invalid_argument("erase: a tangent-space map needs its frames (use the TangentFrames overload)");
    }
    return eraseSpace<1>(tm, map, nullptr, options);
}


EraseResult erase(const TexturedMesh &tm, const bake::NormalMap &map, const TangentFrames &frames,
                  const EraseOptions &options) {
    if (frames.tangents.size() != tm.faceV.size()) {
        throw std::invalid_argument("erase: TangentFrames needs one tangent per face corner");
    }
    for (const int v : tm.faceV) {
        if (static_cast<std::size_t>(v) >= frames.normals.size()) {
            throw std::invalid_argument("erase: TangentFrames needs a normal per surface vertex");
        }
    }
    const detail::FrameField field { frames.tangents, frames.normals };
    return eraseSpace<3>(tm, map, &field, options);
}

} // namespace decimate::seam
