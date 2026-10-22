// This file is part of the decimate library; see decimate.h for version/license details.
/// @file erase.cpp
/// @brief seam::erase -- seam erasure on an object-space normal map (paper section 4).
///
/// Eigen is private to this file: the factorization is the only place that needs it.

#include "erase_domain.h"
#include "erase_system.h"
#include "seam_pieces.h"

#include <Eigen/Sparse>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <decimate/seam/erase.h>
#include <decimate/seam/seam_energy.h>
#include <decimate/seam/seam_topology.h>
#include <stdexcept>
#include <utility>
#include <vector>

namespace decimate::seam {

namespace {

constexpr std::size_t kChannels = 3;

/// Ridge added to K's diagonal only when K is singular to working precision (the 1e10 seam
/// penalty can make it so, e.g. on the straightened bijective bunny atlas).
constexpr double kRidge = 1e-9;

using SparseMat = Eigen::SparseMatrix<double>;
using DenseMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;


/// One factorization of K, reused for every right-hand side (LDLT; LU if LDLT fails; a ridged
/// LDLT if both fail).
class Solver {
public:
    explicit Solver(const SparseMat &k)
        : ldlt_(k) {
        if (ldlt_.info() == Eigen::Success) {
            return;
        }

        // LU tolerates the near-singular K that LDLT rejects.
        lu_.compute(k);
        if (lu_.info() == Eigen::Success) {
            useLu_ = true;
            return;
        }

        // Last resort: a tiny ridge makes the factorization exist at the cost of a negligible bias.
        SparseMat ridge(k.rows(), k.cols());
        ridge.setIdentity();
        ldlt_.compute(k + kRidge * ridge);
        if (ldlt_.info() != Eigen::Success) {
            throw std::runtime_error("erase: linear system could not be factorized");
        }
    }

    [[nodiscard]] DenseMat solve(const DenseMat &b) const {
        return useLu_ ? DenseMat(lu_.solve(b)) : DenseMat(ldlt_.solve(b));
    }

private:
    Eigen::SimplicialLDLT<SparseMat> ldlt_;
    Eigen::SparseLU<SparseMat>       lu_;
    bool                             useLu_ = false;
};

} // namespace


/// Reject maps the solver cannot address (smaller than 2x2 or a pixel buffer of the wrong size).
static void requireMap(const bake::NormalMap &map) {
    const auto count = static_cast<std::size_t>(std::max(map.width, 0)) *
                       static_cast<std::size_t>(std::max(map.height, 0));
    if (map.width < 2 || map.height < 2 || map.pixels.size() != count) {
        throw std::invalid_argument("erase: NormalMap must be at least 2x2 with matching pixels");
    }
}


/// The map as the solver's native texels: bytes scaled to [0,1], channel-interleaved.
static std::vector<double> toTexels(const bake::NormalMap &map) {
    std::vector<double> t;
    t.reserve(map.pixels.size() * kChannels);
    for (const bake::Rgb &p : map.pixels) {
        t.push_back(static_cast<double>(p.r) / 255.0);
        t.push_back(static_cast<double>(p.g) / 255.0);
        t.push_back(static_cast<double>(p.b) / 255.0);
    }
    return t;
}


/// Back to 8 bits: clamp to [0,1], scale, round.
static std::uint8_t quantize(double v) {
    return static_cast<std::uint8_t>(std::lround(std::clamp(v, 0.0, 1.0) * 255.0));
}


/// The assembled triplets as a compressed sparse matrix (duplicates are summed).
static SparseMat toSparse(const detail::EraseSystem &sys) {
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(sys.entries.size());
    for (const detail::SparseEntry &e : sys.entries) {
        triplets.emplace_back(e.row, e.col, e.value);
    }

    SparseMat k(sys.n, sys.n);
    k.setFromTriplets(triplets.begin(), triplets.end());
    return k;
}


/// Solve the normal equations for all three channels with one factorization.
///
/// The seam constraint enters as the penalty w_seam * D (paper section 4.1), exactly as
/// the reference does.  The penalty solution is deliberately NOT sharpened into an exact
/// constraint: the seam operator has near-null directions, along which the exact
/// constrained minimizer moves texels by tens of 8-bit levels for a negligible change of
/// D (measured: up to 26 levels on the golden fixtures); the finite penalty regularizes
/// them.
static DenseMat solveSystem(const detail::EraseSystem &sys) {
    DenseMat b(sys.n, static_cast<int>(kChannels));
    for (int u = 0; u < sys.n; ++u) {
        for (std::size_t ch = 0; ch < kChannels; ++ch) {
            b(u, static_cast<int>(ch)) = sys.rhs[kChannels * static_cast<std::size_t>(u) + ch];
        }
    }

    const SparseMat k = toSparse(sys);
    const Solver    solver(k);
    DenseMat        x = solver.solve(b);

    // A NaN would be quantized to an arbitrary byte and pass for texture.
    if (!x.allFinite()) {
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


/// Write solution values into the double texels and the quantized map.
static void applySolution(const detail::EraseSystem &sys, const DenseMat &x, EraseResult &out) {
    for (int u = 0; u < sys.n; ++u) {
        const auto                          texel = static_cast<std::size_t>(sys.texelOf[static_cast<std::size_t>(u)]);
        bake::Rgb                          &px = out.map.pixels[texel];
        std::array<std::uint8_t, kChannels> bytes {};
        for (std::size_t ch = 0; ch < kChannels; ++ch) {
            const double v = x(u, static_cast<int>(ch));
            out.texels[kChannels * texel + ch] = v;
            bytes[ch] = quantize(v);
        }
        px = { bytes[0], bytes[1], bytes[2] };
    }
}


EraseResult erase(const TexturedMesh &tm, const bake::NormalMap &map, NormalSpace space,
                  const EraseOptions &options) {
    // Tangent-space seams are not well-defined, so reject before doing any work.
    if (space == NormalSpace::Tangent) {
        throw std::invalid_argument(
            "erase: tangent-space seams are not well-defined (Liu et al. 2017, footnote 3); "
            "erase an object-space normal map");
    }
    requireMap(map);

    const detail::TexelGrid grid { map.width, map.height };
    std::vector<double>     texels = toTexels(map);
    const SeamTopology      topo = findSeams(tm);

    // No seam pairs (or only zero-length ones): the texture is already consistent, hand it back unchanged.
    if (!detail::hasSeamLength(topo)) {
        return passthrough(map, std::move(texels));
    }

    // The exact piecewise-polynomial seam mismatch and the texels the solve may touch.
    const auto                valuePieces = detail::seamPieces(tm, topo, grid, detail::PieceKind::Value);
    const auto                gradPieces = detail::seamPieces(tm, topo, grid, detail::PieceKind::Gradient);
    const double              total = detail::totalSeamLength(topo);
    const detail::EraseDomain domain = detail::buildEraseDomain(tm, topo, grid);

    // Global mode trades fidelity to the input for a correction that spreads further.
    detail::EraseWeights weights;
    weights.change = options.global ? 1e2 : 1e4;
    const detail::EraseSystem sys =
        detail::buildEraseSystem({ domain, grid, texels, valuePieces, gradPieces, total, weights });

    EraseResult result = passthrough(map, texels);
    result.dBefore = detail::integrateEnergy(valuePieces, total, texels);
    result.unknownTexels = sys.n - sys.pinned;
    applySolution(sys, solveSystem(sys), result);
    result.dAfter = detail::integrateEnergy(valuePieces, total, result.texels);
    return result;
}

} // namespace decimate::seam
