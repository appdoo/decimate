// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file erase_system.h
/// @brief Internal: the sparse normal equations of the erase energy (paper section 4.1).
///
/// Plain data only (no Eigen): `erase.cpp` hands `blocks` to the factorization.  The
/// system is `K x = rhs` per channel, K symmetric positive definite over the unknown
/// texels, with
///   K = w_seam Q_D + w_C1 Q_G + w_change/|C| diag(trusted) + Q_dirichlet.
/// B = 1 (object space): K is n x n, the three channels are three right-hand sides.  B = 3 (tangent
/// space): Q_D couples the channels through the frames, so the system is one 3n x 3n solve with a single
/// right-hand side.  Either way K is a BxB-block matrix (`blocks`) whose pattern is known before any value;
/// every scalar term is added once per channel on the diagonal of its block.

#include "block_ldlt.h"
#include "erase_domain.h"
#include "seam_pieces.h"

#include <span>
#include <vector>

namespace decimate::seam::detail {

/// @brief Energy weights of paper section 4.1.
struct EraseWeights {
    double seam = 1e10;    ///< w_seam: penalty for D = p^T M p.
    double change = 1e4;   ///< w_change (1e2 in global mode).
    double gradient = 1e2; ///< w_C1: seam-normal derivative continuity.
    double smooth = 1.0;   ///< w_grad: masked Dirichlet (inside + outside gradients).
};


/// @brief Everything `buildEraseSystem` needs.
struct EraseInput {
    const EraseDomain         &domain;
    TexelGrid                  grid;
    std::span<const double>    texels;      ///< width*height*3 input values in [0,1].
    std::span<const SeamPiece> valuePieces; ///< PieceKind::Value (with frames for B = 3).
    std::span<const SeamPiece> gradPieces;  ///< PieceKind::Gradient.
    double                     totalLength = 0.0;
    EraseWeights               weights;
};


/// @brief Assembled normal equations over the unknown texels.
template <int B>
struct EraseSystem {
    int                 n = 0;      ///< Number of unknowns.
    std::vector<int>    unknownOf;  ///< texel -> unknown index, or -1.
    std::vector<int>    texelOf;    ///< unknown index -> texel.
    BlockMatrix<B>      blocks;     ///< K as n x n blocks of BxB, lower triangle.
    std::vector<double> rhs;        ///< n*3, row-major (`3*u + channel`).
    int                 pinned = 0; ///< Unknowns with no energy at all, fixed to input.
};


/// @brief Assemble the system for @p in.
///
/// @pre `in.totalLength > 0` and every piece texel lies inside `in.domain.unknown`.
template <int B>
[[nodiscard]] EraseSystem<B> buildEraseSystem(const EraseInput &in);


extern template EraseSystem<1> buildEraseSystem<1>(const EraseInput &);
extern template EraseSystem<3> buildEraseSystem<3>(const EraseInput &);

} // namespace decimate::seam::detail
