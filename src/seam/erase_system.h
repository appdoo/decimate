// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file erase_system.h
/// @brief Internal: the sparse normal equations of the erase energy (paper section 4.1).
///
/// Plain data only (no Eigen): `erase.cpp` hands `entries` to the factorization.  The
/// system is `K x = rhs` per channel, K symmetric positive definite over the unknown
/// texels, with
///   K = w_seam Q_D + w_C1 Q_G + w_change/|C| diag(trusted) + Q_dirichlet.

#include "erase_domain.h"
#include "seam_pieces.h"

#include <span>
#include <vector>

namespace decimate::seam::detail {

/// @brief One (row, col, value) term; duplicates are summed by the consumer.
struct SparseEntry {
    int    row = 0;
    int    col = 0;
    double value = 0.0;
};


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
    std::span<const SeamPiece> valuePieces; ///< PieceKind::Value.
    std::span<const SeamPiece> gradPieces;  ///< PieceKind::Gradient.
    double                     totalLength = 0.0;
    EraseWeights               weights;
};


/// @brief Assembled normal equations over the unknown texels.
struct EraseSystem {
    int                      n = 0;      ///< Number of unknowns.
    std::vector<int>         unknownOf;  ///< texel -> unknown index, or -1.
    std::vector<int>         texelOf;    ///< unknown index -> texel.
    std::vector<SparseEntry> entries;    ///< Lower+upper triangle, duplicates allowed.
    std::vector<double>      rhs;        ///< n*3, row-major (`3*u + channel`).
    int                      pinned = 0; ///< Unknowns with no energy at all, fixed to input.
};


/// @brief Assemble the system for @p in.
///
/// @pre `in.totalLength > 0` and every piece texel lies inside `in.domain.unknown`.
[[nodiscard]] EraseSystem buildEraseSystem(const EraseInput &in);

} // namespace decimate::seam::detail
