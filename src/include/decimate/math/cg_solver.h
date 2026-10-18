// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file cg_solver.h
/// @brief Jacobi-preconditioned conjugate-gradient solver for sparse symmetric
///        positive-(semi)definite systems — the CG subset of legacy smcr.h.
///
/// Ported from `legacy/math/smcr.h` (`SparseMatrixCR` + `cgiters`, Spec 08-02).
/// Only the CG path is ported; the legacy banded-LU, Cuthill-McKee reordering,
/// and ILU preconditioner are deliberately NOT ported (Beat 2b).
///
/// Used by the LSCM conformal mapper to solve the normal equations AᵀA x = Aᵀb.
/// This is NOT the Spec-05 `SMat` (which is Gauss-Seidel — the wrong method for
/// this system class).
///
/// Boundary contract: pure numeric kernel — no GL/OS, no I/O, no global state.

#include <vector>

namespace decimate::math {

/// @brief Sparse symmetric system A·x = b solved by Jacobi-preconditioned CG.
///
/// Build the matrix by accumulating entries (symmetric — supply the full matrix
/// or the upper/lower consistently as the assembly does), set the right-hand
/// side, then `solve()`.  Entries at the same (row,col) accumulate (`+=`).
class ConjugateGradientSolver {
public:
    /// @param n  Number of unknowns (square n×n system).
    explicit ConjugateGradientSolver(int n);

    /// @brief Accumulate `value` into A(row, col)  (A(row,col) += value).
    void addEntry(int row, int col, double value);

    /// @brief Set b(row) = value.
    void setRhs(int row, double value);

    /// @brief Solve A·x = b by Jacobi-preconditioned conjugate gradient.
    ///
    /// @param tol      Convergence tolerance on the (preconditioned) residual norm.
    /// @param maxIter  Iteration cap; 0 ⇒ default (legacy 15000).
    /// @return The solution vector x (size n).  On a singular/non-converged
    ///         system, returns the best iterate reached at maxIter.
    [[nodiscard]] std::vector<double> solve(double tol = 1e-12, int maxIter = 0) const;

    /// @brief Number of unknowns.
    [[nodiscard]] int size() const noexcept { return n_; }

private:
    int n_;
    // Compressed sparse rows accumulated by addEntry; finalized in solve().
    std::vector<std::vector<std::pair<int, double>>> rows_;  ///< row -> (col,value)
    std::vector<double>                              rhs_;   ///< b
};

}  // namespace decimate::math
