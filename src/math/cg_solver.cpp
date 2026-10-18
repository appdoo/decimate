// This file is part of the decimate library; see decimate.h for version/license details.
/// @file cg_solver.cpp
/// @brief Jacobi-preconditioned conjugate-gradient solver implementation.
///
/// Ported from `legacy/math/smcr.h`:
///   - Sparse matvec: derived from `SMCR::vecmuls` (lines 601-631), but
///     adapted for the `rows_` accumulation storage (full matrix, not
///     SMCR's compact-symmetric format).
///   - CG iteration structure: `SMCR::cgiters` (lines 924-992) — symmetric
///     CG with diagonal initial guess and relative-residual convergence.
///   - Jacobi preconditioner: `SMCR::asolve` (lines 997-1003) — z = b/diag
///     applied each iteration (this lifts the implicit scaling from cgiters'
///     initial guess into a proper per-iteration preconditioner, matching the
///     header contract).
///   - Convergence criterion: `SMCR::result` (lines 911-916) — stops when
///     ‖r‖ < tol * ‖b‖ (relative), or when both norms vanish.
///
/// NOT ported: banded LU, Cuthill-McKee reordering, ILU preconditioner
/// (`bcgiter`), biconjugate-gradient (`bcgiter`), Gauss-Seidel (`gsseidel`),
/// file I/O helpers, `printf` diagnostics, `setmatherrflag`/`ismatherr`
/// (platform-specific FP traps).

#include <decimate/math/cg_solver.h>

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace decimate::math {

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

ConjugateGradientSolver::ConjugateGradientSolver(int n)
    : n_(n), rows_(static_cast<std::size_t>(n)), rhs_(static_cast<std::size_t>(n), 0.0) {
    if (n <= 0) {
        throw std::invalid_argument("ConjugateGradientSolver: n must be > 0");
    }
}

// ---------------------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------------------

void ConjugateGradientSolver::addEntry(int row, int col,  // NOLINT(bugprone-easily-swappable-parameters): row/col are distinct matrix coordinates; value is a scalar — not swappable in practice
                                       double value) {
    assert(row >= 0 && row < n_);
    assert(col >= 0 && col < n_);

    // Scan for an existing (col, value) pair and accumulate.
    for (auto& [c, v] : rows_[static_cast<std::size_t>(row)]) {
        if (c == col) {
            v += value;
            return;
        }
    }
    rows_[static_cast<std::size_t>(row)].emplace_back(col, value);
}

void ConjugateGradientSolver::setRhs(int row, double value) {
    assert(row >= 0 && row < n_);
    rhs_[static_cast<std::size_t>(row)] = value;
}

// ---------------------------------------------------------------------------
// Internal helpers (file-scope, not exposed in header)
// ---------------------------------------------------------------------------

namespace {

/// @brief Compute r = A * x using the accumulated row storage.
///
/// Ported from SMCR::vecmuls (smcr.h line 601). The legacy version exploited
/// compact-symmetric storage (diagonal first, upper triangle only). Here rows_
/// contains the full matrix as assembled by addEntry, so a plain row-dot-product
/// is correct and simpler.
void matvec(const std::vector<std::vector<std::pair<int, double>>>& rows,
            const std::vector<double>& x, std::vector<double>& result) {
    const std::size_t n = rows.size();
    result.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (const auto& [col, val] : rows[i]) {
            result[i] += val * x[static_cast<std::size_t>(col)];
        }
    }
}

/// @brief Euclidean norm of a vector.
///
/// Ported from SMCR::l2norm / SMCR::l2norm2 (smcr.h lines 823-838).
/// The legacy version guarded against FP traps via ismatherr(); we instead rely
/// on IEEE 754 (inf/nan propagate transparently — the solver exits on the
/// iteration cap in degenerate cases).
double l2norm(const std::vector<double>& v) {
    double sum = 0.0;
    for (double x : v) {
        sum += x * x;
    }
    return std::sqrt(sum);
}

/// @brief Extract the diagonal of A from row storage.
///
/// Returns diag[i] = A(i,i), or 1.0 if A(i,i) == 0 (safe reciprocal fallback,
/// matching SMCR::asolve line 999: `sa[i] != 0.0 ? b[i]/sa[i] : b[i]`).
std::vector<double> extractDiag(const std::vector<std::vector<std::pair<int, double>>>& rows) {
    const std::size_t n = rows.size();
    std::vector<double> diag(n, 1.0);  // default: 1.0 → no scaling
    for (std::size_t i = 0; i < n; ++i) {
        for (const auto& [col, val] : rows[i]) {
            if (col == static_cast<int>(i)) {
                diag[i] = (val != 0.0) ? val : 1.0;
                break;
            }
        }
    }
    return diag;
}

/// @brief Apply the Jacobi preconditioner: z[i] = r[i] / diag[i].
///
/// Ported from SMCR::asolve (smcr.h lines 997-1003).
void jacobiApply(const std::vector<double>& diag, const std::vector<double>& r,
                 std::vector<double>& z) {
    const std::size_t n = diag.size();
    z.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        z[i] = r[i] / diag[i];
    }
}

/// @brief CG iteration state bundled for a single step.
struct CGState {
    std::vector<double> x;    ///< Current solution estimate.
    std::vector<double> r;    ///< Current residual r = b − A·x.
    std::vector<double> z;    ///< Preconditioned residual z = M⁻¹r.
    std::vector<double> p;    ///< Search direction.
    std::vector<double> q;    ///< A·p (scratch, updated each step).
    double              rho;  ///< Current ρ = rᵀz.
};

/// @brief Initialise CG state from the system `A x = b`.
///
/// Sets x = M⁻¹b (Jacobi-scaled initial guess, matches cgiters line 941),
/// computes the initial residual and preconditioned direction.
CGState initCG(const std::vector<std::vector<std::pair<int, double>>>& rows,
               const std::vector<double>&                              rhs,
               const std::vector<double>&                              diag) {
    const std::size_t n = rhs.size();
    CGState s;
    s.x.resize(n);
    for (std::size_t i = 0; i < n; ++i) { s.x[i] = rhs[i] / diag[i]; }

    // r = b − A·x
    matvec(rows, s.x, s.r);
    s.r.resize(n);
    for (std::size_t i = 0; i < n; ++i) { s.r[i] = rhs[i] - s.r[i]; }

    jacobiApply(diag, s.r, s.z);
    s.p = s.z;
    s.q.resize(n, 0.0);
    s.rho = 0.0;
    for (std::size_t i = 0; i < n; ++i) { s.rho += s.r[i] * s.z[i]; }
    return s;
}

/// @brief Execute one CG step: update x, r, z, p, rho.
///
/// Ported from cgiters inner loop (smcr.h lines 949-986).
/// @returns false if the denominator pᵀq vanishes (singular / ill-conditioned).
bool cgStep(const std::vector<std::vector<std::pair<int, double>>>& rows,
            const std::vector<double>&                              diag,
            CGState&                                                s) {
    matvec(rows, s.p, s.q);

    double pTq = 0.0;
    for (std::size_t i = 0; i < s.p.size(); ++i) { pTq += s.p[i] * s.q[i]; }
    if (std::abs(pTq) < 1e-300) { return false; }  // singular/ill-conditioned

    const double alpha = s.rho / pTq;
    for (std::size_t i = 0; i < s.x.size(); ++i) {
        s.x[i] += alpha * s.p[i];
        s.r[i] -= alpha * s.q[i];
    }

    jacobiApply(diag, s.r, s.z);

    double rhoNew = 0.0;
    for (std::size_t i = 0; i < s.r.size(); ++i) { rhoNew += s.r[i] * s.z[i]; }
    const double beta = rhoNew / s.rho;
    for (std::size_t i = 0; i < s.p.size(); ++i) { s.p[i] = s.z[i] + beta * s.p[i]; }
    s.rho = rhoNew;
    return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// Solve
// ---------------------------------------------------------------------------

/// Jacobi-preconditioned CG.
///
/// Algorithm structure (ported from cgiters, smcr.h lines 924-992, with the
/// Jacobi preconditioner from asolve applied per-iteration):
///
///   1. diag ← extract diagonal of A
///   2. x₀   ← M⁻¹b  (Jacobi-scaled initial guess; matches cgiters line 941)
///   3. r₀   ← b − A·x₀
///   4. z₀   ← M⁻¹r₀  (preconditioned residual)
///   5. p₀   ← z₀;  ρ₀ ← rᵀz
///   6. loop until ‖r‖ < tol·‖b‖ or iter ≥ maxIter:
///        q     ← A·p
///        α     ← ρ / (pᵀq)
///        x     ← x + α·p
///        r     ← r − α·q
///        z     ← M⁻¹r
///        ρ_new ← rᵀz
///        β     ← ρ_new / ρ
///        p     ← z + β·p
///        ρ     ← ρ_new
///
/// Convergence test (SMCR::result, smcr.h line 911):
///   stop when ‖r‖ < tol·‖b‖, or when ‖b‖ ≈ 0 and ‖r‖ ≈ 0.
std::vector<double> ConjugateGradientSolver::solve(double tol, int maxIter) const {
    const int                  cap   = (maxIter <= 0) ? 15000 : maxIter;
    const std::vector<double>  diag  = extractDiag(rows_);
    const double               bnorm = l2norm(rhs_);

    // Initialise x, r, z, p, rho from the Jacobi-scaled initial guess.
    CGState s = initCG(rows_, rhs_, diag);

    for (int iter = 0; iter < cap; ++iter) {
        // Convergence check (SMCR::result, smcr.h line 911):
        // stop when ‖r‖ < tol·‖b‖ or the system is trivially zero.
        const double rnorm   = l2norm(s.r);
        const bool   trivial = (bnorm < 1e-300 && rnorm < 1e-300);
        if (rnorm < tol * bnorm || trivial) { break; }

        // One CG step (compute alpha, update x/r, apply Jacobi, update p/rho).
        if (!cgStep(rows_, diag, s)) { break; }  // singular / ill-conditioned
    }

    return s.x;
}

}  // namespace decimate::math
