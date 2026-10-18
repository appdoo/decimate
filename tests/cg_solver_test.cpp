// This file is part of the decimate library; see decimate.h for version/license details.
/// @file cg_solver_test.cpp
/// @brief Catch2 tests for ConjugateGradientSolver (Spec 08-02, S1 + S2).
///
/// Test scenarios assert numeric properties of the solver output —
/// component values, residual norms, iteration counts — never text presence.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <decimate/math/cg_solver.h>

#include <cmath>
#include <vector>

using decimate::math::ConjugateGradientSolver;

// ---------------------------------------------------------------------------
// S1 — Known SPD systems with analytic solutions
// ---------------------------------------------------------------------------

TEST_CASE("S1a: 2x2 SPD system [[4,1],[1,3]]*x=[1,2]", "[cg_solver][S1]") {
    // Analytic solution:
    //   det(A) = 4*3 - 1*1 = 11
    //   x0 = (3*1 - 1*2) / 11 = 1/11  ≈  0.090909090909
    //   x1 = (4*2 - 1*1) / 11 = 7/11  ≈  0.636363636364
    constexpr double kX0 = 1.0 / 11.0;
    constexpr double kX1 = 7.0 / 11.0;

    ConjugateGradientSolver solver{2};
    // Row 0: 4*x0 + 1*x1 = 1
    solver.addEntry(0, 0, 4.0);
    solver.addEntry(0, 1, 1.0);
    // Row 1: 1*x0 + 3*x1 = 2
    solver.addEntry(1, 0, 1.0);
    solver.addEntry(1, 1, 3.0);
    solver.setRhs(0, 1.0);
    solver.setRhs(1, 2.0);

    const std::vector<double> x = solver.solve(1e-12);

    REQUIRE(x.size() == 2);
    CHECK_THAT(x[0], Catch::Matchers::WithinAbs(kX0, 1e-6));
    CHECK_THAT(x[1], Catch::Matchers::WithinAbs(kX1, 1e-6));

    // Verify residual ‖Ax−b‖ < 1e-8 using hand-computed Ax.
    // A = [[4,1],[1,3]]
    const double r0 = (4.0 * x[0] + 1.0 * x[1]) - 1.0;
    const double r1 = (1.0 * x[0] + 3.0 * x[1]) - 2.0;
    const double residual = std::sqrt(r0 * r0 + r1 * r1);
    CHECK(residual < 1e-8);
}

TEST_CASE("S1b: 3x3 SPD system [[4,1,0],[1,3,1],[0,1,2]]*x=[1,1,1]", "[cg_solver][S1]") {
    // Analytic solution (by elimination):
    //   Row 0: 4*x0 + x1 = 1   → x1 = 1 - 4*x0
    //   Row 2: x1 + 2*x2 = 1   → x2 = (1 - x1)/2 = 2*x0
    //   Row 1: x0 + 3*(1-4*x0) + 2*x0 = 1
    //          x0 + 3 - 12*x0 + 2*x0 = 1
    //          -9*x0 = -2  →  x0 = 2/9
    //   x1 = 1 - 8/9 = 1/9
    //   x2 = 4/9
    //
    // Verification: row1: 2/9 + 3/9 + 4/9 = 9/9 = 1.  ✓
    constexpr double kX0 = 2.0 / 9.0;
    constexpr double kX1 = 1.0 / 9.0;
    constexpr double kX2 = 4.0 / 9.0;

    ConjugateGradientSolver solver{3};
    // Row 0: 4*x0 + 1*x1
    solver.addEntry(0, 0, 4.0);
    solver.addEntry(0, 1, 1.0);
    // Row 1: 1*x0 + 3*x1 + 1*x2
    solver.addEntry(1, 0, 1.0);
    solver.addEntry(1, 1, 3.0);
    solver.addEntry(1, 2, 1.0);
    // Row 2: 1*x1 + 2*x2
    solver.addEntry(2, 1, 1.0);
    solver.addEntry(2, 2, 2.0);
    solver.setRhs(0, 1.0);
    solver.setRhs(1, 1.0);
    solver.setRhs(2, 1.0);

    const std::vector<double> x = solver.solve(1e-12);

    REQUIRE(x.size() == 3);
    CHECK_THAT(x[0], Catch::Matchers::WithinAbs(kX0, 1e-6));
    CHECK_THAT(x[1], Catch::Matchers::WithinAbs(kX1, 1e-6));
    CHECK_THAT(x[2], Catch::Matchers::WithinAbs(kX2, 1e-6));

    // Residual ‖Ax−b‖ < 1e-8
    const double r0 = (4.0 * x[0] + 1.0 * x[1]) - 1.0;
    const double r1 = (1.0 * x[0] + 3.0 * x[1] + 1.0 * x[2]) - 1.0;
    const double r2 = (1.0 * x[1] + 2.0 * x[2]) - 1.0;
    const double residual = std::sqrt(r0 * r0 + r1 * r1 + r2 * r2);
    CHECK(residual < 1e-8);
}

TEST_CASE("S1c: addEntry accumulates at same (row,col)", "[cg_solver][S1]") {
    // Adding the same entry twice should be equivalent to adding once with sum.
    // [[4,1],[1,3]] assembled via two addEntry(0,0,2.0) calls.
    ConjugateGradientSolver solver{2};
    solver.addEntry(0, 0, 2.0);
    solver.addEntry(0, 0, 2.0);  // cumulative: A(0,0) = 4
    solver.addEntry(0, 1, 1.0);
    solver.addEntry(1, 0, 1.0);
    solver.addEntry(1, 1, 3.0);
    solver.setRhs(0, 1.0);
    solver.setRhs(1, 2.0);

    const std::vector<double> x = solver.solve(1e-12);

    REQUIRE(x.size() == 2);
    CHECK_THAT(x[0], Catch::Matchers::WithinAbs(1.0 / 11.0, 1e-6));
    CHECK_THAT(x[1], Catch::Matchers::WithinAbs(7.0 / 11.0, 1e-6));
}

// ---------------------------------------------------------------------------
// S2 — Preconditioner and convergence properties
// ---------------------------------------------------------------------------

TEST_CASE("S2a: diagonal system converges in ≤ n iterations", "[cg_solver][S2]") {
    // A purely diagonal SPD system A = diag(1, 2, 3, 4) with b = [1, 2, 3, 4].
    // Exact solution: x = [1, 1, 1, 1].
    // For a diagonal system the Jacobi preconditioner makes A→I, so CG
    // converges in exactly 1 iteration (all eigenvalues collapse to 1).
    // We assert convergence within n = 4 iterations.
    constexpr int n = 4;
    ConjugateGradientSolver solver{n};
    for (int i = 0; i < n; ++i) {
        solver.addEntry(i, i, static_cast<double>(i + 1));
        solver.setRhs(i, static_cast<double>(i + 1));
    }

    // Solve with a tight iteration cap = n to verify fast convergence.
    const std::vector<double> x = solver.solve(1e-12, n);

    REQUIRE(x.size() == static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        CHECK_THAT(x[static_cast<std::size_t>(i)], Catch::Matchers::WithinAbs(1.0, 1e-10));
    }
}

TEST_CASE("S2b: diagonally-dominant SPD converges to ‖r‖ < 1e-10", "[cg_solver][S2]") {
    // A = [[10, 1, 0.5],
    //      [ 1, 8, 0.5],
    //      [0.5, 0.5, 6]]
    // b = [1, 1, 1]
    // All diagonal entries dominate their rows → SPD and well-conditioned.
    // Assert: solver converges (residual < 1e-10 after default maxIter).

    ConjugateGradientSolver solver{3};
    // Row 0
    solver.addEntry(0, 0, 10.0);
    solver.addEntry(0, 1, 1.0);
    solver.addEntry(0, 2, 0.5);
    // Row 1
    solver.addEntry(1, 0, 1.0);
    solver.addEntry(1, 1, 8.0);
    solver.addEntry(1, 2, 0.5);
    // Row 2
    solver.addEntry(2, 0, 0.5);
    solver.addEntry(2, 1, 0.5);
    solver.addEntry(2, 2, 6.0);

    solver.setRhs(0, 1.0);
    solver.setRhs(1, 1.0);
    solver.setRhs(2, 1.0);

    const std::vector<double> x = solver.solve(1e-12);

    REQUIRE(x.size() == 3);

    // Compute residual ‖Ax−b‖ directly.
    const double r0 = (10.0 * x[0] + 1.0 * x[1] + 0.5 * x[2]) - 1.0;
    const double r1 = (1.0 * x[0] + 8.0 * x[1] + 0.5 * x[2]) - 1.0;
    const double r2 = (0.5 * x[0] + 0.5 * x[1] + 6.0 * x[2]) - 1.0;
    const double residual = std::sqrt(r0 * r0 + r1 * r1 + r2 * r2);

    CHECK(residual < 1e-10);
}

TEST_CASE("S2c: size() returns the system dimension", "[cg_solver][S2]") {
    // Trivial size check — the constructor stores n correctly.
    ConjugateGradientSolver s5{5};
    CHECK(s5.size() == 5);

    ConjugateGradientSolver s1{1};
    CHECK(s1.size() == 1);
}

TEST_CASE("S2d: solve with maxIter=0 uses legacy default 15000", "[cg_solver][S2]") {
    // Passing maxIter=0 must select the legacy cap of 15000, not stop after 0
    // iterations. Verify by solving a 2x2 system that needs >1 iteration.
    ConjugateGradientSolver solver{2};
    solver.addEntry(0, 0, 4.0);
    solver.addEntry(0, 1, 1.0);
    solver.addEntry(1, 0, 1.0);
    solver.addEntry(1, 1, 3.0);
    solver.setRhs(0, 1.0);
    solver.setRhs(1, 2.0);

    // maxIter = 0 → default 15000; must produce the correct answer.
    const std::vector<double> x = solver.solve(1e-12, 0);

    REQUIRE(x.size() == 2);
    CHECK_THAT(x[0], Catch::Matchers::WithinAbs(1.0 / 11.0, 1e-6));
    CHECK_THAT(x[1], Catch::Matchers::WithinAbs(7.0 / 11.0, 1e-6));
}
