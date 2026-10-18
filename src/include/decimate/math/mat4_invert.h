// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file mat4_invert.h
/// @brief 4x4 double-precision matrix inversion — extracted from legacy math3D.h.
///
/// Pulled from `legacy/math/math3D.h` (`invertMatrixd` function) to satisfy
/// `QError::Optimize`'s need for a 4×4 inverse without dragging all of
/// `math3D.h` (which pulls `windows.h` and OpenGL) into core.
///
/// ## What was pulled
/// The `invertMatrixd(double* out, const double* m)` inline, verbatim — column-major
/// layout, Gaussian elimination with partial pivoting, returning false on singular
/// matrices.  Nothing else from `math3D.h` is included here.
///
/// ## Column-major layout (matching legacy)
/// `MAT(m, r, c) = m[c*4 + r]` — OpenGL convention as used by the original code.
/// `QError::optimize()` builds its matrix in the same layout.

#include <cmath>

namespace decimate::math {

/// @brief Invert a 4×4 double-precision column-major matrix.
///
/// Ported verbatim from `legacy/math/math3D.h::invertMatrixd`.
/// Uses Gaussian elimination with partial pivoting.
///
/// @param out  Receives the inverse (column-major, 16 doubles) on success.
/// @param m    Input matrix (column-major, 16 doubles).
/// @returns    `true` on success; `false` if the matrix is singular.
///
/// Numeric kernel: ported verbatim from the 2003 source. The unrolled
/// Gauss-Jordan body trips function-size / cognitive-complexity and the
/// C-array workspace is intrinsic to the column-major contract — both are
/// deliberately preserved (cpp-conventions permits NOLINT for numeric kernels).
// NOLINTBEGIN(readability-function-size,readability-function-cognitive-complexity,modernize-avoid-c-arrays)
inline bool invertMatrixd(double* out, const double* m) noexcept {
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define SWAP_ROWS_(a, b) \
    do {                 \
        double* _tmp = (a); \
        (a)          = (b); \
        (b)          = _tmp; \
    } while (false)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAT_(m_, r, c) (m_)[(c) * 4 + (r)]

    double  wtmp[4][8];
    double  m0 = 0.0;
    double  m1 = 0.0;
    double  m2 = 0.0;
    double  m3 = 0.0;
    double  s  = 0.0;
    double* r0 = wtmp[0];
    double* r1 = wtmp[1];
    double* r2 = wtmp[2];
    double* r3 = wtmp[3];

    r0[0] = MAT_(m, 0, 0); r0[1] = MAT_(m, 0, 1);
    r0[2] = MAT_(m, 0, 2); r0[3] = MAT_(m, 0, 3);
    r0[4] = 1.0; r0[5] = r0[6] = r0[7] = 0.0;

    r1[0] = MAT_(m, 1, 0); r1[1] = MAT_(m, 1, 1);
    r1[2] = MAT_(m, 1, 2); r1[3] = MAT_(m, 1, 3);
    r1[5] = 1.0; r1[4] = r1[6] = r1[7] = 0.0;

    r2[0] = MAT_(m, 2, 0); r2[1] = MAT_(m, 2, 1);
    r2[2] = MAT_(m, 2, 2); r2[3] = MAT_(m, 2, 3);
    r2[6] = 1.0; r2[4] = r2[5] = r2[7] = 0.0;

    r3[0] = MAT_(m, 3, 0); r3[1] = MAT_(m, 3, 1);
    r3[2] = MAT_(m, 3, 2); r3[3] = MAT_(m, 3, 3);
    r3[7] = 1.0; r3[4] = r3[5] = r3[6] = 0.0;

    /* choose pivot - or die */
    if (std::fabs(r3[0]) > std::fabs(r2[0])) { SWAP_ROWS_(r3, r2); }
    if (std::fabs(r2[0]) > std::fabs(r1[0])) { SWAP_ROWS_(r2, r1); }
    if (std::fabs(r1[0]) > std::fabs(r0[0])) { SWAP_ROWS_(r1, r0); }
    if (0.0 == r0[0]) { return false; }

    /* eliminate first variable */
    m1 = r1[0] / r0[0]; m2 = r2[0] / r0[0]; m3 = r3[0] / r0[0];
    s = r0[1]; r1[1] -= m1 * s; r2[1] -= m2 * s; r3[1] -= m3 * s;
    s = r0[2]; r1[2] -= m1 * s; r2[2] -= m2 * s; r3[2] -= m3 * s;
    s = r0[3]; r1[3] -= m1 * s; r2[3] -= m2 * s; r3[3] -= m3 * s;
    s = r0[4];
    if (s != 0.0) { r1[4] -= m1 * s; r2[4] -= m2 * s; r3[4] -= m3 * s; }
    s = r0[5];
    if (s != 0.0) { r1[5] -= m1 * s; r2[5] -= m2 * s; r3[5] -= m3 * s; }
    s = r0[6];
    if (s != 0.0) { r1[6] -= m1 * s; r2[6] -= m2 * s; r3[6] -= m3 * s; }
    s = r0[7];
    if (s != 0.0) { r1[7] -= m1 * s; r2[7] -= m2 * s; r3[7] -= m3 * s; }

    /* choose pivot - or die */
    if (std::fabs(r3[1]) > std::fabs(r2[1])) { SWAP_ROWS_(r3, r2); }
    if (std::fabs(r2[1]) > std::fabs(r1[1])) { SWAP_ROWS_(r2, r1); }
    if (0.0 == r1[1]) { return false; }

    /* eliminate second variable */
    m2 = r2[1] / r1[1]; m3 = r3[1] / r1[1];
    r2[2] -= m2 * r1[2]; r3[2] -= m3 * r1[2];
    r2[3] -= m2 * r1[3]; r3[3] -= m3 * r1[3];
    s = r1[4]; if (0.0 != s) { r2[4] -= m2 * s; r3[4] -= m3 * s; }
    s = r1[5]; if (0.0 != s) { r2[5] -= m2 * s; r3[5] -= m3 * s; }
    s = r1[6]; if (0.0 != s) { r2[6] -= m2 * s; r3[6] -= m3 * s; }
    s = r1[7]; if (0.0 != s) { r2[7] -= m2 * s; r3[7] -= m3 * s; }

    /* choose pivot - or die */
    if (std::fabs(r3[2]) > std::fabs(r2[2])) { SWAP_ROWS_(r3, r2); }
    if (0.0 == r2[2]) { return false; }

    /* eliminate third variable */
    m3 = r3[2] / r2[2];
    r3[3] -= m3 * r2[3];
    r3[4] -= m3 * r2[4];
    r3[5] -= m3 * r2[5];
    r3[6] -= m3 * r2[6];
    r3[7] -= m3 * r2[7];

    /* last check */
    if (0.0 == r3[3]) { return false; }

    s = 1.0 / r3[3]; /* back substitute row 3 */
    r3[4] *= s; r3[5] *= s; r3[6] *= s; r3[7] *= s;

    m2 = r2[3]; /* back substitute row 2 */
    s  = 1.0 / r2[2];
    r2[4] = s * (r2[4] - r3[4] * m2);
    r2[5] = s * (r2[5] - r3[5] * m2);
    r2[6] = s * (r2[6] - r3[6] * m2);
    r2[7] = s * (r2[7] - r3[7] * m2);
    m1 = r1[3];
    r1[4] -= r3[4] * m1; r1[5] -= r3[5] * m1;
    r1[6] -= r3[6] * m1; r1[7] -= r3[7] * m1;
    m0 = r0[3];
    r0[4] -= r3[4] * m0; r0[5] -= r3[5] * m0;
    r0[6] -= r3[6] * m0; r0[7] -= r3[7] * m0;

    m1 = r1[2]; /* back substitute row 1 */
    s  = 1.0 / r1[1];
    r1[4] = s * (r1[4] - r2[4] * m1);
    r1[5] = s * (r1[5] - r2[5] * m1);
    r1[6] = s * (r1[6] - r2[6] * m1);
    r1[7] = s * (r1[7] - r2[7] * m1);
    m0 = r0[2];
    r0[4] -= r2[4] * m0; r0[5] -= r2[5] * m0;
    r0[6] -= r2[6] * m0; r0[7] -= r2[7] * m0;

    m0 = r0[1]; /* back substitute row 0 */
    s  = 1.0 / r0[0];
    r0[4] = s * (r0[4] - r1[4] * m0);
    r0[5] = s * (r0[5] - r1[5] * m0);
    r0[6] = s * (r0[6] - r1[6] * m0);
    r0[7] = s * (r0[7] - r1[7] * m0);

    MAT_(out, 0, 0) = r0[4]; MAT_(out, 0, 1) = r0[5];
    MAT_(out, 0, 2) = r0[6]; MAT_(out, 0, 3) = r0[7];
    MAT_(out, 1, 0) = r1[4]; MAT_(out, 1, 1) = r1[5];
    MAT_(out, 1, 2) = r1[6]; MAT_(out, 1, 3) = r1[7];
    MAT_(out, 2, 0) = r2[4]; MAT_(out, 2, 1) = r2[5];
    MAT_(out, 2, 2) = r2[6]; MAT_(out, 2, 3) = r2[7];
    MAT_(out, 3, 0) = r3[4]; MAT_(out, 3, 1) = r3[5];
    MAT_(out, 3, 2) = r3[6]; MAT_(out, 3, 3) = r3[7];

    return true;

#undef MAT_
#undef SWAP_ROWS_
}
// NOLINTEND(readability-function-size,readability-function-cognitive-complexity,modernize-avoid-c-arrays)

}  // namespace decimate::math

namespace decimate { using math::invertMatrixd; }  // transitional re-export