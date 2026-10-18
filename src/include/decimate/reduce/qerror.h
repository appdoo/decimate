// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file qerror.h
/// @brief Symmetric 4×4 quadric error matrix — Garland-Heckbert 1997, §3.
///
/// `QError` stores the fundamental error metric for one vertex as the upper
/// triangle of a symmetric 4×4 matrix (10 doubles), exactly matching the
/// legacy `decimate::QError` storage in `legacy/reduce/qerror.h`.
///
/// The mathematics are the K_p matrix from the paper:
///
///     K_p = n·nᵀ  (outer product of the normalised face-plane normal with
///                  the plane offset d folded into the 4th component)
///
/// Quadric error for point p: `e(p) = pᵀ Q p`
///
/// ## Operator mapping (Beat-2a correctness)
///
/// Legacy `Vec3` used non-standard operators that differ from `decimate::math::Vec3`:
///   - Legacy `a ^ b`          = cross product → modern `a.cross(b)`
///   - Legacy `a * b` (Vec3×Vec3) = dot product → modern `a.dot(b)`
///   - Legacy `a * scalar`     = scalar multiply → modern `a * scalar` (same)
///
/// Both critical sites (`setPlane` and `optimize`) are ported with explicit
/// `.cross()` and `.dot()` calls and a side-by-side note in the .cpp.

#include <decimate/math/vec3.h>

namespace decimate {

/// @brief Symmetric 4×4 error quadric, stored as 10 doubles (upper triangle).
///
/// Fields follow the naming convention from the Garland-Heckbert paper:
///
/// @verbatim
/// | aa ab ac ad |
/// | ab bb bc bd |   ← symmetric — only upper triangle stored
/// | ac bc cc cd |
/// | ad bd cd dd |
/// @endverbatim
///
/// This is a pure-computation value type: no virtual functions, no allocation,
/// no OS dependencies.  It is aggregated by `QuadricErrorSystem`.
class QError {
public:
    /// @brief Construct a zero-initialized quadric.
    QError() noexcept;

    /// @brief Construct as the sum of two quadrics.
    ///
    /// Equivalent to `QError q; q.init(q1, q2);` — used by
    /// `QuadricErrorSystem::optimize` to form the combined edge quadric
    /// on the stack without touching the per-vertex storage.
    QError(const QError& q1, const QError& q2) noexcept;

    /// @brief Zero all coefficients.
    void init() noexcept;

    /// @brief Set to the element-wise sum of @p q1 and @p q2.
    void init(const QError& q1, const QError& q2) noexcept;

    /// @brief Set the quadric from the plane through three vertex positions.
    ///
    /// Computes the face-plane normal `n = (p2-p1) × (p3-p1)` (normalised),
    /// then sets the K_p matrix coefficients to `n·nᵀ` with the plane-offset
    /// `d = -n·p1` folded into the last row/column.
    ///
    /// @param p1  First  vertex position.
    /// @param p2  Second vertex position.
    /// @param p3  Third  vertex position.
    void setPlane(const decimate::math::Vec3f& p1,
                  const decimate::math::Vec3f& p2,
                  const decimate::math::Vec3f& p3) noexcept;

    /// @brief Accumulate another quadric into this one (in-place addition).
    ///
    /// Used both during initialisation (each vertex gathers its incident
    /// face quadrics) and after a collapse (`QuadricErrorSystem::update`).
    void add(const QError& other) noexcept;

    /// @brief Evaluate the quadric error `pᵀ Q p` for point @p p.
    ///
    /// Returns the sum of squared distances from @p p to all accumulated
    /// face planes.  Always non-negative for a valid quadric.
    ///
    /// @param p  Query point.
    [[nodiscard]] double getError(const decimate::math::Vec3f& p) const noexcept;

    /// @brief Solve for the cost-minimising position by inverting the quadric.
    ///
    /// Assembles the full 4×4 matrix, replaces the bottom row with
    /// `[0 0 0 1]` (as per Garland-Heckbert §3), and calls `invertMatrixd`.
    /// The optimal position is the last column of the inverse.
    ///
    /// Falls back to `false` if the matrix is singular (caller must use
    /// midpoint fallback — this is the Beat-2a behaviour, F4 is deferred).
    ///
    /// @param[out] p  Receives the optimal position on success.
    /// @returns   `true` if the inverse exists; `false` if singular.
    bool optimize(decimate::math::Vec3f& p) const noexcept;

private:
    // Upper triangle of the symmetric 4×4 error matrix.
    // Row-major naming: aa=Q[0][0], ab=Q[0][1], etc.
    double aa_ = 0.0;
    double ab_ = 0.0;
    double ac_ = 0.0;
    double ad_ = 0.0;
    double bb_ = 0.0;
    double bc_ = 0.0;
    double bd_ = 0.0;
    double cc_ = 0.0;
    double cd_ = 0.0;
    double dd_ = 0.0;
};

}  // namespace decimate
