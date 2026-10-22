// This file is part of the decimate library; see decimate.h for version/license details.
/// @file qerror.cpp
/// @brief `QError` implementation — Garland-Heckbert 1997 symmetric quadric.

#include <decimate/reduce/qerror.h>

#include <cmath> // std::fabs

namespace decimate {

QError::QError() noexcept {
    init();
}


QError::QError(const QError &q1, const QError &q2) noexcept {
    init(q1, q2);
}


void QError::init() noexcept {
    aa_ = ab_ = ac_ = ad_ = 0.0;
    bb_ = bc_ = bd_ = 0.0;
    cc_ = cd_ = 0.0;
    dd_ = 0.0;
}


void QError::init(const QError &q1, const QError &q2) noexcept {
    aa_ = q1.aa_ + q2.aa_;
    ab_ = q1.ab_ + q2.ab_;
    ac_ = q1.ac_ + q2.ac_;
    ad_ = q1.ad_ + q2.ad_;
    bb_ = q1.bb_ + q2.bb_;
    bc_ = q1.bc_ + q2.bc_;
    bd_ = q1.bd_ + q2.bd_;
    cc_ = q1.cc_ + q2.cc_;
    cd_ = q1.cd_ + q2.cd_;
    dd_ = q1.dd_ + q2.dd_;
}


/// Garland-Heckbert K_p matrix from the face plane (legacy2003/reduce/qerror.cpp:72).
///
/// Operator mapping vs legacy:
///   legacy `a ^ b`       → a.cross(b)   (legacy ^ was cross product)
///   legacy `nrm.norm()`  → nrm.normalize()
///   legacy `-(nrm.x*p1.x + ...)` → -static_cast<double>(nrm.dot(p1))
///                                    (legacy * on two Vec3s was dot product)
/// All arithmetic is promoted to double for the matrix coefficients.
void QError::setPlane(const math::Vec3f &p1,
                      const math::Vec3f &p2,
                      const math::Vec3f &p3) noexcept {
    const auto a = p2 - p1;
    const auto b = p3 - p1;

    // A degenerate face keeps a zero normal (normalize() leaves it alone), so it adds no error.
    auto nrm = a.cross(b);
    nrm.normalize();

    // Plane offset d = -(n · p1), the 4th component of the plane vector.
    const double d = -static_cast<double>(nrm.dot(p1));

    // Upper-triangle K_p coefficients: n·nᵀ with d in the 4th component (legacy2003/reduce/qerror.cpp:82-98).
    const double nx = static_cast<double>(nrm.x);
    const double ny = static_cast<double>(nrm.y);
    const double nz = static_cast<double>(nrm.z);

    aa_ = nx * nx;
    ab_ = nx * ny;
    ac_ = nx * nz;
    ad_ = nx * d;
    bb_ = ny * ny;
    bc_ = ny * nz;
    bd_ = ny * d;
    cc_ = nz * nz;
    cd_ = nz * d;
    dd_ = d * d;
}


void QError::add(const QError &other) noexcept {
    aa_ += other.aa_;
    ab_ += other.ab_;
    ac_ += other.ac_;
    ad_ += other.ad_;
    bb_ += other.bb_;
    bc_ += other.bc_;
    bd_ += other.bd_;
    cc_ += other.cc_;
    cd_ += other.cd_;
    dd_ += other.dd_;
}


/// pᵀ Q p, evaluated as a quadratic form (legacy2003/reduce/qerror.cpp:137).
double QError::getError(const math::Vec3f &p) const noexcept {
    const double px = static_cast<double>(p.x);
    const double py = static_cast<double>(p.y);
    const double pz = static_cast<double>(p.z);

    // Q·(p,1) row by row, then dot with (p,1).
    const double x = aa_ * px + ab_ * py + ac_ * pz + ad_;
    const double y = ab_ * px + bb_ * py + bc_ * pz + bd_;
    const double z = ac_ * px + bc_ * py + cc_ * pz + cd_;
    const double w = ad_ * px + bd_ * py + cd_ * pz + dd_;

    return x * px + y * py + z * pz + w;
}


/// Symmetric 3×3 solve for the cost-minimising position.
///
/// Garland-Heckbert §3: the optimal position solves  A·v = -b  where A is the quadric's upper-left
/// 3×3 block and b = (ad, bd, cd). The legacy path (legacy2003/reduce/qerror.cpp:162) built a 4×4 with
/// bottom row [0 0 0 1] and inverted it, taking the last column. That det equals det(A) and that last
/// column is exactly this solution, so this is the same result via the symmetric-cofactor (Cramer)
/// rule — ~10× fewer FLOPs (the 4×4 inverse was the #1 reduce hotspot, ~28% of runtime).
///
/// A is symmetric ⇒ only 6 distinct cofactors. Singular (det≈0) ⇒ return false; the caller falls back
/// to the edge midpoint. NOT bit-identical to the old arithmetic — float order differs — so the reduce
/// golden is re-blessed.
bool QError::optimize(math::Vec3f &p) const noexcept {
    // First-row cofactors give the determinant (adj(A) is symmetric: 6 unique entries in all).
    const double c00 = bb_ * cc_ - bc_ * bc_;
    const double c01 = ac_ * bc_ - ab_ * cc_;
    const double c02 = ab_ * bc_ - ac_ * bb_;

    const double     det = aa_ * c00 + ab_ * c01 + ac_ * c02;
    constexpr double kDetEps = 1e-10; // matches the legacy "fails only when degenerate"
    if (std::fabs(det) < kDetEps) {
        return false; // singular ⇒ caller uses the midpoint
    }

    // Remaining cofactors, only needed once the solve is known to exist.
    const double c11 = aa_ * cc_ - ac_ * ac_;
    const double c12 = ab_ * ac_ - aa_ * bc_;
    const double c22 = aa_ * bb_ - ab_ * ab_;

    const double inv = 1.0 / det;
    // v = -adj(A)·b / det,  with b = (ad, bd, cd).
    p.x = static_cast<float>(-(c00 * ad_ + c01 * bd_ + c02 * cd_) * inv);
    p.y = static_cast<float>(-(c01 * ad_ + c11 * bd_ + c12 * cd_) * inv);
    p.z = static_cast<float>(-(c02 * ad_ + c12 * bd_ + c22 * cd_) * inv);
    return true;
}

} // namespace decimate
