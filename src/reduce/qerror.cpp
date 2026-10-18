// This file is part of the decimate library; see decimate.h for version/license details.
/// @file qerror.cpp
/// @brief `QError` implementation — Garland-Heckbert 1997 symmetric quadric.

// NOLINTBEGIN(modernize-avoid-c-arrays) — 2D C arrays required by invertMatrixd
// which takes `double*` and indexes with the OpenGL column-major MAT(m,r,c)=m[c*4+r]
// macro.  Using std::array<std::array<double,4>,4> would add reinterpret_cast on
// every call and obscure the already-documented layout.  The narrow exception is
// justified: this is a single numeric-kernel private function.

/// ## Operator mapping — side-by-side legacy vs port (Beat-2a requirement)
///
/// ### `setPlane` (was `QError::SetPlane`)
///
/// Legacy (qerror.cpp:72):
/// @code
///   Vec3 a = p2 - p1;
///   Vec3 b = p3 - p1;
///   Vec3 nrm = a ^ b;   // ^ is CROSS PRODUCT in legacy Vec3
///   nrm.norm();
///   double d = -(nrm.x * p1.x + nrm.y * p1.y + nrm.z * p1.z);
/// @endcode
///
/// Port:
/// @code
///   auto a = p2 - p1;
///   auto b = p3 - p1;
///   auto nrm = a.cross(b);    // .cross() replaces ^
///   nrm.normalize();
///   double d = -static_cast<double>(nrm.dot(p1));   // .dot() replaces *
/// @endcode
///
/// ### `optimize` (was `QError::Optimize`)
///
/// Legacy (qerror.cpp:162): fills `double m[4][4]` with the quadric
/// coefficients (symmetric), replaces bottom row with `[0 0 0 1]`, calls
/// `invertMatrixd((double*)inverse, (double*)m)`, then reads column 3 of the
/// inverse as the optimal vertex.
///
/// Port: identical logic — same 4×4 layout, same column-major `invertMatrixd`
/// from `mat4_invert.h`.  No operator mapping issues here (pure doubles).
///
/// `getError` uses only scalar multiplications — no operator mapping needed.

#include <decimate/reduce/qerror.h>

#include <cmath>  // std::fabs

namespace decimate {

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

QError::QError() noexcept {
    init();
}

QError::QError(const QError& q1, const QError& q2) noexcept {
    init(q1, q2);
}

// ---------------------------------------------------------------------------
// init
// ---------------------------------------------------------------------------

void QError::init() noexcept {
    aa_ = ab_ = ac_ = ad_ = 0.0;
    bb_ = bc_ = bd_ = 0.0;
    cc_ = cd_ = 0.0;
    dd_ = 0.0;
}

void QError::init(const QError& q1, const QError& q2) noexcept {
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

// ---------------------------------------------------------------------------
// setPlane — Garland-Heckbert K_p matrix from face plane (qerror.cpp:72)
// ---------------------------------------------------------------------------
// Operator mapping vs legacy:
//   legacy `a ^ b`       → a.cross(b)   (legacy ^ was cross product)
//   legacy `nrm.norm()`  → nrm.normalize()
//   legacy `-(nrm.x*p1.x + ...)` → -static_cast<double>(nrm.dot(p1))
//                                    (legacy * on two Vec3s was dot product)
// All arithmetic is promoted to double for the matrix coefficients.

void QError::setPlane(const decimate::math::Vec3f& p1,
                      const decimate::math::Vec3f& p2,
                      const decimate::math::Vec3f& p3) noexcept {
    // Edge vectors from p1
    const auto a = p2 - p1;
    const auto b = p3 - p1;

    // Normal: cross product (legacy `a ^ b` → a.cross(b))
    auto nrm = a.cross(b);
    nrm.normalize();  // legacy nrm.norm()

    // Plane offset d = -(n · p1)
    // legacy: `-(nrm.x*p1.x + nrm.y*p1.y + nrm.z*p1.z)`
    // = -nrm.dot(p1) in modern form
    const double d = -static_cast<double>(nrm.dot(p1));

    // Fill upper-triangle K_p coefficients: nᵀ·n with d in 4th component
    // Identical to legacy qerror.cpp:82–98
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

// ---------------------------------------------------------------------------
// add — accumulate another quadric (qerror.cpp:105)
// ---------------------------------------------------------------------------

void QError::add(const QError& other) noexcept {
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

// ---------------------------------------------------------------------------
// getError — pᵀ Q p (qerror.cpp:137)
// ---------------------------------------------------------------------------
// Pure scalar arithmetic — no operator mapping needed.

double QError::getError(const decimate::math::Vec3f& p) const noexcept {
    const double px = static_cast<double>(p.x);
    const double py = static_cast<double>(p.y);
    const double pz = static_cast<double>(p.z);

    const double x = aa_ * px + ab_ * py + ac_ * pz + ad_;
    const double y = ab_ * px + bb_ * py + bc_ * pz + bd_;
    const double z = ac_ * px + bc_ * py + cc_ * pz + cd_;
    const double w = ad_ * px + bd_ * py + cd_ * pz + dd_;

    return x * px + y * py + z * pz + w;
}

// ---------------------------------------------------------------------------
// optimize — symmetric 3×3 solve for the cost-minimising position
// ---------------------------------------------------------------------------
// Garland-Heckbert §3: the optimal position solves  A·v = -b  where A is the
// quadric's upper-left 3×3 block and b = (ad, bd, cd). The legacy path built a
// 4×4 with bottom row [0 0 0 1] and inverted it (invertMatrixd), taking the last
// column. That det equals det(A) and that last column is exactly this solution,
// so this is the same result via the symmetric-cofactor (Cramer) rule — ~10×
// fewer FLOPs (the 4×4 inverse was the #1 reduce hotspot, ~28% of runtime).
//
// A is symmetric ⇒ only 6 distinct cofactors. Singular (det≈0) ⇒ return false;
// the caller falls back to the edge midpoint (same contract as the 4×4 path,
// whose pivoting also failed only on a (near-)zero pivot). NOT bit-identical to
// the old arithmetic — float order differs — so the reduce golden is re-blessed.

bool QError::optimize(decimate::math::Vec3f& p) const noexcept {
    // Cofactors of A (symmetric: adj(A) is symmetric, 6 unique entries).
    const double c00 = bb_ * cc_ - bc_ * bc_;
    const double c01 = ac_ * bc_ - ab_ * cc_;
    const double c02 = ab_ * bc_ - ac_ * bb_;

    const double det = aa_ * c00 + ab_ * c01 + ac_ * c02;
    constexpr double kDetEps = 1e-10;  // matches the legacy "fails only when degenerate"
    if (std::fabs(det) < kDetEps) {
        return false;  // singular ⇒ caller uses the midpoint
    }

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
// NOLINTEND(modernize-avoid-c-arrays)

}  // namespace decimate
