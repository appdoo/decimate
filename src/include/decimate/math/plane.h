// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <cstdint>

#include <decimate/math/vec3.h>

namespace decimate::math {

/// Oriented plane (normal + signed distance) with precomputed sign bits for the fast
/// AABB near/far-corner culling test.
///
/// The `signbits` field encodes the sign of each normal component as a 3-bit
/// index (bit2 = x, bit1 = y, bit0 = z; bit set ⇔ that component is negative).
/// Combined with the AABB corner layout (corner index bit set ⇔ min on that
/// axis), this lets AABB::farCorner / AABB::nearCorner pick, in O(1) with no
/// branches, the box corner farthest / nearest along the plane normal — the
/// p-vertex / n-vertex used by the standard fast AABB-vs-plane reject.
///
/// Invariant: `signbits` always matches `nrm`. Assigning `nrm` directly breaks it;
/// use setNrm(). The normal is NOT required to be unit length, but then signedDistance()
/// is scaled by |nrm|.
///
/// Port of the 2003-2019 mudvis/vermessung `vm::Plane` (Apache 2.0, Andor Greissl);
/// German doc-comments translated to English.
template <typename T>
class Plane {
public:
    Vec3<T> nrm;               ///< Plane normal.
    T       dist = T(0);       ///< Signed distance: a point p is on the
                               ///  positive side when nrm·p - dist >= 0.
    std::uint8_t signbits = 0; ///< Per-axis sign of the normal (see class doc).

    /// Degenerate plane: zero normal, distance 0, signbits 0.
    Plane() = default;

    /// Plane from @p normal and signed @p distance (see `dist`); precomputes signbits.
    Plane(const Vec3<T> &normal, T distance)
        : nrm(normal), dist(distance), signbits(calculateSignbits(normal)) {
    }

    /// Set the normal and recompute the sign-bit corner index.
    void setNrm(const Vec3<T> &normal) {
        nrm = normal;
        signbits = calculateSignbits(normal);
    }

    /// Signed distance from point `p` to the plane (>= 0 ⇒ positive side).
    [[nodiscard]] T signedDistance(const Vec3<T> &p) const {
        return nrm.dot(p) - dist;
    }

private:
    /// Pack the sign of each normal component into a 3-bit index.
    /// bit2 = sign(x), bit1 = sign(y), bit0 = sign(z); 1 ⇔ negative.
    static std::uint8_t calculateSignbits(const Vec3<T> &n) {
        return static_cast<std::uint8_t>((sgn(n.x) << 2) | (sgn(n.y) << 1) | sgn(n.z));
    }

    /// 1 for a negative component, else 0 (so -0 and NaN count as non-negative).
    static std::uint8_t sgn(T x) {
        return x < T(0) ? 1 : 0;
    }
};


using Planef = Plane<float>;
using Planed = Plane<double>;

} // namespace decimate::math
