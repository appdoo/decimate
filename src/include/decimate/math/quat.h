// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <decimate/math/vec3.h>
#include <decimate/math/mat4.h>
#include <cmath>
#include <numbers>
#include <algorithm>
#include <type_traits>

namespace decimate::math {

/// Quaternion for rotations, stored as (x, y, z, w) with w the scalar part; contiguous, so
/// `data()` is 4 T in that order.
///
/// Invariant (caller's job): rotation operations (rotate, toMatrix, slerp, ...) assume a UNIT
/// quaternion; the factories return unit quaternions, but arithmetic (+, -, scalar *) does not
/// renormalize. Based on Graphics Gems quaternion articles and the DirectXTK implementation.
template <typename T>
struct Quat {
    static_assert(std::is_floating_point<T>::value, "Quat requires floating point type");

    T x, y, z, w;

    /// Identity rotation (0, 0, 0, 1).
    constexpr Quat() noexcept
        : x(0), y(0), z(0), w(1) {
    }

    /// Raw components, NOT normalized.
    constexpr Quat(T x, T y, T z, T w) noexcept
        : x(x), y(y), z(z), w(w) {
    }

    /// Rotation of @p angle radians around @p axis (normalized here; a zero axis gives NaN).
    static Quat fromAxisAngle(const Vec3<T> &axis, T angle) noexcept {
        T       halfAngle = angle * T(0.5);
        T       s = std::sin(halfAngle);
        Vec3<T> normAxis = axis.normalized();
        return Quat {
            normAxis.x * s,
            normAxis.y * s,
            normAxis.z * s,
            std::cos(halfAngle)
        };
    }

    /// Rotation from Euler angles in radians: @p roll about X, @p pitch about Y, @p yaw about Z
    /// (note the argument order differs from toEuler(), which returns (roll, pitch, yaw)).
    static Quat fromEuler(T pitch, T yaw, T roll) noexcept {
        T cy = std::cos(yaw * T(0.5));
        T sy = std::sin(yaw * T(0.5));
        T cp = std::cos(pitch * T(0.5));
        T sp = std::sin(pitch * T(0.5));
        T cr = std::cos(roll * T(0.5));
        T sr = std::sin(roll * T(0.5));

        return Quat {
            sr * cp * cy - cr * sp * sy, // x
            cr * sp * cy + sr * cp * sy, // y
            cr * cp * sy - sr * sp * cy, // z
            cr * cp * cy + sr * sp * sy  // w
        };
    }

    /// Unit quaternion from the rotation part of @p m (assumed orthonormal; translation ignored).
    /// Picks the numerically largest of w/x/y/z to divide by, so it is stable near 180 degrees.
    ///
    /// Mat4::operator[] is column-major (`m[col][row]`), while the textbook formulas below are
    /// written in R(row, col): every access therefore reads `m[col][row]`.  Reading it as
    /// `m[row][col]` returned the conjugate, i.e. the inverse rotation.
    static Quat fromMatrix(const Mat4<T> &m) noexcept {
        T    trace = m[0][0] + m[1][1] + m[2][2];
        Quat q;

        if (trace > T(0)) {
            T s = std::sqrt(trace + T(1)) * T(2);
            q.w = T(0.25) * s;
            q.x = (m[1][2] - m[2][1]) / s;
            q.y = (m[2][0] - m[0][2]) / s;
            q.z = (m[0][1] - m[1][0]) / s;
        } else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
            T s = std::sqrt(T(1) + m[0][0] - m[1][1] - m[2][2]) * T(2);
            q.w = (m[1][2] - m[2][1]) / s;
            q.x = T(0.25) * s;
            q.y = (m[1][0] + m[0][1]) / s;
            q.z = (m[2][0] + m[0][2]) / s;
        } else if (m[1][1] > m[2][2]) {
            T s = std::sqrt(T(1) + m[1][1] - m[0][0] - m[2][2]) * T(2);
            q.w = (m[2][0] - m[0][2]) / s;
            q.x = (m[1][0] + m[0][1]) / s;
            q.y = T(0.25) * s;
            q.z = (m[2][1] + m[1][2]) / s;
        } else {
            T s = std::sqrt(T(1) + m[2][2] - m[0][0] - m[1][1]) * T(2);
            q.w = (m[0][1] - m[1][0]) / s;
            q.x = (m[2][0] + m[0][2]) / s;
            q.y = (m[2][1] + m[1][2]) / s;
            q.z = T(0.25) * s;
        }

        return q.normalized();
    }

    constexpr Quat(const Quat &) noexcept = default;
    constexpr Quat(Quat &&) noexcept = default;
    constexpr Quat &operator=(const Quat &) noexcept = default;
    constexpr Quat &operator=(Quat &&) noexcept = default;

    /// Component by index 0..3 = x, y, z, w; no bounds check.
    constexpr T &operator[](size_t i) noexcept {
        return (&x)[i];
    }
    constexpr const T &operator[](size_t i) const noexcept {
        return (&x)[i];
    }

    /// Pointer to the 4 contiguous components (x, y, z, w).
    constexpr T *data() noexcept {
        return &x;
    }
    constexpr const T *data() const noexcept {
        return &x;
    }
    constexpr operator T *() noexcept {
        return &x;
    }
    constexpr operator const T *() const noexcept {
        return &x;
    }

    /// Four-component dot product; |dot| = 1 for identical rotations, sign reflects the double cover.
    constexpr T dot(const Quat &q) const noexcept {
        return x * q.x + y * q.y + z * q.z + w * q.w;
    }

    /// Euclidean norm.
    T length() const noexcept {
        return std::sqrt(x * x + y * y + z * z + w * w);
    }

    /// Squared norm; avoids the sqrt.
    constexpr T lengthSquared() const noexcept {
        return x * x + y * y + z * z + w * w;
    }

    /// Unit-length copy; a zero quaternion silently becomes identity.
    Quat normalized() const noexcept {
        T len = length();
        if (len > T(0)) {
            T invLen = T(1) / len;
            return Quat { x * invLen, y * invLen, z * invLen, w * invLen };
        }
        return Quat {}; // Identity if zero length
    }

    /// Normalizes in place; a zero quaternion is left unchanged.
    void normalize() noexcept {
        T len = length();
        if (len > T(0)) {
            T invLen = T(1) / len;
            x *= invLen;
            y *= invLen;
            z *= invLen;
            w *= invLen;
        }
    }

    /// Conjugate; equals the inverse for unit quaternions and is cheaper.
    constexpr Quat conjugate() const noexcept {
        return Quat { -x, -y, -z, w };
    }

    /// True inverse for any non-zero quaternion; a zero quaternion yields identity.
    Quat inverse() const noexcept {
        T lenSq = lengthSquared();
        if (lenSq > T(0)) {
            T invLenSq = T(1) / lenSq;
            return Quat { -x * invLenSq, -y * invLenSq, -z * invLenSq, w * invLenSq };
        }
        return Quat {}; // Identity if zero length
    }

    /// Hamilton product; `a * b` applies @p b first, then a.
    Quat operator*(const Quat &q) const noexcept {
        return Quat {
            w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y - x * q.z + y * q.w + z * q.x,
            w * q.z + x * q.y - y * q.x + z * q.w,
            w * q.w - x * q.x - y * q.y - z * q.z
        };
    }

    /// In-place Hamilton product, see operator*.
    Quat &operator*=(const Quat &q) noexcept {
        *this = *this * q;
        return *this;
    }

    /// Rotates @p v by this unit quaternion (v' = q v q^-1).
    Vec3<T> rotate(const Vec3<T> &v) const noexcept {
        // Cross-product form v + 2 * cross(q.xyz, cross(q.xyz, v) + q.w * v): no full quaternion product.
        Vec3<T> qvec(x, y, z);
        Vec3<T> uv = qvec.cross(v);
        Vec3<T> uuv = qvec.cross(uv);
        uv = uv * (T(2) * w);
        uuv = uuv * T(2);
        return v + uv + uuv;
    }

    /// Componentwise scale; the result is generally not unit length.
    constexpr Quat operator*(T scalar) const noexcept {
        return Quat { x * scalar, y * scalar, z * scalar, w * scalar };
    }

    /// Componentwise sum (not a rotation composition).
    constexpr Quat operator+(const Quat &q) const noexcept {
        return Quat { x + q.x, y + q.y, z + q.z, w + q.w };
    }

    /// Componentwise difference.
    constexpr Quat operator-(const Quat &q) const noexcept {
        return Quat { x - q.x, y - q.y, z - q.z, w - q.w };
    }

    /// Negation; represents the same rotation (q and -q are the double cover).
    constexpr Quat operator-() const noexcept {
        return Quat { -x, -y, -z, -w };
    }

    /// Exact componentwise comparison; q and -q compare unequal though they are the same rotation.
    constexpr bool operator==(const Quat &q) const noexcept {
        return x == q.x && y == q.y && z == q.z && w == q.w;
    }

    constexpr bool operator!=(const Quat &q) const noexcept {
        return !(*this == q);
    }

    /// Normalized linear interpolation toward @p q: cheap, constant direction but not constant angular velocity.
    /// Does NOT take the shortest arc; see slerp() for that.
    constexpr Quat lerp(const Quat &q, T t) const noexcept {
        return (*this * (T(1) - t) + q * t).normalized();
    }

    /// Alias of lerp(), which already normalizes.
    Quat nlerp(const Quat &q, T t) const noexcept {
        return lerp(q, t);
    }

    /// Spherical linear interpolation toward @p q with constant angular velocity, along the
    /// shortest arc (after Graphics Gems / DirectXTK). Both inputs must be unit length.
    Quat slerp(const Quat &q, T t) const noexcept {
        T cosOmega = dot(q);

        // q and -q are the same rotation; pick the one on the near side of the 4D sphere.
        Quat q2 = q;
        if (cosOmega < T(0)) {
            q2 = -q;
            cosOmega = -cosOmega;
        }

        // Use lerp for very small angles to avoid division by zero
        const T slerpEpsilon = T(0.9995);
        if (cosOmega > slerpEpsilon) {
            return lerp(q2, t);
        }

        T omega = std::acos(cosOmega);
        T sinOmega = std::sin(omega);
        T s0 = std::sin((T(1) - t) * omega) / sinOmega;
        T s1 = std::sin(t * omega) / sinOmega;

        return Quat {
            x * s0 + q2.x * s1,
            y * s0 + q2.y * s1,
            z * s0 + q2.z * s1,
            w * s0 + q2.w * s1
        };
    }

    /// Rotation matrix of this unit quaternion (column vectors, right-handed).
    Mat4<T> toMatrix() const noexcept {
        T xx = x * x;
        T yy = y * y;
        T zz = z * z;
        T xy = x * y;
        T xz = x * z;
        T yz = y * z;
        T wx = w * x;
        T wy = w * y;
        T wz = w * z;

        return Mat4<T> {
            T(1) - T(2) * (yy + zz), T(2) * (xy - wz), T(2) * (xz + wy), T(0),
            T(2) * (xy + wz), T(1) - T(2) * (xx + zz), T(2) * (yz - wx), T(0),
            T(2) * (xz - wy), T(2) * (yz + wx), T(1) - T(2) * (xx + yy), T(0),
            T(0), T(0), T(0), T(1)
        };
    }

    /// Extracts axis and angle (radians, 0..2pi); near zero rotation the axis is an arbitrary +X.
    void toAxisAngle(Vec3<T> &axis, T &angle) const noexcept {
        Quat q = *this;
        if (std::abs(q.w) > T(1)) {
            q.normalize();
        }

        angle = T(2) * std::acos(q.w);
        T s = std::sqrt(T(1) - q.w * q.w);

        if (s < T(0.001)) {
            // Near-zero rotation: the axis is numerically meaningless, +X is as good as any.
            axis = Vec3<T> { 1, 0, 0 };
        } else {
            axis = Vec3<T> { q.x / s, q.y / s, q.z / s };
        }
    }

    /// Euler angles in radians as (roll about X, pitch about Y, yaw about Z); order differs from fromEuler().
    Vec3<T> toEuler() const noexcept {
        Vec3<T> euler;

        // Roll (x-axis rotation)
        T sinr_cosp = T(2) * (w * x + y * z);
        T cosr_cosp = T(1) - T(2) * (x * x + y * y);
        euler.x = std::atan2(sinr_cosp, cosr_cosp);

        // Pitch (y-axis rotation); clamp at gimbal lock where asin would get an out-of-range argument.
        T sinp = T(2) * (w * y - z * x);
        if (std::abs(sinp) >= T(1)) {
            euler.y = std::copysign(std::numbers::pi_v<T> / T(2), sinp);
        } else {
            euler.y = std::asin(sinp);
        }

        // Yaw (z-axis rotation)
        T siny_cosp = T(2) * (w * z + x * y);
        T cosy_cosp = T(1) - T(2) * (y * y + z * z);
        euler.z = std::atan2(siny_cosp, cosy_cosp);

        return euler;
    }

    /// Rotation angle in radians, 0..pi (sign of w is ignored).
    T getAngle() const noexcept {
        return T(2) * std::acos(std::min(std::abs(w), T(1)));
    }

    /// Unit rotation axis; +X for a near-identity rotation where the axis is undefined.
    Vec3<T> getAxis() const noexcept {
        T s = std::sqrt(T(1) - w * w);
        if (s < T(0.001)) {
            return Vec3<T> { 1, 0, 0 };
        }
        return Vec3<T> { x / s, y / s, z / s };
    }

    /// Identity rotation.
    static constexpr Quat identity() noexcept {
        return Quat { 0, 0, 0, 1 };
    }

    /// Rotation looking along @p forward with @p up as the vertical hint; @p forward and
    /// @p up must not be parallel (NaN otherwise).
    static Quat lookRotation(const Vec3<T> &forward, const Vec3<T> &up = Vec3<T> { 0, 1, 0 }) noexcept {
        Vec3<T> f = forward.normalized();
        Vec3<T> r = up.cross(f).normalized();
        Vec3<T> u = f.cross(r);

        Mat4<T> m = Mat4<T>::identity();
        // Column-major: column 0 = right, column 1 = up, column 2 = forward.  This used to
        // be stored transposed and only came out right because fromMatrix() read it
        // transposed too.
        m[0][0] = r.x;
        m[0][1] = r.y;
        m[0][2] = r.z;
        m[1][0] = u.x;
        m[1][1] = u.y;
        m[1][2] = u.z;
        m[2][0] = f.x;
        m[2][1] = f.y;
        m[2][2] = f.z;

        return fromMatrix(m);
    }

    /// Shortest-arc rotation taking direction @p from onto direction @p to (both normalized here; zero input gives NaN).
    static Quat fromToRotation(const Vec3<T> &from, const Vec3<T> &to) noexcept {
        Vec3<T> f = from.normalized();
        Vec3<T> t = to.normalized();

        T dot = f.dot(t);

        // Already aligned: avoid the near-zero denominator below.
        if (dot >= T(0.999999)) {
            return identity();
        }

        // Opposite: any perpendicular axis works; try X first, fall back to Y if f is nearly along X.
        if (dot <= T(-0.999999)) {
            Vec3<T> axis = Vec3<T> { 1, 0, 0 }.cross(f);
            if (axis.lengthSquared() < T(0.01)) {
                axis = Vec3<T> { 0, 1, 0 }.cross(f);
            }
            return fromAxisAngle(axis.normalized(), std::numbers::pi_v<T>);
        }

        // Half-angle trick: (cross, 1 + dot) normalized is the shortest-arc quaternion.
        Vec3<T> axis = f.cross(t);
        T       s = std::sqrt((T(1) + dot) * T(2));
        T       invS = T(1) / s;

        return Quat {
            axis.x * invS,
            axis.y * invS,
            axis.z * invS,
            s * T(0.5)
        }
            .normalized();
    }
};


/// Scalar times quaternion; commutes with Quat::operator*(T).
template <typename T>
constexpr Quat<T> operator*(T scalar, const Quat<T> &q) noexcept {
    return q * scalar;
}


using Quatf = Quat<float>;
using Quatd = Quat<double>;


// Layout is part of the contract: 4 tightly packed T, uploadable as-is.
static_assert(std::is_standard_layout<Quatf>::value, "Quatf must be standard layout");
static_assert(sizeof(Quatf) == sizeof(float) * 4, "Quatf must be tightly packed");

} // namespace decimate::math
