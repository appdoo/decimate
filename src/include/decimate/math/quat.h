// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <decimate/math/vec3.h>
#include <decimate/math/mat4.h>
#include <cmath>
#include <numbers>
#include <algorithm>
#include <type_traits>

namespace decimate::math {

/**
 * Quat - Quaternion template class for rotations
 *
 * Stores rotation as (x, y, z, w) where w is the scalar component.
 * Designed for smooth rotation interpolation (slerp) and avoiding gimbal lock.
 *
 * Based on:
 * - Graphics Gems quaternion articles
 * - DirectXTK implementation
 * - Modern game-math library best practices
 *
 * Designed for high performance with:
 * - Contiguous memory layout (castable to T*)
 * - Move semantics and constexpr where possible
 * - SIMD-friendly alignment (16 bytes for float)
 * - Normalized quaternions for rotations
 */
template<typename T>
struct Quat {
    static_assert(std::is_floating_point<T>::value, "Quat requires floating point type");

    T x, y, z, w;

    // Constructors
    constexpr Quat() noexcept : x(0), y(0), z(0), w(1) {}  // Identity quaternion
    constexpr Quat(T x, T y, T z, T w) noexcept : x(x), y(y), z(z), w(w) {}

    // Construct from axis and angle
    static Quat fromAxisAngle(const Vec3<T>& axis, T angle) noexcept {
        T halfAngle = angle * T(0.5);
        T s = std::sin(halfAngle);
        Vec3<T> normAxis = axis.normalized();
        return Quat{
            normAxis.x * s,
            normAxis.y * s,
            normAxis.z * s,
            std::cos(halfAngle)
        };
    }

    // Construct from Euler angles (in radians: pitch, yaw, roll)
    static Quat fromEuler(T pitch, T yaw, T roll) noexcept {
        T cy = std::cos(yaw * T(0.5));
        T sy = std::sin(yaw * T(0.5));
        T cp = std::cos(pitch * T(0.5));
        T sp = std::sin(pitch * T(0.5));
        T cr = std::cos(roll * T(0.5));
        T sr = std::sin(roll * T(0.5));

        return Quat{
            sr * cp * cy - cr * sp * sy,  // x
            cr * sp * cy + sr * cp * sy,  // y
            cr * cp * sy - sr * sp * cy,  // z
            cr * cp * cy + sr * sp * sy   // w
        };
    }

    // Construct from rotation matrix
    static Quat fromMatrix(const Mat4<T>& m) noexcept {
        T trace = m[0][0] + m[1][1] + m[2][2];
        Quat q;

        if (trace > T(0)) {
            T s = std::sqrt(trace + T(1)) * T(2);
            q.w = T(0.25) * s;
            q.x = (m[2][1] - m[1][2]) / s;
            q.y = (m[0][2] - m[2][0]) / s;
            q.z = (m[1][0] - m[0][1]) / s;
        } else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
            T s = std::sqrt(T(1) + m[0][0] - m[1][1] - m[2][2]) * T(2);
            q.w = (m[2][1] - m[1][2]) / s;
            q.x = T(0.25) * s;
            q.y = (m[0][1] + m[1][0]) / s;
            q.z = (m[0][2] + m[2][0]) / s;
        } else if (m[1][1] > m[2][2]) {
            T s = std::sqrt(T(1) + m[1][1] - m[0][0] - m[2][2]) * T(2);
            q.w = (m[0][2] - m[2][0]) / s;
            q.x = (m[0][1] + m[1][0]) / s;
            q.y = T(0.25) * s;
            q.z = (m[1][2] + m[2][1]) / s;
        } else {
            T s = std::sqrt(T(1) + m[2][2] - m[0][0] - m[1][1]) * T(2);
            q.w = (m[1][0] - m[0][1]) / s;
            q.x = (m[0][2] + m[2][0]) / s;
            q.y = (m[1][2] + m[2][1]) / s;
            q.z = T(0.25) * s;
        }

        return q.normalized();
    }

    // Copy and move constructors
    constexpr Quat(const Quat&) noexcept = default;
    constexpr Quat(Quat&&) noexcept = default;
    constexpr Quat& operator=(const Quat&) noexcept = default;
    constexpr Quat& operator=(Quat&&) noexcept = default;

    // Array-like access
    constexpr T& operator[](size_t i) noexcept {
        return (&x)[i];
    }
    constexpr const T& operator[](size_t i) const noexcept {
        return (&x)[i];
    }

    // Direct pointer access
    constexpr T* data() noexcept { return &x; }
    constexpr const T* data() const noexcept { return &x; }
    constexpr operator T*() noexcept { return &x; }
    constexpr operator const T*() const noexcept { return &x; }

    // Quaternion operations
    constexpr T dot(const Quat& q) const noexcept {
        return x * q.x + y * q.y + z * q.z + w * q.w;
    }

    T length() const noexcept {
        return std::sqrt(x * x + y * y + z * z + w * w);
    }

    constexpr T lengthSquared() const noexcept {
        return x * x + y * y + z * z + w * w;
    }

    Quat normalized() const noexcept {
        T len = length();
        if (len > T(0)) {
            T invLen = T(1) / len;
            return Quat{x * invLen, y * invLen, z * invLen, w * invLen};
        }
        return Quat{};  // Identity if zero length
    }

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

    // Conjugate (inverse for unit quaternions)
    constexpr Quat conjugate() const noexcept {
        return Quat{-x, -y, -z, w};
    }

    // Inverse
    Quat inverse() const noexcept {
        T lenSq = lengthSquared();
        if (lenSq > T(0)) {
            T invLenSq = T(1) / lenSq;
            return Quat{-x * invLenSq, -y * invLenSq, -z * invLenSq, w * invLenSq};
        }
        return Quat{};  // Identity if zero length
    }

    // Quaternion multiplication (composition of rotations)
    Quat operator*(const Quat& q) const noexcept {
        return Quat{
            w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y - x * q.z + y * q.w + z * q.x,
            w * q.z + x * q.y - y * q.x + z * q.w,
            w * q.w - x * q.x - y * q.y - z * q.z
        };
    }

    Quat& operator*=(const Quat& q) noexcept {
        *this = *this * q;
        return *this;
    }

    // Rotate a vector by this quaternion
    Vec3<T> rotate(const Vec3<T>& v) const noexcept {
        // Using the formula: v' = q * v * q^-1
        // Optimized version: v' = v + 2 * cross(q.xyz, cross(q.xyz, v) + q.w * v)
        Vec3<T> qvec(x, y, z);
        Vec3<T> uv = qvec.cross(v);
        Vec3<T> uuv = qvec.cross(uv);
        uv = uv * (T(2) * w);
        uuv = uuv * T(2);
        return v + uv + uuv;
    }

    // Scalar operations
    constexpr Quat operator*(T scalar) const noexcept {
        return Quat{x * scalar, y * scalar, z * scalar, w * scalar};
    }

    constexpr Quat operator+(const Quat& q) const noexcept {
        return Quat{x + q.x, y + q.y, z + q.z, w + q.w};
    }

    constexpr Quat operator-(const Quat& q) const noexcept {
        return Quat{x - q.x, y - q.y, z - q.z, w - q.w};
    }

    constexpr Quat operator-() const noexcept {
        return Quat{-x, -y, -z, -w};
    }

    // Comparison
    constexpr bool operator==(const Quat& q) const noexcept {
        return x == q.x && y == q.y && z == q.z && w == q.w;
    }

    constexpr bool operator!=(const Quat& q) const noexcept {
        return !(*this == q);
    }

    // Linear interpolation (fast but not constant angular velocity)
    constexpr Quat lerp(const Quat& q, T t) const noexcept {
        return (*this * (T(1) - t) + q * t).normalized();
    }

    // Normalized linear interpolation (faster than slerp, good approximation)
    Quat nlerp(const Quat& q, T t) const noexcept {
        return lerp(q, t);
    }

    // Spherical linear interpolation (constant angular velocity)
    // Based on Graphics Gems and DirectXTK implementations
    Quat slerp(const Quat& q, T t) const noexcept {
        T cosOmega = dot(q);

        // Ensure shortest path
        Quat q2 = q;
        if (cosOmega < T(0)) {
            q2 = -q;
            cosOmega = -cosOmega;
        }

        // Use lerp for very small angles to avoid division by zero
        const T SLERP_EPSILON = T(0.9995);
        if (cosOmega > SLERP_EPSILON) {
            return lerp(q2, t);
        }

        // Standard slerp
        T omega = std::acos(cosOmega);
        T sinOmega = std::sin(omega);
        T s0 = std::sin((T(1) - t) * omega) / sinOmega;
        T s1 = std::sin(t * omega) / sinOmega;

        return Quat{
            x * s0 + q2.x * s1,
            y * s0 + q2.y * s1,
            z * s0 + q2.z * s1,
            w * s0 + q2.w * s1
        };
    }

    // Convert to rotation matrix
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

        return Mat4<T>{
            T(1) - T(2) * (yy + zz), T(2) * (xy - wz),         T(2) * (xz + wy),         T(0),
            T(2) * (xy + wz),         T(1) - T(2) * (xx + zz), T(2) * (yz - wx),         T(0),
            T(2) * (xz - wy),         T(2) * (yz + wx),         T(1) - T(2) * (xx + yy), T(0),
            T(0),                     T(0),                     T(0),                     T(1)
        };
    }

    // Convert to axis-angle representation
    void toAxisAngle(Vec3<T>& axis, T& angle) const noexcept {
        Quat q = *this;
        if (std::abs(q.w) > T(1)) {
            q.normalize();
        }

        angle = T(2) * std::acos(q.w);
        T s = std::sqrt(T(1) - q.w * q.w);

        if (s < T(0.001)) {
            // Angle is close to 0, axis doesn't matter
            axis = Vec3<T>{1, 0, 0};
        } else {
            axis = Vec3<T>{q.x / s, q.y / s, q.z / s};
        }
    }

    // Convert to Euler angles (pitch, yaw, roll in radians)
    Vec3<T> toEuler() const noexcept {
        Vec3<T> euler;

        // Roll (x-axis rotation)
        T sinr_cosp = T(2) * (w * x + y * z);
        T cosr_cosp = T(1) - T(2) * (x * x + y * y);
        euler.x = std::atan2(sinr_cosp, cosr_cosp);

        // Pitch (y-axis rotation)
        T sinp = T(2) * (w * y - z * x);
        if (std::abs(sinp) >= T(1)) {
            euler.y = std::copysign(std::numbers::pi_v<T> / T(2), sinp);  // Use 90 degrees if out of range
        } else {
            euler.y = std::asin(sinp);
        }

        // Yaw (z-axis rotation)
        T siny_cosp = T(2) * (w * z + x * y);
        T cosy_cosp = T(1) - T(2) * (y * y + z * z);
        euler.z = std::atan2(siny_cosp, cosy_cosp);

        return euler;
    }

    // Get rotation angle
    T getAngle() const noexcept {
        return T(2) * std::acos(std::min(std::abs(w), T(1)));
    }

    // Get rotation axis (returns unit vector)
    Vec3<T> getAxis() const noexcept {
        T s = std::sqrt(T(1) - w * w);
        if (s < T(0.001)) {
            return Vec3<T>{1, 0, 0};  // Arbitrary axis for near-zero rotation
        }
        return Vec3<T>{x / s, y / s, z / s};
    }

    // Static factory methods
    static constexpr Quat identity() noexcept {
        return Quat{0, 0, 0, 1};
    }

    // Look rotation (creates rotation from forward and up vectors)
    static Quat lookRotation(const Vec3<T>& forward, const Vec3<T>& up = Vec3<T>{0, 1, 0}) noexcept {
        Vec3<T> f = forward.normalized();
        Vec3<T> r = up.cross(f).normalized();
        Vec3<T> u = f.cross(r);

        Mat4<T> m = Mat4<T>::identity();
        m[0][0] = r.x; m[1][0] = r.y; m[2][0] = r.z;
        m[0][1] = u.x; m[1][1] = u.y; m[2][1] = u.z;
        m[0][2] = f.x; m[1][2] = f.y; m[2][2] = f.z;

        return fromMatrix(m);
    }

    // Rotation from one vector to another
    static Quat fromToRotation(const Vec3<T>& from, const Vec3<T>& to) noexcept {
        Vec3<T> f = from.normalized();
        Vec3<T> t = to.normalized();

        T dot = f.dot(t);

        // Vectors are parallel
        if (dot >= T(0.999999)) {
            return identity();
        }

        // Vectors are opposite
        if (dot <= T(-0.999999)) {
            Vec3<T> axis = Vec3<T>{1, 0, 0}.cross(f);
            if (axis.lengthSquared() < T(0.01)) {
                axis = Vec3<T>{0, 1, 0}.cross(f);
            }
            return fromAxisAngle(axis.normalized(), std::numbers::pi_v<T>);
        }

        // Standard case
        Vec3<T> axis = f.cross(t);
        T s = std::sqrt((T(1) + dot) * T(2));
        T invS = T(1) / s;

        return Quat{
            axis.x * invS,
            axis.y * invS,
            axis.z * invS,
            s * T(0.5)
        }.normalized();
    }
};

// Non-member scalar multiplication
template<typename T>
constexpr Quat<T> operator*(T scalar, const Quat<T>& q) noexcept {
    return q * scalar;
}

// Common typedefs
using Quatf = Quat<float>;
using Quatd = Quat<double>;

// Ensure standard layout for GPU compatibility
static_assert(std::is_standard_layout<Quatf>::value, "Quatf must be standard layout");
static_assert(sizeof(Quatf) == sizeof(float) * 4, "Quatf must be tightly packed");

} // namespace decimate::math
