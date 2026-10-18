// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <cmath>
#include <type_traits>

namespace decimate::math {

/**
 * Vec3 - 3D vector template class
 *
 * Designed for high performance with:
 * - Contiguous memory layout for GPU upload (castable to T*)
 * - Move semantics and constexpr where possible
 * - Zero-overhead abstractions
 * - SSE/NEON friendly alignment
 */
template<typename T>
struct Vec3 {
    static_assert(std::is_arithmetic<T>::value, "Vec3 requires arithmetic type");

    // Dual access: name the components (.x/.y/.z) or index them (data[i]) —
    // both views alias the same three contiguous T (standard-layout; see the
    // asserts at the bottom of this file). The anonymous struct inside the
    // union is a near-universal compiler extension (glm/Unreal rely on it) and
    // compiles clean on the project's clang/gcc/Emscripten toolchains.
    union {
        T data[3];
        struct { T x, y, z; };
    };

    // Constructors — the default zero-initializes; callers depend on it (e.g.
    // a std::vector<Vec3> normal accumulator that is then +='d into).
    constexpr Vec3() noexcept : x(0), y(0), z(0) {}
    constexpr Vec3(T x, T y, T z) noexcept : x(x), y(y), z(z) {}
    constexpr explicit Vec3(T scalar) noexcept : x(scalar), y(scalar), z(scalar) {}

    // Cross-precision conversion (e.g. Vec3<double> → Vec3<float>); explicit
    // so a change of element type is never silent.
    template<typename U>
    constexpr explicit Vec3(const Vec3<U>& o) noexcept
        : x(static_cast<T>(o.x)), y(static_cast<T>(o.y)), z(static_cast<T>(o.z)) {}

    // Copy and move constructors (trivial for arithmetic T)
    constexpr Vec3(const Vec3&) noexcept = default;
    constexpr Vec3(Vec3&&) noexcept = default;
    constexpr Vec3& operator=(const Vec3&) noexcept = default;
    constexpr Vec3& operator=(Vec3&&) noexcept = default;

    // Array-like access via the union array. NOT constexpr: reading the
    // inactive union member is ill-formed in a constant expression.
    T& operator[](size_t i) noexcept { return data[i]; }
    const T& operator[](size_t i) const noexcept { return data[i]; }

    // Direct pointer access for GPU upload / raw interop (the `data` member is
    // equally accessible directly).
    operator T*() noexcept { return data; }
    operator const T*() const noexcept { return data; }

    // Arithmetic operators
    constexpr Vec3 operator+(const Vec3& v) const noexcept {
        return {x + v.x, y + v.y, z + v.z};
    }

    constexpr Vec3 operator-(const Vec3& v) const noexcept {
        return {x - v.x, y - v.y, z - v.z};
    }

    constexpr Vec3 operator*(T scalar) const noexcept {
        return {x * scalar, y * scalar, z * scalar};
    }

    constexpr Vec3 operator/(T scalar) const noexcept {
        return {x / scalar, y / scalar, z / scalar};
    }

    constexpr Vec3 operator-() const noexcept {
        return {-x, -y, -z};
    }

    // Compound assignment operators
    constexpr Vec3& operator+=(const Vec3& v) noexcept {
        x += v.x; y += v.y; z += v.z;
        return *this;
    }

    constexpr Vec3& operator-=(const Vec3& v) noexcept {
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }

    constexpr Vec3& operator*=(T scalar) noexcept {
        x *= scalar; y *= scalar; z *= scalar;
        return *this;
    }

    constexpr Vec3& operator/=(T scalar) noexcept {
        x /= scalar; y /= scalar; z /= scalar;
        return *this;
    }

    // Comparison operators
    constexpr bool operator==(const Vec3& v) const noexcept {
        return x == v.x && y == v.y && z == v.z;
    }

    constexpr bool operator!=(const Vec3& v) const noexcept {
        return !(*this == v);
    }

    // Vector operations
    constexpr T dot(const Vec3& v) const noexcept {
        return x * v.x + y * v.y + z * v.z;
    }

    constexpr Vec3 cross(const Vec3& v) const noexcept {
        return {
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        };
    }

    // Component-wise multiplication
    constexpr Vec3 operator*(const Vec3& v) const noexcept {
        return {x * v.x, y * v.y, z * v.z};
    }

    constexpr Vec3& operator*=(const Vec3& v) noexcept {
        x *= v.x; y *= v.y; z *= v.z;
        return *this;
    }

    T length() const noexcept {
        return std::sqrt(x * x + y * y + z * z);
    }

    constexpr T lengthSquared() const noexcept {
        return x * x + y * y + z * z;
    }

    T distance(const Vec3& v) const noexcept {
        return (*this - v).length();
    }

    constexpr T distanceSquared(const Vec3& v) const noexcept {
        return (*this - v).lengthSquared();
    }

    Vec3 normalized() const noexcept {
        T len = length();
        return len > T(0) ? (*this / len) : Vec3{};
    }

    void normalize() noexcept {
        T len = length();
        if (len > T(0)) {
            *this /= len;
        }
    }

    // Linear interpolation
    constexpr Vec3 lerp(const Vec3& v, T t) const noexcept {
        return *this + (v - *this) * t;
    }

    // Reflection: reflect this vector across a normal
    Vec3 reflect(const Vec3& normal) const noexcept {
        return *this - normal * (T(2) * dot(normal));
    }

    // Clamp each component between min and max
    constexpr Vec3 clamp(T minVal, T maxVal) const noexcept {
        return {
            x < minVal ? minVal : (x > maxVal ? maxVal : x),
            y < minVal ? minVal : (y > maxVal ? maxVal : y),
            z < minVal ? minVal : (z > maxVal ? maxVal : z)
        };
    }

    // Static factory methods
    static constexpr Vec3 zero() noexcept { return {0, 0, 0}; }
    static constexpr Vec3 one() noexcept { return {1, 1, 1}; }
    static constexpr Vec3 unitX() noexcept { return {1, 0, 0}; }
    static constexpr Vec3 unitY() noexcept { return {0, 1, 0}; }
    static constexpr Vec3 unitZ() noexcept { return {0, 0, 1}; }

    // Min/max component-wise
    static constexpr Vec3 min(const Vec3& a, const Vec3& b) noexcept {
        return {
            a.x < b.x ? a.x : b.x,
            a.y < b.y ? a.y : b.y,
            a.z < b.z ? a.z : b.z
        };
    }

    static constexpr Vec3 max(const Vec3& a, const Vec3& b) noexcept {
        return {
            a.x > b.x ? a.x : b.x,
            a.y > b.y ? a.y : b.y,
            a.z > b.z ? a.z : b.z
        };
    }
};

// Non-member operators for scalar * vector
template<typename T>
constexpr Vec3<T> operator*(T scalar, const Vec3<T>& v) noexcept {
    return v * scalar;
}

// Common typedefs
using Vec3f = Vec3<float>;
using Vec3d = Vec3<double>;
using Vec3i = Vec3<int>;
using Vec3u = Vec3<unsigned int>;

// Ensure standard layout for GPU compatibility
static_assert(std::is_standard_layout<Vec3f>::value, "Vec3f must be standard layout");
static_assert(sizeof(Vec3f) == sizeof(float) * 3, "Vec3f must be tightly packed");

} // namespace decimate::math
