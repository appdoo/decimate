// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <cmath>
#include <type_traits>

namespace decimate::math {

/**
 * Vec2 - 2D vector template class
 *
 * Optimized for texture coordinates and 2D operations
 */
template<typename T>
struct Vec2 {
    static_assert(std::is_arithmetic<T>::value, "Vec2 requires arithmetic type");

    // Dual access: name the components (.x/.y) or index them (data[i]) — both
    // views alias the same two contiguous T (standard-layout; see asserts
    // below). Matches the Vec3/Vec4 union layout for family consistency.
    union {
        T data[2];
        struct { T x, y; };
    };

    // Constructors — the default zero-initializes.
    constexpr Vec2() noexcept : x(0), y(0) {}
    constexpr Vec2(T x, T y) noexcept : x(x), y(y) {}
    constexpr explicit Vec2(T scalar) noexcept : x(scalar), y(scalar) {}

    // Cross-precision conversion (e.g. Vec2<double> → Vec2<float>); explicit
    // so a change of element type is never silent.
    template<typename U>
    constexpr explicit Vec2(const Vec2<U>& o) noexcept
        : x(static_cast<T>(o.x)), y(static_cast<T>(o.y)) {}

    // Copy and move constructors
    constexpr Vec2(const Vec2&) noexcept = default;
    constexpr Vec2(Vec2&&) noexcept = default;
    constexpr Vec2& operator=(const Vec2&) noexcept = default;
    constexpr Vec2& operator=(Vec2&&) noexcept = default;

    // Array-like access via the union array. NOT constexpr: reading the
    // inactive union member is ill-formed in a constant expression.
    T& operator[](size_t i) noexcept { return data[i]; }
    const T& operator[](size_t i) const noexcept { return data[i]; }

    // Direct pointer access for GPU upload / raw interop (the `data` member is
    // equally accessible directly).
    operator T*() noexcept { return data; }
    operator const T*() const noexcept { return data; }

    // Arithmetic operators
    constexpr Vec2 operator+(const Vec2& v) const noexcept {
        return {x + v.x, y + v.y};
    }

    constexpr Vec2 operator-(const Vec2& v) const noexcept {
        return {x - v.x, y - v.y};
    }

    constexpr Vec2 operator*(T scalar) const noexcept {
        return {x * scalar, y * scalar};
    }

    constexpr Vec2 operator/(T scalar) const noexcept {
        return {x / scalar, y / scalar};
    }

    constexpr Vec2 operator-() const noexcept {
        return {-x, -y};
    }

    // Compound assignment operators
    constexpr Vec2& operator+=(const Vec2& v) noexcept {
        x += v.x; y += v.y;
        return *this;
    }

    constexpr Vec2& operator-=(const Vec2& v) noexcept {
        x -= v.x; y -= v.y;
        return *this;
    }

    constexpr Vec2& operator*=(T scalar) noexcept {
        x *= scalar; y *= scalar;
        return *this;
    }

    constexpr Vec2& operator/=(T scalar) noexcept {
        x /= scalar; y /= scalar;
        return *this;
    }

    // Comparison operators
    constexpr bool operator==(const Vec2& v) const noexcept {
        return x == v.x && y == v.y;
    }

    constexpr bool operator!=(const Vec2& v) const noexcept {
        return !(*this == v);
    }

    // Vector operations
    constexpr T dot(const Vec2& v) const noexcept {
        return x * v.x + y * v.y;
    }

    T length() const noexcept {
        return std::sqrt(x * x + y * y);
    }

    constexpr T lengthSquared() const noexcept {
        return x * x + y * y;
    }

    Vec2 normalized() const noexcept {
        T len = length();
        return len > T(0) ? (*this / len) : Vec2{};
    }

    void normalize() noexcept {
        T len = length();
        if (len > T(0)) {
            *this /= len;
        }
    }

    // Static factory methods
    static constexpr Vec2 zero() noexcept { return {0, 0}; }
    static constexpr Vec2 one() noexcept { return {1, 1}; }
    static constexpr Vec2 unitX() noexcept { return {1, 0}; }
    static constexpr Vec2 unitY() noexcept { return {0, 1}; }
};

// Non-member operators
template<typename T>
constexpr Vec2<T> operator*(T scalar, const Vec2<T>& v) noexcept {
    return v * scalar;
}

// Common typedefs
using Vec2f = Vec2<float>;
using Vec2d = Vec2<double>;
using Vec2i = Vec2<int>;
using Vec2u = Vec2<unsigned int>;

// Ensure standard layout for GPU compatibility
static_assert(std::is_standard_layout<Vec2f>::value, "Vec2f must be standard layout");
static_assert(sizeof(Vec2f) == sizeof(float) * 2, "Vec2f must be tightly packed");

} // namespace decimate::math
