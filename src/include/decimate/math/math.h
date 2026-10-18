/**
 * pure-math — a small, fast, header-only C++ vector/matrix math library
 *
 * version 0.6
 *
 * Copyright (C) 2003-2026 Andor Greißl (Greissl)
 * Report bugs and download new versions at https://github.com/appdoo/decimate
 *
 * This library is distributed under the MIT License. See notice at the end of
 * this file.
 *
 * This is decimate's foundational vector-math layer (namespace `decimate::math`),
 * absorbed from the formerly standalone pure-math library. The other core math
 * headers carry a one-line reference back to this header for version/license
 * details.
 */
#pragma once

/// @file PureMath.h
/// @brief Central include for all pure-math types + the canonical version/license
///        notice for the library.
///
/// Usage:
///   #include <decimate/math/math.h>
///   using namespace decimate::math;

// Core math primitives
#include <decimate/math/vec2.h>
#include <decimate/math/vec3.h>
#include <decimate/math/vec4.h>
#include <decimate/math/mat4.h>
#include <decimate/math/quat.h>

// Geometry primitives
#include <decimate/math/plane.h>
#include <decimate/math/aabb.h>

namespace decimate::math {

// Common mathematical constants
namespace constants {
    constexpr float PI = 3.14159265358979323846f;
    constexpr float TWO_PI = 2.0f * PI;
    constexpr float HALF_PI = 0.5f * PI;
    constexpr float DEG_TO_RAD = PI / 180.0f;
    constexpr float RAD_TO_DEG = 180.0f / PI;
    constexpr float EPSILON = 1e-6f;
}

// Utility functions (directly in namespace decimate::math)
    template<typename T>
    inline constexpr T radians(T degrees) noexcept {
        return degrees * static_cast<T>(constants::DEG_TO_RAD);
    }

    template<typename T>
    inline constexpr T degrees(T radians) noexcept {
        return radians * static_cast<T>(constants::RAD_TO_DEG);
    }

    template<typename T>
    inline constexpr T clamp(T value, T min, T max) noexcept {
        return value < min ? min : (value > max ? max : value);
    }

    template<typename T>
    inline constexpr T lerp(T a, T b, T t) noexcept {
        return a + t * (b - a);
    }

} // namespace decimate::math

/*
 * MIT License
 *
 * Copyright (C) 2003-2026 Andor Greißl (Greissl)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
