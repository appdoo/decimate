// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <decimate/math/vec3.h>
#include <decimate/math/vec4.h>
#include <cmath>
#include <type_traits>

namespace decimate::math {

/// 4x4 matrix, column-major (OpenGL/Vulkan upload layout): element (col, row) is `m[col * 4 + row]`.
///
/// Vectors are columns: `M * v` applies M to v, and `A * B` applies B first. Standard layout,
/// tightly packed, so `data()` can be handed to the GPU as 16 contiguous T.
/// The multi-argument constructor takes its 16 arguments in ROW-major reading order
/// (m00 m01 ... m33) and stores them transposed.
template <typename T>
struct Mat4 {
    static_assert(std::is_arithmetic<T>::value, "Mat4 requires arithmetic type");

    // Column 0: m[0], m[1], m[2], m[3]
    // Column 1: m[4], m[5], m[6], m[7]
    // Column 2: m[8], m[9], m[10], m[11]
    // Column 3: m[12], m[13], m[14], m[15]
    /// Column-major storage, see the class comment.
    T m[16];

    /// Identity matrix.
    constexpr Mat4() noexcept
        : m {
              1, 0, 0, 0, // Column 0
              0, 1, 0, 0, // Column 1
              0, 0, 1, 0, // Column 2
              0, 0, 0, 1  // Column 3 (identity matrix)
          } {
    }

    /// @p diagonal on the main diagonal, zero elsewhere (Mat4(1) is identity, Mat4(0) is zero).
    constexpr explicit Mat4(T diagonal) noexcept
        : m {
              diagonal, 0, 0, 0,
              0, diagonal, 0, 0,
              0, 0, diagonal, 0,
              0, 0, 0, diagonal
          } {
    }

    /// Matrix from 16 values in row-major reading order (mRC = row R, column C).
    constexpr Mat4(
        T m00, T m01, T m02, T m03,
        T m10, T m11, T m12, T m13,
        T m20, T m21, T m22, T m23,
        T m30, T m31, T m32, T m33) noexcept
        : m {
              m00, m10, m20, m30, // Column 0 (note: row-major params, column-major storage)
              m01, m11, m21, m31, // Column 1
              m02, m12, m22, m32, // Column 2
              m03, m13, m23, m33  // Column 3
          } {
    }

    constexpr Mat4(const Mat4 &) noexcept = default;
    constexpr Mat4(Mat4 &&) noexcept = default;
    constexpr Mat4 &operator=(const Mat4 &) noexcept = default;
    constexpr Mat4 &operator=(Mat4 &&) noexcept = default;

    /// Pointer to the first element of column @p col (0..3); no bounds check.
    constexpr T *operator[](size_t col) noexcept {
        return &m[col * 4];
    }
    constexpr const T *operator[](size_t col) const noexcept {
        return &m[col * 4];
    }

    /// The 16 contiguous elements, column-major, for GPU upload.
    constexpr T *data() noexcept {
        return m;
    }
    constexpr const T *data() const noexcept {
        return m;
    }
    constexpr operator T *() noexcept {
        return m;
    }
    constexpr operator const T *() const noexcept {
        return m;
    }

    /// Matrix product; the result applies @p other first, then this.
    Mat4 operator*(const Mat4 &other) const noexcept {
        Mat4 result;
        for (size_t col = 0; col < 4; ++col) {
            for (size_t row = 0; row < 4; ++row) {
                T sum = 0;
                for (size_t k = 0; k < 4; ++k) {
                    sum += m[k * 4 + row] * other.m[col * 4 + k];
                }
                result.m[col * 4 + row] = sum;
            }
        }
        return result;
    }

    Mat4 &operator*=(const Mat4 &other) noexcept {
        *this = *this * other;
        return *this;
    }

    /// Transforms @p v as a column vector.
    Vec4<T> operator*(const Vec4<T> &v) const noexcept {
        return {
            m[0] * v.x + m[4] * v.y + m[8] * v.z + m[12] * v.w,
            m[1] * v.x + m[5] * v.y + m[9] * v.z + m[13] * v.w,
            m[2] * v.x + m[6] * v.y + m[10] * v.z + m[14] * v.w,
            m[3] * v.x + m[7] * v.y + m[11] * v.z + m[15] * v.w
        };
    }

    /// Transforms a position (w = 1) and divides by the resulting w; w == 0 yields inf/NaN.
    Vec3<T> transformPoint(const Vec3<T> &v) const noexcept {
        T w = m[3] * v.x + m[7] * v.y + m[11] * v.z + m[15];
        return {
            (m[0] * v.x + m[4] * v.y + m[8] * v.z + m[12]) / w,
            (m[1] * v.x + m[5] * v.y + m[9] * v.z + m[13]) / w,
            (m[2] * v.x + m[6] * v.y + m[10] * v.z + m[14]) / w
        };
    }

    /// Transforms a direction (w = 0): translation is ignored and no perspective divide happens.
    Vec3<T> transformDirection(const Vec3<T> &v) const noexcept {
        return {
            m[0] * v.x + m[4] * v.y + m[8] * v.z,
            m[1] * v.x + m[5] * v.y + m[9] * v.z,
            m[2] * v.x + m[6] * v.y + m[10] * v.z
        };
    }

    /// Every element times @p scalar.
    constexpr Mat4 operator*(T scalar) const noexcept {
        Mat4 result;
        for (size_t i = 0; i < 16; ++i) {
            result.m[i] = m[i] * scalar;
        }
        return result;
    }

    constexpr Mat4 &operator*=(T scalar) noexcept {
        for (size_t i = 0; i < 16; ++i) {
            m[i] *= scalar;
        }
        return *this;
    }

    /// Elementwise sum.
    constexpr Mat4 operator+(const Mat4 &other) const noexcept {
        Mat4 result;
        for (size_t i = 0; i < 16; ++i) {
            result.m[i] = m[i] + other.m[i];
        }
        return result;
    }

    /// Elementwise difference.
    constexpr Mat4 operator-(const Mat4 &other) const noexcept {
        Mat4 result;
        for (size_t i = 0; i < 16; ++i) {
            result.m[i] = m[i] - other.m[i];
        }
        return result;
    }

    /// Exact elementwise comparison (no tolerance).
    constexpr bool operator==(const Mat4 &other) const noexcept {
        for (size_t i = 0; i < 16; ++i) {
            if (m[i] != other.m[i])
                return false;
        }
        return true;
    }

    constexpr bool operator!=(const Mat4 &other) const noexcept {
        return !(*this == other);
    }

    /// Copy with rows and columns swapped.
    constexpr Mat4 transposed() const noexcept {
        return Mat4 {
            m[0], m[4], m[8], m[12],
            m[1], m[5], m[9], m[13],
            m[2], m[6], m[10], m[14],
            m[3], m[7], m[11], m[15]
        };
    }

    /// Transposes in place.
    void transpose() noexcept {
        *this = transposed();
    }

    /// Determinant by Laplace expansion along the first row.
    T determinant() const noexcept {
        T a = m[0], b = m[4], c = m[8], d = m[12];
        T e = m[1], f = m[5], g = m[9], h = m[13];
        T i = m[2], j = m[6], k = m[10], l = m[14];
        T n = m[3], o = m[7], p = m[11], q = m[15];

        T det3_1 = f * (k * q - l * p) - j * (g * q - h * p) + o * (g * l - h * k);
        T det3_2 = e * (k * q - l * p) - i * (g * q - h * p) + n * (g * l - h * k);
        T det3_3 = e * (j * q - l * o) - i * (f * q - h * o) + n * (f * l - h * j);
        T det3_4 = e * (j * p - k * o) - i * (f * p - g * o) + n * (f * k - g * j);

        return a * det3_1 - b * det3_2 + c * det3_3 - d * det3_4;
    }

    /// Inverse by the adjugate method.
    ///
    /// @returns the inverse, or IDENTITY if |det| < 1e-10 (absolute threshold, so a valid matrix
    ///          with tiny scale is also reported as singular). The caller cannot tell the two apart.
    Mat4 inverse() const noexcept {
        T a = m[0], b = m[4], c = m[8], d = m[12];
        T e = m[1], f = m[5], g = m[9], h = m[13];
        T i = m[2], j = m[6], k = m[10], l = m[14];
        T n = m[3], o = m[7], p = m[11], q = m[15];

        // First-row cofactors give the determinant; the rest are only needed if it is usable.
        T c00 = f * (k * q - l * p) - j * (g * q - h * p) + o * (g * l - h * k);
        T c01 = -(e * (k * q - l * p) - i * (g * q - h * p) + n * (g * l - h * k));
        T c02 = e * (j * q - l * o) - i * (f * q - h * o) + n * (f * l - h * j);
        T c03 = -(e * (j * p - k * o) - i * (f * p - g * o) + n * (f * k - g * j));

        T det = a * c00 + b * c01 + c * c02 + d * c03;
        if (std::abs(det) < T(1e-10)) {
            return Mat4 {}; // Return identity if singular
        }

        T invDet = T(1) / det;

        T c10 = -(b * (k * q - l * p) - j * (c * q - d * p) + o * (c * l - d * k));
        T c11 = a * (k * q - l * p) - i * (c * q - d * p) + n * (c * l - d * k);
        T c12 = -(a * (j * q - l * o) - i * (b * q - d * o) + n * (b * l - d * j));
        T c13 = a * (j * p - k * o) - i * (b * p - c * o) + n * (b * k - c * j);

        T c20 = b * (g * q - h * p) - f * (c * q - d * p) + o * (c * h - d * g);
        T c21 = -(a * (g * q - h * p) - e * (c * q - d * p) + n * (c * h - d * g));
        T c22 = a * (f * q - h * o) - e * (b * q - d * o) + n * (b * h - d * f);
        T c23 = -(a * (f * p - g * o) - e * (b * p - c * o) + n * (b * g - c * f));

        T c30 = -(b * (g * l - h * k) - f * (c * l - d * k) + j * (c * h - d * g));
        T c31 = a * (g * l - h * k) - e * (c * l - d * k) + i * (c * h - d * g);
        T c32 = -(a * (f * l - h * j) - e * (b * l - d * j) + i * (b * h - d * f));
        T c33 = a * (f * k - g * j) - e * (b * k - c * j) + i * (b * g - c * f);

        return Mat4 {
            c00 * invDet, c01 * invDet, c02 * invDet, c03 * invDet,
            c10 * invDet, c11 * invDet, c12 * invDet, c13 * invDet,
            c20 * invDet, c21 * invDet, c22 * invDet, c23 * invDet,
            c30 * invDet, c31 * invDet, c32 * invDet, c33 * invDet
        };
    }

    /// Identity matrix.
    static constexpr Mat4 identity() noexcept {
        return Mat4 {};
    }

    /// All-zero matrix.
    static constexpr Mat4 zero() noexcept {
        return Mat4 { T(0) };
    }

    /// Translation by @p v.
    static constexpr Mat4 translate(const Vec3<T> &v) noexcept {
        return Mat4 {
            1, 0, 0, v.x,
            0, 1, 0, v.y,
            0, 0, 1, v.z,
            0, 0, 0, 1
        };
    }

    /// Per-axis scale by @p v.
    static constexpr Mat4 scale(const Vec3<T> &v) noexcept {
        return Mat4 {
            v.x, 0, 0, 0,
            0, v.y, 0, 0,
            0, 0, v.z, 0,
            0, 0, 0, 1
        };
    }

    /// Uniform scale by @p s.
    static constexpr Mat4 scale(T s) noexcept {
        return scale(Vec3<T> { s, s, s });
    }

    /// Right-handed rotation by @p radians around the X axis.
    static Mat4 rotateX(T radians) noexcept {
        T c = std::cos(radians);
        T s = std::sin(radians);
        return Mat4 {
            1, 0, 0, 0,
            0, c, -s, 0,
            0, s, c, 0,
            0, 0, 0, 1
        };
    }

    /// Right-handed rotation by @p radians around the Y axis.
    static Mat4 rotateY(T radians) noexcept {
        T c = std::cos(radians);
        T s = std::sin(radians);
        return Mat4 {
            c, 0, s, 0,
            0, 1, 0, 0,
            -s, 0, c, 0,
            0, 0, 0, 1
        };
    }

    /// Right-handed rotation by @p radians around the Z axis.
    static Mat4 rotateZ(T radians) noexcept {
        T c = std::cos(radians);
        T s = std::sin(radians);
        return Mat4 {
            c, -s, 0, 0,
            s, c, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        };
    }

    /// Right-handed rotation by @p radians around @p axis (normalized here; a zero axis yields NaN).
    static Mat4 rotate(T radians, const Vec3<T> &axis) noexcept {
        Vec3<T> a = axis.normalized();
        T       c = std::cos(radians);
        T       s = std::sin(radians);
        T       t = T(1) - c;

        return Mat4 {
            t * a.x * a.x + c, t * a.x * a.y - s * a.z, t * a.x * a.z + s * a.y, 0,
            t * a.x * a.y + s * a.z, t * a.y * a.y + c, t * a.y * a.z - s * a.x, 0,
            t * a.x * a.z - s * a.y, t * a.y * a.z + s * a.x, t * a.z * a.z + c, 0,
            0, 0, 0, 1
        };
    }

    /// Right-handed view matrix looking from @p eye toward @p center; @p up must not be parallel to the view direction (NaN otherwise).
    static Mat4 lookAt(const Vec3<T> &eye, const Vec3<T> &center, const Vec3<T> &up) noexcept {
        Vec3<T> f = (center - eye).normalized();
        Vec3<T> s = f.cross(up).normalized();
        Vec3<T> u = s.cross(f);

        return Mat4 {
            s.x, s.y, s.z, -s.dot(eye),
            u.x, u.y, u.z, -u.dot(eye),
            -f.x, -f.y, -f.z, f.dot(eye),
            0, 0, 0, 1
        };
    }

    /// Right-handed perspective projection; @p fovY is the vertical field of view in radians,
    /// @p near and @p far are positive distances. Y is flipped for Vulkan clip space (Y down);
    /// depth is mapped to [-1, 1] (OpenGL convention).
    static Mat4 perspective(T fovY, T aspect, T near, T far) noexcept {
        T tanHalfFovy = std::tan(fovY / T(2));
        T range = far - near;

        // Y flip: Vulkan's clip space has Y pointing down.
        return Mat4 {
            T(1) / (aspect * tanHalfFovy), 0, 0, 0,
            0, -T(1) / tanHalfFovy, 0, 0, // Negated for Vulkan
            0, 0, -(far + near) / range, -(T(2) * far * near) / range,
            0, 0, -1, 0
        };
    }

    /// Right-handed orthographic projection of the box [left,right] x [bottom,top] x [-near,-far]; depth maps to [-1, 1].
    static constexpr Mat4 orthographic(T left, T right, T bottom, T top, T near, T far) noexcept {
        T rl = right - left;
        T tb = top - bottom;
        T fn = far - near;

        return Mat4 {
            T(2) / rl, 0, 0, -(right + left) / rl,
            0, T(2) / tb, 0, -(top + bottom) / tb,
            0, 0, -T(2) / fn, -(far + near) / fn,
            0, 0, 0, 1
        };
    }
};


/// Scalar times matrix; commutes with Mat4::operator*(T).
template <typename T>
constexpr Mat4<T> operator*(T scalar, const Mat4<T> &mat) noexcept {
    return mat * scalar;
}


using Mat4f = Mat4<float>;
using Mat4d = Mat4<double>;


// Uploaded to the GPU as-is, so layout is part of the contract.
static_assert(std::is_standard_layout<Mat4f>::value, "Mat4f must be standard layout");
static_assert(sizeof(Mat4f) == sizeof(float) * 16, "Mat4f must be tightly packed");

} // namespace decimate::math
