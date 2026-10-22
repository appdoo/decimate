// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <decimate/math/vec3.h>
#include <decimate/math/mat4.h>
#include <decimate/math/plane.h>
#include <vector>
#include <algorithm>

namespace decimate::math {

/// Axis-aligned bounding box for ray tests, plane culling and spatial queries.
///
/// Based on the mudvis implementation (Apache 2.0, Copyright 2003-2019 Andor Greissl)
/// with optimizations from Graphics Gems (Andrew Woo, Pierre Terdiman).
///
/// Invariant: min <= max componentwise, and `mid` / `corners` match them (every mutator goes
/// through initialize()). The caller must pass ordered corners to the two-corner ctor and
/// initialize(); nothing re-sorts them. A default-constructed box is the zero box at the origin.
///
/// Memory layout:
/// - 8 corner points arranged so the bit-inverted index gives the opposite diagonal corner
/// - Near/far point computation using plane sign bits for fast culling
///
/// Corner arrangement:      5 ---- 1
///                         /|     /|
///                        4 ---- 0 |
///                        | 7 ---| 3
///                        |/     |/
///                        6 ---- 2
template <typename T>
struct AABB {
    static_assert(std::is_floating_point<T>::value, "AABB requires floating point type");

    Vec3<T>              min, max, mid;
    std::vector<Vec3<T>> corners; ///< The 8 corner points, laid out as in the diagram above.

    /// Degenerate box at the origin (min == max == mid == 0); corners are allocated but zero.
    AABB() noexcept {
        corners.resize(8);
    }

    /// Box from two opposite corners; the caller guarantees min <= max componentwise.
    AABB(const Vec3<T> &min, const Vec3<T> &max) noexcept
        : min(min), max(max) {
        initialize(min, max);
    }

    /// Tight box around @p vertices; an empty set yields the zero box at the origin.
    AABB(const std::vector<Vec3<T>> &vertices) noexcept {
        if (vertices.empty()) {
            initialize(Vec3<T>::zero(), Vec3<T>::zero());
            return;
        }

        // Seed with the first vertex so the running min/max need no infinity sentinels.
        min = vertices[0];
        max = vertices[0];

        for (const auto &v : vertices) {
            min.x = std::min(min.x, v.x);
            min.y = std::min(min.y, v.y);
            min.z = std::min(min.z, v.z);

            max.x = std::max(max.x, v.x);
            max.y = std::max(max.y, v.y);
            max.z = std::max(max.z, v.z);
        }

        initialize(min, max);
    }

    /// Sets min/max and recomputes `mid` and the 8 corners; the single place that keeps them in sync.
    void initialize(const Vec3<T> &minPt, const Vec3<T> &maxPt) noexcept {
        min = minPt;
        max = maxPt;
        mid = (min + max) * T(0.5);

        // 8 corners arranged for inverted index = opposite diagonal corner
        corners.resize(8);
        corners[0] = Vec3<T>(max.x, max.y, max.z);
        corners[1] = Vec3<T>(max.x, max.y, min.z);
        corners[2] = Vec3<T>(max.x, min.y, max.z);
        corners[3] = Vec3<T>(max.x, min.y, min.z);
        corners[4] = Vec3<T>(min.x, max.y, max.z);
        corners[5] = Vec3<T>(min.x, max.y, min.z);
        corners[6] = Vec3<T>(min.x, min.y, max.z);
        corners[7] = Vec3<T>(min.x, min.y, min.z);
    }

    /// Corner @p index in 0..7 (see diagram); the index is masked with 7, so out-of-range values wrap.
    const Vec3<T> &getCorner(int index) const noexcept {
        return corners[index & 0x7]; // Wrap to 0-7
    }

    // p-vertex / n-vertex fast plane culling.
    // Port of the 2003-2019 mudvis/vermessung farCorner/nearCorner.
    //
    // The corner FARTHEST along the plane normal (the "p-vertex") is selected
    // in O(1) from the plane's precomputed sign bits: for a negative normal
    // component the min corner is farthest, else the max — and the corner array
    // is laid out so (corner index bit set ⇔ min on that axis), exactly the
    // complement of the plane's signbits (signbit set ⇔ normal negative).

    /// Corner farthest along the plane normal (p-vertex).
    const Vec3<T> &farCorner(const Plane<T> &plane) const noexcept {
        return corners[(~plane.signbits) & 0x7];
    }

    /// Corner nearest along the plane normal (n-vertex); inverse of farCorner.
    const Vec3<T> &nearCorner(const Plane<T> &plane) const noexcept {
        return corners[plane.signbits & 0x7];
    }

    /// True if the whole box is strictly on the NEGATIVE side of `plane`
    /// (i.e. even the p-vertex is behind it) — the fast, conservative reject.
    /// Conservative: never reports "outside" for a box that actually crosses
    /// the plane, so it can over-include but never wrongly exclude.
    bool isOutsidePlane(const Plane<T> &plane) const noexcept {
        return plane.signedDistance(farCorner(plane)) < T(0);
    }

    /// Minimum corner.
    Vec3<T> getMin() const noexcept {
        return min;
    }
    /// Maximum corner.
    Vec3<T> getMax() const noexcept {
        return max;
    }
    /// Center point, (min + max) / 2.
    Vec3<T> getMid() const noexcept {
        return mid;
    }

    /// Full side lengths, max - min.
    Vec3<T> getExtents() const noexcept {
        return max - min;
    }
    /// Half side lengths, (max - min) / 2.
    Vec3<T> getHalfExtents() const noexcept {
        return (max - min) * T(0.5);
    }

    /// Total area of the six faces.
    T surfaceArea() const noexcept {
        Vec3<T> extent = max - min;
        return T(2) * (extent.x * extent.y + extent.x * extent.z + extent.y * extent.z);
    }

    /// Enclosed volume.
    T volume() const noexcept {
        Vec3<T> extent = max - min;
        return extent.x * extent.y * extent.z;
    }

    /// Shortest side length.
    T minExtent() const noexcept {
        Vec3<T> extent = max - min;
        return std::min(std::min(extent.x, extent.y), extent.z);
    }

    /// Longest side length.
    T maxExtent() const noexcept {
        Vec3<T> extent = max - min;
        return std::max(std::max(extent.x, extent.y), extent.z);
    }

    /// True if @p point lies inside or on the boundary of the box.
    bool contains(const Vec3<T> &point) const noexcept {
        return point.x >= min.x && point.x <= max.x &&
               point.y >= min.y && point.y <= max.y &&
               point.z >= min.z && point.z <= max.z;
    }

    /// True if the two boxes overlap or touch (separating-axis test on the three axes).
    bool intersects(const AABB &other) const noexcept {
        if (min.x > other.max.x || min.y > other.max.y || min.z > other.max.z)
            return false;

        if (max.x < other.min.x || max.y < other.min.y || max.z < other.min.z)
            return false;

        return true;
    }

    /// Ray-AABB intersection test, after Andrew Woo (Graphics Gems) with Pierre Terdiman's optimizations.
    ///
    /// @param origin  Ray origin.
    /// @param dir     Ray direction; need not be normalized, t is then in units of |dir|.
    /// @param tMin    Optional out: ray parameter of the entry point; 0 when origin is inside.
    /// @param tMax    Optional out: currently the SAME value as tMin (exit point is not computed).
    /// @return true if the ray (t >= 0) hits the box.
    bool intersectsRay(const Vec3<T> &origin, const Vec3<T> &dir,
                       T *tMin = nullptr, T *tMax = nullptr) const noexcept {
        bool    inside = true;
        Vec3<T> coord;
        Vec3<T> maxT(-1, -1, -1);

        // Per axis: if the origin is outside the slab, the ray can only enter through the near plane.
        for (int i = 0; i < 3; i++) {
            T orig = i == 0 ? origin.x : (i == 1 ? origin.y : origin.z);
            T direction = i == 0 ? dir.x : (i == 1 ? dir.y : dir.z);
            T minVal = i == 0 ? min.x : (i == 1 ? min.y : min.z);
            T maxVal = i == 0 ? max.x : (i == 1 ? max.y : max.z);

            if (orig < minVal) {
                inside = false;
                if (direction != T(0))
                    maxT[i] = (minVal - orig) / direction;
            } else if (orig > maxVal) {
                inside = false;
                if (direction != T(0))
                    maxT[i] = (maxVal - orig) / direction;
            }
        }

        // Origin inside the box: hit at t = 0 by definition.
        if (inside) {
            if (tMin)
                *tMin = T(0);
            if (tMax)
                *tMax = T(0);
            return true;
        }

        // The entry point is on the candidate plane with the largest t (the last slab entered).
        int whichPlane = 0;
        if (maxT.y > maxT.x)
            whichPlane = 1;
        if (maxT.z > maxT[whichPlane])
            whichPlane = 2;

        // All candidates negative (or no candidate): the box is behind the ray.
        if (maxT[whichPlane] < T(0))
            return false;

        // The hit point on the chosen plane must lie within the other two slabs.
        for (int i = 0; i < 3; i++) {
            if (i != whichPlane) {
                T orig = i == 0 ? origin.x : (i == 1 ? origin.y : origin.z);
                T direction = i == 0 ? dir.x : (i == 1 ? dir.y : dir.z);
                T intersect = orig + maxT[whichPlane] * direction;
                T minVal = i == 0 ? min.x : (i == 1 ? min.y : min.z);
                T maxVal = i == 0 ? max.x : (i == 1 ? max.y : max.z);

                if (intersect < minVal || intersect > maxVal)
                    return false;
            }
        }

        if (tMin)
            *tMin = maxT[whichPlane];
        if (tMax) {
            // Exit point is not computed; report the entry t.
            *tMax = maxT[whichPlane];
        }

        return true;
    }

    /// Grows the box to include @p point (never shrinks).
    void expand(const Vec3<T> &point) noexcept {
        min.x = std::min(min.x, point.x);
        min.y = std::min(min.y, point.y);
        min.z = std::min(min.z, point.z);

        max.x = std::max(max.x, point.x);
        max.y = std::max(max.y, point.y);
        max.z = std::max(max.z, point.z);

        initialize(min, max); // keeps mid and corners in sync
    }

    /// Grows the box to include @p other.
    void expand(const AABB &other) noexcept {
        expand(other.min);
        expand(other.max);
    }

    /// Box enclosing this box after transformation by @p m; conservative (not tight) for rotations.
    AABB transformed(const Mat4<T> &m) const noexcept {
        // A rotated box is no longer axis-aligned, so re-fit around all 8 transformed corners.
        std::vector<Vec3<T>> transformedCorners;
        transformedCorners.reserve(8);

        for (const auto &corner : corners) {
            transformedCorners.push_back(m.transformPoint(corner));
        }

        return AABB(transformedCorners);
    }

    /// Smallest box containing both @p a and @p b.
    static AABB merge(const AABB &a, const AABB &b) noexcept {
        Vec3<T> newMin(
            std::min(a.min.x, b.min.x),
            std::min(a.min.y, b.min.y),
            std::min(a.min.z, b.min.z));
        Vec3<T> newMax(
            std::max(a.max.x, b.max.x),
            std::max(a.max.y, b.max.y),
            std::max(a.max.z, b.max.z));
        return AABB(newMin, newMax);
    }
};


using AABBf = AABB<float>;
using AABBd = AABB<double>;


// Vec3 is uploaded to the GPU as-is, so it must stay standard layout.
static_assert(std::is_standard_layout<Vec3f>::value, "Vec3f must be standard layout");

} // namespace decimate::math
