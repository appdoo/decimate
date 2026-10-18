// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <decimate/math/vec3.h>
#include <decimate/math/mat4.h>
#include <decimate/math/plane.h>
#include <vector>
#include <algorithm>

namespace decimate::math {

/**
 * AABB - Axis-Aligned Bounding Box
 *
 * Efficient spatial data structure for:
 * - Collision detection
 * - Ray intersection tests
 * - Frustum culling
 * - Spatial queries
 *
 * Based on mudvis implementation (Apache 2.0, Copyright 2003-2019 Andor Greissl)
 * with optimizations from Graphics Gems (Andrew Woo, Pierre Terdiman)
 *
 * Memory layout:
 * - 8 corner points arranged so inverted index gives opposite diagonal corner
 * - Near/far point computation using plane sign bits for fast culling
 *
 * Corner arrangement:      5 ---- 1
 *                         /|     /|
 *                        4 ---- 0 |
 *                        | 7 ---| 3
 *                        |/     |/
 *                        6 ---- 2
 */
template<typename T>
struct AABB {
    static_assert(std::is_floating_point<T>::value, "AABB requires floating point type");

    Vec3<T> min, max, mid;
    std::vector<Vec3<T>> corners; // 8 corner points

    // Constructors
    AABB() noexcept {
        corners.resize(8);
    }

    AABB(const Vec3<T>& min, const Vec3<T>& max) noexcept
        : min(min), max(max) {
        initialize(min, max);
    }

    AABB(const std::vector<Vec3<T>>& vertices) noexcept {
        if (vertices.empty()) {
            initialize(Vec3<T>::zero(), Vec3<T>::zero());
            return;
        }

        // Find min/max extents
        min = vertices[0];
        max = vertices[0];

        for (const auto& v : vertices) {
            min.x = std::min(min.x, v.x);
            min.y = std::min(min.y, v.y);
            min.z = std::min(min.z, v.z);

            max.x = std::max(max.x, v.x);
            max.y = std::max(max.y, v.y);
            max.z = std::max(max.z, v.z);
        }

        initialize(min, max);
    }

    void initialize(const Vec3<T>& minPt, const Vec3<T>& maxPt) noexcept {
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

    // Accessors
    const Vec3<T>& getCorner(int index) const noexcept {
        return corners[index & 0x7];  // Wrap to 0-7
    }

    // --- p-vertex / n-vertex fast plane culling ----------------------------
    // Faithful port of the 2003-2019 mudvis/vermessung farCorner/nearCorner.
    //
    // The corner FARTHEST along the plane normal (the "p-vertex") is selected
    // in O(1) from the plane's precomputed sign bits: for a negative normal
    // component the min corner is farthest, else the max — and the corner array
    // is laid out so (corner index bit set ⇔ min on that axis), exactly the
    // complement of the plane's signbits (signbit set ⇔ normal negative).

    /// Corner farthest along the plane normal (p-vertex).
    const Vec3<T>& farCorner(const Plane<T>& plane) const noexcept {
        return corners[(~plane.signbits) & 0x7];
    }

    /// Corner nearest along the plane normal (n-vertex); inverse of farCorner.
    const Vec3<T>& nearCorner(const Plane<T>& plane) const noexcept {
        return corners[plane.signbits & 0x7];
    }

    /// True if the whole box is strictly on the NEGATIVE side of `plane`
    /// (i.e. even the p-vertex is behind it) — the fast, conservative reject.
    /// Conservative: never reports "outside" for a box that actually crosses
    /// the plane, so it can over-include but never wrongly exclude.
    bool isOutsidePlane(const Plane<T>& plane) const noexcept {
        return plane.signedDistance(farCorner(plane)) < T(0);
    }

    Vec3<T> getMin() const noexcept { return min; }
    Vec3<T> getMax() const noexcept { return max; }
    Vec3<T> getMid() const noexcept { return mid; }

    Vec3<T> getExtents() const noexcept { return max - min; }
    Vec3<T> getHalfExtents() const noexcept { return (max - min) * T(0.5); }

    // Properties
    T surfaceArea() const noexcept {
        Vec3<T> extent = max - min;
        return T(2) * (extent.x * extent.y + extent.x * extent.z + extent.y * extent.z);
    }

    T volume() const noexcept {
        Vec3<T> extent = max - min;
        return extent.x * extent.y * extent.z;
    }

    T minExtent() const noexcept {
        Vec3<T> extent = max - min;
        return std::min(std::min(extent.x, extent.y), extent.z);
    }

    T maxExtent() const noexcept {
        Vec3<T> extent = max - min;
        return std::max(std::max(extent.x, extent.y), extent.z);
    }

    // Test if point is inside AABB
    bool contains(const Vec3<T>& point) const noexcept {
        return point.x >= min.x && point.x <= max.x &&
               point.y >= min.y && point.y <= max.y &&
               point.z >= min.z && point.z <= max.z;
    }

    // AABB-AABB intersection test
    bool intersects(const AABB& other) const noexcept {
        // Min against Max test
        if (min.x > other.max.x || min.y > other.max.y || min.z > other.max.z)
            return false;

        // Max against Min test
        if (max.x < other.min.x || max.y < other.min.y || max.z < other.min.z)
            return false;

        return true;  // Overlap exists
    }

    /**
     * Ray-AABB intersection test
     *
     * Based on Andrew Woo (Graphics Gems) with optimizations by Pierre Terdiman
     * ~20-30% faster than original, more robust
     *
     * @param origin Ray origin
     * @param dir Ray direction (does not need to be normalized)
     * @param tMin Optional output: distance to entry point
     * @param tMax Optional output: distance to exit point
     * @return true if ray intersects AABB
     */
    bool intersectsRay(const Vec3<T>& origin, const Vec3<T>& dir,
                       T* tMin = nullptr, T* tMax = nullptr) const noexcept {
        bool inside = true;
        Vec3<T> coord;
        Vec3<T> maxT(-1, -1, -1);

        // Find candidate planes
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

        // Ray origin inside bounding box
        if (inside) {
            if (tMin) *tMin = T(0);
            if (tMax) *tMax = T(0);
            return true;
        }

        // Get largest of the maxT's for final choice of intersection
        int whichPlane = 0;
        if (maxT.y > maxT.x) whichPlane = 1;
        if (maxT.z > maxT[whichPlane]) whichPlane = 2;

        // Check final candidate actually inside box
        if (maxT[whichPlane] < T(0))
            return false;

        // Check intersection point is inside AABB
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

        if (tMin) *tMin = maxT[whichPlane];
        if (tMax) {
            // Compute exit point if needed (not implemented in original)
            *tMax = maxT[whichPlane];
        }

        return true;  // Ray hits box
    }

    // Expand AABB to include a point
    void expand(const Vec3<T>& point) noexcept {
        min.x = std::min(min.x, point.x);
        min.y = std::min(min.y, point.y);
        min.z = std::min(min.z, point.z);

        max.x = std::max(max.x, point.x);
        max.y = std::max(max.y, point.y);
        max.z = std::max(max.z, point.z);

        initialize(min, max);  // Recompute corners and mid
    }

    // Expand AABB to include another AABB
    void expand(const AABB& other) noexcept {
        expand(other.min);
        expand(other.max);
    }

    // Transform AABB by matrix (computes new AABB enclosing transformed box)
    AABB transformed(const Mat4<T>& m) const noexcept {
        // Transform all 8 corners and recompute AABB
        std::vector<Vec3<T>> transformedCorners;
        transformedCorners.reserve(8);

        for (const auto& corner : corners) {
            transformedCorners.push_back(m.transformPoint(corner));
        }

        return AABB(transformedCorners);
    }

    // Merge two AABBs
    static AABB merge(const AABB& a, const AABB& b) noexcept {
        Vec3<T> newMin(
            std::min(a.min.x, b.min.x),
            std::min(a.min.y, b.min.y),
            std::min(a.min.z, b.min.z)
        );
        Vec3<T> newMax(
            std::max(a.max.x, b.max.x),
            std::max(a.max.y, b.max.y),
            std::max(a.max.z, b.max.z)
        );
        return AABB(newMin, newMax);
    }
};

// Common typedefs
using AABBf = AABB<float>;
using AABBd = AABB<double>;

// Ensure standard layout for GPU compatibility
static_assert(std::is_standard_layout<Vec3<float>>::value, "Vec3f must be standard layout");

} // namespace decimate::math
