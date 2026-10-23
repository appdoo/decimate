// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file predicates.h
/// @brief Robust geometric predicates: exact sign, adaptive cost.
///
/// After J. R. Shewchuk, "Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric
/// Predicates", Discrete & Computational Geometry 18 (1997); reference code `predicates.c`
/// (public domain).  The fast floating-point filter decides nearly every call; the exact
/// expansion arithmetic only runs near degeneracy.  Float inputs are widened to double
/// (exactly) first.  Every predicate returns a double whose SIGN is exact; its magnitude is
/// only an estimate, so use it for decisions, never as an area.  A non-finite input gives a
/// non-finite result (comparisons against it are false).

#include <decimate/math/vec2.h>
#include <decimate/math/vec3.h>

#include <cstddef>
#include <span>

namespace decimate::math {

/// Orientation of @p c relative to the directed line a->b: > 0 left (counter-clockwise a, b, c),
/// < 0 right, 0 exactly collinear.  Equals the sign of twice the signed area of triangle abc.
[[nodiscard]] double orient2d(Vec2d a, Vec2d b, Vec2d c) noexcept;
[[nodiscard]] double orient2d(Vec2f a, Vec2f b, Vec2f c) noexcept;


/// Side of the plane through a, b, c on which @p d lies: > 0 when d is below the plane, i.e. when a, b, c
/// appear counter-clockwise seen from above (Shewchuk's convention); < 0 above; 0 coplanar.
[[nodiscard]] double orient3d(Vec3d a, Vec3d b, Vec3d c, Vec3d d) noexcept;
[[nodiscard]] double orient3d(Vec3f a, Vec3f b, Vec3f c, Vec3f d) noexcept;


/// > 0 when @p d lies inside the circle through a, b, c, < 0 outside, 0 on it.  Requires a, b, c
/// counter-clockwise (orient2d(a, b, c) > 0); for clockwise input the sign flips.
[[nodiscard]] double incircle(Vec2d a, Vec2d b, Vec2d c, Vec2d d) noexcept;
[[nodiscard]] double incircle(Vec2f a, Vec2f b, Vec2f c, Vec2f d) noexcept;


/// Do segments a-b and c-d cross properly: the interiors meet in exactly one point.  Touching
/// (an end point on the other segment or line), shared end points and collinear overlap are NOT
/// crossings.  Decided by four exact orient2d signs.
[[nodiscard]] bool segmentsCross(Vec2d a, Vec2d b, Vec2d c, Vec2d d) noexcept;
[[nodiscard]] bool segmentsCross(Vec2f a, Vec2f b, Vec2f c, Vec2f d) noexcept;


/// Twice the signed area of the closed polygon @p p (shoelace, double accumulation): > 0 counter-clockwise.
/// A VALUE, not a predicate: rounding error is that of the naive sum.  Use `orient2d` for decisions.
[[nodiscard]] double signedArea2(std::span<const Vec2d> p) noexcept;
[[nodiscard]] double signedArea2(std::span<const Vec2f> p) noexcept;


/// Twice the signed area of triangle abc as a plain double evaluation (the filter's estimate, no
/// exactness): for energies, scales and ratios that need the magnitude.  Its sign can be wrong for
/// slivers; decide with `orient2d`.
[[nodiscard]] constexpr double triangleArea2(Vec2d a, Vec2d b, Vec2d c) noexcept {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

} // namespace decimate::math
