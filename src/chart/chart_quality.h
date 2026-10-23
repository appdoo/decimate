// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file chart_quality.h
/// @brief Private: the O(n + B^2) flattening-quality check (B = boundary edges) shared by the LSCM and the bijective mapper.

#include <decimate/chart/parametrize.h>
#include <decimate/math/predicates.h>

namespace decimate::chart {

/// Twice the signed area of (a, b, c) in double, from float UVs: a VALUE (the sign of a sliver can be wrong; decide
/// with `math::orient2d`).  Bit-identical to the former local determinant.
[[nodiscard]] inline double orientValue(float ax, float ay, float bx, float by, float cx, float cy) noexcept {
    return math::triangleArea2(math::Vec2d { ax, ay }, math::Vec2d { bx, by }, math::Vec2d { cx, cy });
}


/// Fill `ac.flippedTriangles` / `ac.boundarySimple` from `ac.uv` (resets both first): triangles
/// whose UV orientation opposes the chart's (sign of the summed signed area) or is degenerate, and
/// whether two boundary edges cross.  Disk + no flip + simple boundary => injective.
void measureQuality(HalfEdgeChart &ac);

} // namespace decimate::chart
