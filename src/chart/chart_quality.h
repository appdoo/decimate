// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file chart_quality.h
/// @brief Private: the O(n + B^2) flattening-quality check (B = boundary edges) shared by the LSCM and the bijective mapper.

#include <decimate/chart/parametrize.h>

namespace decimate::chart {

/// Orientation determinant of (a, b, c) in double: > 0 when c is left of a->b.
[[nodiscard]] inline double orient2d(float ax, float ay, float bx, float by, float cx, float cy) noexcept {
    return (static_cast<double>(bx) - ax) * (static_cast<double>(cy) - ay) -
           (static_cast<double>(by) - ay) * (static_cast<double>(cx) - ax);
}


/// Do segments a-b and c-d cross properly (interiors intersect at one point, no shared end point)?
/// The one crossing predicate of the O(B^2) simple-boundary checks (`measureQuality`, `hasSelfCrossing`).
[[nodiscard]] inline bool segmentsCross(float ax, float ay, float bx, float by, float cx, float cy, float dx,
                                        float dy) noexcept {
    const double d1 = orient2d(cx, cy, dx, dy, ax, ay);
    const double d2 = orient2d(cx, cy, dx, dy, bx, by);
    const double d3 = orient2d(ax, ay, bx, by, cx, cy);
    const double d4 = orient2d(ax, ay, bx, by, dx, dy);
    return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
}


/// Fill `ac.flippedTriangles` / `ac.boundarySimple` from `ac.uv` (resets both first): triangles
/// whose UV orientation opposes the chart's (sign of the summed signed area) or is degenerate, and
/// whether two boundary edges cross.  Disk + no flip + simple boundary => injective.
void measureQuality(HalfEdgeChart &ac);

} // namespace decimate::chart
