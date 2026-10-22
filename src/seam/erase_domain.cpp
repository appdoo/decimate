// This file is part of the decimate library; see decimate.h for version/license details.
/// @file erase_domain.cpp
/// @brief The unknown / free texel sets of the erase solve.

#include "erase_domain.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

namespace decimate::seam::detail {

namespace {

constexpr double kInsideTolerance = 1e-8; // barycentric slack for "texel centre in triangle"


/// A UV triangle in texel coordinates.
struct TexelTriangle {
    std::array<double, 3> x {};
    std::array<double, 3> y {};
};

} // namespace


/// Mark the four texels of every bilinear cell the segment visits.
static void markSegmentCells(const Segment &seg, TexelGrid grid, std::vector<char> &mask) {
    for (const Cell &c : cellsAlong(seg, grid)) {
        for (const int texel :
             { c.py * grid.width + c.px, c.py * grid.width + c.px + 1,
               (c.py + 1) * grid.width + c.px, (c.py + 1) * grid.width + c.px + 1 }) {
            mask[static_cast<std::size_t>(texel)] = 1;
        }
    }
}


/// Mark all cells of the UV silhouette: both sides of every seam, boundary, foldovers.
static std::vector<char> silhouetteCells(const TexturedMesh &tm, const SeamTopology &topo,
                                         TexelGrid grid) {
    std::vector<char> mask(static_cast<std::size_t>(grid.width) * grid.height, 0);

    // A seam is sampled from both sides, so both half-edges contribute cells.
    for (const SeamPair &pair : topo.seams) {
        markSegmentCells(halfEdgeSegment(tm, pair.a, grid), grid, mask);
        markSegmentCells(halfEdgeSegment(tm, pair.b, grid), grid, mask);
    }

    // Boundary and foldover edges leave the parametrization too: texels beside them are unreliable.
    for (const HalfEdgeRef &he : topo.boundary) {
        markSegmentCells(halfEdgeSegment(tm, he, grid), grid, mask);
    }
    for (const HalfEdgeRef &he : topo.foldovers) {
        markSegmentCells(halfEdgeSegment(tm, he, grid), grid, mask);
    }
    return mask;
}


/// Corner UVs of @p face in texel coordinates (`uv * N - 0.5`).
static TexelTriangle texelTriangle(const TexturedMesh &tm, std::size_t face, TexelGrid grid) {
    TexelTriangle t;
    for (std::size_t k = 0; k < 3; ++k) {
        const math::Vec2f uv = tm.uvs[static_cast<std::size_t>(tm.faceU[3 * face + k])];
        t.x[k] = static_cast<double>(uv.x) * grid.width - 0.5;
        t.y[k] = static_cast<double>(uv.y) * grid.height - 0.5;
    }
    return t;
}


/// True when texel centre (px,py) lies in the closed triangle (barycentrics >= -tol).
static bool containsTexel(const TexelTriangle &t, double px, double py) {
    const double det =
        (t.x[1] - t.x[0]) * (t.y[2] - t.y[0]) - (t.x[2] - t.x[0]) * (t.y[1] - t.y[0]);

    // A zero-area triangle covers no texel centre.
    if (det == 0.0) {
        return false;
    }

    const double l1 = ((px - t.x[0]) * (t.y[2] - t.y[0]) - (t.x[2] - t.x[0]) * (py - t.y[0])) / det;
    const double l2 = ((t.x[1] - t.x[0]) * (py - t.y[0]) - (px - t.x[0]) * (t.y[1] - t.y[0])) / det;
    const double l0 = 1.0 - l1 - l2;
    return l0 >= -kInsideTolerance && l1 >= -kInsideTolerance && l2 >= -kInsideTolerance;
}


/// Texel index range [lo,hi] a triangle extent [minC,maxC] can touch along one axis.
static std::pair<int, int> texelRange(double minC, double maxC, int size) {
    const double limit = static_cast<double>(size);
    const int    lo = std::max(0, static_cast<int>(std::floor(std::max(minC, -1.0))));
    const int    hi = std::min(size - 1, static_cast<int>(std::ceil(std::min(maxC, limit))));
    return { lo, hi };
}


/// Mark the texels of row @p y, columns [x0,x1], whose centre lies in @p t.
static void markRowTexels(const TexelTriangle &t, int y, int x0, int x1, TexelGrid grid,
                          std::vector<char> &mask) {
    for (int x = x0; x <= x1; ++x) {
        if (containsTexel(t, x, y)) {
            mask[static_cast<std::size_t>(y) * grid.width + x] = 1;
        }
    }
}


/// Mark every texel whose centre lies in @p t, scanning only the triangle's bounding box.
static void markTriangleTexels(const TexelTriangle &t, TexelGrid grid, std::vector<char> &mask) {
    const auto [x0, x1] = texelRange(std::min({ t.x[0], t.x[1], t.x[2] }),
                                     std::max({ t.x[0], t.x[1], t.x[2] }), grid.width);
    const auto [y0, y1] = texelRange(std::min({ t.y[0], t.y[1], t.y[2] }),
                                     std::max({ t.y[0], t.y[1], t.y[2] }), grid.height);
    for (int y = y0; y <= y1; ++y) {
        markRowTexels(t, y, x0, x1, grid, mask);
    }
}


/// The texels whose centre lies inside some UV triangle.
static std::vector<char> faceTexels(const TexturedMesh &tm, TexelGrid grid) {
    std::vector<char> mask(static_cast<std::size_t>(grid.width) * grid.height, 0);
    for (std::size_t f = 0; f < tm.faceCount(); ++f) {
        markTriangleTexels(texelTriangle(tm, f, grid), grid, mask);
    }
    return mask;
}


EraseDomain buildEraseDomain(const TexturedMesh &tm, const SeamTopology &topo, TexelGrid grid) {
    const std::vector<char> silhouette = silhouetteCells(tm, topo, grid);
    const std::vector<char> inside = faceTexels(tm, grid);

    // Unknown = anything inside the UV mesh or in a silhouette cell; free = silhouette texels
    // the bake never wrote, which must not be anchored to the input.
    EraseDomain domain;
    domain.unknown.resize(silhouette.size());
    domain.free.resize(silhouette.size());
    for (std::size_t i = 0; i < silhouette.size(); ++i) {
        domain.unknown[i] = static_cast<char>(silhouette[i] != 0 || inside[i] != 0);
        domain.free[i] = static_cast<char>(silhouette[i] != 0 && inside[i] == 0);
        domain.freeCount += domain.free[i];
    }
    return domain;
}

} // namespace decimate::seam::detail
