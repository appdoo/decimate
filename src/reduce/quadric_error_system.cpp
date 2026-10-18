// This file is part of the decimate library; see decimate.h for version/license details.
/// @file quadric_error_system.cpp
/// @brief `QuadricErrorSystem` implementation — port of `CQuadricErrorSystem`.

#include <decimate/reduce/quadric_error_system.h>

#include <decimate/math/vec3.h>

namespace decimate {

void QuadricErrorSystem::initialize(int numPoints, const decimate::math::Vec3f* points,
                                    int numFaces,  const Face* faces) {
    pointRef_ = points;
    errorMetrics_.assign(static_cast<std::size_t>(numPoints), QError{});

    // For each face: compute face-plane quadric and add to each corner vertex.
    // O(F) — matches legacy cquadricerrorsystem.cpp:44
    for (int i = 0; i < numFaces; ++i) {
        const Face& f = faces[i];
        QError faceQuadric;
        faceQuadric.setPlane(points[f.a], points[f.b], points[f.c]);
        errorMetrics_[static_cast<std::size_t>(f.a)].add(faceQuadric);
        errorMetrics_[static_cast<std::size_t>(f.b)].add(faceQuadric);
        errorMetrics_[static_cast<std::size_t>(f.c)].add(faceQuadric);
    }
}

void QuadricErrorSystem::update(int p1, int p2) {
    // Merge p2's quadric into p1 — O(1) scalar additions.
    // Matches legacy cquadricerrorsystem.cpp:81
    errorMetrics_[static_cast<std::size_t>(p1)].add(
        errorMetrics_[static_cast<std::size_t>(p2)]);
}

void QuadricErrorSystem::optimize(int p1, int p2, bool boundary,
                                  double& error, decimate::math::Vec3f& p) const {
    // Combined edge quadric on the stack (no heap allocation).
    // Matches legacy cquadricerrorsystem.cpp:90
    QError combined(errorMetrics_[static_cast<std::size_t>(p1)],
                    errorMetrics_[static_cast<std::size_t>(p2)]);

    // If boundary or inversion fails: midpoint fallback (Beat-2a — F5 deferred).
    if (boundary || !combined.optimize(p)) {
        p = (pointRef_[p1] + pointRef_[p2]) * 0.5F;
    }
    error = combined.getError(p);
}

}  // namespace decimate
