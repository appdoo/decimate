// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file sat_overlap.h
/// @brief Exact 2D triangle-triangle SAT overlap, shared by the atlas tests (antipattern 12).

#include <algorithm>
#include <limits>

#include <decimate/math/vec2.h>

/// 2D SAT overlap test for two triangles p0,p1,p2 and q0,q1,q2.
/// Returns true if the triangles overlap with INTERIOR area (share more than a
/// boundary edge or vertex).  Epsilon ensures floating-point edge-touching is
/// not reported as overlap: the atlas is valid when triangles touch at seams.
/// kEps is in UV units; charts are in [0,1]^2, so 1e-5 is 0.001% of atlas width.
inline bool trianglesOverlapSAT(
    const decimate::math::Vec2f &p0, const decimate::math::Vec2f &p1, const decimate::math::Vec2f &p2,
    const decimate::math::Vec2f &q0, const decimate::math::Vec2f &q1, const decimate::math::Vec2f &q2) {
    constexpr float kEps = 1e-5f; // touching tolerance in UV space

    const decimate::math::Vec2f pVerts[3] = { p0, p1, p2 };
    const decimate::math::Vec2f qVerts[3] = { q0, q1, q2 };

    // Separation test on one axis: true when the two projected intervals are apart.
    auto testAxis = [&](float ax, float ay) -> bool {
        float pMin = std::numeric_limits<float>::max();
        float pMax = std::numeric_limits<float>::lowest();
        float qMin = std::numeric_limits<float>::max();
        float qMax = std::numeric_limits<float>::lowest();
        for (const auto &v : pVerts) {
            float d = v.x * ax + v.y * ay;
            pMin = std::min(pMin, d);
            pMax = std::max(pMax, d);
        }
        for (const auto &v : qVerts) {
            float d = v.x * ax + v.y * ay;
            qMin = std::min(qMin, d);
            qMax = std::max(qMax, d);
        }

        // Separated if intervals are disjoint by more than kEps.
        // Touching within kEps is not counted as interior overlap.
        return (pMax <= qMin + kEps || qMax <= pMin + kEps);
    };

    // The edge normals of both triangles are the only candidate separating axes.
    for (int i = 0; i < 3; ++i) {
        const decimate::math::Vec2f &a = pVerts[i];
        const decimate::math::Vec2f &b = pVerts[(i + 1) % 3];
        const float                  ex = b.x - a.x;
        const float                  ey = b.y - a.y;
        if (testAxis(-ey, ex)) {
            return false;
        } // separated on this axis
    }

    for (int i = 0; i < 3; ++i) {
        const decimate::math::Vec2f &a = qVerts[i];
        const decimate::math::Vec2f &b = qVerts[(i + 1) % 3];
        const float                  ex = b.x - a.x;
        const float                  ey = b.y - a.y;
        if (testAxis(-ey, ex)) {
            return false;
        } // separated on this axis
    }

    // No separating axis found: the interiors overlap.
    return true;
}
