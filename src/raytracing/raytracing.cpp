// This file is part of the decimate library; see decimate.h for version/license details.
/// @file raytracing.cpp
/// @brief Ray/plane and ray/triangle kernels, ported from legacy2003/raytracing/raytracing.cpp onto
/// decimate::math::Vec3f. Contracts are in raytracing.h.
///
/// Operator mapping (legacy Vec3 vs decimate::math::Vec3):
///   legacy  a * b  (Vec3 * Vec3) = dot product  -> a.dot(b)
///   legacy  a ^ b               = cross product -> a.cross(b)
///   legacy  a * s  (scalar)                     -> a * s  (same)
///   legacy  a + b, a - b, +=, -=                -> same
///   legacy  norm()                              -> normalize() (in-place) or normalized()
///
/// Operand order and branch structure are kept exactly as in the 2003 source: the bake's bit-exact
/// equality with the brute-force loop (see bvh.h) depends on the same floats going through the same
/// operations.

#include <decimate/raytracing/raytracing.h>

#include <cmath>

/// File-local shorthand (an alias, not a using-directive): this file is outside namespace decimate.
namespace math = decimate::math;

namespace rtrace {

bool evaluate_t(float &t, const math::Vec3f &nrm, float d,
                const math::Vec3f &origin, const math::Vec3f &dir) {
    float div = nrm.dot(dir);

    // N.D near zero: the ray runs parallel to the plane, there is no intersection to solve for.
    if (std::fabs(div) < kEpsilonF) {
        return false;
    }

    t = -(d + nrm.dot(origin)) / div;
    return true;
}


bool barycentricBadouel(float &alpha, float &beta, const float *p, const float *nrm,
                        const float *p0, const float *p1, const float *p2) {
    int i1 = 0;
    int i2 = 0;

    // Project onto the plane that drops the normal's dominant axis, keeping the projected area largest
    // (the choice depends on the face only, so it could be precomputed).
    float nx = std::fabs(nrm[0]);
    float ny = std::fabs(nrm[1]);
    float nz = std::fabs(nrm[2]);

    if (nx > ny && nx > nz) {
        i1 = 1;
        i2 = 2;
    } else if (ny > nx && ny > nz) {
        i1 = 0;
        i2 = 2;
    } else {
        i1 = 0;
        i2 = 1;
    }

    // Triangle edges and the point, relative to p0, in the projected plane.
    float u0 = p[i1] - p0[i1];
    float v0 = p[i2] - p0[i2];
    float u1 = p1[i1] - p0[i1];
    float v1 = p1[i2] - p0[i2];
    float u2 = p2[i1] - p0[i1];
    float v2 = p2[i2] - p0[i2];

    // Solve p = p0 + alpha*(p1-p0) + beta*(p2-p0); u1 == 0 selects the formulation that avoids dividing by it.
    if (u1 == 0.0f) {
        beta = u0 / u2;
        if ((beta >= 0.0f) && (beta <= 1.0f)) {
            alpha = (v0 - beta * v2) / v1;
            return (alpha >= 0.0f) && ((alpha + beta) <= 1.0f);
        }
    } else {
        beta = (v0 * u1 - u0 * v1) / (v2 * u1 - u2 * v1);
        if ((beta >= 0.0f) && (beta <= 1.0f)) {
            alpha = (u0 - beta * u2) / u1;
            return (alpha >= 0.0f) && ((alpha + beta) <= 1.0f);
        }
    }

    return false;
}


bool intersect_triangle_bfc(math::Vec3f &orig, math::Vec3f &dir,
                            math::Vec3f &v0, math::Vec3f &v1, math::Vec3f &v2,
                            float *t, float *u, float *v) {
    math::Vec3f e1 = v1 - v0;
    math::Vec3f e2 = v2 - v0;
    math::Vec3f pv = dir.cross(e2);
    float       det = e1.dot(pv);

    // Backface culling: the determinant's sign is the facing, so a non-positive one is a backface
    // (or a ray in the triangle's plane).
    if (det < 1e-12f) {
        return false;
    }

    math::Vec3f tv = orig - v0;

    // Reject on u and v before the remaining cross product and the division are paid for.
    *u = tv.dot(pv);
    if (*u < 0.0f || *u > det) {
        return false;
    }

    math::Vec3f qv = tv.cross(e1);
    *v = dir.dot(qv);
    if (*v < 0.0f || *u + *v > det) {
        return false;
    }

    // Hit: scale t, u, v by 1/det (the tests above ran on unscaled values).
    *t = e2.dot(qv);
    float invDet = 1.0f / det;
    *t *= invDet;
    *u *= invDet;
    *v *= invDet;

    return true;
}


bool intersect_triangle(math::Vec3f &orig, math::Vec3f &dir, math::Vec3f &v0,
                        math::Vec3f &v1, math::Vec3f &v2, float *t, float *u,
                        float *v) {
    math::Vec3f e1 = v1 - v0;
    math::Vec3f e2 = v2 - v0;
    math::Vec3f pv = dir.cross(e2);
    float       det = e1.dot(pv);
    math::Vec3f tv = orig - v0;
    math::Vec3f qv = tv.cross(e1);

    // Both signs of det are hits (no culling); a negative det mirrors every inequality.
    if (det > kEpsilonF) {
        *u = tv.dot(pv);
        if (*u < 0.0f || *u > det) {
            return false;
        }

        *v = dir.dot(qv);
        if (*v < 0.0f || *u + *v > det) {
            return false;
        }
    } else if (det < -kEpsilonF) {
        *u = tv.dot(pv);
        if (*u > 0.0f || *u < det) {
            return false;
        }

        *v = dir.dot(qv);
        if (*v > 0.0f || *u + *v < det) {
            return false;
        }
    } else {
        // Ray is parallel to the plane of the triangle
        return false;
    }

    // Hit: scale t, u, v by 1/det. This variant divides; the other two multiply by the reciprocal. The
    // results can differ in the last bit, so the three are not interchangeable.
    *t = e2.dot(qv) / det;
    *u /= det;
    *v /= det;

    return true;
}


bool IntersectTriangle(math::Vec3f &orig, math::Vec3f &dir, math::Vec3f &v1,
                       math::Vec3f &v2, math::Vec3f &v3, float *t, float *u,
                       float *v) {
    math::Vec3f edge1 = v2 - v1;
    math::Vec3f edge2 = v3 - v1;
    math::Vec3f pvec = dir.cross(edge2);
    float       det = edge1.dot(pvec);
    math::Vec3f tvec = orig - v1;
    math::Vec3f qvec = tvec.cross(edge1);

    // Same two-sided test as intersect_triangle.
    if (det > kEpsilonF) {
        *u = tvec.dot(pvec);
        if (*u < 0.0f || *u > det) {
            return false;
        }

        *v = dir.dot(qvec);
        if (*v < 0.0f || *u + *v > det) {
            return false;
        }
    } else if (det < -kEpsilonF) {
        *u = tvec.dot(pvec);
        if (*u > 0.0f || *u < det) {
            return false;
        }

        *v = dir.dot(qvec);
        if (*v > 0.0f || *u + *v < det) {
            return false;
        }
    } else {
        // Ray is parallel to the plane of the triangle
        return false;
    }

    // Hit: scale t, u, v by 1/det.
    float invDet = 1.0f / det;

    *t = edge2.dot(qvec) * invDet;
    (*u) *= invDet;
    (*v) *= invDet;

    return true;
}

} // namespace rtrace
