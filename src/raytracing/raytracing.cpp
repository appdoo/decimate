// This file is part of the decimate library; see decimate.h for version/license details.
// mesh-toolkit — core raytracing implementation
// Ported from legacy/raytracing/raytracing.cpp onto decimate::math::Vec3<float>.
//
// Operator mapping (legacy Vec3 vs decimate::math::Vec3):
//   legacy  a * b  (Vec3 * Vec3) = dot product  → a.dot(b)
//   legacy  a ^ b               = cross product → a.cross(b)
//   legacy  a * s  (scalar)                     → a * s  (same)
//   legacy  a + b, a - b, +=, -=                → same
//   legacy  norm()                              → normalize() (in-place) or normalized()
//
// Algorithm: zero behavior change — operand order and branch structure
// are preserved exactly from the 2003 source.

#include <decimate/raytracing/raytracing.h>

#include <cmath>

namespace rtrace {

// ---------------------------------------------------------------------------
// evaluate_t
// ---------------------------------------------------------------------------
bool evaluate_t(float& t, const decimate::math::Vec3<float>& nrm, float d,
                const decimate::math::Vec3<float>& origin, const decimate::math::Vec3<float>& dir) {
    // N·D: if near zero, ray is parallel to plane
    float div = nrm.dot(dir);

    if (std::fabs(div) < kEpsilonF) {
        return false;
    }

    // t = -(d + N·O) / (N·D)
    t = -(d + nrm.dot(origin)) / div;
    return true;
}

// ---------------------------------------------------------------------------
// barycentricBadouel
// ---------------------------------------------------------------------------
bool barycentricBadouel(float& alpha, float& beta, const float* p, const float* nrm,
                        const float* p0, const float* p1, const float* p2) {
    int i1 = 0;
    int i2 = 0;

    // Select the two dominant axes of the face normal (pre-computable).
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

    float u0 = p[i1] - p0[i1];
    float v0 = p[i2] - p0[i2];
    float u1 = p1[i1] - p0[i1];
    float v1 = p1[i2] - p0[i2];
    float u2 = p2[i1] - p0[i1];
    float v2 = p2[i2] - p0[i2];

    if (u1 == 0.0F) {
        beta = u0 / u2;
        if ((beta >= 0.0F) && (beta <= 1.0F)) {
            alpha = (v0 - beta * v2) / v1;
            return (alpha >= 0.0F) && ((alpha + beta) <= 1.0F);
        }
    } else {
        beta = (v0 * u1 - u0 * v1) / (v2 * u1 - u2 * v1);
        if ((beta >= 0.0F) && (beta <= 1.0F)) {
            alpha = (u0 - beta * u2) / u1;
            return (alpha >= 0.0F) && ((alpha + beta) <= 1.0F);
        }
    }

    return false;
}

// ---------------------------------------------------------------------------
// intersect_triangle_bfc  (backface culling ON)
// ---------------------------------------------------------------------------
bool intersect_triangle_bfc(decimate::math::Vec3<float>& orig, decimate::math::Vec3<float>& dir,
                            decimate::math::Vec3<float>& v0, decimate::math::Vec3<float>& v1, decimate::math::Vec3<float>& v2,
                            float* t, float* u, float* v) {
    // Edge vectors from v0
    decimate::math::Vec3<float> e1 = v1 - v0;
    decimate::math::Vec3<float> e2 = v2 - v0;

    // pv = dir × e2 (used for determinant and u parameter)
    decimate::math::Vec3<float> pv = dir.cross(e2);

    // det = e1 · pv
    float det = e1.dot(pv);

    // Backface culling: reject if determinant is not clearly positive
    if (det < 1e-12F) {
        return false;
    }

    // tv = orig - v0
    decimate::math::Vec3<float> tv = orig - v0;

    // u parameter
    *u = tv.dot(pv);
    if (*u < 0.0F || *u > det) {
        return false;
    }

    // qv = tv × e1
    decimate::math::Vec3<float> qv = tv.cross(e1);

    // v parameter
    *v = dir.dot(qv);
    if (*v < 0.0F || *u + *v > det) {
        return false;
    }

    // Scale to barycentric coordinates
    *t = e2.dot(qv);
    float invDet = 1.0F / det;
    *t *= invDet;
    *u *= invDet;
    *v *= invDet;

    return true;
}

// ---------------------------------------------------------------------------
// intersect_triangle  (backface culling OFF)
// ---------------------------------------------------------------------------
bool intersect_triangle(decimate::math::Vec3<float>& orig, decimate::math::Vec3<float>& dir, decimate::math::Vec3<float>& v0,
                        decimate::math::Vec3<float>& v1, decimate::math::Vec3<float>& v2, float* t, float* u,
                        float* v) {
    // Edge vectors from v0
    decimate::math::Vec3<float> e1 = v1 - v0;
    decimate::math::Vec3<float> e2 = v2 - v0;

    // pv = dir × e2
    decimate::math::Vec3<float> pv = dir.cross(e2);

    // det = e1 · pv
    float det = e1.dot(pv);

    // tv = orig - v0
    decimate::math::Vec3<float> tv = orig - v0;

    // qv = tv × e1
    decimate::math::Vec3<float> qv = tv.cross(e1);

    if (det > kEpsilonF) {
        // u parameter
        *u = tv.dot(pv);
        if (*u < 0.0F || *u > det) {
            return false;
        }

        // v parameter
        *v = dir.dot(qv);
        if (*v < 0.0F || *u + *v > det) {
            return false;
        }
    } else if (det < -kEpsilonF) {
        // u parameter (negative determinant branch)
        *u = tv.dot(pv);
        if (*u > 0.0F || *u < det) {
            return false;
        }

        // v parameter
        *v = dir.dot(qv);
        if (*v > 0.0F || *u + *v < det) {
            return false;
        }
    } else {
        // Ray is parallel to the plane of the triangle
        return false;
    }

    // Scale to barycentric coordinates
    *t = e2.dot(qv) / det;
    *u /= det;
    *v /= det;

    return true;
}

// ---------------------------------------------------------------------------
// IntersectTriangle  (Drew's adaptation from Real Time Rendering)
// ---------------------------------------------------------------------------
bool IntersectTriangle(decimate::math::Vec3<float>& orig, decimate::math::Vec3<float>& dir, decimate::math::Vec3<float>& v1,
                       decimate::math::Vec3<float>& v2, decimate::math::Vec3<float>& v3, float* t, float* u,
                       float* v) {
    // Edge vectors from v1
    decimate::math::Vec3<float> edge1 = v2 - v1;
    decimate::math::Vec3<float> edge2 = v3 - v1;

    // pvec = dir × edge2
    decimate::math::Vec3<float> pvec = dir.cross(edge2);

    // det = edge1 · pvec
    float det = edge1.dot(pvec);

    // tvec = orig - v1
    decimate::math::Vec3<float> tvec = orig - v1;

    // qvec = tvec × edge1
    decimate::math::Vec3<float> qvec = tvec.cross(edge1);

    if (det > kEpsilonF) {
        *u = tvec.dot(pvec);
        if (*u < 0.0F || *u > det) {
            return false;
        }

        *v = dir.dot(qvec);
        if (*v < 0.0F || *u + *v > det) {
            return false;
        }
    } else if (det < -kEpsilonF) {
        *u = tvec.dot(pvec);
        if (*u > 0.0F || *u < det) {
            return false;
        }

        *v = dir.dot(qvec);
        if (*v > 0.0F || *u + *v < det) {
            return false;
        }
    } else {
        // Ray is parallel to the plane of the triangle
        return false;
    }

    float invDet = 1.0F / det;

    *t = edge2.dot(qvec) * invDet;
    (*u) *= invDet;
    (*v) *= invDet;

    return true;
}

}  // namespace rtrace
