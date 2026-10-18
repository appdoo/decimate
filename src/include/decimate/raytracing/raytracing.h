// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

// mesh-toolkit — core raytracing functions
// Ported from legacy/raytracing/ onto decimate::math::Vec3<float>.
// Zero OS/graphics dependencies — pure C++20, Wasm-compilable.

#include <decimate/math/vec3.h>

namespace rtrace {

/// Epsilon for parallel-ray / degenerate-triangle tests.
/// Preserved exactly from legacy/math/math3D.h (EPSILONF = 1e-8f).
inline constexpr float kEpsilonF = 1e-8f;

// ---------------------------------------------------------------------------
// evaluate_t
// ---------------------------------------------------------------------------
/// Compute the ray/plane intersection parameter t.
///
/// Solves:  t = -(d + N·O) / (N·D)
///
/// @param[out] t      Intersection parameter along the ray.
/// @param      nrm    Plane normal vector N.
/// @param      d      Plane equation constant: N·P + d = 0.
/// @param      origin Ray origin O.
/// @param      dir    Ray direction D (need not be unit-length).
/// @return true  if the ray intersects the plane (N·D not near zero).
/// @return false if the ray is parallel to the plane.
bool evaluate_t(float& t, const decimate::math::Vec3<float>& nrm, float d,
                const decimate::math::Vec3<float>& origin, const decimate::math::Vec3<float>& dir);

// ---------------------------------------------------------------------------
// barycentricBadouel
// ---------------------------------------------------------------------------
/// Ray/polygon intersection via Didier Badouel's barycentric method.
///
/// Reference: "An Efficient Ray/Polygon Intersection", Graphics Gems (1990).
///
/// @param[out] alpha  Barycentric coordinate alpha.
/// @param[out] beta   Barycentric coordinate beta.
/// @param      p      Intersection point (float[3]).
/// @param      nrm    Face normal (float[3]).
/// @param      p0     Triangle vertex 0 (float[3]).
/// @param      p1     Triangle vertex 1 (float[3]).
/// @param      p2     Triangle vertex 2 (float[3]).
/// @return true if the point lies inside the triangle.
bool barycentricBadouel(float& alpha, float& beta, const float* p, const float* nrm,
                        const float* p0, const float* p1, const float* p2);

// ---------------------------------------------------------------------------
// intersect_triangle_bfc  (backface culling enabled)
// ---------------------------------------------------------------------------
/// Möller–Trumbore ray/triangle intersection — backface culling ON.
///
/// Reference: "Fast, Minimum Storage Ray/Triangle Intersection",
///   Tomas Möller & Ben Trumbore (1997).
///   http://www.graphics.cornell.edu/pubs/1997/MT97.html
///
/// Rejects triangles whose determinant is < 1e-12 (i.e. backfaces and
/// near-degenerate triangles are culled).
///
/// @param orig  Ray origin.
/// @param dir   Ray direction.
/// @param v0    Triangle vertex 0.
/// @param v1    Triangle vertex 1.
/// @param v2    Triangle vertex 2.
/// @param[out] t  Distance parameter along the ray.
/// @param[out] u  Barycentric u coordinate.
/// @param[out] v  Barycentric v coordinate.
/// @return true if the ray hits the front face of the triangle.
bool intersect_triangle_bfc(decimate::math::Vec3<float>& orig, decimate::math::Vec3<float>& dir,
                            decimate::math::Vec3<float>& v0, decimate::math::Vec3<float>& v1, decimate::math::Vec3<float>& v2,
                            float* t, float* u, float* v);

// ---------------------------------------------------------------------------
// intersect_triangle  (backface culling disabled)
// ---------------------------------------------------------------------------
/// Möller–Trumbore ray/triangle intersection — backface culling OFF.
///
/// Same algorithm as intersect_triangle_bfc but handles both positive and
/// negative determinants, so both front- and back-facing triangles register
/// as hits.
///
/// @param orig  Ray origin.
/// @param dir   Ray direction.
/// @param v0    Triangle vertex 0.
/// @param v1    Triangle vertex 1.
/// @param v2    Triangle vertex 2.
/// @param[out] t  Distance parameter along the ray.
/// @param[out] u  Barycentric u coordinate.
/// @param[out] v  Barycentric v coordinate.
/// @return true if the ray hits the triangle (either face).
bool intersect_triangle(decimate::math::Vec3<float>& orig, decimate::math::Vec3<float>& dir, decimate::math::Vec3<float>& v0,
                        decimate::math::Vec3<float>& v1, decimate::math::Vec3<float>& v2, float* t, float* u,
                        float* v);

// ---------------------------------------------------------------------------
// IntersectTriangle  (Drew's adaption from Real Time Rendering)
// ---------------------------------------------------------------------------
/// Alternative Möller–Trumbore variant attributed to Drew's adaptation from
/// "Real Time Rendering", optimized against Möller's raytri.c.
///
/// Reference: http://www.ce.chalmers.se/staff/tomasm/raytri/raytri.c
///
/// Semantically equivalent to intersect_triangle (no backface culling).
/// Preserved as a separate entry point because downstream code may call
/// either variant.
///
/// @param orig  Ray origin.
/// @param dir   Ray direction.
/// @param v1    Triangle vertex 0 (note: uses v1/v2/v3 naming in legacy).
/// @param v2    Triangle vertex 1.
/// @param v3    Triangle vertex 2.
/// @param[out] t  Distance parameter along the ray.
/// @param[out] u  Barycentric u coordinate.
/// @param[out] v  Barycentric v coordinate.
/// @return true if the ray hits the triangle.
bool IntersectTriangle(decimate::math::Vec3<float>& orig, decimate::math::Vec3<float>& dir, decimate::math::Vec3<float>& v1,
                       decimate::math::Vec3<float>& v2, decimate::math::Vec3<float>& v3, float* t, float* u,
                       float* v);

}  // namespace rtrace
