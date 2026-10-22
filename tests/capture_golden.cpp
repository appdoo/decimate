// This file is part of the decimate library; see decimate.h for version/license details.
// Golden fixture capture for legacy rtrace functions.
//
// Compiles the legacy raytracing code (post-unblock, on clang) and runs
// it over the fixed input set defined in Spec 04, writing the results to
// tests/fixtures/raytracing_golden.txt.
//
// Run once to (re)generate the golden; the characterization test then
// asserts the core port against this fixture.
//
// To regenerate:
//   cmake --build build --target capture_golden
//   build/tests/capture_golden <output_path>

// Pull in legacy raytracing directly (not via mesh-core).
// This is the ONLY consumer of legacy code — it is not linked into mesh-core.
#include "../legacy2003/raytracing/raytracing.h"

namespace rtrace {
// The legacy .cpp defines barycentricBadouelf (with 'f' suffix); the header
// declares barycentricBadouel (without 'f'). The declaration was never matched
// by a definition in the legacy code — only the 'f' variant is live.
// Declare it here so the capture can call the actual defined symbol.
bool barycentricBadouelf(float &alpha, float &beta, const float *p, const float *nrm,
                         const float *p0, const float *p1, const float *p2);
} // namespace rtrace

#include <cstdio>
#include <cstring>

/// Builds a legacy Vec3 from three components.
static Vec3 v3(float x, float y, float z) {
    Vec3 r {};
    r.x = x;
    r.y = y;
    r.z = z;
    return r;
}


/// Writes the raytracing golden fixture (argv[1], default tests/fixtures/raytracing_golden.txt).
///
/// Fixed input set (Spec 04 §Architecture Step 2)
///
/// Triangle T = [v0(0,0,0), v1(1,0,0), v2(0,1,0)] in the z=0 plane.
/// Cases:
///   row 0: interior hit      — orig=(0.25,0.25,1), dir=(0,0,-1)
///   row 1: outside miss      — orig=(2,2,1),       dir=(0,0,-1)
///   row 2: parallel-to-plane — orig=(0.25,0.25,1), dir=(1,0,0)
///   row 3: origin behind     — orig=(0.25,0.25,-1),dir=(0,0,-1)
///   row 4: backface (bfc)    — orig=(0.25,0.25,-1),dir=(0,0,1)
int main(int argc, char *argv[]) {
    const char *outPath = (argc > 1) ? argv[1] : "tests/fixtures/raytracing_golden.txt";

    FILE *fp = std::fopen(outPath, "w");
    if (!fp) {
        std::fprintf(stderr, "Cannot open output file: %s\n", outPath);
        return 1;
    }

    // Unit right triangle T in the z=0 plane.
    Vec3 v0 = v3(0.0F, 0.0F, 0.0F);
    Vec3 v1 = v3(1.0F, 0.0F, 0.0F);
    Vec3 v2 = v3(0.0F, 1.0F, 0.0F);

    // evaluate_t
    // Row 0: standard hit — N=(0,0,1), d=0, O=(0,0,1), D=(0,0,-1) → t=1
    {
        Vec3  nrm = v3(0.0F, 0.0F, 1.0F);
        Vec3  orig = v3(0.0F, 0.0F, 1.0F);
        Vec3  dir = v3(0.0F, 0.0F, -1.0F);
        float t = 0.0F;
        bool  hit = rtrace::evaluate_t(t, nrm, 0.0F, orig, dir);
        std::fprintf(fp, "evaluate_t hit=%d t=%.9g\n", hit ? 1 : 0, static_cast<double>(t));
    }

    // Row 1: parallel — N=(0,0,1), d=0, O=(0,0,1), D=(1,0,0) → no hit
    {
        Vec3  nrm = v3(0.0F, 0.0F, 1.0F);
        Vec3  orig = v3(0.0F, 0.0F, 1.0F);
        Vec3  dir = v3(1.0F, 0.0F, 0.0F);
        float t = 0.0F;
        bool  hit = rtrace::evaluate_t(t, nrm, 0.0F, orig, dir);
        std::fprintf(fp, "evaluate_t_parallel hit=%d t=%.9g\n", hit ? 1 : 0,
                     static_cast<double>(t));
    }

    // intersect_triangle (no backface culling)
    // Row 0: interior hit
    {
        Vec3  orig = v3(0.25F, 0.25F, 1.0F);
        Vec3  dir = v3(0.0F, 0.0F, -1.0F);
        float t = 0.0F;
        float u = 0.0F;
        float vv = 0.0F;
        bool  hit = rtrace::intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &vv);
        std::fprintf(fp, "intersect_triangle_interior hit=%d t=%.9g u=%.9g v=%.9g\n", hit ? 1 : 0,
                     static_cast<double>(t), static_cast<double>(u), static_cast<double>(vv));
    }

    // Row 1: outside miss (u+v > 1)
    {
        Vec3  orig = v3(2.0F, 2.0F, 1.0F);
        Vec3  dir = v3(0.0F, 0.0F, -1.0F);
        float t = 0.0F;
        float u = 0.0F;
        float vv = 0.0F;
        bool  hit = rtrace::intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &vv);
        std::fprintf(fp, "intersect_triangle_outside hit=%d t=%.9g u=%.9g v=%.9g\n", hit ? 1 : 0,
                     static_cast<double>(t), static_cast<double>(u), static_cast<double>(vv));
    }

    // Row 2: ray parallel to plane
    {
        Vec3  orig = v3(0.25F, 0.25F, 1.0F);
        Vec3  dir = v3(1.0F, 0.0F, 0.0F);
        float t = 0.0F;
        float u = 0.0F;
        float vv = 0.0F;
        bool  hit = rtrace::intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &vv);
        std::fprintf(fp, "intersect_triangle_parallel hit=%d t=%.9g u=%.9g v=%.9g\n",
                     hit ? 1 : 0, static_cast<double>(t), static_cast<double>(u),
                     static_cast<double>(vv));
    }

    // Row 3: ray origin behind triangle (t would be negative — still a hit)
    {
        Vec3  orig = v3(0.25F, 0.25F, -1.0F);
        Vec3  dir = v3(0.0F, 0.0F, -1.0F);
        float t = 0.0F;
        float u = 0.0F;
        float vv = 0.0F;
        bool  hit = rtrace::intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &vv);
        std::fprintf(fp, "intersect_triangle_behind hit=%d t=%.9g u=%.9g v=%.9g\n", hit ? 1 : 0,
                     static_cast<double>(t), static_cast<double>(u), static_cast<double>(vv));
    }

    // intersect_triangle_bfc (backface culling ON)
    // Row 4: backface hit — ray from below, hits back face
    {
        Vec3  orig = v3(0.25F, 0.25F, -1.0F);
        Vec3  dir = v3(0.0F, 0.0F, 1.0F);
        float t = 0.0F;
        float u = 0.0F;
        float vv = 0.0F;
        bool  hit = rtrace::intersect_triangle_bfc(orig, dir, v0, v1, v2, &t, &u, &vv);
        std::fprintf(fp, "intersect_triangle_bfc_backface hit=%d t=%.9g u=%.9g v=%.9g\n",
                     hit ? 1 : 0, static_cast<double>(t), static_cast<double>(u),
                     static_cast<double>(vv));
    }

    // Also front face hit for bfc
    {
        Vec3  orig = v3(0.25F, 0.25F, 1.0F);
        Vec3  dir = v3(0.0F, 0.0F, -1.0F);
        float t = 0.0F;
        float u = 0.0F;
        float vv = 0.0F;
        bool  hit = rtrace::intersect_triangle_bfc(orig, dir, v0, v1, v2, &t, &u, &vv);
        std::fprintf(fp, "intersect_triangle_bfc_front hit=%d t=%.9g u=%.9g v=%.9g\n",
                     hit ? 1 : 0, static_cast<double>(t), static_cast<double>(u),
                     static_cast<double>(vv));
    }

    // IntersectTriangle (Drew's adaptation)
    // Interior hit — same geometry as intersect_triangle row 0
    {
        Vec3  orig = v3(0.25F, 0.25F, 1.0F);
        Vec3  dir = v3(0.0F, 0.0F, -1.0F);
        float t = 0.0F;
        float u = 0.0F;
        float vv = 0.0F;
        bool  hit = rtrace::IntersectTriangle(orig, dir, v0, v1, v2, &t, &u, &vv);
        std::fprintf(fp, "IntersectTriangle_interior hit=%d t=%.9g u=%.9g v=%.9g\n", hit ? 1 : 0,
                     static_cast<double>(t), static_cast<double>(u), static_cast<double>(vv));
    }

    // barycentricBadouel
    // Interior point (0.25, 0.25, 0) on triangle T (z=0 plane)
    // Face normal = (0, 0, 1)
    {
        float p[3] = { 0.25F, 0.25F, 0.0F };
        float nrm[3] = { 0.0F, 0.0F, 1.0F };
        float pp0[3] = { 0.0F, 0.0F, 0.0F };
        float pp1[3] = { 1.0F, 0.0F, 0.0F };
        float pp2[3] = { 0.0F, 1.0F, 0.0F };
        float alpha = 0.0F;
        float beta = 0.0F;
        bool  hit = rtrace::barycentricBadouelf(alpha, beta, p, nrm, pp0, pp1, pp2);
        std::fprintf(fp, "barycentricBadouelf_interior hit=%d alpha=%.9g beta=%.9g\n", hit ? 1 : 0,
                     static_cast<double>(alpha), static_cast<double>(beta));
    }

    // Outside point (2, 0, 0) — should be a miss
    {
        float p[3] = { 2.0F, 0.0F, 0.0F };
        float nrm[3] = { 0.0F, 0.0F, 1.0F };
        float pp0[3] = { 0.0F, 0.0F, 0.0F };
        float pp1[3] = { 1.0F, 0.0F, 0.0F };
        float pp2[3] = { 0.0F, 1.0F, 0.0F };
        float alpha = 0.0F;
        float beta = 0.0F;
        bool  hit = rtrace::barycentricBadouelf(alpha, beta, p, nrm, pp0, pp1, pp2);
        std::fprintf(fp, "barycentricBadouelf_outside hit=%d alpha=%.9g beta=%.9g\n", hit ? 1 : 0,
                     static_cast<double>(alpha), static_cast<double>(beta));
    }

    std::fclose(fp);
    std::printf("Golden fixture written to: %s\n", outPath);
    return 0;
}
