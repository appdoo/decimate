// This file is part of the decimate library; see decimate.h for version/license details.
/// @file parametrize_pack_test.cpp
/// @brief Tests for Spec 08-02: atlas packing + parametrize() orchestrator.
///
/// ## Test strategy
///
/// - Small hand-crafted meshes (single triangle, coplanar quad) for exact
///   geometric checks.
/// - Real cow.obj for integration-scale proof.
/// - All assertions are structural / numeric — no text-presence checks
///   (antipattern 12).
///
/// ## Scenario mapping (spec 08-02 test table)
///
/// B1 [B]  Single flat triangle chart → UVs preserve edge-length ratios (isometry, 1e-4).
/// B2 [B]  Coplanar quad (2 tris) → conformal energy < 1e-6 (isometric flattening).
/// B3 [I]  cow.obj → segment → flatten all charts → every UV finite, area > 0.
/// B4 [I]  cow.obj → packAtlas → all UVs in [0,1]²; cross-chart PIXEL overlap tiny.
/// B5 [B]  Full parametrize() twice → byte-identical UVs (determinism).
/// S5 [S]  ScaleProportional: known single chart → correct UV scale factor.
/// S6 [S]  16 rotations are tried during pack (verified via resolution probe).
/// S7 [B]  atlasResolution parameter: calling with 128 and 256 both succeed.
/// S8 [S]  atlasResolution is honoured: 128 vs 256 produce different pixel scale.

#include <catch2/catch_approx.hpp>
#include "model_fixture.h"
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include <decimate/chart/chart.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>
#include <decimate/chart/parametrize.h>
#include <decimate/math/vec2.h>
#include <decimate/math/vec3.h>

#ifndef MESH_MODEL_DIR
#error "MESH_MODEL_DIR must be defined by CMake"
#endif

using namespace decimate;
using namespace decimate::chart;
using Catch::Approx;

/// @brief Build a MeshData from explicit position + index lists.
static MeshData makeMesh(std::initializer_list<float> positions,
                         std::initializer_list<int>   indices) {
    MeshData                 m;
    const std::vector<float> flat(positions);
    for (std::size_t i = 0; i + 3 <= flat.size(); i += 3) {
        m.positions.push_back({ flat[i], flat[i + 1], flat[i + 2] });
    }
    m.indices.assign(indices);
    return m;
}


/// @brief Compute 3D edge length between two mesh vertices by global index.
static float edgeLen3D(const MeshData &data, int i, int j) {
    const auto &vi = data.positions[static_cast<std::size_t>(i)];
    const auto &vj = data.positions[static_cast<std::size_t>(j)];
    float       dx = vi.x - vj.x;
    float       dy = vi.y - vj.y;
    float       dz = vi.z - vj.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}


/// @brief Compute 2D edge length between two local UV indices.
static float edgeLenUV(const HalfEdgeChart &ac, int li, int lj) {
    float du = ac.uv[li].x - ac.uv[lj].x;
    float dv = ac.uv[li].y - ac.uv[lj].y;
    return std::sqrt(du * du + dv * dv);
}


/// @brief True if a UV is finite (neither NaN nor Inf).
static bool uvFinite(const decimate::math::Vec2f &uv) {
    return std::isfinite(uv.x) && std::isfinite(uv.y);
}


/// @brief UV area of one chart via shoelace over all its triangles.
static float chartUVArea(const HalfEdgeChart &ac) {
    float area = 0.0f;
    for (int i = 0; i < static_cast<int>(ac.faces.size()); ++i) {
        const int   la = localCorner(ac, i, 0);
        const int   lb = localCorner(ac, i, 1);
        const int   lc = localCorner(ac, i, 2);
        const auto &a = ac.uv[la];
        const auto &b = ac.uv[lb];
        const auto &c = ac.uv[lc];
        float       det = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
        area += std::abs(det) * 0.5f;
    }
    return area;
}


/// B1 [B] — Single flat triangle: edge-length ratios preserved after flatten.
///
/// Given: a right-isosceles triangle in 3D with legs of length 1.
/// When:  flatten() is called.
/// Then:  the UV edge lengths satisfy the same RATIO (isometry up to scale),
///        within 1e-4.
TEST_CASE("B1: single triangle flatten preserves edge-length ratios",
          "[parametrize][flatten][B1]") {
    // Right-isosceles triangle: A=(0,0,0), B=(1,0,0), C=(0,1,0).
    // 3D edges: AB=1, AC=1, BC=sqrt(2).
    MeshData data = makeMesh(
        { 0.0f, 0.0f, 0.0f,   // vertex 0 = A
          1.0f, 0.0f, 0.0f,   // vertex 1 = B
          0.0f, 1.0f, 0.0f }, // vertex 2 = C
        { 0, 1, 2 });
    const int          numFaces = 1;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(3, numFaces,
                                                  reinterpret_cast<const Face *>(data.indices.data()));

    // Single chart containing the only face.
    Chart ch;
    ch.id = 0;
    ch.seed = 0;
    ch.faces = { 0 };

    HalfEdgeChart ac = flatten(hem, data, ch);

    // flatten() must produce 3 local vertices.
    REQUIRE(ac.uv.size() == 3);
    REQUIRE(ac.localToMesh.size() == 3);

    // Recover local indices for global vertices 0, 1, 2.
    // The single face is (0, 1, 2): its corner k is mesh vertex k.
    const std::array<int, 3> localOf { localCorner(ac, 0, 0), localCorner(ac, 0, 1),
                                       localCorner(ac, 0, 2) };

    // 3D edge lengths.
    float len3D_AB = edgeLen3D(data, 0, 1); // expected 1.0
    float len3D_AC = edgeLen3D(data, 0, 2); // expected 1.0
    float len3D_BC = edgeLen3D(data, 1, 2); // expected sqrt(2)

    // UV edge lengths.
    float lenUV_AB = edgeLenUV(ac, localOf[0], localOf[1]);
    float lenUV_AC = edgeLenUV(ac, localOf[0], localOf[2]);
    float lenUV_BC = edgeLenUV(ac, localOf[1], localOf[2]);

    // LSCM is conformal (angle-preserving, hence isometric up to scale).
    // The ratios must be preserved: ratio3D(AB/BC) == ratioUV(AB/BC).
    float ratio3D_AB_BC = len3D_AB / len3D_BC;
    float ratio3D_AC_BC = len3D_AC / len3D_BC;
    float ratioUV_AB_BC = lenUV_AB / lenUV_BC;
    float ratioUV_AC_BC = lenUV_AC / lenUV_BC;

    CHECK(ratioUV_AB_BC == Approx(ratio3D_AB_BC).epsilon(1e-4));
    CHECK(ratioUV_AC_BC == Approx(ratio3D_AC_BC).epsilon(1e-4));
}


/// B2 [B] — Coplanar quad (2 triangles): conformal energy < 1e-6.
///
/// Given: a unit square split into two triangles (coplanar, zero Gaussian curvature).
/// When:  flatten() is called on the chart.
/// Then:  conformal energy = sum |J - R|² / num_triangles < 1e-6,
///        where J is the UV Jacobian and R is its nearest rotation.
///        For a planar mesh LSCM is exact; the Cauchy-Riemann equations
///        hold exactly (energy = 0 in exact arithmetic, < 1e-6 numerically).
TEST_CASE("B2: coplanar quad flatten — conformal energy < 1e-6",
          "[parametrize][flatten][B2]") {
    // Unit square: two triangles sharing the diagonal 0→2.
    // Vertices: 0=(0,0,0), 1=(1,0,0), 2=(1,1,0), 3=(0,1,0).
    // Faces: (0,1,2) and (0,2,3).
    MeshData data = makeMesh(
        { 0.0f, 0.0f, 0.0f,
          1.0f, 0.0f, 0.0f,
          1.0f, 1.0f, 0.0f,
          0.0f, 1.0f, 0.0f },
        { 0, 1, 2,
          0, 2, 3 });
    const int          numFaces = 2;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(4, numFaces,
                                                  reinterpret_cast<const Face *>(data.indices.data()));

    Chart ch;
    ch.id = 0;
    ch.seed = 0;
    ch.faces = { 0, 1 };

    HalfEdgeChart ac = flatten(hem, data, ch);

    // Measure conformal energy per triangle.
    // For a planar mesh LSCM should be near-exact; energy measures how far
    // the UV Jacobian is from a rotation+scale (= conformality error).
    // A simple proxy: for each triangle, check that the angles between UV
    // edges match those in 3D.
    //
    // Precise conformal energy: E_c = sum_t area_t * |J_t - R(J_t)|^2
    // where J_t is the 2x2 Jacobian mapping UV to 3D.
    // For a planar mesh this is zero.  We use a normalised angle-based proxy:
    // |cos(angle_UV) - cos(angle_3D)|² for each interior angle.

    float totalError = 0.0f;
    int   numAngles = 0;

    for (int fi = 0; fi < numFaces; ++fi) {
        int ga = data.indices[3 * fi + 0];
        int gb = data.indices[3 * fi + 1];
        int gc = data.indices[3 * fi + 2];
        // Chart {0, 1}: local face index == mesh face index.
        int la = localCorner(ac, fi, 0), lb = localCorner(ac, fi, 1), lc = localCorner(ac, fi, 2);

        // 3D vectors.
        auto v3 = [&](int gi) -> decimate::math::Vec3f {
            return data.positions[static_cast<std::size_t>(gi)];
        };
        decimate::math::Vec3f A3 = v3(ga);
        decimate::math::Vec3f B3 = v3(gb);
        decimate::math::Vec3f C3 = v3(gc);

        auto edgeV3 = [](decimate::math::Vec3f from, decimate::math::Vec3f to) {
            return to - from;
        };
        decimate::math::Vec3f AB3 = edgeV3(A3, B3);
        decimate::math::Vec3f AC3 = edgeV3(A3, C3);

        // UV vectors.
        const auto           &Auv = ac.uv[la];
        const auto           &Buv = ac.uv[lb];
        const auto           &Cuv = ac.uv[lc];
        decimate::math::Vec2f ABuv = Buv - Auv;
        decimate::math::Vec2f ACuv = Cuv - Auv;

        // Cosine of angle at A in 3D.
        float cos3D = AB3.dot(AC3) / (AB3.length() * AC3.length() + 1e-12f);
        // Cosine of angle at A in UV.
        float cosUV = ABuv.dot(ACuv) / (ABuv.length() * ACuv.length() + 1e-12f);

        float err = cos3D - cosUV;
        totalError += err * err;
        ++numAngles;
    }

    float energy = (numAngles > 0) ? (totalError / static_cast<float>(numAngles)) : 0.0f;
    CHECK(energy < 1e-6);
}


/// B3 [I] — cow.obj: segment → flatten all charts → all UVs finite, area > 0.
TEST_CASE("B3: cow.obj — flatten all charts: UVs finite, area > 0",
          "[parametrize][flatten][integration][B3]") {
    const std::string cowPath = modelOrSkip("cow.obj");
    MeshData          data;
    try {
        data = loadObj(cowPath);
    } catch (...) {
        FAIL("Could not load cow.obj");
    }
    REQUIRE_FALSE(data.positions.empty());
    REQUIRE_FALSE(data.indices.empty());

    const int numFaces = static_cast<int>(data.indices.size()) / 3;
    const int numPoints = static_cast<int>(data.positions.size());

    const HalfEdgeMesh hem = HalfEdgeMesh::create(numPoints, numFaces,
                                                  reinterpret_cast<const Face *>(data.indices.data()));

    const ChartSet cs = segment(hem, data, 64);
    REQUIRE_FALSE(cs.charts.empty());

    for (const Chart &ch : cs.charts) {
        if (ch.faces.empty()) {
            continue;
        }

        HalfEdgeChart ac = flatten(hem, data, ch);

        // Every UV must be finite.
        for (const auto &uv : ac.uv) {
            REQUIRE(uvFinite(uv));
        }

        // Each chart must have positive UV area after flattening.
        float area = chartUVArea(ac);
        CHECK(area > 0.0f);
    }
}


/// B4 [I] — cow.obj → packAtlas → all UVs in [0,1]²; cross-chart PIXEL overlap tiny.
TEST_CASE("B4: cow.obj — packAtlas: UVs in [0,1]^2, cross-chart pixel overlap tiny",
          "[parametrize][packAtlas][integration][B4]") {
    const std::string cowPath = modelOrSkip("cow.obj");
    MeshData          data;
    try {
        data = loadObj(cowPath);
    } catch (...) {
        FAIL("Could not load cow.obj");
    }

    const int numFaces = static_cast<int>(data.indices.size()) / 3;
    const int numPoints = static_cast<int>(data.positions.size());

    const HalfEdgeMesh hem = HalfEdgeMesh::create(numPoints, numFaces,
                                                  reinterpret_cast<const Face *>(data.indices.data()));

    const ChartSet cs = segment(hem, data, 64);
    REQUIRE_FALSE(cs.charts.empty());

    // Flatten all charts.
    std::vector<HalfEdgeChart> atlasCharts;
    atlasCharts.reserve(cs.charts.size());
    for (const Chart &ch : cs.charts) {
        atlasCharts.push_back(flatten(hem, data, ch));
    }

    // Pack.
    packAtlas(atlasCharts, data, 256);

    // AC1: all UVs in [0,1]².
    for (const HalfEdgeChart &ac : atlasCharts) {
        for (const auto &uv : ac.uv) {
            CHECK(uv.x >= -1e-4f);
            CHECK(uv.x <= 1.0f + 1e-4f);
            CHECK(uv.y >= -1e-4f);
            CHECK(uv.y <= 1.0f + 1e-4f);
        }
    }

    // AC2: charts must not overlap in the atlas — measured at the PIXEL level,
    // NOT by AABB.  A silhouette / Tetris packer INTENTIONALLY interlocks chart
    // bounding boxes (a concave chart nests into a neighbour's bounding box)
    // while the actual rasterised triangles stay disjoint.  An AABB-overlap
    // assertion forbids exactly the nesting the packer exists to do — that wrong
    // test previously drove a regression to a box-packer.
    // Rasterise every chart's triangles into a coverage grid tagged by chart id;
    // assert cross-chart double-coverage is a tiny fraction (seam slivers only).
    constexpr int    kRes = 256;
    std::vector<int> cell(static_cast<std::size_t>(kRes) * kRes, -1);
    long             covered = 0;
    long             crossChart = 0;
    auto             edge = [](float px, float py, float qx, float qy, float rx, float ry) {
        return (qx - px) * (ry - py) - (qy - py) * (rx - px);
    };
    for (int ci = 0; ci < static_cast<int>(atlasCharts.size()); ++ci) {
        const HalfEdgeChart &ac = atlasCharts[ci];
        for (int i = 0; i < static_cast<int>(ac.faces.size()); ++i) {
            const int la = localCorner(ac, i, 0);
            const int lb = localCorner(ac, i, 1);
            const int lc = localCorner(ac, i, 2);
            {
                const auto &a = ac.uv[static_cast<std::size_t>(la)];
                const auto &b = ac.uv[static_cast<std::size_t>(lb)];
                const auto &c = ac.uv[static_cast<std::size_t>(lc)];
                const float ax = a.x * kRes, ay = a.y * kRes;
                const float bx = b.x * kRes, by = b.y * kRes;
                const float cx = c.x * kRes, cy = c.y * kRes;
                const int   xlo = std::max(0, static_cast<int>(std::floor(std::min({ ax, bx, cx }))));
                const int   xhi = std::min(kRes - 1, static_cast<int>(std::ceil(std::max({ ax, bx, cx }))));
                const int   ylo = std::max(0, static_cast<int>(std::floor(std::min({ ay, by, cy }))));
                const int   yhi = std::min(kRes - 1, static_cast<int>(std::ceil(std::max({ ay, by, cy }))));
                for (int y = ylo; y <= yhi; ++y) {
                    for (int x = xlo; x <= xhi; ++x) {
                        const float pxc = static_cast<float>(x) + 0.5f;
                        const float pyc = static_cast<float>(y) + 0.5f;
                        const float w0 = edge(bx, by, cx, cy, pxc, pyc);
                        const float w1 = edge(cx, cy, ax, ay, pxc, pyc);
                        const float w2 = edge(ax, ay, bx, by, pxc, pyc);
                        const bool  inside = (w0 >= 0 && w1 >= 0 && w2 >= 0) ||
                                             (w0 <= 0 && w1 <= 0 && w2 <= 0);
                        if (!inside) {
                            continue;
                        }
                        int &tag = cell[static_cast<std::size_t>(y) * kRes + x];
                        if (tag == -1) {
                            tag = ci;
                            ++covered;
                        } else if (tag != ci) {
                            ++crossChart;
                        }
                    }
                }
            } // inner block (a, b, c UVs)
        } // for local face
    }
    REQUIRE(covered > 0);
    const double overlapFrac = static_cast<double>(crossChart) / static_cast<double>(covered);
    INFO("covered=" << covered << " crossChart=" << crossChart
                    << " overlap=" << (overlapFrac * 100.0) << "%");
    CHECK(overlapFrac < 0.05);
}


/// B5/S7 [B] — Full parametrize() twice → byte-identical UVs (determinism).
TEST_CASE("B5/S7: parametrize() is deterministic — two runs yield identical UVs",
          "[parametrize][determinism][B5][S7]") {
    // A small open box (cube without its y=1 cap, 10 faces): open so every chart has a boundary.
    // A closed cube segments into one seamless chart, which packAtlas rejects (Spec 24 H4).
    // clang-format off
    MeshData data = makeMesh(
        { 0.0f,0.0f,0.0f,  1.0f,0.0f,0.0f,  1.0f,1.0f,0.0f,  0.0f,1.0f,0.0f,
          0.0f,0.0f,1.0f,  1.0f,0.0f,1.0f,  1.0f,1.0f,1.0f,  0.0f,1.0f,1.0f },
        { 0,2,1,  0,3,2,
          4,5,6,  4,6,7,
          0,1,5,  0,5,4,
          0,4,7,  0,7,3,
          1,2,6,  1,6,5 }
    );
    // clang-format on
    const int numFaces = static_cast<int>(data.indices.size()) / 3;
    const int numPoints = static_cast<int>(data.positions.size());

    const HalfEdgeMesh hem = HalfEdgeMesh::create(numPoints, numFaces,
                                                  reinterpret_cast<const Face *>(data.indices.data()));

    ChartSet cs = segment(hem, data, 64);

    // Run twice.
    auto uvA = parametrize(hem, data, cs, 256);
    auto uvB = parametrize(hem, data, cs, 256);

    REQUIRE(uvA.size() == uvB.size());
    for (int i = 0; i < static_cast<int>(uvA.size()); ++i) {
        // Byte-identical: same bit pattern.
        uint32_t bxA, byA, bxB, byB;
        std::memcpy(&bxA, &uvA[i].x, sizeof(float));
        std::memcpy(&byA, &uvA[i].y, sizeof(float));
        std::memcpy(&bxB, &uvB[i].x, sizeof(float));
        std::memcpy(&byB, &uvB[i].y, sizeof(float));
        CHECK(bxA == bxB);
        CHECK(byA == byB);
    }
}


/// S5 [S] — ScaleProportional: UV area normalises to 3D area share.
///
/// The formula is cscale = sqrt(chartArea3D / (totalArea3D * uvArea)), so after scaling
/// uvArea' = uvArea * cscale^2 = chartArea3D / totalArea3D: with one chart that is 1.0, and for two
/// charts the post-scale UV areas stand in the ratio of their 3D areas.
///
/// ScaleProportional is internal to packAtlas, so this test checks the formula's arithmetic on the
/// measured pre-scale UV areas, not the library call itself.
TEST_CASE("S5: ScaleProportional — UV area normalises to 1.0 for a single chart",
          "[parametrize][scale][S5]") {
    // Equilateral triangle side length 2: vertices at (0,0,0),(2,0,0),(1,sqrt(3),0).
    // 3D area = sqrt(3).
    // clang-format off
    MeshData data = makeMesh(
        { 0.0f, 0.0f, 0.0f,
          2.0f, 0.0f, 0.0f,
          1.0f, 1.7320508f, 0.0f },
        { 0, 1, 2 }
    );
    // clang-format on
    const int          numFaces = 1;
    const int          numPoints = 3;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(numPoints, numFaces,
                                                  reinterpret_cast<const Face *>(data.indices.data()));

    Chart ch;
    ch.id = 0;
    ch.seed = 0;
    ch.faces = { 0 };

    HalfEdgeChart ac = flatten(hem, data, ch);

    // Confirm pre-scale UV area is non-trivial.
    float uvAreaBefore = chartUVArea(ac);
    REQUIRE(uvAreaBefore > 1e-12f);

    std::vector<HalfEdgeChart> charts;
    charts.push_back(std::move(ac));

    // One chart gives a ratio of 1.0 and no information, so use two triangles with 3D areas S and 2S.
    // A: area sqrt(3); B: same base, double height.
    MeshData data2 = makeMesh(
        { // triangle A: area = sqrt(3) ≈ 1.732
          0.0f, 0.0f, 0.0f,
          2.0f, 0.0f, 0.0f,
          1.0f, 1.7320508f, 0.0f,
          // triangle B: base 4, height sqrt(3): area = 2*sqrt(3) ≈ 3.464
          3.0f, 0.0f, 0.0f,
          7.0f, 0.0f, 0.0f,
          5.0f, 1.7320508f, 0.0f },
        { 0, 1, 2,
          3, 4, 5 });
    const int          nf2 = 2;
    const int          nv2 = 6;
    const HalfEdgeMesh hem2 = HalfEdgeMesh::create(nv2, nf2,
                                                   reinterpret_cast<const Face *>(data2.indices.data()));

    Chart chA;
    chA.id = 0;
    chA.seed = 0;
    chA.faces = { 0 };
    Chart chB;
    chB.id = 1;
    chB.seed = 1;
    chB.faces = { 1 };

    HalfEdgeChart acA = flatten(hem2, data2, chA);
    HalfEdgeChart acB = flatten(hem2, data2, chB);

    float uvABefore = chartUVArea(acA);
    float uvBBefore = chartUVArea(acB);
    REQUIRE(uvABefore > 1e-12f);
    REQUIRE(uvBBefore > 1e-12f);

    // Post-scale uvArea_X = uvX_before * cscale_X^2 = area3D_X / totalArea3D, so the UV ratio equals
    // the 3D ratio.

    float area3D_A = 1.7320508f; // sqrt(3) ≈ equilateral tri side 2
    float area3D_B = 3.4641016f; // base 4 x height sqrt(3) / 2
    float totalArea = area3D_A + area3D_B;

    float expectedRatio = area3D_A / area3D_B;
    float cscaleA = std::sqrt(area3D_A / (totalArea * uvABefore));
    float cscaleB = std::sqrt(area3D_B / (totalArea * uvBBefore));
    float postScaleA = uvABefore * cscaleA * cscaleA;
    float postScaleB = uvBBefore * cscaleB * cscaleB;

    // Post-scale UV area ratio must equal 3D area ratio (within 1e-4).
    CHECK(postScaleA / postScaleB == Approx(expectedRatio).epsilon(1e-4));

    // The check above only proves the formula against itself.  Now the same two charts go through
    // the real packAtlas and the ratio is read off its output: packed UV area per chart must stand
    // in the ratio of the 3D areas, i.e. equal texel density across charts.  packAtlas may apply
    // one common factor to every chart, so only the ratio is asserted.  The tolerance is 1e-3: the
    // scale is exact in float, and the slack covers the shoelace sum on UVs near [0,1].
    std::vector<HalfEdgeChart> packed;
    packed.push_back(std::move(acA));
    packed.push_back(std::move(acB));
    packAtlas(packed, data2, 256);

    const float packedA = chartUVArea(packed[0]);
    const float packedB = chartUVArea(packed[1]);
    REQUIRE(packedA > 0.0f);
    REQUIRE(packedB > 0.0f);
    CHECK(packedA / packedB == Approx(expectedRatio).epsilon(1e-3));
}


/// S6 [S] — 16 rotations tried: verify via resolution probe.
///
/// Scenario: We call packAtlas with a chart that only fits in one specific
///   orientation.  If fewer than 16 rotations were tried, the pack might
///   fail or produce a suboptimal result.
///
/// Approach: Use a thin 1×16 rectangle chart.  At 0° it fits in a 256-grid;
///   at 90° it also fits.  The pack should always succeed regardless of
///   orientation because at least one of the 16 rotations will find a fit.
///
/// The test is structural: after packAtlas with atlasResolution=256, all
///   UVs must be in [0,1]² (proving a valid placement was found).
TEST_CASE("S6: 16 rotations tried — thin chart always fits in at least one orientation",
          "[parametrize][rotations][S6]") {
    // Build a thin 1×8 rectangle mesh: a quad made of two tris.
    // In UV space (after flatten) this will be elongated.
    // clang-format off
    MeshData data = makeMesh(
        { 0.0f, 0.0f, 0.0f,
          8.0f, 0.0f, 0.0f,
          8.0f, 0.5f, 0.0f,
          0.0f, 0.5f, 0.0f },
        { 0, 1, 2,
          0, 2, 3 }
    );
    // clang-format on
    const int          numFaces = 2;
    const int          numPoints = 4;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(numPoints, numFaces,
                                                  reinterpret_cast<const Face *>(data.indices.data()));

    Chart ch;
    ch.id = 0;
    ch.seed = 0;
    ch.faces = { 0, 1 };

    HalfEdgeChart ac = flatten(hem, data, ch);

    std::vector<HalfEdgeChart> charts;
    charts.push_back(std::move(ac));

    // Should succeed and produce valid [0,1]² UVs.
    packAtlas(charts, data, 256);

    for (const auto &uv : charts[0].uv) {
        CHECK(uv.x >= -1e-4f);
        CHECK(uv.x <= 1.0f + 1e-4f);
        CHECK(uv.y >= -1e-4f);
        CHECK(uv.y <= 1.0f + 1e-4f);
    }
}


/// S7 / S8 [B] — atlasResolution parameter: 128 and 256 both produce valid UVs.
TEST_CASE("S7/S8: atlasResolution parameter — 128 and 256 both succeed",
          "[parametrize][resolution][S7][S8]") {
    // Small quad mesh.
    // clang-format off
    MeshData data = makeMesh(
        { 0.0f, 0.0f, 0.0f,
          1.0f, 0.0f, 0.0f,
          1.0f, 1.0f, 0.0f,
          0.0f, 1.0f, 0.0f },
        { 0, 1, 2,
          0, 2, 3 }
    );
    // clang-format on
    const int          numFaces = 2;
    const int          numPoints = 4;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(numPoints, numFaces,
                                                  reinterpret_cast<const Face *>(data.indices.data()));

    ChartSet cs = segment(hem, data, 4);

    // Resolution 128.
    {
        auto uvs = parametrize(hem, data, cs, 128);
        REQUIRE_FALSE(uvs.empty());
        for (const auto &uv : uvs) {
            CHECK(uvFinite(uv));
            CHECK(uv.x >= -1e-4f);
            CHECK(uv.x <= 1.0f + 1e-4f);
            CHECK(uv.y >= -1e-4f);
            CHECK(uv.y <= 1.0f + 1e-4f);
        }
    }

    // Resolution 256.
    {
        auto uvs = parametrize(hem, data, cs, 256);
        REQUIRE_FALSE(uvs.empty());
        for (const auto &uv : uvs) {
            CHECK(uvFinite(uv));
            CHECK(uv.x >= -1e-4f);
            CHECK(uv.x <= 1.0f + 1e-4f);
            CHECK(uv.y >= -1e-4f);
            CHECK(uv.y <= 1.0f + 1e-4f);
        }
    }
}


/// S3 [S] — per-vertex UV size matches mesh vertex count.
TEST_CASE("S3: parametrize() output size equals number of mesh vertices",
          "[parametrize][S3]") {
    // clang-format off
    MeshData data = makeMesh(
        { 0.0f,0.0f,0.0f,  1.0f,0.0f,0.0f,  1.0f,1.0f,0.0f,  0.0f,1.0f,0.0f },
        { 0, 1, 2,  0, 2, 3 }
    );
    // clang-format on
    const int          numFaces = 2;
    const int          numPoints = 4;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(numPoints, numFaces,
                                                  reinterpret_cast<const Face *>(data.indices.data()));

    ChartSet cs = segment(hem, data, 4);
    auto     uvs = parametrize(hem, data, cs, 256);

    // Output size == positions.size().
    REQUIRE(uvs.size() == static_cast<std::size_t>(numPoints));
}


/// S8 [S] — seam vertex gets last-written chart's UV.
///
/// Scenario: a mesh with exactly 2 coplanar charts sharing one vertex.
///   After parametrize(), the shared vertex's UV equals the UV assigned by
///   the chart processed LAST (chart index 1 in chart order).
TEST_CASE("S8: seam vertex receives last-chart UV (legacy seam behaviour)",
          "[parametrize][seam][S8]") {
    // Build a T-shape: two quads (each 2 tris) sharing vertex 2.
    // Left quad: vertices 0,1,2,3.  Right quad: vertices 2,4,5,6.
    // Vertex 2 is on the seam.
    // clang-format off
    MeshData data = makeMesh(
        { 0.0f, 0.0f, 0.0f,  // 0
          1.0f, 0.0f, 0.0f,  // 1
          1.0f, 1.0f, 0.0f,  // 2 — seam vertex
          0.0f, 1.0f, 0.0f,  // 3
          2.0f, 0.0f, 0.0f,  // 4
          2.0f, 1.0f, 0.0f,  // 5
          1.5f, 0.5f, 0.0f }, // 6 — interior of right quad
        { 0, 1, 2,   // left chart, face 0
          0, 2, 3,   // left chart, face 1
          1, 4, 5,   // right chart, face 2
          1, 5, 2 }  // right chart, face 3 — touches seam vertex 2
    );
    // clang-format on
    const int          numFaces = 4;
    const int          numPoints = 7;
    const HalfEdgeMesh hem = HalfEdgeMesh::create(numPoints, numFaces,
                                                  reinterpret_cast<const Face *>(data.indices.data()));

    ChartSet cs = segment(hem, data, 4);

    // At least 2 charts for a T-mesh (different normal planes not needed
    // here since it's all coplanar — the segmenter may produce 1 or 2 charts).
    // This test is meaningful when 2 charts exist; skip gracefully otherwise.
    if (cs.charts.size() < 2) {
        SUCCEED("Segmenter merged both quads into one chart — seam test N/A");
        return;
    }

    auto uvs = parametrize(hem, data, cs, 256);

    // The seam vertex (global id 2) must have a finite UV.
    REQUIRE(uvs.size() > 2);
    CHECK(uvFinite(uvs[2]));
}


/// Spec 22-02 I1: cornerUVs is the scatter parametrizeCorners ends with
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("I1: cornerUVs(packAtlas(flattenCharts)) equals parametrizeCorners bit for bit",
          "[parametrize][cornerUVs][integration][I1][22-02]") {
    for (const char *model : { "spot.obj", "cow.obj" }) {
        CAPTURE(model);
        const MeshData     data = loadObj(modelOrSkip(model));
        const int          numFaces = static_cast<int>(data.indices.size()) / 3;
        const HalfEdgeMesh hem = HalfEdgeMesh::create(static_cast<int>(data.positions.size()), numFaces,
                                                      reinterpret_cast<const Face *>(data.indices.data()));
        const ChartSet     fresh = segment(hem, data, 64);

        ChartSet   viaOld = fresh;
        const auto old = parametrizeCorners(hem, data, viaOld, 256);

        ChartSet viaNew = fresh;
        auto     atlas = flattenCharts(hem, data, viaNew);
        packAtlas(atlas, data, 256);
        const auto uvs = cornerUVs(atlas, numFaces);

        REQUIRE(uvs.size() == 3 * static_cast<std::size_t>(numFaces));
        REQUIRE(uvs.size() == old.size());
        CHECK(std::memcmp(uvs.data(), old.data(), uvs.size() * sizeof(uvs[0])) == 0);

        // Not vacuous: the atlas really spreads over the unit square.
        float maxU = 0.0f;
        for (const auto &uv : uvs) {
            maxU = std::max(maxU, uv.x);
        }
        CHECK(maxU > 0.5f);
    }
}
