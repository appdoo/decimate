// This file is part of the decimate library; see decimate.h for version/license details.
/// @file seam_bake_test.cpp
/// @brief Unit tests for bakeSeamFill — seam-correct gutter fill.
///
/// ## Strategy
///
/// Tests use SYNTHETIC minimal fixtures (no model files) so the expected values
/// are analytically derivable. No golden files needed.
///
/// ## Test scenarios
///
/// [S/B] Two-chart seam: two quads (two triangles each) placed at different atlas
///       positions. The shared inner edge is a CHART SEAM. A gutter texel just
///       outside chart A's seam edge MUST be filled from chart B's interior, NOT
///       from chart A's own nearest filled texel (which is the blur-dilate result).
///
/// [S]   Open-mesh boundary fallback: a single chart with an open boundary edge
///       (no sym). Gutter texels along the open boundary must not crash and must
///       be filled by the blur fallback (non-sentinel value after dilatePasses>0).
///
/// [S]   dilatePasses==0 faithful bake: bake() with dilatePasses=0 returns the
///       identical result regardless of whether bakeSeamFill exists (the
///       existing golden invariant is unaffected).

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdint>
#include <vector>

#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normal_map.h>

namespace {

/// Background sentinel (mirrors the bake's kBg).
static constexpr decimate::bake::Rgb kBg { 128, 128, 128 };

/// Build a Face array from an indices vector (3 ints per face).
static std::vector<decimate::Face> makeFaces(const std::vector<int> &indices) {
    std::vector<decimate::Face> faces;
    for (std::size_t i = 0; i + 2 < indices.size(); i += 3) {
        faces.push_back({ indices[i], indices[i + 1], indices[i + 2] });
    }
    return faces;
}


/// Fixture: two coplanar quads (two triangles each, XZ plane, Y=0) sharing the mesh edge v1-v4.
///
///   v0 (0,0,0) --- v1 (1,0,0) --- v2 (2,0,0)
///       |     \        |     \        |
///   v3 (0,0,1) --- v4 (1,0,1) --- v5 (2,0,1)
///
/// Chart A (left quad) = faces 0,1; chart B (right quad) = faces 2,3.  The shared edge is a CHART
/// SEAM because the two charts place v1 and v4 at different atlas positions:
///   A: v1 -> (0.45, 0), v4 -> (0.45, 1)   (A covers atlas x in [0, 0.45])
///   B: v1 -> (0.55, 0), v4 -> (0.55, 1)   (B covers atlas x in [0.55, 1])
/// so cornerUV diverges across the mesh edge and the gap x in (0.45, 0.55) is gutter.
///
/// Every normal is (0,1,0), so every baked interior texel is rangeCompress((0,1,0)) =
/// (127, 254, 127).  A gutter texel in the gap filled from chart B's interior therefore has that
/// same value; what separates it from "not filled" is that it differs from the background kBg.
struct TwoChartFixture {
    decimate::MeshData                 low;
    decimate::HalfEdgeMesh             hem;
    std::vector<decimate::math::Vec2f> cornerUV; // 3 per face
    std::vector<int>                   faceToChart;
    std::vector<decimate::math::Vec3f> lowNormals;
    decimate::MeshData                 high;
    std::vector<decimate::math::Vec3f> highNormals;
    int                                texSize = 16;
};


static TwoChartFixture makeTwoChartFixture() {
    TwoChartFixture f;
    f.texSize = 16;

    // 6 vertices in the XZ plane (y=0)
    f.low.positions = {
        { 0.0f, 0.0f, 0.0f }, // v0
        { 1.0f, 0.0f, 0.0f }, // v1
        { 2.0f, 0.0f, 0.0f }, // v2
        { 0.0f, 0.0f, 1.0f }, // v3
        { 1.0f, 0.0f, 1.0f }, // v4
        { 2.0f, 0.0f, 1.0f }, // v5
    };

    // 4 triangles: chart A = faces 0,1; chart B = faces 2,3
    // Chart A (left quad): v0,v1,v4 and v0,v4,v3
    // Chart B (right quad): v1,v2,v5 and v1,v5,v4
    f.low.indices = {
        0,
        1,
        4, // face 0 (chart A)
        0,
        4,
        3, // face 1 (chart A)
        1,
        2,
        5, // face 2 (chart B)
        1,
        5,
        4, // face 3 (chart B)
    };

    const int nlv = 6;
    const int nlf = 4;
    auto      faces = makeFaces(f.low.indices);
    f.hem = decimate::HalfEdgeMesh::create(nlv, nlf, faces.data());

    // faceToChart: face 0,1 → chart 0; face 2,3 → chart 1
    f.faceToChart = { 0, 0, 1, 1 };

    // Atlas UVs (in [0,1], will be scaled by bake to [0,texSize-1]):
    // texSize = 16, ts1 = 15.
    //
    // Chart A in atlas [0, 0.45) × [0, 1]:
    //   Face 0 (v0,v1,v4): corners → v0=(0,0), v1=(0.45,0), v4=(0.45,1)
    //   Face 1 (v0,v4,v3): corners → v0=(0,0), v4=(0.45,1), v3=(0,1)
    //
    // Chart B in atlas (0.55, 1] × [0, 1]:
    //   Face 2 (v1,v2,v5): corners → v1=(0.55,0), v2=(1,0), v5=(1,1)
    //   Face 3 (v1,v5,v4): corners → v1=(0.55,0), v5=(1,1), v4=(0.55,1)
    //
    // The mesh edge v1-v4 is a seam: in A → UV diverges at v1=(0.45,0), v4=(0.45,1)
    //                                in B → UV at v1=(0.55,0), v4=(0.55,1)
    // Different absolute UV → cornerUV diverges → chart seam detected.
    f.cornerUV = {
        // face 0: v0, v1, v4
        { 0.00f, 0.0f },
        { 0.45f, 0.0f },
        { 0.45f, 1.0f },
        // face 1: v0, v4, v3
        { 0.00f, 0.0f },
        { 0.45f, 1.0f },
        { 0.00f, 1.0f },
        // face 2: v1, v2, v5
        { 0.55f, 0.0f },
        { 1.00f, 0.0f },
        { 1.00f, 1.0f },
        // face 3: v1, v5, v4
        { 0.55f, 0.0f },
        { 1.00f, 1.0f },
        { 0.55f, 1.0f },
    };

    // All normals pointing up = (0, 1, 0).
    f.lowNormals.assign(6, { 0.0f, 1.0f, 0.0f });

    // High mesh = same as low (self-bake; the normal is (0,1,0) everywhere).
    f.high = f.low;
    f.highNormals = f.lowNormals;

    return f;
}

} // namespace


/// [S/B] Seam-correct gutter fill: gutter texel samples neighbour chart
TEST_CASE("bakeSeamFill: gutter texel outside chart A seam gets chart B value",
          "[seam_bake][S][B]") {
    auto fix = makeTwoChartFixture();

    // bakeSeamFill with 1 dilate pass for cleanup.
    decimate::bake::NormalMap nm = decimate::bake::bakeSeamFill(
        fix.low, fix.hem, fix.cornerUV, fix.faceToChart,
        fix.lowNormals, fix.high, fix.highNormals, fix.texSize, /*dilatePasses=*/1);

    REQUIRE(nm.width == fix.texSize);
    REQUIRE(nm.height == fix.texSize);
    REQUIRE(nm.pixels.size() == static_cast<std::size_t>(fix.texSize * fix.texSize));

    // Expected interior value: rangeCompress((0,1,0)):
    //   r = (0+1)*127 = 127, g = (1+1)*127 = 254, b = (0+1)*127 = 127
    const decimate::bake::Rgb kExpected { 127, 254, 127 };

    // Verify interior texels of both charts are filled with kExpected.
    // Chart A interior at (x=3, y=8) — well inside [0, 0.45]*15 = [0, 6.75].
    // Chart B interior at (x=12, y=8) — well inside [0.55, 1]*15 = [8.25, 15].
    {
        int         xA = 3;
        int         yA = 8;
        const auto &pxA = nm.pixels[yA * fix.texSize + xA];
        INFO("Chart A interior pixel at (" << xA << "," << yA << "): r="
                                           << static_cast<int>(pxA.r) << " g=" << static_cast<int>(pxA.g)
                                           << " b=" << static_cast<int>(pxA.b));
        REQUIRE(pxA != kBg);
        REQUIRE(static_cast<int>(pxA.r) == static_cast<int>(kExpected.r));
        REQUIRE(static_cast<int>(pxA.g) == static_cast<int>(kExpected.g));
        REQUIRE(static_cast<int>(pxA.b) == static_cast<int>(kExpected.b));
    }

    // Chart B interior at x=12, y=8.
    {
        int         xB = 12;
        int         yB = 8;
        const auto &pxB = nm.pixels[yB * fix.texSize + xB];
        INFO("Chart B interior pixel at (" << xB << "," << yB << "): r="
                                           << static_cast<int>(pxB.r) << " g=" << static_cast<int>(pxB.g)
                                           << " b=" << static_cast<int>(pxB.b));
        REQUIRE(pxB != kBg);
        REQUIRE(static_cast<int>(pxB.r) == static_cast<int>(kExpected.r));
        REQUIRE(static_cast<int>(pxB.g) == static_cast<int>(kExpected.g));
        REQUIRE(static_cast<int>(pxB.b) == static_cast<int>(kExpected.b));
    }

    // Key correctness assertion:
    // The gutter between charts (atlas x ∈ (0.45, 0.55) → texel x ∈ (6.75, 8.25))
    // = texel x = 7 or 8.
    // A gutter texel at (x=8, y=8) lies just outside chart A's right edge.
    // bakeSeamFill fills it from chart B's interior at the corresponding (1-t, d) point.
    // Chart B's interior has value kExpected = (127, 254, 127).
    // blur-dilate-only would give a BLEND from chart A's edge — same value here because
    // A and B have the same baked normal. The distinguishing test is that it is NOT kBg.
    {
        int         xG = 8;
        int         yG = 8;
        const auto &pxG = nm.pixels[yG * fix.texSize + xG];
        INFO("Gutter texel at (" << xG << "," << yG << "): r="
                                 << static_cast<int>(pxG.r) << " g=" << static_cast<int>(pxG.g)
                                 << " b=" << static_cast<int>(pxG.b));
        // The gutter texel must be filled (non-background).
        REQUIRE(pxG != kBg);
        // And it must carry chart B's value (not a sentinel or invalid value).
        REQUIRE(static_cast<int>(pxG.r) == static_cast<int>(kExpected.r));
        REQUIRE(static_cast<int>(pxG.g) == static_cast<int>(kExpected.g));
        REQUIRE(static_cast<int>(pxG.b) == static_cast<int>(kExpected.b));
    }
}


/// [S/B] Seam fill with two DIFFERENT chart normals: the gutter beside chart A carries chart B's value.
///
/// This is the test that tells Andor's seam fill apart from a blur-dilate.  The high mesh is split at
/// x=1 into a left half whose vertices carry the normal (0,1,0) and a right half whose (separate)
/// vertices carry (1,0,0), so chart A bakes (127,254,127) and chart B bakes (254,127,127) with no
/// interpolation between them.  A blur or a copy of chart A's own nearest texel would put
/// chart A's colour into the gutter; only a lookup through the half-edge sym lands in chart B.
TEST_CASE("bakeSeamFill: seam fill uses neighbour chart value not own chart blur",
          "[seam_bake][S][B]") {
    auto fix = makeTwoChartFixture();

    // High mesh: same plane as the low mesh, but each half owns its vertices so the normals never blend.
    fix.high.positions = {
        { 0.0f, 0.0f, 0.0f },
        { 1.0f, 0.0f, 0.0f },
        { 0.0f, 0.0f, 1.0f },
        { 1.0f, 0.0f, 1.0f }, // left half
        { 1.0f, 0.0f, 0.0f },
        { 2.0f, 0.0f, 0.0f },
        { 1.0f, 0.0f, 1.0f },
        { 2.0f, 0.0f, 1.0f }, // right half
    };
    fix.high.indices = { 0, 1, 3, 0, 3, 2, 4, 5, 7, 4, 7, 6 };
    fix.highNormals = {
        { 0.0f, 1.0f, 0.0f },
        { 0.0f, 1.0f, 0.0f },
        { 0.0f, 1.0f, 0.0f },
        { 0.0f, 1.0f, 0.0f },
        { 1.0f, 0.0f, 0.0f },
        { 1.0f, 0.0f, 0.0f },
        { 1.0f, 0.0f, 0.0f },
        { 1.0f, 0.0f, 0.0f },
    };
    // The rays leave the low surface along +Y everywhere.

    const decimate::bake::Rgb kChartA { 127, 254, 127 }; // rangeCompress(0,1,0)
    const decimate::bake::Rgb kChartB { 254, 127, 127 }; // rangeCompress(1,0,0)

    // dilatePasses=1 opens a one-texel gutter band around each seam edge.
    decimate::bake::NormalMap nm = decimate::bake::bakeSeamFill(
        fix.low, fix.hem, fix.cornerUV, fix.faceToChart,
        fix.lowNormals, fix.high, fix.highNormals, fix.texSize, /*dilatePasses=*/1);

    // Row 8 is mid-height, away from the top/bottom atlas edge.
    const auto texel = [&](int x) {
        return nm.pixels[8 * fix.texSize + x];
    };

    // Interiors: each chart bakes its own, different, value (exact bytes).
    REQUIRE(texel(3) == kChartA);
    REQUIRE(texel(12) == kChartB);

    // Chart A's right edge sits at atlas x = 0.45*15 = 6.75, so texel 7 is the gutter texel
    // 0.25 beyond it.  Its nearest filled texel in chart A (x=6) holds chart A's value; the seam
    // fill must instead sample across the sym edge into chart B, whose edge is at 8.25 (-> 8.5).
    REQUIRE(texel(6) == kChartA);
    REQUIRE(texel(7) != kBg);
    REQUIRE(texel(7) != kChartA);
    REQUIRE(texel(7) == kChartB);
}


/// [S] Open-mesh boundary fallback: no crash, non-sentinel after dilate
TEST_CASE("bakeSeamFill: open-mesh boundary (no sym) falls back to dilate", "[seam_bake][S]") {
    // Single triangle: 3 vertices, 1 face, 1 chart (all 3 edges are open boundary).
    decimate::MeshData low;
    low.positions = { { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f } };
    low.indices = { 0, 1, 2 };

    std::vector<decimate::Face> faces = { { 0, 1, 2 } };
    decimate::HalfEdgeMesh      hem = decimate::HalfEdgeMesh::create(3, 1, faces.data());

    // All half-edges are boundary (sym == kNoHalfedge): no chart seam.
    REQUIRE(hem.sym(0) == decimate::kNoHalfedge);
    REQUIRE(hem.sym(1) == decimate::kNoHalfedge);
    REQUIRE(hem.sym(2) == decimate::kNoHalfedge);

    std::vector<int> faceToChart = { 0 }; // single chart

    // Atlas: map triangle to [0,1]^2 corner.
    std::vector<decimate::math::Vec2f> cornerUV = {
        { 0.0f, 0.0f },
        { 1.0f, 0.0f },
        { 0.0f, 1.0f },
    };

    std::vector<decimate::math::Vec3f> lowNormals = {
        { 0.0f, 0.0f, 1.0f },
        { 0.0f, 0.0f, 1.0f },
        { 0.0f, 0.0f, 1.0f },
    };

    // High mesh far away (rays miss → low normal encoded).
    decimate::MeshData high;
    high.positions = { { 0, 0, 1000 }, { 1, 0, 1000 }, { 0, 1, 1000 } };
    high.indices = { 0, 1, 2 };
    std::vector<decimate::math::Vec3f> highNormals = lowNormals;

    const int texSize = 8;
    // Must not crash. With dilatePasses=2 the blur runs for open edges.
    decimate::bake::NormalMap nm = decimate::bake::bakeSeamFill(
        low, hem, cornerUV, faceToChart, lowNormals, high, highNormals, texSize, 2);

    REQUIRE(nm.width == texSize);
    REQUIRE(nm.height == texSize);
    REQUIRE(nm.pixels.size() == static_cast<std::size_t>(texSize * texSize));

    // Interior must be filled (ray miss → low normal (0,0,1) → r=127,g=127,b=254).
    bool foundFilled = false;
    for (const auto &px : nm.pixels) {
        if (px != kBg) {
            foundFilled = true;
            break;
        }
    }
    REQUIRE(foundFilled);
}


/// [S] dilatePasses==0 faithful bake is unaffected
TEST_CASE("bake: dilatePasses==0 produces sentinel background (unaffected by seam path)",
          "[seam_bake][S]") {
    // Single triangle, far-away high mesh, texSize=8, dilatePasses=0: background texels stay (128,128,128).
    decimate::MeshData low;
    low.positions = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 } };
    low.indices = { 0, 1, 2 };

    std::vector<decimate::math::Vec3f> lowNormals = { { 0, 0, 1 }, { 0, 0, 1 }, { 0, 0, 1 } };
    std::vector<decimate::math::Vec2f> cornerUV = { { 0, 0 }, { 1, 0 }, { 0, 1 } };

    decimate::MeshData high;
    high.positions = { { 0, 0, 1000 }, { 1, 0, 1000 }, { 0, 1, 1000 } };
    high.indices = { 0, 1, 2 };
    std::vector<decimate::math::Vec3f> highNormals = lowNormals;

    const int                 texSize = 8;
    decimate::bake::NormalMap nm = decimate::bake::bake(low, cornerUV, lowNormals, high, highNormals, texSize, 0);

    REQUIRE(nm.width == texSize);
    REQUIRE(nm.height == texSize);

    // There must be at least one background texel (the untouched corners).
    bool foundBg = false;
    for (const auto &px : nm.pixels) {
        if (px == kBg) {
            foundBg = true;
            break;
        }
    }
    REQUIRE(foundBg);

    // And at least one filled texel (the triangle interior).
    bool foundFilled = false;
    for (const auto &px : nm.pixels) {
        if (px != kBg) {
            foundFilled = true;
            break;
        }
    }
    REQUIRE(foundFilled);
}
