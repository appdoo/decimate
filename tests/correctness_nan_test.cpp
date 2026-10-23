// This file is part of the decimate library; see decimate.h for version/license details.
/// @file correctness_nan_test.cpp
/// @brief Spec 24-02: NaN / non-finite / out-of-range inputs on the reduce, seam, bake and chart paths.
///
/// N1 [S]  reduce() with a NaN position terminates (the NaN cost used to re-queue forever)
/// N2 [S]  reduce()/reduceHalfEdge() throw std::invalid_argument on an out-of-range face index
/// N3 [B]  a zero-UV-area face does not overwrite a neighbouring chart's texels in bake()
/// N4 [S]  seams of total 3D length 0 are "no seam": erase() is a passthrough, seamEnergy() is 0
/// N5 [S]  non-finite areas / UVs are degenerate in packAtlas() and selfOverlapPairs()
///
/// Every case that could hang runs on a detached worker with a deadline, so a regression fails the
/// test instead of hanging ctest.

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normal_map.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/reduce/halfedge_reduce.h>
#include <decimate/reduce/reducer.h>
#include <decimate/seam/erase.h>
#include <decimate/seam/seam_energy.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/textured_mesh.h>
#include <future>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using decimate::MeshData;
using decimate::math::Vec2f;
using decimate::math::Vec3f;

constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();


/// An n x n vertex grid in the XY plane (gently curved so quadrics are not all coplanar).
MeshData makeGrid(int n) {
    MeshData m;
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            const float fx = static_cast<float>(x);
            const float fy = static_cast<float>(y);
            m.positions.push_back({ fx, fy, 0.05f * fx * fy });
        }
    }
    for (int y = 0; y + 1 < n; ++y) {
        for (int x = 0; x + 1 < n; ++x) {
            const int a = y * n + x;
            const int b = a + 1;
            const int c = a + n;
            const int d = c + 1;
            m.indices.insert(m.indices.end(), { a, b, d, a, d, c });
        }
    }
    return m;
}


/// Runs @p fn on a detached worker; nullopt when it has not finished within @p limit (a hang).
template <class Fn>
std::optional<MeshData> runBounded(Fn fn, std::chrono::seconds limit) {
    auto promise = std::make_shared<std::promise<MeshData>>();
    auto future = promise->get_future();
    std::thread([promise, fn = std::move(fn)] { promise->set_value(fn()); }).detach();
    if (future.wait_for(limit) != std::future_status::ready) {
        return std::nullopt;
    }
    return future.get();
}


/// True when @p what contains @p needle.
bool contains(const char *what, const std::string &needle) {
    return std::string(what).find(needle) != std::string::npos;
}

} // namespace


/// N1 — one NaN position must not hang the lazy-heap loop; the NaN vertex is never a collapse survivor.
TEST_CASE("N1: reduce() with a NaN position terminates and keeps the NaN vertex isolated",
          "[nan][N1]") {
    MeshData in = makeGrid(7);
    in.positions[24] = { kNaN, kNaN, kNaN }; // centre vertex

    const MeshData copy = in;
    const auto     out = runBounded([copy] { return decimate::reduce(copy, 10); }, std::chrono::seconds(5));
    REQUIRE(out.has_value()); // nullopt = the old infinite re-queue loop

    CHECK(out->positions.size() < in.positions.size());
    CHECK(!out->indices.empty());

    // The only non-finite vertex is the original one: no collapse placed a survivor at a NaN point.
    int nonFinite = 0;
    for (const Vec3f &p : out->positions) {
        nonFinite += std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z) ? 0 : 1;
    }
    CHECK(nonFinite == 1);
    for (const int idx : out->indices) {
        CHECK(idx >= 0);
        CHECK(static_cast<std::size_t>(idx) < out->positions.size());
    }
}


/// N2 — an index outside [0, positions.size()) is rejected at the entry, before any table is indexed.
TEST_CASE("N2: reduce()/reduceHalfEdge() throw invalid_argument on out-of-range indices",
          "[nan][N2]") {
    const MeshData good = makeGrid(4);
    const int      outOfRange[] = { -1, static_cast<int>(good.positions.size()) };

    for (const int bad : outOfRange) {
        MeshData in = good;
        in.indices[4] = bad;
        CAPTURE(bad);

        try {
            (void)decimate::reduce(in, 3);
            FAIL("reduce() accepted an out-of-range index");
        } catch (const std::invalid_argument &e) {
            CHECK(contains(e.what(), "reduce:"));
        }

        try {
            (void)decimate::reduceHalfEdge(in, 3);
            FAIL("reduceHalfEdge() accepted an out-of-range index");
        } catch (const std::invalid_argument &e) {
            CHECK(contains(e.what(), "reduceHalfEdge:"));
        }
    }
}


namespace {

/// A valid quad chart (faces 0,1) covering UV [0.25,0.75]^2, plus optionally a zero-UV-area face
/// whose collinear UVs run along the quad's diagonal (every texel centre on it computes 0/0).
struct BakeFixture {
    MeshData             low;
    std::vector<Vec2f>   cornerUV;
    std::vector<Vec3f>   normals;
    static constexpr int kTex = 16;

    explicit BakeFixture(bool withDegenerateFace) {
        low.positions = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 0, 1 }, { 0, 0, 1 } };
        low.indices = { 0, 1, 2, 0, 2, 3 };
        cornerUV = { { 0.25f, 0.25f }, { 0.75f, 0.25f }, { 0.75f, 0.75f }, { 0.25f, 0.25f }, { 0.75f, 0.75f }, { 0.25f, 0.75f } };
        if (withDegenerateFace) {
            low.indices.insert(low.indices.end(), { 0, 1, 2 });
            cornerUV.insert(cornerUV.end(), { { 0.25f, 0.25f }, { 0.75f, 0.75f }, { 0.5f, 0.5f } });
        }
        normals.assign(low.positions.size(), { 0.0f, 1.0f, 0.0f });
    }
};

} // namespace


/// N3 — the valid chart's texels are identical with and without the degenerate face.
TEST_CASE("N3: a zero-UV-area face does not overwrite a valid chart's texels in bake()",
          "[nan][N3]") {
    const BakeFixture clean(false);
    const BakeFixture dirty(true);

    const auto expected = decimate::bake::bake(clean.low, clean.cornerUV, clean.normals, clean.low,
                                               clean.normals, BakeFixture::kTex, 0, {}, 1);
    const auto actual = decimate::bake::bake(dirty.low, dirty.cornerUV, dirty.normals, dirty.low,
                                             dirty.normals, BakeFixture::kTex, 0, {}, 1);

    REQUIRE(actual.pixels.size() == expected.pixels.size());
    int differing = 0;
    for (std::size_t i = 0; i < expected.pixels.size(); ++i) {
        differing += actual.pixels[i] == expected.pixels[i] ? 0 : 1;
    }
    CHECK(differing == 0);

    // The diagonal texel the degenerate face used to claim still carries the quad's up-normal.
    const auto &px = actual.pixels[5 * BakeFixture::kTex + 5];
    CHECK(static_cast<int>(px.g) == 254);
}


namespace {

/// Two triangles sharing the 3D edge (0,1) whose endpoints coincide in space, with different UVs on
/// each side: findSeams reports one seam of 3D length 0.
decimate::TexturedMesh zeroLengthSeamMesh() {
    decimate::TexturedMesh tm;
    tm.positions = { { 0, 0, 0 }, { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 } };
    tm.faceV = { 0, 1, 2, 1, 0, 3 };
    tm.uvs = { { 0.1f, 0.1f }, { 0.4f, 0.1f }, { 0.2f, 0.4f }, // face 0
               { 0.6f, 0.6f },
               { 0.9f, 0.6f },
               { 0.7f, 0.9f } }; // face 1
    tm.faceU = { 0, 1, 2, 3, 4, 5 };
    return tm;
}


decimate::bake::NormalMap gradientMap(int size) {
    decimate::bake::NormalMap map;
    map.width = size;
    map.height = size;
    for (int i = 0; i < size * size; ++i) {
        const auto v = static_cast<std::uint8_t>(i * 3 % 256);
        map.pixels.push_back({ v, static_cast<std::uint8_t>(255 - v), 128 });
    }
    return map;
}

} // namespace


/// N4 — total seam length 0 means no seam, not a NaN system.
TEST_CASE("N4: seams of total 3D length 0 are treated as no seam", "[nan][N4]") {
    const decimate::TexturedMesh tm = zeroLengthSeamMesh();
    const auto                   map = gradientMap(8);

    const auto energy = decimate::seam::seamEnergy(tm, map);
    CHECK(energy[0] == 0.0);
    CHECK(energy[1] == 0.0);
    CHECK(energy[2] == 0.0);

    decimate::seam::EraseResult r;
    REQUIRE_NOTHROW(r = decimate::seam::erase(tm, map, decimate::seam::NormalSpace::Object));
    CHECK(r.map.pixels == map.pixels);
    for (std::size_t c = 0; c < 3; ++c) {
        CHECK(r.dBefore[c] == 0.0);
        CHECK(r.dAfter[c] == 0.0);
    }
    CHECK(r.unknownTexels == 0);
}


namespace {

/// A chart over mesh triangles @p tri with local vertex UVs @p uvs; faces are numbered from @p firstFace.
decimate::chart::HalfEdgeChart handChart(const std::vector<int> &tri, const std::vector<Vec2f> &uvs,
                                         int firstFace = 0) {
    decimate::chart::HalfEdgeChart ac;
    const int                      nf = static_cast<int>(tri.size() / 3);
    ac.mesh = decimate::HalfEdgeMesh::create(static_cast<int>(uvs.size()), nf,
                                             reinterpret_cast<const decimate::Face *>(tri.data()));
    ac.uv = uvs;
    for (int f = 0; f < nf; ++f) {
        ac.faces.push_back(firstFace + f);
    }
    return ac;
}

} // namespace


/// N5a — a NaN 3D area of one chart must not turn every chart's scale into NaN.
TEST_CASE("N5: packAtlas() with a NaN-area chart keeps the other chart finite", "[nan][N5]") {
    // Mesh: quad A (faces 0,1) and quad B (faces 2,3) with one NaN position in B.
    MeshData data;
    data.positions = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 2, 0, 0 }, { 3, 0, 0 }, { 3, 1, 0 }, { kNaN, kNaN, kNaN } };
    data.indices = { 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7 };

    const std::vector<int>                      quad = { 0, 1, 2, 0, 2, 3 };
    const std::vector<Vec2f>                    square = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
    std::vector<decimate::chart::HalfEdgeChart> charts;
    charts.push_back(handChart(quad, square, 0));
    charts.push_back(handChart(quad, square, 2));
    charts[0].localToMesh = { 0, 1, 2, 3 };
    charts[1].localToMesh = { 4, 5, 6, 7 };

    decimate::chart::packAtlas(charts, data, 64);

    for (const Vec2f &uv : charts[0].uv) {
        CHECK(std::isfinite(uv.x));
        CHECK(std::isfinite(uv.y));
    }
}


/// N5b — a NaN-UV triangle is degenerate: it neither hides real overlaps nor invents new ones.
TEST_CASE("N5: selfOverlapPairs() treats NaN-UV triangles as degenerate", "[nan][N5]") {
    // Triangles 0,1 form a clean quad; triangle 2 has a NaN corner.
    const std::vector<int> tri = { 0, 1, 2, 1, 3, 2, 4, 5, 6 };

    auto uvs = std::vector<Vec2f> { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 1, 1 }, { 0.2f, 0.2f }, { 0.4f, 0.2f }, { kNaN, 0.4f } };
    CHECK(decimate::chart::selfOverlapPairs(handChart(tri, uvs)) == 0);

    // Triangle 2 all-NaN: still nothing to report, and no out-of-range grid cell.
    uvs[4] = { kNaN, kNaN };
    uvs[5] = { kNaN, kNaN };
    CHECK(decimate::chart::selfOverlapPairs(handChart(tri, uvs)) == 0);

    // A real overlap next to a NaN triangle is still found.
    const std::vector<int>   folded = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    const std::vector<Vec2f> fuv = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 0.1f, 0.1f }, { 0.6f, 0.1f }, { 0.1f, 0.6f }, { kNaN, 0 }, { kNaN, 1 }, { kNaN, 2 } };
    CHECK(decimate::chart::selfOverlapPairs(handChart(folded, fuv)) == 1);
}
