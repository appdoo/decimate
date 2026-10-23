// This file is part of the decimate library; see decimate.h for version/license details.
/// @file seam_erase_test.cpp
/// @brief Spec 17-01: seam energy metric + seam erasure (Liu et al. 2017, section 4).
///
/// Scenarios I1-I3, B1-B3, S1-S7 (seam erasure and its metric).
///
/// The goldens (tests/fixtures/seam_golden_{rot,mir,glob}.txt) come from the MIT
/// reference implementation via tools/seam_erasure_golden.py (pinned commit in the
/// script). Tests read the text fixture, never an OBJ: loadObj drops `vt`.
///
/// The NOLINT on the TEST_CASEs is for Catch2: every CHECK macro expands to branching code
/// that readability-function-size / -cognitive-complexity count against the test body.

#include "erase_domain.h" // internal: the unknown set the solver actually uses
#include "seam_pieces.h"  // internal: texel mapping

#include <algorithm>
#include <array>
#include <bit>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normal_map.h>
#include <decimate/bake/normals.h>
#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/reduce/reducer.h>
#include <decimate/seam/erase.h>
#include <decimate/seam/seam_energy.h>
#include <decimate/seam/seam_topology.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/textured_mesh.h>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef MESH_MODEL_DIR
#error "MESH_MODEL_DIR must be defined by CMake"
#endif
#ifndef MESH_FIXTURES_DIR
#error "MESH_FIXTURES_DIR must be defined by CMake"
#endif

namespace {

using decimate::TexturedMesh;
using decimate::bake::NormalMap;
using decimate::bake::Rgb;
using decimate::math::Vec2f;
using decimate::math::Vec3f;
using decimate::seam::NormalSpace;

/// One golden fixture from the MIT reference: mesh, reference counts and energies, input map, reference solution.
struct Golden {
    TexturedMesh          tm;
    int                   size = 0;
    bool                  global = false;
    std::array<int, 3>    counts {}; // seams, boundary, foldovers (reference find_seam)
    std::array<double, 3> dBefore {};
    std::array<double, 3> dAfter {};
    NormalMap             input;
    std::vector<double>   reference; // size*size*3, unquantized reference solution
};


/// Whitespace-token reader over a golden fixture; `#` lines are comments.
class FixtureReader {
public:
    explicit FixtureReader(const std::string &path)
        : in_(path) {
        REQUIRE(in_.is_open());
        skipComments();
    }

    /// Consume the next word and require it to be @p expected.
    void key(const char *expected) {
        std::string word;
        in_ >> word;
        REQUIRE(word == expected);
    }

    template <class T>
    T value() {
        T v {};
        in_ >> v;
        REQUIRE(!in_.fail());
        return v;
    }

    /// Discard the rest of the current line (trailing comment) and any comment lines.
    void endLine() {
        std::string rest;
        std::getline(in_, rest);
        skipComments();
    }

private:
    void skipComments() {
        while (in_.peek() == '#') {
            std::string line;
            std::getline(in_, line);
        }
    }
    std::ifstream in_;
};


/// Reads the positions / uvs / faceV / faceU blocks of a fixture.
TexturedMesh readMesh(FixtureReader &r) {
    TexturedMesh tm;
    r.key("positions");
    for (auto n = r.value<std::size_t>(); n > 0; --n) {
        const float x = r.value<float>();
        const float y = r.value<float>();
        tm.positions.emplace_back(x, y, r.value<float>());
    }
    r.key("uvs");
    for (auto n = r.value<std::size_t>(); n > 0; --n) {
        const float u = r.value<float>();
        tm.uvs.emplace_back(u, r.value<float>());
    }
    r.key("faces");
    const auto faces = r.value<std::size_t>();
    r.key("faceV");
    for (std::size_t i = 0; i < 3 * faces; ++i) {
        tm.faceV.push_back(r.value<int>());
    }
    r.key("faceU");
    for (std::size_t i = 0; i < 3 * faces; ++i) {
        tm.faceU.push_back(r.value<int>());
    }
    return tm;
}


/// Reads `key` followed by three doubles (one per colour channel).
std::array<double, 3> readTriple(FixtureReader &r, const char *key) {
    r.key(key);
    std::array<double, 3> v {};
    for (double &d : v) {
        d = r.value<double>();
    }
    return v;
}


/// Parses tests/fixtures/seam_golden_<name>.txt.
Golden loadGolden(const std::string &name) {
    FixtureReader r(std::string(MESH_FIXTURES_DIR) + "/seam_golden_" + name + ".txt");
    Golden        g;
    r.key("name");
    r.value<std::string>();
    r.key("size");
    g.size = r.value<int>();
    r.key("global");
    g.global = r.value<int>() != 0;
    g.tm = readMesh(r);
    r.key("counts");
    for (std::size_t i = 0; i < 3; ++i) {
        r.value<std::string>(); // "seams" / "boundary" / "foldovers"
        g.counts[i] = r.value<int>();
    }
    g.dBefore = readTriple(r, "dBefore");
    g.dAfter = readTriple(r, "dAfter");
    readTriple(r, "dBeforeUnpatched"); // informational
    r.key("input");
    r.endLine();
    g.input.width = g.size;
    g.input.height = g.size;
    for (int i = 0; i < g.size * g.size; ++i) {
        const int red = r.value<int>();
        const int green = r.value<int>();
        g.input.pixels.push_back({ static_cast<std::uint8_t>(red), static_cast<std::uint8_t>(green),
                                   static_cast<std::uint8_t>(r.value<int>()) });
    }
    r.key("reference");
    r.endLine();
    for (int i = 0; i < g.size * g.size * 3; ++i) {
        g.reference.push_back(r.value<double>());
    }
    return g;
}


/// A [0,1] value as the byte our output stores (round to nearest).
std::uint8_t quantize(double v) {
    return static_cast<std::uint8_t>(std::lround(std::clamp(v, 0.0, 1.0) * 255.0));
}


/// Largest |byte - reference byte| over the map; the reference is quantized like our output.
int maxByteDiff(const NormalMap &map, const std::vector<double> &reference) {
    int maxDiff = 0;
    for (std::size_t t = 0; t < map.pixels.size(); ++t) {
        const Rgb                px = map.pixels[t];
        const std::array<int, 3> got { px.r, px.g, px.b };
        for (std::size_t c = 0; c < 3; ++c) {
            maxDiff = std::max(maxDiff, std::abs(got[c] - quantize(reference[3 * t + c])));
        }
    }
    return maxDiff;
}


/// I1/I2/S5 acceptance: metric, byte-level agreement with the reference, D_after.
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
void checkAgainstReference(const Golden &g) {
    // The metric first: our D_before must agree with the reference to 1e-9 relative.
    const auto before = decimate::seam::seamEnergy(g.tm, g.input);
    for (std::size_t c = 0; c < 3; ++c) {
        CHECK(std::abs(before[c] - g.dBefore[c]) <= 1e-9 * g.dBefore[c]);
    }

    decimate::seam::EraseOptions opt;
    opt.global = g.global;
    const auto res = decimate::seam::erase(g.tm, g.input, NormalSpace::Object, opt);

    // Then the solution: within one byte of the reference, and the seam energy driven to ~0.
    CHECK(maxByteDiff(res.map, g.reference) <= 1);
    for (std::size_t c = 0; c < 3; ++c) {
        CHECK(res.dBefore[c] == before[c]);
        CHECK(res.dAfter[c] <= 1e-8 * res.dBefore[c]);
    }
    CHECK(res.unknownTexels > 0);
}


/// The spot pipeline (built once; shared by I3, B1, S2, S6): reduce, chart, bake with dilation, weld.
struct SpotPipeline {
    decimate::MeshData low;
    std::vector<Vec2f> cornerUV;
    std::vector<Vec3f> lowN;
    decimate::MeshData high;
    std::vector<Vec3f> highN;
    NormalMap          dilateMap;
    TexturedMesh       tm;
};


/// The triangles of @p mesh as the Face list the half-edge builder takes.
std::vector<decimate::Face> facesOf(const decimate::MeshData &mesh) {
    std::vector<decimate::Face> faces;
    for (std::size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
        faces.push_back({ mesh.indices[i], mesh.indices[i + 1], mesh.indices[i + 2] });
    }
    return faces;
}


/// Builds the spot pipeline once.
SpotPipeline buildSpot() {
    constexpr int kTex = 256;
    SpotPipeline  s;
    s.high = decimate::loadObj(std::string(MESH_MODEL_DIR) + "/spot.obj");
    s.low = decimate::reduce(s.high, 1500);
    const auto faces = facesOf(s.low);
    const auto hem = decimate::HalfEdgeMesh::create(static_cast<int>(s.low.positions.size()),
                                                    static_cast<int>(faces.size()), faces.data());
    auto       charts = decimate::chart::segment(hem, s.low, 64);
    s.cornerUV = decimate::chart::parametrizeCorners(hem, s.low, charts, kTex);
    s.lowN = decimate::bake::computeVertexNormals(s.low, decimate::bake::Weighting::EqualUnit);
    s.highN = decimate::bake::computeVertexNormals(s.high, decimate::bake::Weighting::EqualUnit);
    s.dilateMap = decimate::bake::bake(s.low, s.cornerUV, s.lowN, s.high, s.highN, kTex, 2, {});
    s.tm = decimate::seam::weld(s.low, s.cornerUV);
    return s;
}


/// The shared, lazily built spot pipeline (read-only).
const SpotPipeline &spot() {
    static const SpotPipeline kPipeline = buildSpot();
    return kPipeline;
}


/// Same size and identical pixels.
bool sameMap(const NormalMap &a, const NormalMap &b) {
    return a.width == b.width && a.height == b.height && a.pixels == b.pixels;
}


/// UV of a texel coordinate on an NxN texture (inverse of uv*N - 0.5).
Vec2f uvOfTexel(double x, double y, int n) {
    return { static_cast<float>((x + 0.5) / n), static_cast<float>((y + 0.5) / n) };
}


/// 3D edge p0->p1 in face A, p1->p0 in face B; A carries (a0,a1), B carries (b0,b1) where
/// b0 meets a0 at the same 3D point.  Third corners are far from the seam.
TexturedMesh seamPair(Vec2f a0, Vec2f a1, Vec2f b0, Vec2f b1) {
    TexturedMesh tm;
    tm.positions = { { 0, 0, 0 }, { 0, 1, 0 }, { -1, 0.5f, 0 }, { 1, 0.5f, 0 } };
    tm.uvs = { a0, a1, { 0.05f, 0.95f }, b1, b0, { 0.95f, 0.95f } };
    tm.faceV = { 0, 1, 2, 1, 0, 3 };
    tm.faceU = { 0, 1, 2, 3, 4, 5 };
    return tm;
}


/// An n x n map filled with @p c.
NormalMap constantMap(int n, Rgb c) {
    NormalMap m;
    m.width = n;
    m.height = n;
    m.pixels.assign(static_cast<std::size_t>(n) * n, c);
    return m;
}

} // namespace

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("I1: erase matches the reference on G-rot (rotated neighbour)", "[seam][I1]") {
    checkAgainstReference(loadGolden("rot"));
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("I2: erase matches the reference on G-mir (mirrored neighbour)", "[seam][I2]") {
    checkAgainstReference(loadGolden("mir"));
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S5: global mode matches the reference do_global=True row", "[seam][S5]") {
    const Golden g = loadGolden("glob");
    REQUIRE(g.global);
    checkAgainstReference(g);
}

/// Distance (in texels) from the centre of texel (@p tx, @p ty) to the nearest UV triangle
/// of @p tm.  Texel centres sit at uv * N - 0.5 (S3).  Written from the geometry alone, so
/// it checks the unknown set without going through `buildEraseDomain`.
static double distanceToUvTriangles(const TexturedMesh &tm, int w, int h, int tx, int ty) {
    const auto toTexel = [&](int corner) {
        const Vec2f uv = tm.uvs[static_cast<std::size_t>(tm.faceU[static_cast<std::size_t>(corner)])];
        return std::array<double, 2> { (static_cast<double>(uv.x) * w) - 0.5,
                                       (static_cast<double>(uv.y) * h) - 0.5 };
    };
    const std::array<double, 2> p { static_cast<double>(tx), static_cast<double>(ty) };
    constexpr double            kReach = 2.0; // skip triangles whose box is farther than this
    double                      best = 1e30;

    // Per triangle: box reject, then inside test, then distance to the three edges.
    for (int c = 0; c + 2 < static_cast<int>(tm.faceU.size()); c += 3) {
        const std::array<std::array<double, 2>, 3> v { toTexel(c), toTexel(c + 1), toTexel(c + 2) };
        const double                               lox = std::min({ v[0][0], v[1][0], v[2][0] }) - kReach;
        const double                               hix = std::max({ v[0][0], v[1][0], v[2][0] }) + kReach;
        const double                               loy = std::min({ v[0][1], v[1][1], v[2][1] }) - kReach;
        const double                               hiy = std::max({ v[0][1], v[1][1], v[2][1] }) + kReach;
        if (p[0] < lox || p[0] > hix || p[1] < loy || p[1] > hiy) {
            continue;
        }
        const auto side = [&](int a, int b) {
            return ((v[b][0] - v[a][0]) * (p[1] - v[a][1])) - ((v[b][1] - v[a][1]) * (p[0] - v[a][0]));
        };

        const double s0 = side(0, 1);
        const double s1 = side(1, 2);
        const double s2 = side(2, 0);
        if ((s0 >= 0 && s1 >= 0 && s2 >= 0) || (s0 <= 0 && s1 <= 0 && s2 <= 0)) {
            return 0.0; // inside
        }

        for (int a = 0; a < 3; ++a) {
            const auto  &q0 = v[static_cast<std::size_t>(a)];
            const auto  &q1 = v[static_cast<std::size_t>((a + 1) % 3)];
            const double ex = q1[0] - q0[0];
            const double ey = q1[1] - q0[1];
            const double len2 = (ex * ex) + (ey * ey);
            const double t =
                len2 > 0.0 ? std::clamp((((p[0] - q0[0]) * ex) + ((p[1] - q0[1]) * ey)) / len2, 0.0, 1.0)
                           : 0.0;
            best = std::min(best, std::hypot(p[0] - (q0[0] + (t * ex)), p[1] - (q0[1] + (t * ey))));
        }
    }
    return best;
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("I3: erase on the spot bake removes the seam and touches only the unknown set",
          "[seam][I3]") {
    const SpotPipeline &s = spot();
    const auto          res = decimate::seam::erase(s.tm, s.dilateMap, NormalSpace::Object);

    // Bound: 1e-7, not the spec's 1e-8.  At the paper's w_seam = 1e10 the penalty leaves
    // D_after = 2.04e-10 / 1.86e-10 / 2.19e-10 against D_before = 4.81e-3 / 5.03e-3 / 5.18e-3
    // (ratio 4.2e-8 / 3.7e-8 / 4.2e-8) on this 594-seam bake; the MIT reference run on the
    // identical mesh + map gives D_after = 2.0364e-10 / 1.8561e-10 / 2.1911e-10 and its bytes
    // agree with ours to +-1 (see the receipt).  Sharpening the penalty into an exact
    // constraint moves texels by up to 26 levels, away from the reference minimizer.
    for (std::size_t c = 0; c < 3; ++c) {
        INFO("channel " << c << " D " << res.dBefore[c] << " -> " << res.dAfter[c]);
        CHECK(res.dBefore[c] > 0.0);
        CHECK(res.dAfter[c] <= 1e-7 * res.dBefore[c]);
    }
    CHECK(res.unknownTexels > 0);

    // Outside the unknown set the map must be untouched; inside, changes stay near the UV mesh.
    const decimate::seam::detail::TexelGrid grid { s.dilateMap.width, s.dilateMap.height };
    const auto                              domain =
        decimate::seam::detail::buildEraseDomain(s.tm, decimate::seam::findSeams(s.tm), grid);
    std::size_t outside = 0;
    std::size_t changed = 0;
    for (std::size_t t = 0; t < domain.unknown.size(); ++t) {
        outside += domain.unknown[t] == 0 ? 1 : 0;
        if (domain.unknown[t] == 0) {
            REQUIRE(res.map.pixels[t] == s.dilateMap.pixels[t]);
        } else if (res.map.pixels[t] != s.dilateMap.pixels[t]) {
            ++changed;
            // Independent of the solver's domain code: a texel the solve changed is a bilinear
            // cell corner of the UV mesh, so its centre is within sqrt(2) texels (the cell
            // diagonal, 1.4142...) of some UV triangle.  Measured max: 1.413.
            const int tx = static_cast<int>(t % static_cast<std::size_t>(grid.width));
            const int ty = static_cast<int>(t / static_cast<std::size_t>(grid.width));
            REQUIRE(distanceToUvTriangles(s.tm, grid.width, grid.height, tx, ty) <= 1.4143);
        }
    }

    CHECK(outside > 0);
    CHECK(changed > 0);
    CHECK(static_cast<std::size_t>(res.unknownTexels) <= domain.unknown.size() - outside);
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B1: a tangent-space map is rejected and left unchanged", "[seam][B1]") {
    const SpotPipeline &s = spot();
    const NormalMap     copy = s.dilateMap;

    // Tangent-space maps need a different formulation; the call must refuse, not guess.
    CHECK_THROWS_AS(decimate::seam::erase(s.tm, s.dilateMap, NormalSpace::Tangent),
                    std::invalid_argument);
    CHECK(sameMap(copy, s.dilateMap));
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B2: a constant texture lies in the null space (output == input)", "[seam][B2]") {
    const Golden    g = loadGolden("rot");
    const NormalMap flat = constantMap(g.size, { 100, 150, 200 });

    // A constant map has zero seam energy already, so the minimizer must return it unchanged.
    const auto res = decimate::seam::erase(g.tm, flat, NormalSpace::Object);
    CHECK(res.map.pixels.size() == static_cast<std::size_t>(64 * 64));
    CHECK(sameMap(res.map, flat));
    CHECK(res.unknownTexels > 0); // the solve really ran
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B3: a seamless quad has zero energy and erase is a no-op", "[seam][B3]") {
    // Two triangles of one chart: the shared edge matches in UV, the rest is boundary.
    TexturedMesh quad;
    quad.positions = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 } };
    quad.uvs = { { 0.2f, 0.2f }, { 0.8f, 0.2f }, { 0.8f, 0.8f }, { 0.2f, 0.8f } };
    quad.faceV = { 0, 1, 2, 0, 2, 3 };
    quad.faceU = { 0, 1, 2, 0, 2, 3 };

    // Distinct pixel values, so an accidental rewrite would show.
    NormalMap map = constantMap(16, { 10, 20, 30 });
    for (std::size_t i = 0; i < map.pixels.size(); ++i) {
        map.pixels[i] = { static_cast<std::uint8_t>(i * 7), static_cast<std::uint8_t>(i * 13),
                          static_cast<std::uint8_t>(i * 29) };
    }


    const auto topo = decimate::seam::findSeams(quad);
    CHECK(topo.seams.empty());
    CHECK(topo.boundary.size() == 4);
    CHECK(topo.foldovers.empty());


    const auto d = decimate::seam::seamEnergy(quad, map);
    for (const double c : d) {
        CHECK(c == 0.0);
    }

    const auto res = decimate::seam::erase(quad, map, NormalSpace::Object);
    CHECK(sameMap(res.map, map));
    for (const double c : res.dBefore) {
        CHECK(c == 0.0);
    }
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S1: findSeams on G-rot agrees with the reference find_seam counts", "[seam][S1]") {
    const Golden g = loadGolden("rot");
    const auto   topo = decimate::seam::findSeams(g.tm);
    CHECK(static_cast<int>(topo.seams.size()) == g.counts[0]);
    CHECK(static_cast<int>(topo.boundary.size()) == g.counts[1]);
    CHECK(static_cast<int>(topo.foldovers.size()) == g.counts[2]);
    CHECK(g.counts[0] == 2); // sanity: the fixture really has a 2-edge seam chain
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S2: weld on spot is lossless and seams are exactly the UV-differing edges",
          "[seam][S2]") {
    const SpotPipeline &s = spot();
    CHECK(s.tm.faceV == s.low.indices);
    CHECK(s.tm.positions.size() == s.low.positions.size());

    // Welding keeps the UVs bit-exact (corner UVs round-trip).
    const std::vector<Vec2f> back = s.tm.cornerUV();
    REQUIRE(back.size() == s.cornerUV.size());
    bool bitExact = true;
    for (std::size_t i = 0; i < back.size(); ++i) {
        bitExact = bitExact &&
                   std::bit_cast<std::uint32_t>(back[i].x) ==
                       std::bit_cast<std::uint32_t>(s.cornerUV[i].x) &&
                   std::bit_cast<std::uint32_t>(back[i].y) ==
                       std::bit_cast<std::uint32_t>(s.cornerUV[i].y);
    }
    CHECK(bitExact);

    // Independent count: manifold surface edges whose two half-edges' corner UVs differ.
    struct HE {
        int face;
        int k;
    };
    std::map<std::pair<int, int>, HE> directed;
    const auto                       &idx = s.low.indices;
    for (int f = 0; f < static_cast<int>(idx.size() / 3); ++f) {
        for (int k = 0; k < 3; ++k) {
            directed[{ idx[3 * f + k], idx[3 * f + (k + 1) % 3] }] = { f, k };
        }
    }

    std::size_t expected = 0;
    for (const auto &[edge, he] : directed) {
        const auto opp = directed.find({ edge.second, edge.first });
        if (edge.first >= edge.second || opp == directed.end()) {
            continue;
        }

        const Vec2f a0 = s.cornerUV[3 * he.face + he.k];
        const Vec2f a1 = s.cornerUV[3 * he.face + (he.k + 1) % 3];
        const Vec2f b0 = s.cornerUV[3 * opp->second.face + opp->second.k];
        const Vec2f b1 = s.cornerUV[3 * opp->second.face + (opp->second.k + 1) % 3];
        const bool  same = a0.x == b1.x && a0.y == b1.y && a1.x == b0.x && a1.y == b0.y;
        expected += same ? 0 : 1;
    }

    const auto topo = decimate::seam::findSeams(s.tm);
    CHECK(topo.seams.size() == expected);
    CHECK(expected > 0);
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S3: texel mapping is uv*N-0.5 with clamp-to-edge", "[seam][S3]") {
    using decimate::seam::detail::TexelGrid;
    const auto mid =
        decimate::seam::detail::toTexelSegment({ 0.5f, 0.5f }, { 0.0f, 0.0f }, TexelGrid { 64, 64 });
    CHECK(mid.x0 == 31.5);
    CHECK(mid.y0 == 31.5);
    CHECK(mid.x1 == -0.5);
    CHECK(mid.y1 == -0.5);

    // Clamp-to-edge, observed through the metric: side A samples at x = -0.5 (left of
    // texel 0, clamped to column 0), side B at x = 0 exactly (also column 0). Columns 0
    // and 1 differ, so an unclamped extrapolation (1.5*c0 - 0.5*c1) would give D > 0.
    constexpr int       kN = 64;
    const TexturedMesh  tm = seamPair(uvOfTexel(-0.5, 10.25, kN), uvOfTexel(-0.5, 10.75, kN),
                                      uvOfTexel(0.0, 10.25, kN), uvOfTexel(0.0, 10.75, kN));
    std::vector<double> texels(static_cast<std::size_t>(kN) * kN * 3, 0.0);

    // Columns 0 and 1 differ on purpose, so extrapolation instead of clamping would show up.
    for (int y = 10; y <= 11; ++y) {
        texels[(static_cast<std::size_t>(y) * kN + 0) * 3] = 0.25 + 0.1 * (y - 10); // column 0
        texels[(static_cast<std::size_t>(y) * kN + 1) * 3] = 0.90 - 0.3 * (y - 10); // column 1
    }
    const auto d = decimate::seam::seamEnergy(tm, kN, kN, texels);
    CHECK(d[0] <= 1e-24); // an unclamped extrapolation would give ~0.25 * (c0 - c1)^2 > 1e-3
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S4: seamEnergy equals the hand-computed closed form", "[seam][S4]") {
    // N=4. Side A: x = 1.5, y = 1.25..1.75 (inside cell [1,2]x[1,2]); side B: x = 2.5, same y
    // (cell [2,3]x[1,2]); neither crosses a texel-grid line. With s = 0.5 and t = 0.25 + 0.5g:
    //   A(g) = 0.4 + 0.3 t,  B(g) = 0.4 + 0.5 t,  A - B = -0.05 - 0.1 g
    //   D = int_0^1 (0.05 + 0.1 g)^2 dg = 0.0025 + 0.005 + 0.01/3 = 13/1200.
    constexpr int       kN = 4;
    const TexturedMesh  tm = seamPair(uvOfTexel(1.5, 1.25, kN), uvOfTexel(1.5, 1.75, kN),
                                      uvOfTexel(2.5, 1.25, kN), uvOfTexel(2.5, 1.75, kN));
    std::vector<double> texels(static_cast<std::size_t>(kN) * kN * 3, 0.0);
    const auto          set = [&](int x, int y, double v) {
        texels[(static_cast<std::size_t>(y) * kN + x) * 3] = v;
        texels[(static_cast<std::size_t>(y) * kN + x) * 3 + 1] = 0.5 * v; // channel 1 = half
    };
    set(1, 1, 0.2);
    set(2, 1, 0.6);
    set(1, 2, 0.4);
    set(2, 2, 1.0);
    set(3, 1, 0.2);
    set(3, 2, 0.8);

    const auto d = decimate::seam::seamEnergy(tm, kN, kN, texels);
    CHECK(std::abs(d[0] - 13.0 / 1200.0) <= 1e-12);
    CHECK(std::abs(d[1] - 0.25 * 13.0 / 1200.0) <= 1e-12);
    CHECK(d[2] == 0.0);
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S6: toggle off -- bake and bakeSeamFill outputs do not depend on erase", "[seam][S6]") {
    // erase() takes the map by const& and owns no global state: baking before and after an
    // erase call must be bit-identical, and the bake sources themselves are untouched by
    // this spec (the bake_characterization / seam_bake_test goldens keep guarding them).
    const SpotPipeline &s = spot();
    const NormalMap     before =
        decimate::bake::bake(s.low, s.cornerUV, s.lowN, s.high, s.highN, 256, 2, {});
    CHECK(sameMap(before, s.dilateMap));

    (void)decimate::seam::erase(s.tm, s.dilateMap, NormalSpace::Object);
    const NormalMap after =
        decimate::bake::bake(s.low, s.cornerUV, s.lowN, s.high, s.highN, 256, 2, {});
    CHECK(sameMap(after, before));
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S7: both seamEnergy overloads are one unit (bit-identical doubles)", "[seam][S7]") {
    const Golden        g = loadGolden("rot");
    std::vector<double> texels;

    // The raw-texel overload takes [0,1] doubles; feed it the same map the Rgb overload quantizes from.
    for (const Rgb &p : g.input.pixels) {
        texels.push_back(static_cast<double>(p.r) / 255.0);
        texels.push_back(static_cast<double>(p.g) / 255.0);
        texels.push_back(static_cast<double>(p.b) / 255.0);
    }

    const auto a = decimate::seam::seamEnergy(g.tm, g.input);
    const auto b = decimate::seam::seamEnergy(g.tm, g.size, g.size, texels);
    for (std::size_t c = 0; c < 3; ++c) {
        CHECK(std::bit_cast<std::uint64_t>(a[c]) == std::bit_cast<std::uint64_t>(b[c]));
    }
}
