// This file is part of the decimate library; see decimate.h for version/license details.
/// @file tangent_space_test.cpp
/// @brief Spec 26: the tangent-space bake is written in MikkTSpace, the basis engines decode with.
///
/// a [S] flat quad, UV = XY: T = (1,0,0), w = +1; U mirrored: T = (-1,0,0), w = -1
/// b [S] cow: corners sharing vertex + UV (+ normal) get identical tangents
/// c [I] cow round-trip: decode every interior texel with the MikkTSpace convention, compare with the
///       object-space bake of the same atlas (max angle bound)
/// Object space is covered by the existing bake goldens, which stay byte-identical.

#include "bake_bvh_common.h"
#include "model_fixture.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <decimate/bake/tangents.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <tuple>

using namespace decimate;
using decimate::bake::NormalMap;
using decimate::bake::Rgb;
using decimate::math::Vec2f;
using decimate::math::Vec3f;
using decimate::math::Vec4f;

namespace {

/// Unit square in the XY plane facing +z, two triangles; `uMirror` flips U.
MeshData quad(std::vector<Vec2f> &uv, std::vector<Vec3f> &nrm, bool uMirror) {
    MeshData m;
    m.positions = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 } };
    m.indices = { 0, 1, 2, 0, 2, 3 };
    nrm.assign(4, { 0.0f, 0.0f, 1.0f });
    for (const int i : m.indices) {
        const Vec3f &p = m.positions[static_cast<std::size_t>(i)];
        uv.push_back({ uMirror ? 1.0f - p.x : p.x, p.y });
    }
    return m;
}


/// The cow reduced to 1500 faces, charted at cap 64, parametrized on a 256 grid (as the bake tests).
bvhtest::Pipeline cowPipeline() {
    return bvhtest::buildPipeline(modelOrSkip("cow.obj"), 1500, 64, 256);
}


/// Decoded unit normal of a texel: inverse of rangeCompress (round to nearest, decode c / 127 - 1).
Vec3f unpack(const Rgb &c) {
    Vec3f n { static_cast<float>(c.r) / 127.0f - 1.0f, static_cast<float>(c.g) / 127.0f - 1.0f,
              static_cast<float>(c.b) / 127.0f - 1.0f };
    n.normalize();
    return n;
}


/// What an engine does with a tangent-space texel `tn` at barycentrics (u,v,w) of corner data (t, n):
/// interpolate T and N unnormalised, B = sign * cross(N, T), n = normalize(tn.x*T + tn.y*B + tn.z*N).
Vec3f engineDecode(const Vec3f &tn, float u, float v, float w, const Vec4f t[3], const Vec3f n[3]) {
    const Vec3f T = Vec3f { t[0].x, t[0].y, t[0].z } * u + Vec3f { t[1].x, t[1].y, t[1].z } * v +
                    Vec3f { t[2].x, t[2].y, t[2].z } * w;
    const Vec3f N = n[0] * u + n[1] * v + n[2] * w;
    const float sign = t[0].w * u + t[1].w * v + t[2].w * w < 0.0f ? -1.0f : 1.0f;
    Vec3f       r = T * tn.x + N.cross(T) * (sign * tn.y) + N * tn.z;
    r.normalize();
    return r;
}


/// Angle in degrees between two unit vectors.
double angleDeg(const Vec3f &a, const Vec3f &b) {
    return std::acos(std::fmin(1.0f, std::fmax(-1.0f, a.dot(b)))) * 180.0 / 3.14159265358979323846;
}

} // namespace


TEST_CASE("Spec 26: flat quad with UV = XY has T = +x and w = +1", "[tangents][spec26][S]") {
    std::vector<Vec2f> uv;
    std::vector<Vec3f> nrm;
    const MeshData     m = quad(uv, nrm, false);
    const auto         t = tangents(m, uv, nrm);
    REQUIRE(t.size() == 6);
    for (const Vec4f &c : t) {
        CHECK(c.x == Catch::Approx(1.0f).margin(1e-6));
        CHECK(c.y == Catch::Approx(0.0f).margin(1e-6));
        CHECK(c.z == Catch::Approx(0.0f).margin(1e-6));
        CHECK(c.w == 1.0f);
    }
}


TEST_CASE("Spec 26: mirrored U flips the sign and T follows +u", "[tangents][spec26][S]") {
    std::vector<Vec2f> uv;
    std::vector<Vec3f> nrm;
    const MeshData     m = quad(uv, nrm, true);
    const auto         t = tangents(m, uv, nrm);
    REQUIRE(t.size() == 6);
    for (const Vec4f &c : t) {
        CHECK(c.x == Catch::Approx(-1.0f).margin(1e-6));
        CHECK(c.y == Catch::Approx(0.0f).margin(1e-6));
        CHECK(c.w == -1.0f);
    }
}


TEST_CASE("Spec 26: cow corners sharing vertex and UV share one tangent", "[tangents][spec26][S]") {
    const auto p = cowPipeline();
    const auto t = tangents(p.low, p.cornerUV, p.lowN);
    REQUIRE(t.size() == p.low.indices.size());

    // Group by (vertex, uv): the same position, normal (per-vertex) and UV, so MikkTSpace welds them.
    std::map<std::tuple<int, float, float>, Vec4f> first;
    int                                            shared = 0;
    int                                            differing = 0;
    for (std::size_t c = 0; c < t.size(); c++) {
        const auto key = std::make_tuple(p.low.indices[c], p.cornerUV[c].x, p.cornerUV[c].y);
        const auto [it, fresh] = first.emplace(key, t[c]);
        if (fresh) {
            continue;
        }
        shared++;
        const Vec4f &a = it->second;
        differing += std::fabs(a.x - t[c].x) + std::fabs(a.y - t[c].y) + std::fabs(a.z - t[c].z) > 1e-5f || a.w != t[c].w;
    }
    CHECK(shared > 1000);
    CHECK(differing == 0);
}


TEST_CASE("Spec 26: tangent-space bake decodes back to the object-space bake", "[tangents][spec26][I]") {
    auto p = cowPipeline();

    // Mirror U on every odd chart so the sign w = -1 is exercised (the charts may then overlap in the
    // atlas; both bakes see the same UVs and the same owner order, so the comparison still holds).
    for (std::size_t c = 0; c < p.cornerUV.size(); c++) {
        if (p.faceToChart[c / 3] % 2 == 1) {
            p.cornerUV[c].x = 1.0f - p.cornerUV[c].x;
        }
    }
    const auto tg = tangents(p.low, p.cornerUV, p.lowN);
    const auto minus = std::count_if(tg.begin(), tg.end(), [](const Vec4f &c) { return c.w < 0.0f; });
    REQUIRE(minus > 100);

    // dilatePasses 0: only texels inside a face are written, so every non-background texel is interior.
    const NormalMap object = bake::bake(p.low, p.cornerUV, p.lowN, p.high, p.highN, 256, 0);
    const NormalMap tangent = bake::bake(p.low, p.cornerUV, p.lowN, p.high, p.highN, 256, 0, tg);

    // Owner of each texel: the last face whose UV triangle covers its centre (the bake's overwrite order).
    const int          tex = 256;
    const float        ts = static_cast<float>(tex);
    std::vector<int>   owner(static_cast<std::size_t>(tex) * tex, -1);
    std::vector<float> bc(owner.size() * 3);
    const int          nf = static_cast<int>(p.low.indices.size() / 3);
    for (int f = 0; f < nf; f++) {
        const Vec2f a = p.cornerUV[3 * f] * ts;
        const Vec2f b = p.cornerUV[3 * f + 1] * ts;
        const Vec2f c = p.cornerUV[3 * f + 2] * ts;
        const float area = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
        if (!(std::fabs(area) > 0.0f)) {
            continue;
        }
        for (int y = 0; y < tex; y++) {
            for (int x = 0; x < tex; x++) {
                const float px = static_cast<float>(x) + 0.5f;
                const float py = static_cast<float>(y) + 0.5f;
                const float u = ((b.x - px) * (c.y - py) - (c.x - px) * (b.y - py)) / area;
                const float v = ((c.x - px) * (a.y - py) - (a.x - px) * (c.y - py)) / area;
                const float w = 1.0f - u - v;
                if (u < 0.0f || v < 0.0f || w < 0.0f) {
                    continue;
                }
                const std::size_t i = static_cast<std::size_t>(y) * tex + x;
                owner[i] = f;
                bc[3 * i] = u;
                bc[3 * i + 1] = v;
                bc[3 * i + 2] = w;
            }
        }
    }

    double maxDeg = 0.0;
    double sumDeg = 0.0;
    int    texels = 0;
    for (std::size_t i = 0; i < owner.size(); i++) {
        const Rgb bg { 128, 128, 128 };
        if (owner[i] < 0 || object.pixels[i] == bg || tangent.pixels[i] == bg) {
            continue;
        }
        const int    f = owner[i];
        const Vec4f  ct[3] = { tg[3 * f], tg[3 * f + 1], tg[3 * f + 2] };
        const Vec3f  cn[3] = { p.lowN[p.low.indices[3 * f]], p.lowN[p.low.indices[3 * f + 1]],
                               p.lowN[p.low.indices[3 * f + 2]] };
        const Vec3f  got = engineDecode(unpack(tangent.pixels[i]), bc[3 * i], bc[3 * i + 1], bc[3 * i + 2], ct, cn);
        const double deg = angleDeg(got, unpack(object.pixels[i]));
        maxDeg = std::fmax(maxDeg, deg);
        sumDeg += deg;
        texels++;
    }
    const double meanDeg = sumDeg / texels;
    CAPTURE(texels, maxDeg, meanDeg);
    REQUIRE(texels > 10000);
    CHECK(maxDeg <= 0.63);
}
