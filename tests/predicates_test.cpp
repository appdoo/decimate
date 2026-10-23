// This file is part of the decimate library; see decimate.h for version/license details.

/// @file predicates_test.cpp
/// @brief Spec 23: the robust predicates return the exact sign.  Reference signs come from integer
///        arithmetic (`__int128`) on lattice coordinates, never from the code under test.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <decimate/math/predicates.h>
#include <limits>
#include <random>
#include <utility>
#include <vector>

using decimate::math::Vec2d;
using decimate::math::Vec2f;
using decimate::math::Vec3d;

namespace {

using i128 = __int128;

/// Sign of an exact integer.
int signOf(i128 v) {
    return (v > 0) - (v < 0);
}


/// Sign of a double.
int signOf(double v) {
    return (v > 0.0) - (v < 0.0);
}


/// The coordinate x as an integer count of 2^-53 (exact for |x| < 2^10 on that lattice).
std::int64_t units53(double x) {
    return static_cast<std::int64_t>(std::ldexp(x, 53));
}


/// The naive double determinant of the same three points (what the old local copies computed).
double naiveOrient2d(Vec2d a, Vec2d b, Vec2d c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}


/// Exact sign of orient2d for points that are integer multiples of 2^-53.
int exactOrient2d(Vec2d a, Vec2d b, Vec2d c) {
    const i128 abx = static_cast<i128>(units53(b.x)) - units53(a.x);
    const i128 aby = static_cast<i128>(units53(b.y)) - units53(a.y);
    const i128 acx = static_cast<i128>(units53(c.x)) - units53(a.x);
    const i128 acy = static_cast<i128>(units53(c.y)) - units53(a.y);
    return signOf(abx * acy - aby * acx);
}

} // namespace


TEST_CASE("Spec 23: orient2d is exact on the near-collinear ulp grid; naive double is not", "[predicates]") {
    // Shewchuk's image: p moves over a 256x256 ulp grid around (0.5, 0.5) while q, r span the line y = x.
    const Vec2d  q { 12.0, 12.0 };
    const Vec2d  r { 24.0, 24.0 };
    const double ulp = std::ldexp(1.0, -53);
    int          naiveWrong = 0;
    int          zeros = 0;
    int          positives = 0;
    int          negatives = 0;
    for (int i = 0; i < 256; ++i) {
        for (int j = 0; j < 256; ++j) {
            const Vec2d p { 0.5 + i * ulp, 0.5 + j * ulp };
            const int   want = exactOrient2d(p, q, r);
            REQUIRE(signOf(decimate::math::orient2d(p, q, r)) == want);
            naiveWrong += signOf(naiveOrient2d(p, q, r)) != want ? 1 : 0;
            zeros += want == 0;
            positives += want > 0;
            negatives += want < 0;
        }
    }
    CHECK(naiveWrong >= 1);
    CHECK(zeros == 256); // the diagonal i == j lies exactly on the line
    CHECK(positives == negatives);
}


TEST_CASE("Spec 23: orient2d is exact on a perturbed float triple", "[predicates]") {
    // Float inputs widen exactly: three nearly collinear floats, the middle one nudged by one float ulp.
    const Vec2f a { 0.1f, 0.1f };
    const Vec2f c { 0.9f, 0.9f };
    const float mid = 0.5f;
    const float up = std::nextafter(mid, 1.0f);
    CHECK(decimate::math::orient2d(a, Vec2f { mid, mid }, c) == 0.0);
    CHECK(decimate::math::orient2d(a, Vec2f { mid, up }, c) < 0.0); // b above the line a->c: a, b, c turn clockwise
    CHECK(decimate::math::orient2d(a, Vec2f { up, mid }, c) > 0.0);
}


TEST_CASE("Spec 23: an FMA-contracted evaluation gets the sign wrong, orient2d does not", "[predicates]") {
    // a = (u, v), b = (2u, 2v), c = origin are exactly collinear, but u*v is not representable: the
    // determinant u*(2v) - v*(2u) is 0 only when both products are rounded alike.  One fused
    // multiply-add leaves the rounding error of the other product behind.
    const double u = 1.0 + std::ldexp(1.0, -30);
    const double v = 1.0 + std::ldexp(1.0, -31);
    const Vec2d  a { u, v };
    const Vec2d  b { 2.0 * u, 2.0 * v };
    const Vec2d  c { 0.0, 0.0 };
    const double contracted = std::fma(a.x - c.x, b.y - c.y, -((a.y - c.y) * (b.x - c.x)));
    CHECK(contracted != 0.0);
    CHECK(decimate::math::orient2d(a, b, c) == 0.0);
    CHECK(decimate::math::orient2d(c, a, b) == 0.0);
}


TEST_CASE("Spec 23: orient2d on non-finite input is non-finite, never a spurious sign", "[predicates]") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double r = decimate::math::orient2d(Vec2d { nan, 0.0 }, Vec2d { 1.0, 0.0 }, Vec2d { 0.0, 1.0 });
    CHECK_FALSE(r > 0.0);
    CHECK_FALSE(r < 0.0);
    CHECK_FALSE(r == 0.0);
}


TEST_CASE("Spec 23: orient3d is exact near coplanarity", "[predicates]") {
    // x, y on a 2^-10 lattice, z on a 2^-53 lattice: scaling each axis by its own lattice gives integer
    // coordinates and an __int128 determinant.  The plane z = 0.5 + (x + y) / 8 is exact; d is moved
    // off it by whole ulps of z.
    std::mt19937_64                    rng(23);
    std::uniform_int_distribution<int> lat(0, 4095);
    std::uniform_int_distribution<int> ulps(-3, 3);
    const double                       xy = std::ldexp(1.0, -10);
    const double                       ulp = std::ldexp(1.0, -53);
    int                                naiveWrong = 0;
    int                                nonzero = 0;
    for (int trial = 0; trial < 20000; ++trial) {
        std::array<Vec3d, 4>        p {};
        std::array<std::int64_t, 4> zi {};
        std::array<std::int64_t, 4> xi {};
        std::array<std::int64_t, 4> yi {};
        for (int k = 0; k < 4; ++k) {
            xi[k] = lat(rng);
            yi[k] = lat(rng);
            const double x = xi[k] * xy;
            const double y = yi[k] * xy;
            const int    off = k == 3 ? ulps(rng) : 0;
            const double z = 0.5 + (x + y) / 8.0 + off * ulp;
            p[k] = { x, y, z };
            zi[k] = units53(z);
        }
        const i128 adx = xi[0] - xi[3], bdx = xi[1] - xi[3], cdx = xi[2] - xi[3];
        const i128 ady = yi[0] - yi[3], bdy = yi[1] - yi[3], cdy = yi[2] - yi[3];
        const i128 adz = zi[0] - zi[3], bdz = zi[1] - zi[3], cdz = zi[2] - zi[3];
        const int  want = signOf(adz * (bdx * cdy - bdy * cdx) + bdz * (cdx * ady - cdy * adx) +
                                 cdz * (adx * bdy - ady * bdx));
        REQUIRE(signOf(decimate::math::orient3d(p[0], p[1], p[2], p[3])) == want);

        const double nx = (p[0].x - p[3].x), mx = (p[1].x - p[3].x), kx = (p[2].x - p[3].x);
        const double ny = (p[0].y - p[3].y), my = (p[1].y - p[3].y), ky = (p[2].y - p[3].y);
        const double nz = (p[0].z - p[3].z), mz = (p[1].z - p[3].z), kz = (p[2].z - p[3].z);
        const double naive = nz * (mx * ky - my * kx) + mz * (kx * ny - ky * nx) + kz * (nx * my - ny * mx);
        naiveWrong += signOf(naive) != want ? 1 : 0;
        nonzero += want != 0 ? 1 : 0;
    }
    CHECK(nonzero > 1000);
    INFO("naive orient3d wrong on " << naiveWrong << " of 20000 near-coplanar quadruples");
    CHECK(naiveWrong >= 1);
}


TEST_CASE("Spec 23: orient3d orientation convention and coplanar zero", "[predicates]") {
    const Vec3d a { 0.0, 0.0, 0.0 }, b { 1.0, 0.0, 0.0 }, c { 0.0, 1.0, 0.0 };
    CHECK(decimate::math::orient3d(a, b, c, Vec3d { 0.2, 0.2, -1.0 }) > 0.0); // below the CCW triangle
    CHECK(decimate::math::orient3d(a, b, c, Vec3d { 0.2, 0.2, 1.0 }) < 0.0);
    CHECK(decimate::math::orient3d(a, b, c, Vec3d { 7.0, -3.0, 0.0 }) == 0.0);
}


/// The lattice point (2 + i)^a (2 - i)^(n - a): every such point has norm 5^n, so they all lie on one circle.
static std::pair<std::int64_t, std::int64_t> gaussianPoint(int a, int n) {
    std::int64_t re = 1;
    std::int64_t im = 0;
    for (int k = 0; k < n; ++k) {
        const std::int64_t f = k < a ? 1 : -1; // imaginary part of the factor (2 + f i)
        const std::int64_t r = re * 2 - im * f;
        im = re * f + im * 2;
        re = r;
    }
    return { re, im };
}


TEST_CASE("Spec 23: incircle is exact on a cocircular lattice", "[predicates]") {
    // Lattice points of norm 5^24 (coordinates ~ 2^28) are exactly cocircular; d moves around one of
    // them by whole lattice units.  The integer coordinates make the __int128 determinant the reference.
    constexpr int n = 24;
    auto          pa = gaussianPoint(24, n);
    auto          pb = gaussianPoint(17, n);
    auto          pc = gaussianPoint(9, n);
    const auto    pd = gaussianPoint(3, n);
    const i128    orient = static_cast<i128>(pb.first - pa.first) * (pc.second - pa.second) -
                           static_cast<i128>(pb.second - pa.second) * (pc.first - pa.first);
    if (orient < 0) {
        std::swap(pb, pc); // incircle wants a, b, c counter-clockwise
    }
    int naiveWrong = 0;
    int zeros = 0;
    for (int i = -8; i <= 8; ++i) {
        for (int j = -8; j <= 8; ++j) {
            const std::int64_t dx = pd.first + i, dy = pd.second + j;
            const i128         adx = pa.first - dx, bdx = pb.first - dx, cdx = pc.first - dx;
            const i128         ady = pa.second - dy, bdy = pb.second - dy, cdy = pc.second - dy;
            const i128         aLift = adx * adx + ady * ady, bLift = bdx * bdx + bdy * bdy, cLift = cdx * cdx + cdy * cdy;
            const int          want = signOf(aLift * (bdx * cdy - cdx * bdy) + bLift * (cdx * ady - adx * cdy) +
                                             cLift * (adx * bdy - bdx * ady));
            const Vec2d        a { double(pa.first), double(pa.second) }, b { double(pb.first), double(pb.second) };
            const Vec2d        c { double(pc.first), double(pc.second) }, d { double(dx), double(dy) };
            REQUIRE(signOf(decimate::math::incircle(a, b, c, d)) == want);
            zeros += want == 0;

            const double fadx = a.x - d.x, fbdx = b.x - d.x, fcdx = c.x - d.x;
            const double fady = a.y - d.y, fbdy = b.y - d.y, fcdy = c.y - d.y;
            const double naive = (fadx * fadx + fady * fady) * (fbdx * fcdy - fcdx * fbdy) +
                                 (fbdx * fbdx + fbdy * fbdy) * (fcdx * fady - fadx * fcdy) +
                                 (fcdx * fcdx + fcdy * fcdy) * (fadx * fbdy - fbdx * fady);
            naiveWrong += signOf(naive) != want ? 1 : 0;
        }
    }
    CHECK(zeros == 1); // d == its lattice point, exactly on the circle
    INFO("naive incircle wrong on " << naiveWrong << " of 289 cocircular-neighbour points");
    CHECK(naiveWrong >= 1);
}


TEST_CASE("Spec 23: incircle inside, outside and clockwise flip", "[predicates]") {
    const Vec2d a { 1.0, 0.0 }, b { 0.0, 1.0 }, c { -1.0, 0.0 };
    CHECK(decimate::math::incircle(a, b, c, Vec2d { 0.0, 0.0 }) > 0.0);
    CHECK(decimate::math::incircle(a, b, c, Vec2d { 2.0, 2.0 }) < 0.0);
    CHECK(decimate::math::incircle(a, b, c, Vec2d { 0.0, -1.0 }) == 0.0);
    CHECK(decimate::math::incircle(a, c, b, Vec2d { 0.0, 0.0 }) < 0.0);
}


TEST_CASE("Spec 23: segmentsCross keeps the proper-crossing semantics", "[predicates]") {
    const Vec2f o { 0.0f, 0.0f };
    const Vec2f p { 2.0f, 2.0f };
    using decimate::math::segmentsCross;
    CHECK(segmentsCross(o, p, Vec2f { 0.0f, 2.0f }, Vec2f { 2.0f, 0.0f }));                          // X
    CHECK_FALSE(segmentsCross(o, p, Vec2f { 1.0f, 1.0f }, Vec2f { 3.0f, 0.0f }));                    // T: end point on the line a-b
    CHECK_FALSE(segmentsCross(o, p, Vec2f { 0.0f, 0.0f }, Vec2f { 2.0f, 0.0f }));                    // shared end point
    CHECK_FALSE(segmentsCross(o, p, Vec2f { 1.0f, 1.0f }, Vec2f { 3.0f, 3.0f }));                    // collinear overlap
    CHECK_FALSE(segmentsCross(o, p, Vec2f { 3.0f, 3.0f }, Vec2f { 4.0f, 4.0f }));                    // collinear disjoint
    CHECK_FALSE(segmentsCross(o, Vec2f { 2.0f, 0.0f }, Vec2f { 0.0f, 1.0f }, Vec2f { 2.0f, 1.0f })); // parallel
    CHECK_FALSE(segmentsCross(o, p, Vec2f { 3.0f, 0.0f }, Vec2f { 4.0f, 0.0f }));                    // lines cross outside the segments
    CHECK(segmentsCross(Vec2d { 0.0, 0.0 }, Vec2d { 2.0, 2.0 }, Vec2d { 0.0, 2.0 }, Vec2d { 2.0, 0.0 }));
}


TEST_CASE("Spec 23: segmentsCross agrees with the naive float test away from degeneracy", "[predicates]") {
    // The old local copy: four double determinants of float inputs, strict sign changes on both lines.
    const auto oldCross = [](Vec2f a, Vec2f b, Vec2f c, Vec2f d) {
        const auto det = [](Vec2f p, Vec2f q, Vec2f r) {
            return (double(q.x) - p.x) * (double(r.y) - p.y) - (double(q.y) - p.y) * (double(r.x) - p.x);
        };
        const double d1 = det(c, d, a), d2 = det(c, d, b), d3 = det(a, b, c), d4 = det(a, b, d);
        return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
    };
    std::mt19937                          rng(7);
    std::uniform_real_distribution<float> u(-1.0f, 1.0f);
    int                                   crossings = 0;
    for (int i = 0; i < 100000; ++i) {
        const Vec2f a { u(rng), u(rng) }, b { u(rng), u(rng) }, c { u(rng), u(rng) }, d { u(rng), u(rng) };
        const bool  want = oldCross(a, b, c, d);
        REQUIRE(decimate::math::segmentsCross(a, b, c, d) == want);
        crossings += want ? 1 : 0;
    }
    CHECK(crossings > 10000);
}


TEST_CASE("Spec 23: signedArea2 of a polygon is twice the area, positive counter-clockwise", "[predicates]") {
    const std::vector<Vec2d> sq { { 0.0, 0.0 }, { 4.0, 0.0 }, { 4.0, 4.0 }, { 0.0, 4.0 } };
    CHECK(decimate::math::signedArea2(sq) == 32.0);
    const std::vector<Vec2d> cw(sq.rbegin(), sq.rend());
    CHECK(decimate::math::signedArea2(cw) == -32.0);
    const std::vector<Vec2f> tri { { 0.0f, 0.0f }, { 3.0f, 0.0f }, { 0.0f, 2.0f } };
    CHECK(decimate::math::signedArea2(tri) == 6.0);
}


TEST_CASE("Spec 23: predicates cost nothing measurable on generic input", "[.bench][predicates]") {
    std::mt19937                           rng(1);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<Vec2d>                     pts(3 * 1000);
    for (Vec2d &p : pts) {
        p = { u(rng), u(rng) };
    }
    double acc = 0.0;
    for (int rep = 0; rep < 20000; ++rep) {
        for (std::size_t i = 0; i < pts.size(); i += 3) {
            acc += decimate::math::orient2d(pts[i], pts[i + 1], pts[i + 2]);
        }
    }
    CHECK(std::isfinite(acc));
}
