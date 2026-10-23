// This file is part of the decimate library; see decimate.h for version/license details.

/// @file predicates.cpp
/// @brief Shewchuk's adaptive robust predicates, ported without global state (see predicates.h).
///
/// The error-free transformations below need every product and sum rounded on its own, so FP
/// contraction is off for this file (also set on the command line in src/CMakeLists.txt).

#include <decimate/math/predicates.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#pragma STDC FP_CONTRACT OFF

namespace decimate::math {

static_assert(std::numeric_limits<double>::is_iec559, "the error bounds assume IEEE-754 double");
static_assert(std::numeric_limits<double>::digits == 53, "the error bounds assume 53-bit doubles");

/// Shewchuk's machine epsilon (half an ulp of 1) and the error bounds derived from it.
constexpr double kEps = 1.0 / 9007199254740992.0; // 2^-53
constexpr double kResultErrBound = (3.0 + 8.0 * kEps) * kEps;
constexpr double kCcwErrBoundA = (3.0 + 16.0 * kEps) * kEps;
constexpr double kCcwErrBoundB = (2.0 + 12.0 * kEps) * kEps;
constexpr double kCcwErrBoundC = (9.0 + 64.0 * kEps) * kEps * kEps;
constexpr double kO3dErrBoundA = (7.0 + 56.0 * kEps) * kEps;
constexpr double kIccErrBoundA = (10.0 + 96.0 * kEps) * kEps;

/// Dekker's splitter 2^27 + 1, used when there is no hardware fused multiply-add.
constexpr double kSplitter = 134217729.0;


/// a + b = x + y exactly (x the rounded sum, y the rounding error).
static void twoSum(double a, double b, double &x, double &y) {
    x = a + b;
    const double bVirt = x - a;
    const double aVirt = x - bVirt;
    y = (a - aVirt) + (b - bVirt);
}


/// As twoSum, valid when |a| >= |b| (or a == 0).
static void fastTwoSum(double a, double b, double &x, double &y) {
    x = a + b;
    y = b - (x - a);
}


/// a - b = x + y exactly.
static void twoDiff(double a, double b, double &x, double &y) {
    x = a - b;
    const double bVirt = a - x;
    const double aVirt = x + bVirt;
    y = (a - aVirt) + (bVirt - b);
}


/// The rounding error of x = a - b.
static double diffTail(double a, double b, double x) {
    const double bVirt = a - x;
    const double aVirt = x + bVirt;
    return (a - aVirt) + (bVirt - b);
}


/// a * b = x + y exactly.  One fused multiply-add where the hardware has it (Apple Silicon, x86 FMA3),
/// Dekker's split otherwise (Wasm).  An explicit fma stays exact under -ffp-contract=off.
static void twoProduct(double a, double b, double &x, double &y) {
    x = a * b;
#if defined(__FP_FAST_FMA) && !defined(DECIMATE_PREDICATES_DEKKER)
    y = std::fma(a, b, -x);
#else
    const double aC = kSplitter * a;
    const double aHi = aC - (aC - a);
    const double aLo = a - aHi;
    const double bC = kSplitter * b;
    const double bHi = bC - (bC - b);
    const double bLo = b - bHi;
    y = aLo * bLo - (((x - aHi * bHi) - aLo * bHi) - aHi * bLo);
#endif
}


/// An expansion: a sum of non-overlapping doubles, smallest magnitude first; the last component
/// carries the sign of the whole.  Only the exact fallback uses these, so a vector is fine.
using Expansion = std::vector<double>;


/// e + b, zero components eliminated (Shewchuk's grow_expansion_zeroelim).
static Expansion growExpansion(const Expansion &e, double b) {
    Expansion h;
    double    q = b;
    for (const double c : e) {
        double qNew = 0.0;
        double hh = 0.0;
        twoSum(q, c, qNew, hh);
        q = qNew;
        if (hh != 0.0) {
            h.push_back(hh);
        }
    }
    if (q != 0.0 || h.empty()) {
        h.push_back(q);
    }
    return h;
}


/// e + f.
static Expansion addExpansions(Expansion e, const Expansion &f) {
    for (const double c : f) {
        e = growExpansion(e, c);
    }
    return e;
}


/// e * b (Shewchuk's scale_expansion_zeroelim).
static Expansion scaleExpansion(const Expansion &e, double b) {
    Expansion h;
    double    q = 0.0;
    double    hh = 0.0;
    twoProduct(e[0], b, q, hh);
    if (hh != 0.0) {
        h.push_back(hh);
    }
    for (std::size_t i = 1; i < e.size(); ++i) {
        double p1 = 0.0;
        double p0 = 0.0;
        double s = 0.0;
        twoProduct(e[i], b, p1, p0);
        twoSum(q, p0, s, hh);
        if (hh != 0.0) {
            h.push_back(hh);
        }
        fastTwoSum(p1, s, q, hh);
        if (hh != 0.0) {
            h.push_back(hh);
        }
    }
    if (q != 0.0 || h.empty()) {
        h.push_back(q);
    }
    return h;
}


/// e * f.
static Expansion mulExpansions(const Expansion &e, const Expansion &f) {
    Expansion r { 0.0 };
    for (const double c : f) {
        r = addExpansions(std::move(r), scaleExpansion(e, c));
    }
    return r;
}


/// -e.
static Expansion negated(Expansion e) {
    for (double &c : e) {
        c = -c;
    }
    return e;
}


/// a - b as an exact two-component expansion.
static Expansion diff(double a, double b) {
    double x = 0.0;
    double y = 0.0;
    twoDiff(a, b, x, y);
    return { y, x };
}


/// ab - cd, exactly.
static Expansion crossTerm(const Expansion &a, const Expansion &b, const Expansion &c, const Expansion &d) {
    return addExpansions(mulExpansions(a, b), negated(mulExpansions(c, d)));
}


/// Exact orient2d of the widened inputs: the determinant of the exact differences.
static double orient2dExact(const double *pa, const double *pb, const double *pc) {
    const Expansion acx = diff(pa[0], pc[0]);
    const Expansion acy = diff(pa[1], pc[1]);
    const Expansion bcx = diff(pb[0], pc[0]);
    const Expansion bcy = diff(pb[1], pc[1]);
    return crossTerm(acx, bcy, acy, bcx).back();
}


/// Stages B and C of orient2d: refine with the rounding tails of the four coordinate differences.
/// @p detSum is the sum of the magnitudes of the two products (the filter's permanent).
static double orient2dAdapt(const double *pa, const double *pb, const double *pc, double detSum) {
    const double acx = pa[0] - pc[0];
    const double bcx = pb[0] - pc[0];
    const double acy = pa[1] - pc[1];
    const double bcy = pb[1] - pc[1];

    // Stage B: the product pair is exact; only the four differences are still rounded.
    double detLeft = 0.0;
    double detLeftTail = 0.0;
    double detRight = 0.0;
    double detRightTail = 0.0;
    twoProduct(acx, bcy, detLeft, detLeftTail);
    twoProduct(acy, bcx, detRight, detRightTail);
    double b0 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0, t = 0.0, u = 0.0, v = 0.0, w = 0.0;
    twoDiff(detLeftTail, detRightTail, t, b0);
    twoSum(detLeft, t, u, v);
    twoDiff(v, detRight, w, b1);
    twoSum(u, w, b3, b2);
    double det = b0 + b1 + b2 + b3;
    double errBound = kCcwErrBoundB * detSum;
    if (det >= errBound || -det >= errBound) {
        return det;
    }

    // Stage C: first-order correction by the difference tails.
    const double acxTail = diffTail(pa[0], pc[0], acx);
    const double bcxTail = diffTail(pb[0], pc[0], bcx);
    const double acyTail = diffTail(pa[1], pc[1], acy);
    const double bcyTail = diffTail(pb[1], pc[1], bcy);
    if (acxTail == 0.0 && acyTail == 0.0 && bcxTail == 0.0 && bcyTail == 0.0) {
        return det;
    }
    errBound = kCcwErrBoundC * detSum + kResultErrBound * std::abs(det);
    det += (acx * bcyTail + bcy * acxTail) - (acy * bcxTail + bcx * acyTail);
    if (det >= errBound || -det >= errBound) {
        return det;
    }
    return orient2dExact(pa, pb, pc);
}


static double orient2dRaw(const double *pa, const double *pb, const double *pc) {
    const double detLeft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
    const double detRight = (pa[1] - pc[1]) * (pb[0] - pc[0]);
    const double det = detLeft - detRight;

    // Opposite signs of the two products: the difference cannot lose its sign to rounding.
    double detSum = 0.0;
    if (detLeft > 0.0) {
        if (detRight <= 0.0) {
            return det;
        }
        detSum = detLeft + detRight;
    } else if (detLeft < 0.0) {
        if (detRight >= 0.0) {
            return det;
        }
        detSum = -detLeft - detRight;
    } else {
        return det;
    }
    const double errBound = kCcwErrBoundA * detSum;
    if (det >= errBound || -det >= errBound) {
        return det;
    }
    return orient2dAdapt(pa, pb, pc, detSum);
}


double orient2d(Vec2d a, Vec2d b, Vec2d c) noexcept {
    const double pa[2] = { a.x, a.y };
    const double pb[2] = { b.x, b.y };
    const double pc[2] = { c.x, c.y };
    return orient2dRaw(pa, pb, pc);
}


double orient2d(Vec2f a, Vec2f b, Vec2f c) noexcept {
    return orient2d(Vec2d { a.x, a.y }, Vec2d { b.x, b.y }, Vec2d { c.x, c.y });
}


/// Exact orient3d from the exact differences to d.
static double orient3dExact(Vec3d a, Vec3d b, Vec3d c, Vec3d d) {
    const Expansion adx = diff(a.x, d.x), ady = diff(a.y, d.y), adz = diff(a.z, d.z);
    const Expansion bdx = diff(b.x, d.x), bdy = diff(b.y, d.y), bdz = diff(b.z, d.z);
    const Expansion cdx = diff(c.x, d.x), cdy = diff(c.y, d.y), cdz = diff(c.z, d.z);
    Expansion       r = mulExpansions(adz, crossTerm(bdx, cdy, bdy, cdx));
    r = addExpansions(std::move(r), mulExpansions(bdz, crossTerm(cdx, ady, cdy, adx)));
    r = addExpansions(std::move(r), mulExpansions(cdz, crossTerm(adx, bdy, ady, bdx)));
    return r.back();
}


double orient3d(Vec3d a, Vec3d b, Vec3d c, Vec3d d) noexcept {
    const double adx = a.x - d.x, bdx = b.x - d.x, cdx = c.x - d.x;
    const double ady = a.y - d.y, bdy = b.y - d.y, cdy = c.y - d.y;
    const double adz = a.z - d.z, bdz = b.z - d.z, cdz = c.z - d.z;

    const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy;
    const double cdxady = cdx * ady, adxcdy = adx * cdy;
    const double adxbdy = adx * bdy, bdxady = bdx * ady;
    const double det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
    if (!std::isfinite(det)) {
        return det;
    }

    const double permanent = (std::abs(bdxcdy) + std::abs(cdxbdy)) * std::abs(adz) +
                             (std::abs(cdxady) + std::abs(adxcdy)) * std::abs(bdz) +
                             (std::abs(adxbdy) + std::abs(bdxady)) * std::abs(cdz);
    const double errBound = kO3dErrBoundA * permanent;
    if (det > errBound || -det > errBound) {
        return det;
    }
    return orient3dExact(a, b, c, d);
}


double orient3d(Vec3f a, Vec3f b, Vec3f c, Vec3f d) noexcept {
    return orient3d(Vec3d { a.x, a.y, a.z }, Vec3d { b.x, b.y, b.z }, Vec3d { c.x, c.y, c.z },
                    Vec3d { d.x, d.y, d.z });
}


/// Exact incircle from the exact differences to d.
static double incircleExact(Vec2d a, Vec2d b, Vec2d c, Vec2d d) {
    const Expansion adx = diff(a.x, d.x), ady = diff(a.y, d.y);
    const Expansion bdx = diff(b.x, d.x), bdy = diff(b.y, d.y);
    const Expansion cdx = diff(c.x, d.x), cdy = diff(c.y, d.y);
    const Expansion aLift = addExpansions(mulExpansions(adx, adx), mulExpansions(ady, ady));
    const Expansion bLift = addExpansions(mulExpansions(bdx, bdx), mulExpansions(bdy, bdy));
    const Expansion cLift = addExpansions(mulExpansions(cdx, cdx), mulExpansions(cdy, cdy));
    Expansion       r = mulExpansions(aLift, crossTerm(bdx, cdy, cdx, bdy));
    r = addExpansions(std::move(r), mulExpansions(bLift, crossTerm(cdx, ady, adx, cdy)));
    r = addExpansions(std::move(r), mulExpansions(cLift, crossTerm(adx, bdy, bdx, ady)));
    return r.back();
}


double incircle(Vec2d a, Vec2d b, Vec2d c, Vec2d d) noexcept {
    const double adx = a.x - d.x, bdx = b.x - d.x, cdx = c.x - d.x;
    const double ady = a.y - d.y, bdy = b.y - d.y, cdy = c.y - d.y;

    const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy;
    const double cdxady = cdx * ady, adxcdy = adx * cdy;
    const double adxbdy = adx * bdy, bdxady = bdx * ady;
    const double aLift = adx * adx + ady * ady;
    const double bLift = bdx * bdx + bdy * bdy;
    const double cLift = cdx * cdx + cdy * cdy;
    const double det = aLift * (bdxcdy - cdxbdy) + bLift * (cdxady - adxcdy) + cLift * (adxbdy - bdxady);
    if (!std::isfinite(det)) {
        return det;
    }

    const double permanent = (std::abs(bdxcdy) + std::abs(cdxbdy)) * aLift +
                             (std::abs(cdxady) + std::abs(adxcdy)) * bLift +
                             (std::abs(adxbdy) + std::abs(bdxady)) * cLift;
    const double errBound = kIccErrBoundA * permanent;
    if (det > errBound || -det > errBound) {
        return det;
    }
    return incircleExact(a, b, c, d);
}


double incircle(Vec2f a, Vec2f b, Vec2f c, Vec2f d) noexcept {
    return incircle(Vec2d { a.x, a.y }, Vec2d { b.x, b.y }, Vec2d { c.x, c.y }, Vec2d { d.x, d.y });
}


bool segmentsCross(Vec2d a, Vec2d b, Vec2d c, Vec2d d) noexcept {
    const double d1 = orient2d(c, d, a);
    const double d2 = orient2d(c, d, b);
    const double d3 = orient2d(a, b, c);
    const double d4 = orient2d(a, b, d);
    return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
}


bool segmentsCross(Vec2f a, Vec2f b, Vec2f c, Vec2f d) noexcept {
    return segmentsCross(Vec2d { a.x, a.y }, Vec2d { b.x, b.y }, Vec2d { c.x, c.y }, Vec2d { d.x, d.y });
}


double signedArea2(std::span<const Vec2d> p) noexcept {
    double s = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const Vec2d &a = p[i];
        const Vec2d &b = p[(i + 1) % p.size()];
        s += a.x * b.y - b.x * a.y;
    }
    return s;
}


double signedArea2(std::span<const Vec2f> p) noexcept {
    double s = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const Vec2f &a = p[i];
        const Vec2f &b = p[(i + 1) % p.size()];
        s += static_cast<double>(a.x) * b.y - static_cast<double>(b.x) * a.y;
    }
    return s;
}

} // namespace decimate::math
