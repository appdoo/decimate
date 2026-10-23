// This file is part of the decimate library; see decimate.h for version/license details.
/// @file normal_bake.cpp
/// @brief Raytraced normal bake - port of normalmap::mapnormals (legacy2003/split/texturesample.cpp).
///
/// Two phases: (1) a serial raster pass records which low face owns each texel (the LAST covering
/// face wins, the 2003 overwrite order, legacy2003/split/texturesample.cpp:172-173,189,365-368);
/// (2) every owned texel is shaded independently (ray against the high mesh through the injected
/// Tracer, normal blend, range-compress) on worker threads. A texel's value depends only on
/// (owner, x, y), so the output is byte-identical for any thread count.

#include <decimate/bake/normal_bake.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <thread>
#include <unordered_map>

#include "bake_detail.h"

namespace decimate::bake {

namespace {

/// Barycentric coordinates (u, v, w) of a point in a triangle; all in [0,1] when inside.
struct Bc3 {
    float u;
    float v;
    float w;
};


/// Coordinate in texel units, [0, texSize]: texel x spans [x, x+1] and its centre is x + 0.5.
struct UV2 {
    float u;
    float v;
};


/// Exclusive |t| bound of the nearest-hit query (the 2003 loop starts at mint = 100000.0,
/// legacy2003/split/texturesample.cpp:84).
constexpr float kMaxRayT = 100000.0f;


/// Per-low-face raster setup, computed once and shared by both bake phases so the two phases
/// evaluate the SAME expressions on the SAME floats (no reordering across the phase split).
/// The 2003 counterpart is the prologue of mapfacert (legacy2003/split/texturesample.cpp:139-170).
struct FaceSetup {
    UV2         t1 {};
    UV2         t2 {};
    UV2         t3 {};
    float       bary = 0.0f; ///< barycentric divisor (t2-t1) x (t3-t1), fixed for the triangle
    int         minu = 0;    ///< integer UV bounding box clamped to the texture
    int         maxu = 0;
    int         minv = 0;
    int         maxv = 0;
    math::Vec3f n1; ///< low vertex normals (ray directions)
    math::Vec3f n2;
    math::Vec3f n3;
    math::Vec3f pa; ///< low vertex positions (ray origins)
    math::Vec3f pb;
    math::Vec3f pc;
    math::Vec4f tg1; ///< tangent-space option: MikkTSpace corner tangents (xyz = T, w = sign)
    math::Vec4f tg2;
    math::Vec4f tg3;
    bool        tValid = false; ///< all three corner tangents exist (else the texel stores the flat normal)
};


/// Everything phase 2 reads. `shadeRow` writes only row y of the map, so rows are independent.
/// Member functions are defined after the file-local helpers they call.
struct ShadeJob {
    const std::vector<FaceSetup>   &setup;
    const std::vector<int>         &owner;
    const detail::Tracer           &tracer;
    const decimate::MeshData       &high;
    const std::vector<math::Vec3f> &highNormals;
    bool                            tangentSpace;

    /// Phase 2: ray + blend + range-compress for one owned texel. Pure function of its arguments.
    /// 2003: mapfacert inner body (legacy2003/split/texturesample.cpp:177-190) and rtraceHighOrderNormal
    /// (:82-122; a miss yields the ray direction, :106-110).
    [[nodiscard]] Rgb shadeTexel(const FaceSetup &f, int x, int y) const;

    /// Shade every owned texel of row `y`; texels without an owner keep the background.
    void shadeRow(int y, NormalMap &nm) const;
};


/// A seam edge in atlas space: two 2D endpoints (in texel coordinates) for
/// each side (A = owning half-edge, B = sym half-edge).
///
/// t in [0,1] parameterizes along the edge from p0 to p1.
/// The sym edge runs p1->p0, so t along A maps to (1-t) along B.
struct SeamEdge {
    math::Vec2f p0;  ///< Atlas start (texel coords), origin vertex of he.
    math::Vec2f p1;  ///< Atlas end   (texel coords), end vertex of he.
    math::Vec2f inN; ///< Inward unit normal (into chart A's interior).

    // Sym edge (B side) — same endpoints but REVERSED order.
    math::Vec2f bp0;  ///< Atlas start for B (= A's p1 direction at t=0 on B).
    math::Vec2f bp1;  ///< Atlas end   for B.
    math::Vec2f bInN; ///< Inward unit normal for B.

    bool hasRefMap = false; ///< true = chart seam (neighbour exists); false = open mesh.
};


/// Production tracer: the exact BVH over the high mesh.
class BvhTracer final : public detail::Tracer {
public:
    explicit BvhTracer(const decimate::MeshData &high)
        : bvh_(high) {
    }

    [[nodiscard]] rtrace::TriangleBvh::Hit nearest(math::Vec3f orig,
                                                   math::Vec3f dir,
                                                   float       maxAbsT) const override {
        return bvh_.nearest(orig, dir, maxAbsT);
    }

private:
    rtrace::TriangleBvh bvh_;
};

} // namespace


/// Barycentric coordinates of (pu, pv) in triangle (t1, t2, t3); `bdiv` is the precomputed
/// divisor (t2-t1) x (t3-t1). 2003: legacy2003/split/texturesample.cpp:55-63.
static Bc3 barycentricUV(float pu, float pv, UV2 t1, UV2 t2, UV2 t3, float bdiv) {
    Bc3 bc;
    bc.u = ((t2.u - pu) * (t3.v - pv) - (t3.u - pu) * (t2.v - pv)) / bdiv;
    bc.v = ((t3.u - pu) * (t1.v - pv) - (t1.u - pu) * (t3.v - pv)) / bdiv;
    bc.w = ((t1.u - pu) * (t2.v - pv) - (t2.u - pu) * (t1.v - pv)) / bdiv;
    return bc;
}


/// All barycentric coordinates in [0,1]? 2003: inside(), legacy2003/split/texturesample.cpp:30.
/// The double-negation form keeps the NaN behavior: NaN comparisons are false, so NaN coordinates
/// pass as "inside". The simplified form would change NaN semantics.
static bool insideTriangle(float u, float v, float w) {
    // NOLINTNEXTLINE(readability-simplify-boolean-expr)
    return !((u < 0.0f || v < 0.0f || w < 0.0f) || (u > 1.0f || v > 1.0f || w > 1.0f));
}


/// Interpolate a Vec3 with barycentric weights. 2003: legacy2003/split/texturesample.cpp:128.
static math::Vec3f interpolateVec3(const Bc3         &bc,
                                   const math::Vec3f &v1,
                                   const math::Vec3f &v2,
                                   const math::Vec3f &v3) {
    return { bc.u * v1.x + bc.v * v2.x + bc.w * v3.x,
             bc.u * v1.y + bc.v * v2.y + bc.w * v3.y,
             bc.u * v1.z + bc.v * v2.z + bc.w * v3.z };
}


/// Round the 8-bit code of `c` in [-1, 1]: round((c + 1) * 127). The decode is `code / 127 - 1`, so
/// code 127 is exactly 0 and 254 is exactly +1; every code is the centre of its bucket, the error
/// is at most 0.5 / 127. 2003 truncated here (legacy2003/split/texturesample.cpp:74-76), which
/// biases every component by half a code, and by a different amount in every tangent frame.
static std::uint8_t encodeComponent(float c) {
    return static_cast<std::uint8_t>(std::floor((c + 1.0f) * 127.0f + 0.5f));
}


namespace detail {

Rgb rangeCompress(math::Vec3f nrm) {
    nrm.normalize();
    return { encodeComponent(nrm.x), encodeComponent(nrm.y), encodeComponent(nrm.z) };
}

} // namespace detail


/// Raster setup of low face `fi`: UVs scaled to texels, clamped bounding box, ray origins/directions
/// and (in tangent space) the corner tangents. Run once per face so both bake phases see identical floats.
static FaceSetup makeFaceSetup(int fi, const decimate::MeshData &low,
                               const std::vector<math::Vec2f> &cornerUV,
                               const std::vector<math::Vec3f> &lowNormals,
                               const std::vector<math::Vec4f> &tangents, int texSize) {
    FaceSetup f;

    // positions is std::vector<Vec3f>; Vec3f is tightly packed, so the flat coord buffer indexes
    // as lp[3*v+k] (the offsets below are the hoisted 3*v bases).
    const float *lp = reinterpret_cast<const float *>(low.positions.data());
    const int   *li = low.indices.data();
    const float  ts = static_cast<float>(texSize);

    // Corner UV in [0,1] -> texels [0, texSize]: texel x is centred at uv = (x + 0.5) / texSize, where
    // GL and the seam erasure read it (the 2003 bake scaled by texSize - 1).
    f.t1 = { cornerUV[3 * fi + 0].x * ts, cornerUV[3 * fi + 0].y * ts };
    f.t2 = { cornerUV[3 * fi + 1].x * ts, cornerUV[3 * fi + 1].y * ts };
    f.t3 = { cornerUV[3 * fi + 2].x * ts, cornerUV[3 * fi + 2].y * ts };

    f.bary = (f.t2.u - f.t1.u) * (f.t3.v - f.t1.v) - (f.t3.u - f.t1.u) * (f.t2.v - f.t1.v);

    // A zero or non-finite UV area has no barycentrics (0/0 = NaN, which insideTriangle accepts) and
    // would claim texels of a valid neighbour: leave it with an empty box so the owner pass skips it.
    if (f.bary == 0.0f || !std::isfinite(f.bary)) {
        f.maxu = -1;
        f.maxv = -1;
        return f;
    }

    // Integer bounding box of the UV triangle, clamped to the texture.
    // NOLINTBEGIN(bugprone-incorrect-roundings): the 2003 bounding-box expression, kept verbatim
    const int minu = static_cast<int>(std::min({ f.t1.u, f.t2.u, f.t3.u }));
    const int maxu = static_cast<int>(std::max({ f.t1.u, f.t2.u, f.t3.u }) + 0.5f);
    const int minv = static_cast<int>(std::min({ f.t1.v, f.t2.v, f.t3.v }));
    const int maxv = static_cast<int>(std::max({ f.t1.v, f.t2.v, f.t3.v }) + 0.5f);
    // NOLINTEND(bugprone-incorrect-roundings)
    f.minu = std::min(std::max(minu, 0), texSize - 1);
    f.maxu = std::min(std::max(maxu, 0), texSize - 1);
    f.minv = std::min(std::max(minv, 0), texSize - 1);
    f.maxv = std::min(std::max(maxv, 0), texSize - 1);

    // Vertex indices and position base offsets (index-hoist: compute once, use three times)
    const int va = li[3 * fi + 0];
    const int vb = li[3 * fi + 1];
    const int vc = li[3 * fi + 2];
    const int oa = va * 3;
    const int ob = vb * 3;
    const int oc = vc * 3;
    f.n1 = lowNormals[va];
    f.n2 = lowNormals[vb];
    f.n3 = lowNormals[vc];
    f.pa = { lp[oa], lp[oa + 1], lp[oa + 2] };
    f.pb = { lp[ob], lp[ob + 1], lp[ob + 2] };
    f.pc = { lp[oc], lp[oc + 1], lp[oc + 2] };

    // Tangent-space option: stored normals become relative to the MikkTSpace frame -> near-flat
    // (0,0,1), so chart values barely differ across a seam. A zero tangent marks a face with no basis.
    if (!tangents.empty()) {
        f.tg1 = tangents[3 * fi + 0];
        f.tg2 = tangents[3 * fi + 1];
        f.tg3 = tangents[3 * fi + 2];
        f.tValid = f.tg1.x * f.tg1.x + f.tg1.y * f.tg1.y + f.tg1.z * f.tg1.z > 0.0f &&
                   f.tg2.x * f.tg2.x + f.tg2.y * f.tg2.y + f.tg2.z * f.tg2.z > 0.0f &&
                   f.tg3.x * f.tg3.x + f.tg3.y * f.tg3.y + f.tg3.z * f.tg3.z > 0.0f;
    }
    return f;
}


/// Barycentric coords of texel centre (x, y) in face `f`.
///
/// Deliberate deviation from 2003: `mapfacert` samples at the texel corner `(float)x, (float)y`
/// (`legacy2003/split/texturesample.cpp:177`); this samples at the centre (`+ 0.5f`), which is the
/// position a bilinear GL lookup reads for texel (x, y).
static Bc3 texelBary(const FaceSetup &f, int x, int y) {
    return barycentricUV(static_cast<float>(x) + 0.5f, static_cast<float>(y) + 0.5f, f.t1, f.t2, f.t3,
                         f.bary);
}


/// Claim the texels of column `x` inside face `fi`'s bounding box whose centre is inside the triangle.
static void claimColumn(const FaceSetup &f, int fi, int x, int texSize, std::vector<int> &owner) {
    for (int y = f.minv; y <= f.maxv; y++) {
        const Bc3 bc = texelBary(f, x, y);
        if (insideTriangle(bc.u, bc.v, bc.w)) {
            owner[static_cast<std::size_t>(y) * texSize + x] = fi;
        }
    }
}


/// Phase 1: mark every texel whose centre lies inside face `fi` as owned by it. Faces are claimed
/// in ascending order, so the last covering face wins - the 2003 overwrite order.
static void claimTexels(const FaceSetup &f, int fi, int texSize, std::vector<int> &owner) {
    for (int x = f.minu; x <= f.maxu; x++) {
        claimColumn(f, fi, x, texSize, owner);
    }
}


/// Express the unit normal `want` in the MikkTSpace frame of the texel, so that the engines' decode
/// `normalize(x*T + y*B + z*N)` (T, N interpolated unnormalised, B = sign * cross(N, T)) returns `want`.
/// The frame is not orthonormal, so this is the inverse of [T B N] (Cramer), not three dot products.
/// A face without a tangent, or a collapsed frame (|det| tiny: T parallel to N), stores the flat normal
/// (0,0,1), which decodes to the interpolated vertex normal: no detail, but no garbage either.
static math::Vec3f toTangentSpace(const FaceSetup &f, const Bc3 &bc, const math::Vec3f &n,
                                  const math::Vec3f &want) {
    const math::Vec3f flat { 0.0f, 0.0f, 1.0f };
    if (!f.tValid) {
        return flat;
    }

    const math::Vec3f t = interpolateVec3(bc, { f.tg1.x, f.tg1.y, f.tg1.z }, { f.tg2.x, f.tg2.y, f.tg2.z },
                                          { f.tg3.x, f.tg3.y, f.tg3.z });
    const float       sign = bc.u * f.tg1.w + bc.v * f.tg2.w + bc.w * f.tg3.w < 0.0f ? -1.0f : 1.0f;
    const math::Vec3f b = n.cross(t) * sign;

    const math::Vec3f bn = b.cross(n);
    const float       det = t.dot(bn);
    if (!(std::fabs(det) > 1e-6f * t.length() * b.length() * n.length())) {
        return flat;
    }
    return { want.dot(bn) / det, want.dot(n.cross(t)) / det, want.dot(t.cross(b)) / det };
}


Rgb ShadeJob::shadeTexel(const FaceSetup &f, int x, int y) const {
    const int *hi = high.indices.data();
    const Bc3  bc = texelBary(f, x, y);

    // Interpolated low-poly normal (ray direction) and surface point
    const math::Vec3f n = interpolateVec3(bc, f.n1, f.n2, f.n3);
    const math::Vec3f p = interpolateVec3(bc, f.pa, f.pb, f.pc);

    // Cast along the low normal; on a hit blend the high vertex normals at the hit point.
    const rtrace::TriangleBvh::Hit hit = tracer.nearest(p, n, kMaxRayT);
    math::Vec3f                    snrm = n; // miss fallback: the low-poly normal
    if (hit.face >= 0) {
        const float u = hit.u;
        const float v = hit.v;
        snrm = highNormals[hi[3 * hit.face + 0]] * (1.0f - (u + v)) + highNormals[hi[3 * hit.face + 1]] * u +
               highNormals[hi[3 * hit.face + 2]] * v;
    }
    snrm.normalize();

    if (tangentSpace) {
        snrm = toTangentSpace(f, bc, n, snrm);
    }
    return detail::rangeCompress(snrm);
}


void ShadeJob::shadeRow(int y, NormalMap &nm) const {
    const std::size_t base = static_cast<std::size_t>(y) * nm.width;
    for (int x = 0; x < nm.width; x++) {
        const int fi = owner[base + x];
        if (fi < 0) {
            continue; // no face covers this texel: background stays
        }
        nm.pixels[base + x] = shadeTexel(setup[fi], x, y);
    }
}


/// Worker count for `threads` (<= 0: all cores) over `rows` rows. 1 when threads are compiled out.
static int workerCount([[maybe_unused]] int threads, [[maybe_unused]] int rows) {
#ifdef DECIMATE_NO_THREADS
    return 1;
#else
    // Never more workers than rows: a worker without a row is pure thread-start cost.
    const int want = threads > 0 ? threads : static_cast<int>(std::max(1U, std::thread::hardware_concurrency()));
    return std::max(1, std::min(want, rows));
#endif
}


/// Run `body(row)` for every row in [0, rows) on `workers` workers (the caller's thread when 1).
/// Rows are handed out by an atomic counter; `body` must write only data owned by its row.
template <typename Body>
static void forEachRow(int rows, int workers, const Body &body) {
    if (workers <= 1) {
        for (int r = 0; r < rows; r++) {
            body(r);
        }
        return;
    }
#ifndef DECIMATE_NO_THREADS
    // Workers pull the next unclaimed row, so uneven row costs balance themselves; the caller works too.
    std::atomic<int> next { 0 };
    const auto       drain = [&] {
        for (int r = next.fetch_add(1); r < rows; r = next.fetch_add(1)) {
            body(r);
        }
    };
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(workers - 1));
    for (int w = 1; w < workers; w++) {
        pool.emplace_back(drain);
    }
    drain();
    for (auto &t : pool) {
        t.join();
    }
#endif
}


/// One dilation pass: read from src, write into dst.
/// Background texels receive the integer average of non-background 8-neighbours.
/// Filled texels and isolated background texels copy through unchanged.
/// 2003: the inner loop of dilate(), legacy2003/split/texturesample.cpp:222-275 (double buffer, 8 neighbours, % wrap).
// NOLINTNEXTLINE(readability-function-cognitive-complexity,readability-function-size)
static void dilateOnce(int texsize, const std::vector<Rgb> &src, std::vector<Rgb> &dst) {
    const Rgb background { 128, 128, 128 };

    for (int y = 0; y < texsize; y++) {
        for (int x = 0; x < texsize; x++) {
            const int idx = y * texsize + x;

            const auto &cur = src[idx];
            auto       &out = dst[idx];

            // Filled texel: copy through unchanged.
            if (cur != background) {
                out = cur;
                continue;
            }

            // Background texel: integer average of non-background 8-neighbours. The % wrap and the
            // negative-index fixup are kept verbatim from 2003 (legacy2003/split/texturesample.cpp:241-244).
            int accu[3] = { 0, 0, 0 }; // NOLINT(modernize-avoid-c-arrays)
            int acnt = 0;

            for (int i = x - 1; i < x + 2; i++) {
                for (int j = y - 1; j < y + 2; j++) {
                    if (i == x && j == y) {
                        continue;
                    }

                    int xn = i % texsize;
                    int yn = j % texsize;
                    if (xn < 0) {
                        xn = texsize - 1;
                    }
                    if (yn < 0) {
                        yn = texsize - 1;
                    }

                    const auto &nb = src[yn * texsize + xn];
                    if (nb != background) {
                        accu[0] += nb.r;
                        accu[1] += nb.g;
                        accu[2] += nb.b;
                        acnt++;
                    }
                }
            }

            // Average of the filled neighbours, or stay background if there are none.
            if (acnt > 0) {
                out.r = static_cast<std::uint8_t>(accu[0] / acnt);
                out.g = static_cast<std::uint8_t>(accu[1] / acnt);
                out.b = static_cast<std::uint8_t>(accu[2] / acnt);
            } else {
                out = cur; // stay background
            }
        }
    }
}


/// `passes`-pass dilation: fill background texels from non-background neighbours.
/// Double-buffer, 8-neighbour average, % wrap; each pass grows the filled region by 1 texel.
/// The pass count is the caller's choice (see bake()'s `dilatePasses`); the 2003 dilate() runs a
/// fixed 6 (legacy2003/split/texturesample.cpp:222).
static void dilate(NormalMap &nm, int passes) {
    const int texsize = nm.width;

    // Ping-pong between the map and a scratch buffer so a pass never reads its own output.
    std::vector<Rgb>  scratch(texsize * texsize);
    std::vector<Rgb> *img = &nm.pixels;
    std::vector<Rgb> *scr = &scratch;

    for (int pass = 0; pass < passes; pass++) {
        dilateOnce(texsize, *img, *scr);
        std::swap(img, scr);
    }

    // Odd pass count leaves the result in `scratch`; copy it back.
    if (img != &nm.pixels) {
        nm.pixels = *img;
    }
}


namespace detail {

// NOLINTNEXTLINE(readability-function-size): the public bake() argument list plus the tracer
NormalMap bakeWith(const Tracer &tracer, const decimate::MeshData &low,
                   const std::vector<math::Vec2f> &lowCornerUV,
                   const std::vector<math::Vec3f> &lowNormals, const decimate::MeshData &high,
                   const std::vector<math::Vec3f> &highNormals, int texSize, int dilatePasses,
                   const std::vector<math::Vec4f> &tangents, int threads) {
    // Everything starts as the background sentinel; only owned texels are overwritten.
    NormalMap nm;
    nm.width = texSize;
    nm.height = texSize;
    nm.pixels.assign(static_cast<std::size_t>(texSize) * texSize, Rgb { 128, 128, 128 });

    const int nlf = static_cast<int>(low.indices.size()) / 3;

    // Phase 1 (serial, cheap): who owns each texel. Last covering face wins.
    std::vector<FaceSetup> setup;
    setup.reserve(static_cast<std::size_t>(nlf));
    std::vector<int> owner(static_cast<std::size_t>(texSize) * texSize, -1);
    for (int fi = 0; fi < nlf; fi++) {
        setup.push_back(makeFaceSetup(fi, low, lowCornerUV, lowNormals, tangents, texSize));
        claimTexels(setup.back(), fi, texSize, owner);
    }

    // Phase 2 (parallel): shade each owned texel; each row is written by exactly one worker.
    const ShadeJob job { setup, owner, tracer, high, highNormals, !tangents.empty() };
    forEachRow(texSize, workerCount(threads, texSize), [&](int y) { job.shadeRow(y, nm); });

    // Dilation is opt-in (dilatePasses == 0 leaves crisp charts, no normals outside them). The 2003
    // mapnormals always dilates, 6 passes (legacy2003/split/texturesample.cpp:382 calling dilate(),
    // :214, loop at :222); the default 0 is deliberate: the seam fill (bakeSeamFill) replaces dilate
    // at chart edges by baking the real adjacent triangle into the gutter. The packer must reserve
    // matching clearance for the grown gutter (see packAtlas).
    if (dilatePasses > 0) {
        dilate(nm, dilatePasses);
    }
    return nm;
}

} // namespace detail


// NOLINTNEXTLINE(readability-function-size): public signature, the 2003 argument list plus `threads`
NormalMap bake(const decimate::MeshData       &low,
               const std::vector<math::Vec2f> &lowCornerUV,
               const std::vector<math::Vec3f> &lowNormals,
               const decimate::MeshData       &high,
               const std::vector<math::Vec3f> &highNormals,
               int                             texSize,
               int                             dilatePasses,
               const std::vector<math::Vec4f> &tangents,
               int                             threads) {
    // The BVH is built once per bake and shared read-only by all workers.
    const BvhTracer tracer(high);
    return detail::bakeWith(tracer, low, lowCornerUV, lowNormals, high, highNormals, texSize, dilatePasses,
                            tangents, threads);
}


/// Project point `g` onto the line segment [p0, p1].
/// Returns (t, d²) where t = clamped param in [0,1] and d² = squared distance.
static std::pair<float, float> projectOnSegment(math::Vec2f g,
                                                math::Vec2f p0,
                                                math::Vec2f p1) {
    float ex = p1.x - p0.x;
    float ey = p1.y - p0.y;
    float len2 = ex * ex + ey * ey;

    // Degenerate segment: the nearest point is p0.
    if (len2 < 1e-12f) {
        float dx = g.x - p0.x;
        float dy = g.y - p0.y;
        return { 0.0f, dx * dx + dy * dy };
    }

    float t = ((g.x - p0.x) * ex + (g.y - p0.y) * ey) / len2;
    t = std::max(0.0f, std::min(1.0f, t));
    float fx = p0.x + t * ex - g.x;
    float fy = p0.y + t * ey - g.y;
    return { t, fx * fx + fy * fy };
}


/// Bilinear sample the atlas at floating-point texel coordinates (u, v), skipping background taps
/// so a chart edge does not bleed the sentinel colour. (u, v) are in texel units with texel x centred
/// at x + 0.5, like every coordinate of the bake; they clamp to the outermost centres. Returns the
/// background sentinel when no tap holds data.
static Rgb bilinearSample(const std::vector<Rgb> &pixels, int texSize, float u, float v) {
    const Rgb kBg { 128, 128, 128 };

    // Centre convention: the tap at index i sits at i + 0.5. Clamp to the outermost centres.
    u = std::max(0.0f, std::min(static_cast<float>(texSize - 1), u - 0.5f));
    v = std::max(0.0f, std::min(static_cast<float>(texSize - 1), v - 0.5f));

    int x0 = static_cast<int>(u);
    int y0 = static_cast<int>(v);
    int x1 = std::min(x0 + 1, texSize - 1);
    int y1 = std::min(y0 + 1, texSize - 1);

    float fu = u - static_cast<float>(x0);
    float fv = v - static_cast<float>(y0);

    const Rgb &c00 = pixels[y0 * texSize + x0];
    const Rgb &c10 = pixels[y0 * texSize + x1];
    const Rgb &c01 = pixels[y1 * texSize + x0];
    const Rgb &c11 = pixels[y1 * texSize + x1];

    // Nothing to blend: stay background.
    if (c00 == kBg && c10 == kBg && c01 == kBg && c11 == kBg) {
        return kBg;
    }

    // Background taps are transparent: weight only the filled ones and renormalise below.
    float r = 0.0f;
    float g = 0.0f;
    float b = 0.0f;
    float w = 0.0f;

    auto addTap = [&](const Rgb &c, float weight) {
        if (c == kBg) {
            return;
        }
        r += weight * static_cast<float>(c.r);
        g += weight * static_cast<float>(c.g);
        b += weight * static_cast<float>(c.b);
        w += weight;
    };

    addTap(c00, (1.0f - fu) * (1.0f - fv));
    addTap(c10, fu * (1.0f - fv));
    addTap(c01, (1.0f - fu) * fv);
    addTap(c11, fu * fv);

    // Weight vanished (e.g. only a zero-weight tap is filled): treat as no data.
    if (w < 1e-6f) {
        return kBg;
    }
    return { static_cast<std::uint8_t>(r / w + 0.5f),
             static_cast<std::uint8_t>(g / w + 0.5f),
             static_cast<std::uint8_t>(b / w + 0.5f) };
}


/// Edge geometry of half-edge `he` (= 3*f + k) in texel-scale coordinates.
/// Returns {p0 = origin UV, p1 = end UV, opposite vertex UV} scaled to texels [0, texSize].
/// Origin vertex UV, writeFace convention:
///   k=0 (ac): origin=a (cornerUV[3f+0]), end=c (cornerUV[3f+2])
///   k=1 (cb): origin=c (cornerUV[3f+2]), end=b (cornerUV[3f+1])
///   k=2 (ba): origin=b (cornerUV[3f+1]), end=a (cornerUV[3f+0])
/// Third vertex (opposite, inside the face):
///   k=0: opposite = b (cornerUV[3f+1])
///   k=1: opposite = a (cornerUV[3f+0])
///   k=2: opposite = c (cornerUV[3f+2])
static std::tuple<math::Vec2f, math::Vec2f, math::Vec2f>
heUV(decimate::HeIndex he, const std::vector<math::Vec2f> &cornerUV, float ts) {
    int f = static_cast<int>(he / 3u);
    int k = static_cast<int>(he % 3u);

    // Corner slot (0..2) of the origin, end and opposite vertex for each k.
    //   k=0: origin=cornerUV[3f+0], end=cornerUV[3f+2], opp=cornerUV[3f+1]
    //   k=1: origin=cornerUV[3f+2], end=cornerUV[3f+1], opp=cornerUV[3f+0]
    //   k=2: origin=cornerUV[3f+1], end=cornerUV[3f+0], opp=cornerUV[3f+2]
    int origIdx = (k == 0) ? 0 : (k == 1) ? 2
                                          : 1;
    int endIdx = (k == 0) ? 2 : (k == 1) ? 1
                                         : 0;
    int oppIdx = (k == 0) ? 1 : (k == 1) ? 0
                                         : 2;

    // Corner UVs are in [0,1]; the seam math works in texels.
    auto scale = [&](math::Vec2f uv) {
        return math::Vec2f { uv.x * ts, uv.y * ts };
    };

    return { scale(cornerUV[3 * f + origIdx]),
             scale(cornerUV[3 * f + endIdx]),
             scale(cornerUV[3 * f + oppIdx]) };
}


/// 2D inward unit normal of edge [p0, p1]: the perpendicular pointing toward `opp`.
/// Zero for a degenerate edge.
static math::Vec2f inwardNormal2D(math::Vec2f p0,
                                  math::Vec2f p1,
                                  math::Vec2f opp) {
    float ex = p1.x - p0.x;
    float ey = p1.y - p0.y;
    float len = std::sqrt(ex * ex + ey * ey);

    if (len < 1e-12f) {
        return { 0.0f, 0.0f };
    }

    // Of the two perpendiculars (-ey, ex) and (ey, -ex), keep the one on opp's side.
    float nx = -ey / len;
    float ny = ex / len;
    float dot = nx * (opp.x - p0.x) + ny * (opp.y - p0.y);
    if (dot < 0.0f) {
        nx = -nx;
        ny = -ny;
    }
    return { nx, ny };
}


/// Seam-edge catalog: every half-edge on a chart boundary, i.e. with no sym (open mesh boundary)
/// or with a sym in a DIFFERENT chart, keyed by half-edge index. Chart seams carry the sym side's
/// atlas edge too, which is what the gutter fill samples from.
static std::unordered_map<decimate::HeIndex, SeamEdge>
buildSeamEdges(const decimate::HalfEdgeMesh   &hem,
               const std::vector<math::Vec2f> &cornerUV,
               const std::vector<int>         &faceToChart,
               float                           ts) {
    std::unordered_map<decimate::HeIndex, SeamEdge> seams;
    const auto                                     &hes = hem.halfedges();

    for (decimate::HeIndex i = 0; i < static_cast<decimate::HeIndex>(hes.size()); ++i) {
        const auto       &h = hes[i];
        decimate::HeIndex s = h.sym;

        // Only half-edges on a CHART BOUNDARY matter (open-mesh boundary or a different chart).
        bool isOpenBoundary = (s == decimate::kNoHalfedge);
        bool isChartSeam = false;

        if (!isOpenBoundary && s < static_cast<decimate::HeIndex>(hes.size())) {
            int chartA = faceToChart[h.face];
            int chartB = faceToChart[hes[s].face];
            isChartSeam = (chartA != chartB);
        }

        if (!isOpenBoundary && !isChartSeam) {
            continue; // interior edge: both sides are in the same chart
        }

        // The A side: this half-edge's own atlas edge.
        auto [p0, p1, opp] = heUV(i, cornerUV, ts);
        math::Vec2f inN = inwardNormal2D(p0, p1, opp);

        SeamEdge se;
        se.p0 = p0;
        se.p1 = p1;
        se.inN = inN;
        se.hasRefMap = isChartSeam;

        // The B side, only across a chart seam: the sym edge runs p1_b -> p0_b in B's atlas
        // (B's p0 matches A's p1, the sym reversal).
        if (isChartSeam) {
            auto [bp0, bp1, bopp] = heUV(s, cornerUV, ts);
            se.bp0 = bp0;
            se.bp1 = bp1;
            se.bInN = inwardNormal2D(bp0, bp1, bopp);
        }

        seams.emplace(i, se);
    }
    return seams;
}


NormalMap bakeSeamFill(const decimate::MeshData       &low,
                       const decimate::HalfEdgeMesh   &hem,
                       const std::vector<math::Vec2f> &lowCornerUV,
                       const std::vector<int>         &faceToChart,
                       const std::vector<math::Vec3f> &lowNormals,
                       const decimate::MeshData       &high,
                       const std::vector<math::Vec3f> &highNormals,
                       int                             texSize,
                       int                             dilatePasses,
                       const std::vector<math::Vec4f> &tangents,
                       int                             threads) {
    // Pass 1: the interior bake, without dilation.
    NormalMap nm = bake(low, lowCornerUV, lowNormals, high, highNormals, texSize, 0, tangents, threads);

    // Seam catalog. A mesh edge is a CHART SEAM iff sym(he) != kNoHalfedge AND the two faces' charts
    // differ (same 3D edge, different atlas locations). Idea: Ray et al. 2010 "Invisible Seams"
    // indirection map, NOT the MS patent US10,839,590 polygon-projection mechanism.
    const float ts = static_cast<float>(texSize);
    auto        seams = buildSeamEdges(hem, lowCornerUV, faceToChart, ts);

    // No chart seams (single-chart mesh or trivial case): the plain blur dilate is all there is.
    if (seams.empty()) {
        if (dilatePasses > 0) {
            dilate(nm, dilatePasses);
        }
        return nm;
    }

    // Flatten the seam map to a vector for O(E) nearest-edge search per texel.
    std::vector<std::pair<decimate::HeIndex, SeamEdge>> seamVec(seams.begin(), seams.end());

    // Pass 2: seam-correct gutter fill. Pass 1 stays frozen (all samples read nm), results go to pass2.
    const Rgb        kBg { 128, 128, 128 };
    std::vector<Rgb> pass2 = nm.pixels;

    // The gutter is `dilatePasses` texels wide (same width as the blur dilate). Every background texel
    // inside it takes its nearest seam edge: a chart seam samples the neighbour chart, an open boundary
    // is left to the blur dilate. One O(texels x seam edges) pass - the simple correctness path.
    for (int y = 0; y < texSize; ++y) {
        for (int x = 0; x < texSize; ++x) {
            const int idx = y * texSize + x;
            if (nm.pixels[idx] != kBg) {
                continue; // already filled in pass 1
            }

            math::Vec2f g { static_cast<float>(x) + 0.5f, static_cast<float>(y) + 0.5f };

            // Nearest seam edge and the foot parameter on it.
            float           bestDist2 = std::numeric_limits<float>::max();
            float           bestT = 0.0f;
            bool            bestHasRef = false;
            const SeamEdge *bestEdge = nullptr;

            for (const auto &[he, se] : seamVec) {
                auto [t, d2] = projectOnSegment(g, se.p0, se.p1);
                if (d2 < bestDist2) {
                    bestDist2 = d2;
                    bestT = t;
                    bestHasRef = se.hasRefMap;
                    bestEdge = &se;
                }
            }

            if (bestEdge == nullptr) {
                continue; // no seam edge at all — leave background
            }

            const float gDist = std::sqrt(bestDist2); // perpendicular distance d

            // Bound the fill to the gutter band: only texels within `dilatePasses`
            // of a seam edge are filled. Without this the loop floods the ENTIRE
            // empty atlas (every background texel finds *some* nearest seam edge),
            // smearing far regions and producing seam-like artefacts near borders.
            if (gDist > static_cast<float>(dilatePasses)) {
                continue; // beyond the gutter — leave background
            }

            if (!bestHasRef) {
                // Open-mesh boundary: classic dilate will handle this texel.
                // Leave it as background for the blur pass below.
                continue;
            }

            // Chart seam: sample neighbour chart B at the sym point.
            // Invariants (concept-seamless-bake.md):
            //   1. Param reversal: t along e_a -> (1-t) along e_b.
            //   2. Perpendicular flip: d outside A = d INSIDE B (bInN points in).
            float tB = 1.0f - bestT;
            float srcU = bestEdge->bp0.x + tB * (bestEdge->bp1.x - bestEdge->bp0.x) + bestEdge->bInN.x * gDist;
            float srcV = bestEdge->bp0.y + tB * (bestEdge->bp1.y - bestEdge->bp0.y) + bestEdge->bInN.y * gDist;

            // A sample that lands in B's background (very thin chart) stays for the dilate fallback below.
            Rgb sampled = bilinearSample(nm.pixels, texSize, srcU, srcV);
            if (sampled != kBg) {
                pass2[idx] = sampled;
            }
        }
    }

    nm.pixels = std::move(pass2);

    // Blur-dilate for open-mesh boundaries and any remaining background.
    if (dilatePasses > 0) {
        dilate(nm, dilatePasses);
    }
    return nm;
}

} // namespace decimate::bake
