// This file is part of the decimate library; see decimate.h for version/license details.
/// @file normal_bake.cpp
/// @brief Raytraced normal bake — port of normalmap::mapnormals.
///
/// Behavior-preserving port of legacy/split/texturesample.cpp.
/// Arithmetic is identical to the 2003 original; only the container
/// types and control-flow structure are modernized.
///
/// getTestFaces is RESTORED (faithful port): per low face, cull the high mesh to
/// the faces inside the low triangle's 1.5×-scaled swept prism, then raytrace
/// only those — the dense→sparse normal-transfer optimization (Greissl 2003).
/// It is NOT a general RT accel and NOT droppable "perf code": it encodes the
/// correspondence the bake transfers.  Output is identical when high is not
/// denser than low (every dense face already sits in its own sparse triangle's
/// prism, so nothing is culled); the speedup applies when high ≫ low.

#include <decimate/bake/normal_bake.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>

#include <decimate/raytracing/raytracing.h>  // rtrace::intersect_triangle

namespace decimate::bake {

// ---------------------------------------------------------------------------
// Internal helpers — mirror the legacy static functions exactly
// ---------------------------------------------------------------------------

/// Barycentric coordinates of point (pu, pv) inside triangle (t1, t2, t3).
struct Bc3 {
    float u;
    float v;
    float w;
};

/// UV pair — corner texture coordinate in [0, texSize-1].
struct UV2 {
    float u;
    float v;
};

/// Compute barycentric coords of point (pu, pv) in triangle (t1, t2, t3).
/// Precomputed divisor bdiv = (t2-t1) × (t3-t1).
/// Mirrors legacy barycentric(float u, float v, Tex2&, Tex2&, Tex2&, float).
static Bc3 barycentricUV(float pu, float pv,
                          UV2 t1, UV2 t2, UV2 t3,
                          float bdiv) {
    Bc3 bc;
    bc.u = ((t2.u - pu) * (t3.v - pv) - (t3.u - pu) * (t2.v - pv)) / bdiv;
    bc.v = ((t3.u - pu) * (t1.v - pv) - (t1.u - pu) * (t3.v - pv)) / bdiv;
    bc.w = ((t1.u - pu) * (t2.v - pv) - (t2.u - pu) * (t1.v - pv)) / bdiv;
    return bc;
}

/// All barycentric coordinates in [0,1]?  Mirrors legacy inside().
/// The double-negation form preserves NaN behavior: NaN comparisons return
/// false, so !false = true — NaN barycentric coords would pass as "inside",
/// matching the legacy behavior.  The simplified form changes NaN semantics.
static bool insideTriangle(float u, float v, float w) {
    // NOLINTNEXTLINE(readability-simplify-boolean-expr)
    return !((u < 0.0F || v < 0.0F || w < 0.0F) ||
             (u > 1.0F || v > 1.0F || w > 1.0F));
}

/// Interpolate Vec3 using barycentric coordinates.
/// Mirrors legacy interpolate_vec3().
static decimate::math::Vec3<float> interpolateVec3(const Bc3&              bc,
                                          const decimate::math::Vec3<float>& v1,
                                          const decimate::math::Vec3<float>& v2,
                                          const decimate::math::Vec3<float>& v3) {
    return {bc.u * v1.x + bc.v * v2.x + bc.w * v3.x,
            bc.u * v1.y + bc.v * v2.y + bc.w * v3.y,
            bc.u * v1.z + bc.v * v2.z + bc.w * v3.z};
}

/// Range-compress a normalized normal into 8-bit RGB.
/// Mirrors legacy rangecompress(): nrm.norm() then (c+1)*127.
static Rgb rangeCompress(decimate::math::Vec3<float> nrm) {
    nrm.normalize();
    return {static_cast<std::uint8_t>((nrm.x + 1.0F) * 127.0F),
            static_cast<std::uint8_t>((nrm.y + 1.0F) * 127.0F),
            static_cast<std::uint8_t>((nrm.z + 1.0F) * 127.0F)};
}

/// True if point `p` lies inside the low triangle's swept prism (6 side planes).
///
/// Faithful port of legacy `contains()` (texturesample.cpp).  The OR-of-pairs
/// per edge admits the wedge swept by the triangle's edges along its vertex
/// normals.  `nrms`/`d` define the 6 planes built by `getTestFaces`.
static bool prismContains(const decimate::math::Vec3<float>& p,
                          const std::array<decimate::math::Vec3<float>, 6>& nrms,
                          const std::array<float, 6>& d) {
    const bool ba = (nrms[0].dot(p) - d[0] >= 0) || (nrms[1].dot(p) - d[1] >= 0);
    const bool cb = (nrms[2].dot(p) - d[2] >= 0) || (nrms[3].dot(p) - d[3] >= 0);
    const bool ac = (nrms[4].dot(p) - d[4] >= 0) || (nrms[5].dot(p) - d[5] >= 0);
    return ba && cb && ac;
}

/// Cull the high faces to those near low face `fi` — the dense-onto-sparse
/// normal-transfer optimization (Greissl 2003).  NOT a general RT accel: it asks
/// "which dense faces sit within THIS sparse triangle's swept envelope", which
/// is exactly the correspondence the bake transfers.  Scale the low triangle
/// 1.5× about its centroid, build 6 side planes from each edge × the adjacent
/// vertex normals, and keep any high face with ≥1 vertex inside.  `mapFaceRt`
/// then raytraces only those.  Faithful port of legacy `getTestFaces()`.
///
/// @param fi        Low face index.
/// @param lp,li     Low positions / indices.
/// @param lowNormals Per-low-vertex normals (define the sweep direction).
/// @param hp,hi,nhf High positions / indices / face count.
/// @param[out] out  Culled high-face indices (0-based face ids into hi).
static void getTestFaces(int fi,
                         const float* lp, const int* li,
                         const std::vector<decimate::math::Vec3<float>>& lowNormals,
                         const float* hp, const int* hi, int nhf,
                         std::vector<int>& out) {
    out.clear();

    const int la = li[3 * fi + 0];
    const int lb = li[3 * fi + 1];
    const int lc = li[3 * fi + 2];
    const int oa = la * 3;
    const int ob = lb * 3;
    const int oc = lc * 3;

    decimate::math::Vec3<float> a{lp[oa], lp[oa + 1], lp[oa + 2]};
    decimate::math::Vec3<float> b{lp[ob], lp[ob + 1], lp[ob + 2]};
    decimate::math::Vec3<float> c{lp[oc], lp[oc + 1], lp[oc + 2]};

    // Scale the triangle 1.5× about its centroid (legacy `scale = 1.5f`).
    const decimate::math::Vec3<float> m = (a + b + c) / 3.0F;
    constexpr float scale = 1.5F;
    a = ((a - m) * scale) + m;
    b = ((b - m) * scale) + m;
    c = ((c - m) * scale) + m;

    const decimate::math::Vec3<float> ab = a - b;
    const decimate::math::Vec3<float> bc = b - c;
    const decimate::math::Vec3<float> ca = c - a;

    // 6 side-plane normals = edge × adjacent vertex normal (legacy order).
    std::array<decimate::math::Vec3<float>, 6> nrms{
        ab.cross(lowNormals[la]).normalized(),
        ab.cross(lowNormals[lb]).normalized(),
        bc.cross(lowNormals[lb]).normalized(),
        bc.cross(lowNormals[lc]).normalized(),
        ca.cross(lowNormals[lc]).normalized(),
        ca.cross(lowNormals[la]).normalized()};
    std::array<float, 6> d{
        nrms[0].dot(a), nrms[1].dot(b), nrms[2].dot(b),
        nrms[3].dot(c), nrms[4].dot(c), nrms[5].dot(a)};

    for (int i = 0; i < nhf; i++) {
        const int ha = hi[3 * i + 0] * 3;
        const int hb = hi[3 * i + 1] * 3;
        const int hc = hi[3 * i + 2] * 3;
        const decimate::math::Vec3<float> pa{hp[ha], hp[ha + 1], hp[ha + 2]};
        const decimate::math::Vec3<float> pb{hp[hb], hp[hb + 1], hp[hb + 2]};
        const decimate::math::Vec3<float> pc{hp[hc], hp[hc + 1], hp[hc + 2]};
        if (prismContains(pa, nrms, d) || prismContains(pb, nrms, d) ||
            prismContains(pc, nrms, d)) {
            out.push_back(i);
        }
    }
}

/// Cast a ray from `orig` in direction `dir` against all high faces.
/// Returns the nearest hit normal via Möller-Trumbore (using core
/// rtrace::intersect_triangle), or falls back to `dir` if no hit.
/// Mirrors legacy rtraceHighOrderNormal() exactly.
static decimate::math::Vec3<float> rtraceHighOrderNormal(  // NOLINT(readability-function-size)
    decimate::math::Vec3<float>                       orig,
    decimate::math::Vec3<float>                       dir,
    const float*                            hp,
    const int*                              hi,
    int                                     nhf,
    const std::vector<decimate::math::Vec3<float>>&   highNormals) {
    float mint    = 100000.0F;
    float minu    = 0.0F;
    float minv    = 0.0F;
    int   minface = -1;

    for (int i = 0; i < nhf; i++) {
        int ia = hi[3 * i + 0] * 3;
        int ib = hi[3 * i + 1] * 3;
        int ic = hi[3 * i + 2] * 3;

        // intersect_triangle takes non-const refs — copy vertices into mutable locals.
        decimate::math::Vec3<float> v0{hp[ia], hp[ia + 1], hp[ia + 2]};
        decimate::math::Vec3<float> v1{hp[ib], hp[ib + 1], hp[ib + 2]};
        decimate::math::Vec3<float> v2{hp[ic], hp[ic + 1], hp[ic + 2]};

        float t = 0.0F;
        float u = 0.0F;
        float v = 0.0F;

        if (rtrace::intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &v)) {
            if (std::fabs(t) < std::fabs(mint)) {
                mint    = t;
                minface = i;
                minu    = u;
                minv    = v;
            }
        }
    }

    // Miss fallback: return the low-poly normal (the ray direction).
    if (minface < 0) {
        return dir;
    }

    // Barycentric interpolation of the high-poly normal at the hit point.
    float u = minu;
    float v = minv;
    int   a = hi[3 * minface + 0];
    int   b = hi[3 * minface + 1];
    int   c = hi[3 * minface + 2];

    return highNormals[a] * (1.0F - (u + v)) +
           highNormals[b] * u +
           highNormals[c] * v;
}

/// Rasterize one low face into the atlas and raytrace each covered texel.
/// Mirrors legacy mapfacert().
// NOLINT(readability-function-size): 10 params — 5 distinct data streams
// (low positions, low indices, corner UVs, low normals, high mesh triplet).
// Each stream is a separate array; bundling them would obscure the parallel
// low/high data-flow this function exists to connect.
static void mapFaceRt(NormalMap&                              nm,        // NOLINT(readability-function-size)
                       int                                     fi,
                       const float*                            lp,
                       const int*                              li,
                       const std::vector<decimate::math::Vec2<float>>&   cornerUV,
                       const std::vector<decimate::math::Vec3<float>>&   lowNormals,
                       const float*                            hp,
                       const int*                              hi,
                       int                                     nhf,
                       const std::vector<decimate::math::Vec3<float>>&   highNormals,
                       bool                                    tangentSpace) {
    const int texSize = nm.width;
    const float ts1   = static_cast<float>(texSize - 1);

    // Corner UV in [0,1] → scaled to [0, texSize-1]
    UV2 t1{cornerUV[3 * fi + 0].x * ts1, cornerUV[3 * fi + 0].y * ts1};
    UV2 t2{cornerUV[3 * fi + 1].x * ts1, cornerUV[3 * fi + 1].y * ts1};
    UV2 t3{cornerUV[3 * fi + 2].x * ts1, cornerUV[3 * fi + 2].y * ts1};

    // Barycentric divisor (fixed for the whole triangle)
    float bary = (t2.u - t1.u) * (t3.v - t1.v) - (t3.u - t1.u) * (t2.v - t1.v);

    // Integer UV bounding box, clamped to texture extent
    int minu = static_cast<int>(std::min({t1.u, t2.u, t3.u}));
    int maxu = static_cast<int>(std::max({t1.u, t2.u, t3.u}) + 0.5F);
    int minv = static_cast<int>(std::min({t1.v, t2.v, t3.v}));
    int maxv = static_cast<int>(std::max({t1.v, t2.v, t3.v}) + 0.5F);

    minu = std::max(minu, 0);  minu = std::min(minu, texSize - 1);
    maxu = std::max(maxu, 0);  maxu = std::min(maxu, texSize - 1);
    minv = std::max(minv, 0);  minv = std::min(minv, texSize - 1);
    maxv = std::max(maxv, 0);  maxv = std::min(maxv, texSize - 1);

    // Vertex indices and position base offsets (index-hoist: compute once, use three times)
    int va = li[3 * fi + 0];
    int vb = li[3 * fi + 1];
    int vc = li[3 * fi + 2];

    int oa = va * 3;
    int ob = vb * 3;
    int oc = vc * 3;

    const decimate::math::Vec3<float> n1 = lowNormals[va];
    const decimate::math::Vec3<float> n2 = lowNormals[vb];
    const decimate::math::Vec3<float> n3 = lowNormals[vc];

    const decimate::math::Vec3<float> pa{lp[oa], lp[oa + 1], lp[oa + 2]};
    const decimate::math::Vec3<float> pb{lp[ob], lp[ob + 1], lp[ob + 2]};
    const decimate::math::Vec3<float> pc{lp[oc], lp[oc + 1], lp[oc + 2]};

    // Tangent-space option: per-face tangent from the UV gradient (constant over
    // the face). Stored normals become relative to the surface frame → near-flat
    // ≈(0,0,1), so chart values barely differ across a seam.
    decimate::math::Vec3<float> tFace{0.0F, 0.0F, 0.0F};
    bool              tValid = false;
    if (tangentSpace) {
        const decimate::math::Vec2<float>& uA = cornerUV[3 * fi + 0];
        const decimate::math::Vec2<float>& uB = cornerUV[3 * fi + 1];
        const decimate::math::Vec2<float>& uC = cornerUV[3 * fi + 2];
        const float du1 = uB.x - uA.x; const float dv1 = uB.y - uA.y;
        const float du2 = uC.x - uA.x; const float dv2 = uC.y - uA.y;
        const float det = du1 * dv2 - du2 * dv1;
        if (std::fabs(det) > 1e-12F) {
            const float r = 1.0F / det;
            tFace = (pb - pa) * (dv2 * r) - (pc - pa) * (dv1 * r);
            tValid = tFace.lengthSquared() > 1e-12F;
        }
    }

    for (int x = minu; x <= maxu; x++) {
        for (int y = minv; y <= maxv; y++) {
            Bc3 bc = barycentricUV(static_cast<float>(x) + 0.5F, static_cast<float>(y) + 0.5F,
                                    t1, t2, t3, bary);

            if (!insideTriangle(bc.u, bc.v, bc.w)) {
                continue;
            }

            // Interpolated low-poly normal (ray direction) and surface point
            decimate::math::Vec3<float> n = interpolateVec3(bc, n1, n2, n3);
            decimate::math::Vec3<float> p = interpolateVec3(bc, pa, pb, pc);

            // Cast the ray against the high mesh; range-compress the result
            decimate::math::Vec3<float> snrm = rtraceHighOrderNormal(p, n, hp, hi, nhf, highNormals);
            snrm.normalize();  // mirrors mapfacert snrm.norm() before rangecompress

            if (tangentSpace && tValid) {
                // Express snrm in the per-texel tangent frame (N = interpolated low
                // normal, T = Gram-Schmidt'd face tangent, B = N×T).
                decimate::math::Vec3<float> nN = n;       nN.normalize();
                decimate::math::Vec3<float> tT = tFace - nN * nN.dot(tFace);
                tT.normalize();
                const decimate::math::Vec3<float> bB = nN.cross(tT);
                snrm = decimate::math::Vec3<float>{snrm.dot(tT), snrm.dot(bB), snrm.dot(nN)};
            }

            nm.pixels[y * texSize + x] = rangeCompress(snrm);
        }
    }
}

/// One dilation pass: read from src, write into dst.
/// Background texels receive the integer average of non-background 8-neighbours.
/// Filled texels and isolated background texels copy through unchanged.
/// Mirrors legacy dilate() inner loop: double-buffer, 8-neighbour, % wrap.
// NOLINTNEXTLINE(readability-function-cognitive-complexity,readability-function-size)
static void dilateOnce(int texsize,
                        const std::vector<Rgb>& src, std::vector<Rgb>& dst) {
    const Rgb background{128, 128, 128};

    for (int y = 0; y < texsize; y++) {
        for (int x = 0; x < texsize; x++) {
            int idx = y * texsize + x;

            const auto& cur = src[idx];
            auto&       out = dst[idx];

            // Filled texel: copy through unchanged.
            if (cur != background) {
                out = cur;
                continue;
            }

            // Background texel: integer average of non-background 8-neighbours.
            // Mirrors legacy dilate() % wrap + negative-index fixup exactly.
            int accu[3] = {0, 0, 0};  // NOLINT(modernize-avoid-c-arrays)
            int acnt    = 0;

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

                    const auto& nb = src[yn * texsize + xn];
                    if (nb != background) {
                        accu[0] += nb.r;
                        accu[1] += nb.g;
                        accu[2] += nb.b;
                        acnt++;
                    }
                }
            }

            if (acnt > 0) {
                out.r = static_cast<std::uint8_t>(accu[0] / acnt);
                out.g = static_cast<std::uint8_t>(accu[1] / acnt);
                out.b = static_cast<std::uint8_t>(accu[2] / acnt);
            } else {
                out = cur;  // stay background
            }
        }
    }
}

/// `passes`-pass dilation: fill background texels from non-background neighbours.
/// Double-buffer, 8-neighbour average, % wrap. NOT part of the 2003 pipeline —
/// this is the opt-in gutter extension (see bake()'s `dilatePasses`); each pass
/// grows the filled region by 1 texel.
static void dilate(NormalMap& nm, int passes) {
    const int texsize = nm.width;

    std::vector<Rgb> scratch(texsize * texsize);
    std::vector<Rgb>* img = &nm.pixels;
    std::vector<Rgb>* scr = &scratch;

    for (int pass = 0; pass < passes; pass++) {
        dilateOnce(texsize, *img, *scr);
        std::swap(img, scr);
    }

    // Odd pass count leaves the result in `scratch`; copy it back.
    if (img != &nm.pixels) {
        nm.pixels = *img;
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

NormalMap bake(const decimate::MeshData&                   low,
               const std::vector<decimate::math::Vec2<float>>&   lowCornerUV,
               const std::vector<decimate::math::Vec3<float>>&   lowNormals,
               const decimate::MeshData&                   high,
               const std::vector<decimate::math::Vec3<float>>&   highNormals,
               int                                      texSize,
               int                                      dilatePasses,
               bool                                     tangentSpace) {
    NormalMap nm;
    nm.width  = texSize;
    nm.height = texSize;
    nm.pixels.assign(static_cast<std::size_t>(texSize) * texSize, Rgb{128, 128, 128});

    // positions is std::vector<Vec3f>; Vec3f is standard-layout and tightly
    // packed, so reinterpret to the flat coord buffer the cull/rtrace machinery
    // (getTestFaces / mapFaceRt / rtrace) indexes as lp[3*v+k] — unchanged.
    const float* lp = reinterpret_cast<const float*>(low.positions.data());
    const int*   li = low.indices.data();
    const int    nlf = static_cast<int>(low.indices.size()) / 3;

    const float* hp = reinterpret_cast<const float*>(high.positions.data());
    const int*   hi = high.indices.data();
    const int    nhf = static_cast<int>(high.indices.size()) / 3;

    // Per-low-face: cull the high mesh to the faces within this triangle's swept
    // prism (getTestFaces), then raytrace only those — the dense-onto-sparse
    // optimization.  `testFace` collects culled high-face ids; `testHi` packs
    // their 3 indices into the contiguous form mapFaceRt/rtrace consume (it is
    // the legacy `testFaces` index array).  Both reused across faces.
    // Pure RT: raytrace EVERY high face per low face (no prism cull). The prism
    // optimization (getTestFaces) is LOSSY on sparse/decimated low meshes — the
    // 1.5x swept envelope can exclude the correct high face, producing speckles
    // (measured: speckles vanish with the cull off). It is kept below, gated off,
    // until the Spec-13 SAH-BVH makes brute-force RT fast enough to retire it.
    // Flip to false to restore the (faster, lossy) culled path on dense meshes.
    constexpr bool kBakeDisableCulling = true;

    std::vector<int> testFace;
    std::vector<int> testHi;
    for (int fi = 0; fi < nlf; fi++) {
        if (kBakeDisableCulling) {
            mapFaceRt(nm, fi, lp, li, lowCornerUV, lowNormals, hp, hi, nhf, highNormals, tangentSpace);
            continue;
        }
        getTestFaces(fi, lp, li, lowNormals, hp, hi, nhf, testFace);
        testHi.clear();
        testHi.reserve(testFace.size() * 3);
        for (const int hf : testFace) {
            testHi.push_back(hi[3 * hf + 0]);
            testHi.push_back(hi[3 * hf + 1]);
            testHi.push_back(hi[3 * hf + 2]);
        }
        mapFaceRt(nm, fi, lp, li, lowCornerUV, lowNormals, hp,
                  testHi.data(), static_cast<int>(testFace.size()), highNormals, tangentSpace);
    }

    // Original 2003 output has NO dilation (charts crisp; no normals outside a
    // chart). dilatePasses > 0 is the opt-in gutter extension for seam-free mesh
    // texturing — the packer must reserve matching clearance (see packAtlas).
    if (dilatePasses > 0) {
        dilate(nm, dilatePasses);
    }
    return nm;
}

// ---------------------------------------------------------------------------
// Seam-correct gutter fill helpers
// ---------------------------------------------------------------------------

/// A seam edge in atlas space: two 2D endpoints (in texel coordinates) for
/// each side (A = owning half-edge, B = sym half-edge).
///
/// t ∈ [0,1] parameterizes along the edge from p0 to p1.
/// The sym edge runs p1→p0, so t along A maps to (1-t) along B.
struct SeamEdge {
    decimate::math::Vec2<float> p0;   ///< Atlas start (texel coords), origin vertex of he.
    decimate::math::Vec2<float> p1;   ///< Atlas end   (texel coords), end vertex of he.
    decimate::math::Vec2<float> inN;  ///< Inward unit normal (into chart A's interior).

    // Sym edge (B side) — same endpoints but REVERSED order.
    decimate::math::Vec2<float> bp0;  ///< Atlas start for B (= A's p1 direction at t=0 on B).
    decimate::math::Vec2<float> bp1;  ///< Atlas end   for B.
    decimate::math::Vec2<float> bInN; ///< Inward unit normal for B.

    bool hasRefMap = false;  ///< true = chart seam (neighbour exists); false = open mesh.
};

/// Project point `g` onto the line segment [p0, p1].
/// Returns (t, d²) where t = clamped param in [0,1] and d² = squared distance.
static std::pair<float, float> projectOnSegment(decimate::math::Vec2<float> g,
                                                 decimate::math::Vec2<float> p0,
                                                 decimate::math::Vec2<float> p1) {
    float ex = p1.x - p0.x;
    float ey = p1.y - p0.y;
    float len2 = ex * ex + ey * ey;
    if (len2 < 1e-12F) {
        float dx = g.x - p0.x;
        float dy = g.y - p0.y;
        return {0.0F, dx * dx + dy * dy};
    }
    float t = ((g.x - p0.x) * ex + (g.y - p0.y) * ey) / len2;
    t = std::max(0.0F, std::min(1.0F, t));
    float fx = p0.x + t * ex - g.x;
    float fy = p0.y + t * ey - g.y;
    return {t, fx * fx + fy * fy};
}

/// Bilinear sample the atlas at floating-point texel coordinates (u, v).
/// Clamps to [0, texSize-1]². Returns background sentinel on out-of-range.
static Rgb bilinearSample(const std::vector<Rgb>& pixels, int texSize, float u, float v) {
    const Rgb kBg{128, 128, 128};
    // Clamp to valid range
    u = std::max(0.0F, std::min(static_cast<float>(texSize - 1), u));
    v = std::max(0.0F, std::min(static_cast<float>(texSize - 1), v));

    int x0 = static_cast<int>(u);
    int y0 = static_cast<int>(v);
    int x1 = std::min(x0 + 1, texSize - 1);
    int y1 = std::min(y0 + 1, texSize - 1);

    float fu = u - static_cast<float>(x0);
    float fv = v - static_cast<float>(y0);

    const Rgb& c00 = pixels[y0 * texSize + x0];
    const Rgb& c10 = pixels[y0 * texSize + x1];
    const Rgb& c01 = pixels[y1 * texSize + x0];
    const Rgb& c11 = pixels[y1 * texSize + x1];

    // If all four taps are background, return background.
    if (c00 == kBg && c10 == kBg && c01 == kBg && c11 == kBg) {
        return kBg;
    }

    // Bilinear blend — skip background taps by treating them as transparent.
    // Simple approach: blend only non-background taps.
    float r = 0.0F;
    float g = 0.0F;
    float b = 0.0F;
    float w = 0.0F;

    auto addTap = [&](const Rgb& c, float weight) {
        if (c == kBg) {
            return;
        }
        r += weight * static_cast<float>(c.r);
        g += weight * static_cast<float>(c.g);
        b += weight * static_cast<float>(c.b);
        w += weight;
    };

    addTap(c00, (1.0F - fu) * (1.0F - fv));
    addTap(c10, fu * (1.0F - fv));
    addTap(c01, (1.0F - fu) * fv);
    addTap(c11, fu * fv);

    if (w < 1e-6F) {
        return kBg;
    }
    return {static_cast<std::uint8_t>(r / w),
            static_cast<std::uint8_t>(g / w),
            static_cast<std::uint8_t>(b / w)};
}

/// Build the UV pair for a half-edge (texel-scale coordinates).
/// he = 3*f + k.
/// Returns {p0 = origin UV, p1 = end UV} scaled to [0, texSize-1].
/// Origin vertex UV: writeFace convention:
///   k=0 (ac): origin=a (cornerUV[3f+0]), end=c (cornerUV[3f+2])
///   k=1 (cb): origin=c (cornerUV[3f+2]), end=b (cornerUV[3f+1])
///   k=2 (ba): origin=b (cornerUV[3f+1]), end=a (cornerUV[3f+0])
/// Third vertex (opposite, inside the face):
///   k=0: opposite = b (cornerUV[3f+1])
///   k=1: opposite = a (cornerUV[3f+0])
///   k=2: opposite = c (cornerUV[3f+2])
static std::tuple<decimate::math::Vec2<float>, decimate::math::Vec2<float>, decimate::math::Vec2<float>>
heUV(decimate::HeIndex he, const std::vector<decimate::math::Vec2<float>>& cornerUV, float ts1) {
    int f = static_cast<int>(he / 3u);
    int k = static_cast<int>(he % 3u);

    // UV indices for each slot
    //   k=0: origin=cornerUV[3f+0], end=cornerUV[3f+2], opp=cornerUV[3f+1]
    //   k=1: origin=cornerUV[3f+2], end=cornerUV[3f+1], opp=cornerUV[3f+0]
    //   k=2: origin=cornerUV[3f+1], end=cornerUV[3f+0], opp=cornerUV[3f+2]
    int origIdx = (k == 0) ? 0 : (k == 1) ? 2 : 1;
    int endIdx  = (k == 0) ? 2 : (k == 1) ? 1 : 0;
    int oppIdx  = (k == 0) ? 1 : (k == 1) ? 0 : 2;

    auto scale = [&](decimate::math::Vec2<float> uv) {
        return decimate::math::Vec2<float>{uv.x * ts1, uv.y * ts1};
    };

    return {scale(cornerUV[3 * f + origIdx]),
            scale(cornerUV[3 * f + endIdx]),
            scale(cornerUV[3 * f + oppIdx])};
}

/// Compute 2D inward unit normal for edge [p0, p1] with opposite vertex `opp`.
/// Inward = perpendicular pointing toward `opp`.
static decimate::math::Vec2<float> inwardNormal2D(decimate::math::Vec2<float> p0,
                                         decimate::math::Vec2<float> p1,
                                         decimate::math::Vec2<float> opp) {
    float ex = p1.x - p0.x;
    float ey = p1.y - p0.y;
    float len = std::sqrt(ex * ex + ey * ey);
    if (len < 1e-12F) {
        return {0.0F, 0.0F};
    }
    // Two candidate perpendiculars: (−ey, ex) and (ey, −ex).
    float nx = -ey / len;
    float ny =  ex / len;
    // Dot with (opp - p0): choose the side that matches opp.
    float dot = nx * (opp.x - p0.x) + ny * (opp.y - p0.y);
    if (dot < 0.0F) {
        nx = -nx;
        ny = -ny;
    }
    return {nx, ny};
}

/// Build the seam-edge catalog from the half-edge mesh and corner UVs.
/// Returns a flat list of SeamEdge (only boundary half-edges, i.e. those
/// where the mesh boundary OR chart seam lives — the set of he with no sym
/// OR with a sym in a DIFFERENT chart).
///
/// Also returns the atlas-space seam segment per half-edge for fast lookup.
/// refMap: he → SeamEdge index (only for he that are chart-boundary).
static std::unordered_map<decimate::HeIndex, SeamEdge>
buildSeamEdges(const decimate::HalfEdgeMesh& hem,
               const std::vector<decimate::math::Vec2<float>>& cornerUV,
               const std::vector<int>& faceToChart,
               float ts1) {
    std::unordered_map<decimate::HeIndex, SeamEdge> seams;
    const auto& hes = hem.halfedges();

    for (decimate::HeIndex i = 0; i < static_cast<decimate::HeIndex>(hes.size()); ++i) {
        const auto& h = hes[i];
        decimate::HeIndex s = h.sym;

        // Only process half-edges that lie on a CHART BOUNDARY
        // (either open-mesh boundary OR different chart).
        bool isOpenBoundary = (s == decimate::kNoHalfedge);
        bool isChartSeam = false;

        if (!isOpenBoundary && s < static_cast<decimate::HeIndex>(hes.size())) {
            int chartA = faceToChart[h.face];
            int chartB = faceToChart[hes[s].face];
            isChartSeam = (chartA != chartB);
        }

        if (!isOpenBoundary && !isChartSeam) {
            continue;  // interior edge — skip
        }

        // Build the SeamEdge for half-edge i.
        auto [p0, p1, opp] = heUV(i, cornerUV, ts1);
        decimate::math::Vec2<float> inN = inwardNormal2D(p0, p1, opp);

        SeamEdge se;
        se.p0   = p0;
        se.p1   = p1;
        se.inN  = inN;
        se.hasRefMap = isChartSeam;

        if (isChartSeam) {
            // sym edge: runs p1_b → p0_b in B's atlas.
            // B's p0 = A's p1 direction-wise (sym reversal).
            auto [bp0, bp1, bopp] = heUV(s, cornerUV, ts1);
            se.bp0  = bp0;
            se.bp1  = bp1;
            se.bInN = inwardNormal2D(bp0, bp1, bopp);
        }

        seams.emplace(i, se);
    }
    return seams;
}

// ---------------------------------------------------------------------------
// bakeSeamFill — public API
// ---------------------------------------------------------------------------

NormalMap bakeSeamFill(const decimate::MeshData&                   low,
                       const decimate::HalfEdgeMesh&               hem,
                       const std::vector<decimate::math::Vec2<float>>&   lowCornerUV,
                       const std::vector<int>&                 faceToChart,
                       const std::vector<decimate::math::Vec3<float>>&   lowNormals,
                       const decimate::MeshData&                   high,
                       const std::vector<decimate::math::Vec3<float>>&   highNormals,
                       int                                      texSize,
                       int                                      dilatePasses,
                       bool                                     tangentSpace) {
    // --- Pass 1: interior bake (unchanged behavior, no dilation) -------
    NormalMap nm = bake(low, lowCornerUV, lowNormals, high, highNormals, texSize, 0, tangentSpace);

    // --- Build seam-edge catalog from HEM + cornerUV + faceToChart -----
    // Seam detection:
    // A mesh edge is a CHART SEAM iff sym(he) != kNoHalfedge AND the two
    // faces' charts differ (same 3D edge, different atlas locations).
    // Implementation cites: Ray et al. 2010 "Invisible Seams" indirection map.
    // NOT the MS patent US10,839,590 polygon-projection mechanism.
    const float ts1 = static_cast<float>(texSize - 1);
    auto seams = buildSeamEdges(hem, lowCornerUV, faceToChart, ts1);

    if (seams.empty()) {
        // No chart seams (single-chart mesh or trivial case): fall back to dilate.
        if (dilatePasses > 0) {
            dilate(nm, dilatePasses);
        }
        return nm;
    }

    // Flatten the seam map to a vector for O(E) nearest-edge search per texel.
    std::vector<std::pair<decimate::HeIndex, SeamEdge>> seamVec(seams.begin(), seams.end());

    // --- Pass 2: seam-correct gutter fill ------------------------------
    const Rgb kBg{128, 128, 128};
    std::vector<Rgb> pass2 = nm.pixels;  // mutable output; pass1 stays frozen

    // Gutter distance: how many texels out from the seam edge we fill.
    // Match dilatePasses as the gutter width (same semantics as blur dilate).
    // For the seam path, we fill texels within `dilatePasses` perpendicular
    // distance of a seam edge, regardless of 8-neighbour fill status.
    // SIMPLE CORRECTNESS PATH: fill every background texel whose nearest seam
    // edge is a chart seam → sample neighbour; nearest is open → classic dilate.
    // One O(texels × seam_edges) pass — KISS for this prototype.
    for (int y = 0; y < texSize; ++y) {
        for (int x = 0; x < texSize; ++x) {
            const int idx = y * texSize + x;
            if (nm.pixels[idx] != kBg) {
                continue;  // already filled in pass 1 — copy through
            }

            decimate::math::Vec2<float> g{static_cast<float>(x), static_cast<float>(y)};

            // Find the nearest seam edge and foot parameters.
            float bestDist2 = std::numeric_limits<float>::max();
            float bestT     = 0.0F;
            bool  bestHasRef = false;
            const SeamEdge* bestEdge = nullptr;

            for (const auto& [he, se] : seamVec) {
                auto [t, d2] = projectOnSegment(g, se.p0, se.p1);
                if (d2 < bestDist2) {
                    bestDist2  = d2;
                    bestT      = t;
                    bestHasRef = se.hasRefMap;
                    bestEdge   = &se;
                }
            }

            if (bestEdge == nullptr) {
                continue;  // no seam edge at all — leave background
            }

            const float gDist = std::sqrt(bestDist2);  // perpendicular distance d

            // Bound the fill to the gutter band: only texels within `dilatePasses`
            // of a seam edge are filled. Without this the loop floods the ENTIRE
            // empty atlas (every background texel finds *some* nearest seam edge),
            // smearing far regions and producing seam-like artefacts near borders.
            if (gDist > static_cast<float>(dilatePasses)) {
                continue;  // beyond the gutter — leave background
            }

            if (!bestHasRef) {
                // Open-mesh boundary: classic dilate will handle this texel.
                // Leave it as background for the blur pass below.
                continue;
            }

            // Seam edge: sample neighbour chart B at the sym point.
            // Correctness invariants (concept-seamless-bake.md):
            //   1. Param reversal: t along e_a → (1-t) along e_b.
            //   2. Perpendicular flip: d outside A = d INSIDE B (bInN points in).
            float tB = 1.0F - bestT;
            float srcU = bestEdge->bp0.x + tB * (bestEdge->bp1.x - bestEdge->bp0.x)
                         + bestEdge->bInN.x * gDist;
            float srcV = bestEdge->bp0.y + tB * (bestEdge->bp1.y - bestEdge->bp0.y)
                         + bestEdge->bInN.y * gDist;

            Rgb sampled = bilinearSample(nm.pixels, texSize, srcU, srcV);
            if (sampled != kBg) {
                pass2[idx] = sampled;
            }
            // else: sample landed in background of B (very thin chart) — leave
            // for the classic dilate fallback below.
        }
    }

    nm.pixels = std::move(pass2);

    // Blur-dilate for open-mesh boundaries and any remaining background.
    if (dilatePasses > 0) {
        dilate(nm, dilatePasses);
    }
    return nm;
}

}  // namespace decimate::bake
