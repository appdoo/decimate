// This file is part of the decimate library; see decimate.h for version/license details.
/// @file seam_straighten_test.cpp
/// @brief Spec 17-02: seam straightening (Liu et al. 2017, section 6) + the cutter's
///        straighten-threshold parameter.
///
/// Scenarios I1, I2, B1, B2, S1-S8 (seam straightening).
///
/// The NOLINT on the TEST_CASEs is for Catch2: every CHECK macro expands to branching code
/// that readability-function-size / -cognitive-complexity count against the test body.

#include "model_fixture.h"
#include "sat_overlap.h"

#include <algorithm>
#include <array>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normals.h>
#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/reduce/reducer.h>
#include <decimate/seam/erase.h>
#include <decimate/seam/seam_topology.h>
#include <decimate/seam/straighten.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/textured_mesh.h>
#include <functional>
#include <map>
#include <numbers>
#include <string>
#include <vector>

namespace {

using decimate::MeshData;
using decimate::TexturedMesh;
using decimate::math::Vec2f;
using decimate::math::Vec3f;
using decimate::seam::StraightenOptions;
using decimate::seam::StraightenResult;

/// Signed int to container index.
std::size_t idx(int i) {
    return static_cast<std::size_t>(i);
}


/// Twice the signed UV area of @p face (positive = counter-clockwise, i.e. not flipped).
double signedArea(const TexturedMesh &tm, std::size_t face) {
    const Vec2f a = tm.uvs[idx(tm.faceU[3 * face])];
    const Vec2f b = tm.uvs[idx(tm.faceU[(3 * face) + 1])];
    const Vec2f c = tm.uvs[idx(tm.faceU[(3 * face) + 2])];
    return ((static_cast<double>(b.x) - a.x) * (static_cast<double>(c.y) - a.y)) -
           ((static_cast<double>(b.y) - a.y) * (static_cast<double>(c.x) - a.x));
}


/// Bit-identical UV (a protected vertex must not move by even one ulp).
bool sameBits(Vec2f a, Vec2f b) {
    return std::bit_cast<std::uint32_t>(a.x) == std::bit_cast<std::uint32_t>(b.x) &&
           std::bit_cast<std::uint32_t>(a.y) == std::bit_cast<std::uint32_t>(b.y);
}


/// Largest distance any UV vertex moved.
double maxMove(const TexturedMesh &before, const TexturedMesh &after) {
    double worst = 0.0;
    for (std::size_t u = 0; u < before.uvs.size(); ++u) {
        worst = std::max(worst, std::hypot(static_cast<double>(after.uvs[u].x) - before.uvs[u].x,
                                           static_cast<double>(after.uvs[u].y) - before.uvs[u].y));
    }
    return worst;
}


/// UV vertices sitting on surface vertex @p sv, ordered by (u, v).
std::vector<int> copiesOf(const TexturedMesh &tm, int sv) {
    std::vector<int> out;
    for (std::size_t i = 0; i < tm.faceV.size(); ++i) {
        if (tm.faceV[i] == sv && std::ranges::find(out, tm.faceU[i]) == out.end()) {
            out.push_back(tm.faceU[i]);
        }
    }
    std::ranges::sort(out, [&](int a, int b) {
        return std::pair(tm.uvs[idx(a)].x, tm.uvs[idx(a)].y) <
               std::pair(tm.uvs[idx(b)].x, tm.uvs[idx(b)].y);
    });
    return out;
}


/// Are all UV copies of surface vertex @p sv bit-identical before and after?
bool copiesUnchanged(const TexturedMesh &before, const TexturedMesh &after, int sv) {
    return std::ranges::all_of(copiesOf(before, sv), [&](int u) {
        return sameBits(before.uvs[idx(u)], after.uvs[idx(u)]);
    });
}


/// Planar grid fixture: unit-grid triangles, a chart label per face, one UV offset
/// per chart.  UV = 0.05 + 0.1 * (x, y) + offset(chart).  Cells listed in `flipped` use the
/// other diagonal; that is how fixtures keep a triangle from touching two silhouette edges
/// (which would make it an "ear" and protect its vertices).
struct Fixture {
    MeshData           mesh;
    std::vector<Vec2f> cornerUV;
    std::vector<int>   chartOfFace;
    TexturedMesh       tm;
};


/// Maps (face index, centroid x, centroid y) to the chart that face belongs to.
using ChartFn = std::function<int(int face, double cx, double cy)>;
using Cells = std::vector<std::pair<int, int>>;

/// The two triangles of cell (a,b,c,d), split along a-c or, flipped, b-d.
std::array<std::array<int, 3>, 2> splitCell(int a, int b, int c, int d, bool flip) {
    if (flip) {
        return { { { a, b, d }, { b, c, d } } };
    }
    return { { { a, b, c }, { a, c, d } } };
}


/// Appends triangle @p t to @p fx, labelled by @p chartOf and placed at its chart's UV offset.
void addTriangle(Fixture &fx, const std::array<int, 3> &t, const ChartFn &chartOf,
                 const std::vector<Vec2f> &offset) {
    double cx = 0.0;
    double cy = 0.0;
    for (const int v : t) {
        cx += fx.mesh.positions[idx(v)].x / 3.0;
        cy += fx.mesh.positions[idx(v)].y / 3.0;
    }
    const int chart = chartOf(static_cast<int>(fx.chartOfFace.size()), cx, cy);
    fx.chartOfFace.push_back(chart);
    for (const int v : t) {
        fx.mesh.indices.push_back(v);
        const Vec3f p = fx.mesh.positions[idx(v)];
        const Vec2f o = offset[idx(chart)];
        fx.cornerUV.emplace_back(0.05f + (0.1f * p.x) + o.x, 0.05f + (0.1f * p.y) + o.y);
    }
}


/// A w x h unit grid of cells, two triangles each; the chart layout is given by @p chartOf.
Fixture gridFixture(int w, int h, const ChartFn &chartOf, const std::vector<Vec2f> &offset,
                    const Cells &flipped) {
    Fixture fx;
    for (int j = 0; j <= h; ++j) {
        for (int i = 0; i <= w; ++i) {
            fx.mesh.positions.emplace_back(static_cast<float>(i), static_cast<float>(j), 0.0f);
        }
    }
    const auto vid = [&](int i, int j) {
        return (j * (w + 1)) + i;
    };
    for (int cell = 0; cell < w * h; ++cell) {
        const int  i = cell % w;
        const int  j = cell / w;
        const bool flip = std::ranges::find(flipped, std::pair(i, j)) != flipped.end();
        const auto tris =
            splitCell(vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1), flip);
        addTriangle(fx, tris[0], chartOf, offset);
        addTriangle(fx, tris[1], chartOf, offset);
    }
    fx.tm = decimate::seam::weld(fx.mesh, fx.cornerUV);
    return fx;
}


/// Pushes the interior chain vertices (2,1)..(2,3) of both copies @p noise away from the chain.
void bendChain(Fixture &fx, float noise) {
    for (std::size_t i = 0; i < fx.mesh.indices.size(); ++i) {
        const Vec3f p = fx.mesh.positions[idx(fx.mesh.indices[i])];
        if (p.x == 2.0f && p.y >= 1.0f && p.y <= 3.0f) {
            fx.cornerUV[i].x += fx.chartOfFace[i / 3] == 0 ? -noise : noise;
        }
    }
    fx.tm = decimate::seam::weld(fx.mesh, fx.cornerUV);
}


/// Two charts, a vertical 4-edge seam at x = 2 (surface vertices (2,0)..(2,4)); the cells
/// are split so that no triangle touches two silhouette edges (no ears).  With @p noise
/// the three interior chain vertices of both copies are pushed @p noise perpendicular to
/// the chain.
Fixture stripFixture(float noise) {
    Fixture fx = gridFixture(4, 4, [](int, double cx, double) { return cx < 2.0 ? 0 : 1; },
                             { { 0.0f, 0.05f }, { 0.3f, 0.05f } }, { { 3, 0 }, { 0, 3 }, { 1, 0 }, { 2, 3 } });
    if (noise != 0.0f) {
        bendChain(fx, noise);
    }
    return fx;
}


/// Surface vertex id of chain vertex (2, j) in the 5x5 strip fixture.
int stripVertex(int j) {
    return (j * 5) + 2;
}


/// Chart A has a spike (t1,t2,t3 with the tip t2 pointing up) holding the interior vertex
/// w; chart B sits above.  Flattening the chain pulls the tip below w.
Fixture foldFixture(float wy) {
    Fixture fx;

    // 0..4 chain t0..t4, 5 b0, 6 b1, 7 w, 8..10 u0..u2
    fx.mesh.positions = { { 0, 1, 0 }, { 0.25f, 1, 0 }, { 0.5f, 1.4f, 0 }, { 0.75f, 1, 0 }, { 1, 1, 0 }, { 0, 0, 0 }, { 1, 0, 0 }, { 0.5f, wy, 0 }, { 0, 2, 0 }, { 0.5f, 2, 0 }, { 1, 2, 0 } };
    const std::vector<std::array<int, 3>> a = { { 7, 3, 2 }, { 7, 2, 1 }, { 7, 1, 3 }, { 5, 6, 3 }, { 6, 4, 3 }, { 5, 3, 1 }, { 5, 1, 0 } };
    const std::vector<std::array<int, 3>> b = { { 0, 1, 8 }, { 1, 9, 8 }, { 1, 2, 9 }, { 2, 3, 9 }, { 3, 10, 9 }, { 3, 4, 10 } };
    const auto                            add = [&](const std::array<int, 3> &t, int chart) {
        fx.chartOfFace.push_back(chart);
        for (const int v : t) {
            fx.mesh.indices.push_back(v);
            const Vec3f p = fx.mesh.positions[idx(v)];
            fx.cornerUV.emplace_back(0.05f + (0.3f * p.x) + (chart == 0 ? 0.0f : 0.5f),
                                                                0.05f + (0.3f * p.y));
        }
    };

    for (const auto &t : a) {
        add(t, 0);
    }
    for (const auto &t : b) {
        add(t, 1);
    }
    fx.tm = decimate::seam::weld(fx.mesh, fx.cornerUV);
    return fx;
}


using decimate::bake::NormalMap;

/// The reduce -> segment -> parametrize -> weld pipeline of one model.
struct Pipeline {
    MeshData                  high;
    MeshData                  low;
    decimate::HalfEdgeMesh    hem;
    decimate::chart::ChartSet charts;
    std::vector<Vec2f>        cornerUV;
    TexturedMesh              tm;
};


constexpr int kTex = 256;
constexpr int kLowVertices = 1500;

/// Runs the pipeline on the model at @p path.
Pipeline buildPipeline(const std::string &path) {
    Pipeline p;
    p.high = decimate::load(path);
    p.low = decimate::reduce(p.high, kLowVertices);
    std::vector<decimate::Face> faces;
    for (std::size_t i = 0; i + 2 < p.low.indices.size(); i += 3) {
        faces.push_back({ p.low.indices[i], p.low.indices[i + 1], p.low.indices[i + 2] });
    }
    p.hem = decimate::HalfEdgeMesh::create(static_cast<int>(p.low.positions.size()),
                                           static_cast<int>(faces.size()), faces.data());
    p.charts = decimate::chart::segment(p.hem, p.low, 64);
    p.cornerUV = decimate::chart::parametrizeCorners(p.hem, p.low, p.charts, kTex);
    p.tm = decimate::seam::weld(p.low, p.cornerUV);
    return p;
}


/// The shared, lazily built spot pipeline (read-only).
const Pipeline &spot() {
    static const Pipeline kSpot = buildPipeline(std::string(MESH_MODEL_DIR) + "/spot.obj");
    return kSpot;
}


/// FNV-1a over a list of ints (golden fingerprint of the cutter's output).
std::uint64_t fingerprint(const std::vector<int> &values) {
    std::uint64_t h = 1469598103934665603ULL;
    for (const int v : values) {
        h = (h ^ static_cast<std::uint32_t>(v)) * 1099511628211ULL;
    }
    return h;
}


/// Per chart: seed face and face count -- the rest of what `segment()` returns.
std::vector<int> chartSummary(const decimate::chart::ChartSet &cs) {
    std::vector<int> out;
    for (const auto &c : cs.charts) {
        out.push_back(c.seed);
        out.push_back(static_cast<int>(c.faces.size()));
    }
    return out;
}


/// Edge-length-ratio agreement of the two copies of a seam, over all unprotected chain vertices.
struct RatioStats {
    int    pairs = 0;         ///< Adjacent seam-edge pairs at unprotected chain vertices.
    int    parallel = 0;      ///< ... of which the 4D (both copies) directions are anti-parallel.
    double worstExcess = 0.0; ///< Largest (mismatch - tolerance) over the parallel pairs.
    double worstAll = 0.0;    ///< Largest relative mismatch over all pairs.
};


/// One copy's UV edge leaving a seam vertex: the UV vertex at v and the edge vector.
struct Side {
    int    uvAtV = 0;
    double dx = 0.0;
    double dy = 0.0;
};


/// The two UV edges of seam @p pair leaving surface vertex @p v, keyed by the UV vertex at v.
std::array<Side, 2> sidesAt(const TexturedMesh &tm, const decimate::seam::SeamPair &pair, int v) {
    std::array<Side, 2> out {};
    std::size_t         k = 0;
    for (const auto &he : { pair.a, pair.b }) {
        const std::size_t base = 3 * idx(he.face);
        const std::size_t c0 = base + idx(he.corner);
        const std::size_t c1 = base + ((idx(he.corner) + 1) % 3);
        const bool        startsAtV = tm.faceV[c0] == v;
        const int         atV = startsAtV ? tm.faceU[c0] : tm.faceU[c1];
        const int         far = startsAtV ? tm.faceU[c1] : tm.faceU[c0];
        out[k++] = { atV, static_cast<double>(tm.uvs[idx(far)].x) - tm.uvs[idx(atV)].x,
                     static_cast<double>(tm.uvs[idx(far)].y) - tm.uvs[idx(atV)].y };
    }
    return out;
}


/// Seam pairs and boundary edges incident to each surface vertex.
struct VertexEdges {
    std::map<int, std::vector<const decimate::seam::SeamPair *>> seams;
    std::map<int, int>                                           boundary;
};


/// Indexes the seams and boundary edges of @p tm by surface vertex.
VertexEdges edgesPerVertex(const TexturedMesh &tm, const decimate::seam::SeamTopology &topo) {
    const auto vertexAt = [&](decimate::seam::HalfEdgeRef he, std::size_t k) {
        return tm.faceV[(3 * idx(he.face)) + ((idx(he.corner) + k) % 3)];
    };
    VertexEdges out;
    for (const auto &s : topo.seams) {
        out.seams[vertexAt(s.a, 0)].push_back(&s);
        out.seams[vertexAt(s.a, 1)].push_back(&s);
    }
    for (const auto &b : topo.boundary) {
        ++out.boundary[vertexAt(b, 0)];
        ++out.boundary[vertexAt(b, 1)];
    }
    return out;
}


/// Edge length of one side.
double length(const Side &s) {
    return std::hypot(s.dx, s.dy);
}


/// Folds the adjacent seam edges e, f (sides matched at their shared vertex) into @p st.
void addRatioPair(const std::array<Side, 2> &e, const std::array<Side, 2> &f, RatioStats &st) {
    const double r0 = length(e[0]) / length(f[0]);
    const double r1 = length(e[1]) / length(f[1]);
    const double dot =
        (e[0].dx * f[0].dx) + (e[0].dy * f[0].dy) + (e[1].dx * f[1].dx) + (e[1].dy * f[1].dy);
    const double cosine =
        dot / (std::hypot(length(e[0]), length(e[1])) * std::hypot(length(f[0]), length(f[1])));
    const double mismatch = std::abs(r0 - r1) / std::max(r0, r1);
    const double shortest = std::min({ length(e[0]), length(e[1]), length(f[0]), length(f[1]) });
    ++st.pairs;
    st.worstAll = std::max(st.worstAll, mismatch);
    if (cosine < -(1.0 - 1e-6)) {
        ++st.parallel;
        // 1e-5 relative, plus what float UVs of this edge length can resolve.
        st.worstExcess = std::max(st.worstExcess, mismatch - (1e-5 + (4e-7 / shortest)));
    }
}


/// Edge-length ratio agreement along every unprotected, interior chain vertex of @p tm.
RatioStats seamRatios(const TexturedMesh &tm, const std::vector<std::uint8_t> &isProtected) {
    const VertexEdges edges = edgesPerVertex(tm, decimate::seam::findSeams(tm));
    RatioStats        st;
    for (const auto &[v, pairs] : edges.seams) {
        const auto copies = copiesOf(tm, v);
        const auto boundary = edges.boundary.find(v);
        const bool chainInterior =
            pairs.size() == 2 && boundary == edges.boundary.end() &&
            std::ranges::all_of(copies, [&](int u) { return isProtected[idx(u)] == 0; });
        if (!chainInterior) {
            continue;
        }
        const auto e = sidesAt(tm, *pairs[0], v);
        auto       f = sidesAt(tm, *pairs[1], v);
        if (f[0].uvAtV != e[0].uvAtV) {
            std::swap(f[0], f[1]);
        }
        if (f[0].uvAtV == e[0].uvAtV && f[1].uvAtV == e[1].uvAtV) {
            addRatioPair(e, f, st);
        }
    }
    return st;
}


/// Axis-aligned UV box of one face.
struct Box {
    float x0, y0, x1, y1;
};


/// Do the boxes share no point?
bool disjoint(const Box &a, const Box &b) {
    return a.x1 < b.x0 || b.x1 < a.x0 || a.y1 < b.y0 || b.y1 < a.y0;
}


/// The UV box of every face (for cheap pair rejection).
std::vector<Box> faceBoxes(const TexturedMesh &tm) {
    std::vector<Box> box(tm.faceCount());
    for (std::size_t f = 0; f < box.size(); ++f) {
        const Vec2f a = tm.uvs[idx(tm.faceU[3 * f])];
        const Vec2f b = tm.uvs[idx(tm.faceU[(3 * f) + 1])];
        const Vec2f c = tm.uvs[idx(tm.faceU[(3 * f) + 2])];
        box[f] = { std::min({ a.x, b.x, c.x }), std::min({ a.y, b.y, c.y }), std::max({ a.x, b.x, c.x }),
                   std::max({ a.y, b.y, c.y }) };
    }
    return box;
}


/// Number of triangles after @p f (from other charts) whose interior overlaps triangle @p f.
int overlapsAfter(const TexturedMesh &tm, const std::vector<Box> &box,
                  const std::vector<int> &faceToChart, std::size_t f) {
    const auto corner = [&](std::size_t face, std::size_t k) {
        return tm.uvs[idx(tm.faceU[(3 * face) + k])];
    };
    int overlaps = 0;
    for (std::size_t g = f + 1; g < box.size(); ++g) {
        if (faceToChart[f] == faceToChart[g] || disjoint(box[f], box[g])) {
            continue;
        }
        overlaps += trianglesOverlapSAT(corner(f, 0), corner(f, 1), corner(f, 2), corner(g, 0),
                                        corner(g, 1), corner(g, 2))
                        ? 1
                        : 0;
    }
    return overlaps;
}


/// Pairs of triangles from different charts whose interiors overlap (exact SAT).
int overlappingTrianglePairs(const TexturedMesh &tm, const std::vector<int> &faceToChart) {
    const std::vector<Box> box = faceBoxes(tm);
    int                    overlaps = 0;
    for (std::size_t f = 0; f < box.size(); ++f) {
        overlaps += overlapsAfter(tm, box, faceToChart, f);
    }
    return overlaps;
}


/// Triangles whose signed area changed sign or became zero.
int flippedTriangles(const TexturedMesh &before, const TexturedMesh &after) {
    int flipped = 0;
    for (std::size_t f = 0; f < before.faceCount(); ++f) {
        const double a = signedArea(before, f);
        const double b = signedArea(after, f);
        flipped += (a != 0.0 && (b == 0.0 || (a > 0.0) != (b > 0.0))) ? 1 : 0;
    }
    return flipped;
}


// The NOLINT on the check helpers is for Catch2: every CHECK macro expands to branching code.

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
/// Every protected UV vertex is bit-identical after straightening, and the count matches the stats.
void checkProtectedBitIdentical(const TexturedMesh &before, const StraightenResult &res) {
    int protectedChecked = 0;
    for (std::size_t u = 0; u < before.uvs.size(); ++u) {
        if (res.protectedUV[u] != 0) {
            ++protectedChecked;
            CHECK(sameBits(before.uvs[u], res.mesh.uvs[u]));
        }
    }
    CHECK(protectedChecked == res.stats.protectedVertices);
}


// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
/// Bakes with the straightened UVs and checks that seam erasure still removes the seam energy.
void checkErasure(const Pipeline &p, const TexturedMesh &straightened) {
    const auto lowN =
        decimate::bake::computeVertexNormals(p.low, decimate::bake::Weighting::EqualUnit);
    const auto highN =
        decimate::bake::computeVertexNormals(p.high, decimate::bake::Weighting::EqualUnit);
    const NormalMap map =
        decimate::bake::bake(p.low, straightened.cornerUV(), lowN, p.high, highN, kTex, 2, false);
    const auto erased =
        decimate::seam::erase(straightened, map, decimate::seam::NormalSpace::Object);
    // 1e-7, not the spec's 1e-8: at the paper's w_seam = 1e10 the seam penalty leaves
    // D_after at ~3e-8 x D_before on these bakes (the same limit 17-01's I3 documents).
    for (std::size_t c = 0; c < 3; ++c) {
        CHECK(erased.dBefore[c] > 0.0);
        CHECK(erased.dAfter[c] <= 1e-7 * erased.dBefore[c]);
    }
}


/// What a pipeline run measured, for the printed summary.
struct RunReport {
    const char *name;
    RatioStats  before;
    RatioStats  after;
    int         overlapBefore;
    int         overlapAfter;
};


/// Prints one summary line per pipeline run (for the ctest log, not asserted).
void printRun(const Pipeline &p, const StraightenResult &res, const RunReport &r) {
    const auto &st = res.stats;
    std::printf("[straighten] %s: chains seam=%d boundary=%d protected=%d (rules %d/%d/%d/%d) "
                "ears=%d (same-neighbour %d) "
                "reverted=%d maxMove=%.5f ratio-pairs=%d parallel %d->%d worst-mismatch %.3e->%.3e "
                "overlap %d->%d\n",
                r.name, st.seamChains, st.boundaryChains, st.protectedVertices,
                st.protectedByRule[0], st.protectedByRule[1], st.protectedByRule[2],
                st.protectedByRule[3], st.earFaces, st.sameNeighbourEars, st.revertedCharts,
                maxMove(p.tm, res.mesh), r.after.pairs, r.before.parallel, r.after.parallel,
                r.before.worstAll, r.after.worstAll, r.overlapBefore, r.overlapAfter);
}


/// The I1/I2 contract on one pipeline; returns the straightened result for further checks.
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
StraightenResult checkStraightenedAtlas(const Pipeline &p, const char *name) {
    const StraightenResult res = decimate::seam::straighten(p.tm, StraightenOptions { .tau = 0.01 });

    // Topology and protected vertices are untouched; no triangle flips.
    CHECK(res.stats.seamChains > 0);
    CHECK(res.mesh.faceV == p.tm.faceV);
    CHECK(res.mesh.faceU == p.tm.faceU);
    checkProtectedBitIdentical(p.tm, res);
    CHECK(flippedTriangles(p.tm, res.mesh) == 0);

    // Charts stay separated (exact SAT) -- the input is the reference.
    const int overlapBefore = overlappingTrianglePairs(p.tm, p.charts.faceToChart);
    const int overlapAfter = overlappingTrianglePairs(res.mesh, p.charts.faceToChart);
    CHECK(overlapBefore == 0);
    CHECK(overlapAfter == 0);

    // Edge-length ratios of straightened seam runs agree across the two copies.
    const RatioStats before = seamRatios(p.tm, res.protectedUV);
    const RatioStats after = seamRatios(res.mesh, res.protectedUV);
    CHECK(after.parallel > before.parallel);
    CHECK(after.worstExcess <= 0.0);

    // Erasure on top of the straightened atlas still removes the seam energy.
    checkErasure(p, res.mesh);
    printRun(p, res, { name, before, after, overlapBefore, overlapAfter });
    return res;
}

} // namespace

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S1: an already straight, ratio-matched chain does not move", "[seam][straighten][S1]") {
    const Fixture          fx = stripFixture(0.0f);
    const StraightenResult res = decimate::seam::straighten(fx.tm, StraightenOptions { .tau = 0.01 });

    CHECK(res.stats.earFaces == 0); // fixture precondition: the whole chain is straightenable
    CHECK(res.stats.seamChains == 1);
    CHECK(res.stats.revertedCharts == 0);
    CHECK(maxMove(fx.tm, res.mesh) <= 1e-6);
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S2: a noisy chain is straightened and its edge-length ratios match",
          "[seam][straighten][S2]") {
    const Fixture          fx = stripFixture(0.002f);
    const StraightenResult res = decimate::seam::straighten(fx.tm, StraightenOptions { .tau = 0.01 });
    CHECK(res.stats.earFaces == 0);
    CHECK(res.stats.seamChains == 1);
    // The input chain really is bent: the interior copies sit 0.002 off the end-to-end line.
    CHECK(maxMove(fx.tm, res.mesh) >= 0.0019);

    // Collect each copy's chain, ordered by u, to compare the two sides.
    std::array<std::array<Vec2f, 5>, 2> side {}; // [copy][j]
    for (int j = 0; j <= 4; ++j) {
        const std::vector<int> copies = copiesOf(res.mesh, stripVertex(j));
        REQUIRE(copies.size() == 2);
        side[0][idx(j)] = res.mesh.uvs[idx(copies[0])];
        side[1][idx(j)] = res.mesh.uvs[idx(copies[1])];
    }

    // The interior vertices must lie exactly on the end-to-end line of their copy.
    const std::array<float, 2> lineX { 0.25f, 0.55f }; // end-to-end segment of each copy
    for (std::size_t s = 0; s < 2; ++s) {
        for (std::size_t j = 1; j <= 3; ++j) {
            CHECK(std::abs(side[s][j].x - lineX[s]) <= 1e-6);
        }
    }

    // Equal edge-length ratios on both copies, so the seam samples the same surface points from both sides.
    for (std::size_t j = 1; j <= 3; ++j) {
        const auto len = [&](std::size_t s, std::size_t a, std::size_t b) {
            return std::hypot(static_cast<double>(side[s][a].x) - side[s][b].x,
                              static_cast<double>(side[s][a].y) - side[s][b].y);
        };
        const double ratio0 = len(0, j - 1, j) / len(0, j, j + 1);
        const double ratio1 = len(1, j - 1, j) / len(1, j, j + 1);
        INFO("vertex " << j << " ratios " << ratio0 << " vs " << ratio1);
        CHECK(std::abs(ratio0 - ratio1) <= 1e-6);
    }
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S3: three charts meeting at one vertex protect it (rule 1)", "[seam][straighten][S3]") {
    // Pie fixture: three charts of equal angle around the grid centre.
    constexpr int kSide = 6;
    const double  twoPi = 2.0 * std::numbers::pi;
    const Fixture fx =
        gridFixture(kSide, kSide,
                    [&](int, double cx, double cy) {
                        double angle = std::atan2(cy - 3.0, cx - 3.0);
                        angle = angle < 0.0 ? angle + twoPi : angle;
                        return static_cast<int>(angle / (twoPi / 3.0));
                    },
                    { { 0.0f, 0.0f }, { 0.35f, 0.0f }, { 0.0f, 0.35f } }, { { kSide - 1, 0 }, { 0, kSide - 1 } });
    const int center = (3 * (kSide + 1)) + 3;
    REQUIRE(copiesOf(fx.tm, center).size() == 3);

    const StraightenResult res = decimate::seam::straighten(fx.tm, StraightenOptions { .tau = 0.01 });
    CHECK(res.stats.protectedByRule[0] >= 3);
    CHECK(copiesUnchanged(fx.tm, res.mesh, center));
    for (const int u : copiesOf(fx.tm, center)) {
        CHECK(res.protectedUV[idx(u)] != 0);
    }
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S4: a seam ending on a 3D boundary keeps the junction fixed (rule 2)",
          "[seam][straighten][S4]") {
    const Fixture          fx = stripFixture(0.002f);
    const StraightenResult res = decimate::seam::straighten(fx.tm, StraightenOptions { .tau = 0.01 });
    CHECK(res.stats.protectedByRule[1] >= 2);

    // Both chain ends (on the 3D boundary) stay put ...
    for (const int j : { 0, 4 }) {
        CHECK(copiesOf(fx.tm, stripVertex(j)).size() == 2);
        CHECK(copiesUnchanged(fx.tm, res.mesh, stripVertex(j)));
    }
    // ... while the interior of the chain did move.
    CHECK_FALSE(copiesUnchanged(fx.tm, res.mesh, stripVertex(2)));
}


// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S5: a foldover reverts the chart and its protection propagates",
          "[seam][straighten][S5]") {
    const Fixture fx = foldFixture(1.15f);

    // Control: a tolerance below the spike height keeps the tip as a corner, nothing folds.
    const StraightenResult keep = decimate::seam::straighten(fx.tm, StraightenOptions { .tau = 0.01 });
    CHECK(keep.stats.revertedCharts == 0);
    CHECK(keep.stats.seamChains == 1);

    // tau = 0.5 flattens the spike, which drags its tip through the interior vertex.
    const StraightenResult res = decimate::seam::straighten(fx.tm, StraightenOptions { .tau = 0.5 });
    CHECK(res.stats.revertedCharts == 1);
    REQUIRE(res.mesh.uvs.size() == fx.tm.uvs.size());

    // A reverted chart restores its UVs bit for bit; the other chart follows via rule 4.
    for (std::size_t u = 0; u < fx.tm.uvs.size(); ++u) {
        CHECK(sameBits(fx.tm.uvs[u], res.mesh.uvs[u])); // chart A, and chart B via rule 4
    }
    for (int sv = 0; sv <= 7; ++sv) { // chart A's vertices, every copy (the chain's included)
        for (const int u : copiesOf(fx.tm, sv)) {
            CHECK(res.protectedUV[idx(u)] != 0);
        }
    }
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S6: a closed boundary loop never coarsens below three edges", "[seam][straighten][S6]") {
    // A fan of 8 triangles around a centre vertex: the rim is one closed boundary loop.
    constexpr int kRim = 8;
    MeshData      disc;
    disc.positions.emplace_back(0.0f, 0.0f, 0.0f);
    std::vector<Vec2f> corner;
    const auto         rim = [&](int k) {
        const double a = 2.0 * std::numbers::pi * k / kRim;
        return Vec2f { static_cast<float>(0.5 + (0.3 * std::cos(a))),
                       static_cast<float>(0.5 + (0.3 * std::sin(a))) };
    };
    for (int k = 0; k < kRim; ++k) {
        const Vec2f p = rim(k);
        disc.positions.emplace_back((p.x - 0.5f) / 0.3f, (p.y - 0.5f) / 0.3f, 0.0f);
    }
    for (int k = 0; k < kRim; ++k) {
        for (const int v : { 0, 1 + k, 1 + ((k + 1) % kRim) }) {
            disc.indices.push_back(v);
            corner.push_back(v == 0 ? Vec2f { 0.5f, 0.5f } : rim(v - 1));
        }
    }

    // tau = 1.0 is loose enough to coarsen the whole rim, which the three-edge floor must stop.
    const TexturedMesh     tm = decimate::seam::weld(disc, corner);
    const StraightenResult res = decimate::seam::straighten(tm, StraightenOptions { .tau = 1.0 });
    CHECK(res.stats.boundaryChains == 1);
    CHECK(res.stats.revertedCharts == 0);

    int corners = 0; // rim vertices where the loop turns
    for (int k = 1; k <= kRim; ++k) {
        const Vec2f  prev = res.mesh.uvs[idx(copiesOf(res.mesh, 1 + ((k + kRim - 2) % kRim))[0])];
        const Vec2f  cur = res.mesh.uvs[idx(copiesOf(res.mesh, k)[0])];
        const Vec2f  next = res.mesh.uvs[idx(copiesOf(res.mesh, 1 + (k % kRim))[0])];
        const double cross =
            ((static_cast<double>(cur.x) - prev.x) * (static_cast<double>(next.y) - cur.y)) -
            ((static_cast<double>(cur.y) - prev.y) * (static_cast<double>(next.x) - cur.x));
        corners += std::abs(cross) > 1e-6 ? 1 : 0;
    }
    CHECK(corners >= 3);
    CHECK(maxMove(tm, res.mesh) > 0.01); // it really coarsened

    for (std::size_t f = 0; f < res.mesh.faceCount(); ++f) {
        CHECK(signedArea(res.mesh, f) > 1e-6);
    }
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S7: ears are counted, classified and protected (rule 3)", "[seam][straighten][S7]") {
    // Faces 14 and 17 form chart 1 inside chart 0, face 16 is chart 2.  Face 14's two
    // other edges both face chart 0; face 17's face chart 0 and chart 2.
    const Fixture          fx = gridFixture(6, 6,
                                            [](int face, double, double) {
                                       return face == 14 || face == 17 ? 1 : (face == 16 ? 2 : 0);
                                            },
                                            { { 0.0f, 0.0f }, { 0.4f, 0.0f }, { 0.4f, 0.3f } }, { { 5, 0 }, { 0, 5 } });
    const StraightenResult res = decimate::seam::straighten(fx.tm, StraightenOptions { .tau = 0.01 });
    CHECK(res.stats.earFaces == 2);
    CHECK(res.stats.sameNeighbourEars == 1);
    CHECK(res.stats.protectedByRule[2] > 0);

    // Every corner of an ear is protected.
    for (const std::size_t f : { 14U, 17U }) {
        for (std::size_t k = 0; k < 3; ++k) {
            CHECK(res.protectedUV[idx(fx.tm.faceU[(3 * f) + k])] != 0);
        }
    }
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S8: straighten never touches faces or UV indexing", "[seam][straighten][S8]") {
    for (const float noise : { 0.0f, 0.002f }) {
        const Fixture          fx = stripFixture(noise);
        const StraightenResult res =
            decimate::seam::straighten(fx.tm, StraightenOptions { .tau = 0.01 });
        CHECK(res.mesh.faceV == fx.tm.faceV);
        CHECK(res.mesh.faceU == fx.tm.faceU);
        CHECK(res.mesh.uvs.size() == fx.tm.uvs.size());
        CHECK(res.mesh.positions == fx.tm.positions);
    }
}


// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("I1: straighten -> bake -> erase on the closed spot atlas", "[seam][straighten][I1]") {
    checkStraightenedAtlas(spot(), "spot");
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("I2: straighten -> bake -> erase on the open bunny atlas", "[seam][straighten][I2]") {
    const Pipeline         bunny = buildPipeline(modelOrSkip("bunny.iv"));
    const StraightenResult res = checkStraightenedAtlas(bunny, "bunny");
    CHECK(res.stats.boundaryChains > 0);
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B1: the default threshold is byte-stable", "[seam][straighten][B1]") {
    const Pipeline &p = spot();
    // The cutter's own output: p.charts has since been repaired in place by the 19-03 injectivity
    // guard (parametrizeCorners), and the guard is deliberately not part of segment().
    const auto cut = decimate::chart::segment(p.hem, p.low, 64);
    const auto explicit03 = decimate::chart::segment(p.hem, p.low, 64, 0.3);
    CHECK(explicit03.faceToChart == cut.faceToChart);
    CHECK(chartSummary(explicit03) == chartSummary(cut));
    CHECK(decimate::chart::kStraightenThreshold == 0.3);

    // Fingerprints re-captured by Spec 19-01 (mode CentroidMode::TrueCentroid, the default) from the
    // restored 2003 partition (empty seed history, float edge-cost denominator, no final regrow).
    // The previous values pinned the pre-19-01 cutter (cap never bound).
    CHECK(fingerprint(cut.faceToChart) == 6303312352358204886ULL);
    CHECK(fingerprint(chartSummary(cut)) == 1573457813196365394ULL);
    const Pipeline bunny = buildPipeline(modelOrSkip("bunny.iv"));
    const auto     bunnyCut = decimate::chart::segment(bunny.hem, bunny.low, 64);
    CHECK(fingerprint(bunnyCut.faceToChart) == 16017822343820123847ULL);
    CHECK(fingerprint(chartSummary(bunnyCut)) == 17478674128315925162ULL);
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B2: the threshold reaches the cutter's comparison", "[seam][straighten][B2]") {
    const Pipeline &p = spot();

    // Threshold 0 never folds a face, 2 always does: the two partitions must differ.
    const auto never = decimate::chart::segment(p.hem, p.low, 64, 0.0);
    const auto always = decimate::chart::segment(p.hem, p.low, 64, 2.0);
    int        differing = 0;
    for (std::size_t f = 0; f < never.faceToChart.size(); ++f) {
        differing += never.faceToChart[f] != always.faceToChart[f] ? 1 : 0;
    }
    std::printf("[segment] spot: threshold 0.0 vs 2.0 differ on %d of %zu faces\n", differing,
                never.faceToChart.size());
    CHECK(differing >= 1);
}
