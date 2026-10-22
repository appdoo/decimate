// This file is part of the decimate library; see decimate.h for version/license details.
/// @file self_seam_cut_test.cpp
/// @brief Spec 19-02: the 2003 self-seam cut - intra-chart seams are exported by the segmenter
///        and `flatten()` builds the chart's local mesh by corner, cut along them.
///
/// B1 [B]  bunny/16 base plate: one boundary loop, no flipped triangle
/// B2 [B]  open tube, one chart: cut into a disk (1 loop instead of 2), no flip
/// B3 [B]  bunny/16 parametrizeCorners: every corner UV is its local corner's UV
/// S1 [S]  cow/24 (Legacy2003Middle): 4 seam half-edges (2003 keeps 3 edges = 6)
/// S2 [S]  seamHalfEdges contract on segment() results (both faces in the chart, sym listed, sorted)
/// S3 [S]  spot/64 has no intra-chart seam, so the cut is a no-op (UV golden: halfedge_chart_test B1)
/// S4 [S]  seamHalfEdges contract holds for every chart after the injectivity guard split charts

#include "model_fixture.h"

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>
#include <set>
#include <string>
#include <vector>

using namespace decimate;
using namespace decimate::chart;

namespace {

/// A mesh plus its half-edge mesh.
struct Loaded {
    MeshData     data;
    HalfEdgeMesh hem;
};


/// Builds the half-edge mesh over @p d.
Loaded fromData(MeshData d) {
    Loaded l;
    l.data = std::move(d);
    l.hem = HalfEdgeMesh::create(static_cast<int>(l.data.positions.size()),
                                 static_cast<int>(l.data.indices.size() / 3),
                                 reinterpret_cast<const Face *>(l.data.indices.data()));
    return l;
}


/// Loads a model file and builds its half-edge mesh.
Loaded loadMesh(const std::string &path) {
    return fromData(load(path));
}


/// Number of boundary loops of a local mesh: a walk over firstBoundaryEdge/nextBoundaryEdge
/// (`boundaryPointIndices()` returns one loop only).
int boundaryLoops(const HalfEdgeMesh &m) {
    std::vector<char> seen(m.halfEdgeCount(), 0);
    int               loops = 0;
    for (HeIndex h = 0; h < m.halfEdgeCount(); ++h) {
        if (m.sym(h) != kNoHalfedge || seen[h] != 0) {
            continue;
        }
        ++loops;
        HeIndex x = h;
        while (x != kNoHalfedge && seen[x] == 0) {
            seen[x] = 1;
            x = m.nextBoundaryEdge(x);
        }
    }
    return loops;
}


/// Local triangles whose UV signed area has the sign opposite to the chart's majority.
int flippedTriangles(const HalfEdgeChart &ac) {
    int positive = 0;
    int negative = 0;
    for (int i = 0; i < static_cast<int>(ac.faces.size()); ++i) {
        const auto &a = ac.uv[static_cast<std::size_t>(localCorner(ac, i, 0))];
        const auto &b = ac.uv[static_cast<std::size_t>(localCorner(ac, i, 1))];
        const auto &c = ac.uv[static_cast<std::size_t>(localCorner(ac, i, 2))];
        const float area = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
        positive += area > 0.0f ? 1 : 0;
        negative += area < 0.0f ? 1 : 0;
    }
    return std::min(positive, negative);
}


/// True when a mesh vertex occurs as more than one local vertex (the chart was cut through it).
bool hasSplitVertex(const HalfEdgeChart &ac) {
    const std::set<int> distinct(ac.localToMesh.begin(), ac.localToMesh.end());
    return distinct.size() < ac.localToMesh.size();
}


/// The seamHalfEdges contract (S2): sorted, both faces in the chart, sym listed too.
// NOLINTNEXTLINE(readability-function-size): one predicate per clause of the contract
bool seamContractHolds(const HalfEdgeMesh &hem, const ChartSet &cs) {
    for (const Chart &c : cs.charts) {
        if (!std::ranges::is_sorted(c.seamHalfEdges)) {
            return false;
        }
        for (const int h : c.seamHalfEdges) {
            const HeIndex he = static_cast<HeIndex>(h);
            const HeIndex s = hem.sym(he);
            if (s == kNoHalfedge) {
                return false;
            }
            const bool inChart = cs.faceToChart[hem.halfedge(he).face] == c.id &&
                                 cs.faceToChart[hem.halfedge(s).face] == c.id;
            const bool symListed = std::ranges::binary_search(c.seamHalfEdges, static_cast<int>(s));
            if (!inChart || !symListed) {
                return false;
            }
        }
    }
    return true;
}


/// Total number of seam half-edges over all charts.
std::size_t totalSeams(const ChartSet &cs) {
    std::size_t n = 0;
    for (const Chart &c : cs.charts) {
        n += c.seamHalfEdges.size();
    }
    return n;
}


/// Open tube: 16 segments x 4 rings of quads (5 vertex rings), radius 1, height 2.
MeshData makeTube() {
    constexpr int kSeg = 16;
    constexpr int kRings = 5;
    MeshData      d;
    for (int r = 0; r < kRings; ++r) {
        for (int s = 0; s < kSeg; ++s) {
            const float a = 6.2831853f * static_cast<float>(s) / kSeg;
            d.positions.emplace_back(std::cos(a), std::sin(a), 2.0f * static_cast<float>(r) / (kRings - 1));
        }
    }
    for (int r = 0; r + 1 < kRings; ++r) {
        for (int s = 0; s < kSeg; ++s) {
            const int s1 = (s + 1) % kSeg;
            const int a = r * kSeg + s;
            const int b = r * kSeg + s1;
            const int c = (r + 1) * kSeg + s;
            const int e = (r + 1) * kSeg + s1;
            d.indices.insert(d.indices.end(), { a, b, e, a, e, c });
        }
    }
    return d;
}

} // namespace

/// B1 - the bunny base plate is a disk after the cut
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B1: bunny/16 base plate flattens to one boundary loop and no flipped triangle",
          "[chart][selfseam][B1]") {
    const auto     l = loadMesh(modelOrSkip("bunny.iv"));
    const ChartSet cs = segment(l.hem, l.data, 16);

    // The base plate owns the face with the lowest centroid y.
    // (Sum of the corner y values; the 3x factor does not change the arg-min.)
    int   lowest = 0;
    float lowestY = 1e30f;
    for (std::size_t f = 0; f < l.data.indices.size() / 3; ++f) {
        float y = 0.0f;
        for (std::size_t k = 0; k < 3; ++k) {
            y += l.data.positions[static_cast<std::size_t>(l.data.indices[3 * f + k])].y;
        }
        if (y < lowestY) {
            lowestY = y;
            lowest = static_cast<int>(f);
        }
    }

    const Chart        &plate = cs.charts[static_cast<std::size_t>(cs.faceToChart[static_cast<std::size_t>(lowest)])];
    const HalfEdgeChart ac = flatten(l.hem, l.data, plate);
    INFO("base plate chart " << plate.id << ": faces " << plate.faces.size() << ", seam half-edges "
                             << plate.seamHalfEdges.size() << ", local vertices " << ac.localToMesh.size());
    CHECK(boundaryLoops(ac.mesh) == 1);
    CHECK(flippedTriangles(ac) == 0);
    CHECK(hasSplitVertex(ac));
}

/// B2 - an open tube in one chart is cut open into a disk
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B2: an open tube in one chart is cut into a disk", "[chart][selfseam][B2]") {
    const auto l = fromData(makeTube());
    REQUIRE(l.data.positions.size() == 80);
    REQUIRE(l.data.indices.size() / 3 == 128);

    // Cap 1 forces the whole tube into a single chart.
    const ChartSet cs = segment(l.hem, l.data, 1);
    REQUIRE(cs.charts.size() == 1);
    const Chart &chart = cs.charts[0];
    CHECK_FALSE(chart.seamHalfEdges.empty());

    // Cutting along the exported seams turns the tube into a disk.
    const HalfEdgeChart ac = flatten(l.hem, l.data, chart);
    CHECK(boundaryLoops(ac.mesh) == 1); // an uncut tube has 2 (top and bottom ring)
    CHECK(ac.localToMesh.size() > 80);
    CHECK(flippedTriangles(ac) == 0);
}

/// B3 - corner UVs come from the local triangulation
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B3: bunny/16 parametrizeCorners reads every corner through its local triangulation",
          "[chart][selfseam][B3]") {
    const auto l = loadMesh(modelOrSkip("bunny.iv"));
    ChartSet   viaCorners = segment(l.hem, l.data, 16);
    const auto cornerUV = parametrizeCorners(l.hem, l.data, viaCorners, 256);

    // The reference: the same partition flattened and packed chart by chart.
    ChartSet                   again = segment(l.hem, l.data, 16);
    std::vector<HalfEdgeChart> atlas = flattenCharts(l.hem, l.data, again);
    packAtlas(atlas, l.data, 256);
    REQUIRE(again.charts.size() == viaCorners.charts.size());

    // Every corner of every face must carry the UV of its local (possibly split) vertex.
    bool anySplit = false;
    int  compared = 0;
    int  mismatches = 0;
    for (const HalfEdgeChart &ac : atlas) {
        anySplit = anySplit || hasSplitVertex(ac);
        for (int i = 0; i < static_cast<int>(ac.faces.size()); ++i) {
            for (int k = 0; k < 3; ++k) {
                const auto &want = ac.uv[static_cast<std::size_t>(localCorner(ac, i, k))];
                const auto  face = static_cast<std::size_t>(ac.faces[static_cast<std::size_t>(i)]);
                const auto &got = cornerUV[3 * face + static_cast<std::size_t>(k)];
                mismatches += (want.x != got.x || want.y != got.y) ? 1 : 0;
                ++compared;
            }
        }
    }

    CHECK(anySplit); // otherwise the lookup-by-vertex-id shortcut would also pass
    CHECK(compared == static_cast<int>(l.data.indices.size()));
    CHECK(mismatches == 0);
}

/// S1 - the stale-edge fix: cow/24 keeps 2 intra-chart seam edges where 2003 keeps 3
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S1: cow/24 (Legacy2003Middle) exports 4 seam half-edges", "[chart][selfseam][S1]") {
    const auto     l = loadMesh(modelOrSkip("cow.obj"));
    const ChartSet cs =
        segment(l.hem, l.data, 24, kStraightenThreshold, CentroidMode::Legacy2003Middle);
    CHECK(totalSeams(cs) == 4);
}

/// S2 - the export contract
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S2: seamHalfEdges lists both halves of intra-chart edges, sorted", "[chart][selfseam][S2]") {
    for (const char *model : { "cow.obj", "bunny.iv" }) {
        const auto l = loadMesh(modelOrSkip(model));
        for (const int cap : { 16, 64 }) {
            CAPTURE(model, cap);
            const ChartSet cs = segment(l.hem, l.data, cap);
            CHECK(seamContractHolds(l.hem, cs));
        }
    }
    const auto spot = loadMesh(std::string(MESH_MODEL_DIR) + "/spot.obj");
    CHECK(seamContractHolds(spot.hem, segment(spot.hem, spot.data, 64)));
}

/// S3 - no seams, no cut
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S3: spot/64 has no intra-chart seam, so the cut is a no-op", "[chart][selfseam][S3]") {
    const auto l = loadMesh(std::string(MESH_MODEL_DIR) + "/spot.obj");
    ChartSet   cs = segment(l.hem, l.data, 64);
    CHECK(totalSeams(cs) == 0);

    const auto atlas = flattenCharts(l.hem, l.data, cs);
    for (const HalfEdgeChart &ac : atlas) {
        CHECK_FALSE(hasSplitVertex(ac));
    }
    CHECK(totalSeams(cs) == 0);
}

/// S4 (19-03 A1) - split pieces keep only the seams wholly inside them
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S4: after the guard every chart's seamHalfEdges satisfies S2 and none separates two charts",
          "[chart][selfseam][S4]") {
    for (const char *model : { "spot.obj", "cow.obj", "bunny.iv" }) {
        const std::string path = std::string(model) == "spot.obj"
                                     ? std::string(MESH_MODEL_DIR) + "/spot.obj"
                                     : modelOrSkip(model);
        const auto        l = loadMesh(path);
        ChartSet          cs = segment(l.hem, l.data, 64);
        InjectivityStats  st;
        const auto        atlas = flattenCharts(l.hem, l.data, cs, {}, &st);
        CAPTURE(model, st.splits);
        CHECK(seamContractHolds(l.hem, cs));
        CHECK(st.selfOverlappingAfter == 0);
    }
}
