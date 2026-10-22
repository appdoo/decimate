// This file is part of the decimate library; see decimate.h for version/license details.
/// @file correctness_chart_test.cpp
/// @brief Spec 24-01 — hangs and stack depth in the chart module (H1-H5).
///
/// H1 [S]  The straighten fixpoint loop is bounded: a step that never settles stops at the cap.
/// H2 [S]  BoundaryFlags::erase peels a 200 000-edge stub chain without recursing per edge.
/// H3 [S]  nextBoundaryEdge on a corrupt (cyclic) fan returns kNoHalfedge instead of spinning.
/// H4 [B]  packAtlas rejects a boundary-less (closed) chart before the ratio loop.
/// H5 [S]  packAtlas rejects atlasResolution < 2.

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>

#include "segmenter_detail.h"

using namespace decimate;
using namespace decimate::chart;


namespace {

/// Triangle strip of `n` quads: top row t_0..t_n (ids 0..n), bottom row b_0..b_n (ids n+1..2n+1).
/// The top edges t_i - t_{i+1} are the mesh boundary; everything else is shared.
HalfEdgeMesh makeStrip(int n) {
    std::vector<int> idx;
    idx.reserve(static_cast<std::size_t>(n) * 6);
    for (int i = 0; i < n; ++i) {
        const int t0 = i;
        const int t1 = i + 1;
        const int b0 = n + 1 + i;
        const int b1 = n + 2 + i;
        idx.insert(idx.end(), { t0, b0, t1, b0, b1, t1 });
    }
    return HalfEdgeMesh::create(2 * n + 2, 2 * n, reinterpret_cast<const Face *>(idx.data()));
}


/// Closed tetrahedron as one chart: every edge is shared, so the chart has no boundary loop.
HalfEdgeChart makeTetraChart(MeshData &data) {
    data.positions = { { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, { 0.0f, 0.0f, 1.0f } };
    data.indices = { 0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3 };

    HalfEdgeChart ac;
    ac.mesh = HalfEdgeMesh::create(4, 4, reinterpret_cast<const Face *>(data.indices.data()));
    ac.uv = { { 0.0f, 0.0f }, { 1.0f, 0.0f }, { 0.0f, 1.0f }, { 1.0f, 1.0f } };
    ac.localToMesh = { 0, 1, 2, 3 };
    ac.faces = { 0, 1, 2, 3 };
    return ac;
}


/// One flat triangle chart: a valid input whose only defect can be the resolution.
HalfEdgeChart makeTriangleChart(MeshData &data) {
    data.positions = { { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f } };
    data.indices = { 0, 1, 2 };

    HalfEdgeChart ac;
    ac.mesh = HalfEdgeMesh::create(3, 1, reinterpret_cast<const Face *>(data.indices.data()));
    ac.uv = { { 0.0f, 0.0f }, { 1.0f, 0.0f }, { 0.0f, 1.0f } };
    ac.localToMesh = { 0, 1, 2 };
    ac.faces = { 0 };
    return ac;
}

} // namespace


TEST_CASE("H1: runBounded stops a step that never settles at the cap", "[chart][correctness][H1]") {
    int  calls = 0;
    auto pingPong = [&calls] {
        ++calls;
        return true; // two faces trading places forever: always "changed"
    };

    const int passes = detail::runBounded(pingPong, 64);

    CHECK(passes == 64);
    CHECK(calls == 64);
}


TEST_CASE("H1: runBounded returns at the fixpoint, before the cap", "[chart][correctness][H1]") {
    int  calls = 0;
    auto settles = [&calls] {
        return ++calls < 3;
    }; // true, true, false

    CHECK(detail::runBounded(settles, 64) == 3);
    CHECK(calls == 3);
}


TEST_CASE("H1: segment with straightenThreshold 2.0 terminates and covers every face",
          "[chart][correctness][H1]") {
    // Cube, 12 triangles: threshold 2.0 lets straighten() move any face with two foreign edges.
    const std::vector<int> idx = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4,
                                   2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5 };
    MeshData               data;
    for (int i = 0; i < 8; ++i) {
        data.positions.push_back({ float(i & 1), float((i >> 1) & 1), float((i >> 2) & 1) });
    }
    data.indices = idx;
    const HalfEdgeMesh mesh = HalfEdgeMesh::create(8, 12, reinterpret_cast<const Face *>(idx.data()));

    const ChartSet set = segment(mesh, data, 4, 2.0);

    REQUIRE(set.faceToChart.size() == 12);
    for (const int c : set.faceToChart) {
        CHECK(c >= 0);
        CHECK(c < static_cast<int>(set.charts.size()));
    }
}


TEST_CASE("H2: erase peels a 200000-edge stub chain iteratively", "[chart][correctness][H2]") {
    constexpr int kChain = 200000;
    HalfEdgeMesh  mesh = makeStrip(kChain);

    detail::BoundaryFlags bf;
    bf.init(mesh);
    bf.reset();

    // Leave only the top boundary path flagged: a chain of kChain edges whose inner vertices have degree 2.
    const std::uint32_t topCount = kChain + 1;
    HeIndex             first = kNoHalfedge;
    for (HeIndex h = 0; h < mesh.halfEdgeCount(); ++h) {
        const std::uint32_t endV = mesh.halfedge(h).end;
        const std::uint32_t startV = mesh.halfedge(mesh.next(mesh.next(h))).end;
        if (endV < topCount && startV < topCount) {
            if (startV + 1 == endV && startV == 0) {
                first = h; // t_0 -> t_1 : the chain's tip
            }
            continue;
        }
        bf.setFlag(h, false);
    }
    REQUIRE(first != kNoHalfedge);
    REQUIRE(bf.isBoundary(first));

    bf.erase(first);

    // One erase dissolved the whole chain: every edge unflagged, every degree back to zero.
    std::size_t flagged = 0;
    for (const bool f : bf.flags) {
        flagged += f ? 1 : 0;
    }
    std::size_t degreeSum = 0;
    for (const int d : bf.degrees) {
        degreeSum += static_cast<std::size_t>(d);
    }
    CHECK(flagged == 0);
    CHECK(degreeSum == 0);
}


TEST_CASE("H3: nextBoundaryEdge on a cyclic non-manifold fan returns the sentinel",
          "[chart][correctness][H3]") {
    // Fan of four triangles around vertex 0; outer edges 1-2, 2-3, 3-4, 4-1 are the boundary.
    const std::vector<int> idx = { 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 1 };
    HalfEdgeMesh           mesh = HalfEdgeMesh::create(5, 4, reinterpret_cast<const Face *>(idx.data()));

    // The honest walk: from the first boundary edge, rotate once, and stop at the next unpaired edge.
    const HeIndex h = mesh.firstBoundaryEdge();
    REQUIRE(h != kNoHalfedge);
    const HeIndex n0 = mesh.next(h);
    const HeIndex n1 = mesh.next(mesh.sym(n0));
    REQUIRE(mesh.sym(n1) == kNoHalfedge);
    REQUIRE(mesh.nextBoundaryEdge(h) == n1);

    // Corrupt the fan: n1 gets the twin of n0, so the rotation n0 -> next(sym) = n1 -> next(sym) = n1
    // revisits n1 forever and never returns to h.
    mesh.halfedge(n1).sym = mesh.sym(n0);

    CHECK(mesh.nextBoundaryEdge(h) == kNoHalfedge);
}


TEST_CASE("H4: packAtlas rejects a boundary-less chart and names it", "[chart][correctness][H4]") {
    MeshData                   data;
    std::vector<HalfEdgeChart> charts;
    charts.push_back(makeTetraChart(data));

    REQUIRE(charts[0].mesh.firstBoundaryEdge() == kNoHalfedge);
    REQUIRE_THROWS_AS(packAtlas(charts, data, 256), std::invalid_argument);
    try {
        packAtlas(charts, data, 256);
    } catch (const std::invalid_argument &e) {
        const std::string what = e.what();
        CHECK(what.find("packAtlas") != std::string::npos);
        CHECK(what.find("chart 0") != std::string::npos);
    }
}


TEST_CASE("H5: packAtlas rejects atlasResolution below 2", "[chart][correctness][H5]") {
    MeshData                   data;
    std::vector<HalfEdgeChart> charts;
    charts.push_back(makeTriangleChart(data));

    for (const int res : { 1, 0, -7 }) {
        CAPTURE(res);
        CHECK_THROWS_AS(packAtlas(charts, data, res), std::invalid_argument);
    }
    CHECK_NOTHROW(packAtlas(charts, data, 2));
}
