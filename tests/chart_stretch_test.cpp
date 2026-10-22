// This file is part of the decimate library; see decimate.h for version/license details.
/// @file chart_stretch_test.cpp
/// @brief Spec 20: closed charts are opened (Gu 2002, 3.1) and over-stretched charts are cut open.
///
/// a [S]  bound off: charts, seams and UV are those of the plain `flatten`, no cuts
/// b [B]  closed cube: opened to a disk, injective, packs, no overlap
/// c [B]  synthetic torus (genus 1): opened to a disk, injective, packs, no overlap
/// d [I]  cow as one chart, bound on: cut open, bound reported honestly, injective, no overlap
/// e [S]  the cut path: from a corner of the worst face to the boundary, interior edges only, both halves seamed
/// f [S]  two runs are bit-identical
/// [.chartcut] hidden: charts, cuts, stretch, density and time of cow, bunny, horse (start 1).

#include "model_fixture.h"
#include "sat_overlap.h"
#include "chart_cut.h"

#include <algorithm>
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <numbers>
#include <string>
#include <vector>

using namespace decimate;
using namespace decimate::chart;

namespace {

/// A mesh and its half-edge structure.
struct Loaded {
    MeshData     data;
    HalfEdgeMesh hem;
};


Loaded makeLoaded(MeshData data) {
    Loaded l;
    l.data = std::move(data);
    l.hem = HalfEdgeMesh::create(static_cast<int>(l.data.positions.size()), static_cast<int>(l.data.indices.size() / 3),
                                 reinterpret_cast<const Face *>(l.data.indices.data()));
    return l;
}


/// Unit cube, 12 outward-wound triangles.
Loaded makeCube() {
    MeshData d;
    for (int i = 0; i < 8; ++i) {
        d.positions.push_back({ static_cast<float>(i & 1), static_cast<float>((i >> 1) & 1), static_cast<float>(i >> 2) });
    }
    d.indices = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5 };
    return makeLoaded(std::move(d));
}


/// Torus (genus 1) of 16 x 8 quads, tube radius 0.7 around a circle of radius 2.
Loaded makeTorus() {
    constexpr int around = 16;
    constexpr int tube = 8;
    MeshData      d;
    for (int i = 0; i < around; ++i) {
        for (int j = 0; j < tube; ++j) {
            const double a = 2.0 * std::numbers::pi * i / around;
            const double b = 2.0 * std::numbers::pi * j / tube;
            const double r = 2.0 + 0.7 * std::cos(b);
            d.positions.push_back({ static_cast<float>(r * std::cos(a)), static_cast<float>(r * std::sin(a)),
                                    static_cast<float>(0.7 * std::sin(b)) });
        }
    }
    const auto at = [](int i, int j) {
        return (i % around) * tube + (j % tube);
    };
    for (int i = 0; i < around; ++i) {
        for (int j = 0; j < tube; ++j) {
            d.indices.insert(d.indices.end(), { at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j), at(i + 1, j + 1), at(i, j + 1) });
        }
    }
    return makeLoaded(std::move(d));
}


/// Interior-overlapping triangle pairs inside one packed chart (exact SAT).
long selfOverlaps(const HalfEdgeChart &ac) {
    const auto uvOf = [&](std::size_t i, int k) {
        return ac.uv[static_cast<std::size_t>(localCorner(ac, static_cast<int>(i), k))];
    };
    long overlaps = 0;
    for (std::size_t a = 0; a < ac.faces.size(); ++a) {
        for (std::size_t b = a + 1; b < ac.faces.size(); ++b) {
            overlaps += trianglesOverlapSAT(uvOf(a, 0), uvOf(a, 1), uvOf(a, 2), uvOf(b, 0), uvOf(b, 1), uvOf(b, 2)) ? 1 : 0;
        }
    }
    return overlaps;
}


/// V - E + F of a local mesh.  Connected, so chi = 1 means a disk: genus 0 with one boundary loop.
long eulerCharacteristic(const HalfEdgeMesh &m) {
    long boundary = 0;
    for (HeIndex h = 0; h < m.halfEdgeCount(); ++h) {
        boundary += m.sym(h) == kNoHalfedge ? 1 : 0;
    }
    return static_cast<long>(m.pointCount()) - (static_cast<long>(m.halfEdgeCount()) + boundary) / 2 +
           static_cast<long>(m.faceCount());
}


/// One chart of @p l, flattened (closed charts opened first), packed: must be a disk that packs without overlap.
void requireOpenedDisk(const Loaded &l, ChartSet cs) {
    InjectivityStats st;
    auto             atlas = flattenCharts(l.hem, l.data, cs, {}, &st);

    REQUIRE(atlas.size() == 1);
    CHECK(st.closedOpened == 1);
    CHECK(atlas[0].injective());
    CHECK(eulerCharacteristic(atlas[0].mesh) == 1);
    REQUIRE_NOTHROW(packAtlas(atlas, l.data, 128));
    CHECK(selfOverlaps(atlas[0]) == 0);
}


InjectivityOptions boundOptions() {
    InjectivityOptions o;
    o.maxStretchL2 = 1.05;
    o.maxStretchLinf = 4.0;
    return o;
}

} // namespace


TEST_CASE("Spec 20 a: bound off leaves the plain flatten untouched", "[chartstretch]") {
    const Loaded l = makeLoaded(load(modelOrSkip("cow.obj")));
    ChartSet     cs = segment(l.hem, l.data, 4);

    InjectivityOptions o;
    o.flattener = Flattener::LscmSplit; // no bijective heal: each final chart is its plain LSCM map
    InjectivityStats st;
    const auto       atlas = flattenCharts(l.hem, l.data, cs, o, &st);

    CHECK(st.cuts == 0);
    CHECK(st.closedOpened == 0);
    REQUIRE(atlas.size() == cs.charts.size());
    for (std::size_t c = 0; c < atlas.size(); ++c) {
        const HalfEdgeChart plain = flatten(l.hem, l.data, cs.charts[c]);
        REQUIRE(atlas[c].uv.size() == plain.uv.size());
        for (std::size_t i = 0; i < plain.uv.size(); ++i) {
            REQUIRE((atlas[c].uv[i].x == plain.uv[i].x && atlas[c].uv[i].y == plain.uv[i].y));
        }
    }
}


TEST_CASE("Spec 20 b: a closed cube chart opens to a disk", "[chartstretch]") {
    const Loaded cube = makeCube();
    ChartSet     cs = segment(cube.hem, cube.data, 1);
    REQUIRE(cs.charts[0].seamHalfEdges.empty());
    requireOpenedDisk(cube, cs);
}


TEST_CASE("Spec 20 c: a torus chart opens to a disk", "[chartstretch]") {
    const Loaded torus = makeTorus();
    ChartSet     cs = segment(torus.hem, torus.data, 1);
    cs.charts[0].seamHalfEdges.clear(); // the segmenter flags its own cut edges; the chart must be truly closed
    requireOpenedDisk(torus, cs);
}


TEST_CASE("Spec 20 d: cow as one chart is cut open under the stretch bound", "[chartstretch]") {
    const Loaded     l = makeLoaded(load(modelOrSkip("cow.obj")));
    ChartSet         cs = segment(l.hem, l.data, 1);
    InjectivityStats st;
    auto             atlas = flattenCharts(l.hem, l.data, cs, boundOptions(), &st);

    REQUIRE(atlas.size() == 1);
    CHECK(st.cuts > 0);
    CHECK(st.nonInjectiveAfter == 0);

    // The report is the measurement; the bound itself is met on cow (measured), and a shortfall would print.
    const ChartStretch measured = chartStretch(l.data, atlas[0]);
    CHECK(st.maxChartL2 == Catch::Approx(measured.l2).epsilon(1e-9));
    CHECK(st.maxChartLinf == Catch::Approx(measured.linf).epsilon(1e-9));
    INFO("L2 " << st.maxChartL2 << " Linf " << st.maxChartLinf);
    CHECK(st.maxChartL2 <= 1.05);
    CHECK(st.maxChartLinf <= 4.0);

    REQUIRE_NOTHROW(packAtlas(atlas, l.data, 256));
    CHECK(selfOverlaps(atlas[0]) == 0);
}


TEST_CASE("Spec 20 e: the cut path runs from the worst face to the boundary through interior edges", "[chartstretch]") {
    const Loaded       l = makeLoaded(load(modelOrSkip("cow.obj")));
    ChartSet           cs = segment(l.hem, l.data, 1);
    InjectivityOptions open;
    open.flattener = Flattener::LscmSplit;
    const auto atlas = flattenCharts(l.hem, l.data, cs, open, nullptr); // opens the closed chart, nothing else
    Chart     &chart = cs.charts[0];

    std::vector<double> gamma;
    const ChartStretch  s = chartStretch(l.data, atlas[0], &gamma);
    const int           face = atlas[0].faces[static_cast<std::size_t>(s.worstFace)];
    const auto          path = augmentationPaths(l.hem, l.data, cs, 0, { face }, 1);
    REQUIRE(!path.empty());

    // Vertices of a half-edge; the outline is every vertex of an edge that is not interior to the chart.
    const auto ends = [&](int e) {
        const auto h = static_cast<HeIndex>(e);
        return std::array<int, 2> { static_cast<int>(l.hem.halfedge(l.hem.next(l.hem.next(h))).end),
                                    static_cast<int>(l.hem.halfedge(h).end) };
    };
    const auto interior = [&](HeIndex h) {
        const HeIndex t = l.hem.sym(h);
        return t != kNoHalfedge && cs.faceToChart[l.hem.halfedge(t).face] == 0 &&
               !std::ranges::binary_search(chart.seamHalfEdges, static_cast<int>(h));
    };
    std::vector<char> outline(l.hem.pointCount(), 0);
    for (const int f : chart.faces) {
        for (HeIndex h = 3 * static_cast<HeIndex>(f); h < 3 * static_cast<HeIndex>(f) + 3; ++h) {
            outline[ends(static_cast<int>(h))[0]] |= interior(h) ? 0 : 1;
            outline[ends(static_cast<int>(h))[1]] |= interior(h) ? 0 : 1;
        }
    }

    for (const int e : path) {
        CHECK(interior(static_cast<HeIndex>(e)));
    }
    const auto first = ends(path.front());
    const auto last = ends(path.back());
    const auto corner = [&](int v) {
        return std::ranges::any_of(std::array { 0, 1, 2 }, [&](int k) { return l.hem.halfedge(3 * face + k).end == static_cast<HeIndex>(v); });
    };
    CHECK((corner(first[0]) || corner(first[1])));
    CHECK((outline[last[0]] != 0 || outline[last[1]] != 0));

    addCutEdges(l.hem, chart, path);
    for (const int e : path) {
        CHECK(std::ranges::binary_search(chart.seamHalfEdges, e));
        CHECK(std::ranges::binary_search(chart.seamHalfEdges, static_cast<int>(l.hem.sym(static_cast<HeIndex>(e)))));
    }
}


TEST_CASE("Spec 20 f: two bounded runs are bit-identical", "[chartstretch]") {
    const Loaded l = makeLoaded(load(modelOrSkip("cow.obj")));
    ChartSet     a = segment(l.hem, l.data, 1);
    ChartSet     b = segment(l.hem, l.data, 1);
    const auto   atlasA = flattenCharts(l.hem, l.data, a, boundOptions(), nullptr);
    const auto   atlasB = flattenCharts(l.hem, l.data, b, boundOptions(), nullptr);

    CHECK(a.charts[0].seamHalfEdges == b.charts[0].seamHalfEdges);
    REQUIRE(atlasA[0].uv.size() == atlasB[0].uv.size());
    for (std::size_t i = 0; i < atlasA[0].uv.size(); ++i) {
        REQUIRE((atlasA[0].uv[i].x == atlasB[0].uv[i].x && atlasA[0].uv[i].y == atlasB[0].uv[i].y));
    }
}


TEST_CASE("Spec 20: measurement, start 1, bound 1.05 / 4, pack 1024", "[.chartcut]") {
    for (const char *name : { "cow.obj", "bunny.iv", "horse.iv" }) {
        const Loaded     l = makeLoaded(load(modelOrSkip(name)));
        const auto       t0 = std::chrono::steady_clock::now();
        ChartSet         cs = segment(l.hem, l.data, 1);
        InjectivityStats st;
        auto             atlas = flattenCharts(l.hem, l.data, cs, boundOptions(), &st);
        packAtlas(atlas, l.data, 1024);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

        double density = 0.0; // summed UV area in the [0,1]^2 atlas
        for (const HalfEdgeChart &c : atlas) {
            for (std::size_t i = 0; i < c.faces.size(); ++i) {
                const auto p = [&](int k) {
                    return c.uv[static_cast<std::size_t>(localCorner(c, static_cast<int>(i), k))];
                };
                density += 0.5 * std::abs((p(1).x - p(0).x) * (p(2).y - p(0).y) - (p(2).x - p(0).x) * (p(1).y - p(0).y));
            }
        }
        std::printf("%-9s %zu charts %2d cuts L2 %.3f Linf %.3f density %.3f effective %.3f %6.2fs\n", name, atlas.size(),
                    st.cuts, st.maxChartL2, st.maxChartLinf, density, density / (st.maxChartL2 * st.maxChartL2), seconds);
    }
}
