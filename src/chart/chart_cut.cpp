// This file is part of the decimate library; see decimate.h for version/license details.
/// @file chart_cut.cpp
/// @brief Spec 20 phase 2 — cut graphs inside a chart (see chart_cut.h).

#include "chart_cut.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <queue>
#include <utility>

namespace decimate::chart {

namespace {

/// One direction of an interior chart edge.
struct Arc {
    int    to = 0;
    int    edge = 0; ///< Lower half-edge id of the edge.
    double len = 0.0;
};


/// The chart's interior edges as a vertex graph, indexed by mesh vertex.
struct ChartGraph {
    std::vector<std::vector<Arc>> arcs;
    std::vector<char>             onBoundary; ///< Vertex touches an outline edge or a seam.
    std::vector<int>              vertices;   ///< Vertices of the chart, ascending.
    std::vector<int>              edges;      ///< Interior edge ids, ascending.
};


/// Shortest-path tree of a multi-source Dijkstra.
struct PathTree {
    std::vector<double> dist;
    std::vector<int>    pred;     ///< Previous vertex on the path, -1 at a source.
    std::vector<int>    predEdge; ///< Edge to the previous vertex.
};


constexpr double kUnreached = std::numeric_limits<double>::infinity();

} // namespace


/// Origin vertex of half-edge @p h (the tip of its CCW predecessor).
static int originOf(const HalfEdgeMesh &mesh, HeIndex h) {
    return static_cast<int>(mesh.halfedge(mesh.next(mesh.next(h))).end);
}


/// 3D length of the edge between mesh vertices @p a and @p b.
static double edgeLength(const MeshData &data, int a, int b) {
    const auto d = data.positions[static_cast<std::size_t>(a)] - data.positions[static_cast<std::size_t>(b)];
    return std::sqrt(static_cast<double>(d.lengthSquared()));
}


/// True when half-edge @p h of a face of @p chartId is an interior, uncut edge of the chart.
static bool isInterior(const HalfEdgeMesh &mesh, const ChartSet &charts, int chartId, HeIndex h) {
    const HeIndex s = mesh.sym(h);
    if (s == kNoHalfedge || charts.faceToChart[mesh.halfedge(s).face] != chartId) {
        return false;
    }
    const auto &seams = charts.charts[static_cast<std::size_t>(chartId)].seamHalfEdges;
    return !std::ranges::binary_search(seams, static_cast<int>(h));
}


/// The interior-edge graph of @p chartId, with its outline vertices marked.
static ChartGraph buildGraph(const HalfEdgeMesh &mesh, const MeshData &data, const ChartSet &charts, int chartId) {
    const std::size_t n = mesh.pointCount();
    ChartGraph        g;
    g.arcs.resize(n);
    g.onBoundary.assign(n, 0);
    std::vector<char> used(n, 0);

    // Every half-edge of the chart's faces: interior ones become arcs, the rest mark the outline.
    for (const int f : charts.charts[static_cast<std::size_t>(chartId)].faces) {
        for (HeIndex h = static_cast<HeIndex>(f) * 3; h < static_cast<HeIndex>(f) * 3 + 3; ++h) {
            const int from = originOf(mesh, h);
            const int to = static_cast<int>(mesh.halfedge(h).end);
            used[static_cast<std::size_t>(from)] = 1;
            used[static_cast<std::size_t>(to)] = 1;
            if (!isInterior(mesh, charts, chartId, h)) {
                g.onBoundary[static_cast<std::size_t>(from)] = 1;
                g.onBoundary[static_cast<std::size_t>(to)] = 1;
                continue;
            }
            const int edge = static_cast<int>(std::min(h, mesh.sym(h)));
            g.arcs[static_cast<std::size_t>(from)].push_back({ to, edge, edgeLength(data, from, to) });
            g.edges.push_back(edge);
        }
    }
    for (std::size_t v = 0; v < n; ++v) {
        if (used[v] != 0) {
            g.vertices.push_back(static_cast<int>(v));
        }
    }

    // Each interior edge was met from both halves.
    std::ranges::sort(g.edges);
    g.edges.erase(std::unique(g.edges.begin(), g.edges.end()), g.edges.end());
    return g;
}


/// Multi-source Dijkstra over the interior arcs; equal distances go to the lower predecessor id,
/// so the tree does not depend on arc order.
static PathTree dijkstra(const ChartGraph &g, const std::vector<int> &sources) {
    const std::size_t n = g.arcs.size();
    PathTree          t;
    t.dist.assign(n, kUnreached);
    t.pred.assign(n, -1);
    t.predEdge.assign(n, -1);

    using Item = std::pair<double, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<>> queue;
    for (const int s : sources) {
        t.dist[static_cast<std::size_t>(s)] = 0.0;
        queue.emplace(0.0, s);
    }

    while (!queue.empty()) {
        const auto [d, u] = queue.top();
        queue.pop();
        if (d > t.dist[static_cast<std::size_t>(u)]) {
            continue; // stale entry
        }
        for (const Arc &a : g.arcs[static_cast<std::size_t>(u)]) {
            const auto   v = static_cast<std::size_t>(a.to);
            const double nd = d + a.len;
            const bool   better = nd < t.dist[v] || (nd == t.dist[v] && t.pred[v] >= 0 && u < t.pred[v]);
            if (!better) {
                continue;
            }
            t.dist[v] = nd;
            t.pred[v] = u;
            t.predEdge[v] = a.edge;
            queue.emplace(nd, a.to);
        }
    }
    return t;
}


/// The reached vertex with the largest distance (ties: lowest id); -1 if none.
static int farthest(const ChartGraph &g, const PathTree &t) {
    int best = -1;
    for (const int v : g.vertices) {
        const double d = t.dist[static_cast<std::size_t>(v)];
        if (d != kUnreached && (best < 0 || d > t.dist[static_cast<std::size_t>(best)])) {
            best = v;
        }
    }
    return best;
}


/// Edges from @p v back to its source along the tree.
static std::vector<int> pathBack(const PathTree &t, int v) {
    std::vector<int> edges;
    for (; t.pred[static_cast<std::size_t>(v)] >= 0; v = t.pred[static_cast<std::size_t>(v)]) {
        edges.push_back(t.predEdge[static_cast<std::size_t>(v)]);
    }
    return edges;
}


/// Genus 0: the double-sweep path, at least two edges so the opened boundary has three vertices.
static std::vector<int> sphereCut(const ChartGraph &g) {
    const int        a = farthest(g, dijkstra(g, { g.vertices.front() }));
    const PathTree   fromA = dijkstra(g, { a });
    const int        b = farthest(g, fromA);
    std::vector<int> path = pathBack(fromA, b);
    if (path.size() != 1) {
        return path;
    }

    // a - b alone would open a digon; continue from b along its shortest edge to a third vertex.
    const Arc *next = nullptr;
    for (const Arc &arc : g.arcs[static_cast<std::size_t>(b)]) {
        const bool shorter = next == nullptr || arc.len < next->len || (arc.len == next->len && arc.to < next->to);
        if (arc.to != a && shorter) {
            next = &arc;
        }
    }
    if (next != nullptr) {
        path.insert(path.begin(), next->edge);
    }
    return path;
}


/// Visit the unseen interior neighbours of face @p f: mark them seen, mark the edge to them crossed, queue them.
static void crossToNeighbours(const HalfEdgeMesh &mesh, const ChartSet &charts, int chartId, int f,
                              std::vector<char> &seen, std::vector<char> &crossed,
                              std::priority_queue<int, std::vector<int>, std::greater<>> &frontier) {
    for (HeIndex h = static_cast<HeIndex>(f) * 3; h < static_cast<HeIndex>(f) * 3 + 3; ++h) {
        if (!isInterior(mesh, charts, chartId, h)) {
            continue;
        }
        const int other = static_cast<int>(mesh.halfedge(mesh.sym(h)).face);
        if (seen[static_cast<std::size_t>(other)] != 0) {
            continue;
        }
        seen[static_cast<std::size_t>(other)] = 1;
        crossed[std::min(h, mesh.sym(h))] = 1;
        frontier.push(other);
    }
}


/// Dual spanning tree from the seed face, lowest face id first; returns the interior edges it does not cross.
static std::vector<int> cotreeEdges(const HalfEdgeMesh &mesh, const ChartSet &charts, int chartId,
                                    const ChartGraph &g) {
    const Chart                                               &chart = charts.charts[static_cast<std::size_t>(chartId)];
    std::vector<char>                                          seen(mesh.faceCount(), 0);
    std::vector<char>                                          crossed(mesh.halfEdgeCount(), 0);
    std::priority_queue<int, std::vector<int>, std::greater<>> frontier;
    seen[static_cast<std::size_t>(chart.seed)] = 1;
    frontier.push(chart.seed);

    while (!frontier.empty()) {
        const int f = frontier.top();
        frontier.pop();
        crossToNeighbours(mesh, charts, chartId, f, seen, crossed, frontier);
    }

    std::vector<int> cut;
    for (const int e : g.edges) {
        if (crossed[static_cast<std::size_t>(e)] == 0) {
            cut.push_back(e);
        }
    }
    return cut;
}


/// Remove edges with a degree-1 end vertex, repeatedly (the dangling branches of a cut graph).
static std::vector<int> pruned(const HalfEdgeMesh &mesh, const std::vector<int> &edges) {
    const std::size_t             n = mesh.pointCount();
    std::vector<std::vector<int>> incident(n); // edge slots per vertex
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto h = static_cast<HeIndex>(edges[i]);
        incident[static_cast<std::size_t>(originOf(mesh, h))].push_back(static_cast<int>(i));
        incident[mesh.halfedge(h).end].push_back(static_cast<int>(i));
    }

    std::vector<char> alive(edges.size(), 1);
    std::vector<int>  degree(n);
    std::vector<int>  leaves;
    for (std::size_t v = 0; v < n; ++v) {
        degree[v] = static_cast<int>(incident[v].size());
        if (degree[v] == 1) {
            leaves.push_back(static_cast<int>(v));
        }
    }

    for (std::size_t head = 0; head < leaves.size(); ++head) { // NOLINT(modernize-loop-convert): grows while iterated
        const auto v = static_cast<std::size_t>(leaves[head]);
        for (const int slot : incident[v]) {
            if (alive[static_cast<std::size_t>(slot)] == 0) {
                continue;
            }
            alive[static_cast<std::size_t>(slot)] = 0;
            const auto h = static_cast<HeIndex>(edges[static_cast<std::size_t>(slot)]);
            const auto a = static_cast<std::size_t>(originOf(mesh, h));
            const auto other = a == v ? static_cast<std::size_t>(mesh.halfedge(h).end) : a;
            --degree[v];
            if (--degree[other] == 1) {
                leaves.push_back(static_cast<int>(other));
            }
        }
    }

    std::vector<int> kept;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (alive[i] != 0) {
            kept.push_back(edges[i]);
        }
    }
    return kept;
}


bool isClosedChart(const HalfEdgeMesh &mesh, const ChartSet &charts, int chartId) {
    const Chart &chart = charts.charts[static_cast<std::size_t>(chartId)];
    if (chart.faces.empty() || !chart.seamHalfEdges.empty()) {
        return false;
    }
    for (const int f : chart.faces) {
        for (HeIndex h = static_cast<HeIndex>(f) * 3; h < static_cast<HeIndex>(f) * 3 + 3; ++h) {
            if (!isInterior(mesh, charts, chartId, h)) {
                return false;
            }
        }
    }
    return true;
}


std::vector<int> closedChartCut(const HalfEdgeMesh &mesh, const MeshData &data, const ChartSet &charts, int chartId) {
    const ChartGraph g = buildGraph(mesh, data, charts, chartId);
    if (g.vertices.size() < 4) {
        return {};
    }

    // Euler characteristic of the closed chart: V - E + F, with E = 3F / 2.
    const auto faces = static_cast<long>(charts.charts[static_cast<std::size_t>(chartId)].faces.size());
    const long chi = static_cast<long>(g.vertices.size()) - static_cast<long>(g.edges.size()) + faces;

    std::vector<int> cut = chi == 2 ? sphereCut(g) : pruned(mesh, cotreeEdges(mesh, charts, chartId, g));
    std::ranges::sort(cut);
    return cut;
}


std::vector<int> augmentationPaths(const HalfEdgeMesh &mesh, const MeshData &data, const ChartSet &charts, int chartId,
                                   const std::vector<int> &faces, int maxPaths) {
    const ChartGraph g = buildGraph(mesh, data, charts, chartId);
    std::vector<int> sources;
    for (const int v : g.vertices) {
        if (g.onBoundary[static_cast<std::size_t>(v)] != 0) {
            sources.push_back(v);
        }
    }
    const PathTree tree = dijkstra(g, sources);

    std::vector<char> used(mesh.pointCount(), 0);
    std::vector<int>  cut;
    for (const int face : faces) {
        if (maxPaths == 0) {
            break;
        }

        // The corner deepest inside the chart: where a cut relieves the most stretched region.
        int corner = -1;
        for (HeIndex h = static_cast<HeIndex>(face) * 3; h < static_cast<HeIndex>(face) * 3 + 3; ++h) {
            const int    v = static_cast<int>(mesh.halfedge(h).end);
            const double d = tree.dist[static_cast<std::size_t>(v)];
            const bool   deeper = corner < 0 || d > tree.dist[static_cast<std::size_t>(corner)] ||
                                  (d == tree.dist[static_cast<std::size_t>(corner)] && v < corner);
            corner = d != kUnreached && deeper ? v : corner;
        }
        if (corner < 0 || g.onBoundary[static_cast<std::size_t>(corner)] != 0) {
            continue;
        }

        // Vertex-disjoint from the paths already taken, else this face waits for the next round.
        std::vector<int> vertices;
        std::vector<int> edges;
        for (int v = corner; v >= 0; v = tree.pred[static_cast<std::size_t>(v)]) {
            vertices.push_back(v);
            edges.push_back(tree.predEdge[static_cast<std::size_t>(v)]);
        }
        edges.pop_back(); // the source has no edge
        if (std::ranges::any_of(vertices, [&](int v) { return used[static_cast<std::size_t>(v)] != 0; })) {
            continue;
        }
        for (const int v : vertices) {
            used[static_cast<std::size_t>(v)] = 1;
        }
        cut.insert(cut.end(), edges.begin(), edges.end());
        --maxPaths;
    }
    return cut;
}


void addCutEdges(const HalfEdgeMesh &mesh, Chart &chart, const std::vector<int> &edges) {
    for (const int e : edges) {
        chart.seamHalfEdges.push_back(e);
        chart.seamHalfEdges.push_back(static_cast<int>(mesh.sym(static_cast<HeIndex>(e))));
    }
    std::ranges::sort(chart.seamHalfEdges);
    chart.seamHalfEdges.erase(std::unique(chart.seamHalfEdges.begin(), chart.seamHalfEdges.end()),
                              chart.seamHalfEdges.end());
}

} // namespace decimate::chart
