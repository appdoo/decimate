// This file is part of the decimate library; see decimate.h for version/license details.
/// @file find_seams.cpp
/// @brief seam::findSeams -- the single UV-based seam predicate (paper section 3).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <decimate/seam/seam_topology.h>
#include <unordered_map>
#include <utility>

namespace decimate::seam {

namespace {

using EdgeMap = std::unordered_map<std::uint64_t, HalfEdgeRef>;

} // namespace


/// Pack a directed vertex pair into one hash key.
static std::uint64_t edgeKey(int from, int to) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(from)) << 32U) |
           static_cast<std::uint32_t>(to);
}


/// Index of corner @p corner (taken mod 3) of face @p face in `faceV` / `faceU`.
static std::size_t cornerSlot(int face, int corner) {
    return 3 * static_cast<std::size_t>(face) + static_cast<std::size_t>(corner % 3);
}


/// Surface vertex id at corner @p corner of face @p face.
static int vertexAt(const TexturedMesh &tm, int face, int corner) {
    return tm.faceV[cornerSlot(face, corner)];
}


/// UV vertex id at corner @p corner of face @p face.
static int uvIdAt(const TexturedMesh &tm, int face, int corner) {
    return tm.faceU[cornerSlot(face, corner)];
}


/// UV position at corner @p corner of face @p face.
static math::Vec2f uvAt(const TexturedMesh &tm, int face, int corner) {
    return tm.uvs[static_cast<std::size_t>(uvIdAt(tm, face, corner))];
}


/// Directed surface edge -> half-edge; on a non-manifold edge the last face wins.
static EdgeMap directedEdges(const TexturedMesh &tm) {
    EdgeMap edges;
    edges.reserve(tm.faceV.size());
    for (int f = 0; f < static_cast<int>(tm.faceCount()); ++f) {
        for (int k = 0; k < 3; ++k) {
            edges[edgeKey(vertexAt(tm, f, k), vertexAt(tm, f, k + 1))] = { f, k };
        }
    }
    return edges;
}


/// Undirected edges as sorted (smaller, larger) vertex pairs.
static std::vector<std::pair<int, int>> undirectedEdges(const EdgeMap &edges) {
    std::vector<std::pair<int, int>> out;
    out.reserve(edges.size());
    for (const auto &entry : edges) {
        const auto from = static_cast<int>(entry.first >> 32U);
        const auto to = static_cast<int>(entry.first & 0xFFFFFFFFU);
        out.emplace_back(std::min(from, to), std::max(from, to));
    }

    // Sorted, deduplicated output makes the pair order (and so the seam order) deterministic.
    std::ranges::sort(out);
    out.erase(std::ranges::unique(out).begin(), out.end());
    return out;
}


/// True when the two half-edges reference the same UV edge (by id or by exact value).
static bool uvEdgesMatch(const TexturedMesh &tm, HalfEdgeRef fwd, HalfEdgeRef bwd) {
    const int fs = uvIdAt(tm, fwd.face, fwd.corner);
    const int fe = uvIdAt(tm, fwd.face, fwd.corner + 1);
    const int bs = uvIdAt(tm, bwd.face, bwd.corner);
    const int be = uvIdAt(tm, bwd.face, bwd.corner + 1);
    if (fs == be && fe == bs) {
        return true;
    }

    // Distinct UV ids can still sit at identical positions; the value test makes such edges match.
    const math::Vec2f f0 = tm.uvs[static_cast<std::size_t>(fs)];
    const math::Vec2f f1 = tm.uvs[static_cast<std::size_t>(fe)];
    const math::Vec2f b0 = tm.uvs[static_cast<std::size_t>(bs)];
    const math::Vec2f b1 = tm.uvs[static_cast<std::size_t>(be)];
    return f0.x == b1.x && f0.y == b1.y && f1.x == b0.x && f1.y == b0.y;
}


/// Orientation of @p c relative to the line a->b (twice the signed triangle area).
static double orientation(math::Vec2f a, math::Vec2f b, math::Vec2f c) {
    const double ax = static_cast<double>(a.x) - static_cast<double>(c.x);
    const double ay = static_cast<double>(a.y) - static_cast<double>(c.y);
    const double bx = static_cast<double>(b.x) - static_cast<double>(c.x);
    const double by = static_cast<double>(b.y) - static_cast<double>(c.y);
    return ax * by - ay * bx;
}


/// -1, 0 or +1.
static int signOf(double v) {
    return (v > 0.0) - (v < 0.0);
}


/// A matching edge is a foldover when the opposite corners lie on the same side of it.
static bool isFoldover(const TexturedMesh &tm, HalfEdgeRef fwd, HalfEdgeRef bwd) {
    const math::Vec2f a = uvAt(tm, fwd.face, fwd.corner);
    const math::Vec2f b = uvAt(tm, fwd.face, fwd.corner + 1);
    const math::Vec2f c = uvAt(tm, fwd.face, fwd.corner + 2);
    const math::Vec2f d = uvAt(tm, bwd.face, bwd.corner + 2);
    return signOf(orientation(a, b, c)) == signOf(orientation(a, b, d));
}


/// Euclidean length of the surface edge @p he, in object space.
static double edgeLength(const TexturedMesh &tm, HalfEdgeRef he) {
    const math::Vec3f p = tm.positions[static_cast<std::size_t>(vertexAt(tm, he.face, he.corner))];
    const math::Vec3f q =
        tm.positions[static_cast<std::size_t>(vertexAt(tm, he.face, he.corner + 1))];
    const double dx = static_cast<double>(q.x) - static_cast<double>(p.x);
    const double dy = static_cast<double>(q.y) - static_cast<double>(p.y);
    const double dz = static_cast<double>(q.z) - static_cast<double>(p.z);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}


/// Classify the manifold edge (@p fwd, @p bwd) into @p topo.
static void classifyPair(const TexturedMesh &tm, HalfEdgeRef fwd, HalfEdgeRef bwd,
                         SeamTopology &topo) {
    if (!uvEdgesMatch(tm, fwd, bwd)) {
        topo.seams.push_back({ fwd, bwd, edgeLength(tm, fwd) });
        return;
    }
    if (isFoldover(tm, fwd, bwd)) {
        topo.foldovers.push_back(fwd);
    }
}


SeamTopology findSeams(const TexturedMesh &tm) {
    const EdgeMap edges = directedEdges(tm);
    SeamTopology  topo;
    for (const auto &[lo, hi] : undirectedEdges(edges)) {
        if (lo == hi) {
            continue; // degenerate face edge: no direction, no seam
        }

        // Both directions present = interior edge; exactly one = open boundary.
        const auto fwd = edges.find(edgeKey(lo, hi));
        const auto bwd = edges.find(edgeKey(hi, lo));
        if (fwd != edges.end() && bwd != edges.end()) {
            classifyPair(tm, fwd->second, bwd->second, topo);
            continue;
        }
        topo.boundary.push_back(fwd != edges.end() ? fwd->second : bwd->second);
    }
    return topo;
}

} // namespace decimate::seam
