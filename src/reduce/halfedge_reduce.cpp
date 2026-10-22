// This file is part of the decimate library; see decimate.h for version/license details.
/// @file halfedge_reduce.cpp
/// @brief Half-edge edge-collapse decimator (see halfedge_reduce.h).
///
/// Modernization of the author's 2003 `CHalfedgeCollapse`. The connectivity +
/// guarded collapse operator live in the shared `corner_mesh` module (the 8-byte
/// Directed-Edges / Corner-Table representation); this file is only the
/// *scheduler*: order edges by initial QEM cost (lowest first), then collapse in
/// active-set passes, recomputing each edge's optimum at visit. Quadrics reuse
/// `QuadricErrorSystem`.

#include <decimate/reduce/halfedge_reduce.h>

#include <decimate/math/vec3.h>
#include <decimate/reduce/corner_mesh.h>
#include <decimate/reduce/quadric_error_system.h>
#include <decimate/types/mesh_types.h>

#include <algorithm>
#include <vector>

namespace decimate {

using corner::He8;
using corner::HeIndex;
using corner::kNoHE;
using math::Vec3f;

namespace {

/// A canonical half-edge with its initial collapse cost, for the one-time sort.
struct EdgeCost {
    double  err = 0.0;
    HeIndex he = 0;
};

} // namespace


/// One canonical half-edge per undirected edge, ordered by initial collapse cost.
static std::vector<HeIndex> buildEdgeOrder(const std::vector<He8> &he, const QuadricErrorSystem &es) {
    const int             numHE = static_cast<int>(he.size());
    std::vector<EdgeCost> edges;
    edges.reserve(static_cast<std::size_t>(numHE));
    for (int i = 0; i < numHE; ++i) {
        const HeIndex h = static_cast<HeIndex>(i);
        const HeIndex s = he[h].sym;
        if (s != kNoHE && h > s) {
            continue; // the twin is the canonical half-edge of this edge
        }

        double err = 0.0;
        Vec3f  p;
        es.optimize(static_cast<int>(corner::originOf(he, h)), static_cast<int>(corner::tipOf(he, h)),
                    s == kNoHE, err, p);
        edges.push_back({ err, h });
    }

    std::ranges::sort(edges, [](const EdgeCost &a, const EdgeCost &b) { return a.err < b.err; });
    std::vector<HeIndex> order;
    order.reserve(edges.size());
    for (const EdgeCost &e : edges) {
        order.push_back(e.he);
    }
    return order;
}


/// One sweep: recompute optimum at visit, collapse if legal; compact `order` to the still-live
/// worklist (drop dead/collapsed, keep guarded-out for retry). Returns whether any edge collapsed.
static bool sweep(std::vector<He8> &he, QuadricErrorSystem &es, std::vector<Vec3f> &pos,
                  std::vector<char> &deadFace, std::vector<char> &deadVert, int &numValid, int target,
                  std::vector<HeIndex> &order) {
    bool        progressed = false;
    std::size_t w = 0;
    for (std::size_t r = 0; r < order.size(); ++r) {
        const HeIndex h = order[r];
        if (deadFace[h / 3u] != 0) {
            continue; // its face is gone
        }
        if (numValid <= target) {
            order[w++] = h; // target reached mid-sweep: keep the rest untouched
            continue;
        }

        const HeIndex p2 = corner::tipOf(he, h);
        const HeIndex p1 = corner::originOf(he, h);
        if (p1 == p2 || deadVert[p1] != 0 || deadVert[p2] != 0) {
            continue; // an endpoint was merged away by an earlier collapse
        }

        // The optimum depends on quadrics that earlier collapses changed, so recompute at visit.
        double err = 0.0;
        Vec3f  p;
        es.optimize(static_cast<int>(p1), static_cast<int>(p2), he[h].sym == kNoHE, err, p);
        if (!corner::collapseLegal(he, pos, h, p1, p2, p)) {
            order[w++] = h; // guarded-out — retry
            continue;
        }

        corner::collapse(he, h, p1, deadFace);
        pos[p1] = p;
        es.update(static_cast<int>(p1), static_cast<int>(p2));
        deadVert[p2] = 1;
        --numValid;
        progressed = true;
    }
    order.resize(w);
    return progressed;
}


/// Repeats sweeps until the vertex target is met or a sweep makes no progress.
static void collapseToTarget(std::vector<He8> &he, QuadricErrorSystem &es, std::vector<Vec3f> &pos,
                             std::vector<char> &deadFace, std::vector<char> &deadVert, int &numValid,
                             int target, std::vector<HeIndex> order) {
    while (numValid > target && !order.empty()) {
        if (!sweep(he, es, pos, deadFace, deadVert, numValid, target, order)) {
            break;
        }
    }
}


MeshData reduceHalfEdge(const MeshData &input, int targetVertexCount) {
    corner::requireValidIndices(input, "reduceHalfEdge");

    const int n = static_cast<int>(input.positions.size());
    if (n == 0 || targetVertexCount >= n) {
        return input;
    }

    const int numFaces = static_cast<int>(input.indices.size()) / 3;
    // Face is three consecutive ints, so the flat index array is read as faces without a copy.
    static_assert(sizeof(Face) == 3 * sizeof(int), "Face must be three packed ints");
    const auto *faces = reinterpret_cast<const Face *>(input.indices.data());

    std::vector<He8>   he = corner::build(input); // 8-byte slab, hash-paired
    std::vector<Vec3f> pos = input.positions;
    QuadricErrorSystem es;
    es.initialize(n, pos.data(), numFaces, faces);


    // Tombstones. Vertices no face references (loader phantoms) start dead so they never count as live.
    std::vector<char> deadFace(static_cast<std::size_t>(numFaces), 0);
    std::vector<char> deadVert(static_cast<std::size_t>(n), 0);
    std::vector<char> referenced(static_cast<std::size_t>(n), 0);
    for (int idx : input.indices) {
        referenced[static_cast<std::size_t>(idx)] = 1;
    }

    int numValid = 0;
    for (int v = 0; v < n; ++v) {
        if (referenced[static_cast<std::size_t>(v)] == 0) {
            deadVert[static_cast<std::size_t>(v)] = 1;
        } else {
            ++numValid;
        }
    }


    std::vector<HeIndex> order = buildEdgeOrder(he, es);
    collapseToTarget(he, es, pos, deadFace, deadVert, numValid, targetVertexCount, std::move(order));

    return corner::reconstruct(he, pos, deadFace, deadVert, numFaces);
}

} // namespace decimate
