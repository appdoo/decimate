// This file is part of the decimate library; see decimate.h for version/license details.
/// @file reducer.cpp
/// @brief `reduce()` — public QEM simplification entry point.
///
/// Heap-driven greedy QEM edge-collapse on the shared 8-byte corner-table
/// (`decimate::corner`): always collapse the current lowest-error edge,
/// guard each collapse (link condition + flip) so the result stays 2-manifold,
/// and place the survivor at the QEM-optimal point. The heap is **lazy**: a
/// stale priority is detected and re-queued at pop time (recompute-and-validate),
/// which gives current-min order without an explicit decrease-key.
///
/// This supersedes the 2003 `CEdgeCollapse` (CSet vector-of-vectors) port: same
/// QEM metric and greedy heap order, on a denser representation, plus the
/// link/flip guard the 2003 reducer lacked.

#include <decimate/reduce/reducer.h>

#include <decimate/math/vec3.h>
#include <decimate/reduce/corner_mesh.h>
#include <decimate/reduce/quadric_error_system.h>
#include <decimate/types/mesh_types.h>

#include <cstddef>
#include <queue>
#include <vector>

namespace decimate {

using math::Vec3f;
using corner::He8;
using corner::HeIndex;
using corner::kNoHE;

namespace {

// A heap candidate: the cached collapse cost + optimal point for half-edge `he`.
// Stale candidates (cost no longer current) are detected and re-queued at pop.
struct Cand {
    double      cost = 0.0;
    HeIndex     he   = 0;
    math::Vec3f p;
};
struct CostGreater {
    bool operator()(const Cand& a, const Cand& b) const noexcept { return a.cost > b.cost; }
};
using CandHeap = std::priority_queue<Cand, std::vector<Cand>, CostGreater>;

Cand evaluate(const std::vector<He8>& he, const QuadricErrorSystem& es, HeIndex h) {
    Cand   c;
    c.he = h;
    es.optimize(static_cast<int>(corner::originOf(he, h)), static_cast<int>(corner::tipOf(he, h)),
                he[h].sym == kNoHE, c.cost, c.p);
    return c;
}

// Mark vertices referenced by no face (e.g. loader phantom verts) dead up front;
// returns the live vertex count.
int markPhantoms(const MeshData& input, std::vector<char>& deadVert) {
    std::vector<char> referenced(deadVert.size(), 0);
    for (int idx : input.indices) {
        referenced[static_cast<std::size_t>(idx)] = 1;
    }
    int live = 0;
    for (std::size_t v = 0; v < deadVert.size(); ++v) {
        if (referenced[v] == 0) {
            deadVert[v] = 1;
        } else {
            ++live;
        }
    }
    return live;
}

}  // namespace

MeshData reduce(const MeshData& input, int targetVertexCount) {
    const int numVerts = static_cast<int>(input.positions.size());
    const int numFaces = static_cast<int>(input.indices.size()) / 3;
    if (numVerts <= targetVertexCount || numFaces == 0) {
        return input;
    }
    const auto* faces = reinterpret_cast<const Face*>(input.indices.data());

    std::vector<He8>   he  = corner::build(input);
    std::vector<Vec3f> pos = input.positions;
    QuadricErrorSystem es;
    es.initialize(numVerts, pos.data(), numFaces, faces);

    std::vector<char> deadFace(static_cast<std::size_t>(numFaces), 0);
    std::vector<char> deadVert(static_cast<std::size_t>(numVerts), 0);
    int               numValid = markPhantoms(input, deadVert);

    CandHeap  heap;
    const int numHE = static_cast<int>(he.size());
    for (int i = 0; i < numHE; ++i) {
        const HeIndex h = static_cast<HeIndex>(i);
        const HeIndex s = he[h].sym;
        if (s == kNoHE || h < s) {  // one candidate per undirected edge
            heap.push(evaluate(he, es, h));
        }
    }

    while (numValid > targetVertexCount && !heap.empty()) {
        const Cand c = heap.top();
        heap.pop();
        const HeIndex h = c.he;
        if (deadFace[h / 3u] != 0) {
            continue;  // edge's face already collapsed
        }
        const HeIndex p2 = corner::tipOf(he, h);
        const HeIndex p1 = corner::originOf(he, h);
        if (p1 == p2 || deadVert[p1] != 0 || deadVert[p2] != 0) {
            continue;  // stale
        }
        const Cand cur = evaluate(he, es, h);
        if (cur.cost != c.cost) {
            heap.push(cur);  // priority went stale -> re-queue with current cost
            continue;
        }
        if (!corner::collapseLegal(he, pos, h, p1, p2, cur.p)) {
            continue;  // link/flip guard rejects; re-queued later via a neighbour collapse
        }
        corner::collapse(he, h, p1, deadFace);
        pos[p1] = cur.p;
        es.update(static_cast<int>(p1), static_cast<int>(p2));
        deadVert[p2] = 1;
        --numValid;
        // No proactive re-queue: incident edges keep their (now stale-low) entries
        // and are demoted on pop via the recompute-and-validate check above —
        // collapses raise neighbour cost, so this is correct greedy order.
    }

    return corner::reconstruct(he, pos, deadFace, deadVert, numFaces);
}

}  // namespace decimate
