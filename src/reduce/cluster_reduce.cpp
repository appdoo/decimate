// This file is part of the decimate library; see decimate.h for version/license details.
/// @file cluster_reduce.cpp
/// @brief Octree vertex-clustering decimator (see cluster_reduce.h).
///
/// Pipeline: vertex AABB → per-vertex quadrics (QError) → adaptive octree
/// (subdivide while a node holds > maxLeaf vertices) → per-cell quadric-optimal
/// representative (clamped to the cell box) → face remap + degenerate drop.
/// The octree is internal here (a flat counting-sort partition over a vertex
/// index array, not a pointer tree) — extract it only if a second caller appears.

#include <decimate/reduce/cluster_reduce.h>

#include <decimate/math/vec3.h>
#include <decimate/reduce/qerror.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace decimate {

using math::Vec3f;

namespace {

/// Recursion bound: 8^20 cells is far finer than any real mesh needs; it exists only so a blob of
/// (near-)coincident points terminates instead of recursing forever — they merge into one cell at
/// the bound, which is the correct result.
constexpr int kMaxDepth = 20;


/// Axis-aligned box, lean (just two corners) — unlike math::AABB, which carries an 8-corner vector
/// unsuited to one-per-octree-node use.
struct Box {
    Vec3f               lo;
    Vec3f               hi;
    [[nodiscard]] Vec3f mid() const noexcept {
        return (lo + hi) * 0.5f;
    }
};


/// Adaptive octree clustering: partitions a vertex-index array in place and stamps each vertex with
/// its leaf cell id. No pointer tree — each recursion owns a [begin,end) slice of `order`,
/// counting-sorted into 8 octant groups. `partitionOctants` and `subdivide` are defined below, after the
/// static helpers they call.
struct OctreeClusterer {
    const std::vector<Vec3f> &pos;
    int                       maxLeaf;
    std::vector<int>          order;     // vertex indices, partitioned in place
    std::vector<int>          scratch;   // counting-sort destination buffer
    std::vector<int>          clusterId; // result: vertex -> leaf cell id
    int                       clusterCount = 0;

    OctreeClusterer(const std::vector<Vec3f> &p, int maxLeafVerts)
        : pos(p),
          maxLeaf(maxLeafVerts),
          order(p.size()),
          scratch(p.size()),
          clusterId(p.size(), -1) {
        for (int i = 0; i < static_cast<int>(order.size()); ++i) {
            order[i] = i;
        }
    }

    void run(const Box &root) {
        subdivide(0, static_cast<int>(order.size()), root, 0);
    }

    /// Gives every vertex of `order[begin,end)` one fresh cell id.
    void assignLeaf(int begin, int end) {
        const int id = clusterCount++;
        for (int i = begin; i < end; ++i) {
            clusterId[order[i]] = id;
        }
    }

    /// Counting-sort order[begin,end) by octant; return the 9 group boundaries.
    std::array<int, 9> partitionOctants(int begin, int end, const Vec3f &mid);

    /// Splits the slice into octants until it is small enough (or deep enough) to be one cell.
    void subdivide(int begin, int end, const Box &box, int depth);
};

} // namespace


/// Component-wise minimum.
static Vec3f minv(const Vec3f &a, const Vec3f &b) noexcept {
    return { std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z) };
}


/// Component-wise maximum.
static Vec3f maxv(const Vec3f &a, const Vec3f &b) noexcept {
    return { std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z) };
}


/// Keeps a representative inside its cell's vertex box.
static Vec3f clampToBox(const Vec3f &p, const Vec3f &lo, const Vec3f &hi) noexcept {
    return { std::clamp(p.x, lo.x, hi.x), std::clamp(p.y, lo.y, hi.y), std::clamp(p.z, lo.z, hi.z) };
}


/// Octant 0..7 from the per-axis sign relative to the split point.
static int octantOf(const Vec3f &p, const Vec3f &mid) noexcept {
    return (p.x >= mid.x ? 1 : 0) | (p.y >= mid.y ? 2 : 0) | (p.z >= mid.z ? 4 : 0);
}


/// The box of one octant of @p b (bit 0/1/2 of @p octant = upper half in x/y/z).
static Box childBox(const Box &b, int octant) noexcept {
    const Vec3f m = b.mid();
    Box         c;
    c.lo.x = (octant & 1) ? m.x : b.lo.x;
    c.hi.x = (octant & 1) ? b.hi.x : m.x;
    c.lo.y = (octant & 2) ? m.y : b.lo.y;
    c.hi.y = (octant & 2) ? b.hi.y : m.y;
    c.lo.z = (octant & 4) ? m.z : b.lo.z;
    c.hi.z = (octant & 4) ? b.hi.z : m.z;
    return c;
}


/// Bounding box of a non-empty point set — the octree root.
static Box boundsOf(const std::vector<Vec3f> &pts) noexcept {
    Box b { pts[0], pts[0] };
    for (const Vec3f &p : pts) {
        b.lo = minv(b.lo, p);
        b.hi = maxv(b.hi, p);
    }
    return b;
}


std::array<int, 9> OctreeClusterer::partitionOctants(int begin, int end, const Vec3f &mid) {
    // Histogram, so every octant knows where its group starts.
    std::array<int, 8> count {};
    for (int i = begin; i < end; ++i) {
        ++count[octantOf(pos[order[i]], mid)];
    }
    std::array<int, 9> bound {};
    bound[0] = begin;
    for (int o = 0; o < 8; ++o) {
        bound[o + 1] = bound[o] + count[o];
    }

    // Scatter into scratch, then copy back: stable, and `order` keeps its slice in place.
    std::array<int, 8> cursor {};
    for (int o = 0; o < 8; ++o) {
        cursor[o] = bound[o];
    }
    for (int i = begin; i < end; ++i) {
        const int v = order[i];
        scratch[cursor[octantOf(pos[v], mid)]++] = v;
    }
    for (int i = begin; i < end; ++i) {
        order[i] = scratch[i];
    }
    return bound;
}


void OctreeClusterer::subdivide(int begin, int end, const Box &box, int depth) {
    if (end - begin <= maxLeaf || depth >= kMaxDepth) {
        assignLeaf(begin, end);
        return;
    }

    const std::array<int, 9> bound = partitionOctants(begin, end, box.mid());
    for (int o = 0; o < 8; ++o) {
        if (bound[o] == bound[o + 1]) {
            continue; // empty octant
        }
        subdivide(bound[o], bound[o + 1], childBox(box, o), depth + 1);
    }
}


/// Accumulates one normalized face-plane quadric onto each of its three vertices, skipping degenerate
/// (zero-area) faces whose normal would be undefined.
static std::vector<QError> vertexQuadrics(const MeshData &mesh) {
    std::vector<QError> vq(mesh.positions.size());
    const int           triCount = static_cast<int>(mesh.indices.size()) / 3;
    for (int t = 0; t < triCount; ++t) {
        const std::size_t base = static_cast<std::size_t>(t) * 3;
        const int         a = mesh.indices[base];
        const int         b = mesh.indices[base + 1];
        const int         c = mesh.indices[base + 2];
        const Vec3f       n = (mesh.positions[b] - mesh.positions[a])
                            .cross(mesh.positions[c] - mesh.positions[a]);
        if (n.lengthSquared() < 1e-20f) {
            continue; // degenerate triangle — no plane
        }

        QError q;
        q.setPlane(mesh.positions[a], mesh.positions[b], mesh.positions[c]);
        vq[a].add(q);
        vq[b].add(q);
        vq[c].add(q);
    }
    return vq;
}


/// One representative per cell: quadric-optimal point (centroid fallback on a singular solve),
/// clamped to the cell's vertex AABB.
static std::vector<Vec3f> cellRepresentatives(const MeshData         &mesh,
                                              const std::vector<int> &clusterId,
                                              int                     cellCount) {
    const std::vector<QError> vq = vertexQuadrics(mesh);

    std::vector<QError> q(cellCount);
    std::vector<Vec3f>  lo(cellCount);
    std::vector<Vec3f>  hi(cellCount);
    std::vector<Vec3f>  sum(cellCount);
    std::vector<int>    cnt(cellCount, 0);

    // Fold every vertex into its cell: summed quadric, vertex box, centroid accumulator.
    const int n = static_cast<int>(mesh.positions.size());
    for (int v = 0; v < n; ++v) {
        const int    c = clusterId[v];
        const Vec3f &p = mesh.positions[v];
        q[c].add(vq[v]);
        lo[c] = (cnt[c] == 0) ? p : minv(lo[c], p);
        hi[c] = (cnt[c] == 0) ? p : maxv(hi[c], p);
        sum[c] += p;
        ++cnt[c];
    }

    // The box clamp stops a near-singular solve from flinging a representative into a spike.
    std::vector<Vec3f> reps(cellCount);
    for (int c = 0; c < cellCount; ++c) {
        Vec3f p;
        if (!q[c].optimize(p)) {
            p = sum[c] / static_cast<float>(cnt[c]); // centroid fallback
        }
        reps[c] = clampToBox(p, lo[c], hi[c]);
    }
    return reps;
}


/// Remaps every face through clusterId, dropping any triangle whose three vertices no longer map to
/// three distinct cells.
static std::vector<int> remapFaces(const MeshData &mesh, const std::vector<int> &clusterId) {
    std::vector<int> out;
    out.reserve(mesh.indices.size());
    const int triCount = static_cast<int>(mesh.indices.size()) / 3;
    for (int t = 0; t < triCount; ++t) {
        const std::size_t base = static_cast<std::size_t>(t) * 3;
        const int         a = clusterId[mesh.indices[base]];
        const int         b = clusterId[mesh.indices[base + 1]];
        const int         c = clusterId[mesh.indices[base + 2]];
        if (a == b || b == c || a == c) {
            continue; // degenerate after collapse
        }

        out.push_back(a);
        out.push_back(b);
        out.push_back(c);
    }
    return out;
}


MeshData reduceClustering(const MeshData &input, int targetVertexCount) {
    const int n = static_cast<int>(input.positions.size());
    if (n == 0 || targetVertexCount >= n) {
        return input;
    }

    // Leaf capacity n/target makes the cell count land near the target.
    const int target = std::max(1, targetVertexCount);
    const int maxLeaf = std::max(1, static_cast<int>(std::llround(
                                        static_cast<double>(n) / target)));

    OctreeClusterer clusterer(input.positions, maxLeaf);
    clusterer.run(boundsOf(input.positions));

    MeshData out;
    out.positions = cellRepresentatives(input, clusterer.clusterId, clusterer.clusterCount);
    out.indices = remapFaces(input, clusterer.clusterId);
    return out;
}

} // namespace decimate
