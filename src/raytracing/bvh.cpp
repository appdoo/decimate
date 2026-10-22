// This file is part of the decimate library; see decimate.h for version/license details.
/// @file bvh.cpp
/// @brief Binned-SAH BVH with an exact nearest-|t| line query (see bvh.h for the contract).

#include <decimate/raytracing/bvh.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <decimate/raytracing/raytracing.h>

/// File-local shorthand (an alias, not a using-directive): this file is outside namespace decimate.
namespace math = decimate::math;

namespace rtrace {

namespace {

// Tuning constants of the build and of the conservative traversal bound.
constexpr int   kBins = 16;
constexpr int   kMaxLeaf = 4;
constexpr int   kMaxDepth = 48; // beyond this, split by median (bounds the traversal stack)
constexpr int   kStackSize = 128;
constexpr float kBoxPad = 1e-5f;            // x scene diagonal: float-rounding slack in intersect_triangle
constexpr float kBoundSlack = 1.0f - 1e-4f; // relative slack on the |t| lower bound


/// Axis-aligned box; the default-constructed one is "empty" (inverted), so the first grow() adopts the operand.
struct Box {
    std::array<float, 3> lo { std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                              std::numeric_limits<float>::max() };
    std::array<float, 3> hi { std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(),
                              std::numeric_limits<float>::lowest() };

    /// Union with another box.
    void grow(const Box &o) {
        for (int k = 0; k < 3; k++) {
            lo[k] = std::min(lo[k], o.lo[k]);
            hi[k] = std::max(hi[k], o.hi[k]);
        }
    }

    /// Extend to contain a point.
    void addPoint(const math::Vec3f &p) {
        for (int k = 0; k < 3; k++) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    }

    /// Surface area; an empty box has none (its inverted extent would otherwise give a bogus positive value).
    [[nodiscard]] float area() const {
        const float dx = hi[0] - lo[0];
        const float dy = hi[1] - lo[1];
        const float dz = hi[2] - lo[2];
        return dx < 0.0f ? 0.0f : 2.0f * (dx * dy + dy * dz + dz * dx);
    }
};


/// Faces of one axis binned by centroid: per-bin bounds and counts.
struct Binning {
    std::array<Box, kBins> box;
    std::array<int, kBins> cnt {};
};


/// Best SAH plane found so far: faces in bins < `bin` of `axis` go left.
struct Plane {
    float cost = std::numeric_limits<float>::max();
    int   axis = -1;
    int   bin = 0;
};

} // namespace


/// Prefix/suffix sweeps give left/right (area, count) for each of the kBins-1 planes of `axis`.
static void bestPlaneOnAxis(const Binning &bins, int axis, Plane &best) {
    std::array<float, kBins> leftArea {};
    std::array<int, kBins>   leftCnt {};
    Box                      acc;
    int                      n = 0;

    // Left sweep: bounds and counts of everything left of each plane.
    for (int b = 0; b < kBins - 1; b++) {
        acc.grow(bins.box[b]);
        n += bins.cnt[b];
        leftArea[b] = acc.area();
        leftCnt[b] = n;
    }

    // Right sweep: pair each plane with its right side and keep the cheapest.
    acc = Box {};
    n = 0;
    for (int b = kBins - 1; b > 0; b--) {
        acc.grow(bins.box[b]);
        n += bins.cnt[b];
        const int l = leftCnt[b - 1];
        if (l == 0 || n == 0) {
            continue; // a plane with an empty side splits nothing
        }
        const float cost = leftArea[b - 1] * static_cast<float>(l) + acc.area() * static_cast<float>(n);
        if (cost < best.cost) {
            best = { cost, axis, b };
        }
    }
}


/// Top-down builder; friend of TriangleBvh so it can fill the node array.
struct BvhBuilder {
    TriangleBvh       &bvh;
    std::vector<Box>   box;   // padded per-face box
    std::vector<float> cen;   // 3 per face
    std::vector<int>   order; // face ids, partitioned in place

    /// Bin of face `f` along `axis`; the clamp keeps the centroid at the box's far edge in the last bin.
    [[nodiscard]] int binOf(int f, int axis, const Box &cenBox, float scale) const {
        return std::min(kBins - 1, static_cast<int>((cen[3 * f + axis] - cenBox.lo[axis]) * scale));
    }
    [[nodiscard]] Binning binFaces(int begin, int end, int axis, const Box &cenBox) const;
    int                   build(int begin, int end, int depth);
    bool                  sahSplit(int begin, int end, const Box &cenBox, int &mid);
    int                   medianSplit(int begin, int end, const Box &cenBox);
    void                  computeBoxes(const decimate::MeshData &mesh);
};


/// Bin the faces order[begin,end) by centroid along `axis` (per-bin bounds and counts).
Binning BvhBuilder::binFaces(int begin, int end, int axis, const Box &cenBox) const {
    const float scale = static_cast<float>(kBins) / (cenBox.hi[axis] - cenBox.lo[axis]);
    Binning     bins;
    for (int i = begin; i < end; i++) {
        const int f = order[i];
        const int b = binOf(f, axis, cenBox, scale);
        bins.box[b].grow(box[f]);
        bins.cnt[b]++;
    }
    return bins;
}


// Binned SAH over all three axes. Partitions order[begin,end) and returns the split position in
// `mid`, or false if every candidate leaves one side empty (caller falls back to a median split).
bool BvhBuilder::sahSplit(int begin, int end, const Box &cenBox, int &mid) {
    Plane best;

    // An axis with zero centroid extent cannot separate anything; skip it.
    for (int axis = 0; axis < 3; axis++) {
        if (cenBox.hi[axis] - cenBox.lo[axis] > 0.0f) {
            bestPlaneOnAxis(binFaces(begin, end, axis, cenBox), axis, best);
        }
    }
    if (best.axis < 0) {
        return false;
    }

    const float scale = static_cast<float>(kBins) / (cenBox.hi[best.axis] - cenBox.lo[best.axis]);
    const auto  it = std::partition(order.begin() + begin, order.begin() + end,
                                    [&](int f) { return binOf(f, best.axis, cenBox, scale) < best.bin; });
    mid = static_cast<int>(it - order.begin());
    return mid > begin && mid < end;
}


// Median split on the widest centroid axis (also the NaN/degenerate-centroid fallback).
int BvhBuilder::medianSplit(int begin, int end, const Box &cenBox) {
    int   axis = 0;
    float ext = cenBox.hi[0] - cenBox.lo[0];

    // Widest axis: the best chance to separate the faces.
    for (int k = 1; k < 3; k++) {
        if (cenBox.hi[k] - cenBox.lo[k] > ext) {
            ext = cenBox.hi[k] - cenBox.lo[k];
            axis = k;
        }
    }

    const int mid = begin + (end - begin) / 2;
    std::nth_element(order.begin() + begin, order.begin() + mid, order.begin() + end,
                     [&](int a, int b) { return cen[3 * a + axis] < cen[3 * b + axis]; });
    return mid;
}


/// Recursively build the subtree over order[begin,end) in depth-first order; returns its node index.
/// Layout contract: the left child is always `nodeIdx + 1`, the right child's index is stored in `first`.
int BvhBuilder::build(int begin, int end, int depth) {
    const int nodeIdx = static_cast<int>(bvh.nodes_.size());
    bvh.nodes_.push_back({});

    // Node bounds from the face boxes; centroid bounds drive the split.
    Box nodeBox;
    Box cenBox;
    for (int i = begin; i < end; i++) {
        const int f = order[i];
        nodeBox.grow(box[f]);
        Box c;
        for (int k = 0; k < 3; k++) {
            c.lo[k] = c.hi[k] = cen[3 * f + k];
        }
        cenBox.grow(c);
    }
    for (int k = 0; k < 3; k++) {
        bvh.nodes_[nodeIdx].lo[k] = nodeBox.lo[k];
        bvh.nodes_[nodeIdx].hi[k] = nodeBox.hi[k];
    }

    const int count = end - begin;

    // Small enough: make a leaf.
    if (count <= kMaxLeaf) {
        bvh.nodes_[nodeIdx].first = begin;
        bvh.nodes_[nodeIdx].count = count;
        return nodeIdx;
    }

    // SAH first; the median split is the fallback and also bounds the depth (hence the traversal stack).
    int mid = 0;
    if (depth >= kMaxDepth || !sahSplit(begin, end, cenBox, mid)) {
        mid = medianSplit(begin, end, cenBox);
    }
    build(begin, mid, depth + 1); // left child == nodeIdx + 1
    const int right = build(mid, end, depth + 1);
    bvh.nodes_[nodeIdx].first = right;
    bvh.nodes_[nodeIdx].count = 0;
    return nodeIdx;
}


/// Padded per-face boxes and centroids; the pad (1e-5 x scene diagonal) absorbs float rounding
/// inside intersect_triangle so a box test can never reject a face that kernel would accept.
void BvhBuilder::computeBoxes(const decimate::MeshData &mesh) {
    const int nf = static_cast<int>(mesh.indices.size()) / 3;

    box.resize(nf);
    cen.resize(3 * static_cast<std::size_t>(nf));
    order.resize(nf);
    Box scene;

    // Per-face box and centroid; the scene box only sizes the pad below.
    for (int f = 0; f < nf; f++) {
        order[f] = f;
        Box &bx = box[f];
        for (int c = 0; c < 3; c++) {
            bx.addPoint(mesh.positions[mesh.indices[3 * f + c]]);
        }
        for (int k = 0; k < 3; k++) {
            cen[3 * f + k] = 0.5f * (bx.lo[k] + bx.hi[k]);
        }
        scene.grow(bx);
    }

    // A non-finite scene gets no pad rather than a NaN one.
    const float dx = scene.hi[0] - scene.lo[0];
    const float dy = scene.hi[1] - scene.lo[1];
    const float dz = scene.hi[2] - scene.lo[2];
    const float diag = std::sqrt(dx * dx + dy * dy + dz * dz);
    const float pad = std::isfinite(diag) ? kBoxPad * diag : 0.0f;
    for (auto &bx : box) {
        for (int k = 0; k < 3; k++) {
            bx.lo[k] -= pad;
            bx.hi[k] += pad;
        }
    }
}


TriangleBvh::TriangleBvh(const decimate::MeshData &mesh) {
    const int nf = static_cast<int>(mesh.indices.size()) / 3;

    // Empty mesh: no nodes, every query misses.
    if (nf == 0) {
        return;
    }

    BvhBuilder b { *this, {}, {}, {} };
    b.computeBoxes(mesh);
    nodes_.reserve(static_cast<std::size_t>(nf));
    b.build(0, nf, 0);

    // Slot order = partitioned order; copy vertices in index order v0, v1, v2 (the brute-force order).
    faceId_ = b.order;
    tris_.resize(9 * static_cast<std::size_t>(nf));
    for (int s = 0; s < nf; s++) {
        const int f = faceId_[s];
        for (int c = 0; c < 3; c++) {
            const auto &p = mesh.positions[mesh.indices[3 * f + c]];
            tris_[9 * s + 3 * c + 0] = p.x;
            tris_[9 * s + 3 * c + 1] = p.y;
            tris_[9 * s + 3 * c + 2] = p.z;
        }
    }
}


/// Run the exact kernel over a leaf's triangles and keep the argmin over (|t|, face index).
void TriangleBvh::testLeaf(const Node &n, math::Vec3f &orig,
                           math::Vec3f &dir, float &bestAbs, Hit &best) const {
    for (int s = n.first; s < n.first + n.count; s++) {
        const float *q = &tris_[9 * static_cast<std::size_t>(s)];
        math::Vec3f  v0 { q[0], q[1], q[2] };
        math::Vec3f  v1 { q[3], q[4], q[5] };
        math::Vec3f  v2 { q[6], q[7], q[8] };
        float        t = 0.0f;
        float        u = 0.0f;
        float        v = 0.0f;
        if (!intersect_triangle(orig, dir, v0, v1, v2, &t, &u, &v)) {
            continue;
        }
        const float at = std::fabs(t);
        const int   face = faceId_[s];

        // Argmin over (|t|, face): strict < like the brute loop, plus the lower index on an exact
        // tie (a miss has face == -1, so a tie at the bound maxAbsT itself never wins).
        if (at < bestAbs || (at == bestAbs && best.face >= 0 && face < best.face)) {
            bestAbs = at;
            best = { face, t, u, v };
        }
    }
}


/// Intersect the running [tmin, tmax] with one slab of the line. False if the line misses it.
/// NaN (inf/inf, 0/0) gives no constraint: the bound stays conservative.
static bool clipSlab(float lo, float hi, float o, float d, float &tmin, float &tmax) {
    if (d == 0.0f) {
        return o >= lo && o <= hi;
    }
    float t1 = (lo - o) / d;
    float t2 = (hi - o) / d;
    if (t1 != t1 || t2 != t2) {
        return true;
    }
    if (t1 > t2) {
        std::swap(t1, t2);
    }
    tmin = std::max(tmin, t1);
    tmax = std::min(tmax, t2);
    return true;
}


/// Conservative lower bound on |t| over the part of the LINE inside box [lo, hi].
/// Returns false if the line misses the box.
static bool lineBoxBound(const float *lo, const float *hi, const float *o, const float *d, float &lb) {
    float tmin = -std::numeric_limits<float>::infinity();
    float tmax = std::numeric_limits<float>::infinity();
    for (int k = 0; k < 3; k++) {
        if (!clipSlab(lo[k], hi[k], o[k], d[k], tmin, tmax)) {
            return false;
        }
    }
    if (tmin > tmax) {
        return false;
    }
    const bool straddles = tmin <= 0.0f && tmax >= 0.0f;
    lb = straddles ? 0.0f : std::min(std::fabs(tmin), std::fabs(tmax)) * kBoundSlack;
    return true;
}


/// Traversal stack: node index + the |t| lower bound it was pushed with.
struct TriangleBvh::Stack {
    std::array<int, kStackSize>   node {};
    std::array<float, kStackSize> lb {};
    int                           sp = 0;

    void push(int n, float bound) {
        node[sp] = n;
        lb[sp] = bound;
        sp++;
    }
};


/// Push the children of inner node `cur` that the line can reach, nearer one on top.
void TriangleBvh::pushChildren(int cur, const math::Vec3f &orig,
                               const math::Vec3f &dir, Stack &st) const {
    const Node &n = nodes_[cur];
    const Node &l = nodes_[cur + 1];
    const Node &r = nodes_[n.first];
    float       ll = 0.0f;
    float       rl = 0.0f;
    const bool  hl = lineBoxBound(l.lo, l.hi, orig.data, dir.data, ll);
    const bool  hr = lineBoxBound(r.lo, r.hi, orig.data, dir.data, rl);

    // Nearer child on top so it is visited first and tightens the bound for the other.
    if (hl && hr) {
        const bool leftNear = ll <= rl;
        st.push(leftNear ? n.first : cur + 1, leftNear ? rl : ll); // farther first
        st.push(leftNear ? cur + 1 : n.first, leftNear ? ll : rl);
    } else if (hl) {
        st.push(cur + 1, ll);
    } else if (hr) {
        st.push(n.first, rl);
    }
}


TriangleBvh::Hit TriangleBvh::nearest(math::Vec3f orig, math::Vec3f dir,
                                      float maxAbsT) const {
    Hit best;

    if (nodes_.empty()) {
        return best;
    }

    // Brute force: a NaN or zero direction fails every det test; a NaN origin yields NaN t, which
    // never wins. All three are misses.
    const bool nan = orig.x != orig.x || orig.y != orig.y || orig.z != orig.z || dir.x != dir.x ||
                     dir.y != dir.y || dir.z != dir.z;
    if (nan || (dir.x == 0.0f && dir.y == 0.0f && dir.z == 0.0f)) {
        return best;
    }

    // Seed the traversal with the root; a line that misses the root box misses everything.
    float bestAbs = maxAbsT;
    Stack st;
    float rootLb = 0.0f;
    if (!lineBoxBound(nodes_[0].lo, nodes_[0].hi, orig.data, dir.data, rootLb)) {
        return best;
    }
    st.push(0, rootLb);

    // Depth-first with pruning: skip subtrees whose |t| lower bound cannot beat the current best.
    while (st.sp > 0) {
        st.sp--;
        if (st.lb[st.sp] > bestAbs) {
            continue; // strictly farther than the current best; equal bounds may still tie
        }
        const int   cur = st.node[st.sp];
        const Node &n = nodes_[cur];
        if (n.count > 0) {
            testLeaf(n, orig, dir, bestAbs, best);
            continue;
        }
        pushChildren(cur, orig, dir, st);
    }
    return best;
}

} // namespace rtrace
