// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file bvh.h
/// @brief Exact nearest-hit query over a triangle soup (binned-SAH bounding volume hierarchy).
///
/// `TriangleBvh::nearest` returns EXACTLY what the 2003 brute-force loop returns
/// (legacy2003/split/texturesample.cpp:84,97): the argmin over faces of |t|, where every face is
/// tested with `rtrace::intersect_triangle` (a LINE test: negative t counts) and a face only wins
/// if |t| < |best| (so the lowest face index wins an exact |t| tie). It is an accelerator, never
/// an approximation: node culling is conservative and the leaf test is the same kernel call on the
/// same floats, so face, t, u and v are bit-identical to the brute-force result.

#include <vector>

#include <decimate/math/vec3.h>
#include <decimate/types/mesh_data.h>

namespace rtrace {

/// Immutable acceleration structure over a triangle soup; built once, queried concurrently.
class TriangleBvh {
public:
    /// Nearest-hit record. `face == -1` means miss (u, v, t are then unspecified zeros).
    struct Hit {
        int   face = -1;
        float t = 0.0f;
        float u = 0.0f;
        float v = 0.0f;
    };

    /// Build over all triangles of `mesh` (positions + indices, 3 per face).
    /// Copies the triangle vertices; the mesh need not outlive the BVH.
    /// An empty mesh is valid: every query misses.
    explicit TriangleBvh(const decimate::MeshData &mesh);

    /// Nearest face along the LINE `orig + t * dir` (t of either sign) by smallest |t|.
    ///
    /// @param maxAbsT  Exclusive bound: a hit needs |t| < maxAbsT (the bake passes 100000).
    /// @return the hit whose (|t|, face index) is lexicographically smallest, or a miss.
    ///         NaN/zero directions and NaN origins miss, as the brute-force loop does.
    /// Thread-safe: the query only reads the BVH.
    [[nodiscard]] Hit nearest(decimate::math::Vec3f orig, decimate::math::Vec3f dir,
                              float maxAbsT) const;

    /// Number of triangles indexed.
    [[nodiscard]] int faceCount() const {
        return static_cast<int>(faceId_.size());
    }

private:
    friend struct BvhBuilder;

    /// Flattened tree node in depth-first order: the left child of an inner node is the next node,
    /// the right child's index is stored in `first`.
    struct Node {
        float lo[3]; // NOLINT(modernize-avoid-c-arrays)
        float hi[3]; // NOLINT(modernize-avoid-c-arrays)
        int   first; ///< leaf: first slot in tris_/faceId_; inner: right-child node index
        int   count; ///< leaf: triangle count (> 0); inner: 0 (left child is node + 1)
    };

    struct Stack;

    /// Traversal step: push the reachable children of inner node `cur`, nearer one on top.
    void pushChildren(int cur, const decimate::math::Vec3f &orig, const decimate::math::Vec3f &dir,
                      Stack &st) const;

    /// Leaf step: run the exact triangle kernel over the leaf and update the running argmin
    /// (`bestAbs`, `best`). Takes `orig`/`dir` by non-const reference only because the kernel does.
    void testLeaf(const Node &n, decimate::math::Vec3f &orig, decimate::math::Vec3f &dir,
                  float &bestAbs, Hit &best) const;

    std::vector<Node>  nodes_;
    std::vector<float> tris_;   ///< 9 floats per slot: v0, v1, v2 in face index order
    std::vector<int>   faceId_; ///< slot -> original face index
};

} // namespace rtrace
