// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file bijective_detail.h
/// @brief Private: the energy / line-search kernels of `flattenBijective` (Smith & Schaefer 2015),
///        exposed for white-box tests (gradient by finite differences, exactness of the flip bound).
///
/// State vector layout everywhere: `x[2 i], x[2 i + 1]` = (u, v) of local vertex `i`.

#include <memory>
#include <vector>

#include <decimate/chart/parametrize.h>
#include <decimate/types/mesh_data.h>

namespace decimate::chart::detail {

/// One triangle's precomputed 3D side: the inverse of its isometrically flattened reference
/// (paper Figure 3, `R`) and its 3D area `A_P`.
struct RefTri {
    double m00 = 0.0;
    double m01 = 0.0;
    double m10 = 0.0;
    double m11 = 0.0;
    double area = 0.0;
};


/// The chart as the optimizer sees it.
struct Problem {
    int                 n = 0;      ///< Local vertices.
    std::vector<int>    tri;        ///< 3 local corners (a, b, c) per face.
    std::vector<RefTri> ref;        ///< Parallel to the faces.
    std::vector<int>    bEdge;      ///< Boundary edges as (origin, end) local vertex pairs.
    std::vector<int>    bVert;      ///< Boundary vertices, each once.
    double              eps = 0.0;  ///< Barrier support: mean 3D boundary edge length / 4 (section 3.2).
    double              area = 0.0; ///< Total 3D area.
};


/// Build the optimizer's view of `ac`: local vertices, reference triangles, boundary, eps and total area.
[[nodiscard]] Problem buildProblem(const MeshData &data, const HalfEdgeChart &ac);


/// E_D (Eq. 3 summed over the triangles); +inf when a triangle is degenerate or inverted.
/// `grad` (size 2n, accumulated into) is optional; `maxTri` receives the largest per-triangle
/// value when given.
double distortion(const Problem &p, const std::vector<double> &x, std::vector<double> *grad,
                  double *maxTri = nullptr);


/// The barrier's spatial hash (section 3.2, Figure 6), kept across energy evaluations.
///
/// Each boundary vertex is listed ONCE, in the cell of its position `x0` when the hash was built.
/// A later evaluation at `x` pads every query box by `delta = max_v |x_v - x0_v|_inf`, so the
/// candidates are a superset of all pairs within eps (a pair farther than eps contributes exactly 0),
/// and the hash is rebuilt at `x` only once `delta` exceeds one cell.  One cache belongs to one
/// chart's optimisation; it holds no globals and is not thread-safe.
class BarrierCache {
public:
    /// Build the hash at the boundary positions of `x`.
    BarrierCache(const Problem &p, const std::vector<double> &x);
    ~BarrierCache();
    BarrierCache(BarrierCache &&) noexcept;
    BarrierCache &operator=(BarrierCache &&) noexcept;
    BarrierCache(const BarrierCache &) = delete;
    BarrierCache &operator=(const BarrierCache &) = delete;

    /// Times the hash was built (the constructor counts as the first).
    [[nodiscard]] int builds() const;
    /// Cell side of the current hash.
    [[nodiscard]] double cell() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    friend double         barrier(const Problem &, const std::vector<double> &, std::vector<double> *, BarrierCache &);
};


/// E_B (section 3.2): for every boundary edge and every boundary vertex off that edge,
/// `max(0, eps / dist - 1)^2`; +inf when a vertex touches an edge.  `grad` is accumulated into.
/// Pairs are summed per edge in ascending vertex id, so the result does not depend on the hash.
/// This overload builds a fresh hash at `x`.
double barrier(const Problem &p, const std::vector<double> &x, std::vector<double> *grad);


/// Same value and gradient, bit for bit, as the overload above, with the hash taken from `cache`
/// (rebuilt at `x` first when the boundary moved more than one cell since it was built).
double barrier(const Problem &p, const std::vector<double> &x, std::vector<double> *grad,
               BarrierCache &cache);


/// Largest step `t <= cap` such that `x + s d` for every `s` in `[0, t)` inverts no triangle and
/// lets no boundary vertex cross a boundary edge (section 3.3, Eq. 4 for triangles and for
/// (boundary edge, boundary vertex) pairs).  Exactly the first singularity when one lies below
/// `cap`, else `cap`.
double maxStep(const Problem &p, const std::vector<double> &x, const std::vector<double> &d,
               double cap);

} // namespace decimate::chart::detail
