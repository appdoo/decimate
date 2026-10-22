// This file is part of the decimate library; see decimate.h for version/license details.
/// @file bijective_mapper.cpp
/// @brief Bijective free-boundary flattening — Smith & Schaefer 2015, "Bijective Parameterization
///        with Free Boundaries", ACM TOG 34(4):70.  Implemented from the paper.
///
/// For one disk chart:
///   1. Tutte embedding (section 3.3, the paper's start point): boundary on the unit circle by 3D
///      arc length, interior by mean-value weights, one sparse LU solve (Eigen, private to this TU).
///   2. L-BFGS on f = E_D + E_B (section 3.3):
///        E_D  isometric distortion, Eq. 1-3:  sum_T A_T (1 + A_T^2 / A_uv^2) ||J_T||_F^2
///        E_B  boundary barrier, section 3.2:  sum over boundary edges (U1,U2) and boundary vertices
///             U_i off the edge of max(0, eps / dist(U1,U2,U_i) - 1)^2, eps = mean 3D boundary
///             edge length / 4, candidates from a spatial hash (Figure 6)
///   3. Flip-avoiding line search (section 3.3): t_max is the smallest positive root of Eq. 4 over
///      the triangles and over the (boundary edge, boundary vertex) pairs the hash yields; the
///      backtracking starts below t_max, so every iterate is bijective.
///
/// Deviations from the paper, both measured on bunny/16 chart 3 (7199 faces, 3964 vertices):
///   - L-BFGS initial matrix.  The paper's plain L-BFGS (H0 = gamma I) crawled: E_D/area 27.7
///     after 1500 iterations, 21.4 after 3000 (history 8..40 and uncapped steps changed nothing),
///     where LSCM reaches 4.6 locally injective.  With H0 = the inverse cotangent Laplacian of the
///     3D chart (`Laplacian` below; the Hessian E_D tends to near the isometric optimum) the same
///     chart converges to 4.01 in 1556 iterations.  The energies, the barrier and the line search
///     are the paper's; only the quasi-Newton scaling differs.
///   - Convergence test: the last 10 iterations together must lower f by less than `tolerance`
///     of it (a single iteration can be starved by the flip bound).  The paper does not state its
///     criterion.
///   - The Tutte start is scaled to the chart's 3D area (E_D has a definite scale).
///   - Line search: starts at min(1, 0.8 t_max) and halves under an Armijo test, the shape of
///     libigl's flip_avoiding_line_search (MPL 2.0; no code taken).  The paper starts at t_max.
///
/// Boundary contract: no GL/OS/I/O/global state; Eigen never reaches a public header.

#include "bijective_detail.h"
#include "chart_quality.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <tuple>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

#include <decimate/chart/parametrize.h>
#include <decimate/math/vec3.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>

namespace decimate::chart {

namespace {

using Vec3 = math::Vec3f;

/// Positive infinity: the energy of an inverted or touching configuration.
constexpr double kInf = std::numeric_limits<double>::infinity();


/// L-BFGS memory: the number of (s, y) pairs kept.
constexpr int kHistory = 8;


/// Three UV points: a triangle (a, b, c), or a boundary edge (a, b) and a boundary vertex c.
struct Tri3 {
    const double *a;
    const double *b;
    const double *c;
};


/// Gradient slots of three points; all null when no gradient is wanted.
struct Grad3 {
    double *a;
    double *b;
    double *c;
};


/// 2x2 Jacobian J = D M of one triangle plus its determinant and squared Frobenius norm.
struct Jac {
    double j00, j01, j10, j11, det, frob;
};


/// Boundary vertices bucketed by uniform cells; a vertex may be listed in several cells (its
/// sweep box during the line search).  Entries are sorted by (column, row) with a direct offset
/// per column, so a query costs one short search per column of its box, not one per cell.
class BoundaryGrid {
public:
    BoundaryGrid(double cell, double x0, double y0)
        : cell_(cell), x0_(x0), y0_(y0) {
    }

    /// List vertex @p v in every cell the box [lo, hi] touches.
    void insert(int v, const double *lo, const double *hi) {
        const auto first = cellOf(lo);
        const auto last = cellOf(hi);
        for (std::int64_t ix = first.first; ix <= last.first; ++ix) {
            for (std::int64_t iy = first.second; iy <= last.second; ++iy) {
                entries_.push_back({ static_cast<std::int32_t>(ix), static_cast<std::int32_t>(iy), v });
            }
        }
    }

    /// Sort and index the columns; call once after the last insert.
    void finish() {
        std::ranges::sort(entries_, [](const Entry &a, const Entry &b) {
            return std::tie(a.ix, a.iy, a.v) < std::tie(b.ix, b.iy, b.v);
        });
        const std::int32_t columns = entries_.empty() ? 0 : entries_.back().ix + 1;
        colStart_.assign(static_cast<std::size_t>(columns) + 1, 0);
        for (const Entry &e : entries_) {
            ++colStart_[static_cast<std::size_t>(e.ix) + 1];
        }
        std::partial_sum(colStart_.begin(), colStart_.end(), colStart_.begin());
    }

    /// Call f(v) for every entry whose cell touches the box [lo, hi] (a vertex spanning several
    /// cells is reported once per cell), column by column, rows ascending.
    template <class F>
    void query(const double *lo, const double *hi, F &&f) const {
        const auto [ix0, iy0] = cellOf(lo);
        const auto [ix1, iy1] = cellOf(hi);
        // No entry lies below cell 0 or past the last column: clamping skips empty space.
        const std::int64_t lastColumn = static_cast<std::int64_t>(colStart_.size()) - 2;
        for (std::int64_t ix = std::max<std::int64_t>(ix0, 0); ix <= std::min(ix1, lastColumn); ++ix) {
            const auto first = entries_.begin() + colStart_[static_cast<std::size_t>(ix)];
            const auto last = entries_.begin() + colStart_[static_cast<std::size_t>(ix) + 1];
            auto       it = std::lower_bound(first, last, iy0, [](const Entry &e, std::int64_t row) { return e.iy < row; });
            for (; it != last && it->iy <= iy1; ++it) {
                f(it->v);
            }
        }
    }

    [[nodiscard]] double cell() const {
        return cell_;
    }

private:
    struct Entry {
        std::int32_t ix;
        std::int32_t iy;
        int          v;
    };

    [[nodiscard]] std::pair<std::int64_t, std::int64_t> cellOf(const double *p) const {
        return { static_cast<std::int64_t>(std::floor((p[0] - x0_) / cell_)),
                 static_cast<std::int64_t>(std::floor((p[1] - y0_) / cell_)) };
    }

    double                    cell_;
    double                    x0_;
    double                    y0_;
    std::vector<Entry>        entries_;
    std::vector<std::int32_t> colStart_;
};


/// Bounding box of the boundary vertices' sweep segments [x, x + cap d] (cap may be 0).
struct Box {
    std::array<double, 2> lo = { kInf, kInf };
    std::array<double, 2> hi = { -kInf, -kInf };
    void                  add(const double *p) {
        lo[0] = std::min(lo[0], p[0]);
        lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]);
        hi[1] = std::max(hi[1], p[1]);
    }
    [[nodiscard]] double extent() const {
        return std::max(hi[0] - lo[0], hi[1] - lo[1]);
    }
};


/// Closest point of edge (a, b) to vertex c: parameter `t` along the edge, `(wx, wy)` = c minus that point.
struct Closest {
    double t = 0.0;
    double wx = 0.0;
    double wy = 0.0;
    double dist = 0.0;
};


/// What edgeBarrier needs besides the chart state: the hash, how far any vertex may have moved since
/// it was built, and a scratch list for the candidates.
struct BarrierQuery {
    const BoundaryGrid &grid;
    double              slack;
    std::vector<int>   &near;
};


using detail::Problem;

/// Initial inverse Hessian of the L-BFGS: the inverse cotangent Laplacian of the 3D chart.  Near
/// the isometric optimum E_D tends to twice the Dirichlet energy, whose Hessian this is, so the
/// quasi-Newton steps are well scaled from the first iteration on (the paper's plain L-BFGS with
/// H0 = gamma I crawls on stiff charts; see the file header for the measurement).  Prefactored once.
class Laplacian {
public:
    Laplacian(const MeshData &data, const HalfEdgeChart &ac, int n);

    /// r = H0 q on the interleaved (u, v) state.
    void apply(const std::vector<double> &q, std::vector<double> &r) const {
        solve<1>({ &q }, { &r });
    }

    /// r = H0 q and rb = H0 qb in one walk over the factor (the walk over the factor is the cost, so two right-hand sides cost little more than one).
    void apply(const std::vector<double> &q, std::vector<double> &r, const std::vector<double> &qb,
               std::vector<double> &rb) const {
        solve<2>({ &q, &qb }, { &r, &rb });
    }

private:
    /// K interleaved (u, v) right-hand sides solved together.  Eigen's solve runs the coordinates one
    /// after the other and each walks the whole factor; this walks it once.  Per coordinate the operations
    /// and their order are those of `ldlt_.solve`, so the result is bit-identical.
    template <std::size_t K>
    void solve(const std::array<const std::vector<double> *, K> &in, const std::array<std::vector<double> *, K> &out) const {
        constexpr std::size_t w = 2 * K; // doubles per vertex in the work vector
        const auto            n = static_cast<std::size_t>(n_);
        const auto           &lower = ldlt_.matrixL().nestedExpression();
        const int            *outer = lower.outerIndexPtr();
        const int            *inner = lower.innerIndexPtr();
        const double         *val = lower.valuePtr();
        const int            *perm = ldlt_.permutationP().indices().data();
        const auto           &diag = ldlt_.vectorD(); // by value in Eigen: bound once

        // x = P q.
        std::vector<double> &x = scratch_;
        x.resize(w * n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = 0; k < K; ++k) {
                x[w * static_cast<std::size_t>(perm[i]) + 2 * k] = (*in[k])[2 * i];
                x[w * static_cast<std::size_t>(perm[i]) + 2 * k + 1] = (*in[k])[2 * i + 1];
            }
        }

        // L x = x (unit lower, column by column).
        for (std::size_t i = 0; i < n; ++i) {
            for (int e = outer[i]; e < outer[i + 1]; ++e) {
                const auto j = w * static_cast<std::size_t>(inner[e]);
                for (std::size_t c = 0; c < w; ++c) {
                    x[j + c] -= val[e] * x[w * i + c];
                }
            }
        }

        // x = D^-1 x (Eigen multiplies by the inverse).
        for (std::size_t i = 0; i < n; ++i) {
            const double inverse = 1.0 / diag[static_cast<Eigen::Index>(i)];
            for (std::size_t c = 0; c < w; ++c) {
                x[w * i + c] *= inverse;
            }
        }

        // L^T x = x, last column first, each product subtracted in turn as Eigen's transposed solve does.
        for (std::size_t i = n; i-- > 0;) {
            std::array<double, w> t {};
            for (std::size_t c = 0; c < w; ++c) {
                t[c] = x[w * i + c];
            }
            for (int e = outer[i]; e < outer[i + 1]; ++e) {
                const auto j = w * static_cast<std::size_t>(inner[e]);
                for (std::size_t c = 0; c < w; ++c) {
                    t[c] -= val[e] * x[j + c];
                }
            }
            for (std::size_t c = 0; c < w; ++c) {
                x[w * i + c] = t[c];
            }
        }

        // r = P^-1 x.
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = 0; k < K; ++k) {
                (*out[k])[2 * i] = x[w * static_cast<std::size_t>(perm[i]) + 2 * k];
                (*out[k])[2 * i + 1] = x[w * static_cast<std::size_t>(perm[i]) + 2 * k + 1];
            }
        }
    }

    Eigen::Index                                       n_;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt_;
    mutable std::vector<double>                        scratch_; ///< Work vector of `solve`, kept across calls.
};


/// An accepted line-search step and the energy there (`t == 0` when nothing decreased f).
struct Step {
    double t = 0.0;
    double f = 0.0;
};


/// The objective at the line-search origin and its slope along the search direction.
struct Descent {
    double f0;
    double slope;
};


/// Euclidean inner product of two equally long vectors.
double dot(const std::vector<double> &a, const std::vector<double> &b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        s += a[i] * b[i];
    }
    return s;
}


/// L-BFGS history: the (s, y) pairs and y.H0 y of each.  The newest pair's y.H0 y is left pending: the
/// next direction solves H0 for it and for its own vector in one walk over the factor.
struct History {
    std::deque<std::vector<double>> s;
    std::deque<std::vector<double>> y;
    std::deque<double>              yH0y;
    bool                            pending = false; ///< The newest yH0y is not computed yet.

    void push(std::vector<double> sk, std::vector<double> yk) {
        yH0y.push_back(0.0);
        pending = true;
        s.push_back(std::move(sk));
        y.push_back(std::move(yk));
        if (static_cast<int>(s.size()) > kHistory) {
            s.pop_front();
            y.pop_front();
            yH0y.pop_front();
        }
    }
    void clear() {
        s.clear();
        y.clear();
        yH0y.clear();
        pending = false;
    }
};

} // namespace


namespace detail {


/// 3D position of local vertex @p li.
static Vec3 pos3(const MeshData &data, const HalfEdgeChart &ac, int li) {
    return data.positions[static_cast<std::size_t>(ac.localToMesh[static_cast<std::size_t>(li)])];
}


/// Isometric reference of one 3D triangle (Figure 3, R): p1 = (0,0), p2 = (|e1|, 0), p3 above the
/// x axis.  Returns the inverse of [p2 - p1, p3 - p1] and the 3D area.
static RefTri referenceTriangle(const Vec3 &p1, const Vec3 &p2, const Vec3 &p3) {
    const Vec3   e1 = p2 - p1;
    const Vec3   e2 = p3 - p1;
    const double x2 = e1.length();
    const double x3 = x2 > 0.0 ? e1.dot(e2) / x2 : 0.0;
    const double y3 = std::sqrt(std::max(0.0, static_cast<double>(e2.lengthSquared()) - x3 * x3));
    const double det = x2 * y3; // 2 * area
    RefTri       r;
    r.area = 0.5 * det;
    if (det <= 0.0) {
        return r; // degenerate 3D triangle: contributes nothing
    }
    // inverse of [[x2, x3], [0, y3]]
    r.m00 = y3 / det;
    r.m01 = -x3 / det;
    r.m10 = 0.0;
    r.m11 = x2 / det;
    return r;
}


/// Jacobian of the map from the reference triangle to the UV triangle (ua, ub, uc).
static Jac jacobian(const RefTri &r, const double *ua, const double *ub, const double *uc) {
    const double d1x = ub[0] - ua[0];
    const double d1y = ub[1] - ua[1];
    const double d2x = uc[0] - ua[0];
    const double d2y = uc[1] - ua[1];
    Jac          j;
    j.j00 = d1x * r.m00 + d2x * r.m10;
    j.j01 = d1x * r.m01 + d2x * r.m11;
    j.j10 = d1y * r.m00 + d2y * r.m10;
    j.j11 = d1y * r.m01 + d2y * r.m11;
    j.det = j.j00 * j.j11 - j.j01 * j.j10;
    j.frob = j.j00 * j.j00 + j.j01 * j.j01 + j.j10 * j.j10 + j.j11 * j.j11;
    return j;
}


/// Eq. 2 times the 3D area for one triangle, and its gradient w.r.t. the three UV corners:
/// dE/dJ = 2 A [(1 + det^-2) J - (||J||^2 / det^2) J^-T], then dE/dD = dE/dJ M^T.
static double triangleEnergy(const RefTri &r, const Jac &j, const Grad3 &g) {
    const double inv2 = 1.0 / (j.det * j.det);
    const double k = 1.0 + inv2;
    const double e = r.area * k * j.frob;
    if (g.a == nullptr) {
        return e;
    }
    const double s = j.frob * inv2 / j.det; // ||J||^2 / det^3, scales adj(J)^T
    const double g00 = 2.0 * r.area * (k * j.j00 - s * j.j11);
    const double g01 = 2.0 * r.area * (k * j.j01 + s * j.j10);
    const double g10 = 2.0 * r.area * (k * j.j10 + s * j.j01);
    const double g11 = 2.0 * r.area * (k * j.j11 - s * j.j00);
    const double bx = g00 * r.m00 + g01 * r.m01; // d/dU2 = (G M^T) column 0
    const double by = g10 * r.m00 + g11 * r.m01;
    const double cx = g00 * r.m10 + g01 * r.m11; // d/dU3 = (G M^T) column 1
    const double cy = g10 * r.m10 + g11 * r.m11;
    g.b[0] += bx;
    g.b[1] += by;
    g.c[0] += cx;
    g.c[1] += cy;
    g.a[0] -= bx + cx;
    g.a[1] -= by + cy;
    return e;
}


/// Grid of the boundary vertices, each listed over its sweep [x_v, x_v + cap d_v].  The cell is
/// 2 eps, coarsened so that no sweep spans more than ~8 cells per axis (early iterations move far).
static BoundaryGrid sweepGrid(const Problem &p, const std::vector<double> &x, const std::vector<double> &d,
                              double cap) {
    Box              all;
    double           widest = 0.0;
    std::vector<Box> sweep(p.bVert.size());
    for (std::size_t k = 0; k < p.bVert.size(); ++k) {
        const std::size_t           i = 2 * static_cast<std::size_t>(p.bVert[k]);
        const std::array<double, 2> end = { x[i] + cap * d[i], x[i + 1] + cap * d[i + 1] };
        sweep[k].add(&x[i]);
        sweep[k].add(end.data());
        all.add(sweep[k].lo.data());
        all.add(sweep[k].hi.data());
        widest = std::max(widest, sweep[k].extent());
    }
    // Cap the columns at ~4 per boundary vertex: a state flung far from a tiny eps must not allocate
    // (extent / eps) column offsets.  A coarser cell only widens the candidate sets.
    const double cell = std::max({ 2.0 * p.eps, widest / 8.0, all.extent() / (4.0 * static_cast<double>(p.bVert.size()) + 64.0) });
    BoundaryGrid grid(cell > 0.0 ? cell : 1.0, all.lo[0], all.lo[1]);
    for (std::size_t k = 0; k < p.bVert.size(); ++k) {
        grid.insert(p.bVert[k], sweep[k].lo.data(), sweep[k].hi.data());
    }
    grid.finish();
    return grid;
}


/// Closest point of the edge `u.a`-`u.b` to the vertex `u.c`.
static Closest closestOnEdge(const Tri3 &u) {
    const double ex = u.b[0] - u.a[0];
    const double ey = u.b[1] - u.a[1];
    const double l2 = ex * ex + ey * ey;
    Closest      c;
    c.t = l2 > 0.0 ? std::clamp(((u.c[0] - u.a[0]) * ex + (u.c[1] - u.a[1]) * ey) / l2, 0.0, 1.0) : 0.0;
    c.wx = u.c[0] - u.a[0] - c.t * ex;
    c.wy = u.c[1] - u.a[1] - c.t * ey;
    c.dist = std::sqrt(c.wx * c.wx + c.wy * c.wy);
    return c;
}


/// Barrier of one (edge a-b, vertex c) pair; gradient accumulated when `g` has slots.
static double pairBarrier(const Tri3 &u, double eps, const Grad3 &g) {
    const auto [t, wx, wy, dist] = closestOnEdge(u);
    if (dist >= eps) {
        return 0.0;
    }
    if (dist <= 0.0) {
        return kInf;
    }
    const double r = eps / dist - 1.0;
    if (g.a == nullptr) {
        return r * r;
    }
    // dE/dd = -2 r eps / d^2; dd/dU_c = w / d, dd/dU_a = -(1 - t) w / d, dd/dU_b = -t w / d
    // (envelope theorem: t is the minimizer of the point-segment distance).
    const double s = -2.0 * r * eps / (dist * dist) / dist;
    const double fx = s * wx;
    const double fy = s * wy;
    g.c[0] += fx;
    g.c[1] += fy;
    g.a[0] -= (1.0 - t) * fx;
    g.a[1] -= (1.0 - t) * fy;
    g.b[0] -= t * fx;
    g.b[1] -= t * fy;
    return r * r;
}


/// Smallest root of a t^2 + b t + c in (0, cap], or cap when there is none; `second` receives
/// the other root (or cap).
static double firstRoot(double a, double b, double c, double cap, double *second) {
    *second = cap;
    if (std::fabs(a) < 1e-300) {
        const double t = b != 0.0 ? -c / b : cap;
        return (t > 0.0 && t < cap) ? t : cap;
    }
    const double disc = b * b - 4.0 * a * c;
    if (disc < 0.0) {
        return cap;
    }
    const double sq = std::sqrt(disc);
    const double q = -0.5 * (b + (b < 0.0 ? -sq : sq)); // stable: no cancellation
    double       r1 = q / a;
    double       r2 = q != 0.0 ? c / q : cap;
    if (r1 > r2) {
        std::swap(r1, r2);
    }
    const auto valid = [cap](double t) {
        return (t > 0.0 && t < cap) ? t : cap;
    };
    *second = valid(r2);
    const double first = valid(r1);
    return first < cap ? first : *second;
}


/// Roots of Eq. 4 for the triangle `u` moving along `d`: the first in (0, cap) is returned, the
/// second through @p second (cap when absent).  The orientation determinant is quadratic in t.
static double orientRoots(const Tri3 &u, const Tri3 &d, double cap, double *second) {
    const double e1x = u.b[0] - u.a[0];
    const double e1y = u.b[1] - u.a[1];
    const double e2x = u.c[0] - u.a[0];
    const double e2y = u.c[1] - u.a[1];
    const double v1x = d.b[0] - d.a[0];
    const double v1y = d.b[1] - d.a[1];
    const double v2x = d.c[0] - d.a[0];
    const double v2y = d.c[1] - d.a[1];
    const double a = v1x * v2y - v1y * v2x;
    const double b = e1x * v2y - e1y * v2x + v1x * e2y - v1y * e2x;
    const double c = e1x * e2y - e1y * e2x;
    return firstRoot(a, b, c, cap, second);
}


/// Is vertex c within the extent of edge (a, b) at the collinear instant @p t along `d`?
static bool withinEdge(const Tri3 &u, const Tri3 &d, double t) {
    const double ex = u.b[0] + t * d.b[0] - u.a[0] - t * d.a[0];
    const double ey = u.b[1] + t * d.b[1] - u.a[1] - t * d.a[1];
    const double wx = u.c[0] + t * d.c[0] - u.a[0] - t * d.a[0];
    const double wy = u.c[1] + t * d.c[1] - u.a[1] - t * d.a[1];
    const double l2 = ex * ex + ey * ey;
    const double s = wx * ex + wy * ey;
    return l2 > 0.0 && s >= -1e-12 * l2 && s <= l2 * (1.0 + 1e-12);
}


/// First t in (0, cap) at which boundary vertex c crosses boundary edge (a, b): the collinearity
/// quadratic of Eq. 4, keeping only roots where the vertex lies within the edge's extent.
static double crossingStep(const Tri3 &u, const Tri3 &d, double cap) {
    double       second = cap;
    const double first = orientRoots(u, d, cap, &second);
    if (first < cap && withinEdge(u, d, first)) {
        return first;
    }
    if (second < cap && withinEdge(u, d, second)) {
        return second;
    }
    return cap;
}


/// Positions (or directions) of boundary edge @p e and boundary vertex @p v in `x`.
static Tri3 edgeVertex(const Problem &p, const std::vector<double> &x, std::size_t e, int v) {
    return { &x[2 * static_cast<std::size_t>(p.bEdge[2 * e])],
             &x[2 * static_cast<std::size_t>(p.bEdge[2 * e + 1])],
             &x[2 * static_cast<std::size_t>(v)] };
}


/// E_B of boundary edge @p e over the boundary vertices the hash yields near it.  The hash holds
/// each vertex at its position when it was built, so the query box grows by `q.slack`.  Candidates
/// within eps are summed in ascending vertex id, so the result is the same for any hash that yields
/// a superset of them (cached or fresh).
static double edgeBarrier(const Problem &p, const std::vector<double> &x, std::vector<double> *grad,
                          const BarrierQuery &q, std::size_t e) {
    const int                   pv = p.bEdge[2 * e];
    const int                   qv = p.bEdge[2 * e + 1];
    const double               *up = &x[2 * static_cast<std::size_t>(pv)];
    const double               *uq = &x[2 * static_cast<std::size_t>(qv)];
    const double                pad = (p.eps + q.slack) * (1.0 + 1e-9); // relative cushion for the rounding of the box
    const std::array<double, 2> lo = { std::min(up[0], uq[0]) - pad, std::min(up[1], uq[1]) - pad };
    const std::array<double, 2> hi = { std::max(up[0], uq[0]) + pad, std::max(up[1], uq[1]) + pad };
    q.near.clear();
    q.grid.query(lo.data(), hi.data(), [&](int v) {
        if (v != pv && v != qv && closestOnEdge(edgeVertex(p, x, e, v)).dist < p.eps) {
            q.near.push_back(v);
        }
    });
    std::ranges::sort(q.near);
    double total = 0.0;
    for (const int v : q.near) {
        Grad3 g { nullptr, nullptr, nullptr };
        if (grad != nullptr) {
            g = { &(*grad)[2 * static_cast<std::size_t>(pv)], &(*grad)[2 * static_cast<std::size_t>(qv)],
                  &(*grad)[2 * static_cast<std::size_t>(v)] };
        }
        total += pairBarrier(edgeVertex(p, x, e, v), p.eps, g);
    }
    return total;
}


/// Lower @p t to the first crossing of boundary edge @p e by any boundary vertex the hash yields
/// in the edge's sweep box.
static double edgeCrossing(const Problem &p, const std::vector<double> &x, const std::vector<double> &d,
                           const BoundaryGrid &grid, std::size_t e, double t) {
    const int                   pv = p.bEdge[2 * e];
    const int                   qv = p.bEdge[2 * e + 1];
    const std::size_t           ip = 2 * static_cast<std::size_t>(pv);
    const std::size_t           iq = 2 * static_cast<std::size_t>(qv);
    Box                         box;
    const std::array<double, 2> pEnd = { x[ip] + t * d[ip], x[ip + 1] + t * d[ip + 1] };
    const std::array<double, 2> qEnd = { x[iq] + t * d[iq], x[iq + 1] + t * d[iq + 1] };
    box.add(&x[ip]);
    box.add(&x[iq]);
    box.add(pEnd.data());
    box.add(qEnd.data());
    grid.query(box.lo.data(), box.hi.data(), [&](int v) {
        if (v == pv || v == qv) {
            return;
        }
        t = crossingStep(edgeVertex(p, x, e, v), edgeVertex(p, d, e, v), t);
    });
    return t;
}


Problem buildProblem(const MeshData &data, const HalfEdgeChart &ac) {
    Problem p;
    p.n = static_cast<int>(ac.uv.size());

    // One reference triangle per face; the 3D area is the weight of its energy.
    const int nf = static_cast<int>(ac.faces.size());
    p.tri.resize(static_cast<std::size_t>(nf) * 3);
    p.ref.resize(static_cast<std::size_t>(nf));
    for (int f = 0; f < nf; ++f) {
        const std::size_t base = 3 * static_cast<std::size_t>(f);
        p.tri[base] = localCorner(ac, f, 0);
        p.tri[base + 1] = localCorner(ac, f, 1);
        p.tri[base + 2] = localCorner(ac, f, 2);
        p.ref[static_cast<std::size_t>(f)] = referenceTriangle(
            pos3(data, ac, p.tri[base]), pos3(data, ac, p.tri[base + 1]), pos3(data, ac, p.tri[base + 2]));
        p.area += p.ref[static_cast<std::size_t>(f)].area;
    }

    // Twin-less half-edges are the chart boundary; collect them as (origin, end) pairs.
    const HalfEdgeMesh &m = ac.mesh;
    std::vector<char>   onBoundary(static_cast<std::size_t>(p.n), 0);
    double              seamLength = 0.0;
    for (HeIndex h = 0; h < m.halfEdgeCount(); ++h) {
        if (m.sym(h) != kNoHalfedge) {
            continue;
        }
        const int o = static_cast<int>(m.halfedge(m.next(m.next(h))).end);
        const int e = static_cast<int>(m.halfedge(h).end);
        p.bEdge.push_back(o);
        p.bEdge.push_back(e);
        onBoundary[static_cast<std::size_t>(o)] = 1;
        onBoundary[static_cast<std::size_t>(e)] = 1;
        seamLength += (pos3(data, ac, e) - pos3(data, ac, o)).length();
    }

    for (int v = 0; v < p.n; ++v) {
        if (onBoundary[static_cast<std::size_t>(v)] != 0) {
            p.bVert.push_back(v);
        }
    }
    // "we choose eps to be the average length of the chart's 3D seam edges divided by 4" (3.2)
    p.eps = p.bEdge.empty() ? 0.0 : seamLength / (0.5 * static_cast<double>(p.bEdge.size())) / 4.0;
    return p;
}


double distortion(const Problem &p, const std::vector<double> &x, std::vector<double> *grad,
                  double *maxTri) {
    double            total = 0.0;
    double            worst = 0.0;
    const std::size_t nf = p.ref.size();
    for (std::size_t f = 0; f < nf; ++f) {
        const RefTri &r = p.ref[f];
        if (r.area <= 0.0) {
            continue; // degenerate in 3D: no reference to measure distortion against
        }
        const std::size_t ia = 2 * static_cast<std::size_t>(p.tri[3 * f]);
        const std::size_t ib = 2 * static_cast<std::size_t>(p.tri[3 * f + 1]);
        const std::size_t ic = 2 * static_cast<std::size_t>(p.tri[3 * f + 2]);
        const Jac         j = jacobian(r, &x[ia], &x[ib], &x[ic]);
        if (!(j.det > 0.0)) {
            return kInf; // flipped or collapsed (also NaN): outside the energy's domain
        }
        Grad3 g { nullptr, nullptr, nullptr };
        if (grad != nullptr) {
            g = { &(*grad)[ia], &(*grad)[ib], &(*grad)[ic] };
        }
        const double e = triangleEnergy(r, j, g);
        total += e;
        worst = std::max(worst, e / r.area);
    }
    if (maxTri != nullptr) {
        *maxTri = worst;
    }
    return total;
}


/// The cached hash, its build-time snapshot and the query scratch behind `BarrierCache`.
struct BarrierCache::Impl {
    BoundaryGrid        grid;
    std::vector<double> x0; ///< Boundary positions the grid was built at, (u, v) per `bVert` entry.
    double              cell;
    int                 builds = 1;
    std::vector<int>    near; ///< Scratch of edgeBarrier.

    Impl(const Problem &p, const std::vector<double> &x)
        : grid(sweepGrid(p, x, x, 0.0)), x0(2 * p.bVert.size()), cell(grid.cell()) {
        snapshot(p, x);
    }
    /// Remember the boundary positions the grid is (re)built at.
    void snapshot(const Problem &p, const std::vector<double> &x) {
        for (std::size_t k = 0; k < p.bVert.size(); ++k) {
            const std::size_t i = 2 * static_cast<std::size_t>(p.bVert[k]);
            x0[2 * k] = x[i];
            x0[2 * k + 1] = x[i + 1];
        }
    }
    /// Largest coordinate move of any boundary vertex since the grid was built.
    [[nodiscard]] double displacement(const Problem &p, const std::vector<double> &x) const {
        double delta = 0.0;
        for (std::size_t k = 0; k < p.bVert.size(); ++k) {
            const std::size_t i = 2 * static_cast<std::size_t>(p.bVert[k]);
            delta = std::max({ delta, std::fabs(x[i] - x0[2 * k]), std::fabs(x[i + 1] - x0[2 * k + 1]) });
        }
        return delta;
    }
    /// Rebuild the grid at `x` once the boundary moved more than a cell.
    void rebuild(const Problem &p, const std::vector<double> &x) {
        grid = sweepGrid(p, x, x, 0.0);
        cell = grid.cell();
        snapshot(p, x);
        ++builds;
    }
};


BarrierCache::BarrierCache(const Problem &p, const std::vector<double> &x)
    : impl_(std::make_unique<Impl>(p, x)) {
}
BarrierCache::~BarrierCache() = default;
BarrierCache::BarrierCache(BarrierCache &&) noexcept = default;
BarrierCache &BarrierCache::operator=(BarrierCache &&) noexcept = default;
int           BarrierCache::builds() const {
    return impl_->builds;
}
double BarrierCache::cell() const {
    return impl_->cell;
}


double barrier(const Problem &p, const std::vector<double> &x, std::vector<double> *grad,
               BarrierCache &cache) {
    if (p.eps <= 0.0) {
        return 0.0; // no boundary, no barrier
    }

    // The hash holds the boundary where it was built; queries pad their box by how far it moved since.
    BarrierCache::Impl &c = *cache.impl_;
    double              slack = c.displacement(p, x);
    if (!(slack <= c.cell)) { // also rebuilds on NaN
        c.rebuild(p, x);
        slack = 0.0;
    }
    const BarrierQuery query { c.grid, slack, c.near };
    double             total = 0.0;
    const std::size_t  ne = p.bEdge.size() / 2;
    for (std::size_t e = 0; e < ne && total < kInf; ++e) {
        total += edgeBarrier(p, x, grad, query, e);
    }
    return total;
}


double barrier(const Problem &p, const std::vector<double> &x, std::vector<double> *grad) {
    BarrierCache fresh(p, x);
    return barrier(p, x, grad, fresh);
}


double maxStep(const Problem &p, const std::vector<double> &x, const std::vector<double> &d,
               double cap) {
    double            t = cap;
    const std::size_t nf = p.ref.size();
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t ia = 2 * static_cast<std::size_t>(p.tri[3 * f]);
        const std::size_t ib = 2 * static_cast<std::size_t>(p.tri[3 * f + 1]);
        const std::size_t ic = 2 * static_cast<std::size_t>(p.tri[3 * f + 2]);
        double            second = t;
        t = orientRoots({ &x[ia], &x[ib], &x[ic] }, { &d[ia], &d[ib], &d[ic] }, t, &second);
    }
    if (p.bVert.empty()) {
        return t; // no boundary to cross
    }

    // Boundary: vertices listed over their sweep up to the triangle bound; each edge queries its
    // own sweep box, which shrinks as t does (section 3.3, last paragraph).
    const BoundaryGrid grid = sweepGrid(p, x, d, t);
    const std::size_t  ne = p.bEdge.size() / 2;
    for (std::size_t e = 0; e < ne; ++e) {
        t = edgeCrossing(p, x, d, grid, e, t);
    }
    return t;
}

} // namespace detail


/// Is the chart's local mesh a topological disk: one boundary loop and Euler characteristic 1?
static bool isDisk(const HalfEdgeChart &ac) {
    const HalfEdgeMesh &m = ac.mesh;
    const std::uint32_t nb =
        static_cast<std::uint32_t>(std::ranges::count_if(m.halfedges(), [](const Halfedge &h) {
            return h.sym == kNoHalfedge;
        }));
    if (nb == 0 || m.boundaryPointIndices().size() != nb) {
        return false; // closed, or more than one boundary loop
    }

    // Euler characteristic V - E + F == 1 for a disk; E follows from F and the boundary count.
    const long long v = m.pointCount();
    const long long f = m.faceCount();
    const long long e = (3 * f + nb) / 2;
    return v - e + f == 1;
}


/// Mean-value weights (Floater 2003) of every directed edge, as sparse triplets of the Tutte
/// system: row i holds sum_j w_ij on the diagonal and -w_ij off it.  Positive weights, so Tutte's
/// theorem applies.
static void meanValueRows(const MeshData &data, const HalfEdgeChart &ac,
                          std::vector<Eigen::Triplet<double>> &rows) {
    const int nf = static_cast<int>(ac.faces.size());
    for (int f = 0; f < nf; ++f) {
        for (int k = 0; k < 3; ++k) {
            const int         i = localCorner(ac, f, k);
            const int         j = localCorner(ac, f, (k + 1) % 3);
            const int         l = localCorner(ac, f, (k + 2) % 3);
            const math::Vec3f a = data.positions[static_cast<std::size_t>(ac.localToMesh[static_cast<std::size_t>(j)])] -
                                  data.positions[static_cast<std::size_t>(ac.localToMesh[static_cast<std::size_t>(i)])];
            const math::Vec3f b = data.positions[static_cast<std::size_t>(ac.localToMesh[static_cast<std::size_t>(l)])] -
                                  data.positions[static_cast<std::size_t>(ac.localToMesh[static_cast<std::size_t>(i)])];
            const double      la = std::max(static_cast<double>(a.length()), 1e-300);
            const double      lb = std::max(static_cast<double>(b.length()), 1e-300);
            const double      theta = std::atan2(static_cast<double>(a.cross(b).length()), static_cast<double>(a.dot(b)));
            const double      tanHalf = std::min(std::tan(0.5 * std::min(theta, std::numbers::pi - 1e-4)), 1e6);
            rows.emplace_back(i, i, tanHalf / la + tanHalf / lb);
            rows.emplace_back(i, j, -tanHalf / la);
            rows.emplace_back(i, l, -tanHalf / lb);
        }
    }
}


/// Boundary loop placed on the unit circle by 3D arc length; `rhs` gets the rows for those vertices.
static void circleBoundary(const MeshData &data, const HalfEdgeChart &ac, std::vector<char> &fixed,
                           Eigen::MatrixX2d &rhs) {
    const std::vector<std::uint32_t> loop = ac.mesh.boundaryPointIndices();
    const std::size_t                nb = loop.size();
    std::vector<double>              arc(nb + 1, 0.0);
    for (std::size_t k = 0; k < nb; ++k) {
        const auto a = static_cast<int>(loop[k]);
        const auto b = static_cast<int>(loop[(k + 1) % nb]);
        const auto pa = data.positions[static_cast<std::size_t>(ac.localToMesh[static_cast<std::size_t>(a)])];
        const auto pb = data.positions[static_cast<std::size_t>(ac.localToMesh[static_cast<std::size_t>(b)])];
        arc[k + 1] = arc[k] + (pb - pa).length();
    }
    const double total = arc[nb] > 0.0 ? arc[nb] : 1.0;
    for (std::size_t k = 0; k < nb; ++k) {
        const double phi = 2.0 * std::numbers::pi * arc[k] / total;
        fixed[loop[k]] = 1;
        rhs(static_cast<Eigen::Index>(loop[k]), 0) = std::cos(phi);
        rhs(static_cast<Eigen::Index>(loop[k]), 1) = std::sin(phi);
    }
}


/// Signed UV area of the chart (twice, summed).
static double signedArea2(const Problem &p, const std::vector<double> &x) {
    double s = 0.0;
    for (std::size_t f = 0; f < p.ref.size(); ++f) {
        const double *a = &x[2 * static_cast<std::size_t>(p.tri[3 * f])];
        const double *b = &x[2 * static_cast<std::size_t>(p.tri[3 * f + 1])];
        const double *c = &x[2 * static_cast<std::size_t>(p.tri[3 * f + 2])];
        s += (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    }
    return s;
}


/// Solve the Tutte system; writes the interleaved state `x`.  Orientation is made positive and the
/// embedding scaled to the 3D area (the isometric energy has a definite scale).
static std::vector<double> tutteState(const MeshData &data, const HalfEdgeChart &ac, const Problem &p) {
    const auto        n = static_cast<Eigen::Index>(p.n);
    std::vector<char> fixed(static_cast<std::size_t>(p.n), 0);
    Eigen::MatrixX2d  rhs = Eigen::MatrixX2d::Zero(n, 2);
    circleBoundary(data, ac, fixed, rhs);

    std::vector<Eigen::Triplet<double>> rows;
    rows.reserve(ac.faces.size() * 9 + static_cast<std::size_t>(p.n));
    meanValueRows(data, ac, rows);
    std::erase_if(rows, [&](const Eigen::Triplet<double> &t) { return fixed[static_cast<std::size_t>(t.row())] != 0; });
    for (Eigen::Index i = 0; i < n; ++i) {
        if (fixed[static_cast<std::size_t>(i)] != 0) {
            rows.emplace_back(static_cast<int>(i), static_cast<int>(i), 1.0);
        }
    }
    Eigen::SparseMatrix<double> sys(n, n);
    sys.setFromTriplets(rows.begin(), rows.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu(sys);
    const Eigen::MatrixX2d                       sol = lu.solve(rhs);

    std::vector<double> x(2 * static_cast<std::size_t>(p.n));
    for (Eigen::Index i = 0; i < n; ++i) {
        x[2 * static_cast<std::size_t>(i)] = sol(i, 0);
        x[2 * static_cast<std::size_t>(i) + 1] = sol(i, 1);
    }

    // E_D needs a positively oriented map at the chart's 3D scale.
    const double area2 = signedArea2(p, x);
    const double flip = area2 < 0.0 ? -1.0 : 1.0;
    const double scale = area2 != 0.0 ? std::sqrt(2.0 * p.area / std::fabs(area2)) : 1.0;
    for (std::size_t i = 0; i < x.size(); i += 2) {
        x[i] *= scale;
        x[i + 1] *= flip * scale;
    }
    return x;
}


/// The chart's current UV as the interleaved state: orientation made positive and scaled to the 3D
/// area, exactly as `tutteState` leaves its embedding (E_D is scale dependent).
static std::vector<double> currentState(const HalfEdgeChart &ac, const Problem &p) {
    std::vector<double> x(2 * ac.uv.size());
    for (std::size_t i = 0; i < ac.uv.size(); ++i) {
        x[2 * i] = ac.uv[i].x;
        x[2 * i + 1] = ac.uv[i].y;
    }
    const double area2 = signedArea2(p, x);
    const double flip = area2 < 0.0 ? -1.0 : 1.0;
    const double scale = std::sqrt(2.0 * p.area / std::fabs(area2));
    for (std::size_t i = 0; i < x.size(); i += 2) {
        x[i] *= scale;
        x[i + 1] *= flip * scale;
    }
    return x;
}


/// Write the interleaved state back to the chart's UVs and re-measure its quality.
static void storeState(const std::vector<double> &x, HalfEdgeChart &ac) {
    for (std::size_t i = 0; i < ac.uv.size(); ++i) {
        ac.uv[i] = { static_cast<float>(x[2 * i]), static_cast<float>(x[2 * i + 1]) };
    }
    measureQuality(ac);
}


Laplacian::Laplacian(const MeshData &data, const HalfEdgeChart &ac, int n)
    : n_(n) {
    std::vector<Eigen::Triplet<double>> rows;
    rows.reserve(ac.faces.size() * 12);
    const int nf = static_cast<int>(ac.faces.size());
    double    trace = 0.0;
    for (int f = 0; f < nf; ++f) {
        for (int k = 0; k < 3; ++k) { // corner i opposite the edge (j, l): weight cot(theta_i) / 2
            const int         i = localCorner(ac, f, k);
            const int         j = localCorner(ac, f, (k + 1) % 3);
            const int         l = localCorner(ac, f, (k + 2) % 3);
            const math::Vec3f a = detail::pos3(data, ac, j) - detail::pos3(data, ac, i);
            const math::Vec3f b = detail::pos3(data, ac, l) - detail::pos3(data, ac, i);
            const double      cross = a.cross(b).length();
            const double      w = cross > 0.0 ? 0.5 * std::max(0.0, static_cast<double>(a.dot(b)) / cross) : 0.0;
            rows.emplace_back(j, j, w);
            rows.emplace_back(l, l, w);
            rows.emplace_back(j, l, -w);
            rows.emplace_back(l, j, -w);
            trace += 2.0 * w;
        }
    }

    // The Laplacian is singular (constants); a tiny diagonal shift makes it invertible.
    const double shift = 1e-6 * trace / static_cast<double>(n);
    for (int i = 0; i < n; ++i) {
        rows.emplace_back(i, i, shift > 0.0 ? shift : 1.0);
    }
    Eigen::SparseMatrix<double> lap(n, n);
    lap.setFromTriplets(rows.begin(), rows.end());
    ldlt_.compute(lap);
}


/// f(x) = E_D + E_B and its gradient.
static double objective(const Problem &p, const std::vector<double> &x, std::vector<double> *grad,
                        detail::BarrierCache &cache) {
    if (grad != nullptr) {
        std::ranges::fill(*grad, 0.0);
    }

    const double ed = detail::distortion(p, x, grad);
    return ed < kInf ? ed + detail::barrier(p, x, grad, cache) : ed;
}


/// Two-loop recursion: d = -H g from the (s, y) history, H0 = gamma L^-1 with the usual
/// gamma = s.y / y.H0 y scaling.
static std::vector<double> lbfgsDirection(History &h, const Laplacian &h0, const std::vector<double> &g) {
    const auto &s = h.s;
    const auto &y = h.y;
    std::vector<double> q = g;
    const std::size_t   m = s.size();
    std::vector<double> alpha(m);
    std::vector<double> rho(m);

    // First loop of the two-loop recursion, newest pair first.
    for (std::size_t k = m; k-- > 0;) {
        rho[k] = 1.0 / dot(y[k], s[k]);
        alpha[k] = rho[k] * dot(s[k], q);
        for (std::size_t i = 0; i < q.size(); ++i) {
            q[i] -= alpha[k] * y[k][i];
        }
    }

    std::vector<double> d(q.size());
    if (h.pending) {
        std::vector<double> h0y(q.size());
        h0.apply(q, d, y[m - 1], h0y);
        h.yH0y[m - 1] = dot(y[m - 1], h0y);
        h.pending = false;
    } else {
        h0.apply(q, d);
    }
    const double gamma = m > 0 ? dot(s[m - 1], y[m - 1]) / h.yH0y[m - 1] : 1.0;
    for (double &v : d) {
        v *= gamma;
    }

    // Second loop, oldest pair first.
    for (std::size_t k = 0; k < m; ++k) {
        const double beta = rho[k] * dot(y[k], d);
        for (std::size_t i = 0; i < d.size(); ++i) {
            d[i] += (alpha[k] - beta) * s[k][i];
        }
    }
    for (double &v : d) {
        v = -v;
    }
    return d;
}


constexpr double kMaxMove = 64.0; ///< Longest line-search move, in extents of the chart's UV box.


/// Backtracking from below the first singularity (Armijo, halving); the trial point is left in
/// `xNew`.
static Step lineSearch(const Problem &p, const std::vector<double> &x, const std::vector<double> &d,
                       const Descent &from, detail::BarrierCache &cache, std::vector<double> &xNew) {
    Step step;
    step.t = std::min(1.0, 0.8 * detail::maxStep(p, x, d, 1.0 / 0.8));
    // A chart that only grows inverts no triangle, so maxStep does not bound the step.  A stiff
    // barrier gradient (a one-triangle sliver) then throws it out to 1e5 times its own size, where
    // the spatial hash walks (extent / eps)^2 empty cells per query.  Bound the move to kMaxMove
    // extents; healthy descents stay far below that.
    const auto [lo, hi] = std::ranges::minmax_element(x);
    double dMax = 0.0;
    for (const double v : d) {
        dMax = std::max(dMax, std::fabs(v));
    }
    step.t = dMax > 0.0 ? std::min(step.t, kMaxMove * (*hi - *lo) / dMax) : step.t;

    // Halve until the Armijo condition holds; 60 halvings reach below any representable gain.
    for (int k = 0; k < 60; ++k) {
        for (std::size_t i = 0; i < x.size(); ++i) {
            xNew[i] = x[i] + step.t * d[i];
        }
        step.f = objective(p, xNew, nullptr, cache);
        if (step.f <= from.f0 + 1e-4 * step.t * from.slope) {
            return step;
        }
        step.t *= 0.5;
    }

    step.t = 0.0; // no decrease found: the caller stops
    return step;
}


constexpr int kWindow = 10; ///< Iterations over which the convergence test measures the decrease.

/// Minimise f from `x` in place; returns the iterations run.  Converged when the last `kWindow`
/// iterations together lowered f by less than `opts.tolerance` of it (a single iteration can be
/// starved by the flip bound without the descent being over).
static int optimise(const Problem &p, const Laplacian &h0, std::vector<double> &x,
                    const BijectiveOptions &opts) {
    std::vector<double>  g(x.size());
    std::vector<double>  gNew(x.size());
    std::vector<double>  xNew(x.size());
    std::deque<double>   past; // f of the last kWindow iterations
    History              h;
    detail::BarrierCache cache(p, x);
    double               f = objective(p, x, &g, cache);
    int                  it = 0;
    for (; it < opts.maxIterations; ++it) {
        std::vector<double> d = lbfgsDirection(h, h0, g);
        double              slope = dot(g, d);
        if (!(slope < 0.0)) { // not a descent direction: restart from the preconditioned gradient
            h.clear();
            d = lbfgsDirection(h, h0, g);
            slope = dot(g, d);
        }
        const Step step = lineSearch(p, x, d, { f, slope }, cache, xNew);
        if (step.t == 0.0) {
            break; // no step lowers f: stuck at a (local) optimum or the flip bound
        }


        // The (s, y) pair needs the gradient at the accepted point.
        const double fNew = step.f;
        objective(p, xNew, &gNew, cache);
        std::vector<double> s(x.size());
        std::vector<double> y(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            s[i] = xNew[i] - x[i];
            y[i] = gNew[i] - g[i];
        }

        if (dot(s, y) > 1e-12 * dot(y, y)) { // curvature condition keeps H positive definite
            h.push(std::move(s), std::move(y));
        }

        x.swap(xNew);
        g.swap(gNew);
        past.push_back(f);
        f = fNew;
        if (static_cast<int>(past.size()) > kWindow) {
            past.pop_front();
        }

        // Converged: a full window of iterations barely lowered f.
        if (static_cast<int>(past.size()) == kWindow && past.front() - f < opts.tolerance * std::fabs(f)) {
            ++it;
            break;
        }
    }
    return it;
}


bool flattenTutte(const MeshData &data, HalfEdgeChart &ac) {
    if (ac.faces.empty() || !isDisk(ac)) {
        return false; // Tutte's theorem needs a disk
    }

    const Problem p = detail::buildProblem(data, ac);
    storeState(tutteState(data, ac, p), ac);
    return true;
}


bool flattenBijective(const MeshData &data, HalfEdgeChart &ac, const BijectiveOptions &opts,
                      BijectiveStats *stats) {
    const auto t0 = std::chrono::steady_clock::now();
    if (ac.faces.empty() || !isDisk(ac)) {
        return false; // the paper's section-4 hole filling is not ported
    }

    // A bijective start is required; the chart's own map only qualifies when it is injective.
    const Problem       p = detail::buildProblem(data, ac);
    const bool          fromCurrent = opts.start == BijectiveStart::Current && ac.injective();
    std::vector<double> x = fromCurrent ? currentState(ac, p) : tutteState(data, ac, p);
    const Laplacian     h0(data, ac, p.n);
    const int           iterations = optimise(p, h0, x, opts);
    storeState(x, ac);
    if (stats != nullptr) {
        stats->iterations = iterations;
        stats->energy = isometricEnergy(data, ac, &stats->maxEnergy);
        stats->seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
    return true;
}


double isometricEnergy(const MeshData &data, const HalfEdgeChart &ac, double *maxTriangle) {
    const Problem       p = detail::buildProblem(data, ac);
    std::vector<double> x(2 * ac.uv.size());
    for (std::size_t i = 0; i < ac.uv.size(); ++i) {
        x[2 * i] = ac.uv[i].x;
        x[2 * i + 1] = ac.uv[i].y;
    }
    double       worst = 0.0;
    const double e = detail::distortion(p, x, nullptr, &worst);
    if (maxTriangle != nullptr) {
        *maxTriangle = e < kInf ? worst : kInf;
    }
    return p.area > 0.0 ? e / p.area : e;
}

} // namespace decimate::chart
