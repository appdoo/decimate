// This file is part of the decimate library; see decimate.h for version/license details.
/// @file straighten_warp.cpp
/// @brief Harmonic interior displacement per chart + foldover detection (section 6).

#include "straighten_warp.h"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace decimate::seam::detail {

namespace {

using Vec2d = std::array<double, 2>;
using SparseMat = Eigen::SparseMatrix<double>;

constexpr double kRidge = 1e-12;      // keeps rows of fully degenerate fans non-singular
constexpr double kDegenerate = 1e-12; // |cross| / (|e1||e2|) below this: no cotangent


/// Displacement field: chain displacements fixed, interiors to be solved.
struct Field {
    std::vector<double> dx;
    std::vector<double> dy;
};


/// Interior unknowns of one chart plus the system over them.
struct ChartSystem {
    std::vector<int>                    unknowns;
    std::vector<int>                   *local = nullptr; // UV vertex -> row, -1 outside
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::MatrixXd                     rhs;
    std::vector<double>                 diag;
};

} // namespace


/// Vertex id as a container index.
static std::size_t idx(int v) {
    return static_cast<std::size_t>(v);
}


/// UV of vertex @p u of the mesh, in double.
static Vec2d uvOf(const TexturedMesh &tm, int u) {
    const auto &p = tm.uvs[idx(u)];
    return { static_cast<double>(p.x), static_cast<double>(p.y) };
}


/// Cotangent of the angle at @p a in triangle (a, b, c); 0 for a degenerate triangle.
static double cotangent(const Vec2d &a, const Vec2d &b, const Vec2d &c) {
    const double e1x = b[0] - a[0];
    const double e1y = b[1] - a[1];
    const double e2x = c[0] - a[0];
    const double e2y = c[1] - a[1];
    const double cross = std::abs((e1x * e2y) - (e1y * e2x));
    const double scale = std::sqrt((e1x * e1x + e1y * e1y) * (e2x * e2x + e2y * e2y));

    // A sliver would give a huge weight; treating it as flat keeps the system well conditioned.
    if (scale <= 0.0 || cross <= kDegenerate * scale) {
        return 0.0;
    }
    return ((e1x * e2x) + (e1y * e2y)) / cross;
}


/// Adds the spring of weight @p w between @p uj and @p uk to the system.
static void addSpring(ChartSystem &sys, const Field &field, int uj, int uk, double w) {
    const int                               rj = (*sys.local)[idx(uj)];
    const int                               rk = (*sys.local)[idx(uk)];
    const std::array<std::array<int, 2>, 2> ends { { { rj, uk }, { rk, uj } } };
    for (std::size_t s = 0; s < 2; ++s) {
        // Only unknown endpoints get a row; the other side is either coupled or Dirichlet data.
        const int row = ends[s][0];
        if (row < 0) {
            continue;
        }
        const int other = ends[s][1];
        const int otherRow = s == 0 ? rk : rj;
        sys.diag[idx(row)] += w;
        if (otherRow >= 0) {
            sys.triplets.emplace_back(row, otherRow, -w);
            continue;
        }
        sys.rhs(row, 0) += w * field.dx[idx(other)];
        sys.rhs(row, 1) += w * field.dy[idx(other)];
    }
}


/// Add the three cotangent springs of triangle @p face.
static void assembleFace(const TexturedMesh &tm, int face, const Field &field, ChartSystem &sys) {
    std::array<int, 3>   u {};
    std::array<Vec2d, 3> p {};
    for (std::size_t k = 0; k < 3; ++k) {
        u[k] = tm.faceU[(3 * idx(face)) + k];
        p[k] = uvOf(tm, u[k]);
    }

    // The spring between two vertices is weighted by the angle opposite their edge.
    for (std::size_t k = 0; k < 3; ++k) {
        const std::size_t j = (k + 1) % 3;
        const std::size_t l = (k + 2) % 3;
        addSpring(sys, field, u[j], u[l], 0.5 * cotangent(p[k], p[j], p[l]));
    }
}


/// Solve the chart system into @p x; false when no factorization gives a finite solution.
static bool solve(const ChartSystem &sys, Eigen::MatrixXd &x) {
    const auto                          n = static_cast<Eigen::Index>(sys.unknowns.size());
    SparseMat                           a(n, n);
    std::vector<Eigen::Triplet<double>> all = sys.triplets;
    for (Eigen::Index i = 0; i < n; ++i) {
        all.emplace_back(i, i, sys.diag[static_cast<std::size_t>(i)] + kRidge);
    }
    a.setFromTriplets(all.begin(), all.end());

    // Negative cotangent weights can make the matrix indefinite: LDLT first, LU as the fallback.
    Eigen::SimplicialLDLT<SparseMat> ldlt(a);
    if (ldlt.info() == Eigen::Success) {
        x = ldlt.solve(sys.rhs);
        if (ldlt.info() == Eigen::Success && x.allFinite()) {
            return true;
        }
    }

    Eigen::SparseLU<SparseMat> lu(a);
    if (lu.info() != Eigen::Success) {
        return false;
    }
    x = lu.solve(sys.rhs);
    return lu.info() == Eigen::Success && x.allFinite();
}


/// Numbers UV vertex @p u as an unknown if no silhouette edge touches it.
static void claimInterior(const Atlas &atlas, int u, std::vector<int> &local,
                          std::vector<int> &unknowns) {
    if (!atlas.incident[idx(u)].empty() || local[idx(u)] >= 0) {
        return;
    }
    local[idx(u)] = static_cast<int>(unknowns.size());
    unknowns.push_back(u);
}


/// The chart's UV vertices that touch no silhouette edge, numbered through @p local.
static std::vector<int> interiorVertices(const Atlas &atlas, const std::vector<int> &faces,
                                         std::vector<int> &local) {
    std::vector<int> unknowns;
    for (const int f : faces) {
        for (std::size_t k = 0; k < 3; ++k) {
            claimInterior(atlas, atlas.tm->faceU[(3 * idx(f)) + k], local, unknowns);
        }
    }
    return unknowns;
}


/// Solves the interior displacement of one chart in place; false when the solve failed.
static bool warpChart(const Atlas &atlas, const std::vector<int> &faces, std::vector<int> &local,
                      Field &field) {
    ChartSystem sys;
    sys.local = &local;
    sys.unknowns = interiorVertices(atlas, faces, local);
    bool ok = true;

    // A chart without interior vertices has nothing to solve: its silhouette fully determines it.
    if (!sys.unknowns.empty()) {
        sys.diag.assign(sys.unknowns.size(), 0.0);
        sys.rhs = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(sys.unknowns.size()), 2);
        for (const int f : faces) {
            assembleFace(*atlas.tm, f, field, sys);
        }
        Eigen::MatrixXd x;
        ok = solve(sys, x);
        for (std::size_t i = 0; ok && i < sys.unknowns.size(); ++i) {
            field.dx[idx(sys.unknowns[i])] = x(static_cast<Eigen::Index>(i), 0);
            field.dy[idx(sys.unknowns[i])] = x(static_cast<Eigen::Index>(i), 1);
        }
    }

    // `local` is shared across charts: reset only the entries this chart set.
    for (const int u : sys.unknowns) {
        local[idx(u)] = -1;
    }
    return ok;
}


/// Twice the signed area of triangle (a, b, c).
static double signedArea(const Vec2d &a, const Vec2d &b, const Vec2d &c) {
    return ((b[0] - a[0]) * (c[1] - a[1])) - ((b[1] - a[1]) * (c[0] - a[0]));
}


/// UV of vertex @p u in the displaced UV array, in double.
static Vec2d uvOf(const std::vector<math::Vec2f> &uvs, int u) {
    return { static_cast<double>(uvs[idx(u)].x), static_cast<double>(uvs[idx(u)].y) };
}


/// True if a triangle of @p faces changed orientation or became degenerate.
static bool folds(const TexturedMesh &tm, const std::vector<math::Vec2f> &after,
                  const std::vector<int> &faces) {
    for (const int f : faces) {
        std::array<int, 3> u {};
        for (std::size_t k = 0; k < 3; ++k) {
            u[k] = tm.faceU[(3 * idx(f)) + k];
        }
        const double before = signedArea(uvOf(tm, u[0]), uvOf(tm, u[1]), uvOf(tm, u[2]));
        const double now = signedArea(uvOf(after, u[0]), uvOf(after, u[1]), uvOf(after, u[2]));

        // Triangles already degenerate in the input cannot be judged and are skipped.
        if (before != 0.0 && (!std::isfinite(now) || now == 0.0 || (before > 0.0) != (now > 0.0))) {
            return true;
        }
    }
    return false;
}


WarpOutput warp(const Atlas &atlas, const ChainTargets &targets) {
    // Displacement of the moved chain vertices is the Dirichlet data; everything else starts at 0.
    const TexturedMesh &tm = *atlas.tm;
    Field               field { std::vector<double>(tm.uvs.size(), 0.0), std::vector<double>(tm.uvs.size(), 0.0) };
    for (std::size_t u = 0; u < tm.uvs.size(); ++u) {
        if (targets.moved[u] != 0) {
            field.dx[u] = targets.pos[u][0] - static_cast<double>(tm.uvs[u].x);
            field.dy[u] = targets.pos[u][1] - static_cast<double>(tm.uvs[u].y);
        }
    }

    // Charts are independent: each interior solve sees only its own silhouette displacements.
    WarpOutput                out;
    std::vector<int>          local(tm.uvs.size(), -1);
    std::vector<std::uint8_t> solved(atlas.chartFaces.size(), 1);
    for (std::size_t c = 0; c < atlas.chartFaces.size(); ++c) {
        solved[c] = warpChart(atlas, atlas.chartFaces[c], local, field) ? 1 : 0;
    }

    // Untouched vertices keep their input bits; only displaced ones are rewritten.
    out.uvs = tm.uvs;
    for (std::size_t u = 0; u < tm.uvs.size(); ++u) {
        if (field.dx[u] != 0.0 || field.dy[u] != 0.0) {
            out.uvs[u] = { static_cast<float>(static_cast<double>(tm.uvs[u].x) + field.dx[u]),
                           static_cast<float>(static_cast<double>(tm.uvs[u].y) + field.dy[u]) };
        }
    }

    for (std::size_t c = 0; c < atlas.chartFaces.size(); ++c) {
        if (solved[c] == 0 || folds(tm, out.uvs, atlas.chartFaces[c])) {
            out.failedCharts.push_back(static_cast<int>(c));
        }
    }
    return out;
}

} // namespace decimate::seam::detail
