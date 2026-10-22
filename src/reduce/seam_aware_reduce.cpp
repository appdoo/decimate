// This file is part of the decimate library; see decimate.h for version/license details.
/// @file seam_aware_reduce.cpp
/// @brief reduceSeamAware / classifySeamEdges -- 5D QEM edge-collapse that keeps the UV
///        parametrization seamless (Liu et al. 2017, section 5).  See seam_aware_reduce.h.
///
/// Layout: (1) topology predicates on the two corner tables, (2) the 5D quadrics and the
/// closed-form placement, (3) edge evaluation, legality and the collapse itself, (4) the
/// lazy-heap driver.  Eigen lives only in this file.

#include <decimate/reduce/seam_aware_reduce.h>

#include <decimate/reduce/corner_mesh.h>
#include <decimate/types/mesh_data.h>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace decimate {

namespace {

using corner::He8;
using corner::HeIndex;
using corner::kNoHE;
using math::Vec2f;
using math::Vec3f;
using Mat6 = Eigen::Matrix<double, 6, 6>;
using Vec2d = Eigen::Vector2d;
using Vec5 = Eigen::Matrix<double, 5, 1>;

/// Cost of a candidate that cannot be collapsed.
constexpr double kInf = std::numeric_limits<double>::infinity();

/// Collinearity (|sin|) and edge-length-ratio tolerance of an unifiable edge pair.
constexpr double kUnifyTol = 1e-5;

/// Read-only view of the two corner tables and the UV positions (z = 0).
struct View {
    const std::vector<He8>   &s;  ///< Surface table.
    const std::vector<He8>   &u;  ///< UV table; half-edge h is the same corner as in `s`.
    const std::vector<Vec3f> &uv; ///< UV vertex positions as (u, v, 0).
};


/// What kind of edge a surface half-edge is, by the twins it has in the surface and UV tables.
enum class Kind {
    Interior, ///< Surface twin and UV twin exist.
    Seam,     ///< Surface twin exists, UV twin does not (either side).
    Boundary, ///< No surface twin: a 3D boundary edge.
};


/// The boundary half-edges meeting at a vertex: those ending at it and those leaving it.
struct BoundaryLinks {
    HeIndex in = kNoHE;  ///< Last boundary half-edge ending at the vertex.
    HeIndex out = kNoHE; ///< Last boundary half-edge leaving the vertex.
    int     nIn = 0;     ///< How many end at it.
    int     nOut = 0;    ///< How many leave it.
};


/// One UV copy of an edge with its two neighbours along the UV boundary, in chain order:
/// c - a - b - g, where (a, b) is the edge (a = copy of the edge's origin p1).
struct Side {
    HeIndex a = 0;        ///< Origin copy of the edge.
    HeIndex b = 0;        ///< Tip copy of the edge.
    HeIndex c = 0;        ///< Far end of the adjacent boundary edge at `a` (valid iff `loOk`).
    HeIndex g = 0;        ///< Far end of the adjacent boundary edge at `b` (valid iff `hiOk`).
    HeIndex dLo = kNoHE;  ///< The adjacent boundary half-edge at `a`.
    HeIndex dHi = kNoHE;  ///< The adjacent boundary half-edge at `b`.
    bool    loOk = false; ///< `a` has exactly one in- and one out-boundary edge.
    bool    hiOk = false; ///< Same for `b`.
};


/// Two consecutive UV edges p->q, q->r: collinear in the same direction, and their
/// length ratio |qr| / |pq|.
struct Pair {
    bool   collinear = false;
    double ratio = 0.0;
};


/// An edge together with its UV copies and which ends may move ("free", paper fig. 9).
struct Chain {
    std::array<Side, 2> side;
    int                 count = 1; ///< 2 for a seam (two copies), 1 for a 3D boundary edge.
    bool                freeLo = false;
    bool                freeHi = false;
};


/// The two corner tables of a mesh plus its UV positions.
struct Tables {
    std::vector<He8>   s;
    std::vector<He8>   u;
    std::vector<Vec3f> uv;
};


/// One UV copy's contribution: its quadric and the line its new UV is confined to,
/// `uv = base + t * dir` (dir = 0 for a fixed UV).
struct SideTerm {
    Mat6  q;
    Vec2d base;
    Vec2d dir;
};


/// Result of the side-sum minimization: shared position x, line parameter t, and the cost there.
struct Solved {
    Eigen::Vector3d x = Eigen::Vector3d::Zero();
    double          t = 0.0;
    double          cost = kInf;
};


/// Free 5D placement (x, y, z, u, v) of a non-silhouette interior edge.
struct Solved5 {
    Eigen::Vector3d x = Eigen::Vector3d::Zero();
    Vec2d           uv = Vec2d::Zero();
    double          cost = kInf;
};


/// The cached cost and placement of one candidate edge.
struct Placement {
    double               cost = kInf;
    Kind                 kind = Kind::Interior;
    Vec3f                xyz;
    std::array<Vec2f, 2> w; ///< New UV of the surviving vertex of each UV copy (A, B).
};


/// Everything a collapse mutates, plus the read-only views of it.
struct State {
    std::vector<He8>   s;
    std::vector<He8>   u;
    std::vector<Vec3f> pos;
    std::vector<Vec3f> uv;
    std::vector<Mat6>  quad;
    std::vector<char>  deadFace;
    std::vector<char>  deadV;
    [[nodiscard]] View view() const {
        return { s, u, uv };
    }
};


/// A heap candidate: the cost cached for half-edge `he`; stale ones are re-queued at pop.
struct Cand {
    double  cost = 0.0;
    HeIndex he = 0;
};


/// Orders the heap so the lowest cost pops first.
struct CostGreater {
    bool operator()(const Cand &a, const Cand &b) const noexcept {
        return a.cost > b.cost;
    }
};


using CandHeap = std::priority_queue<Cand, std::vector<Cand>, CostGreater>;

} // namespace


/// CCW predecessor of half-edge @p h in its face loop.
static HeIndex prevOf(HeIndex h) noexcept {
    return corner::nextA(corner::nextA(h));
}


/// Widens a half-edge/vertex id to a container index.
static std::size_t idx(HeIndex i) noexcept {
    return static_cast<std::size_t>(i);
}


/// Classifies surface half-edge @p h by which of its twins exist.
static Kind kindOf(const View &v, HeIndex h) noexcept {
    if (v.s[h].sym == kNoHE) {
        return Kind::Boundary;
    }
    return v.u[h].sym != kNoHE ? Kind::Interior : Kind::Seam;
}


/// Calls @p visit once for every half-edge ending at the vertex @p incoming ends at.
/// (`corner::forEachIncident` visits `incoming` twice on a boundary fan.)
template <class F>
static void forEachIncidentOnce(const std::vector<He8> &he, HeIndex incoming, F visit) {
    bool first = true;
    corner::forEachIncident(he, incoming, [&](HeIndex x) {
        if (x == incoming && !first) {
            return;
        }
        first = x != incoming && first;
        visit(x);
    });
}


/// True iff the vertex @p incoming ends at touches a boundary half-edge of @p he.
static bool onBoundary(const std::vector<He8> &he, HeIndex incoming) {
    bool found = false;
    forEachIncidentOnce(he, incoming, [&](HeIndex x) {
        found = found || he[x].sym == kNoHE || he[corner::nextA(x)].sym == kNoHE;
    });
    return found;
}


/// Collects the boundary half-edges around the vertex @p incoming ends at.
static BoundaryLinks boundaryLinks(const std::vector<He8> &he, HeIndex incoming) {
    BoundaryLinks l;
    forEachIncidentOnce(he, incoming, [&](HeIndex x) {
        const HeIndex o = corner::nextA(x);
        if (he[x].sym == kNoHE) {
            l.in = x;
            ++l.nIn;
        }
        if (he[o].sym == kNoHE) {
            l.out = o;
            ++l.nOut;
        }
    });
    return l;
}


/// UV position of UV vertex @p q, widened to double for the solver.
static Vec2d uvOf(const View &v, HeIndex q) noexcept {
    return { static_cast<double>(v.uv[q].x), static_cast<double>(v.uv[q].y) };
}


/// Chain-order description of UV half-edge @p h.  @p forward: `h` runs a -> b; otherwise
/// `h` runs b -> a (the second copy of a seam).
static Side makeSide(const View &v, HeIndex h, bool forward) {
    Side sd;
    sd.a = forward ? corner::originOf(v.u, h) : corner::tipOf(v.u, h);
    sd.b = forward ? corner::tipOf(v.u, h) : corner::originOf(v.u, h);

    // The boundary neighbours at both ends of the edge, in chain order.
    const BoundaryLinks atA = boundaryLinks(v.u, forward ? prevOf(h) : h);
    const BoundaryLinks atB = boundaryLinks(v.u, forward ? h : prevOf(h));
    sd.loOk = atA.nIn == 1 && atA.nOut == 1;
    sd.hiOk = atB.nIn == 1 && atB.nOut == 1;
    sd.dLo = forward ? atA.in : atA.out;
    sd.dHi = forward ? atB.out : atB.in;

    // The far ends are only meaningful when the end has exactly one neighbour on each side.
    if (sd.loOk) {
        sd.c = forward ? corner::originOf(v.u, sd.dLo) : corner::tipOf(v.u, sd.dLo);
    }
    if (sd.hiOk) {
        sd.g = forward ? corner::tipOf(v.u, sd.dHi) : corner::originOf(v.u, sd.dHi);
    }
    return sd;
}


/// Collinearity and length ratio of the edge pair p->q, q->r.
static Pair pairOf(const Vec2d &p, const Vec2d &q, const Vec2d &r) noexcept {
    const Vec2d  d = q - p;
    const Vec2d  e = r - q;
    const double len = d.norm() * e.norm();
    Pair         out;
    if (len <= 0.0) {
        return out; // a zero-length edge has no direction to compare
    }

    const double cross = (d.x() * e.y()) - (d.y() * e.x());
    out.collinear = std::abs(cross) <= kUnifyTol * len && d.dot(e) > 0.0;
    out.ratio = e.norm() / d.norm();
    return out;
}


/// The adjacent-edge pair at the lo end (c,a,b) or the hi end (a,b,g) of a side.
static Pair endPair(const View &v, const Side &sd, bool lo) noexcept {
    return lo ? pairOf(uvOf(v, sd.c), uvOf(v, sd.a), uvOf(v, sd.b))
              : pairOf(uvOf(v, sd.a), uvOf(v, sd.b), uvOf(v, sd.g));
}


/// Do the two edge-length ratios agree within kUnifyTol (relative)?
static bool ratiosMatch(double r0, double r1) noexcept {
    return std::abs(r0 - r1) <= kUnifyTol * std::max(r0, r1);
}


/// Is the end of the edge unifiable with its neighbour on every copy?
static bool endFree(const View &v, const Chain &ch, bool lo) {
    const Side   &s0 = ch.side[0];
    const bool    ok0 = lo ? s0.loOk : s0.hiOk;
    const HeIndex d0 = lo ? s0.dLo : s0.dHi;
    const Pair    p0 = ok0 ? endPair(v, s0, lo) : Pair {};
    if (!ok0 || !p0.collinear) {
        return false;
    }

    // A 3D boundary edge has one copy: the neighbour must be a boundary edge too.
    if (ch.count == 1) {
        return v.s[d0].sym == kNoHE;
    }

    // A seam has two copies: both neighbours must be the two sides of one seam edge, and must
    // continue the line with the same length ratio, so both copies keep a common parametrization.
    const Side   &s1 = ch.side[1];
    const HeIndex d1 = lo ? s1.dLo : s1.dHi;
    if (!(lo ? s1.loOk : s1.hiOk) || v.s[d0].sym != d1) {
        return false; // neighbours are not the two sides of one seam edge
    }
    const Pair p1 = endPair(v, s1, lo);
    return p1.collinear && ratiosMatch(p0.ratio, p1.ratio);
}


/// Builds the chain of @p h (one copy for a boundary edge, two for a seam) and tests both ends.
static Chain analyse(const View &v, HeIndex h, Kind kind) {
    Chain ch;
    ch.side[0] = makeSide(v, h, true);
    if (kind == Kind::Seam) {
        ch.count = 2;
        ch.side[1] = makeSide(v, v.s[h].sym, false);
    }
    ch.freeLo = endFree(v, ch, true);
    ch.freeHi = endFree(v, ch, false);
    return ch;
}


/// Throws std::invalid_argument unless @p tm satisfies the `TexturedMesh` contract this file relies on.
static void validate(const TexturedMesh &tm) {
    if (tm.faceV.size() != tm.faceU.size() || tm.faceV.size() % 3 != 0) {
        throw std::invalid_argument("TexturedMesh: faceV/faceU size mismatch");
    }

    // Each UV vertex must belong to exactly one surface vertex: the flat per-UV quadric array relies on it.
    std::vector<int> owner(tm.uvs.size(), -1);
    for (std::size_t i = 0; i < tm.faceV.size(); ++i) {
        const int sv = tm.faceV[i];
        const int uvx = tm.faceU[i];
        if (sv < 0 || uvx < 0 || idx(static_cast<HeIndex>(sv)) >= tm.positions.size() ||
            idx(static_cast<HeIndex>(uvx)) >= tm.uvs.size()) {
            throw std::invalid_argument("TexturedMesh: index out of range");
        }
        int &o = owner[static_cast<std::size_t>(uvx)];
        if (o >= 0 && o != sv) {
            throw std::invalid_argument("TexturedMesh: a UV vertex is shared by two vertices");
        }
        o = sv;
    }
}


/// Builds the surface and UV corner tables of @p tm; half-edge h names the same corner in both.
static Tables buildTables(const TexturedMesh &tm) {
    Tables   t;
    MeshData surface;
    surface.indices = tm.faceV;
    MeshData uvMesh;
    uvMesh.indices = tm.faceU;
    t.s = corner::build(surface);
    t.u = corner::build(uvMesh);
    t.uv.reserve(tm.uvs.size());
    for (const Vec2f &p : tm.uvs) {
        t.uv.emplace_back(p.x, p.y, 0.0f);
    }
    return t;
}


/// Counts the seam edges of the live faces that have an unifiable end vs none.
static SeamCollapsibility classifyLive(const View &v, const std::vector<char> &deadFace) {
    SeamCollapsibility out;
    for (std::size_t i = 0; i < v.s.size(); ++i) {
        const HeIndex h = static_cast<HeIndex>(i);
        if (deadFace[i / 3] != 0 || kindOf(v, h) != Kind::Seam || h > v.s[h].sym) {
            continue; // dead, not a seam, or the lower-indexed twin counts this seam
        }

        const Chain ch = analyse(v, h, Kind::Seam);
        ++(ch.freeLo || ch.freeHi ? out.collapsible : out.uncollapsible);
    }
    return out;
}


/// Garland-Heckbert 1998 quadric of the 5D triangle (p, q, r), weighted by @p weight, in
/// the homogeneous 6x6 form [[A, b], [b^T, c]] so that cost(x) = (x,1)^T Q (x,1).
static Mat6 faceQuadric(const Vec5 &p, const Vec5 &q, const Vec5 &r, double weight) {
    Mat6         out = Mat6::Zero();
    const Vec5   d1 = q - p;
    const double l1 = d1.norm();
    if (l1 <= 0.0 || weight <= 0.0) {
        return out; // degenerate triangle: no plane to measure distance to
    }

    // Orthonormal frame (e1, e2) of the triangle's 2D plane inside 5D, by Gram-Schmidt.
    const Vec5   e1 = d1 / l1;
    const Vec5   t = (r - p) - (e1.dot(r - p) * e1);
    const double l2 = t.norm();
    if (l2 <= 0.0) {
        return out;
    }
    const Vec5 e2 = t / l2;

    // Squared distance to the plane = |x - p|^2 minus its in-plane part, expanded into A, b, c.
    const double pe1 = p.dot(e1);
    const double pe2 = p.dot(e2);
    const Vec5   b = (pe1 * e1) + (pe2 * e2) - p;
    out.topLeftCorner<5, 5>() = Eigen::Matrix<double, 5, 5>::Identity() - (e1 * e1.transpose()) -
                                (e2 * e2.transpose());
    out.topRightCorner<5, 1>() = b;
    out.bottomLeftCorner<1, 5>() = b.transpose();
    out(5, 5) = p.squaredNorm() - (pe1 * pe1) - (pe2 * pe2);
    return out * weight;
}


/// The 5D point (x, y, z, u, v) of one face corner.
static Vec5 wedgePoint(const TexturedMesh &tm, std::size_t corner) {
    const Vec3f &p = tm.positions[idx(static_cast<HeIndex>(tm.faceV[corner]))];
    const Vec2f &t = tm.uvs[idx(static_cast<HeIndex>(tm.faceU[corner]))];
    return Vec5 { p.x, p.y, p.z, t.x, t.y };
}


/// One quadric per UV vertex.  Every UV vertex belongs to exactly one surface vertex
/// (validated), so the paper's hash map keyed by the (vertex, UV) pair is this flat array.
static std::vector<Mat6> buildQuadrics(const TexturedMesh &tm) {
    std::vector<Mat6> quad(tm.uvs.size(), Mat6::Zero());

    // Each face contributes its area-weighted quadric to its three UV vertices.
    for (std::size_t f = 0; f < tm.faceCount(); ++f) {
        const std::size_t     c = 3 * f;
        const Vec5            p = wedgePoint(tm, c);
        const Vec5            q = wedgePoint(tm, c + 1);
        const Vec5            r = wedgePoint(tm, c + 2);
        const Eigen::Vector3d e1 = p.head<3>() - q.head<3>();
        const Eigen::Vector3d e2 = p.head<3>() - r.head<3>();
        const Mat6            fq = faceQuadric(p, q, r, 0.5 * e1.cross(e2).norm());
        for (std::size_t k = 0; k < 3; ++k) {
            quad[idx(static_cast<HeIndex>(tm.faceU[c + k]))] += fq;
        }
    }
    return quad;
}


/// Minimizer of y^T H y + 2 g^T y by (H + lambda I) y = -g + lambda y0.  The tiny ridge
/// keeps a rank-deficient (flat) quadric well posed and pulls the free directions to y0.
static Eigen::VectorXd ridgeSolve(const Eigen::MatrixXd &hess, const Eigen::VectorXd &lin,
                                  const Eigen::VectorXd &y0) {
    const auto            k = hess.rows();
    const double          scale = hess.trace() / static_cast<double>(k);
    const double          lambda = scale > 0.0 ? 1e-9 * scale : 1e-12;
    const Eigen::MatrixXd a = hess + (lambda * Eigen::MatrixXd::Identity(k, k));
    return a.ldlt().solve((-lin) + (lambda * y0));
}


/// Minimizes the sum of the side costs over x (shared by all sides) and, when @p withT,
/// the common line parameter t in [0,1].  Unknowns (x, y, z, t, 1); every equality of the
/// paper's QP is the substitution, so what is left is one box constraint on t, resolved
/// by clamping t and re-minimizing x (exact for a convex quadratic with one boxed variable).
static Solved minimizeSides(const std::vector<SideTerm> &terms, bool withT, const Eigen::Vector3d &x0) {
    // Substitute each side's line into its 6x6 quadric, giving one 5x5 form over (x, y, z, t, 1).
    Eigen::Matrix<double, 5, 5> m = Eigen::Matrix<double, 5, 5>::Zero();
    for (const SideTerm &term : terms) {
        Eigen::Matrix<double, 6, 5> p = Eigen::Matrix<double, 6, 5>::Zero();
        p.topLeftCorner<3, 3>().setIdentity();
        p.block<2, 1>(3, 3) = term.dir;
        p.block<2, 1>(3, 4) = term.base;
        p(5, 4) = 1.0;
        m += p.transpose() * term.q * p;
    }

    const Eigen::Index    k = withT ? 4 : 3;
    const Eigen::MatrixXd hess = m.topLeftCorner(k, k);
    const Eigen::VectorXd lin = m.block(0, 4, k, 1);
    Eigen::VectorXd       y0(k);
    y0.head(3) = x0;
    if (withT) {
        y0(3) = 0.5;
    }
    Eigen::VectorXd y = ridgeSolve(hess, lin, y0);

    // t outside [0,1]: pin it at the nearest end and re-solve for x alone.
    if (withT && (y(3) < 0.0 || y(3) > 1.0)) {
        const double          t = std::clamp(y(3), 0.0, 1.0);
        const Eigen::VectorXd lin3 = lin.head(3) + (hess.block(0, 3, 3, 1) * t);
        y.head(3) = ridgeSolve(hess.topLeftCorner(3, 3), lin3, x0);
        y(3) = t;
    }

    Solved out;
    out.x = y.head(3);
    out.t = withT ? y(3) : 0.0;
    const double cost = y.dot(hess * y) + (2.0 * lin.dot(y)) + m(4, 4);
    out.cost = std::isfinite(cost) ? std::max(cost, 0.0) : kInf;
    return out;
}


/// Minimizer of the unconstrained 5D quadric @p q (a non-silhouette interior edge), regularized toward @p y0.
static Solved5 minimizeFree(const Mat6 &q, const Vec5 &y0) {
    const Eigen::MatrixXd hess = q.topLeftCorner<5, 5>();
    const Eigen::VectorXd lin = q.topRightCorner<5, 1>();
    const Eigen::VectorXd y = ridgeSolve(hess, lin, y0);
    Solved5               out;
    out.x = y.head(3);
    out.uv = y.tail(2);
    const double cost = y.dot(hess * y) + (2.0 * lin.dot(y)) + q(5, 5);
    out.cost = std::isfinite(cost) ? std::max(cost, 0.0) : kInf;
    return out;
}


/// Narrows a solver UV to float.
static Vec2f toVec2f(const Vec2d &p) {
    return { static_cast<float>(p.x()), static_cast<float>(p.y()) };
}


/// Widens a position to double for the solver.
static Eigen::Vector3d toEigen(const Vec3f &p) {
    return { p.x, p.y, p.z };
}


/// Narrows a solver position to float.
static Vec3f toVec3f(const Eigen::Vector3d &p) {
    return { static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()) };
}


/// Cost and placement of an interior edge (surface and UV twins both exist).
static Placement evalInterior(const State &st, HeIndex h) {
    const HeIndex q1 = corner::originOf(st.u, h);
    const HeIndex q2 = corner::tipOf(st.u, h);
    const bool    sil1 = onBoundary(st.u, prevOf(h));
    const bool    sil2 = onBoundary(st.u, h);
    Placement     pl;
    if (sil1 && sil2) {
        return pl; // paper 5.1: an interior edge between two silhouette vertices
    }

    // Both vertices' quadrics together, starting from the edge midpoint.
    const Mat6            q = st.quad[q1] + st.quad[q2];
    const Eigen::Vector3d x0 =
        0.5 * (toEigen(st.pos[corner::originOf(st.s, h)]) + toEigen(st.pos[corner::tipOf(st.s, h)]));
    // No silhouette end: nothing pins the UV, so place the full 5D optimum.
    if (!sil1 && !sil2) {
        const Vec2d mid = 0.5 * (uvOf(st.view(), q1) + uvOf(st.view(), q2));
        Vec5        y0;
        y0 << x0, mid;
        const Solved5 sol = minimizeFree(q, y0);
        pl.cost = sol.cost;
        pl.xyz = toVec3f(sol.x);
        pl.w[0] = toVec2f(sol.uv);
        return pl;
    }

    // One silhouette end: its UV is the only legal one, so only x moves.
    const HeIndex               fixedUV = sil1 ? q1 : q2;
    const std::vector<SideTerm> terms { { q, uvOf(st.view(), fixedUV), Vec2d::Zero() } };
    const Solved                sol = minimizeSides(terms, false, x0);
    pl.cost = sol.cost;
    pl.xyz = toVec3f(sol.x);
    pl.w[0] = { st.uv[fixedUV].x, st.uv[fixedUV].y }; // bit-exact copy of the silhouette UV
    return pl;
}


/// UV that copy @p sd gets: on the line when both ends slide, else the fixed end's UV.
static SideTerm termOf(const State &st, const Side &sd, const Chain &ch) {
    const View v = st.view();
    SideTerm   t { st.quad[sd.a] + st.quad[sd.b], Vec2d::Zero(), Vec2d::Zero() };
    if (ch.freeLo && ch.freeHi) {
        t.base = uvOf(v, sd.c);
        t.dir = uvOf(v, sd.g) - uvOf(v, sd.c);
    } else {
        t.base = uvOf(v, ch.freeLo ? sd.b : sd.a);
    }
    return t;
}


/// Cost and placement of a seam or 3D-boundary edge: the UVs slide along the existing straight line.
static Placement evalChain(const State &st, HeIndex h, Kind kind) {
    const Chain ch = analyse(st.view(), h, kind);
    Placement   pl;
    pl.kind = kind;
    if (!ch.freeLo && !ch.freeHi) {
        return pl; // neither end unifiable: never collapsible
    }

    const bool            slide = ch.freeLo && ch.freeHi;
    std::vector<SideTerm> terms;
    for (int i = 0; i < ch.count; ++i) {
        terms.push_back(termOf(st, ch.side[static_cast<std::size_t>(i)], ch));
    }
    const Eigen::Vector3d x0 =
        0.5 * (toEigen(st.pos[corner::originOf(st.s, h)]) + toEigen(st.pos[corner::tipOf(st.s, h)]));
    const Solved sol = minimizeSides(terms, slide, x0);
    pl.cost = sol.cost;
    pl.xyz = toVec3f(sol.x);

    // Every copy lands on its own line at the common t (or on its fixed end's UV).
    for (std::size_t i = 0; i < terms.size(); ++i) {
        const Side   &sd = ch.side[i];
        const HeIndex fixedVertex = ch.freeLo ? sd.b : sd.a;
        pl.w[i] = slide ? toVec2f(terms[i].base + (sol.t * terms[i].dir))
                        : Vec2f { st.uv[fixedVertex].x, st.uv[fixedVertex].y };
    }
    return pl;
}


/// Cost and placement of the collapse of surface half-edge @p h.
static Placement evaluate(const State &st, HeIndex h) {
    const Kind kind = kindOf(st.view(), h);
    return kind == Kind::Interior ? evalInterior(st, h) : evalChain(st, h, kind);
}


/// UV-table guard for the collapse of UV half-edge @p hu to @p w.
static bool uvLegal(const State &st, HeIndex hu, const Vec2f &w) {
    return corner::collapseLegal(st.u, st.uv, hu, corner::originOf(st.u, hu),
                                 corner::tipOf(st.u, hu), Vec3f { w.x, w.y, 0.0f });
}


/// Every collapse goes through corner::collapseLegal: on the surface, and on the UV table
/// for each UV half-edge it removes.  No separate pinch test is needed for seams: a seam
/// edge is only a candidate when an end is unifiable, and a unifiable end has two seam
/// edges whose two UV copies both span the whole vertex fan (`endFree`), i.e. the vertex is
/// interior, so the edge never joins two 3D-boundary vertices.
static bool legal(const State &st, HeIndex h, const Placement &pl) {
    const HeIndex p1 = corner::originOf(st.s, h);
    const HeIndex p2 = corner::tipOf(st.s, h);
    if (!corner::collapseLegal(st.s, st.pos, h, p1, p2, pl.xyz)) {
        return false;
    }
    if (pl.kind != Kind::Seam) {
        return uvLegal(st, h, pl.w[0]);
    }
    return uvLegal(st, h, pl.w[0]) && uvLegal(st, st.s[h].sym, pl.w[1]);
}


/// Collapses UV half-edge @p hu: its tip merges into its origin, which takes @p w.
static void collapseUV(State &st, HeIndex hu, const Vec2f &w) {
    const HeIndex keep = corner::originOf(st.u, hu);
    const HeIndex gone = corner::tipOf(st.u, hu);
    corner::collapse(st.u, hu, keep, st.deadFace);
    st.uv[keep] = Vec3f { w.x, w.y, 0.0f };
    st.quad[keep] += st.quad[gone];
}


/// Performs the collapse of @p h with placement @p pl on the surface table and both UV copies.
static void applyCollapse(State &st, HeIndex h, const Placement &pl) {
    const HeIndex p1 = corner::originOf(st.s, h);
    const HeIndex p2 = corner::tipOf(st.s, h);
    const HeIndex twin = st.s[h].sym;

    collapseUV(st, h, pl.w[0]);
    if (pl.kind == Kind::Seam) {
        collapseUV(st, twin, pl.w[1]);
    }

    corner::collapse(st.s, h, p1, st.deadFace);
    st.pos[p1] = pl.xyz;
    st.deadV[p2] = 1;
}


/// Builds the mutable collapse state of @p tm; vertices no face uses start dead.
static State makeState(const TexturedMesh &tm) {
    Tables t = buildTables(tm);
    State  st;
    st.s = std::move(t.s);
    st.u = std::move(t.u);
    st.uv = std::move(t.uv);
    st.pos = tm.positions;
    st.quad = buildQuadrics(tm);

    st.deadFace.assign(tm.faceCount(), 0);
    st.deadV.assign(tm.positions.size(), 1);
    for (const int v : tm.faceV) {
        st.deadV[static_cast<std::size_t>(v)] = 0; // phantom vertices start dead
    }
    return st;
}


/// Queues every undirected surface edge with its initial cost.
static CandHeap seedHeap(const State &st) {
    CandHeap heap;
    for (std::size_t i = 0; i < st.s.size(); ++i) {
        const HeIndex h = static_cast<HeIndex>(i);
        const HeIndex sym = st.s[h].sym;
        if (sym == kNoHE || h < sym) { // one candidate per undirected edge
            heap.push({ evaluate(st, h).cost, h });
        }
    }
    return heap;
}


/// Pops until one legal, current-cost collapse is found and applies it.  Returns -1 when
/// the heap is exhausted, else 1 if the collapse was a seam collapse and 0 otherwise.
static int collapseNext(State &st, CandHeap &heap) {
    while (!heap.empty()) {
        const Cand c = heap.top();
        heap.pop();
        const HeIndex h = c.he;
        if (st.deadFace[h / 3] != 0) {
            continue; // the edge's face is already collapsed away
        }

        // A neighbour collapse may have changed this edge's cost since it was queued.
        const Placement cur = evaluate(st, h);
        if (cur.cost != c.cost) {
            heap.push({ cur.cost, h }); // stale priority: re-queue with the current cost
            continue;
        }
        if (std::isinf(cur.cost) || !legal(st, h, cur)) {
            continue; // not collapsible now; the edge is dropped, not re-queued
        }

        applyCollapse(st, h, cur);
        return cur.kind == Kind::Seam ? 1 : 0;
    }
    return -1;
}


/// Numbers the entries of @p id marked 0 (used) 0, 1, 2, ... in ascending order and appends
/// `convert(src[i])` to @p dst for each; unused entries stay -1.
template <class Src, class Dst, class Convert>
static void compactUsed(std::vector<int> &id, const std::vector<Src> &src, std::vector<Dst> &dst,
                        Convert convert) {
    for (std::size_t i = 0; i < id.size(); ++i) {
        if (id[i] != 0) {
            continue;
        }
        id[i] = static_cast<int>(dst.size());
        dst.push_back(convert(src[i]));
    }
}


/// Compacts the live faces into a fresh TexturedMesh (vertices and UVs in ascending id
/// order, as `corner::reconstruct`), faces in their original corner order.
static TexturedMesh assemble(const State &st) {
    std::vector<int>     newV(st.pos.size(), -1);
    std::vector<int>     newU(st.uv.size(), -1);
    std::vector<HeIndex> liveBase; // first half-edge of every live face
    for (std::size_t f = 0; f < st.deadFace.size(); ++f) {
        if (st.deadFace[f] == 0) {
            liveBase.push_back(static_cast<HeIndex>(3 * f));
        }
    }

    // Mark the vertices and UVs that live faces still use (0 = used, see compactUsed).
    for (const HeIndex base : liveBase) {
        const auto sv = corner::faceVertices(st.s, base);
        const auto uvx = corner::faceVertices(st.u, base);
        for (std::size_t k = 0; k < 3; ++k) {
            newV[sv[k]] = 0;
            newU[uvx[k]] = 0;
        }
    }

    TexturedMesh out;
    compactUsed(newV, st.pos, out.positions, [](const Vec3f &p) { return p; });
    compactUsed(newU, st.uv, out.uvs, [](const Vec3f &p) { return Vec2f { p.x, p.y }; });

    // Faces keep their original corner order, re-indexed to the compacted vertices.
    for (const HeIndex base : liveBase) {
        const auto sv = corner::faceVertices(st.s, base);
        const auto uvx = corner::faceVertices(st.u, base);
        for (std::size_t k = 0; k < 3; ++k) {
            out.faceV.push_back(newV[sv[k]]);
            out.faceU.push_back(newU[uvx[k]]);
        }
    }
    return out;
}


/// Total seam edges in a classification.
static int countSeams(const SeamCollapsibility &c) noexcept {
    return c.collapsible + c.uncollapsible;
}


SeamCollapsibility classifySeamEdges(const TexturedMesh &tm) {
    validate(tm);
    const Tables            t = buildTables(tm);
    const std::vector<char> noneDead(tm.faceCount(), 0);
    return classifyLive({ t.s, t.u, t.uv }, noneDead);
}


SeamAwareResult reduceSeamAware(const TexturedMesh &in, int targetVertexCount) {
    const SeamCollapsibility before = classifySeamEdges(in); // validates
    SeamAwareResult          res;
    res.stats.seamEdgesBefore = countSeams(before);

    State       st = makeState(in);
    std::size_t numValid = 0;
    for (const char d : st.deadV) {
        numValid += d == 0 ? 1 : 0;
    }

    // Nothing to do: hand the input back untouched, counters at zero.
    if (static_cast<std::size_t>(std::max(targetVertexCount, 0)) >= numValid ||
        in.faceCount() == 0) {
        res.mesh = in;
        res.stats.seamEdgesAfter = res.stats.seamEdgesBefore;
        return res;
    }

    res.stats.collapsibleSeamEdges = before.collapsible;
    res.stats.uncollapsibleSeamEdges = before.uncollapsible;
    // Greedy loop: cheapest legal edge first, until the target is met or no legal collapse remains.
    CandHeap heap = seedHeap(st);
    while (numValid > static_cast<std::size_t>(std::max(targetVertexCount, 0))) {
        const int step = collapseNext(st, heap);
        if (step < 0) {
            break;
        }
        res.stats.collapsedSeamEdges += step;
        --numValid;
    }

    res.mesh = assemble(st);
    res.stats.seamEdgesAfter = countSeams(classifySeamEdges(res.mesh));
    return res;
}

} // namespace decimate
