// This file is part of the decimate library; see decimate.h for version/license details.
/// @file conformal_mapper.cpp
/// @brief LSCM chart flattening — a port of 2003 `CConformalMapper` (`legacy2003/chartification/cconformalmapper.cpp`)
/// with a different solver (sparse Cholesky for CG) and a robustness layer; UVs equal 2003 only to the golden tolerance.
///
/// ## Algorithm (Lévy et al. 2002 — Least-Squares Conformal Maps)
///
/// For each chart:
///   1. Build the chart's local vertex set by CORNER, merging corners only across non-seam edges
///      (2003 `CreateVertexStructure` after `CutBoundaryPath`,
///      `legacy2003/chartification/cchartification.cpp:395-407`), so a chart that touches itself is cut into a disk.
///   2. Identify boundary vertices; pin the farthest-apart pair to fixed UVs.
///   3. For each triangle, compute complex Cauchy-Riemann gradients via
///      Gram-Schmidt orthonormalization in the triangle's 3D tangent plane.
///   4. Directly assemble the normal equations AᵀA x = Aᵀb for free
///      (unpinned) vertices and solve them by sparse Cholesky (Eigen SimplicialLDLT).
///      2003 solved by CG, tol 1e-12, cap 15000 (legacy2003/chartification/cconformalmapper.cpp:306);
///      on 6k-face charts that took seconds — the direct solve is exact and takes milliseconds.
///   5. Back-substitute to fill HalfEdgeChart::uv for every local vertex.
///
/// ## Operator mapping from legacy Vec3 to math::Vec3
///
///   Legacy `(ba) * e1` — Vec3::operator*(Vec3) was the dot product.
///   → `ba.dot(e1)` — math::Vec3::operator*(Vec3) is component-wise; must use .dot().
///
///   Legacy Gram-Schmidt: `w2 = w2 - (w1 * (w2 * w1))`
///   Here `w2 * w1` is a dot product (scalar), then `w1 * scalar` scales w1.
///   → `e2 = ca - e1 * e1.dot(ca)` (extract dot as a scalar, then scalar-multiply e1).
///
///   No legacy `^` (cross product) is used in this file; no mapping needed for cross.
///
/// ## Pinning (legacy cconformalmapper.cpp:48-193)
///
///   Diameter() (lines 48–82) finds the pair of boundary vertices (start of edge i,
///   end of edge j) with the maximum 3D Euclidean distance.
///
///   InitPinned() (lines 135–193) swaps those two vertices to the last two positions
///   in the chart-local vertex list and assigns:
///     vertices[size-2].u = 0.0f, vertices[size-2].v = 0.0f   → UV (0, 0)
///     vertices[size-1].u = 1.0f, vertices[size-1].v = 1.0f   → UV (1, 1)
///
///   The modern port avoids physically swapping vertices; it simply records
///   pin0Local (→ UV (0,0)) and pin1Local (→ UV (1,1)) and handles them
///   explicitly in the assembly and back-substitution steps.
///
/// ## Solution vector layout
///
///   x = [u_free0, …, u_free{N-1}, v_free0, …, v_free{N-1}]  (u-block then v-block)
///
///   Mirrors 2003 CopyResult (legacy2003/chartification/cconformalmapper.cpp:35-42):
///     vertices[offset + i].u = (float) vec[i];
///     vertices[offset + i].v = (float) vec[num + i];
///
/// Boundary contract: no GL/OS/I/O/global state; C++20; RAII throughout.

#include <decimate/chart/parametrize.h>

#include "chart_quality.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/SparseCholesky>

#include <decimate/math/vec2.h>
#include <decimate/math/vec3.h>
#include <decimate/chart/chart.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>

namespace decimate::chart {

namespace {

using UV = math::Vec2f;


/// One entry in a sparse matrix row: (column, value).
struct RowEntry {
    int    col;
    double val;
};


/// Union-find over small integer ids; the root is always the smallest id of its set.
class DisjointSet {
public:
    explicit DisjointSet(std::size_t n)
        : parent_(n) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    [[nodiscard]] int find(int x) {
        while (parent_[static_cast<std::size_t>(x)] != x) {
            int &p = parent_[static_cast<std::size_t>(x)];
            p = parent_[static_cast<std::size_t>(p)]; // path halving
            x = p;
        }
        return x;
    }

    void unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a != b) {
            parent_[static_cast<std::size_t>(std::max(a, b))] = std::min(a, b);
        }
    }

private:
    std::vector<int> parent_;
};


/// The chart's local triangulation: vertices by corner, cut along its seams.
struct CornerMesh {
    std::vector<int> localToMesh; ///< local vertex -> mesh vertex id (many-to-one across seams).
    std::vector<int> localTri;    ///< 3 local vertices (corners a, b, c) per chart face.
};


/// Pinning state: local indices + UV values for the two pinned vertices.
struct PinState {
    int    pin0Local;   ///< Local index of pin0 (→ UV 0,0).
    int    pin1Local;   ///< Local index of pin1 (→ UV 1,1).
    int    numFree;     ///< Number of free (unpinned) local vertices.
    double uPin0 = 0.0; ///< u for pin0.
    double vPin0 = 0.0; ///< v for pin0.
    double uPin1 = 1.0; ///< u for pin1.
    double vPin1 = 1.0; ///< v for pin1.
};

} // namespace


/// 3D position of a mesh vertex.
[[nodiscard]] static math::Vec3f vertexPos(const MeshData &data, int vid) noexcept {
    return data.positions[static_cast<std::size_t>(vid)];
}


/// Accumulate one sparse row's outer product into AtA and Atb:
/// AtA[ea.col][eb.col] += ea.val * eb.val for every entry pair, atb[ea.col] += ea.val * rhs for every entry.
static void accumRow(std::vector<Eigen::Triplet<double>> &ata,
                     std::vector<double>                 &atb,
                     const std::vector<RowEntry>         &row,
                     double                               rhs) {
    // every entry pair feeds the matrix; the right-hand side takes one term per entry
    for (const auto &ea : row) {
        atb[static_cast<std::size_t>(ea.col)] += ea.val * rhs;
        for (const auto &eb : row) {
            ata.emplace_back(ea.col, eb.col, ea.val * eb.val);
        }
    }
}


/// Merge the corners of the chart across its non-seam intra-chart edges.
/// A corner is a half-edge slot: half-edge `3f + e` stands for the corner at its `end` vertex;
/// slot `3 * i + e` is that corner for local face `i` (face `chart.faces[i]`).
static DisjointSet mergeCorners(const HalfEdgeMesh &mesh, const Chart &chart) {
    const std::size_t numFaces = chart.faces.size();

    // mesh face -> local face, to tell whether a neighbour belongs to this chart
    std::unordered_map<int, int> localFace;
    localFace.reserve(numFaces * 2);
    for (std::size_t i = 0; i < numFaces; ++i) {
        localFace.emplace(chart.faces[i], static_cast<int>(i));
    }

    // sorted seams, so the per-edge test is a binary search
    std::vector<int> seams = chart.seamHalfEdges;
    std::ranges::sort(seams);

    const auto slotOf = [&](HeIndex h, int face) {
        return 3 * face + static_cast<int>(h % 3);
    };

    // walk every half-edge once; an interior non-seam edge glues the two corners at each of its ends
    DisjointSet corners(numFaces * 3);
    for (std::size_t slot = 0; slot < numFaces * 3; ++slot) {
        const std::size_t i = slot / 3;
        const HeIndex     h = static_cast<HeIndex>(chart.faces[i]) * 3 + static_cast<HeIndex>(slot % 3);
        const HeIndex     s = mesh.sym(h);
        if (s == kNoHalfedge || s < h) {
            continue; // each edge once, from its lower half
        }
        const auto other = localFace.find(static_cast<int>(mesh.halfedge(s).face));
        if (other == localFace.end() || std::ranges::binary_search(seams, static_cast<int>(h))) {
            continue; // the neighbour is outside the chart, or the edge is cut
        }

        // h runs a->b, s runs b->a: tip(h) ~ origin(s), origin(h) ~ tip(s)
        const int fi = static_cast<int>(i);
        const int fj = other->second;
        corners.unite(slotOf(h, fi), slotOf(mesh.next(mesh.next(s)), fj));
        corners.unite(slotOf(mesh.next(mesh.next(h)), fi), slotOf(s, fj));
    }
    return corners;
}


/// Build the local vertex structure of the chart by corner.  Half-edge `3f + e` is corner
/// `2 - e` (a, b, c) of face `f` (`HalfEdgeMesh::writeFace`); each merge class of mergeCorners is one
/// local vertex.  Inter-chart and seam edges merge nothing, so they surface as the boundary of the
/// local mesh: the corner-wise equivalent of `CChartification::CutBoundaryPath` +
/// `CreateVertexStructure` (`legacy2003/chartification/cchartification.cpp:395-500`).
/// Local vertices are numbered in first-encounter order of the face scan (corner a, b, c).
static CornerMesh buildCornerMesh(const HalfEdgeMesh &mesh, const Chart &chart) {
    const std::size_t numFaces = chart.faces.size();

    // which corners are one vertex
    DisjointSet corners = mergeCorners(mesh, chart);

    CornerMesh out;
    out.localTri.resize(numFaces * 3);

    // one local label per merge class, keyed by the class root
    std::vector<int> labelOfRoot(numFaces * 3, -1);

    // scan the corners in face order; corner k of face i is slot 2 - k
    for (std::size_t corner = 0; corner < numFaces * 3; ++corner) {
        const std::size_t i = corner / 3;
        const std::size_t k = corner % 3;
        const int         slot = 3 * static_cast<int>(i) + 2 - static_cast<int>(k);
        int              &label = labelOfRoot[static_cast<std::size_t>(corners.find(slot))];

        // first corner of a class opens a new local vertex
        if (label < 0) {
            label = static_cast<int>(out.localToMesh.size());
            const HeIndex h = static_cast<HeIndex>(chart.faces[i]) * 3 + static_cast<HeIndex>(2 - k);
            out.localToMesh.push_back(static_cast<int>(mesh.halfedge(h).end));
        }
        out.localTri[corner] = label;
    }
    return out;
}


/// The two boundary vertices farthest apart in 3D space, as LOCAL indices.  The boundary is the cut
/// local mesh's own (`boundaryPointIndices()`): the implicit `sym == NULL` boundary left by
/// `CutBoundaryPath` (`cchartification.cpp:484-500`).  The O(n^2) farthest-pair search is `Diameter()`
/// (`cconformalmapper.cpp:48-82`); the pair is pinned as in `InitPinned()` (`:135-193`).
/// Throws std::runtime_error if a boundary index lies outside the local vertex range.
/// A chart without two boundary vertices (closed) falls back to local vertices 0 and 1.
static std::pair<int, int> findPinVertices(const HalfEdgeMesh     &localMesh,
                                           const std::vector<int> &localToMesh,
                                           const MeshData         &data) {
    const int numLocal = static_cast<int>(localToMesh.size());

    // the cut mesh's own boundary, guarded against indices outside the local range
    std::vector<int> bverts;
    for (const std::uint32_t li : localMesh.boundaryPointIndices()) {
        if (li >= localToMesh.size()) {
            throw std::runtime_error("findPinVertices: boundary vertex " + std::to_string(li) +
                                     " outside the chart's local vertices");
        }
        bverts.push_back(static_cast<int>(li));
    }

    // closed chart: no boundary to pick from, any two vertices will do
    if (bverts.size() < 2) {
        bverts.clear();
        bverts.push_back(0);
        if (numLocal >= 2) {
            bverts.push_back(1);
        }
    }

    const auto pos = [&](int local) {
        return vertexPos(data, localToMesh[static_cast<std::size_t>(local)]);
    };

    // O(n^2) scan for the farthest pair; strict '>' keeps the first pair on ties
    int   pin0 = bverts[0];
    int   pin1 = bverts.size() > 1 ? bverts[1] : bverts[0];
    float maxDist = 0.0f;
    for (std::size_t ii = 0; ii < bverts.size(); ++ii) {
        const math::Vec3f pi = pos(bverts[ii]);
        for (std::size_t jj = ii + 1; jj < bverts.size(); ++jj) {
            const float d = (pos(bverts[jj]) - pi).length();
            if (d > maxDist) {
                maxDist = d;
                pin0 = bverts[ii];
                pin1 = bverts[jj];
            }
        }
    }
    return { pin0, pin1 };
}


/// Assemble the LSCM normal equations AᵀAx = Aᵀb for one chart: per triangle, the complex
/// Cauchy-Riemann gradients (Gram-Schmidt in the tangent plane) are accumulated into `ata` and `atb`.
/// Pinned-vertex contributions are folded into `atb` immediately.
/// Ported from CConformalMapper::Parametrize (cconformalmapper.cpp:268–370).
static void assembleLSCMSystem(const Chart                         &chart,
                               const MeshData                      &data,
                               const std::vector<int>              &localTri,
                               const std::vector<int>              &freeIdx,
                               const PinState                      &pins,
                               std::vector<Eigen::Triplet<double>> &ata,
                               std::vector<double>                 &atb) {
    for (std::size_t i = 0; i < chart.faces.size(); ++i) {
        // mesh vertex ids and local ids of the triangle's corners
        const int base = chart.faces[i] * 3;
        const int va = data.indices[static_cast<std::size_t>(base + 0)];
        const int vb = data.indices[static_cast<std::size_t>(base + 1)];
        const int vc = data.indices[static_cast<std::size_t>(base + 2)];
        const int la = localTri[3 * i + 0];
        const int lb = localTri[3 * i + 1];
        const int lc = localTri[3 * i + 2];

        const math::Vec3f posA = vertexPos(data, va);
        const math::Vec3f posB = vertexPos(data, vb);
        const math::Vec3f posC = vertexPos(data, vc);

        // Gram-Schmidt orthonormalization (ports Gradient() + orthonormal()).
        const math::Vec3f ba = posB - posA;
        const math::Vec3f ca = posC - posA;
        const math::Vec3f e1 = ba.normalized();
        math::Vec3f       e2 = ca - e1 * e1.dot(ca);
        e2.normalize();

        // Triangle 2D coordinates in the orthonormal tangent basis (corner a sits at the origin).
        constexpr float x1 = 0.0f;
        constexpr float y1 = 0.0f;
        const float     x2 = ba.dot(e1);
        const float     y2 = ba.dot(e2);
        const float     x3 = ca.dot(e1);
        const float     y3 = ca.dot(e2);

        // Degenerate triangle: no usable gradient, contributes no row.
        const float dT = (x1 * y2 - y1 * x2) + (x2 * y3 - y2 * x3) + (x3 * y1 - y3 * x1);
        if (std::abs(dT) < 1.0e-12f) {
            continue;
        }
        const float sqrtDT = std::sqrt(std::abs(dT));

        // Complex Cauchy-Riemann gradients w = (wr, wi) / sqrtDT.
        const int    lids[3] = { la, lb, lc };
        const double wrs[3] = { static_cast<double>((x3 - x2) / sqrtDT),
                                static_cast<double>((x1 - x3) / sqrtDT),
                                static_cast<double>((x2 - x1) / sqrtDT) };
        const double wis[3] = { static_cast<double>((y3 - y2) / sqrtDT),
                                static_cast<double>((y1 - y3) / sqrtDT),
                                static_cast<double>((y2 - y1) / sqrtDT) };

        // One real and one imaginary equation per triangle, as sparse rows over the free unknowns.
        std::vector<RowEntry> rowReal;
        std::vector<RowEntry> rowImag;
        rowReal.reserve(6);
        rowImag.reserve(6);
        double rhsReal = 0.0;
        double rhsImag = 0.0;

        // A pinned corner has a known UV, so it moves to the right-hand side; a free corner becomes row entries.
        for (int k = 0; k < 3; ++k) {
            const int    lid = lids[k];
            const double wr = wrs[k];
            const double wi = wis[k];

            if (lid == pins.pin0Local) {
                rhsReal -= wr * pins.uPin0 + (-wi) * pins.vPin0;
                rhsImag -= wi * pins.uPin0 + wr * pins.vPin0;
            } else if (lid == pins.pin1Local) {
                rhsReal -= wr * pins.uPin1 + (-wi) * pins.vPin1;
                rhsImag -= wi * pins.uPin1 + wr * pins.vPin1;
            } else {
                const int fi = freeIdx[static_cast<std::size_t>(lid)];
                rowReal.push_back({ fi, wr });
                rowReal.push_back({ fi + pins.numFree, -wi });
                rowImag.push_back({ fi, wi });
                rowImag.push_back({ fi + pins.numFree, wr });
            }
        }

        accumRow(ata, atb, rowReal, rhsReal);
        accumRow(ata, atb, rowImag, rhsImag);
    }
}


/// Write the solution `x` back into `uv` for every local vertex: free vertices get
/// u = x[fi], v = x[fi + numFree]; pinned vertices get their fixed UV from `pins`.
/// Ports CConformalMapper::CopyResult (cconformalmapper.cpp:35-42).
static void backSubstituteUVs(const std::vector<double> &x,
                              const std::vector<int>    &freeIdx,
                              const PinState            &pins,
                              int                        numLocal,
                              std::vector<math::Vec2f>  &uv) {
    for (int li = 0; li < numLocal; ++li) {
        // pinned vertices were never unknowns: they take their fixed UV, the rest read the solution
        float u = 0.0f;
        float v = 0.0f;
        if (li == pins.pin0Local) {
            u = static_cast<float>(pins.uPin0);
            v = static_cast<float>(pins.vPin0);
        } else if (li == pins.pin1Local) {
            u = static_cast<float>(pins.uPin1);
            v = static_cast<float>(pins.vPin1);
        } else {
            const int fi = freeIdx[static_cast<std::size_t>(li)];
            u = static_cast<float>(x[static_cast<std::size_t>(fi)]);
            v = static_cast<float>(x[static_cast<std::size_t>(fi + pins.numFree)]);
        }

        uv[static_cast<std::size_t>(li)] = math::Vec2f { u, v };
    }
}


/// Signed doubled area of triangle abc in UV space, accumulated in double; positive when abc is counter-clockwise.
[[nodiscard]] static double orient(const UV &a, const UV &b, const UV &c) noexcept {
    return static_cast<double>(b.x - a.x) * (c.y - a.y) - static_cast<double>(b.y - a.y) * (c.x - a.x);
}


/// Flattening quality in O(n + B^2) (n triangles, B boundary edges; the B^2 is the boundary
/// self-crossing loop): flipped/degenerate triangles against the chart's orientation (sign of the summed
/// signed area), and boundary edges crossing each other.  Disk + no flip + simple boundary => the map
/// is injective.  Shared with the bijective mapper via chart_quality.h.
void measureQuality(HalfEdgeChart &ac) {
    const HalfEdgeMesh &m = ac.mesh;
    const int           nf = static_cast<int>(m.faceCount());
    ac.flippedTriangles = 0;
    ac.boundarySimple = true;

    // The summed signed area fixes the chart's orientation, whichever way the flattening came out.
    double total = 0.0;
    for (int f = 0; f < nf; ++f) {
        total += orient(ac.uv[localCorner(ac, f, 0)], ac.uv[localCorner(ac, f, 1)], ac.uv[localCorner(ac, f, 2)]);
    }

    // A triangle against that orientation, or with no area (NaN included), is flipped.
    const double sign = total < 0.0 ? -1.0 : 1.0;
    for (int f = 0; f < nf; ++f) {
        const double a = sign * math::orient2d(ac.uv[localCorner(ac, f, 0)], ac.uv[localCorner(ac, f, 1)],
                                               ac.uv[localCorner(ac, f, 2)]);
        if (!(a > 0.0)) {
            ++ac.flippedTriangles;
        }
    }

    // Boundary edges (no sym) as (origin, end) pairs.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
    for (HeIndex h = 0; h < m.halfEdgeCount(); ++h) {
        if (m.sym(h) == kNoHalfedge) {
            edges.emplace_back(m.halfedge(m.next(m.next(h))).end, m.halfedge(h).end);
        }
    }

    // Every edge pair once; the first proper crossing makes the boundary non-simple.
    for (std::size_t i = 0; i < edges.size() && ac.boundarySimple; ++i) {
        const auto [a, b] = edges[i];
        for (std::size_t j = i + 1; j < edges.size(); ++j) {
            const auto [c, d] = edges[j];
            if (a == c || a == d || b == c || b == d) {
                continue; // neighbours share an endpoint and touch there by construction
            }
            if (math::segmentsCross(ac.uv[a], ac.uv[b], ac.uv[c], ac.uv[d])) {
                ac.boundarySimple = false;
                break;
            }
        }
    }
}


/// Ridge added to K's diagonal when the plain factorization fails, relative to K's largest diagonal
/// entry (the retry of `src/seam/erase.cpp`'s Solver).  K is singular to working precision when a
/// vertex's incident triangles were all skipped as degenerate (`|dT| < 1.0e-12f`, zero rows); the
/// 2003 CG tolerated that, a direct solve does not.
constexpr double kRidge = 1e-9;

/// Solve K x = atb by sparse LDLT; on failure retry once with `K + ridge * I`.  The ridge leaves a
/// vertex without any equation at its right-hand side (0), so the solution stays finite.
/// Returns false when both factorizations fail; x is then left as given (zeros).
static bool solveNormalEquations(const Eigen::SparseMatrix<double> &k, const std::vector<double> &atb,
                                 std::vector<double> &x) {
    const auto                                         n = static_cast<Eigen::Index>(atb.size());
    const Eigen::Map<const Eigen::VectorXd>            rhs(atb.data(), n);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(k);

    // Plain factorization failed: regularize the diagonal relative to its largest entry and retry once.
    if (ldlt.info() != Eigen::Success) {
        double maxDiag = 0.0;
        for (Eigen::Index i = 0; i < n; ++i) {
            maxDiag = std::max(maxDiag, k.coeff(i, i));
        }
        Eigen::SparseMatrix<double> ridge(n, n);
        ridge.setIdentity();
        ldlt.compute(k + kRidge * (maxDiag > 0.0 ? maxDiag : 1.0) * ridge);
    }

    if (ldlt.info() != Eigen::Success) {
        return false;
    }

    Eigen::Map<Eigen::VectorXd>(x.data(), n) = ldlt.solve(rhs);
    return true;
}


/// LSCM-flatten a single chart into 2D UV space.  Port of CConformalMapper::Parametrize(CChart*)
/// (cconformalmapper.cpp:247-322) with its helpers Gradient() (:89-126), Diameter() (:48-82) and
/// InitPinned() (:135-193).  Instead of the 2003 two-matrix system (am free, bm pinned) the normal
/// equations AᵀA are accumulated as sparse triplets and factorized by SimplicialLDLT.
HalfEdgeChart flatten(const HalfEdgeMesh &mesh,
                      const MeshData     &data,
                      const Chart        &chart) {
    HalfEdgeChart result;
    result.faces = chart.faces;

    // nothing to flatten
    if (chart.faces.empty()) {
        return result;
    }

    // Local vertices by corner (cut along chart.seamHalfEdges) and the owned local half-edge mesh
    // built from them; boundary and adjacency are queries on it.
    // reinterpret_cast is safe: Face is {int a, b, c} with the layout of three consecutive ints
    // (mesh_types.h; used identically in tests/half_edge_boundary_test.cpp::makeHem).
    CornerMesh cm = buildCornerMesh(mesh, chart);
    const int  numLocal = static_cast<int>(cm.localToMesh.size());
    result.localToMesh = std::move(cm.localToMesh);
    result.uv.resize(static_cast<std::size_t>(numLocal));
    result.mesh = HalfEdgeMesh::create(
        numLocal, static_cast<int>(chart.faces.size()),
        reinterpret_cast<const decimate::Face *>(cm.localTri.data()));

    // Pins come from the cut mesh boundary (Spec 09-01 amendment).
    const auto [pin0Local, pin1Local] = findPinVertices(result.mesh, result.localToMesh, data);

    // Pinned vertices are not unknowns: remap the rest densely (freeIdx[li] = free index, or -1).
    std::vector<int> freeIdx(static_cast<std::size_t>(numLocal), -1);
    int              numFree = 0;
    for (int li = 0; li < numLocal; ++li) {
        if (li != pin0Local && li != pin1Local) {
            freeIdx[static_cast<std::size_t>(li)] = numFree++;
        }
    }

    // Assemble the LSCM system.
    // Pinned UV values (legacy cconformalmapper.cpp:160-163): pin0 -> UV (0, 0), pin1 -> UV (1, 1).
    const PinState                      pins { pin0Local, pin1Local, numFree };
    const int                           sysSize = 2 * numFree;
    std::vector<Eigen::Triplet<double>> ata;
    std::vector<double>                 atb(static_cast<std::size_t>(sysSize), 0.0);
    assembleLSCMSystem(chart, data, cm.localTri, freeIdx, pins, ata, atb);

    // Solve by sparse Cholesky.
    Eigen::SparseMatrix<double> k(sysSize, sysSize);
    k.setFromTriplets(ata.begin(), ata.end()); // sums duplicates, as addEntry accumulated
    std::vector<double> x(static_cast<std::size_t>(sysSize), 0.0);
    const bool          solved = solveNormalEquations(k, atb, x);

    // Back-substitute the UVs and grade the result.
    backSubstituteUVs(x, freeIdx, pins, numLocal, result.uv);
    measureQuality(result);

    if (!solved) {
        // Defined result: free vertices at (0,0), and the chart is never reported injective.
        result.flippedTriangles = std::max(result.flippedTriangles, 1);
    }
    return result;
}


} // namespace decimate::chart
