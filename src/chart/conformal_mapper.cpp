// This file is part of the decimate library; see decimate.h for version/license details.
/// @file conformal_mapper.cpp
/// @brief LSCM chart flattening — behavior-preserving port of CConformalMapper.
///
/// ## Algorithm (Lévy et al. 2002 — Least-Squares Conformal Maps)
///
/// For each chart:
///   1. Collect the unique local vertex set over all chart faces (face-scan order).
///   2. Identify boundary vertices; pin the two farthest-apart pair to fixed UVs.
///   3. For each triangle, compute complex Cauchy-Riemann gradients via
///      Gram-Schmidt orthonormalization in the triangle's 3D tangent plane.
///   4. Directly assemble the normal equations AᵀA x = Aᵀb for free
///      (unpinned) vertices and solve via ConjugateGradientSolver.
///   6. Back-substitute to fill HalfEdgeChart::uv for every local vertex.
///
/// ## Operator mapping from legacy Vec3 to decimate::math::Vec3
///
///   Legacy `(ba) * e1` — Vec3::operator*(Vec3) was the dot product.
///   → `ba.dot(e1)` — decimate::math::Vec3::operator*(Vec3) is component-wise; must use .dot().
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
///   Matches legacy CopyResult (cconformalmapper.cpp:36-42):
///     vertices[offset + i].u = (float) vec[i];
///     vertices[offset + i].v = (float) vec[num + i];
///
/// Boundary contract: no GL/OS/I/O/global state; C++20; RAII throughout.

#include <decimate/chart/parametrize.h>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <decimate/math/vec2.h>
#include <decimate/math/vec3.h>
#include <decimate/math/cg_solver.h>
#include <decimate/chart/chart.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>

namespace decimate::chart {

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

namespace {

/// @brief Extract the 3D position of a mesh vertex.
[[nodiscard]] decimate::math::Vec3<float> vertexPos(const MeshData& data, int vid) noexcept {
    return data.positions[static_cast<std::size_t>(vid)];
}

/// @brief One entry in a sparse matrix row: (column, value).
struct RowEntry {
    int    col;
    double val;
};

/// @brief Accumulate one sparse row's outer product into AtA and Atb.
///
/// For each pair of entries (ea, eb) in row: AtA[ea.col][eb.col] += ea.val * eb.val.
/// For each entry ea: atb[ea.col] += ea.val * rhs.
void accumRow(decimate::math::ConjugateGradientSolver& cg,
              std::vector<double>&                  atb,
              const std::vector<RowEntry>&          row,
              double                                rhs) {
    for (const auto& ea : row) {
        atb[static_cast<std::size_t>(ea.col)] += ea.val * rhs;
        for (const auto& eb : row) {
            cg.addEntry(ea.col, eb.col, ea.val * eb.val);
        }
    }
}

/// @brief Local vertex set for a chart.
struct LocalVertexSet {
    std::vector<int>             meshId;   ///< local index → mesh vertex id.
    std::unordered_map<int, int> localOf;  ///< mesh vertex id → local index.
};

/// @brief Collect the unique local vertex set from chart faces.
///
/// Scans MeshData::indices in face-list order; first encounter determines
/// each vertex's local index.  Ports the Step-1 block of Parametrize().
LocalVertexSet collectLocalVertices(const Chart& chart, const MeshData& data) {
    LocalVertexSet vs;
    vs.meshId.reserve(chart.faces.size() * 3);
    for (int fid : chart.faces) {
        const int base = fid * 3;
        for (int k = 0; k < 3; ++k) {
            const int vid = data.indices[static_cast<std::size_t>(base + k)];
            if (vs.localOf.find(vid) == vs.localOf.end()) {
                vs.localOf[vid] = static_cast<int>(vs.meshId.size());
                vs.meshId.push_back(vid);
            }
        }
    }
    return vs;
}

/// @brief Find the two boundary vertices farthest apart in 3D space.
///
/// Derives the boundary vertex set from the cut chart's own
/// `localMesh.boundaryPointIndices()` — the faithful transform of
/// `CChartification::CreateVertexStructure` (`cchartification.cpp:395`), which
/// cuts the mesh (`CutBoundaryPath` nulls `sym` on boundary half-edges,
/// `cchartification.cpp:484+`) making the boundary implicit (`sym==NULL`).
/// `boundaryPointIndices()` walks exactly that implicit boundary, matching the
/// 2003 `Boundaries()→packer` hand-off exactly.
///
/// Each local boundary index is mapped to its mesh vertex id via `localToMesh`,
/// then the O(n²) farthest-pair search (Diameter) finds the pin pair.
///
/// The 08-03 checked local-lookup (D1 guard) is preserved via the
/// `boundaryPointIndices()` local indices which are guaranteed to be in
/// `[0, localToMesh.size())` by the cut mesh invariant.
///
/// Ports Diameter() + InitPinned() (cconformalmapper.cpp:48-193).
std::pair<int, int> findPinVertices(const HalfEdgeMesh&      localMesh,
                                    const std::vector<int>&  localToMesh,
                                    const MeshData&          data,
                                    const LocalVertexSet&    vs) {
    // Collect boundary vertex mesh-ids from the cut mesh boundary loop.
    // boundaryPointIndices() walks sym==kNoHalfedge half-edges — exactly the
    // boundary left implicit by the cut (cchartification.cpp:484+).
    const std::vector<std::uint32_t> bndLocal = localMesh.boundaryPointIndices();

    std::vector<int> bverts;
    bverts.reserve(bndLocal.size());
    for (std::uint32_t li : bndLocal) {
        bverts.push_back(localToMesh[static_cast<std::size_t>(li)]);
    }

    const int numLocal = static_cast<int>(vs.meshId.size());
    if (bverts.size() < 2) {
        bverts.clear();
        if (numLocal >= 1) { bverts.push_back(vs.meshId[0]); }
        if (numLocal >= 2) { bverts.push_back(vs.meshId[1]); }
    }

    // O(|bverts|²) farthest-pair search (Diameter, cconformalmapper.cpp:48-82).
    int   pin0MeshId = bverts[0];
    int   pin1MeshId = bverts.size() > 1 ? bverts[1] : bverts[0];
    float maxDist    = 0.0f;
    for (std::size_t ii = 0; ii < bverts.size(); ++ii) {
        const decimate::math::Vec3<float> pi = vertexPos(data, bverts[ii]);
        for (std::size_t jj = ii + 1; jj < bverts.size(); ++jj) {
            const float d = (vertexPos(data, bverts[jj]) - pi).length();
            if (d > maxDist) {
                maxDist    = d;
                pin0MeshId = bverts[ii];
                pin1MeshId = bverts[jj];
            }
        }
    }
    // Checked lookup: if a pin vertex is absent from localOf, fail loudly with
    // context rather than throwing a bare std::out_of_range.  (08-03 Fix 2 / D1 guard)
    auto checkedLocal = [&](int meshVid) -> int {
        const auto it = vs.localOf.find(meshVid);
        if (it == vs.localOf.end()) {
            throw std::runtime_error(
                "findPinVertices: vertex " + std::to_string(meshVid) +
                " absent from localOf — boundary vertex not in chart.faces?");
        }
        return it->second;
    };

    const int pin0Local = checkedLocal(pin0MeshId);
    const int pin1Local = checkedLocal(pin1MeshId);
    return {pin0Local, pin1Local};
}

/// @brief Pinning state: local indices + UV values for the two pinned vertices.
struct PinState {
    int    pin0Local;   ///< Local index of pin0 (→ UV 0,0).
    int    pin1Local;   ///< Local index of pin1 (→ UV 1,1).
    int    numFree;     ///< Number of free (unpinned) local vertices.
    double uPin0 = 0.0; ///< u for pin0.
    double vPin0 = 0.0; ///< v for pin0.
    double uPin1 = 1.0; ///< u for pin1.
    double vPin1 = 1.0; ///< v for pin1.
};

/// @brief Assemble the LSCM normal equations AᵀAx = Aᵀb for one chart.
///
/// Iterates over the chart's triangles, computes complex Cauchy-Riemann
/// gradients via Gram-Schmidt orthonormalization, and accumulates them
/// directly into `cg` (AtA matrix) and `atb` (Aᵀb vector).
/// Pinned-vertex contributions are folded into `atb` immediately.
///
/// Ported from CConformalMapper::Parametrize (cconformalmapper.cpp:268–370).
void assembleLSCMSystem(const Chart&                          chart,
                               const MeshData&                       data,
                               const std::unordered_map<int, int>&   localOf,
                               const std::vector<int>&               freeIdx,
                               const PinState&                       pins,
                               decimate::math::ConjugateGradientSolver&  cg,
                               std::vector<double>&                  atb) {
    for (int fid : chart.faces) {
        const int base = fid * 3;
        const int va   = data.indices[static_cast<std::size_t>(base + 0)];
        const int vb   = data.indices[static_cast<std::size_t>(base + 1)];
        const int vc   = data.indices[static_cast<std::size_t>(base + 2)];

        // Checked local lookup — throw with context if any face vertex is absent
        // from localOf.  After 08-03 Fix 1 this cannot fire (collectLocalVertices
        // and assembleLSCMSystem both iterate the same chart.faces), but a bare
        // .at() would reproduce the D1 failure mode silently — prefer loud failure.
        auto localAt = [&](int vid) -> int {
            const auto it = localOf.find(vid);
            if (it == localOf.end()) {
                throw std::runtime_error(
                    "assembleLSCMSystem: vertex " + std::to_string(vid) +
                    " absent from localOf in chart " + std::to_string(chart.id) +
                    " — invariant violation (Fix 1 not applied?)");
            }
            return it->second;
        };
        const int la = localAt(va);
        const int lb = localAt(vb);
        const int lc = localAt(vc);

        const decimate::math::Vec3<float> posA = vertexPos(data, va);
        const decimate::math::Vec3<float> posB = vertexPos(data, vb);
        const decimate::math::Vec3<float> posC = vertexPos(data, vc);

        // Gram-Schmidt orthonormalization (ports Gradient() + orthonormal()).
        const decimate::math::Vec3<float> ba = posB - posA;
        const decimate::math::Vec3<float> ca = posC - posA;
        decimate::math::Vec3<float> e1 = ba.normalized();
        decimate::math::Vec3<float> e2 = ca - e1 * e1.dot(ca);
        e2.normalize();

        // Triangle 2D coordinates in the orthonormal tangent basis.
        constexpr float x1 = 0.0f;
        constexpr float y1 = 0.0f;
        const float     x2 = ba.dot(e1);
        const float     y2 = ba.dot(e2);
        const float     x3 = ca.dot(e1);
        const float     y3 = ca.dot(e2);

        const float dT = (x1 * y2 - y1 * x2) + (x2 * y3 - y2 * x3)
                       + (x3 * y1 - y3 * x1);
        if (std::abs(dT) < 1e-12f) { continue; }  // degenerate — skip
        const float sqrtDT = std::sqrt(std::abs(dT));

        // Complex Cauchy-Riemann gradients w = (wr, wi) / sqrtDT.
        const int    lids[3] = {la, lb, lc};
        const double wrs[3]  = {static_cast<double>((x3 - x2) / sqrtDT),
                                 static_cast<double>((x1 - x3) / sqrtDT),
                                 static_cast<double>((x2 - x1) / sqrtDT)};
        const double wis[3]  = {static_cast<double>((y3 - y2) / sqrtDT),
                                 static_cast<double>((y1 - y3) / sqrtDT),
                                 static_cast<double>((y2 - y1) / sqrtDT)};

        std::vector<RowEntry> rowReal;
        std::vector<RowEntry> rowImag;
        rowReal.reserve(6);
        rowImag.reserve(6);
        double rhsReal = 0.0;
        double rhsImag = 0.0;

        for (int k = 0; k < 3; ++k) {
            const int lid = lids[k];
            const double wr = wrs[k];
            const double wi = wis[k];

            if (lid == pins.pin0Local) {
                rhsReal -= wr * pins.uPin0 + (-wi) * pins.vPin0;
                rhsImag -= wi * pins.uPin0 +   wr  * pins.vPin0;
            } else if (lid == pins.pin1Local) {
                rhsReal -= wr * pins.uPin1 + (-wi) * pins.vPin1;
                rhsImag -= wi * pins.uPin1 +   wr  * pins.vPin1;
            } else {
                const int fi = freeIdx[static_cast<std::size_t>(lid)];
                rowReal.push_back({fi,              wr});
                rowReal.push_back({fi + pins.numFree, -wi});
                rowImag.push_back({fi,              wi});
                rowImag.push_back({fi + pins.numFree,  wr});
            }
        }
        accumRow(cg, atb, rowReal, rhsReal);
        accumRow(cg, atb, rowImag, rhsImag);
    }
}

/// @brief Write CG solution `x` back into `uv` for every local vertex.
///
/// Free vertices: u = x[fi], v = x[fi + numFree].
/// Pinned vertices: hardcoded UV from `pins`.
///
/// Ports CConformalMapper::CopyResult (cconformalmapper.cpp:36-42).
void backSubstituteUVs(const std::vector<double>&    x,
                               const std::vector<int>&       freeIdx,
                               const PinState&               pins,
                               int                           numLocal,
                               std::vector<decimate::math::Vec2<float>>& uv) {
    for (int li = 0; li < numLocal; ++li) {
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
        uv[static_cast<std::size_t>(li)] = decimate::math::Vec2<float>{u, v};
    }
}

}  // anonymous namespace

// ---------------------------------------------------------------------------
// flatten — public entry point
// ---------------------------------------------------------------------------

/// @brief LSCM-flatten a single chart into 2D UV space.
///
/// Behavior-preserving port of CConformalMapper::Parametrize(CChart*)
/// (cconformalmapper.cpp:247–322) combined with its helper methods
/// Gradient() (lines 89–126), Diameter() (lines 48–82), and
/// InitPinned() (lines 135–193).
///
/// The normal-equation assembly replaces the legacy SMCR two-matrix approach
/// (am for free, bm for pinned) with direct accumulation into the
/// ConjugateGradientSolver, which is the correct solver for this system
/// class (NOT the Spec-05 SMat / Gauss-Seidel path).
HalfEdgeChart flatten([[maybe_unused]] const HalfEdgeMesh& mesh,
                      const MeshData&                       data,
                      const Chart&                          chart) {
    HalfEdgeChart result;
    result.faces = chart.faces;

    if (chart.faces.empty()) {
        return result;
    }

    // Step 1: Collect the unique local vertex set over all chart faces.
    const LocalVertexSet vs       = collectLocalVertices(chart, data);
    const int            numLocal = static_cast<int>(vs.meshId.size());
    result.localToMesh = vs.meshId;
    result.uv.resize(static_cast<std::size_t>(numLocal));

    // Step 2: Build the per-face local-vertex index array (3 ints per face).
    // This is the face array fed to HalfEdgeMesh::create — the local triangulation
    // that makes boundary and adjacency derivable as queries.
    // `Face` is 3 contiguous ints (mesh_types.h:18), same layout as int[3].
    std::vector<int> localTri;
    localTri.reserve(chart.faces.size() * 3);
    for (int fid : chart.faces) {
        const int base = fid * 3;
        for (int k = 0; k < 3; ++k) {
            const int vid = data.indices[static_cast<std::size_t>(base + k)];
            localTri.push_back(vs.localOf.at(vid));
        }
    }

    // Step 3: Build the owned local HalfEdgeMesh from the local triangulation.
    // reinterpret_cast is safe: Face is {int a, b, c} with the same layout as
    // three consecutive ints (confirmed in mesh_types.h and used identically in
    // tests/half_edge_boundary_test.cpp::makeHem).
    const int numLocalFaces = static_cast<int>(chart.faces.size());
    result.mesh = HalfEdgeMesh::create(
        numLocal, numLocalFaces,
        reinterpret_cast<const decimate::Face*>(localTri.data()));

    // Step 4: Find the two pin vertices from the cut mesh boundary.
    // Uses result.mesh.boundaryPointIndices() — the faithful port of
    // CChartification::CreateVertexStructure (cchartification.cpp:395), where
    // CutBoundaryPath nulls sym on boundary half-edges (cchartification.cpp:484+)
    // making the boundary implicit (sym==NULL).  boundaryPointIndices() walks
    // exactly that implicit boundary, matching the 2003 Boundaries()→packer
    // hand-off.  This is the correct boundary for pin selection (Spec 09-01 AMENDMENT).
    const auto [pin0Local, pin1Local] =
        findPinVertices(result.mesh, result.localToMesh, data, vs);

    // Step 5: Build free-vertex index remapping (freeIdx[li] = free index, or -1).
    std::vector<int> freeIdx(static_cast<std::size_t>(numLocal), -1);
    int numFree = 0;
    for (int li = 0; li < numLocal; ++li) {
        if (li != pin0Local && li != pin1Local) {
            freeIdx[static_cast<std::size_t>(li)] = numFree++;
        }
    }

    // Steps 6-7: Assemble LSCM system, solve via CG, back-substitute UVs.
    // Pinned UV values (legacy cconformalmapper.cpp:160-163):
    //   pin0 → UV (0, 0),  pin1 → UV (1, 1)
    const PinState pins{pin0Local, pin1Local, numFree};
    const int      sysSize = 2 * numFree;
    decimate::math::ConjugateGradientSolver cg(sysSize);
    std::vector<double> atb(static_cast<std::size_t>(sysSize), 0.0);
    assembleLSCMSystem(chart, data, vs.localOf, freeIdx, pins, cg, atb);
    for (int col = 0; col < sysSize; ++col) {
        cg.setRhs(col, atb[static_cast<std::size_t>(col)]);
    }
    const std::vector<double> x = cg.solve(1e-12, 15000);
    backSubstituteUVs(x, freeIdx, pins, numLocal, result.uv);
    return result;
}

}  // namespace decimate::chart
