// This file is part of the decimate library; see decimate.h for version/license details.
/// @file chart_legacy_golden.cpp
/// @brief Spec 19-01: live A/B of segment() against the verbatim-compiled 2003 CChartification.
///
/// For each mesh x cap the 2003 cutter (legacy_side.cpp) and `segment(..., Legacy2003Middle)` run
/// on the same input.  2003 charts carry no ids (CChartSet is a pointer-ordered std::set<CChart*>,
/// legacy2003/chartification/cchartset.h:24), so charts are matched by MAXIMAL OVERLAP: each
/// legacy chart is paired with the port chart sharing most faces; a face "differs" when its port
/// chart is not the match of its legacy chart.  Exit status != 0 on any mismatch.
///
/// spot.obj is NOT in the list: on a closed mesh the 2003 ApplySeeds reads the uninitialised
/// CChart::searchIndex (legacy2003/chartification/cchart.h:42, ctor legacy2003/chartification/cchart.cpp:23-31), so the
/// legacy output there is undefined behaviour that the port deliberately does not emulate (D-d).
///
/// Spec 19-02 extends the A/B to the chart-local vertex structure: the summed `CChart::VertexCount()`
/// of 2003 (`CreateVertexStructure`, legacy2003/chartification/cchartification.cpp:395-480) against the summed
/// `HalfEdgeChart::localToMesh.size()` of `flatten()`.  They are equal except where the port erases
/// the stale Straighten edge (a documented divergence, segmenter.cpp): `localDelta` pins the exact
/// port-minus-2003 difference per case, so any other drift fails.  Measured with the stale fix
/// disabled, every case equals 2003 exactly (delta 0), which attributes all nonzero deltas
/// (cow/24 -1, cow/32 -1, cow/64 -1, horse/128 -4, bunny/8 -21, bunny/32 -18) to it.
///
/// Only this test uses CentroidMode::Legacy2003Middle (2003 middle(), legacy2003/math/algorithm.h:70-73).
#include <array>
#include <chrono>
#include <cstdio>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <filesystem>

#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>

#ifndef MESH_MODEL_DIR
#error "MESH_MODEL_DIR must be defined by CMake"
#endif

/// The 2003 side, defined in legacy_side.cpp (compiled verbatim as C++14; see that file).
int legacyChartify(int np, const float *pts, int nf, const int *idx, int k,
                   int *faceToChart, int *localVerts);

namespace {

/// One mesh x cap case with its pinned expectations.
struct Case {
    const char *model;
    int         cap;
    int         expectedCharts; ///< Same on both sides (compile-time measurement, Spec 19-01).
    int         localDelta;     ///< Summed local vertices, port minus 2003 (Spec 19-02; stale-edge fix).
};

/// The measured cases; chart counts are equal on both sides (0 = not pinned), localDelta is the pinned port-minus-2003.
constexpr std::array<Case, 14> kCases = { {
    { "cow.obj", 16, 16, 0 },
    { "cow.obj", 24, 0, -1 },
    { "cow.obj", 32, 31, -1 },
    { "cow.obj", 64, 61, -1 },
    { "cow.obj", 128, 113, 0 },
    { "horse.iv", 16, 16, 0 },
    { "horse.iv", 32, 30, 0 },
    { "horse.iv", 64, 55, 0 },
    { "horse.iv", 128, 100, -4 },
    { "bunny.iv", 8, 0, -21 },
    { "bunny.iv", 16, 16, 0 },
    { "bunny.iv", 32, 32, -18 },
    { "bunny.iv", 64, 60, 0 },
    { "bunny.iv", 128, 114, 0 },
} };


/// Faces whose port chart is not the max-overlap match of their legacy chart.
int differingFaces(const std::vector<int> &legacy, const std::vector<int> &port) {
    // Overlap histogram, then the best port chart per legacy chart.
    std::map<std::pair<int, int>, int> overlap;
    for (std::size_t f = 0; f < legacy.size(); ++f) {
        ++overlap[{ legacy[f], port[f] }];
    }
    std::map<int, std::pair<int, int>> best; // legacy chart -> (port chart, overlap)
    for (const auto &[key, n] : overlap) {
        auto &b = best[key.first];
        if (n > b.second) {
            b = { key.second, n };
        }
    }

    int diff = 0;
    for (std::size_t f = 0; f < legacy.size(); ++f) {
        if (best[legacy[f]].first != port[f]) {
            ++diff;
        }
    }
    return diff;
}


/// Summed chart-local vertex count of the port (`flatten()`) and total exported seam half-edges.
std::pair<int, int> portLocalVertices(const decimate::HalfEdgeMesh &hem, const decimate::MeshData &d,
                                      const decimate::chart::ChartSet &cs) {
    int local = 0;
    int seams = 0;
    for (const decimate::chart::Chart &ch : cs.charts) {
        local += static_cast<int>(decimate::chart::flatten(hem, d, ch).localToMesh.size());
        seams += static_cast<int>(ch.seamHalfEdges.size());
    }
    return { local, seams };
}


/// Outcome of one mesh x cap case.
enum class Result { Pass,
                    Fail,
                    Skipped };


/// Runs both cutters on @p c and prints one row.
Result runCase(const Case &c, std::map<std::string, decimate::MeshData> &meshes) {
    using namespace decimate;
    const std::string path = std::string(MESH_MODEL_DIR) + "/" + c.model;

    // A missing model is a skip, not a failure (the models are optional fixtures).
    if (!std::filesystem::exists(path)) {
        std::printf("SKIP %s k=%d (model not present)\n", c.model, c.cap);
        return Result::Skipped;
    }

    // Meshes are loaded once and shared between the caps of one model.
    if (meshes.find(c.model) == meshes.end()) {
        meshes[c.model] = load(path);
    }

    // The legacy API takes flat arrays.
    const MeshData    &d = meshes[c.model];
    const std::size_t  nf = d.indices.size() / 3;
    const std::size_t  np = d.positions.size();
    std::vector<float> pts(3 * np);
    for (std::size_t i = 0; i < np; ++i) {
        pts[3 * i] = d.positions[i].x;
        pts[3 * i + 1] = d.positions[i].y;
        pts[3 * i + 2] = d.positions[i].z;
    }
    const std::vector<int> idx(d.indices.begin(), d.indices.end());
    std::vector<int>       legacyF(nf);
    std::vector<int>       localVerts(static_cast<std::size_t>(c.cap) + 1, 0);


    // Run both cutters on the same input and time them.
    const auto            t0 = std::chrono::steady_clock::now();
    const int             legacyCharts = legacyChartify(static_cast<int>(np), pts.data(), static_cast<int>(nf),
                                                        idx.data(), c.cap, legacyF.data(), localVerts.data());
    const auto            t1 = std::chrono::steady_clock::now();
    const HalfEdgeMesh    hem = HalfEdgeMesh::create(static_cast<int>(np), static_cast<int>(nf),
                                                     reinterpret_cast<const Face *>(d.indices.data()));
    const chart::ChartSet cs = chart::segment(hem, d, c.cap, chart::kStraightenThreshold,
                                              chart::CentroidMode::Legacy2003Middle);
    const auto            t2 = std::chrono::steady_clock::now();


    // A negative count is the legacy cutter's error return; using it as an offset below
    // (`begin() + legacyCharts`) would be undefined behaviour, so it fails the case right here.
    if (legacyCharts < 0) {
        std::printf("FAIL %s k=%d legacyChartify returned %d (hard error)\n", c.model, c.cap, legacyCharts);
        return Result::Fail;
    }

    // A face the 2003 cutter never reached (-1) makes the overlap match meaningless, so it fails the case.
    int dropped = 0;
    for (const int x : legacyF) {
        dropped += (x < 0) ? 1 : 0;
    }
    const int  portCharts = static_cast<int>(cs.charts.size());
    const int  diff = (dropped == 0) ? differingFaces(legacyF, cs.faceToChart) : -1;
    const bool chartsOk = legacyCharts == portCharts &&
                          (c.expectedCharts == 0 || legacyCharts == c.expectedCharts);

    const int legacyLocal =
        std::accumulate(localVerts.begin(), localVerts.begin() + legacyCharts, 0);
    const auto [portLocal, seams] = portLocalVertices(hem, d, cs);

    const bool ok = chartsOk && dropped == 0 && diff == 0 && portLocal - legacyLocal == c.localDelta;
    std::printf("%s %s k=%d faces=%zu legacyCharts=%d portCharts=%d expected=%d legacyDropped=%d "
                "facesDiffering=%d localVerts legacy=%d port=%d delta=%d (pinned %d) seamHalfEdges=%d "
                "tLegacy=%.2fs tPort=%.2fs\n",
                ok ? "PASS" : "FAIL", c.model, c.cap, nf, legacyCharts, portCharts,
                c.expectedCharts, dropped, diff, legacyLocal, portLocal, portLocal - legacyLocal,
                c.localDelta, seams,
                std::chrono::duration<double>(t1 - t0).count(),
                std::chrono::duration<double>(t2 - t1).count());
    return ok ? Result::Pass : Result::Fail;
}

} // namespace


/// Test driver: exit 0 = all cases pass, 77 = no model present (CTest SKIP), else 1.
// NOLINTNEXTLINE(bugprone-exception-escape): a test driver; an escaping exception IS the failure
int main() {
    std::map<std::string, decimate::MeshData> meshes;
    int                                       failures = 0;
    int                                       ran = 0;

    // Run every case so one failure still prints the whole table.
    for (const Case &c : kCases) {
        const Result r = runCase(c, meshes);
        failures += (r == Result::Fail) ? 1 : 0;
        ran += (r == Result::Skipped) ? 0 : 1;
    }
    if (ran == 0) {
        std::printf("no model present: skipped\n");
        return 77; // CTest SKIP_RETURN_CODE
    }
    return failures == 0 ? 0 : 1;
}
