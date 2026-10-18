// This file is part of the decimate library; see decimate.h for version/license details.
/// @file reduce_chartify_timing.cpp
/// @brief Diagnostic: time the reduce→chartify pipeline at several target
///        resolutions, to show that chartification cost tracks the REDUCED
///        mesh size (the 2003 interactive pipeline decimated first), not the
///        raw mesh size.  NOT a test.
///
/// Usage: reduce_chartify_timing <model> [targetVerts1 targetVerts2 ...]
/// Default targets: raw (0 = no reduce), 10000, 5000, 1000 vertices.

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <decimate/chart/chart.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>
#include <decimate/chart/parametrize.h>
#include <decimate/reduce/reducer.h>

using namespace decimate;
using namespace decimate::chart;

namespace {

double msSince(std::chrono::steady_clock::time_point t0) {
    const auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

struct Row {
    int    targetVerts = 0;
    int    faces       = 0;
    int    charts      = 0;
    double reduceMs    = 0.0;
    double segMs       = 0.0;
    double paramMs     = 0.0;
};

Row runOnce(const MeshData& base, int targetVerts) {
    Row r;
    r.targetVerts = targetVerts;

    MeshData data = base;
    if (targetVerts > 0) {
        const auto t0 = std::chrono::steady_clock::now();
        data = decimate::reduce(base, targetVerts);
        r.reduceMs = msSince(t0);
    }
    r.faces = static_cast<int>(data.indices.size()) / 3;
    if (r.faces == 0) { return r; }

    const int np = static_cast<int>(data.positions.size());
    const HalfEdgeMesh hem = HalfEdgeMesh::create(
        np, r.faces, reinterpret_cast<const Face*>(data.indices.data()));

    const auto t1 = std::chrono::steady_clock::now();
    const ChartSet cs = segment(hem, data, 64);
    r.segMs = msSince(t1);
    r.charts = static_cast<int>(cs.charts.size());

    const auto t2 = std::chrono::steady_clock::now();
    const std::vector<decimate::math::Vec2<float>> uv = parametrizeCorners(hem, data, cs, 256);
    r.paramMs = msSince(t2);
    (void)uv;

    return r;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <model> [targetVerts...]\n", argv[0]);
        return 2;
    }
    const std::string model = argv[1];

    std::vector<int> targets;
    if (argc > 2) {
        for (int i = 2; i < argc; ++i) { targets.push_back(std::atoi(argv[i])); }
    } else {
        targets = {0, 10000, 5000, 1000};  // 0 = raw, no reduction
    }

    const MeshData base = load(model);
    const int baseFaces = static_cast<int>(base.indices.size()) / 3;
    const int baseVerts = static_cast<int>(base.positions.size());

    std::printf("\n=== %s  (raw: %d verts, %d faces) ===\n",
                model.c_str(), baseVerts, baseFaces);
    std::printf("%-12s %-8s %-7s %-10s %-10s %-10s %-10s\n",
                "targetVerts", "faces", "charts", "reduce_ms", "seg_ms",
                "param_ms", "total_ms");

    std::FILE* sf = std::fopen("/tmp/reduce_timing.txt", "ab");
    for (int t : targets) {
        const Row r = runOnce(base, t);
        const double total = r.reduceMs + r.segMs + r.paramMs;
        std::printf("%-12d %-8d %-7d %-10.1f %-10.1f %-10.1f %-10.1f\n",
                    r.targetVerts, r.faces, r.charts, r.reduceMs, r.segMs,
                    r.paramMs, total);
        if (sf != nullptr) {
            std::fprintf(sf,
                "%s target=%d faces=%d charts=%d reduce=%.1f seg=%.1f param=%.1f total=%.1f\n",
                model.c_str(), r.targetVerts, r.faces, r.charts, r.reduceMs,
                r.segMs, r.paramMs, total);
        }
    }
    if (sf != nullptr) { std::fclose(sf); }
    return 0;
}
