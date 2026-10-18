// This file is part of the decimate library; see decimate.h for version/license details.
/// @file reduce_bench.cpp
/// @brief Head-to-head benchmark: QEM edge-collapse vs octree clustering. NOT a test.
///
/// Reports best-of-N wall time, throughput, and a fair quality proxy for both
/// reducers across target ratios and meshes, so the speed/quality trade-off is
/// measured, not guessed.
///
/// ## Fair comparison
/// Clustering's target is only approximate (cells, not an exact count), so at a
/// given ratio it overshoots QEM's vertex count. To compare quality honestly we
/// **match the triangle budget**: run QEM at the ratio, then binary-search the
/// clustering target until its output vertex count matches QEM's, and compare at
/// that matched size.
///
/// ## Quality proxy
/// Mean & max distance from a sample of input vertices to the nearest output
/// **triangle** (point-to-surface, not point-to-vertex — vertex distance would
/// unfairly penalize QEM for correctly deleting redundant coplanar vertices),
/// normalized by the input bbox diagonal and printed as a percentage.
///
/// Usage:
///   reduce_bench [--reps N] [model.obj ...]
/// With no model args it benches model/spot.obj plus any of cow/dragon/Armadillo
/// present under model/.
///
/// ## Profiling QEM (sampler)
///   macOS:  xcrun xctrace record --template "Time Profiler" --launch -- \
///             ./reduce_bench --reps 200 model/dragon.obj
///   Linux:  perf record -g ./reduce_bench --reps 200 model/dragon.obj && perf report

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include <decimate/types/mesh_data.h>
#include <decimate/math/vec3.h>
#include <decimate/reduce/reducer.h>
#include <decimate/reduce/cluster_reduce.h>
#include <decimate/reduce/halfedge_reduce.h>

using decimate::MeshData;
using decimate::math::Vec3f;

// Both reducers share this signature, so the bench drives them through one path.
using ReduceFn = MeshData (*)(const MeshData&, int);

namespace {

double msSince(std::chrono::steady_clock::time_point t0) {
    const auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

// Squared distance from point p to triangle (a,b,c) — Ericson, Real-Time
// Collision Detection §5.1.5. Guard-clause Voronoi-region tests keep it flat.
double pointTriDist2(const Vec3f& p, const Vec3f& a, const Vec3f& b, const Vec3f& c) noexcept {
    const Vec3f ab = b - a;
    const Vec3f ac = c - a;
    const Vec3f ap = p - a;
    const float d1 = ab.dot(ap);
    const float d2 = ac.dot(ap);
    if (d1 <= 0.0F && d2 <= 0.0F) {
        return ap.lengthSquared();  // vertex a
    }
    const Vec3f bp = p - b;
    const float d3 = ab.dot(bp);
    const float d4 = ac.dot(bp);
    if (d3 >= 0.0F && d4 <= d3) {
        return bp.lengthSquared();  // vertex b
    }
    const Vec3f cp = p - c;
    const float d5 = ab.dot(cp);
    const float d6 = ac.dot(cp);
    if (d6 >= 0.0F && d5 <= d6) {
        return cp.lengthSquared();  // vertex c
    }
    const float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0F && d1 >= 0.0F && d3 <= 0.0F) {
        const float v = d1 / (d1 - d3);
        return (a + ab * v - p).lengthSquared();  // edge ab
    }
    const float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0F && d2 >= 0.0F && d6 <= 0.0F) {
        const float w = d2 / (d2 - d6);
        return (a + ac * w - p).lengthSquared();  // edge ac
    }
    const float va = d3 * d6 - d5 * d4;
    if (va <= 0.0F && (d4 - d3) >= 0.0F && (d5 - d6) >= 0.0F) {
        const float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (b + (c - b) * w - p).lengthSquared();  // edge bc
    }
    const float denom = 1.0F / (va + vb + vc);
    return (a + ab * (vb * denom) + ac * (vc * denom) - p).lengthSquared();  // interior
}

struct Quality {
    double meanFrac = 0.0;
    double maxFrac  = 0.0;
};

// Point-to-surface error of `out` relative to `in`, sampled and normalized by
// the input bbox diagonal. Brute force over a capped sample — a measurement,
// off the timed path.
Quality measureQuality(const MeshData& in, const MeshData& out) {
    const int outTris = static_cast<int>(out.indices.size()) / 3;
    if (outTris == 0 || in.positions.empty()) {
        return {};
    }
    Vec3f lo = in.positions[0];
    Vec3f hi = in.positions[0];
    for (const Vec3f& p : in.positions) {
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
    const double diag   = std::max(1e-12, static_cast<double>((hi - lo).length()));
    const int    inV    = static_cast<int>(in.positions.size());
    const int    stride = std::max(1, inV / 500);  // cap the sample at ~500

    double sum = 0.0;
    double mx  = 0.0;
    int    cnt = 0;
    for (int i = 0; i < inV; i += stride) {
        const Vec3f q    = in.positions[i];
        double      best = std::numeric_limits<double>::max();
        for (int t = 0; t < outTris; ++t) {
            const std::size_t b  = static_cast<std::size_t>(t) * 3;
            const double      d2 = pointTriDist2(q, out.positions[out.indices[b]],
                                                 out.positions[out.indices[b + 1]],
                                                 out.positions[out.indices[b + 2]]);
            best = std::min(best, d2);
        }
        const double d = std::sqrt(best) / diag;
        sum += d;
        mx = std::max(mx, d);
        ++cnt;
    }
    return {sum / cnt, mx};
}

struct Row {
    int     outVerts = 0;
    int     outFaces = 0;
    double  bestMs   = 0.0;
    Quality quality;
};

// Best-of-N (min = least OS noise) after a warm-up, for one strategy at a target.
Row benchOne(ReduceFn fn, const MeshData& base, int target, int reps) {
    Row r;
    const MeshData warm = fn(base, target);  // warm-up (caches, branch predictor)
    r.outVerts = static_cast<int>(warm.positions.size());
    r.outFaces = static_cast<int>(warm.indices.size()) / 3;
    r.quality  = measureQuality(base, warm);

    r.bestMs = std::numeric_limits<double>::max();
    for (int i = 0; i < reps; ++i) {
        const auto     t0  = std::chrono::steady_clock::now();
        const MeshData out = fn(base, target);
        r.bestMs = std::min(r.bestMs, msSince(t0));
        if (out.indices.empty()) { std::fprintf(stderr, "empty result\n"); }
    }
    return r;
}

// Binary-search the clustering target so its output vertex count matches QEM's,
// for an equal-triangle-budget comparison. (outV rises with target.)
int matchedClusterTarget(const MeshData& base, int desiredOutV) {
    int lo      = 1;
    int hi      = static_cast<int>(base.positions.size());
    int bestT   = 1;
    int bestErr = std::numeric_limits<int>::max();
    for (int it = 0; it < 20 && lo <= hi; ++it) {
        const int mid = (lo + hi) / 2;
        const int got = static_cast<int>(decimate::reduceClustering(base, mid).positions.size());
        const int err = std::abs(got - desiredOutV);
        if (err < bestErr) {
            bestErr = err;
            bestT   = mid;
        }
        if (got > desiredOutV) {
            hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }
    return bestT;
}

void printRow(const char* tag, int inF, double ratio, const Row& r) {
    const double secs  = r.bestMs / 1000.0;
    const double mtris = secs > 0 ? (static_cast<double>(inF) / secs) / 1e6 : 0.0;
    std::printf("%-8s %-7.2f %-8d %-8d %-10.2f %-10.2f %-9.3f %-9.3f\n",
                tag, ratio, r.outVerts, r.outFaces, r.bestMs, mtris,
                r.quality.meanFrac * 100.0, r.quality.maxFrac * 100.0);
}

void benchModel(const std::string& path, int reps) {
    MeshData base;
    try {
        base = decimate::load(path);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "skip %s: %s\n", path.c_str(), e.what());
        return;
    }
    const int inV = static_cast<int>(base.positions.size());
    const int inF = static_cast<int>(base.indices.size()) / 3;
    if (inF == 0) { std::fprintf(stderr, "skip %s: no faces\n", path.c_str()); return; }

    std::printf("\n=== %s  (%d verts, %d faces, reps=%d) ===\n", path.c_str(), inV, inF, reps);
    std::printf("(CLUSTER rows are matched to QEM's vertex count for equal-budget quality)\n");
    std::printf("%-8s %-7s %-8s %-8s %-10s %-10s %-9s %-9s\n",
                "strategy", "ratio", "outV", "outF", "best_ms", "Mtris/s", "qMean%", "qMax%");

    for (double ratio : {0.50, 0.25, 0.10, 0.05, 0.01}) {
        const int qemTarget = std::max(4, static_cast<int>(static_cast<double>(inV) * ratio));
        const Row q         = benchOne(decimate::reduce, base, qemTarget, reps);
        printRow("QEM", inF, ratio, q);

        const int clTarget = matchedClusterTarget(base, q.outVerts);
        const Row c        = benchOne(decimate::reduceClustering, base, clTarget, reps);
        printRow("CLUSTER", inF, ratio, c);

        // Half-edge in-place collapse (modernized 2003 CHalfedgeCollapse).
        const Row he = benchOne(decimate::reduceHalfEdge, base, qemTarget, reps);
        printRow("HALFEDGE", inF, ratio, he);
    }
}

}  // namespace

int main(int argc, char** argv) {
    int reps = 7;
    std::vector<std::string> models;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--reps") == 0 && i + 1 < argc) {
            reps = std::max(1, std::atoi(argv[++i]));
        } else {
            models.emplace_back(argv[i]);
        }
    }

    if (models.empty()) {
        const std::string dir = MESH_MODEL_DIR;
        models.push_back(dir + "/spot.obj");  // always present (CC0)
        for (const char* m : {"cow.obj", "dragon.obj", "Armadillo.obj"}) {
            models.push_back(dir + "/" + m);  // benchModel skips if absent
        }
    }

    for (const std::string& m : models) { benchModel(m, reps); }
    return 0;
}
