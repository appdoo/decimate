// This file is part of the decimate library; see decimate.h for version/license details.
// Golden UV capture + check for the LSCM + atlas pipeline (Spec 09-01 Task 1).
//
// Loads model/spot.obj and runs:
//   load -> HalfEdgeMesh::create -> segment(64) -> parametrizeCorners(256)
//
// The result is compared by TOLERANCE, never by byte hash: the same source compiled at
// -O0 and -O2 reorders float arithmetic, so the raw UV bytes differ by a few ulps while
// the geometry is identical.  A byte hash is therefore not reproducible across build types.
//
// Fixture format (text, tests/fixtures/spot_uv_golden.txt):
//   '#' lines   comments
//   one line per sampled face:  <face> <u0> <v0> <u1> <v1> <u2> <v2>   (%.9g, round-trips float)
// Only every `stride`-th face is stored (default 8, ~700 faces) to keep the fixture small;
// the sampling is deterministic (face % stride == 0); a chart smaller than the
// stride may fall between samples.
//
// Usage:
//   capture_uv_golden <out_path> [stride]            (re)generate the fixture
//   capture_uv_golden --check <golden> <eps>         compare a fresh run against the fixture;
//                                                    exit 0 iff every |dUV| <= eps
//
// Regenerate (use a Debug build so the fixture is the unoptimized reference):
//   cmake --build build --target capture_uv_golden
//   build/tests/capture_uv_golden tests/fixtures/spot_uv_golden.txt

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <decimate/chart/chart.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>
#include <decimate/chart/parametrize.h>

#ifndef MESH_MODEL_DIR
#  error "MESH_MODEL_DIR must be defined by CMake"
#endif

namespace {

/// @brief Run the pipeline on spot.obj; returns per-corner UVs (3 per face) or empty on failure.
std::vector<decimate::math::Vec2<float>> computeSpotUV() {
    const std::string spotPath = std::string(MESH_MODEL_DIR) + "/spot.obj";
    decimate::MeshData data;
    try {
        data = decimate::load(spotPath);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Failed to load %s: %s\n", spotPath.c_str(), e.what());
        return {};
    }
    if (data.indices.empty()) {
        std::fprintf(stderr, "spot.obj produced an empty mesh\n");
        return {};
    }
    const int numFaces  = static_cast<int>(data.indices.size()) / 3;
    const int numPoints = static_cast<int>(data.positions.size());
    const decimate::HalfEdgeMesh hem = decimate::HalfEdgeMesh::create(
        numPoints, numFaces,
        reinterpret_cast<const decimate::Face*>(data.indices.data()));
    try {
        const auto chartSet = decimate::chart::segment(hem, data, 64);
        auto uv = decimate::chart::parametrizeCorners(hem, data, chartSet, 256);
        if (uv.size() != static_cast<std::size_t>(3 * numFaces)) {
            std::fprintf(stderr, "cornerUV size mismatch: expected %d, got %zu\n",
                         3 * numFaces, uv.size());
            return {};
        }
        return uv;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "pipeline failed: %s\n", e.what());
        return {};
    }
}

int writeFixture(const char* outPath, int stride) {
    const auto uv = computeSpotUV();
    if (uv.empty()) {
        return 1;
    }
    FILE* fp = std::fopen(outPath, "w");
    if (fp == nullptr) {
        std::fprintf(stderr, "Cannot open output file: %s\n", outPath);
        return 1;
    }
    const int numFaces = static_cast<int>(uv.size()) / 3;
    std::fprintf(fp, "# spot.obj per-corner UV golden: segment(64) + parametrizeCorners(256)\n");
    std::fprintf(fp, "# faces=%d stride=%d ; line = face u0 v0 u1 v1 u2 v2\n", numFaces, stride);
    for (int f = 0; f < numFaces; f += stride) {
        const auto* c = &uv[static_cast<std::size_t>(3 * f)];
        std::fprintf(fp, "%d %.9g %.9g %.9g %.9g %.9g %.9g\n", f,
                     static_cast<double>(c[0].x), static_cast<double>(c[0].y),
                     static_cast<double>(c[1].x), static_cast<double>(c[1].y),
                     static_cast<double>(c[2].x), static_cast<double>(c[2].y));
    }
    std::fclose(fp);
    std::printf("Wrote %s (faces=%d, stride=%d)\n", outPath, numFaces, stride);
    return 0;
}

int checkFixture(const char* goldenPath, double eps) {
    const auto uv = computeSpotUV();
    if (uv.empty()) {
        return 1;
    }
    FILE* fp = std::fopen(goldenPath, "r");
    if (fp == nullptr) {
        std::fprintf(stderr, "Cannot open golden: %s\n", goldenPath);
        return 1;
    }
    char   line[512];
    int    rows = 0;
    double maxDelta = 0.0;
    while (std::fgets(line, sizeof line, fp) != nullptr) {
        if (line[0] == '#') {
            continue;
        }
        int    f = 0;
        double g[6];
        if (std::sscanf(line, "%d %lf %lf %lf %lf %lf %lf", &f,  // NOLINT(cert-err34-c)
                        &g[0], &g[1], &g[2], &g[3], &g[4], &g[5]) != 7 ||
            f < 0 || static_cast<std::size_t>(3 * f + 2) >= uv.size()) {
            std::fprintf(stderr, "Malformed golden line: %s", line);
            std::fclose(fp);
            return 1;
        }
        for (int k = 0; k < 3; ++k) {
            const auto& c = uv[static_cast<std::size_t>(3 * f + k)];
            maxDelta = std::max({maxDelta, std::fabs(static_cast<double>(c.x) - g[2 * k]),
                                 std::fabs(static_cast<double>(c.y) - g[2 * k + 1])});
        }
        ++rows;
    }
    std::fclose(fp);
    std::printf("UV golden: %d faces compared, max |dUV| = %.3e (eps = %.1e)\n", rows, maxDelta, eps);
    if (rows == 0 || !(maxDelta <= eps)) {
        std::fprintf(stderr, "UV golden MISMATCH (max |dUV| %.3e > eps %.1e, or empty fixture)\n",
                     maxDelta, eps);
        return 1;
    }
    return 0;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc >= 4 && std::strcmp(argv[1], "--check") == 0) {
        return checkFixture(argv[2], std::atof(argv[3]));
    }
    const char* outPath = (argc > 1) ? argv[1] : "tests/fixtures/spot_uv_golden.txt";
    const int   stride  = (argc > 2) ? std::max(1, std::atoi(argv[2])) : 8;
    return writeFixture(outPath, stride);
}
