// This file is part of the decimate library; see decimate.h for version/license details.
/// @file raster_bench.cpp
/// @brief Real-workload benchmark of the atlas packer, driven from a fixture.
///        Diagnostic, NOT a unit test.
///
/// The slow pipeline (load cow -> reduce -> segment -> flatten) is run ONCE and
/// its output — the reduced MeshData + the flattened HalfEdgeCharts that feed
/// `packAtlas` — is cached to `tests/fixtures/pack_bench.fixture`.  On every run
/// the bench just loads that fixture and times ONLY `packAtlas` in a tight loop.
/// packAtlas internally drives the chart-silhouette rasterizer (Andor's
/// `rasterLine`) across 16 rotations x all charts x all edges — the packer's
/// real hot path.  Delete the fixture to regenerate it.
///
/// packAtlas mutates each chart's `uv` in place, so a pristine copy is restored
/// before every timed rep.

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include <decimate/chart/chart.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>
#include <decimate/chart/parametrize.h>
#include <decimate/reduce/reducer.h>

#ifndef MESH_MODEL_DIR
#error "MESH_MODEL_DIR must be defined by CMake"
#endif
#ifndef MESH_FIXTURES_DIR
#error "MESH_FIXTURES_DIR must be defined by CMake"
#endif

using namespace decimate;
using namespace decimate::chart;

namespace {

/// Low-poly target: the pipeline decimates before it unwraps.
constexpr int kReduceTarget = 1000;

/// Cached pipeline output the bench times against; delete it to regenerate.
const std::string kFixture = std::string(MESH_FIXTURES_DIR) + "/pack_bench.fixture";

/// Milliseconds elapsed since @p t0.
double msSince(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now() - t0)
        .count();
}


/// Writes a vector as `count v0 v1 ...` (fixture format: plain text, deterministic, readable).
template <typename T>
void writeVec(std::ofstream &os, const std::vector<T> &v) {
    os << v.size();
    for (const T &x : v) {
        os << ' ' << x;
    }
    os << '\n';
}


/// Writes UVs as `count u0 v0 u1 v1 ...`.
void writeUVs(std::ofstream &os, const std::vector<decimate::math::Vec2f> &v) {
    os << v.size();
    for (const auto &p : v) {
        os << ' ' << p.x << ' ' << p.y;
    }
    os << '\n';
}


/// Writes positions in the flat-float on-disk format: the leading count is the FLOAT count (3 per
/// vertex), so fixtures written before the Vec3 storage change still load.
void writePositions(std::ofstream &os, const std::vector<decimate::math::Vec3f> &v) {
    os << v.size() * 3;
    for (const auto &p : v) {
        os << ' ' << p.x << ' ' << p.y << ' ' << p.z;
    }
    os << '\n';
}


/// Reads what writeVec wrote.
template <typename T>
std::vector<T> readVec(std::ifstream &is) {
    std::size_t n = 0;
    is >> n;
    std::vector<T> v(n);
    for (std::size_t i = 0; i < n; ++i) {
        is >> v[i];
    }
    return v;
}


/// Reads what writePositions wrote.
std::vector<decimate::math::Vec3f> readPositions(std::ifstream &is) {
    std::size_t n = 0; // float count (3 per vertex)
    is >> n;
    std::vector<decimate::math::Vec3f> v(n / 3);
    for (auto &p : v) {
        is >> p.x >> p.y >> p.z;
    }
    return v;
}


/// Reads what writeUVs wrote.
std::vector<decimate::math::Vec2f> readUVs(std::ifstream &is) {
    std::size_t n = 0;
    is >> n;
    std::vector<decimate::math::Vec2f> v(n);
    for (std::size_t i = 0; i < n; ++i) {
        is >> v[i].x >> v[i].y;
    }
    return v;
}


/// The cached input of packAtlas: the reduced mesh and its flattened charts.
struct Fixture {
    MeshData                   data;
    std::vector<HalfEdgeChart> charts;
};


/// @brief Rebuild the local HalfEdgeMesh for a HalfEdgeChart from serialized state.
///
/// @p localTri is the chart's local triangulation: 3 local vertices per chart face, read
/// from the flattened chart with `localCorner`.  A self-seam cut chart cannot be rebuilt from
/// `localToMesh` (many-to-one), so the triangulation itself is stored in the fixture.
void rebuildChartMesh(HalfEdgeChart &c, const std::vector<int> &localTri) {
    c.mesh = HalfEdgeMesh::create(
        static_cast<int>(c.localToMesh.size()), static_cast<int>(c.faces.size()),
        reinterpret_cast<const Face *>(localTri.data()));
}


/// The chart's local triangulation: 3 local vertices per chart face (see `rebuildChartMesh`).
std::vector<int> localTriangulation(const HalfEdgeChart &c) {
    std::vector<int> localTri;
    for (int i = 0; i < static_cast<int>(c.faces.size()); ++i) {
        for (int k = 0; k < 3; ++k) {
            localTri.push_back(localCorner(c, i, k));
        }
    }
    return localTri;
}


/// Serializes the fixture to @p path.
void writeFixture(const Fixture &f, const std::string &path) {
    std::ofstream os(path);
    writePositions(os, f.data.positions);
    writeVec(os, f.data.indices);
    os << f.charts.size() << '\n';
    for (const HalfEdgeChart &c : f.charts) {
        writeVec(os, c.faces);
        writeUVs(os, c.uv);
        writeVec(os, c.localToMesh);
        writeVec(os, localTriangulation(c));
    }
}


/// Loads the fixture from @p path; false if it is missing or truncated.
bool readFixture(Fixture &f, const std::string &path) {
    std::ifstream is(path);
    if (!is) {
        return false;
    }
    f.data.positions = readPositions(is);
    f.data.indices = readVec<int>(is);
    std::size_t nc = 0;
    is >> nc;
    if (!is) {
        return false;
    }
    f.charts.resize(nc);
    for (HalfEdgeChart &c : f.charts) {
        c.faces = readVec<int>(is);
        c.uv = readUVs(is);
        c.localToMesh = readVec<int>(is);
        rebuildChartMesh(c, readVec<int>(is));
    }
    return static_cast<bool>(is);
}


/// Runs the slow pipeline (load, reduce, segment, flatten) once to produce the fixture.
Fixture generate() {
    Fixture f;
    f.data = load(std::string(MESH_MODEL_DIR) + "/cow.obj");

    // Reduce first, then chartify the low-poly mesh, as the pipeline does.
    f.data = decimate::reduce(f.data, kReduceTarget);
    const int          numFaces = static_cast<int>(f.data.indices.size()) / 3;
    const int          np = static_cast<int>(f.data.positions.size());
    const HalfEdgeMesh hem = HalfEdgeMesh::create(
        np, numFaces, reinterpret_cast<const Face *>(f.data.indices.data()));

    const ChartSet cs = segment(hem, f.data, 64);
    f.charts.reserve(cs.charts.size());
    for (const Chart &c : cs.charts) {
        f.charts.push_back(flatten(hem, f.data, c));
    }
    return f;
}

} // namespace


/// Times packAtlas on the cached fixture at three atlas resolutions.
int main() {
    Fixture fx;
    if (!readFixture(fx, kFixture)) {
        std::fprintf(stderr, "[bench] fixture not found — generating %s ...\n", kFixture.c_str());
        fx = generate();
        writeFixture(fx, kFixture);
        std::fprintf(stderr, "[bench] fixture written.\n");
    }

    const int numFaces = static_cast<int>(fx.data.indices.size()) / 3;
    std::printf("\n=== packAtlas() benchmark (from fixture) ===\n");
    std::printf("cow reduced->%d : %d faces, %zu charts\n",
                kReduceTarget, numFaces, fx.charts.size());
    std::printf("(times ONLY packAtlas: 16 rotations x charts x edges -> rasterLine)\n\n");

    for (const int res : { 128, 256, 512 }) {
        constexpr int kReps = 15;
        double        total = 0.0;
        std::int64_t  sink = 0;
        for (int r = 0; r < kReps; ++r) {
            // Fresh copy per rep: packAtlas rewrites the UVs in place.
            std::vector<HalfEdgeChart> charts = fx.charts; // fresh copy (packAtlas mutates uv)
            const auto                 t0 = std::chrono::steady_clock::now();
            packAtlas(charts, fx.data, res);
            total += msSince(t0);
            // Consume the result so the optimizer cannot drop the timed call.
            for (const HalfEdgeChart &c : charts) {
                if (!c.uv.empty()) {
                    sink += static_cast<std::int64_t>(c.uv[0].x * 1000.0f);
                }
            }
        }
        std::printf("  res=%-4d : packAtlas mean %8.3f ms  (x%d reps)  [sink %lld]\n",
                    res, total / kReps, kReps, static_cast<long long>(sink));
        std::fflush(stdout);
    }
    return 0;
}
