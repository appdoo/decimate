// Shared by the Wasm bindings and the native reference binary (W2): one definition of the
// web pipeline, so wasm vs native differ only by floating-point behaviour of the target.
#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <string>
#include <vector>

#include <decimate/bake/normals.h>
#include <decimate/pipeline.h>
#include <decimate/reduce/reducer.h>
#include <decimate/types/mesh_data.h>

namespace decimate_web {

struct Options {
    double targetRatio = 0.25;
    int    maxCharts = 16;
    int    resolution = 512;
    bool   tangentSpace = false;
};

struct Result {
    std::vector<float>         positions; ///< 3 per output vertex
    std::vector<unsigned>      indices;   ///< 3 per face
    std::vector<float>         uv;        ///< 2 per corner
    std::vector<float>         normals;   ///< 3 per corner
    std::vector<float>         tangents;  ///< 4 per corner
    int                        width = 0, height = 0;
    std::vector<unsigned char> rgb;
    std::size_t                inVerts = 0, outVerts = 0, charts = 0;
    bool                       erased = false; ///< false: erasure failed (e.g. out of memory), dilate map kept
    double                     secLoad = 0, secReduce = 0, secChartify = 0, secBake = 0;
};

class Stopwatch {
public:
    double lap() {
        const auto now = std::chrono::steady_clock::now();
        const double s = std::chrono::duration<double>(now - t_).count();
        t_ = now;
        return s;
    }

private:
    std::chrono::steady_clock::time_point t_ = std::chrono::steady_clock::now();
};

/// load -> reduce(high, ratio*V) -> chartify(low) -> bakeAtlas(atlas, high). Throws std::runtime_error.
inline Result run(const std::string &path, const Options &o) {
    Result      r;
    Stopwatch   sw;
    const auto  high = decimate::load(path);
    r.secLoad = sw.lap();
    r.inVerts = high.positions.size();

    const int target = std::max(4, static_cast<int>(static_cast<double>(r.inVerts) * o.targetRatio));
    const auto low = decimate::reduce(high, target);
    r.secReduce = sw.lap();

    decimate::ChartOptions co;
    co.maxCharts = o.maxCharts;
    co.atlasResolution = o.resolution;
    const auto atlas = decimate::chartify(low, co);
    r.secChartify = sw.lap();

    decimate::BakeOptions bo;
    bo.space = o.tangentSpace ? decimate::BakeSpace::Tangent : decimate::BakeSpace::Object;
    bo.gutter = decimate::Gutter::Seamless;
    bo.resolution = o.resolution;
    const auto baked = decimate::bakeAtlas(atlas, high, bo);
    r.secBake = sw.lap();
    r.erased = baked.gutterUsed == decimate::Gutter::Erased;

    const auto &m = atlas.mesh;
    r.outVerts = m.positions.size();
    r.charts = atlas.charts.charts.size();
    for (const auto &p : m.positions) {
        r.positions.insert(r.positions.end(), { p.x, p.y, p.z });
    }
    r.indices.assign(m.indices.begin(), m.indices.end());

    const auto vn = decimate::bake::computeVertexNormals(m, decimate::bake::Weighting::EqualUnit);
    r.normals.reserve(m.indices.size() * 3);
    for (const int i : m.indices) {
        const auto &n = vn[static_cast<std::size_t>(i)];
        r.normals.insert(r.normals.end(), { n.x, n.y, n.z });
    }
    for (const auto &t : atlas.uv) {
        r.uv.insert(r.uv.end(), { t.x, t.y });
    }
    for (const auto &t : atlas.tangents) {
        r.tangents.insert(r.tangents.end(), { t.x, t.y, t.z, t.w });
    }

    r.width = baked.map.width;
    r.height = baked.map.height;
    r.rgb.reserve(baked.map.pixels.size() * 3);
    for (const auto &px : baked.map.pixels) {
        r.rgb.insert(r.rgb.end(), { px.r, px.g, px.b });
    }
    return r;
}

} // namespace decimate_web
