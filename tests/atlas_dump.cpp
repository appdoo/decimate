// This file is part of the decimate library; see decimate.h for version/license details.
/// @file atlas_dump.cpp
/// @brief Diagnostic: render the packed UV atlas of a model as a colored
///        silhouette PPM, and print packing density — for visual comparison
///        against the 2003 thesis atlas figures.  NOT a test.
///
/// Usage: atlas_dump <model.obj> <out.ppm> [resolution] [maxChartFaces]

#include <algorithm>
#include <array>
#include <cstdint>
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

// Distinct pastel palette (same flavor as the thesis chart figures).
std::array<std::uint8_t, 3> chartColor(int id) {
    static const std::array<std::array<std::uint8_t, 3>, 16> pal = {{
        {{232, 120, 120}}, {{120, 180, 232}}, {{140, 210, 140}}, {{230, 200, 120}},
        {{200, 140, 220}}, {{120, 210, 210}}, {{235, 160, 200}}, {{170, 200, 120}},
        {{160, 160, 235}}, {{220, 170, 130}}, {{120, 220, 170}}, {{210, 130, 160}},
        {{180, 220, 235}}, {{235, 220, 150}}, {{150, 190, 170}}, {{200, 180, 220}},
    }};
    return pal[static_cast<std::size_t>(id) % pal.size()];
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr,
            "usage: %s <model> <out.ppm> [res] [maxChartFaces] [reduceToVerts]\n",
            argv[0]);
        return 2;
    }
    const std::string model = argv[1];
    const std::string out = argv[2];
    const int res = (argc > 3) ? std::atoi(argv[3]) : 512;
    const int maxChartFaces = (argc > 4) ? std::atoi(argv[4]) : 64;
    const int reduceToVerts = (argc > 5) ? std::atoi(argv[5]) : 0;  // 0 = no reduce

    // The pipeline unwraps the LOW-POLY mesh: decimate first, then chartify the
    // reduced model (the 2003 driver did exactly this before CChartification).
    MeshData data = load(model);
    if (reduceToVerts > 0) {
        data = decimate::reduce(data, reduceToVerts);
    }
    const int numFaces = static_cast<int>(data.indices.size()) / 3;
    const int np = static_cast<int>(data.positions.size());
    const HalfEdgeMesh hem = HalfEdgeMesh::create(
        np, numFaces, reinterpret_cast<const Face*>(data.indices.data()));
    const ChartSet cs = segment(hem, data, maxChartFaces);
    const std::vector<decimate::math::Vec2<float>> uv = parametrizeCorners(hem, data, cs, 256);

    // Rasterise every UV triangle, tagged by chart color; count coverage.
    std::vector<std::array<std::uint8_t, 3>> img(
        static_cast<std::size_t>(res) * res, {{255, 255, 255}});
    std::vector<std::uint8_t> covered(static_cast<std::size_t>(res) * res, 0);
    auto edge = [](float px, float py, float qx, float qy, float rx, float ry) {
        return (qx - px) * (ry - py) - (qy - py) * (rx - px);
    };
    for (int f = 0; f < numFaces; ++f) {
        const int chart = cs.faceToChart[static_cast<std::size_t>(f)];
        const auto col = chartColor(chart);
        const auto& a = uv[static_cast<std::size_t>(3 * f)];
        const auto& b = uv[static_cast<std::size_t>(3 * f + 1)];
        const auto& c = uv[static_cast<std::size_t>(3 * f + 2)];
        const float ax = a.x * res, ay = a.y * res;
        const float bx = b.x * res, by = b.y * res;
        const float cx = c.x * res, cy = c.y * res;
        const int xlo = std::max(0, static_cast<int>(std::min({ax, bx, cx})));
        const int xhi = std::min(res - 1, static_cast<int>(std::max({ax, bx, cx})) + 1);
        const int ylo = std::max(0, static_cast<int>(std::min({ay, by, cy})));
        const int yhi = std::min(res - 1, static_cast<int>(std::max({ay, by, cy})) + 1);
        for (int y = ylo; y <= yhi; ++y) {
            for (int x = xlo; x <= xhi; ++x) {
                const float px = static_cast<float>(x) + 0.5F;
                const float py = static_cast<float>(y) + 0.5F;
                const float w0 = edge(bx, by, cx, cy, px, py);
                const float w1 = edge(cx, cy, ax, ay, px, py);
                const float w2 = edge(ax, ay, bx, by, px, py);
                const bool inside = (w0 >= 0 && w1 >= 0 && w2 >= 0) ||
                                    (w0 <= 0 && w1 <= 0 && w2 <= 0);
                if (!inside) { continue; }
                const std::size_t idx = static_cast<std::size_t>(y) * res + x;
                img[idx] = col;
                covered[idx] = 1;
            }
        }
    }

    // Flip vertically so v=0 is at the bottom (atlas convention) and write PPM.
    std::FILE* fp = std::fopen(out.c_str(), "wb");
    if (!fp) { std::fprintf(stderr, "cannot open %s\n", out.c_str()); return 1; }
    std::fprintf(fp, "P6\n%d %d\n255\n", res, res);
    for (int y = res - 1; y >= 0; --y) {
        for (int x = 0; x < res; ++x) {
            std::fwrite(img[static_cast<std::size_t>(y) * res + x].data(), 1, 3, fp);
        }
    }
    std::fclose(fp);

    long cov = 0;
    for (const auto v : covered) { cov += v; }
    const double density = 100.0 * static_cast<double>(cov) /
                           (static_cast<double>(res) * res);
    std::printf("%-24s faces=%-6d charts=%-3zu density=%.1f%%  -> %s\n",
                model.c_str(), numFaces, cs.charts.size(), density, out.c_str());

    // Also write a sidecar summary (shell stdout is unreliable in this env).
    if (std::FILE* sf = std::fopen((out + ".txt").c_str(), "wb")) {
        std::fprintf(sf, "%s faces=%d charts=%zu density=%.1f%%\n",
                     model.c_str(), numFaces, cs.charts.size(), density);
        std::fclose(sf);
    }
    return 0;
}
