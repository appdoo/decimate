// wasm_ref <model> <out.rgb> : native run with the web defaults; prints counts and timings as JSON,
// dumps the baked map bytes to <out.rgb> for the W2 byte comparison.
#include <cstdio>
#include <fstream>
#include <string>

#include "pipeline_run.h"

int main(int argc, char **argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: wasm_ref <model> <out.rgb>\n");
        return 2;
    }
    const auto r = decimate_web::run(argv[1], {});
    std::ofstream(argv[2], std::ios::binary)
        .write(reinterpret_cast<const char *>(r.rgb.data()), static_cast<std::streamsize>(r.rgb.size()));
    // positions then uv (raw float32), for diagnosing a W2 mismatch.
    std::ofstream geo(std::string(argv[2]) + ".geo", std::ios::binary);
    geo.write(reinterpret_cast<const char *>(r.positions.data()), static_cast<std::streamsize>(r.positions.size() * 4));
    geo.write(reinterpret_cast<const char *>(r.uv.data()), static_cast<std::streamsize>(r.uv.size() * 4));
    std::printf("{\"inVerts\":%zu,\"outVerts\":%zu,\"faces\":%zu,\"charts\":%zu,\"width\":%d,\"height\":%d,"
                "\"load\":%.4f,\"reduce\":%.4f,\"chartify\":%.4f,\"bake\":%.4f}\n",
                r.inVerts, r.outVerts, r.indices.size() / 3, r.charts, r.width, r.height, r.secLoad,
                r.secReduce, r.secChartify, r.secBake);
    return 0;
}
