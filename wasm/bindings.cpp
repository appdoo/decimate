// Embind surface of the Wasm build: one function, plain data out (spec 34).
//   process(fileName, bytes, {targetRatio, maxCharts, resolution, tangentSpace}) -> object
//   loadMesh(fileName, bytes) -> {positions, indices}  (the input as the core reads it)
#include <cctype>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include "pipeline_run.h"

using emscripten::val;

namespace {

// typed_memory_view aliases the Wasm heap (and picks the JS array type from T); slice() copies it
// into an array the JS side owns, safe across heap growth.
template <class T>
val copyOut(const std::vector<T> &v) {
    return val(emscripten::typed_memory_view(v.size(), v.data())).call<val>("slice");
}

template <class T>
T opt(const val &o, const char *key, T fallback) {
    if (o.isUndefined() || o.isNull() || !o.hasOwnProperty(key)) {
        return fallback;
    }
    return o[key].as<T>();
}

// The loaders dispatch on the file suffix; keep only the (lower-cased) extension, never trust the name.
std::string memfsPath(const std::string &fileName) {
    const auto dot = fileName.rfind('.');
    std::string ext = dot == std::string::npos ? "" : fileName.substr(dot);
    for (auto &c : ext) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return "/input" + ext;
}

// Write the JS bytes to MEMFS, run @p f on the path, remove the file again (also when f throws).
template <class F>
val withStagedFile(const std::string &fileName, const val &bytes, F f) {
    const std::string path = memfsPath(fileName);
    const auto        n = bytes["length"].as<std::size_t>();
    std::vector<std::uint8_t> buf(n);
    val(emscripten::typed_memory_view(n, buf.data())).call<void>("set", bytes);
    {
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        out.write(reinterpret_cast<const char *>(buf.data()), static_cast<std::streamsize>(n));
        if (!out) {
            throw std::runtime_error("cannot stage input file");
        }
    }
    try {
        val r = f(path);
        std::remove(path.c_str());
        return r;
    } catch (...) {
        std::remove(path.c_str());
        throw;
    }
}

val toJs(const decimate_web::Result &r) {
    val out = val::object();
    out.set("positions", copyOut(r.positions));
    out.set("indices", copyOut(r.indices));
    out.set("uv", copyOut(r.uv));
    out.set("normals", copyOut(r.normals));
    out.set("tangents", copyOut(r.tangents));

    val map = val::object();
    map.set("width", r.width);
    map.set("height", r.height);
    map.set("rgb", copyOut(r.rgb));
    out.set("map", map);

    val seconds = val::object();
    seconds.set("load", r.secLoad);
    seconds.set("reduce", r.secReduce);
    seconds.set("chartify", r.secChartify);
    seconds.set("bake", r.secBake);
    val stats = val::object();
    stats.set("inVerts", static_cast<unsigned>(r.inVerts));
    stats.set("outVerts", static_cast<unsigned>(r.outVerts));
    stats.set("charts", static_cast<unsigned>(r.charts));
    stats.set("erased", r.erased);
    stats.set("seconds", seconds);
    out.set("stats", stats);
    return out;
}

val process(const std::string &fileName, const val &bytes, const val &options) {
    try {
        decimate_web::Options o;
        o.targetRatio = opt(options, "targetRatio", o.targetRatio);
        o.maxCharts = opt(options, "maxCharts", o.maxCharts);
        o.resolution = opt(options, "resolution", o.resolution);
        o.tangentSpace = opt(options, "tangentSpace", o.tangentSpace);
        if (o.targetRatio <= 0.0 || o.targetRatio > 1.0 || o.maxCharts < 1 || o.resolution < 2) {
            throw std::runtime_error("invalid options");
        }

        return withStagedFile(fileName, bytes, [&](const std::string &path) { return toJs(decimate_web::run(path, o)); });
    } catch (const std::exception &e) {
        val::global("Error").new_(std::string(e.what())).throw_();
    }
    return val::undefined();
}


// The input mesh exactly as the core loads it (any supported format), for the original pane.
val loadMesh(const std::string &fileName, const val &bytes) {
    try {
        return withStagedFile(fileName, bytes, [](const std::string &path) {
            const auto       m = decimate::load(path);
            std::vector<float> positions;
            positions.reserve(m.positions.size() * 3);
            for (const auto &p : m.positions) {
                positions.insert(positions.end(), { p.x, p.y, p.z });
            }
            const std::vector<unsigned> indices(m.indices.begin(), m.indices.end());
            val out = val::object();
            out.set("positions", copyOut(positions));
            out.set("indices", copyOut(indices));
            return out;
        });
    } catch (const std::exception &e) {
        val::global("Error").new_(std::string(e.what())).throw_();
    }
    return val::undefined();
}

} // namespace

EMSCRIPTEN_BINDINGS(decimate_web) {
    emscripten::function("process", &process);
    emscripten::function("loadMesh", &loadMesh);
}
