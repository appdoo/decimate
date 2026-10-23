// This file is part of the decimate library; see decimate.h for version/license details.
/// @file seamless_sweep.cpp
/// @brief Diagnostic tool (NOT a test): how chartification drives seam cost.
///
/// For every `maxCharts` x cutter-threshold setting it runs the whole seamless pipeline
/// (reduce -> segment -> parametrize -> weld -> straighten -> bake -> erase -> seam-aware
/// LODs) and writes one CSV row per stage -- the shape of Table 1/2 of Liu et al. 2017.
/// Spec 17-03.  The 12-setting run belongs in a Release build (`build-rel`); the canonical
/// -O0 build runs a single setting from ctest (`seamless_sweep_csv`).
///
/// Usage: seamless_sweep [--mesh <file>] [--maxcharts 16,32,64,128]
///                       [--thresholds 0.3,0.5,1.0] [--texsizes 256] [--tangent]
///                       [--out <file.csv>]
///
/// Rows per setting and texture size (object space): chartify, straighten, dilate,
/// seamfill, erased, lod50, lod10.  With `--tangent`, and for every texture size after
/// the first, dilate / seamfill rows in tangent space are added (a tangent-space seam is
/// not erasable, so there is no erase or LOD row for it).  Fields that do not apply to a
/// stage are empty.

#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normals.h>
#include <decimate/bake/tangents.h>
#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/reduce/reducer.h>
#include <decimate/reduce/seam_aware_reduce.h>
#include <decimate/seam/erase.h>
#include <decimate/seam/seam_energy.h>
#include <decimate/seam/seam_topology.h>
#include <decimate/seam/straighten.h>
#include <decimate/types/mesh_data.h>
#include <exception>
#include <optional>
#include <utility>
#include <string>
#include <vector>

#ifndef MESH_MODEL_DIR
#error "MESH_MODEL_DIR must be defined by CMake"
#endif

namespace {

using decimate::MeshData;
using decimate::TexturedMesh;
using decimate::bake::NormalMap;

/// Low-poly vertex target of the reduce stage, and the atlas resolution of the first texture size.
constexpr int    kLowVertices = 1500;
constexpr int    kAtlasResolution = 256;
constexpr double kStraightenTau = 0.01;

/// CSV header; the column order is the field order of writeRow.
constexpr const char *kHeader =
    "maxCharts,cutThreshold,texSize,space,stage,targetFraction,achievedVertices,earFaces,"
    "sameNeighbourEars,"
    "protectedVertices,seamEdges,seamLength3D,seamTexels,dR,dG,dB,collapsibleSeamEdges,"
    "uncollapsibleSeamEdges,collapsedSeamEdges,revertedCharts";

/// A CSV field that may be absent (printed empty).
using Opt = std::optional<double>;

/// One CSV row; an unset field prints empty.
struct Row {
    int                                  maxCharts = 0;
    double                               threshold = 0.0;
    int                                  texSize = 0;
    const char                          *space = "object";
    const char                          *stage = "";
    Opt                                  targetFraction, achievedVertices, earFaces, sameNeighbourEars, protectedVertices, seamEdges, seamLength3D, seamTexels;
    std::optional<std::array<double, 3>> d;
    Opt                                  collapsible, uncollapsible, collapsed, reverted;
};


/// Formats one optional field with nine significant digits; empty if unset.
std::string field(const Opt &v) {
    if (!v) {
        return {};
    }
    std::array<char, 64> buf {};
    std::snprintf(buf.data(), buf.size(), "%.9g", *v);
    return buf.data();
}


/// Writes @p r as one CSV line in kHeader's column order.
void writeRow(std::FILE *out, const Row &r) {
    const auto dField = [&](std::size_t c) {
        return r.d ? field((*r.d)[c]) : std::string();
    };
    std::fprintf(out, "%d,%.9g,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", r.maxCharts,
                 r.threshold, r.texSize, r.space, r.stage, field(r.targetFraction).c_str(),
                 field(r.achievedVertices).c_str(),
                 field(r.earFaces).c_str(), field(r.sameNeighbourEars).c_str(),
                 field(r.protectedVertices).c_str(), field(r.seamEdges).c_str(),
                 field(r.seamLength3D).c_str(), field(r.seamTexels).c_str(), dField(0).c_str(),
                 dField(1).c_str(), dField(2).c_str(), field(r.collapsible).c_str(),
                 field(r.uncollapsible).c_str(), field(r.collapsed).c_str(),
                 field(r.reverted).c_str());
}


/// Everything that does not depend on the texture size.
struct Atlas {
    MeshData                           high;
    MeshData                           low;
    decimate::HalfEdgeMesh             hem;
    std::vector<int>                   faceToChart;
    TexturedMesh                       chartified;
    TexturedMesh                       straightened;
    std::vector<decimate::math::Vec3f> lowNormals;
    std::vector<decimate::math::Vec3f> highNormals;
    decimate::seam::StraightenStats    straightenStats;
};


/// Runs reduce -> segment -> parametrize -> weld -> straighten and the normals for one setting.
Atlas buildAtlas(const MeshData &high, int maxCharts, double threshold) {
    Atlas a;
    a.high = high;
    a.low = decimate::reduce(high, kLowVertices);
    std::vector<decimate::Face> faces;
    for (std::size_t i = 0; i + 2 < a.low.indices.size(); i += 3) {
        faces.push_back({ a.low.indices[i], a.low.indices[i + 1], a.low.indices[i + 2] });
    }

    a.hem = decimate::HalfEdgeMesh::create(static_cast<int>(a.low.positions.size()),
                                           static_cast<int>(faces.size()), faces.data());
    auto       charts = decimate::chart::segment(a.hem, a.low, maxCharts, threshold);
    const auto cornerUV = decimate::chart::parametrizeCorners(a.hem, a.low, charts, kAtlasResolution);

    a.faceToChart = charts.faceToChart; // after the guard may have split charts (19-03)
    a.chartified = decimate::seam::weld(a.low, cornerUV);
    auto st = decimate::seam::straighten(a.chartified, { .tau = kStraightenTau });
    a.straightened = std::move(st.mesh);
    a.straightenStats = st.stats;

    a.lowNormals = decimate::bake::computeVertexNormals(a.low, decimate::bake::Weighting::EqualUnit);
    a.highNormals =
        decimate::bake::computeVertexNormals(a.high, decimate::bake::Weighting::EqualUnit);
    return a;
}


/// One chartification setting of the sweep.
struct Setting {
    int    maxCharts = 0;
    double threshold = 0.0;
};


/// A row carrying only the setting, texture size and stage (all other fields unset).
Row keyOf(const Setting &setting, int texSize, const char *stage) {
    Row r;
    r.maxCharts = setting.maxCharts;
    r.threshold = setting.threshold;
    r.texSize = texSize;
    r.stage = stage;
    return r;
}


/// Seam count and total 3D seam length of @p tm, by `findSeams`.
void seamSummary(const TexturedMesh &tm, Row &r) {
    const auto seams = decimate::seam::findSeams(tm);
    double     length = 0.0;
    for (const auto &s : seams.seams) {
        length += s.length3D;
    }
    r.seamEdges = static_cast<double>(seams.seams.size());
    r.seamLength3D = length;
}


/// Fills the collapsible / uncollapsible seam-edge counts of @p r from @p tm.
void classifyInto(const TexturedMesh &tm, Row &r) {
    const auto c = decimate::classifySeamEdges(tm);
    r.collapsible = c.collapsible;
    r.uncollapsible = c.uncollapsible;
}


/// Writes the chartify and straighten rows.
void emitAtlasRows(std::FILE *out, const Atlas &a, const Row &key) {
    Row chart = key;
    chart.stage = "chartify";
    chart.earFaces = a.straightenStats.earFaces;
    chart.sameNeighbourEars = a.straightenStats.sameNeighbourEars;
    seamSummary(a.chartified, chart);
    classifyInto(a.chartified, chart);
    writeRow(out, chart);

    Row straight = key;
    straight.stage = "straighten";
    straight.protectedVertices = a.straightenStats.protectedVertices;
    straight.reverted = a.straightenStats.revertedCharts;
    seamSummary(a.straightened, straight);
    classifyInto(a.straightened, straight);
    writeRow(out, straight);
}


/// MikkTSpace tangents of the straightened atlas, or none for object space.
std::vector<decimate::math::Vec4f> tangentsFor(const Atlas &a, bool tangent) {
    if (!tangent) {
        return {};
    }
    return decimate::tangents(a.low, a.straightened.cornerUV(), a.lowNormals);
}


/// Normal map with blur dilation (2 passes) on the straightened atlas.
NormalMap bakeDilate(const Atlas &a, int texSize, bool tangent) {
    return decimate::bake::bake(a.low, a.straightened.cornerUV(), a.lowNormals, a.high,
                                a.highNormals, texSize, 2, tangentsFor(a, tangent));
}


/// Normal map with seam-correct gutter fill (2 passes) on the straightened atlas.
NormalMap bakeFill(const Atlas &a, int texSize, bool tangent) {
    return decimate::bake::bakeSeamFill(a.low, a.hem, a.straightened.cornerUV(), a.faceToChart,
                                        a.lowNormals, a.high, a.highNormals, texSize, 2, tangentsFor(a, tangent));
}


/// Writes the dilate and seamfill rows: seam energy of each gutter strategy.
void emitBakeRows(std::FILE *out, const Atlas &a, const Row &key, bool tangent) {
    const char *space = tangent ? "tangent" : "object";
    Row         dilate = key;
    dilate.stage = "dilate";
    dilate.space = space;
    dilate.d = decimate::seam::seamEnergy(a.straightened, bakeDilate(a, key.texSize, tangent));
    writeRow(out, dilate);

    Row fill = key;
    fill.stage = "seamfill";
    fill.space = space;
    fill.d = decimate::seam::seamEnergy(a.straightened, bakeFill(a, key.texSize, tangent));
    writeRow(out, fill);
}


/// The lod50 / lod10 row of one seam-aware reduction.
Row lodRow(const Row &key, int percent, const decimate::SeamAwareResult &res,
           const std::array<double, 3> &d) {
    Row r = key;
    r.stage = percent == 50 ? "lod50" : "lod10";
    // The reducer stops above the target when no legal collapse remains: the target is what
    // was asked for, achievedVertices what the output has.

    r.targetFraction = percent / 100.0;
    r.achievedVertices = static_cast<double>(res.mesh.positions.size());
    seamSummary(res.mesh, r);
    // The reducer's own count is authoritative for the LOD; the summary above only supplies the length.
    r.seamEdges = res.stats.seamEdgesAfter;
    r.d = d;
    r.collapsible = res.stats.collapsibleSeamEdges;
    r.uncollapsible = res.stats.uncollapsibleSeamEdges;
    r.collapsed = res.stats.collapsedSeamEdges;
    return r;
}


/// erased + lod rows (object space).
void emitErasedRows(std::FILE *out, const Atlas &a, const Row &key) {
    const NormalMap map = bakeDilate(a, key.texSize, false);
    const auto      erased =
        decimate::seam::erase(a.straightened, map, decimate::seam::NormalSpace::Object);
    Row er = key;
    er.stage = "erased";
    er.d = erased.dAfter;
    er.seamTexels = erased.unknownTexels;
    writeRow(out, er);

    // LODs are reduced from the straightened atlas, and measured against the same erased texels.
    const int nV = static_cast<int>(a.straightened.positions.size());
    for (const int percent : { 50, 10 }) {
        const auto lod = decimate::reduceSeamAware(a.straightened, nV * percent / 100);
        const auto d = decimate::seam::seamEnergy(lod.mesh, key.texSize, key.texSize, erased.texels);
        writeRow(out, lodRow(key, percent, lod, d));
    }
}


/// Command-line options with the sweep's defaults.
struct Options {
    std::string         mesh = std::string(MESH_MODEL_DIR) + "/bunny.iv";
    std::vector<int>    maxCharts { 16, 32, 64, 128 };
    std::vector<double> thresholds { 0.3, 0.5, 1.0 };
    std::vector<int>    texSizes { kAtlasResolution };
    bool                tangent = false;
    std::string         out;
};


/// Splits a comma-separated list (an empty field stays an empty string).
std::vector<std::string> split(const std::string &s) {
    std::vector<std::string> parts;
    std::size_t              begin = 0;
    while (begin <= s.size()) {
        const std::size_t end = s.find(',', begin);
        parts.push_back(s.substr(begin, end == std::string::npos ? end : end - begin));
        if (end == std::string::npos) {
            break;
        }
        begin = end + 1;
    }
    return parts;
}


/// Parses a comma-separated int list; throws std::invalid_argument on a bad field.
std::vector<int> intList(const std::string &s) {
    std::vector<int> out;
    for (const std::string &part : split(s)) {
        out.push_back(std::stoi(part));
    }
    return out;
}


/// Parses a comma-separated double list; throws std::invalid_argument on a bad field.
std::vector<double> doubleList(const std::string &s) {
    std::vector<double> out;
    for (const std::string &part : split(s)) {
        out.push_back(std::stod(part));
    }
    return out;
}


/// Applies `--name value`; false for an unknown option.
bool applyOption(Options &o, const std::string &name, const std::string &value) {
    if (name == "--mesh") {
        o.mesh = value;
    } else if (name == "--out") {
        o.out = value;
    } else if (name == "--maxcharts") {
        o.maxCharts = intList(value);
    } else if (name == "--thresholds") {
        o.thresholds = doubleList(value);
    } else if (name == "--texsizes") {
        o.texSizes = intList(value);
    } else {
        return false;
    }
    return true;
}


/// False (after printing usage) on a malformed command line.
bool parseArgs(int argc, char **argv, Options &o) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--tangent") {
            o.tangent = true;
        } else if (i + 1 < argc && applyOption(o, arg, argv[i + 1])) {
            ++i;
        } else {
            std::fprintf(stderr,
                         "usage: %s [--mesh f] [--maxcharts a,b] [--thresholds a,b] [--texsizes a,b] "
                         "[--tangent] [--out f.csv]\n",
                         argv[0]);
            return false;
        }
    }
    return true;
}


/// Builds the atlas once for @p setting, then emits all rows for every texture size.
void runSetting(std::FILE *out, const MeshData &high, const Options &o, const Setting &setting) {
    const Atlas atlas = buildAtlas(high, setting.maxCharts, setting.threshold);
    for (std::size_t t = 0; t < o.texSizes.size(); ++t) {
        const Row key = keyOf(setting, o.texSizes[t], "");
        emitAtlasRows(out, atlas, key);
        emitBakeRows(out, atlas, key, false);
        emitErasedRows(out, atlas, key);
        if (o.tangent || t > 0) {
            emitBakeRows(out, atlas, key, true);
        }
        std::fflush(out);
    }
}


/// The whole sweep; exit code 0 on success, 2 for a usage or I/O error.
int run(int argc, char **argv) {
    Options o;
    if (!parseArgs(argc, argv, o)) {
        return 2;
    }
    const MeshData high = decimate::load(o.mesh);
    std::FILE     *out = o.out.empty() ? stdout : std::fopen(o.out.c_str(), "wb");
    if (out == nullptr) {
        std::fprintf(stderr, "cannot open %s\n", o.out.c_str());
        return 2;
    }

    std::fprintf(out, "%s\n", kHeader);
    for (const int maxCharts : o.maxCharts) {
        for (const double threshold : o.thresholds) {
            std::fprintf(stderr, "[sweep] maxCharts %d threshold %g\n", maxCharts, threshold);
            runSetting(out, high, o, { maxCharts, threshold });
        }
    }
    if (out != stdout) {
        std::fclose(out);
    }
    return 0;
}

} // namespace


/// Entry point: any exception from the pipeline becomes exit code 1.
int main(int argc, char **argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception &e) {
        std::fprintf(stderr, "seamless_sweep: %s\n", e.what());
        return 1;
    }
}
