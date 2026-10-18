// This file is part of the decimate library; see decimate.h for version/license details.
/// @file mesh_defects.cpp
/// @brief Topological defect detector for reduced meshes. NOT a test — a diagnostic.
///
/// The point-to-surface quality metric (reduce_bench) is shape- and
/// topology-BLIND: a pinched ear or a non-manifold fold sits right on the
/// original surface, so its distance error is ~0. This tool measures the things
/// that metric misses, on the OUTPUT mesh:
///   - non-manifold edges  (an undirected edge used by != 2 faces; 1 = boundary,
///                          >=3 = non-manifold)
///   - pinch / bowtie verts (incident faces form more than one fan around v —
///                          "closed but unusual": locally two cones meeting at a
///                          point, which keeps every edge manifold yet is wrong)
///   - degenerate faces     (repeated index)
///   - spikes               (an output vertex flung far from the input surface)
/// Each defect prints its world position and normalized height y' in [0,1]
/// (0 = bbox bottom, 1 = top) so a cluster can be located on the model by eye —
/// e.g. bunny ears ~ y'>0.8, neck ~ y'~0.55.
///
/// Usage: mesh_defects [--target N | --ratio R] model.{obj,iv} [more...]

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#include <decimate/math/vec3.h>
#include <decimate/reduce/reducer.h>
#include <decimate/types/mesh_data.h>

using decimate::MeshData;
using decimate::math::Vec3f;

namespace {

// Undirected edge key (min,max) packed into 64 bits.
std::uint64_t edgeKey(int a, int b) {
    const std::uint32_t lo = static_cast<std::uint32_t>(std::min(a, b));
    const std::uint32_t hi = static_cast<std::uint32_t>(std::max(a, b));
    return (static_cast<std::uint64_t>(hi) << 32) | lo;
}

struct DefectReport {
    int   nonManifoldEdges = 0;
    int   boundaryEdges    = 0;
    int   pinchVerts       = 0;
    int   degenFaces       = 0;
    int   spikes           = 0;
    int   flippedEdges     = 0;  ///< interior edges traversed the SAME way by both faces
    int   foldEdges        = 0;  ///< interior edges where adjacent faces crease sharply back
    int   slivers          = 0;  ///< near-degenerate (thin) triangles
    Vec3f lo;                    ///< output bbox, for normalized-height reporting
    Vec3f hi;
};

// Disjoint-set over a small local index space (incident-face count).
struct DSU {
    std::vector<int> parent;
    explicit DSU(int n) : parent(static_cast<std::size_t>(n)) {
        std::iota(parent.begin(), parent.end(), 0);
    }
    int find(int x) {
        while (parent[static_cast<std::size_t>(x)] != x) {
            parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
            x = parent[static_cast<std::size_t>(x)];
        }
        return x;
    }
    void join(int a, int b) { parent[static_cast<std::size_t>(find(a))] = find(b); }
};

void printDefect(const char* kind, const Vec3f& p, const Vec3f& lo, const Vec3f& hi, int& shown) {
    if (shown >= 12) {
        return;  // cap the per-kind listing; the count is the headline
    }
    const float yspan = std::max(1e-9F, hi.y - lo.y);
    const float yn    = (p.y - lo.y) / yspan;
    std::printf("    %-14s (% .4f % .4f % .4f)  y'=%.2f\n", kind, p.x, p.y, p.z, yn);
    ++shown;
}

struct Bounds {
    Vec3f lo;
    Vec3f hi;
};
Bounds bboxOf(const std::vector<Vec3f>& pts) {
    Bounds b{pts.empty() ? Vec3f{} : pts[0], pts.empty() ? Vec3f{} : pts[0]};
    for (const Vec3f& p : pts) {
        b.lo = {std::min(b.lo.x, p.x), std::min(b.lo.y, p.y), std::min(b.lo.z, p.z)};
        b.hi = {std::max(b.hi.x, p.x), std::max(b.hi.y, p.y), std::max(b.hi.z, p.z)};
    }
    return b;
}

// Edge-incidence + per-vertex incident-face lists, in one O(F) pass.
void scanEdges(const MeshData& out, DefectReport& rep,
               std::unordered_map<std::uint64_t, int>& edgeCount,
               std::vector<std::vector<int>>& vFaces) {
    const int outF = static_cast<int>(out.indices.size()) / 3;
    for (int f = 0; f < outF; ++f) {
        const std::size_t b  = static_cast<std::size_t>(f) * 3;
        const int         a0 = out.indices[b], a1 = out.indices[b + 1], a2 = out.indices[b + 2];
        if (a0 == a1 || a1 == a2 || a0 == a2) {
            ++rep.degenFaces;
            continue;
        }
        ++edgeCount[edgeKey(a0, a1)];
        ++edgeCount[edgeKey(a1, a2)];
        ++edgeCount[edgeKey(a2, a0)];
        vFaces[static_cast<std::size_t>(a0)].push_back(f);
        vFaces[static_cast<std::size_t>(a1)].push_back(f);
        vFaces[static_cast<std::size_t>(a2)].push_back(f);
    }
    int shown = 0;
    for (const auto& [key, count] : edgeCount) {
        if (count == 1) {
            ++rep.boundaryEdges;
        } else if (count >= 3) {
            ++rep.nonManifoldEdges;
            const int v0 = static_cast<int>(key & 0xffffffffU);
            printDefect("nonmanif-edge", out.positions[static_cast<std::size_t>(v0)], rep.lo, rep.hi,
                        shown);
        }
    }
}

// A vertex is a pinch (bowtie) when its incident faces form more than one fan:
// group faces by the spoke edge (v,x) they share, union, count components.
bool isPinch(const MeshData& out, int v, const std::vector<int>& fs) {
    DSU                          dsu(static_cast<int>(fs.size()));
    std::unordered_map<int, int> spokeOwner;  // other-vertex x -> local face idx
    for (int li = 0; li < static_cast<int>(fs.size()); ++li) {
        const std::size_t b = static_cast<std::size_t>(fs[static_cast<std::size_t>(li)]) * 3;
        for (int k = 0; k < 3; ++k) {
            const int x = out.indices[b + static_cast<std::size_t>(k)];
            if (x == v) {
                continue;
            }
            const auto [it, inserted] = spokeOwner.try_emplace(x, li);
            if (!inserted) {
                dsu.join(li, it->second);
            }
        }
    }
    int roots = 0;
    for (int li = 0; li < static_cast<int>(fs.size()); ++li) {
        roots += (dsu.find(li) == li) ? 1 : 0;
    }
    return roots > 1;
}

void scanPinches(const MeshData& out, DefectReport& rep,
                 const std::vector<std::vector<int>>& vFaces) {
    int shown = 0;
    for (int v = 0; v < static_cast<int>(out.positions.size()); ++v) {
        const std::vector<int>& fs = vFaces[static_cast<std::size_t>(v)];
        if (fs.size() >= 2 && isPinch(out, v, fs)) {
            ++rep.pinchVerts;
            printDefect("pinch-vert", out.positions[static_cast<std::size_t>(v)], rep.lo, rep.hi,
                        shown);
        }
    }
}

// An output vertex whose nearest INPUT vertex is > 2% of the bbox diagonal away
// has been flung off the surface — a spike.
void scanSpikes(const MeshData& in, const MeshData& out, DefectReport& rep) {
    const Bounds  ib      = bboxOf(in.positions);
    const double  diag    = std::max(1e-9, static_cast<double>((ib.hi - ib.lo).length()));
    const double  thresh2 = (0.02 * diag) * (0.02 * diag);
    int           shown   = 0;
    for (const Vec3f& q : out.positions) {
        double best = std::numeric_limits<double>::max();
        for (const Vec3f& ip : in.positions) {
            best = std::min(best, static_cast<double>((q - ip).lengthSquared()));
            if (best <= thresh2) {
                break;
            }
        }
        if (best > thresh2) {
            ++rep.spikes;
            printDefect("spike", q, rep.lo, rep.hi, shown);
        }
    }
}

// Winding consistency: in a coherently-oriented mesh every directed edge (a->b)
// appears at most once — the adjacent face traverses it as (b->a). A directed
// edge seen twice means two neighbours are wound the same way (one is flipped),
// which shades dark / reads as a gash even though the topology is intact.
void scanWinding(const MeshData& out, DefectReport& rep) {
    const int                              outF = static_cast<int>(out.indices.size()) / 3;
    std::unordered_map<std::uint64_t, int> seen;  // packed (a<<32|b) -> face
    seen.reserve(out.indices.size());
    int shown = 0;
    for (int f = 0; f < outF; ++f) {
        const std::size_t          b = static_cast<std::size_t>(f) * 3;
        const std::array<int, 3>   v = {out.indices[b], out.indices[b + 1], out.indices[b + 2]};
        for (int k = 0; k < 3; ++k) {
            const int a0 = v[static_cast<std::size_t>(k)], a1 = v[static_cast<std::size_t>((k + 1) % 3)];
            if (a0 == a1) {
                continue;
            }
            const std::uint64_t dir = (static_cast<std::uint64_t>(a0) << 32) | static_cast<std::uint32_t>(a1);
            if (!seen.try_emplace(dir, f).second) {
                ++rep.flippedEdges;
                printDefect("flipped-wind", out.positions[static_cast<std::size_t>(a0)], rep.lo,
                            rep.hi, shown);
            }
        }
    }
}

Vec3f faceNormal(const MeshData& m, int f) {
    const std::size_t b = static_cast<std::size_t>(f) * 3;
    const Vec3f&      a = m.positions[static_cast<std::size_t>(m.indices[b])];
    const Vec3f&      p = m.positions[static_cast<std::size_t>(m.indices[b + 1])];
    const Vec3f&      q = m.positions[static_cast<std::size_t>(m.indices[b + 2])];
    return (p - a).cross(q - a);
}

// Fold-over: an interior edge whose two faces crease sharply back toward each
// other (normals dot < cos(threshold)). A near-180° crease is the dark "slit"
// the eye sees and the red gash in the normal view — geometry, not topology.
void scanFolds(const MeshData& out, DefectReport& rep,
               const std::vector<std::vector<int>>& vFaces) {
    std::unordered_map<std::uint64_t, int> edgeFace;  // undirected edge -> first face
    edgeFace.reserve(out.indices.size());
    int       shown    = 0;
    const int outF     = static_cast<int>(out.indices.size()) / 3;
    for (int f = 0; f < outF; ++f) {
        const std::size_t b    = static_cast<std::size_t>(f) * 3;
        const std::array<int, 3> v = {out.indices[b], out.indices[b + 1], out.indices[b + 2]};
        if (v[0] == v[1] || v[1] == v[2] || v[0] == v[2]) {
            continue;
        }
        for (int k = 0; k < 3; ++k) {
            const std::uint64_t e = edgeKey(v[static_cast<std::size_t>(k)],
                                            v[static_cast<std::size_t>((k + 1) % 3)]);
            const auto [it, inserted] = edgeFace.try_emplace(e, f);
            if (inserted) {
                continue;
            }
            const Vec3f n0 = faceNormal(out, it->second);
            const Vec3f n1 = faceNormal(out, f);
            const float l0 = n0.length();
            const float l1 = n1.length();
            if (l0 < 1e-12F || l1 < 1e-12F) {
                continue;
            }
            const float cosang = n0.dot(n1) / (l0 * l1);
            if (cosang < -0.5F) {  // > 120 degrees between normals = severe crease
                ++rep.foldEdges;
                printDefect("fold-crease", out.positions[static_cast<std::size_t>(v[static_cast<std::size_t>(k)])],
                            rep.lo, rep.hi, shown);
            }
        }
    }
    (void)vFaces;
}

// Sliver: a triangle whose smallest altitude is a tiny fraction of its longest
// edge (a needle). Thin shards shade as dark lines and pollute the normal field.
void scanSlivers(const MeshData& out, DefectReport& rep) {
    const int outF  = static_cast<int>(out.indices.size()) / 3;
    int       shown = 0;
    for (int f = 0; f < outF; ++f) {
        const std::size_t b = static_cast<std::size_t>(f) * 3;
        const Vec3f&      a = out.positions[static_cast<std::size_t>(out.indices[b])];
        const Vec3f&      p = out.positions[static_cast<std::size_t>(out.indices[b + 1])];
        const Vec3f&      q = out.positions[static_cast<std::size_t>(out.indices[b + 2])];
        const float       e0 = (p - a).length();
        const float       e1 = (q - p).length();
        const float       e2 = (a - q).length();
        const float       longest = std::max({e0, e1, e2});
        if (longest < 1e-12F) {
            continue;
        }
        const float area   = 0.5F * (p - a).cross(q - a).length();
        const float minAlt = 2.0F * area / longest;  // shortest altitude
        if (minAlt / longest < 0.02F) {              // aspect worse than 1:50
            ++rep.slivers;
            printDefect("sliver", a, rep.lo, rep.hi, shown);
        }
    }
}

DefectReport analyze(const MeshData& in, const MeshData& out) {
    DefectReport rep;
    const Bounds ob = bboxOf(out.positions);
    rep.lo          = ob.lo;
    rep.hi          = ob.hi;

    std::unordered_map<std::uint64_t, int> edgeCount;
    edgeCount.reserve(out.indices.size());
    std::vector<std::vector<int>> vFaces(out.positions.size());

    scanEdges(out, rep, edgeCount, vFaces);
    scanPinches(out, rep, vFaces);
    scanWinding(out, rep);
    scanFolds(out, rep, vFaces);
    scanSlivers(out, rep);
    scanSpikes(in, out, rep);
    return rep;
}

void run(const std::string& path, int target, double ratio) {
    MeshData in;
    try {
        in = decimate::load(path);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "skip %s: %s\n", path.c_str(), e.what());
        return;
    }
    const int inV = static_cast<int>(in.positions.size());
    const int inF = static_cast<int>(in.indices.size()) / 3;
    if (inF == 0) {
        std::fprintf(stderr, "skip %s: no faces\n", path.c_str());
        return;
    }
    const int tgt = target > 0 ? target : std::max(4, static_cast<int>(inV * ratio));

    const MeshData out = decimate::reduce(in, tgt);
    std::printf("\n=== %s  in=%dV/%dF  target=%d  out=%dV/%dF ===\n", path.c_str(), inV, inF, tgt,
                static_cast<int>(out.positions.size()), static_cast<int>(out.indices.size()) / 3);

    const DefectReport r = analyze(in, out);
    std::printf("  non-manifold edges : %d\n", r.nonManifoldEdges);
    std::printf("  pinch/bowtie verts : %d\n", r.pinchVerts);
    std::printf("  boundary edges     : %d   (original mesh boundary -> stays; growth = tearing)\n",
                r.boundaryEdges);
    std::printf("  flipped-winding edg: %d   (face wound opposite its neighbour -> dark shading)\n",
                r.flippedEdges);
    std::printf("  fold-crease edges  : %d   (>120 deg between adjacent faces -> dark slit/gash)\n",
                r.foldEdges);
    std::printf("  sliver triangles   : %d   (aspect worse than 1:50)\n", r.slivers);
    std::printf("  degenerate faces   : %d\n", r.degenFaces);
    std::printf("  spikes (>2%% diag)  : %d\n", r.spikes);
}

}  // namespace

int main(int argc, char** argv) {
    int    target = 0;
    double ratio  = 0.03;  // ~ the bunny guard floor, where defects showed
    std::vector<std::string> models;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--target") == 0 && i + 1 < argc) {
            target = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--ratio") == 0 && i + 1 < argc) {
            ratio = std::atof(argv[++i]);
        } else {
            models.emplace_back(argv[i]);
        }
    }
    if (models.empty()) {
        std::fprintf(stderr, "usage: mesh_defects [--target N | --ratio R] model...\n");
        return 1;
    }
    for (const std::string& m : models) {
        run(m, target, ratio);
    }
    return 0;
}
