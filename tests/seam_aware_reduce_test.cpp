// This file is part of the decimate library; see decimate.h for version/license details.
/// @file seam_aware_reduce_test.cpp
/// @brief Spec 17-03: seam-aware decimation (Liu et al. 2017, section 5).
///
/// Scenarios I1, I2, B1, S1-S8 (seam-aware decimation; I3, the
/// sweep subprocess, lives in CMake: `seamless_sweep_csv`).  Fixtures are planar ladders
/// whose seam geometry is known in closed form, so every bound below is an exact value.
///
/// The NOLINT on the TEST_CASEs is for Catch2: every CHECK macro expands to branching code
/// that readability-function-size / -cognitive-complexity count against the test body.

#include "model_fixture.h"

#include <algorithm>
#include <array>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normals.h>
#include <decimate/chart/chart.h>
#include <decimate/chart/parametrize.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/reduce/cluster_reduce.h>
#include <decimate/reduce/halfedge_reduce.h>
#include <decimate/reduce/reducer.h>
#include <decimate/reduce/seam_aware_reduce.h>
#include <decimate/seam/erase.h>
#include <decimate/seam/seam_energy.h>
#include <decimate/seam/seam_topology.h>
#include <decimate/seam/straighten.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/textured_mesh.h>
#include <functional>
#include <map>
#include <numbers>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

using decimate::MeshData;
using decimate::SeamAwareResult;
using decimate::TexturedMesh;
using decimate::math::Vec2f;
using decimate::math::Vec3f;

/// int -> size_t for indexing (indices are ints throughout the mesh types).
std::size_t idx(int i) {
    return static_cast<std::size_t>(i);
}


/// Bit-identical UVs (== would treat -0 and +0 as equal and hide a sign flip).
bool sameBits(Vec2f a, Vec2f b) {
    return std::bit_cast<std::uint32_t>(a.x) == std::bit_cast<std::uint32_t>(b.x) &&
           std::bit_cast<std::uint32_t>(a.y) == std::bit_cast<std::uint32_t>(b.y);
}


/// Builds a TexturedMesh from a triangle soup with per-corner UVs through the real `weld`.
TexturedMesh weldSoup(const std::vector<Vec3f> &pos, const std::vector<std::array<int, 3>> &tris,
                      const std::vector<std::array<Vec2f, 3>> &cornerUV) {
    MeshData mesh;
    mesh.positions = pos;
    std::vector<Vec2f> uv;
    for (std::size_t f = 0; f < tris.size(); ++f) {
        for (std::size_t k = 0; k < 3; ++k) {
            mesh.indices.push_back(tris[f][k]);
            uv.push_back(cornerUV[f][k]);
        }
    }
    return decimate::seam::weld(mesh, uv);
}


/// A ladder: two charts (left, right) glued along a vertical seam of `n` edges.  Outer
/// columns zigzag (in 3D and UV) so no 3D-boundary edge is collinear with its neighbour
/// and no interior edge (every vertex is on a UV silhouette) is collapsible: the seam is
/// the only thing that can collapse.  The seam copies default to a straight, uniformly
/// spaced vertical line; the callbacks override individual copies.
struct Ladder {
    int                       n = 4;
    float                     rightYScale = 1.0f;
    std::function<Vec2f(int)> leftSeam;  ///< UV of seam vertex j on the left chart.
    std::function<Vec2f(int)> rightSeam; ///< UV of seam vertex j on the right chart.
};


/// Builds the ladder described by @p o as a welded TexturedMesh.
TexturedMesh buildLadder(const Ladder &o) {
    // Three columns of n+1 vertices: left outer, seam (surface), right outer.
    const int          n = o.n;
    std::vector<Vec3f> pos;
    for (int j = 0; j <= n; ++j) {
        pos.emplace_back(-1.0f - (0.3f * static_cast<float>(j % 2)), static_cast<float>(j), 0.0f);
    }
    for (int j = 0; j <= n; ++j) {
        pos.emplace_back(0.0f, static_cast<float>(j), 0.0f);
    }
    for (int j = 0; j <= n; ++j) {
        pos.emplace_back(1.0f + (0.3f * static_cast<float>(j % 2)), static_cast<float>(j), 0.0f);
    }
    // Vertex ids per column.
    const auto left = [](int j) {
        return j;
    };
    const auto mid = [n](int j) {
        return n + 1 + j;
    };
    const auto right = [n](int j) {
        return (2 * n) + 2 + j;
    };

    // UVs: outer columns follow the 3D layout; seam copies default to a straight line.
    const auto uvLeft = [&](int j) {
        return Vec2f { 0.25f + (0.1f * pos[idx(j)].x), 0.05f + (0.1f * pos[idx(j)].y) };
    };
    const auto uvRight = [&](int j) {
        return Vec2f { 0.55f + (0.1f * pos[idx(right(j))].x),
                       0.05f + (0.1f * o.rightYScale * pos[idx(right(j))].y) };
    };
    const auto seamL = [&](int j) {
        return o.leftSeam ? o.leftSeam(j) : Vec2f { 0.25f, 0.05f + (0.1f * static_cast<float>(j)) };
    };
    const auto seamR = [&](int j) {
        return o.rightSeam ? o.rightSeam(j)
                           : Vec2f { 0.55f, 0.05f + (0.1f * o.rightYScale * static_cast<float>(j)) };
    };

    // Each quad strip becomes two triangles per chart; the seam copies get their own UVs.
    std::vector<std::array<int, 3>>   tris;
    std::vector<std::array<Vec2f, 3>> uv;
    for (int j = 0; j < n; ++j) {
        tris.push_back({ left(j), mid(j), mid(j + 1) });
        uv.push_back({ uvLeft(left(j)), seamL(j), seamL(j + 1) });
        tris.push_back({ left(j), mid(j + 1), left(j + 1) });
        uv.push_back({ uvLeft(left(j)), seamL(j + 1), uvLeft(left(j + 1)) });
        tris.push_back({ mid(j), right(j), right(j + 1) });
        uv.push_back({ seamR(j), uvRight(j), uvRight(j + 1) });
        tris.push_back({ mid(j), right(j + 1), mid(j + 1) });
        uv.push_back({ seamR(j), uvRight(j + 1), seamR(j + 1) });
    }
    return weldSoup(pos, tris, uv);
}


/// Straight seam, uniform spacing, right copy stretched 1.5x: ratio-matched.
TexturedMesh ladderStraight(int n = 4) {
    return buildLadder({ .n = n, .rightYScale = 1.5f });
}


/// Zigzag seam copies: no two adjacent seam edges are collinear anywhere.
TexturedMesh ladderBent() {
    const auto dx = [](int j) {
        return std::array<float, 5> { 0.0f, -0.02f, 0.02f, -0.02f, 0.0f }[idx(j)];
    };
    return buildLadder({ .n = 4,
                         .leftSeam =
                             [&](int j) {
                                 return Vec2f { 0.25f + dx(j), 0.05f + (0.1f * static_cast<float>(j)) };
                             },
                         .rightSeam = [&](int j) {
                             return Vec2f { 0.55f + dx(j), 0.05f + (0.1f * static_cast<float>(j)) };
                         } });
}


/// Straight but ratio-mismatched: the right copy's vertices are spaced non-uniformly.
TexturedMesh ladderMismatched() {
    const auto ry = [](int j) {
        return std::array<float, 5> { 0.0f, 0.1f, 0.3f, 0.4f, 0.8f }[idx(j)];
    };
    return buildLadder(
        { .n = 4, .rightSeam = [&](int j) {
             return Vec2f { 0.55f, 0.05f + ry(j) };
         } });
}


/// d,e unifiable (straight, uniform), f bent: chain v0..v3, f = (v2, v3).
TexturedMesh ladderOneNeighbour() {
    const auto dx = [](int j) {
        return j == 3 ? 0.04f : 0.0f;
    };
    return buildLadder({ .n = 3,
                         .leftSeam =
                             [&](int j) {
                                 return Vec2f { 0.25f + dx(j), 0.05f + (0.1f * static_cast<float>(j)) };
                             },
                         .rightSeam = [&](int j) {
                             return Vec2f { 0.55f + dx(j), 0.05f + (0.1f * static_cast<float>(j)) };
                         } });
}


/// A convex hexagon fan: one interior vertex, six boundary vertices, one chart.
TexturedMesh hexFan() {
    std::vector<Vec3f> pos { { 0.0f, 0.0f, 0.0f } };
    for (int k = 0; k < 6; ++k) {
        const float a = static_cast<float>(k) * std::numbers::pi_v<float> / 3.0f;
        pos.emplace_back(std::cos(a), std::sin(a), 0.0f);
    }
    std::vector<std::array<int, 3>>   tris;
    std::vector<std::array<Vec2f, 3>> uv;
    const auto                        uvOf = [&](int v) {
        return Vec2f { 0.5f + (0.2f * pos[idx(v)].x), 0.5f + (0.2f * pos[idx(v)].y) };
    };
    for (int k = 0; k < 6; ++k) {
        const int a = 1 + k;
        const int b = 1 + ((k + 1) % 6);
        tris.push_back({ 0, a, b });
        uv.push_back({ uvOf(0), uvOf(a), uvOf(b) });
    }
    return weldSoup(pos, tris, uv);
}


/// A 2x1 ribbon whose every vertex is on the boundary; the boundary zigzags so no
/// boundary edge is collinear with its neighbour.  Every interior edge joins two
/// silhouette vertices.
TexturedMesh ribbon() {
    const std::vector<Vec3f>              pos { { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.3f, 0.0f }, { 2.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, { 1.0f, 0.7f, 0.0f }, { 2.0f, 1.0f, 0.0f } };
    const std::vector<std::array<int, 3>> tris { { 0, 1, 4 }, { 0, 4, 3 }, { 1, 2, 5 }, { 1, 5, 4 } };
    std::vector<std::array<Vec2f, 3>>     uv;
    for (const auto &t : tris) {
        std::array<Vec2f, 3> c;
        for (std::size_t k = 0; k < 3; ++k) {
            c[k] = { 0.2f + (0.2f * pos[idx(t[k])].x), 0.2f + (0.2f * pos[idx(t[k])].y) };
        }
        uv.push_back(c);
    }
    return weldSoup(pos, tris, uv);
}


/// An n x n vertex grid mapped to the unit disc (square -> disc), one chart: the
/// boundary is a closed polygon on the unit circle.  With @p square the grid stays a
/// square (collinear boundary runs).
TexturedMesh gridPatch(int n, bool square) {
    // Square -> disc mapping keeps the grid regular but makes the boundary a circle.
    std::vector<Vec3f> pos;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            const double x = (2.0 * i / (n - 1)) - 1.0;
            const double y = (2.0 * j / (n - 1)) - 1.0;
            const double px = square ? x : x * std::sqrt(1.0 - (0.5 * y * y));
            const double py = square ? y : y * std::sqrt(1.0 - (0.5 * x * x));
            pos.emplace_back(static_cast<float>(px), static_cast<float>(py), 0.0f);
        }
    }
    std::vector<std::array<int, 3>>   tris;
    std::vector<std::array<Vec2f, 3>> uv;
    const auto                        uvOf = [&](int v) {
        return Vec2f { 0.5f + (0.4f * pos[idx(v)].x), 0.5f + (0.4f * pos[idx(v)].y) };
    };
    for (int j = 0; j + 1 < n; ++j) {
        for (int i = 0; i + 1 < n; ++i) {
            const int a = (j * n) + i;
            const int b = a + 1;
            const int c = a + n + 1;
            const int d = a + n;
            tris.push_back({ a, b, c });
            uv.push_back({ uvOf(a), uvOf(b), uvOf(c) });
            tris.push_back({ a, c, d });
            uv.push_back({ uvOf(a), uvOf(c), uvOf(d) });
        }
    }
    return weldSoup(pos, tris, uv);
}


/// Link of every vertex: v -> (a -> b) for each face (v, a, b).  Reports directed edges used
/// twice and vertices two faces leave through one neighbour through @p note.
using LinkMap = std::map<int, std::map<int, int>>;

LinkMap vertexLinks(const TexturedMesh &tm, const std::function<void(const std::string &)> &note) {
    std::set<std::pair<int, int>> directed;
    LinkMap                       link;
    for (std::size_t c = 0; c < tm.faceV.size(); ++c) {
        const std::size_t f = c / 3;
        const std::size_t k = c % 3;
        const int         v = tm.faceV[c];
        const int         a = tm.faceV[(3 * f) + ((k + 1) % 3)];
        const int         b = tm.faceV[(3 * f) + ((k + 2) % 3)];
        if (!directed.emplace(v, a).second) {
            note("edge " + std::to_string(v) + "->" + std::to_string(a) +
                 " has two faces on one side");
        }
        if (!link[v].emplace(a, b).second) {
            note("vertex " + std::to_string(v) + ": two faces leave through neighbour " +
                 std::to_string(a));
        }
    }
    return link;
}


/// Is the successor map @p succ (a link) one path or one cycle?
bool isOneFan(const std::map<int, int> &succ) {
    std::set<int> heads;
    for (const auto &[a, b] : succ) {
        heads.insert(b);
    }
    int start = succ.begin()->first;
    for (const auto &[a, b] : succ) {
        start = heads.count(a) == 0 ? a : start;
    }
    std::size_t walked = 1;
    for (int cur = succ.at(start); succ.count(cur) != 0 && cur != start && walked <= succ.size();
         cur = succ.at(cur)) {
        ++walked;
    }
    return walked == succ.size();
}


/// Number of manifoldness defects: directed edges used by two faces, vertices two faces
/// leave through one neighbour, vertices whose link is not one fan.  0 = a clean 2-manifold
/// (with boundary).  @p first receives a description of the first defect.
int manifoldDefects(const TexturedMesh &tm, std::string *first = nullptr) {
    int        defects = 0;
    const auto note = [&](const std::string &what) {
        if (defects++ == 0 && first != nullptr) {
            *first = what;
        }
    };
    for (const auto &[v, succ] : vertexLinks(tm, note)) {
        if (!isOneFan(succ)) {
            note("vertex " + std::to_string(v) + ": link is not one fan");
        }
    }
    return defects;
}


/// True if the mesh has no manifoldness defect (see manifoldDefects).
bool isManifold(const TexturedMesh &tm) {
    return manifoldDefects(tm) == 0;
}


/// Twice the 3D area of face @p f.
double surfaceArea2(const TexturedMesh &tm, std::size_t f) {
    const Vec3f a = tm.positions[idx(tm.faceV[3 * f])];
    const Vec3f b = tm.positions[idx(tm.faceV[(3 * f) + 1])];
    const Vec3f c = tm.positions[idx(tm.faceV[(3 * f) + 2])];
    return (b - a).cross(c - a).length();
}


/// Twice the signed UV area of face @p f (sign = orientation).
double uvArea2(const TexturedMesh &tm, std::size_t f) {
    const Vec2f a = tm.uvs[idx(tm.faceU[3 * f])];
    const Vec2f b = tm.uvs[idx(tm.faceU[(3 * f) + 1])];
    const Vec2f c = tm.uvs[idx(tm.faceU[(3 * f) + 2])];
    return (static_cast<double>(b.x - a.x) * static_cast<double>(c.y - a.y)) -
           (static_cast<double>(b.y - a.y) * static_cast<double>(c.x - a.x));
}


/// True if no face has (near-)zero 3D area.
bool noZeroAreaFace(const TexturedMesh &tm) {
    for (std::size_t f = 0; f < tm.faceCount(); ++f) {
        if (surfaceArea2(tm, f) <= 1e-12) {
            return false;
        }
    }
    return true;
}


/// Number of UV triangles with (near-)zero area.
int degenerateUvFaces(const TexturedMesh &tm) {
    int bad = 0;
    for (std::size_t f = 0; f < tm.faceCount(); ++f) {
        bad += std::abs(uvArea2(tm, f)) <= 1e-14 ? 1 : 0;
    }
    return bad;
}


/// Face -> connected component of the UV mesh (faces sharing a UV vertex).
std::vector<int> uvComponents(const TexturedMesh &tm, int &count) {
    std::vector<int> parent(tm.uvs.size());
    for (std::size_t i = 0; i < parent.size(); ++i) {
        parent[i] = static_cast<int>(i);
    }
    const std::function<int(int)> find = [&](int x) {
        while (parent[idx(x)] != x) {
            parent[idx(x)] = parent[idx(parent[idx(x)])];
            x = parent[idx(x)];
        }
        return x;
    };
    for (std::size_t f = 0; f < tm.faceCount(); ++f) {
        parent[idx(find(tm.faceU[3 * f]))] = find(tm.faceU[(3 * f) + 1]);
        parent[idx(find(tm.faceU[(3 * f) + 1]))] = find(tm.faceU[(3 * f) + 2]);
    }
    std::map<int, int> label;
    std::vector<int>   faceComp(tm.faceCount());
    for (std::size_t f = 0; f < tm.faceCount(); ++f) {
        const int root = find(tm.faceU[3 * f]);
        faceComp[f] = label.emplace(root, static_cast<int>(label.size())).first->second;
    }
    count = static_cast<int>(label.size());
    return faceComp;
}


/// UV triangles whose orientation disagrees with the majority of their UV component: a
/// fold-over inside a chart.  (Charts may be mirrored as a whole; only a mix flags.)
int uvOrientationConflicts(const TexturedMesh &tm) {
    int                             count = 0;
    const std::vector<int>          comp = uvComponents(tm, count);
    std::vector<std::array<int, 2>> signs(idx(count), { 0, 0 });
    for (std::size_t f = 0; f < tm.faceCount(); ++f) {
        const double a = uvArea2(tm, f);
        signs[idx(comp[f])][a > 0.0 ? 1 : 0] += a == 0.0 ? 0 : 1;
    }
    int conflicts = 0;
    for (const auto &[neg, pos] : signs) {
        conflicts += std::min(neg, pos);
    }
    return conflicts;
}


/// Is @p p inside triangle (a,b,c) of either winding, within distance @p tol of its edges?
bool insideTriangle(Vec2f p, Vec2f a, Vec2f b, Vec2f c, double tol) {
    const auto edge = [&](Vec2f s, Vec2f e) {
        const double ex = static_cast<double>(e.x) - s.x;
        const double ey = static_cast<double>(e.y) - s.y;
        const double cross = (ex * (static_cast<double>(p.y) - s.y)) -
                             (ey * (static_cast<double>(p.x) - s.x));
        return cross / std::hypot(ex, ey); // signed distance, left positive
    };
    const double area = ((static_cast<double>(b.x) - a.x) * (static_cast<double>(c.y) - a.y)) -
                        ((static_cast<double>(b.y) - a.y) * (static_cast<double>(c.x) - a.x));
    const double s = area >= 0.0 ? 1.0 : -1.0;
    return s * edge(a, b) >= -tol && s * edge(b, c) >= -tol && s * edge(c, a) >= -tol;
}


/// A UV mesh with its connected components (faces sharing a UV vertex).
struct Components {
    const TexturedMesh &mesh;
    std::vector<int>    faceComp;
    int                 count = 0;
};


/// Labels the UV components of @p tm.
Components componentsOf(const TexturedMesh &tm) {
    Components c { tm, {}, 0 };
    c.faceComp = uvComponents(tm, c.count);
    return c;
}


/// Is @p p inside a triangle of component @p comp of @p in (tolerance @p tol)?
bool pointInComponent(const Components &in, int comp, Vec2f p, double tol) {
    const TexturedMesh &m = in.mesh;
    for (std::size_t g = 0; g < m.faceCount(); ++g) {
        if (in.faceComp[g] == comp &&
            insideTriangle(p, m.uvs[idx(m.faceU[3 * g])], m.uvs[idx(m.faceU[(3 * g) + 1])],
                           m.uvs[idx(m.faceU[(3 * g) + 2])], tol)) {
            return true;
        }
    }
    return false;
}


/// Is every UV vertex of component @p comp of @p lod inside the union of the triangles of
/// the single input component @p inComp?
bool componentInside(const Components &in, int inComp, const Components &lod, int comp,
                     double tol) {
    for (std::size_t c = 0; c < lod.mesh.faceU.size(); ++c) {
        const Vec2f p = lod.mesh.uvs[idx(lod.mesh.faceU[c])];
        if (lod.faceComp[c / 3] == comp && !pointInComponent(in, inComp, p, tol)) {
            return false;
        }
    }
    return true;
}


/// Number of LOD UV components NOT contained in the union of one input component.
int componentsOutside(const TexturedMesh &in, const TexturedMesh &lod, double tol) {
    const Components inC = componentsOf(in);
    const Components lodC = componentsOf(lod);
    int              outside = 0;
    for (int c = 0; c < lodC.count; ++c) {
        bool ok = false;
        for (int i = 0; i < inC.count && !ok; ++i) {
            ok = componentInside(inC, i, lodC, c, tol);
        }
        outside += ok ? 0 : 1;
    }
    return outside;
}


/// Seam vertices of a ladder output, each with the UV of its left copy (u < 0.4) and its
/// right copy, ordered by surface y.
struct ChainPoint {
    float y;
    Vec2f left;
    Vec2f right;
};


/// Records the UV of both endpoints of seam half-edge @p side, by surface vertex and copy.
void collectCopies(const TexturedMesh &tm, const decimate::seam::HalfEdgeRef &side,
                   std::map<int, std::pair<Vec2f, Vec2f>> &byVertex) {
    for (std::size_t k = 0; k < 2; ++k) {
        const std::size_t corner = (3 * idx(side.face)) + ((idx(side.corner) + k) % 3);
        const int         sv = tm.faceV[corner];
        const Vec2f       uv = tm.uvs[idx(tm.faceU[corner])];
        (uv.x < 0.4f ? byVertex[sv].first : byVertex[sv].second) = uv;
    }
}


/// The seam chain of a ladder, ordered by surface y.
std::vector<ChainPoint> chainOf(const TexturedMesh &tm) {
    std::map<int, std::pair<Vec2f, Vec2f>> byVertex;
    for (const auto &pair : decimate::seam::findSeams(tm).seams) {
        collectCopies(tm, pair.a, byVertex);
        collectCopies(tm, pair.b, byVertex);
    }
    // Sort by surface y so chains from different reductions compare index by index.
    std::vector<ChainPoint> out;
    for (const auto &[sv, copies] : byVertex) {
        out.push_back({ tm.positions[idx(sv)].y, copies.first, copies.second });
    }
    std::ranges::sort(out, {}, &ChainPoint::y);
    return out;
}


/// FNV-1a style fingerprint of vertex count and index list (positions are not hashed).
std::uint64_t topologyHash(const MeshData &m) {
    std::uint64_t h = 1469598103934665603ULL;
    h = (h ^ m.positions.size()) * 1099511628211ULL;
    for (const int i : m.indices) {
        h = (h ^ static_cast<std::uint32_t>(i)) * 1099511628211ULL;
    }
    return h;
}


/// Boundary edges of the UV mesh, as UV-vertex pairs (a half-edge whose reverse is absent).
std::vector<std::pair<int, int>> uvBoundary(const TexturedMesh &tm) {
    std::set<std::pair<int, int>> all;
    for (std::size_t f = 0; f < tm.faceCount(); ++f) {
        for (std::size_t k = 0; k < 3; ++k) {
            all.emplace(tm.faceU[(3 * f) + k], tm.faceU[(3 * f) + ((k + 1) % 3)]);
        }
    }
    std::vector<std::pair<int, int>> out;
    for (const auto &[a, b] : all) {
        if (all.count({ b, a }) == 0) {
            out.emplace_back(a, b);
        }
    }
    return out;
}


/// Distance from @p p to segment (a,b).
double distToSegment(Vec2f p, Vec2f a, Vec2f b) {
    const double ex = static_cast<double>(b.x) - a.x;
    const double ey = static_cast<double>(b.y) - a.y;
    const double len2 = (ex * ex) + (ey * ey);
    double       t = len2 > 0.0 ? (((static_cast<double>(p.x) - a.x) * ex) +
                             ((static_cast<double>(p.y) - a.y) * ey)) /
                                len2
                                : 0.0;
    t = std::clamp(t, 0.0, 1.0);
    return std::hypot(static_cast<double>(p.x) - (a.x + (t * ex)),
                      static_cast<double>(p.y) - (a.y + (t * ey)));
}


/// Is the boundary one closed loop (every vertex has one outgoing and one incoming
/// boundary edge, and walking from any vertex visits them all)?
bool oneClosedLoop(const std::vector<std::pair<int, int>> &edges) {
    std::map<int, int> next;
    for (const auto &[a, b] : edges) {
        if (!next.emplace(a, b).second) {
            return false;
        }
    }
    if (next.empty()) {
        return false;
    }
    const int   start = next.begin()->first;
    std::size_t steps = 0;
    int         cur = start;
    do {
        const auto it = next.find(cur);
        if (it == next.end()) {
            return false;
        }
        cur = it->second;
        ++steps;
    } while (cur != start && steps <= next.size());
    return steps == next.size();
}


/// B1 body: reducing to @p target keeps the UV boundary one closed loop on the input's boundary.
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
void checkSilhouettePreserved(const TexturedMesh &tm, int target) {
    const SeamAwareResult res = decimate::reduceSeamAware(tm, target);
    REQUIRE(res.mesh.positions.size() < tm.positions.size());
    const auto inBoundary = uvBoundary(tm);
    const auto outBoundary = uvBoundary(res.mesh);
    CHECK(oneClosedLoop(inBoundary));
    CHECK(oneClosedLoop(outBoundary));
    CHECK(res.stats.seamEdgesBefore == 0);
    CHECK(res.stats.collapsedSeamEdges == 0);
    for (const auto &[a, b] : outBoundary) {
        for (const int q : { a, b }) {
            const Vec2f p = res.mesh.uvs[idx(q)];
            double      best = 1e9;
            for (const auto &[ia, ib] : inBoundary) {
                best = std::min(best, distToSegment(p, tm.uvs[idx(ia)], tm.uvs[idx(ib)]));
            }
            CHECK(best <= 1e-6);
        }
    }
    CHECK(isManifold(res.mesh));
    CHECK(noZeroAreaFace(res.mesh));
}


using decimate::bake::NormalMap;

constexpr int kTex = 256;
constexpr int kLowVertices = 1500;

/// 17-02's pipeline: reduce -> segment -> parametrize -> weld -> straighten -> bake -> erase.
struct Full {
    TexturedMesh                tm;
    decimate::seam::EraseResult erased;
};


/// Builds the full atlas of @p path (optionally straightened) with its erased-seam statistics.
Full buildFull(const std::string &path, bool straighten) {
    const MeshData              high = decimate::load(path);
    const MeshData              low = decimate::reduce(high, kLowVertices);
    std::vector<decimate::Face> faces;
    for (std::size_t i = 0; i + 2 < low.indices.size(); i += 3) {
        faces.push_back({ low.indices[i], low.indices[i + 1], low.indices[i + 2] });
    }
    const auto hem = decimate::HalfEdgeMesh::create(static_cast<int>(low.positions.size()),
                                                    static_cast<int>(faces.size()), faces.data());
    auto       charts = decimate::chart::segment(hem, low, 64);
    const auto cornerUV = decimate::chart::parametrizeCorners(hem, low, charts, kTex);
    Full       full;
    full.tm = decimate::seam::weld(low, cornerUV);
    if (straighten) {
        full.tm = decimate::seam::straighten(full.tm, { .tau = 0.01 }).mesh;
    }
    const auto lowN = decimate::bake::computeVertexNormals(low, decimate::bake::Weighting::EqualUnit);
    const auto highN =
        decimate::bake::computeVertexNormals(high, decimate::bake::Weighting::EqualUnit);
    const NormalMap map =
        decimate::bake::bake(low, full.tm.cornerUV(), lowN, high, highN, kTex, 2, false);
    full.erased = decimate::seam::erase(full.tm, map, decimate::seam::NormalSpace::Object);
    return full;
}


/// A reduced LOD with its seam energy per channel.
struct LodReport {
    SeamAwareResult       res;
    std::array<double, 3> dTotal {};
};


/// Reduces @p full to @p target vertices and measures the seam energy of the result.
LodReport lodOf(const Full &full, int target) {
    LodReport r;
    r.res = decimate::reduceSeamAware(full.tm, target);
    r.dTotal = decimate::seam::seamEnergy(r.res.mesh, kTex, kTex, full.erased.texels);
    return r;
}


/// Asserts the I1/I2 invariants of one LOD against the input atlas.
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity,bugprone-easily-swappable-parameters)
void checkLod(const Full &full, const LodReport &lod, const char *name, int percent,
              int inputDefects, std::size_t achievedVertices) {
    const TexturedMesh &m = lod.res.mesh;
    INFO(name << " LOD target " << percent << "%");
    // The achieved count, exactly.  The reducer is deterministic (equal input, bit-equal
    // output) and stops above the target when no legal collapse remains, so these numbers
    // are the measured result of THIS atlas, not the target: a regression that removes fewer
    // vertices fails here.  They legitimately move only if the atlas pipeline or the
    // collapse criteria change -- then re-measure and update, do not loosen.
    CHECK(m.positions.size() == achievedVertices);
    // The reducer must not create a defect.  (The open bunny's 1500-vertex input already
    // carries a pinched vertex from reduce(); the spot input is clean, so its LODs must be.)
    CHECK(manifoldDefects(m) <= inputDefects);
    CHECK(noZeroAreaFace(m));
    CHECK(degenerateUvFaces(m) <= degenerateUvFaces(full.tm));
    CHECK(uvOrientationConflicts(m) <= uvOrientationConflicts(full.tm));
    for (std::size_t c = 0; c < 3; ++c) {
        CHECK(lod.dTotal[c] <= 1e-6 * full.erased.dBefore[c]);
    }
    CHECK(componentsOutside(full.tm, m, 1e-6) == 0);
}


/// Prints the LOD statistics (diagnostic only, never asserted).
void printLod(const char *name, const char *atlas, int percent, const Full &full,
              const LodReport &lod) {
    const auto &st = lod.res.stats;
    std::printf("[reduceSeamAware] %s %s %d%%: faces %zu->%zu verts %zu->%zu | D_total "
                "%.3e %.3e %.3e (D0 %.3e %.3e %.3e) | seams %d->%d collapsible %d "
                "uncollapsible %d collapsed %d\n",
                name, atlas, percent, full.tm.faceCount(), lod.res.mesh.faceCount(),
                full.tm.positions.size(), lod.res.mesh.positions.size(), lod.dTotal[0],
                lod.dTotal[1], lod.dTotal[2], full.erased.dBefore[0], full.erased.dBefore[1],
                full.erased.dBefore[2], st.seamEdgesBefore, st.seamEdgesAfter,
                st.collapsibleSeamEdges, st.uncollapsibleSeamEdges, st.collapsedSeamEdges);
}


void printPlain(const std::string &path, const char *name);

/// I1/I2 on one model: the straightened atlas is asserted; the unstraightened one is
/// measured and printed for comparison.
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
void checkModel(const std::string &path, const char *name, bool expectSeamCollapses,
                std::array<std::size_t, 2> achievedVertices) {
    const Full  full = buildFull(path, true);
    std::string firstDefect;
    const int   inputDefects = manifoldDefects(full.tm, &firstDefect);
    std::printf("[input] %s manifold defects %d (%s) zeroArea %d degenerateUV %d\n", name,
                inputDefects, firstDefect.c_str(), noZeroAreaFace(full.tm) ? 1 : 0,
                degenerateUvFaces(full.tm));
    const int nV = static_cast<int>(full.tm.positions.size());
    for (const int percent : { 50, 10 }) {
        const LodReport lod = lodOf(full, percent == 50 ? nV / 2 : nV / 10);
        checkLod(full, lod, name, percent, inputDefects, achievedVertices[percent == 50 ? 0 : 1]);
        printLod(name, "straightened", percent, full, lod);
        if (expectSeamCollapses) {
            CHECK(lod.res.stats.collapsedSeamEdges > 0);
        }
    }
    // The unstraightened comparison costs a second bake + erase: always on the small spot,
    // on request (DECIMATE_LOD_PLAIN=1) elsewhere.  Measured and printed, not asserted.
    if (!expectSeamCollapses || std::getenv("DECIMATE_LOD_PLAIN") != nullptr) {
        printPlain(path, name);
    }
}


void printPlain(const std::string &path, const char *name) {
    const Full plain = buildFull(path, false);
    const int  nP = static_cast<int>(plain.tm.positions.size());
    for (const int percent : { 50, 10 }) {
        const LodReport lod = lodOf(plain, percent == 50 ? nP / 2 : nP / 10);
        printLod(name, "plain", percent, plain, lod);
    }
}

} // namespace


/// S1, S6: classification of seam edges
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S1: bent seam copies never collapse", "[seam][reduce][S1]") {
    const TexturedMesh tm = ladderBent();
    REQUIRE(tm.positions.size() == 15);
    const SeamAwareResult res = decimate::reduceSeamAware(tm, 8);
    CHECK(res.stats.collapsedSeamEdges == 0);
    CHECK(res.stats.uncollapsibleSeamEdges == 4);
    CHECK(res.stats.collapsibleSeamEdges == 0);
    CHECK(res.stats.seamEdgesBefore == 4);
    CHECK(res.stats.seamEdgesAfter == 4);
    CHECK(chainOf(res.mesh).size() == 5);
}


// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S6: classifySeamEdges on the fixtures", "[seam][reduce][S6]") {
    const auto bent = decimate::classifySeamEdges(ladderBent());
    CHECK(bent.collapsible == 0);
    CHECK(bent.uncollapsible == 4);

    const auto straight = decimate::classifySeamEdges(ladderStraight());
    CHECK(straight.collapsible == 4);
    CHECK(straight.uncollapsible == 0);

    // Straight, but the right copy's edge-length ratios differ: collinearity alone is not enough.
    const auto mismatched = decimate::classifySeamEdges(ladderMismatched());
    CHECK(mismatched.collapsible == 0);
    CHECK(mismatched.uncollapsible == 4);

    // d,e unifiable, f not: d and e collapse, f does not.
    const auto one = decimate::classifySeamEdges(ladderOneNeighbour());
    CHECK(one.collapsible == 2);
    CHECK(one.uncollapsible == 1);

    // The UV-defined seam count equals findSeams' on a fixture.
    const TexturedMesh tm = ladderStraight();
    CHECK(static_cast<std::size_t>(straight.collapsible + straight.uncollapsible) ==
          decimate::seam::findSeams(tm).seams.size());
}


/// S2, S3: unifiable neighbours
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S2: one unifiable neighbour -- the new UV pair is the shared endpoint", "[seam][reduce][S2]") {
    const TexturedMesh tm = ladderOneNeighbour();
    const auto         before = chainOf(tm);
    REQUIRE(before.size() == 4);
    const SeamAwareResult res = decimate::reduceSeamAware(tm, static_cast<int>(tm.positions.size()) - 1);
    REQUIRE(res.stats.collapsedSeamEdges == 1);
    const auto after = chainOf(res.mesh);
    REQUIRE(after.size() == 3);
    // v1 is gone; v0, v2, v3 keep their input UVs bit-for-bit on both copies, so the
    // merged vertex sits exactly on the shared endpoint e1(1) / e2(1) (or d1(0) / d2(0)).
    const std::array<std::size_t, 3> kept { 0, 2, 3 };
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(sameBits(after[i].left, before[kept[i]].left));
        CHECK(sameBits(after[i].right, before[kept[i]].right));
    }
}


// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S3: two unifiable neighbours -- slide on the line, same parameter", "[seam][reduce][S3]") {
    const TexturedMesh    tm = ladderStraight();
    const SeamAwareResult res = decimate::reduceSeamAware(tm, 6);
    CHECK(res.stats.collapsedSeamEdges >= 1);
    CHECK(res.stats.collapsibleSeamEdges == 4);
    const auto chain = chainOf(res.mesh);
    REQUIRE(chain.size() >= 2);
    REQUIRE(chain.size() < 5);
    for (const auto &p : chain) {
        CHECK(std::abs(p.left.x - 0.25f) <= 1e-6f); // collinear: the input line
        CHECK(std::abs(p.right.x - 0.55f) <= 1e-6f);
    }
    // Adjacent edge-length ratios match across copies (the right copy is 1.5x longer).
    for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
        const double lenL = std::hypot(chain[i + 1].left.x - chain[i].left.x,
                                       chain[i + 1].left.y - chain[i].left.y);
        const double lenR = std::hypot(chain[i + 1].right.x - chain[i].right.x,
                                       chain[i + 1].right.y - chain[i].right.y);
        CHECK(std::abs((lenR / lenL) - 1.5) <= 1e-5 * 1.5);
    }
}


/// S4, S5: silhouette rules on interior edges
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S4: one silhouette endpoint -- the new wedge takes that vertex's UV", "[seam][reduce][S4]") {
    const TexturedMesh tm = hexFan();
    REQUIRE(tm.positions.size() == 7);
    const SeamAwareResult res = decimate::reduceSeamAware(tm, 6);
    REQUIRE(res.mesh.positions.size() == 6);
    REQUIRE(res.mesh.uvs.size() == 6);
    // The interior vertex (UV vertex 0) merged into a rim vertex: every surviving UV is
    // one of the input rim UVs, bit-for-bit.
    for (const Vec2f &uv : res.mesh.uvs) {
        bool found = false;
        for (std::size_t q = 1; q < tm.uvs.size(); ++q) {
            found = found || sameBits(uv, tm.uvs[q]);
        }
        CHECK(found);
    }
    CHECK(res.mesh.faceCount() == 4);
    CHECK(isManifold(res.mesh));
}


// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S5: an interior edge between two silhouette vertices never collapses", "[seam][reduce][S5]") {
    const TexturedMesh    tm = ribbon();
    const SeamAwareResult res = decimate::reduceSeamAware(tm, 3);
    CHECK(res.mesh.positions.size() == 6);
    CHECK(res.mesh.faceCount() == 4);
    CHECK(res.mesh.faceV == tm.faceV);
    CHECK(res.mesh.faceU == tm.faceU);
    for (std::size_t v = 0; v < tm.positions.size(); ++v) {
        CHECK(res.mesh.positions[v] == tm.positions[v]);
    }
}


/// S7, S8: contracts
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S7: target >= vertex count returns the input untouched", "[seam][reduce][S7]") {
    const TexturedMesh tm = ladderStraight();
    const int          nV = static_cast<int>(tm.positions.size());
    for (const int target : { nV, nV + 10 }) {
        const SeamAwareResult res = decimate::reduceSeamAware(tm, target);
        CHECK(res.mesh.positions == tm.positions);
        CHECK(res.mesh.faceV == tm.faceV);
        CHECK(res.mesh.faceU == tm.faceU);
        REQUIRE(res.mesh.uvs.size() == tm.uvs.size());
        for (std::size_t q = 0; q < tm.uvs.size(); ++q) {
            CHECK(sameBits(res.mesh.uvs[q], tm.uvs[q]));
        }
        CHECK(res.stats.seamEdgesBefore == 4);
        CHECK(res.stats.seamEdgesAfter == 4);
        CHECK(res.stats.collapsibleSeamEdges == 0);
        CHECK(res.stats.uncollapsibleSeamEdges == 0);
        CHECK(res.stats.collapsedSeamEdges == 0);
    }
}


// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S8: the untextured reducer family is unchanged", "[seam][reduce][S8]") {
    const MeshData spot = decimate::loadObj(std::string(MESH_MODEL_DIR) + "/spot.obj");
    // Topology fingerprints captured from HEAD before reduceSeamAware existed.
    CHECK(topologyHash(decimate::reduce(spot, 1500)) == 6197302785880904982ULL);
    CHECK(topologyHash(decimate::reduceHalfEdge(spot, 1500)) == 2018509701008911115ULL);
    CHECK(topologyHash(decimate::reduceClustering(spot, 1500)) == 4749293316192915343ULL);
}


/// B1: the silhouette shape is preserved
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("B1: a single-chart disc keeps its UV silhouette", "[seam][reduce][B1]") {
    const TexturedMesh disc = gridPatch(8, false);
    REQUIRE(disc.positions.size() == 64);
    const SeamAwareResult res = decimate::reduceSeamAware(disc, 20);
    std::printf("[reduceSeamAware] disc 64 -> %zu vertices, %zu faces\n", res.mesh.positions.size(),
                res.mesh.faceCount());
    checkSilhouettePreserved(disc, 20);

    // The collinear-boundary variant: a square grid slides along its straight sides.
    const TexturedMesh    square = gridPatch(8, true);
    const SeamAwareResult sq = decimate::reduceSeamAware(square, 20);
    std::printf("[reduceSeamAware] square 64 -> %zu vertices, %zu faces\n", sq.mesh.positions.size(),
                sq.mesh.faceCount());
    CHECK(sq.mesh.positions.size() <= 20);
    checkSilhouettePreserved(square, 20);
}


/// I1, I2: the real atlases, original texture on LOD UVs
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("I1: seam-aware LODs of the closed spot atlas keep the erased seam erased", "[seam][reduce][I1][release-pinned]") {
    // Pinned vertex counts, 50% and 10% of the input.  Re-pin history of the 10% count:
    // pre-19: 358; bijective repair (chart 4 re-flattened, not split): 19-01 382, 19-03 391;
    // Spec 19 A3 (bijective flattener is the default): 381; Spec 21 (spike-shift footprint
    // changes packing): 377; Spec 22-02 (barrier candidates summed in ascending vertex id, a
    // one-time float-order change of the bijective flattener; spot_uv_golden unchanged at
    // 5e-10): 379.  TrueCentroid, build-rel.
    // The count is float-order sensitive across build types (19-01 saw 382 in build/, 381 in
    // build-rel); 50% stays 750.
    checkModel(std::string(MESH_MODEL_DIR) + "/spot.obj", "spot", false, { 750, 379 });
}


/// I2: as I1, on the open bunny.
// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("I2: seam-aware LODs of the open bunny atlas", "[seam][reduce][I2][release-pinned]") {
    // Spec 19-01 re-pinned under the TrueCentroid default: 50% -> 749 (was 750), 10% -> 388 (was 289).
    // Spec 19-03 re-pin (injectivity guard ON + conservative packer silhouette): 50% -> 749, 10% -> 393.
    // Bijective-repair re-pin (folded charts re-flattened, not split; erase needs the ridged retry):
    // 50% -> 749, 10% -> 387, build-rel.
    // Spec 19 A3 re-pin (bijective flattener is the default for every chart): 50% 749 -> 750,
    // 10% 387 -> 381, build-rel.
    // reduceSeamAware decimates to AT MOST the target (seam_aware_reduce.h:70, loop at
    // seam_aware_reduce.cpp:730), so 750 was a measured value, not a contract.
    // Spec 21 re-pin (spike-shift footprint changes packing): 10% 381 -> 389.
    checkModel(modelOrSkip("bunny.iv"), "bunny", true, { 750, 389 });
}
