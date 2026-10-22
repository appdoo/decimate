// This file is part of the decimate library; see decimate.h for version/license details.
/// @file straighten_atlas.cpp
/// @brief Atlas analysis + protection rules (Liu et al. 2017, section 6, rules 1-4).

#include "straighten_atlas.h"

#include <algorithm>
#include <cstddef>
#include <decimate/seam/seam_topology.h>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace decimate::seam::detail {

/// Index of corner @p corner (taken mod 3) of face @p face in `faceV` / `faceU`.
static std::size_t slot(int face, int corner) {
    return (3 * static_cast<std::size_t>(face)) + static_cast<std::size_t>(corner % 3);
}


/// The silhouette edge for half-edge @p he, linked to @p twin (-1 for a 3D boundary edge).
static SilEdge edgeOf(const TexturedMesh &tm, HalfEdgeRef he, int twin) {
    return { tm.faceU[slot(he.face, he.corner)], tm.faceU[slot(he.face, he.corner + 1)], twin,
             he.face };
}


/// Surface vertex at corner `he.corner + offset` of the half-edge's face.
static int surfaceVertex(const TexturedMesh &tm, HalfEdgeRef he, int offset) {
    return tm.faceV[slot(he.face, he.corner + offset)];
}


/// Fills `svOf` and `copies`; throws when a UV vertex belongs to two surface vertices.
static void mapUvToSurface(const TexturedMesh &tm, Atlas &atlas) {
    atlas.svOf.assign(tm.uvs.size(), -1);
    atlas.copies.assign(tm.positions.size(), {});
    for (std::size_t i = 0; i < tm.faceU.size(); ++i) {
        const auto u = static_cast<std::size_t>(tm.faceU[i]);
        if (atlas.svOf[u] == -1) {
            atlas.svOf[u] = tm.faceV[i];
            atlas.copies[static_cast<std::size_t>(tm.faceV[i])].push_back(tm.faceU[i]);
        } else if (atlas.svOf[u] != tm.faceV[i]) {
            throw std::invalid_argument(
                "straighten: a UV vertex is shared by two surface vertices");
        }
    }
}


/// Fills `edges`, `incident` and the per-surface-vertex seam / boundary counts.
static void collectSilhouette(const TexturedMesh &tm, const SeamTopology &topo, Atlas &atlas) {
    atlas.seamCount.assign(tm.positions.size(), 0);
    atlas.boundaryCount.assign(tm.positions.size(), 0);

    // Seam half-edges are pushed pairwise so that `twin` is simply the neighbouring index.
    for (const SeamPair &pair : topo.seams) {
        const int ia = static_cast<int>(atlas.edges.size());
        atlas.edges.push_back(edgeOf(tm, pair.a, ia + 1));
        atlas.edges.push_back(edgeOf(tm, pair.b, ia));
        for (int k = 0; k < 2; ++k) {
            ++atlas.seamCount[static_cast<std::size_t>(surfaceVertex(tm, pair.a, k))];
        }
    }
    for (const HalfEdgeRef &he : topo.boundary) {
        atlas.edges.push_back(edgeOf(tm, he, -1));
        for (int k = 0; k < 2; ++k) {
            ++atlas.boundaryCount[static_cast<std::size_t>(surfaceVertex(tm, he, k))];
        }
    }

    // Per-vertex adjacency lets the chain walk step from an edge to the next one.
    atlas.incident.assign(tm.uvs.size(), {});
    for (std::size_t e = 0; e < atlas.edges.size(); ++e) {
        atlas.incident[static_cast<std::size_t>(atlas.edges[e].u0)].push_back(static_cast<int>(e));
        atlas.incident[static_cast<std::size_t>(atlas.edges[e].u1)].push_back(static_cast<int>(e));
    }
}


/// Union-find root of @p x with path halving.
static int findRoot(std::vector<int> &parent, int x) {
    while (parent[static_cast<std::size_t>(x)] != x) {
        parent[static_cast<std::size_t>(x)] =
            parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
        x = parent[static_cast<std::size_t>(x)];
    }
    return x;
}


/// Pack a directed vertex pair into one hash key.
static std::uint64_t edgeKey(int from, int to) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(from)) << 32U) |
           static_cast<std::uint32_t>(to);
}


/// Charts = faces connected across surface edges that are not seams.
static void buildCharts(const TexturedMesh &tm, const SeamTopology &topo, Atlas &atlas) {
    const std::size_t         faces = tm.faceCount();
    std::vector<std::uint8_t> isSeamHalf(3 * faces, 0);
    for (const SeamPair &pair : topo.seams) {
        isSeamHalf[slot(pair.a.face, pair.a.corner)] = 1;
        isSeamHalf[slot(pair.b.face, pair.b.corner)] = 1;
    }

    // Index every directed edge by its vertex pair to find each half-edge's opposite.
    std::unordered_map<std::uint64_t, int> directed;
    directed.reserve(3 * faces);
    for (std::size_t h = 0; h < 3 * faces; ++h) {
        directed[edgeKey(tm.faceV[h], tm.faceV[(h / 3 * 3) + ((h + 1) % 3)])] = static_cast<int>(h);
    }

    // Merge the two faces on either side of every non-seam interior edge.
    std::vector<int> parent(faces);
    std::iota(parent.begin(), parent.end(), 0);
    for (std::size_t h = 0; h < 3 * faces; ++h) {
        const auto opp = directed.find(edgeKey(tm.faceV[(h / 3 * 3) + ((h + 1) % 3)], tm.faceV[h]));
        if (opp == directed.end() || isSeamHalf[h] != 0) {
            continue;
        }
        parent[static_cast<std::size_t>(findRoot(parent, static_cast<int>(h / 3)))] =
            findRoot(parent, opp->second / 3);
    }

    // Number the components in order of their first face.
    atlas.chartOfFace.assign(faces, -1);
    std::vector<int> chartOfRoot(faces, -1);
    for (std::size_t f = 0; f < faces; ++f) {
        int &id = chartOfRoot[static_cast<std::size_t>(findRoot(parent, static_cast<int>(f)))];
        if (id < 0) {
            id = static_cast<int>(atlas.chartFaces.size());
            atlas.chartFaces.emplace_back();
        }
        atlas.chartOfFace[f] = id;
        atlas.chartFaces[static_cast<std::size_t>(id)].push_back(static_cast<int>(f));
    }
}


/// Ears (rule 3): faces with exactly two silhouette half-edges.
static void findEars(const TexturedMesh &tm, Atlas &atlas) {
    // Remember the first three silhouette edges of each face; more than two disqualifies it anyway.
    std::vector<std::array<int, 3>> perFace(tm.faceCount(), { -1, -1, -1 });
    std::vector<int>                count(tm.faceCount(), 0);
    for (std::size_t e = 0; e < atlas.edges.size(); ++e) {
        const auto f = static_cast<std::size_t>(atlas.edges[e].face);
        if (count[f] < 3) {
            perFace[f][static_cast<std::size_t>(count[f])] = static_cast<int>(e);
        }
        ++count[f];
    }

    for (std::size_t f = 0; f < tm.faceCount(); ++f) {
        if (count[f] != 2) {
            continue;
        }
        atlas.earFaces.push_back(static_cast<int>(f));

        // Only ears whose both edges are seams have neighbour charts to compare.
        const SilEdge &e0 = atlas.edges[static_cast<std::size_t>(perFace[f][0])];
        const SilEdge &e1 = atlas.edges[static_cast<std::size_t>(perFace[f][1])];
        if (e0.twin < 0 || e1.twin < 0) {
            continue;
        }
        const int c0 = atlas.chartOfFace[static_cast<std::size_t>(
            atlas.edges[static_cast<std::size_t>(e0.twin)].face)];
        const int c1 = atlas.chartOfFace[static_cast<std::size_t>(
            atlas.edges[static_cast<std::size_t>(e1.twin)].face)];
        atlas.sameNeighbourEars += c0 == c1 ? 1 : 0;
    }
}


/// OR @p bit into the mask of every UV vertex in @p vertices.
static void setBits(std::vector<std::uint8_t> &mask, const std::vector<int> &vertices,
                    std::uint8_t bit) {
    for (const int u : vertices) {
        mask[static_cast<std::size_t>(u)] |= bit;
    }
}


/// Rules (1) and (2), which look at the surface vertex.
static void surfaceRules(const Atlas &atlas, std::vector<std::uint8_t> &mask) {
    for (std::size_t sv = 0; sv < atlas.copies.size(); ++sv) {
        const int seams = atlas.seamCount[sv];
        if (seams != 0 && seams != 2) {
            setBits(mask, atlas.copies[sv], 1); // rule 1
        }
        if (seams >= 1 && atlas.boundaryCount[sv] >= 1) {
            setBits(mask, atlas.copies[sv], 2); // rule 2
        }
    }
}


/// A chain needs exactly two silhouette edges through a vertex; anything else (e.g. a fan
/// whose two sides share one UV vertex) cannot be walked, so it joins rule (1).
static void degreeGuard(const Atlas &atlas, std::vector<std::uint8_t> &mask) {
    for (std::size_t u = 0; u < mask.size(); ++u) {
        const std::size_t degree = atlas.incident[u].size();
        mask[u] |= degree != 0 && degree != 2 ? 1 : 0;
    }
}


/// Rule (3): the corners of every ear face.
static void earRule(const Atlas &atlas, std::vector<std::uint8_t> &mask) {
    for (const int f : atlas.earFaces) {
        for (int k = 0; k < 3; ++k) {
            mask[static_cast<std::size_t>(atlas.tm->faceU[slot(f, k)])] |= 4;
        }
    }
}


/// Rule (4): protection travels to every UV copy of a protected surface vertex.
static void propagateRule(const Atlas &atlas, const std::vector<std::uint8_t> &forced,
                          std::vector<std::uint8_t> &mask) {
    for (const auto &copies : atlas.copies) {
        const bool any = std::ranges::any_of(copies, [&](int u) {
            return (mask[static_cast<std::size_t>(u)] & kRuleBits) != 0 ||
                   forced[static_cast<std::size_t>(u)] != 0;
        });

        // Only copies that no other rule protected get the propagation bit.
        for (const int u : copies) {
            auto &m = mask[static_cast<std::size_t>(u)];
            m |= any && (m & kRuleBits) == 0 ? 8 : 0;
        }
    }
}


Atlas analyse(const TexturedMesh &tm) {
    Atlas atlas;
    atlas.tm = &tm;
    mapUvToSurface(tm, atlas);

    // `findSeams` is the single seam predicate; everything else derives from its topology.
    const SeamTopology topo = findSeams(tm);
    collectSilhouette(tm, topo, atlas);
    buildCharts(tm, topo, atlas);
    findEars(tm, atlas);
    return atlas;
}


std::vector<std::uint8_t> protectionMasks(const Atlas                     &atlas,
                                          const std::vector<std::uint8_t> &forced) {
    std::vector<std::uint8_t> mask(atlas.tm->uvs.size(), 0);

    // Rule (4) must run last: it only tops up vertices the other rules left unprotected.
    surfaceRules(atlas, mask);
    degreeGuard(atlas, mask);
    earRule(atlas, mask);
    propagateRule(atlas, forced, mask);
    return mask;
}

} // namespace decimate::seam::detail
