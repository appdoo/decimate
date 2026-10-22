// This file is part of the decimate library; see decimate.h for version/license details.
/// @file straighten_chains.cpp
/// @brief Edge-chain components, RDP coarsening and arc-length reparametrization (section 6).

#include "straighten_chains.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>

namespace decimate::seam::detail {

namespace {

constexpr std::size_t kMinLoopVertices = 3;


/// Index range [lo, hi] of a polyline.
struct Span {
    std::size_t lo;
    std::size_t hi;
};


/// One chain of silhouette edges through non-protected degree-2 vertices.
struct Chain {
    std::vector<int> verts; ///< UV vertices; closed loops repeat the first vertex last.
    std::vector<int> edges; ///< Silhouette edges, `edges[i]` joins `verts[i]` and `verts[i+1]`.
    bool             closed = false;
};


/// All chains of an atlas plus the chain each silhouette edge belongs to.
struct ChainBuild {
    std::vector<Chain> chains;
    std::vector<int>   chainOfEdge;
};

} // namespace


/// Euclidean distance in all four components.
static double dist(const Pt4 &a, const Pt4 &b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        sum += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return std::sqrt(sum);
}


/// Distance from @p p to the segment a-b (to `a` when the segment is a point).
static double segmentDistance(const Pt4 &p, const Pt4 &a, const Pt4 &b) {
    double len2 = 0.0;
    double dot = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        len2 += (b[i] - a[i]) * (b[i] - a[i]);
        dot += (p[i] - a[i]) * (b[i] - a[i]);
    }
    if (len2 <= 0.0) {
        return dist(p, a);
    }

    // Project onto the segment and clamp, so the distance is to the segment, not its line.
    const double t = std::clamp(dot / len2, 0.0, 1.0);
    Pt4          foot {};
    for (std::size_t i = 0; i < 4; ++i) {
        foot[i] = a[i] + (t * (b[i] - a[i]));
    }
    return dist(p, foot);
}


/// Farthest point strictly inside (lo, hi) from the segment lo-hi.
static std::pair<std::size_t, double> farthest(const std::vector<Pt4> &pts, Span s) {
    std::size_t best = s.lo;
    double      bestDist = -1.0;
    for (std::size_t i = s.lo + 1; i < s.hi; ++i) {
        const double d = segmentDistance(pts[i], pts[s.lo], pts[s.hi]);
        if (d > bestDist) {
            bestDist = d;
            best = i;
        }
    }
    return { best, bestDist };
}


/// RDP on [lo, hi]: both ends are already flagged in @p keep.
static void refine(const std::vector<Pt4> &pts, Span first, double tau,
                   std::vector<std::uint8_t> &keep) {
    // An explicit stack instead of recursion: long chains would otherwise overflow it.
    std::vector<Span> stack { first };
    while (!stack.empty()) {
        const Span s = stack.back();
        stack.pop_back();
        if (s.hi <= s.lo + 1) {
            continue;
        }
        const auto [idx, d] = farthest(pts, s);
        if (d <= tau) {
            continue;
        }
        keep[idx] = 1;
        stack.push_back({ s.lo, idx });
        stack.push_back({ idx, s.hi });
    }
}


/// Kept vertices of a closed loop, not counting the repeated last one.
static std::size_t keptLoopVertices(const std::vector<std::uint8_t> &keep) {
    return static_cast<std::size_t>(std::count(keep.begin(), keep.end() - 1, 1));
}


/// Adds the unkept loop vertex farthest from its enclosing kept segment.
static void keepOneMore(const std::vector<Pt4> &pts, std::vector<std::uint8_t> &keep) {
    std::size_t lo = 0;
    std::size_t bestIdx = 0;
    double      bestDist = -1.0;
    for (std::size_t i = 1; i < pts.size(); ++i) {
        if (keep[i] == 0) {
            continue;
        }
        const auto [idx, d] = farthest(pts, { lo, i });
        if (d > bestDist) {
            bestDist = d;
            bestIdx = idx;
        }
        lo = i;
    }
    keep[bestIdx] = 1;
}


/// Re-spaces the points strictly between kept points @p a and @p b along the segment a-b.
static void respaceSpan(const std::vector<Pt4> &pts, const std::vector<double> &arc, std::size_t a,
                        std::size_t b, std::vector<Pt4> &out) {
    const double span = arc[b] - arc[a];
    for (std::size_t i = a + 1; i < b; ++i) {
        // Zero-length span: fall back to uniform spacing instead of dividing by zero.
        const double t = span > 0.0 ? (arc[i] - arc[a]) / span
                                    : static_cast<double>(i - a) / static_cast<double>(b - a);
        for (std::size_t c = 0; c < 4; ++c) {
            out[i][c] = pts[a][c] + (t * (pts[b][c] - pts[a][c]));
        }
    }
}


/// The other end of silhouette edge @p e as seen from vertex @p v.
static int otherEnd(const SilEdge &e, int v) {
    return e.u0 == v ? e.u1 : e.u0;
}


/// The second silhouette edge at a degree-2 @p vertex, given the first one.
static int otherEdgeAt(const Atlas &atlas, int vertex, int edge) {
    const auto &inc = atlas.incident[static_cast<std::size_t>(vertex)];
    return inc[0] == edge ? inc[1] : inc[0];
}


/// Extends the chain from its current end through non-protected vertices.
/// Returns true if it came back to @p stop (a closed loop).
static bool extend(const Atlas &atlas, const std::vector<std::uint8_t> &isProtected, int stop,
                   std::deque<int> &verts, std::deque<int> &edges, bool forward) {
    int v = forward ? verts.back() : verts.front();
    int e = forward ? edges.back() : edges.front();
    while (isProtected[static_cast<std::size_t>(v)] == 0) {
        const int next = otherEdgeAt(atlas, v, e);
        if (next == stop) {
            return true;
        }
        v = otherEnd(atlas.edges[static_cast<std::size_t>(next)], v);
        e = next;
        if (forward) {
            edges.push_back(e);
            verts.push_back(v);
        } else {
            edges.push_front(e);
            verts.push_front(v);
        }
    }
    return false;
}


/// The maximal chain containing silhouette edge @p seed.
static Chain chainThrough(const Atlas &atlas, const std::vector<std::uint8_t> &isProtected,
                          int seed) {
    const SilEdge  &s = atlas.edges[static_cast<std::size_t>(seed)];
    std::deque<int> verts { s.u0, s.u1 };
    std::deque<int> edges { seed };
    Chain           chain;

    // Walk forward; if the walk returns to the seed the chain is a loop, else walk the other way.
    chain.closed = extend(atlas, isProtected, seed, verts, edges, true);
    if (chain.closed) {
        verts.push_back(verts.front());
    } else {
        extend(atlas, isProtected, seed, verts, edges, false);
    }
    chain.verts.assign(verts.begin(), verts.end());
    chain.edges.assign(edges.begin(), edges.end());
    return chain;
}


/// All edges seams, or all edges boundary: the protection rules guarantee it; checked anyway.
static bool uniformKind(const Atlas &atlas, const Chain &c) {
    const bool seam = atlas.edges[static_cast<std::size_t>(c.edges.front())].twin >= 0;
    return std::ranges::all_of(c.edges, [&](int e) {
        return (atlas.edges[static_cast<std::size_t>(e)].twin >= 0) == seam;
    });
}


/// True when the chain has at least one vertex the straightening may move.
static bool hasInterior(const Chain &c) {
    return c.closed ? c.edges.size() >= kMinLoopVertices : c.edges.size() >= 2;
}


/// Partition every silhouette edge into chains.
static ChainBuild buildChains(const Atlas &atlas, const std::vector<std::uint8_t> &isProtected) {
    ChainBuild out;
    out.chainOfEdge.assign(atlas.edges.size(), -1);
    for (std::size_t e = 0; e < atlas.edges.size(); ++e) {
        if (out.chainOfEdge[e] >= 0) {
            continue;
        }
        out.chains.push_back(chainThrough(atlas, isProtected, static_cast<int>(e)));
        for (const int member : out.chains.back().edges) {
            out.chainOfEdge[static_cast<std::size_t>(member)] =
                static_cast<int>(out.chains.size()) - 1;
        }
    }
    return out;
}


/// The UV vertex @p u as a 4D point with zero third and fourth components.
static Pt4 point2(const TexturedMesh &tm, int u) {
    const auto &p = tm.uvs[static_cast<std::size_t>(u)];
    return { static_cast<double>(p.x), static_cast<double>(p.y), 0.0, 0.0 };
}


/// Record the new position of UV vertex @p u.
static void store(ChainTargets &t, int u, double x, double y) {
    t.pos[static_cast<std::size_t>(u)] = { x, y };
    t.moved[static_cast<std::size_t>(u)] = 1;
}


/// A boundary chain: 2D RDP, vertices re-spaced by normalized original edge length.
static void straightenBoundary(const Atlas &atlas, const Chain &chain, double tau,
                               ChainTargets &t) {
    std::vector<Pt4> pts;
    for (const int u : chain.verts) {
        pts.push_back(point2(*atlas.tm, u));
    }
    const auto keep = coarsen(pts, chain.closed, tau);
    const auto moved = reparametrize(pts, keep);

    // A closed loop repeats its first vertex last; do not store it twice.
    const std::size_t count = chain.closed ? chain.verts.size() - 1 : chain.verts.size();
    for (std::size_t i = 0; i < count; ++i) {
        if (keep[i] == 0) {
            store(t, chain.verts[i], moved[i][0], moved[i][1]);
        }
    }
    ++t.boundaryChains;
}


/// The vertex on the twin side that sits on the same surface vertex as @p u.
static int copyOn(const Atlas &atlas, const SilEdge &twin, int u) {
    const int sv = atlas.svOf[static_cast<std::size_t>(u)];
    return atlas.svOf[static_cast<std::size_t>(twin.u0)] == sv ? twin.u0 : twin.u1;
}


/// The vertices of the copy chain, aligned with `chain.verts`.
static std::vector<int> copyVerts(const Atlas &atlas, const Chain &chain) {
    std::vector<int> other;
    for (std::size_t i = 0; i < chain.edges.size(); ++i) {
        const SilEdge &twin = atlas.edges[static_cast<std::size_t>(
            atlas.edges[static_cast<std::size_t>(chain.edges[i])].twin)];
        other.push_back(copyOn(atlas, twin, chain.verts[i]));
    }

    // The last vertex has no outgoing edge of its own; take it from the last edge's twin.
    const SilEdge &last = atlas.edges[static_cast<std::size_t>(
        atlas.edges[static_cast<std::size_t>(chain.edges.back())].twin)];
    other.push_back(copyOn(atlas, last, chain.verts.back()));
    return other;
}


/// A seam chain and its copy: 4D RDP, both sides re-spaced by the same 4D arc length.
static void straightenSeam(const Atlas &atlas, const Chain &chain, double tau, ChainTargets &t) {
    // Concatenating both copies makes one RDP decision (and one arc length) serve both sides.
    const std::vector<int> other = copyVerts(atlas, chain);
    std::vector<Pt4>       pts;
    for (std::size_t i = 0; i < chain.verts.size(); ++i) {
        const Pt4 a = point2(*atlas.tm, chain.verts[i]);
        const Pt4 b = point2(*atlas.tm, other[i]);
        pts.push_back({ a[0], a[1], b[0], b[1] });
    }
    const auto keep = coarsen(pts, chain.closed, tau);
    const auto moved = reparametrize(pts, keep);

    const std::size_t count = chain.closed ? chain.verts.size() - 1 : chain.verts.size();
    for (std::size_t i = 0; i < count; ++i) {
        if (keep[i] != 0) {
            continue;
        }
        store(t, chain.verts[i], moved[i][0], moved[i][1]);
        store(t, other[i], moved[i][2], moved[i][3]);
    }
    ++t.seamChains;
}


std::vector<std::uint8_t> coarsen(const std::vector<Pt4> &pts, bool closed, double tau) {
    const std::size_t         n = pts.size();
    std::vector<std::uint8_t> keep(n, 0);

    // Fewer than two points: nothing to coarsen.
    if (n < 2) {
        keep.assign(n, 1);
        return keep;
    }
    keep.front() = 1;
    keep.back() = 1;
    if (!closed) {
        refine(pts, { 0, n - 1 }, tau, keep);
        return keep;
    }

    // A loop has no natural end points: split it at the vertex farthest from the start.
    std::size_t far = 0;
    for (std::size_t i = 1; i + 1 < n; ++i) {
        far = dist(pts[i], pts[0]) > dist(pts[far], pts[0]) ? i : far;
    }
    keep[far] = 1;
    refine(pts, { 0, far }, tau, keep);
    refine(pts, { far, n - 1 }, tau, keep);

    // Fewer than three vertices would collapse the loop to a line or a point.
    while (keptLoopVertices(keep) < std::min(kMinLoopVertices, n - 1)) {
        keepOneMore(pts, keep);
    }
    return keep;
}


std::vector<Pt4> reparametrize(const std::vector<Pt4> &pts, const std::vector<std::uint8_t> &keep) {
    // Cumulative original arc length: the parameter along which dropped points are re-spaced.
    std::vector<Pt4>    out = pts;
    std::vector<double> arc(pts.size(), 0.0);
    for (std::size_t i = 1; i < pts.size(); ++i) {
        arc[i] = arc[i - 1] + dist(pts[i - 1], pts[i]);
    }

    std::size_t a = 0;
    for (std::size_t b = 1; b < pts.size(); ++b) {
        if (keep[b] == 0) {
            continue;
        }
        respaceSpan(pts, arc, a, b, out);
        a = b;
    }
    return out;
}


ChainTargets straightenChains(const Atlas &atlas, const std::vector<std::uint8_t> &isProtected,
                              double tau) {
    ChainTargets t;
    t.pos.assign(atlas.tm->uvs.size(), { 0.0, 0.0 });
    t.moved.assign(atlas.tm->uvs.size(), 0);

    const ChainBuild          build = buildChains(atlas, isProtected);
    std::vector<std::uint8_t> done(build.chains.size(), 0);
    for (std::size_t c = 0; c < build.chains.size(); ++c) {
        const Chain &chain = build.chains[c];
        if (done[c] != 0 || !hasInterior(chain) || !uniformKind(atlas, chain)) {
            continue;
        }
        done[c] = 1;

        // A seam chain is handled together with its copy, so mark the copy done as well.
        const int twin = atlas.edges[static_cast<std::size_t>(chain.edges.front())].twin;
        if (twin < 0) {
            straightenBoundary(atlas, chain, tau, t);
            continue;
        }
        done[static_cast<std::size_t>(build.chainOfEdge[static_cast<std::size_t>(twin)])] = 1;
        straightenSeam(atlas, chain, tau, t);
    }
    return t;
}

} // namespace decimate::seam::detail
