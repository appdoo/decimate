// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file segmenter_detail.h
/// @brief Private: the segmenter's boundary bookkeeping and its bounded fixpoint driver, exposed
///        for white-box tests (stub-chain depth, iteration cap).

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <decimate/mesh/half_edge_mesh.h>

namespace decimate::chart::detail {

/// @brief Run @p step until it returns false (fixpoint) or @p cap passes have been executed.
/// @return The number of passes executed; equals @p cap exactly when the cap stopped the loop.
///
/// A bound on `while (step()) {}`: a step that keeps flipping the same elements back and forth
/// never reaches a fixpoint, and an unbounded loop would hang the caller.
template <class Step>
int runBounded(Step &&step, int cap) {
    int passes = 0;
    while (passes < cap) {
        ++passes;
        if (!step()) {
            break;
        }
    }
    return passes;
}


/// Boundary flag per half-edge plus vertex boundary degrees (port of 2003 `CChartBoundary`,
/// `legacy2003/chartification/cchartboundary.cpp`).  The full reset per `expandCharts` (O(numHalfedges)) is
/// preserved (Beat 2a); a generation-counter optimisation is a Beat-2b forbidden path.
struct BoundaryFlags {
    std::vector<bool>             flags;       ///< One flag per half-edge: true = boundary.
    std::vector<int>              degrees;     ///< Per-vertex boundary degree count.
    std::vector<int>              baseDegrees; ///< Full-mesh degrees (all edges boundary); constant per mesh.
    const decimate::HalfEdgeMesh *mesh = nullptr;

    void init(const decimate::HalfEdgeMesh &m) {
        mesh = &m;
        const int nhe = static_cast<int>(m.halfEdgeCount());
        const int np = static_cast<int>(m.pointCount());
        flags.assign(static_cast<std::size_t>(nhe), false);
        degrees.assign(static_cast<std::size_t>(np), 0);
        computeBaseDegrees();
    }

    /// @brief Reset to full-boundary state — O(numHalfedges).
    ///
    /// Mirrors CChartBoundary::Reset(): all flags = true; degrees computed
    /// from unique undirected edges (each edge processed once via the "lower" rule).
    void reset() {
        std::fill(flags.begin(), flags.end(), true);
        degrees = baseDegrees;
    }

    /// Full-mesh degree table, computed once in init().  Each undirected edge is
    /// processed once via the "lower face" rule.
    void computeBaseDegrees() {
        const int nhe = static_cast<int>(mesh->halfEdgeCount());
        const int np = static_cast<int>(mesh->pointCount());
        baseDegrees.assign(static_cast<std::size_t>(np), 0);

        std::vector<bool> visited(static_cast<std::size_t>(nhe), false);

        for (int i = 0; i < nhe; ++i) {
            if (visited[static_cast<std::size_t>(i)]) {
                continue;
            }
            decimate::HeIndex       lo = static_cast<decimate::HeIndex>(i);
            const decimate::HeIndex s = mesh->sym(lo);
            if (s != decimate::kNoHalfedge &&
                mesh->halfedge(s).face < mesh->halfedge(lo).face) {
                lo = s;
            }
            visited[static_cast<std::size_t>(lo)] = true;
            if (s != decimate::kNoHalfedge) {
                visited[static_cast<std::size_t>(s)] = true;
            }

            const auto         &he = mesh->halfedge(lo);
            const std::uint32_t endV =
                mesh->halfedge(mesh->next(mesh->next(lo))).end;
            baseDegrees[static_cast<std::size_t>(he.end)]++;
            baseDegrees[static_cast<std::size_t>(endV)]++;
        }
    }

    /// True when @p h (or its twin) is still flagged boundary.
    [[nodiscard]] bool isBoundary(decimate::HeIndex h) const noexcept {
        return flags[static_cast<std::size_t>(lower(h))];
    }

    /// @brief Erase half-edge h from boundary — mirrors CChartBoundary::Erase.
    ///
    /// 2003 recurses once per peeled stub edge; here the recursion is an explicit worklist with
    /// the same depth-first order and the same state at every decision, so a stub chain of any
    /// length costs heap, not stack.
    void erase(decimate::HeIndex h) {
        // A frame is one pending erase; `resume` = its end-vertex subtree is done, the start-vertex check is next.
        struct Frame {
            decimate::HeIndex h;
            bool              resume;
        };
        std::vector<Frame> work;
        work.push_back({ h, false });

        while (!work.empty()) {
            const Frame f = work.back();
            work.pop_back();

            if (f.resume) {
                const decimate::HeIndex nn = mesh->next(mesh->next(f.h));
                const std::uint32_t     startV = mesh->halfedge(nn).end;
                if (degrees[static_cast<std::size_t>(startV)] < 2) {
                    pushStub(searchBoundary(nn), work);
                }
                continue;
            }

            const decimate::HeIndex lo = lower(f.h);
            if (!flags[static_cast<std::size_t>(lo)]) {
                continue;
            }
            flags[static_cast<std::size_t>(lo)] = false;
            const decimate::HeIndex s = mesh->sym(lo);
            if (s != decimate::kNoHalfedge) {
                flags[static_cast<std::size_t>(s)] = false;
            }

            const std::uint32_t endV = mesh->halfedge(lo).end;
            const std::uint32_t startV = mesh->halfedge(mesh->next(mesh->next(lo))).end;

            // The edge is gone from both end vertices' boundary degree.
            if (degrees[static_cast<std::size_t>(endV)] > 0) {
                degrees[static_cast<std::size_t>(endV)]--;
            }
            if (degrees[static_cast<std::size_t>(startV)] > 0) {
                degrees[static_cast<std::size_t>(startV)]--;
            }

            // A vertex left with fewer than two boundary edges would end a dangling stub: peel the
            // next boundary edge around it as well, so boundaries stay closed loops.  The start-vertex
            // check runs after the end-vertex subtree, as the recursive original ordered it.
            work.push_back({ lo, true });
            if (degrees[static_cast<std::size_t>(endV)] < 2) {
                pushStub(searchBoundary(lo), work);
            }
        }
    }

    /// @brief Direct set — mirrors CChartBoundary::SetBoundaryFlag.
    void setFlag(decimate::HeIndex h, bool boundary) {
        if (boundary == isBoundary(h)) {
            return;
        }
        const decimate::HeIndex lo = lower(h);
        flags[static_cast<std::size_t>(lo)] = boundary;
        const decimate::HeIndex s = mesh->sym(lo);
        if (s != decimate::kNoHalfedge) {
            flags[static_cast<std::size_t>(s)] = boundary;
        }
        const auto         &he = mesh->halfedge(lo);
        const int           delta = boundary ? 1 : -1;
        const std::uint32_t endV = he.end;
        const std::uint32_t startV = mesh->halfedge(mesh->next(mesh->next(lo))).end;
        degrees[static_cast<std::size_t>(endV)] += delta;
        degrees[static_cast<std::size_t>(startV)] += delta;
    }

private:
    /// The half-edge of the pair {h, sym(h)} that lies in the lower face; flags live on it alone.
    [[nodiscard]] decimate::HeIndex lower(decimate::HeIndex h) const noexcept {
        const decimate::HeIndex s = mesh->sym(h);
        if (s != decimate::kNoHalfedge &&
            mesh->halfedge(s).face < mesh->halfedge(h).face) {
            return s;
        }
        return h;
    }

    /// Queue the stub edge @p inv found by `searchBoundary` for erasure; `kNoHalfedge` = nothing to peel.
    template <class Work>
    static void pushStub(decimate::HeIndex inv, Work &work) {
        if (inv != decimate::kNoHalfedge) {
            work.push_back({ inv, false });
        }
    }

    /// Find a boundary half-edge in the 1-ring around `start`.
    /// Mirrors CChartBoundary::SearchBoundary.
    [[nodiscard]] decimate::HeIndex searchBoundary(decimate::HeIndex start) const {
        // Walk the 1-ring one way (next->sym); a mesh boundary ends the walk early.
        decimate::HeIndex loop = start;
        bool              done1 = false;
        do {
            if (isBoundary(loop)) {
                return loop;
            }
            const decimate::HeIndex nxt = mesh->next(loop);
            if (isBoundary(nxt)) {
                return nxt;
            }
            const decimate::HeIndex nxtSym = mesh->sym(nxt);
            if (nxtSym == decimate::kNoHalfedge) {
                done1 = true;
                break;
            }
            loop = nxtSym;
        } while (loop != start);

        if (!done1 && loop == start) {
            return decimate::kNoHalfedge;
        }

        // The first walk stopped at a mesh boundary: cover the other side of the ring too.
        loop = start;
        do {
            if (isBoundary(loop)) {
                return loop;
            }
            const decimate::HeIndex s = mesh->sym(loop);
            if (s != decimate::kNoHalfedge && isBoundary(s)) {
                return s;
            }
            if (s != decimate::kNoHalfedge) {
                loop = mesh->next(mesh->next(s));
            } else {
                break;
            }
        } while (loop != start);

        return decimate::kNoHalfedge;
    }
};


} // namespace decimate::chart::detail
