// This file is part of the decimate library; see decimate.h for version/license details.
/// @file injectivity.cpp
/// @brief Spec 19-03 — detect LSCM self-overlap per chart and repair it.
///
/// The 2003 thesis names the defect (`diploma_tex/chartification.tex:98`) and only mitigates it
/// heuristically (`Straighten`, `:100`).  Here it is detected exactly (triangle-pair interior
/// overlap in UV) and repaired: by default the chart is re-flattened with the bijective
/// free-boundary method (Smith & Schaefer 2015, `flattenBijective`); a chart that method refuses
/// (not a disk) — or every folded chart under `Flattener::LscmSplit` — is split at the fold by
/// `splitChart`, and only the two pieces are re-flattened.

#include "chart_cut.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <decimate/chart/parametrize.h>
#include <limits>
#include <utility>
#include <vector>

namespace decimate::chart {

namespace {

/// Vertex-disjoint cut paths added per stretch round.
constexpr int kPathsPerRound = 4;


/// One UV triangle in double precision, so the overlap test is not limited by float rounding.
struct Tri {
    std::array<double, 3> x {};
    std::array<double, 3> y {};
};


/// True when the projections of @p a and @p b onto the normal of one edge of @p t leave a
/// gap >= -eps (eps in UV units).
bool separatedByEdgeOf(const Tri &t, const Tri &a, const Tri &b, double eps) {
    for (int e = 0; e < 3; ++e) {
        const auto n = static_cast<std::size_t>((e + 1) % 3);
        const auto i = static_cast<std::size_t>(e);

        // Project both triangles onto the edge normal; a gap between the two intervals separates them.
        const double nx = -(t.y[n] - t.y[i]);
        const double ny = t.x[n] - t.x[i];
        double       minA = std::numeric_limits<double>::max();
        double       maxA = std::numeric_limits<double>::lowest();
        double       minB = minA;
        double       maxB = maxA;
        for (std::size_t k = 0; k < 3; ++k) {
            const double pa = nx * a.x[k] + ny * a.y[k];
            const double pb = nx * b.x[k] + ny * b.y[k];
            minA = std::min(minA, pa);
            maxA = std::max(maxA, pa);
            minB = std::min(minB, pb);
            maxB = std::max(maxB, pb);
        }

        // The normal is unnormalised, so the tolerance scales with its length.
        const double slack = eps * std::sqrt(nx * nx + ny * ny);
        if (maxA <= minB + slack || maxB <= minA + slack) {
            return true;
        }
    }
    return false;
}


/// Separating-axis test over the six edge normals of the two triangles.
bool separated(const Tri &a, const Tri &b, double eps) {
    return separatedByEdgeOf(a, a, b, eps) || separatedByEdgeOf(b, a, b, eps);
}


/// Signed area of @p t; the sign is the UV orientation, ~0 means degenerate.
double signedArea(const Tri &t) {
    return 0.5 * ((t.x[1] - t.x[0]) * (t.y[2] - t.y[0]) - (t.x[2] - t.x[0]) * (t.y[1] - t.y[0]));
}


/// UV triangles of a chart (local face order) and the UV bounding box.
struct UvTris {
    std::vector<Tri>     tri;
    double               x0 = std::numeric_limits<double>::max();
    double               y0 = std::numeric_limits<double>::max();
    double               x1 = std::numeric_limits<double>::lowest();
    double               y1 = std::numeric_limits<double>::lowest();
    [[nodiscard]] double diag() const {
        return std::hypot(x1 - x0, y1 - y0);
    }
};


/// Read the chart's UV triangles per corner (so a self-seam cut resolves) and their bounding box.
UvTris loadTriangles(const HalfEdgeChart &ac) {
    UvTris out;
    out.tri.resize(ac.faces.size());
    for (std::size_t i = 0; i < ac.faces.size(); ++i) {
        for (std::size_t k = 0; k < 3; ++k) {
            const auto  local = localCorner(ac, static_cast<int>(i), static_cast<int>(k));
            const auto &uv = ac.uv[static_cast<std::size_t>(local)];
            out.tri[i].x[k] = uv.x;
            out.tri[i].y[k] = uv.y;
            out.x0 = std::min(out.x0, static_cast<double>(uv.x));
            out.y0 = std::min(out.y0, static_cast<double>(uv.y));
            out.x1 = std::max(out.x1, static_cast<double>(uv.x));
            out.y1 = std::max(out.y1, static_cast<double>(uv.y));
        }
    }
    return out;
}


/// Uniform g x g grid over the UV box; each non-degenerate triangle is listed in every cell its
/// bounding box touches, and remembers the first cell (box[0], box[1]) it touches.
struct Grid {
    int                             g = 1;
    std::vector<std::array<int, 4>> box;  ///< minCx, minCy, maxCx, maxCy per triangle
    std::vector<std::vector<int>>   cell; ///< g*g lists of triangle indices, row-major

    [[nodiscard]] std::vector<int> &at(int u, int v) {
        return cell[static_cast<std::size_t>(v) * static_cast<std::size_t>(g) +
                    static_cast<std::size_t>(u)];
    }

    /// List triangle @p i in every cell its bounding box touches.
    void insert(int i) {
        const auto &b = box[static_cast<std::size_t>(i)];
        for (int v = b[1]; v <= b[3]; ++v) {
            for (int u = b[0]; u <= b[2]; ++u) {
                at(u, v).push_back(i);
            }
        }
    }

    /// True when cell (u, v) is the first cell triangles @p i and @p j both touch.
    [[nodiscard]] bool firstCommonCell(std::size_t i, std::size_t j, int u, int v) const {
        return std::max(box[i][0], box[j][0]) == u && std::max(box[i][1], box[j][1]) == v;
    }
};


/// Bin the non-degenerate triangles into a grid with about one triangle per cell, to avoid the O(n^2) pair scan.
Grid buildGrid(const UvTris &t) {
    const std::size_t n = t.tri.size();
    Grid              gr;
    gr.g = std::max(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n)))));
    gr.box.resize(n);
    gr.cell.resize(static_cast<std::size_t>(gr.g) * static_cast<std::size_t>(gr.g));

    // The 1e-300 keeps a flat box (zero width or height) from dividing by zero.
    const double diag = t.diag();
    const auto   cx = [&](double v) {
        return std::clamp(static_cast<int>((v - t.x0) / (t.x1 - t.x0 + 1e-300) * gr.g), 0,
                            gr.g - 1);
    };
    const auto cy = [&](double v) {
        return std::clamp(static_cast<int>((v - t.y0) / (t.y1 - t.y0 + 1e-300) * gr.g), 0,
                          gr.g - 1);
    };

    for (std::size_t i = 0; i < n; ++i) {
        const Tri &tr = t.tri[i];
        if (!(std::fabs(signedArea(tr)) > 1e-12 * diag * diag)) {
            continue; // zero-area (or NaN-UV) triangles never overlap anything
        }
        const auto [ax, bx] = std::minmax({ tr.x[0], tr.x[1], tr.x[2] });
        const auto [ay, by] = std::minmax({ tr.y[0], tr.y[1], tr.y[2] });
        gr.box[i] = { cx(ax), cy(ay), cx(bx), cy(by) };
        gr.insert(static_cast<int>(i));
    }
    return gr;
}


/// Overlap census of one flattened chart.  `involved[i]` is set for local faces in a pair.
struct Census {
    int               pairs = 0;
    std::vector<char> involved;
};


/// Test every pair of cell (u, v) that has this cell as its FIRST common cell (each pair once).
void countCell(Grid &gr, const UvTris &t, double eps, int u, int v, Census &out) {
    const std::vector<int> &members = gr.at(u, v);
    for (std::size_t p = 0; p < members.size(); ++p) {
        const auto i = static_cast<std::size_t>(members[p]);
        for (std::size_t q = p + 1; q < members.size(); ++q) {
            const auto j = static_cast<std::size_t>(members[q]);
            const bool hit = gr.firstCommonCell(i, j, u, v) && !separated(t.tri[i], t.tri[j], eps);
            out.pairs += hit ? 1 : 0;
            out.involved[i] |= hit ? 1 : 0;
            out.involved[j] |= hit ? 1 : 0;
        }
    }
}


/// Count the overlapping triangle pairs of @p ac and mark the faces in them.
Census census(const HalfEdgeChart &ac) {
    Census out;
    out.involved.assign(ac.faces.size(), 0);
    if (ac.faces.size() < 2) {
        return out; // a single triangle cannot overlap itself
    }

    const UvTris t = loadTriangles(ac);
    if (!(t.diag() > 0.0) || !std::isfinite(t.diag())) {
        return out; // all UVs coincide, or none is finite: no usable box
    }

    // Eps is relative to the box, so the verdict is scale invariant.
    Grid gr = buildGrid(t);
    for (int v = 0; v < gr.g; ++v) {
        for (int u = 0; u < gr.g; ++u) {
            countCell(gr, t, 1e-6 * t.diag(), u, v, out);
        }
    }
    return out;
}


/// Face across edge @p e of face @p f, or -1 on a mesh boundary.
int faceAcross(const HalfEdgeMesh &mesh, int f, int e) {
    const HeIndex s = mesh.sym(static_cast<HeIndex>(f * 3 + e));
    return s == kNoHalfedge ? -1 : static_cast<int>(mesh.halfedge(s).face);
}


/// Dual-BFS distance (inside chart @p chartId) of every face from @p src; -1 outside.
std::vector<int> chartDistance(const HalfEdgeMesh &mesh, const ChartSet &cs, int chartId, int src) {
    std::vector<int> dist(cs.faceToChart.size(), -1);
    std::vector<int> queue { src };
    dist[static_cast<std::size_t>(src)] = 0;
    const auto claim = [&](int nb, int d) {
        if (nb < 0 || cs.faceToChart[static_cast<std::size_t>(nb)] != chartId ||
            dist[static_cast<std::size_t>(nb)] >= 0) {
            return;
        }
        dist[static_cast<std::size_t>(nb)] = d;
        queue.push_back(nb);
    };
    for (std::size_t head = 0; head < queue.size();
         ++head) { // NOLINT(modernize-loop-convert): queue grows while iterated
        const int f = queue[head];
        for (int e = 0; e < 3; ++e) {
            claim(faceAcross(mesh, f, e), dist[static_cast<std::size_t>(f)] + 1);
        }
    }
    return dist;
}


/// The involved face farthest from the chart's seed (ties: lowest face index); -1 if none.
int pickFoldFace(const HalfEdgeMesh &mesh, const ChartSet &cs, int chartId, const HalfEdgeChart &ac,
                 const Census &c) {
    const int              seed = cs.charts[static_cast<std::size_t>(chartId)].seed;
    const std::vector<int> dist = chartDistance(mesh, cs, chartId, seed);
    int                    best = -1;
    for (std::size_t i = 0; i < ac.faces.size(); ++i) {
        const int f = ac.faces[i];
        if (c.involved[i] == 0 || f == seed) {
            continue;
        }

        // Farthest from the seed wins; ties go to the lowest face index so the split is deterministic.
        const bool farther =
            best < 0 || dist[static_cast<std::size_t>(f)] > dist[static_cast<std::size_t>(best)];
        const bool tie =
            best >= 0 &&
            dist[static_cast<std::size_t>(f)] == dist[static_cast<std::size_t>(best)] && f < best;
        if (farther || tie) {
            best = f;
        }
    }
    return best;
}


/// The flattened charts of a ChartSet plus their overlap-pair counts, kept in step with the set.
struct Guard {
    const HalfEdgeMesh        &mesh;
    const MeshData            &data;
    ChartSet                  &charts;
    const InjectivityOptions  &opts;
    InjectivityStats          &st;
    std::vector<HalfEdgeChart> atlas;
    std::vector<int>           pairs;

    /// Flatten the chart just added to the set and record its overlap count.
    void append(std::size_t chart) {
        atlas.push_back(flatten(mesh, data, charts.charts[chart]));
        pairs.push_back(atlas.back().injective() ? 0 : selfOverlapPairs(atlas.back()));
    }
    /// Re-flatten a chart whose faces changed (after a split) and refresh its overlap count.
    void refresh(std::size_t chart) {
        atlas[chart] = flatten(mesh, data, charts.charts[chart]);
        // flatten's O(n + B^2) quality check decides; the pair search runs only to locate a fold.
        pairs[chart] = atlas[chart].injective() ? 0 : selfOverlapPairs(atlas[chart]);
    }
    /// Number of charts that still fold.
    [[nodiscard]] int foldedCount() const {
        return static_cast<int>(
            std::count_if(pairs.begin(), pairs.end(), [](int p) { return p > 0; }));
    }

    /// Re-flatten chart @p c bijectively (Smith & Schaefer 2015) as the strategy asks.  Bijective:
    /// every disk chart, started from `Current` (its LSCM map when injective, else Tutte; the
    /// caller's `bijective.start` is not consulted).  Lscm: only a folded chart, started from
    /// `bijective.start`.  LscmSplit / LscmRaw: never.  A refusal (non-disk) leaves it to the split.
    void heal(std::size_t c) {
        const bool       folded = pairs[c] > 0;
        BijectiveOptions bo = opts.bijective;
        switch (opts.flattener) {
            case Flattener::Bijective:
                bo.start = BijectiveStart::Current;
                break;
            case Flattener::Lscm:
                if (!folded) {
                    return;
                }
                break;
            case Flattener::LscmSplit:
            case Flattener::LscmRaw:
                return;
        }

        const auto t0 = std::chrono::steady_clock::now();
        const bool ran = flattenBijective(data, atlas[c], bo);
        st.bijectiveSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

        // A refusal (non-disk) leaves the LSCM map in place; the caller's split round handles it.
        const bool ok = ran && atlas[c].injective();
        st.flattenedBijectively += ran ? 1 : 0;
        st.repairedBijectively += folded && ok ? 1 : 0;
        st.bijectiveFailed += folded && !ok ? 1 : 0;
        if (ran) {
            pairs[c] = ok ? 0 : std::max(selfOverlapPairs(atlas[c]), 1);
        }
    }

    /// Split chart @p c at its fold face and re-flatten the two pieces; false if there is no fold
    /// face.
    bool repair(std::size_t c) {
        const int fold =
            pickFoldFace(mesh, charts, static_cast<int>(c), atlas[c], census(atlas[c]));
        if (fold < 0) {
            return false;
        }
        splitChart(mesh, data, charts, static_cast<int>(c), fold, charts.centroidMode);
        refresh(c);
        append(charts.charts.size() - 1);
        heal(c);
        heal(charts.charts.size() - 1);
        return true;
    }

    /// Cut chart @p c open for the stretch bound (Gu, Gortler, Hoppe 2002, section 3): add up to `kPathsPerRound`
    /// vertex-disjoint paths from its worst faces to the boundary and re-flatten (healed), while the chart is over the
    /// bound and a round still lowers its L2 stretch.  A round that does not is undone.
    void cutToBound(std::size_t c) {
        while (pairs[c] == 0) {
            std::vector<double> gamma;
            const ChartStretch  s = chartStretch(data, atlas[c], &gamma);
            const bool          over = (opts.maxStretchL2 > 0.0 && s.l2 > opts.maxStretchL2) ||
                                       (opts.maxStretchLinf > 0.0 && s.linf > opts.maxStretchLinf);
            if (!over) {
                return;
            }

            // Worst faces first, in mesh face ids; a skipped (flipped) triangle has gamma < 0.
            std::vector<int> worst;
            for (std::size_t i = 0; i < gamma.size(); ++i) {
                if (gamma[i] >= 0.0) {
                    worst.push_back(static_cast<int>(i));
                }
            }
            std::ranges::stable_sort(worst, [&](int a, int b) { return gamma[static_cast<std::size_t>(a)] > gamma[static_cast<std::size_t>(b)]; });
            for (int &i : worst) {
                i = atlas[c].faces[static_cast<std::size_t>(i)];
            }
            const std::vector<int> paths = augmentationPaths(mesh, data, charts, static_cast<int>(c), worst, kPathsPerRound);
            if (paths.empty()) {
                return;
            }

            std::vector<int> seams = charts.charts[c].seamHalfEdges;
            HalfEdgeChart    before = std::move(atlas[c]);
            addCutEdges(mesh, charts.charts[c], paths);
            refresh(c);
            heal(c);
            if (pairs[c] == 0 && chartStretch(data, atlas[c]).l2 < s.l2) {
                ++st.cuts;
                continue;
            }
            charts.charts[c].seamHalfEdges = std::move(seams);
            atlas[c] = std::move(before);
            pairs[c] = 0;
            return;
        }
    }

    /// One round over the charts folded at its start; returns the number of splits.
    int round() {
        const std::size_t count = charts.charts.size();
        int               done = 0;
        for (std::size_t c = 0; c < count; ++c) {
            if (pairs[c] > 0 && charts.charts[c].faces.size() >= 2 && repair(c)) {
                ++done;
            }
        }
        return done;
    }
};

} // namespace

int selfOverlapPairs(const HalfEdgeChart &ac) {
    return census(ac).pairs;
}


/// Gamma >= gamma of one triangle's uv -> 3D map (Sander 2001), unnormalised; false for a flipped or
/// degenerate UV triangle.  @p p are the three 3D corners, @p t the three UV corners.
static bool triangleSingularValues(const std::array<math::Vec3f, 3> &p, const std::array<math::Vec2f, 3> &t,
                                   double &big, double &small) {
    const double s1 = t[0].x, t1 = t[0].y, s2 = t[1].x, t2 = t[1].y, s3 = t[2].x, t3 = t[2].y;
    const double a2 = (s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1);
    if (!(a2 > 0.0)) {
        return false;
    }

    // Partial derivatives of the 3D position along u and v; a2 is twice the UV area.
    double ss = 0.0, st = 0.0, tt = 0.0;
    for (int k = 0; k < 3; ++k) {
        const auto   i = static_cast<std::size_t>(k);
        const double su = (p[0][i] * (t2 - t3) + p[1][i] * (t3 - t1) + p[2][i] * (t1 - t2)) / a2;
        const double sv = (p[0][i] * (s3 - s2) + p[1][i] * (s1 - s3) + p[2][i] * (s2 - s1)) / a2;
        ss += su * su;
        st += su * sv;
        tt += sv * sv;
    }
    const double r = std::sqrt((ss - tt) * (ss - tt) + 4.0 * st * st);
    big = std::sqrt((ss + tt + r) / 2.0);
    small = std::sqrt(std::max(0.0, (ss + tt - r) / 2.0));
    return true;
}


ChartStretch chartStretch(const MeshData &data, const HalfEdgeChart &ac, std::vector<double> *perFace) {
    const std::size_t   n = ac.faces.size();
    std::vector<double> gamma(n, -1.0);
    double              sumA3 = 0.0;
    double              sumAuv = 0.0;
    double              sumL2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const auto                 f = static_cast<std::size_t>(ac.faces[i]);
        std::array<math::Vec3f, 3> p {};
        std::array<math::Vec2f, 3> t {};
        for (std::size_t k = 0; k < 3; ++k) {
            p[k] = data.positions[static_cast<std::size_t>(data.indices[3 * f + k])];
            t[k] = ac.uv[static_cast<std::size_t>(localCorner(ac, static_cast<int>(i), static_cast<int>(k)))];
        }

        double big = 0.0, small = 0.0;
        if (!triangleSingularValues(p, t, big, small)) {
            continue;
        }
        const double a3 = 0.5 * std::sqrt(static_cast<double>((p[1] - p[0]).cross(p[2] - p[0]).lengthSquared()));
        const double auv = 0.5 * (static_cast<double>(t[1].x - t[0].x) * (t[2].y - t[0].y) -
                                  static_cast<double>(t[2].x - t[0].x) * (t[1].y - t[0].y));
        gamma[i] = big;
        sumA3 += a3;
        sumAuv += auv;
        sumL2 += 0.5 * (big * big + small * small) * a3;
    }

    // Normalise by the chart's own area ratio, as packAtlas scales it.
    ChartStretch out;
    if (!(sumA3 > 0.0) || !(sumAuv > 0.0)) {
        return out;
    }
    const double k = std::sqrt(sumAuv / sumA3);
    out.l2 = k * std::sqrt(sumL2 / sumA3);
    for (std::size_t i = 0; i < n; ++i) {
        if (gamma[i] < 0.0) {
            continue;
        }
        gamma[i] *= k;
        if (out.worstFace < 0 || gamma[i] > out.linf) {
            out.linf = gamma[i];
            out.worstFace = static_cast<int>(i);
        }
    }
    if (perFace != nullptr) {
        *perFace = std::move(gamma);
    }
    return out;
}


/// Cut every closed chart open (see `closedChartCut`); returns how many were opened.
static int openClosedCharts(const HalfEdgeMesh &mesh, const MeshData &data, ChartSet &charts) {
    int opened = 0;
    for (std::size_t c = 0; c < charts.charts.size(); ++c) {
        if (!isClosedChart(mesh, charts, static_cast<int>(c))) {
            continue;
        }
        const std::vector<int> cut = closedChartCut(mesh, data, charts, static_cast<int>(c));
        addCutEdges(mesh, charts.charts[c], cut);
        opened += cut.empty() ? 0 : 1;
    }
    return opened;
}


std::vector<HalfEdgeChart> flattenCharts(const HalfEdgeMesh &mesh, const MeshData &data,
                                         ChartSet &charts, const InjectivityOptions &opts,
                                         InjectivityStats *stats) {
    InjectivityStats st;
    Guard            guard { mesh, data, charts, opts, st, {}, {} };
    guard.atlas.reserve(charts.charts.size());

    // A chart without a boundary cannot be flattened or packed: open it before anything flattens.
    if (opts.flattener != Flattener::LscmRaw) {
        st.closedOpened = openClosedCharts(mesh, data, charts);
    }

    // First flatten everything, so the "before" statistic reflects plain LSCM.
    for (std::size_t c = 0; c < charts.charts.size(); ++c) {
        guard.append(c);
    }
    st.selfOverlappingBefore = guard.foldedCount();

    // Heal in place first; only what stays folded costs a split.
    for (std::size_t c = 0; c < charts.charts.size(); ++c) {
        guard.heal(c);
    }

    // Stretch bound (Spec 20): cut over-bound charts open; folds the cuts leave are split below.
    const bool bounded = opts.maxStretchL2 > 0.0 || opts.maxStretchLinf > 0.0;
    for (std::size_t c = 0; bounded && opts.flattener != Flattener::LscmRaw && c < charts.charts.size(); ++c) {
        guard.cutToBound(c);
    }

    // LscmRaw detects only, so it gets no split rounds.
    const int maxRounds = opts.flattener == Flattener::LscmRaw ? 0 : opts.maxRounds;
    while (st.rounds < maxRounds) {
        const int splitNow = guard.round();
        if (splitNow == 0) {
            break; // nothing left that a split can fix
        }
        st.splits += splitNow;
        ++st.rounds;
    }

    st.selfOverlappingAfter = guard.foldedCount();
    st.nonInjectiveAfter = static_cast<int>(std::count_if(
        guard.atlas.begin(), guard.atlas.end(), [](const HalfEdgeChart &a) { return !a.injective(); }));

    for (const HalfEdgeChart &a : guard.atlas) {
        const ChartStretch cs = chartStretch(data, a);
        st.maxChartL2 = std::max(st.maxChartL2, cs.l2);
        st.maxChartLinf = std::max(st.maxChartLinf, cs.linf);
    }

    if (stats != nullptr) {
        *stats = st;
    }
    return std::move(guard.atlas); // a member of Guard, so no implicit move
}

} // namespace decimate::chart
