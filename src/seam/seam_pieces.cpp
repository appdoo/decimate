// This file is part of the decimate library; see decimate.h for version/license details.
/// @file seam_pieces.cpp
/// @brief Exact piecewise-polynomial seam mismatch (paper eq. 3-8) -- see seam_pieces.h.

#include "seam_pieces.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace decimate::seam::detail {

namespace {

constexpr double kMinPiece = 1e-12; // intervals shorter than this carry no energy


/// s(gamma) = s0 + ds * gamma: one axis of a bilinear cell, relative to the cell origin.
struct AxisLine {
    double s0 = 0.0;
    double ds = 0.0;
};


/// Local (s,t) of a segment inside one cell, as linear functions of gamma.
struct LocalLine {
    AxisLine s;
    AxisLine t;
};


/// Coefficients over the cell's texels (p00, p10, p01, p11) of a degree-2 polynomial.
struct CellPolynomial {
    std::array<double, 4> c2 {};
    std::array<double, 4> c1 {};
    std::array<double, 4> c0 {};
};


/// Both sides of one seam pair, in gamma-aligned texel segments.
struct PairSegments {
    Segment                  a;
    Segment                  b;
    math::Vec2d              nA; // inward normals (Gradient pieces only)
    math::Vec2d              nB;
    std::array<EdgeFrame, 2> frame; // tangent frames of A and B (identity unless framed)
};

} // namespace


/// Parameters in (0,1) where the axis c0 -> c1 crosses an integer texel line in [0,maxCoord].
static void addAxisCrossings(double c0, double c1, int maxCoord, std::vector<double> &out) {
    const double d = c1 - c0;

    // A constant axis never crosses a grid line; a non-finite one has no meaningful crossing.
    if (d == 0.0 || !std::isfinite(c0) || !std::isfinite(c1)) {
        return;
    }

    // Clamping the range bounds the loop when a UV lies far outside the texture.
    const double lo = std::clamp(std::min(c0, c1), -1.0, static_cast<double>(maxCoord) + 1.0);
    const double hi = std::clamp(std::max(c0, c1), -1.0, static_cast<double>(maxCoord) + 1.0);
    const int    kLo = std::max(0, static_cast<int>(std::ceil(lo)));
    const int    kHi = std::min(maxCoord, static_cast<int>(std::floor(hi)));
    for (int k = kLo; k <= kHi; ++k) {
        const double g = (static_cast<double>(k) - c0) / d;
        if (g > 0.0 && g < 1.0) {
            out.push_back(g);
        }
    }
}


/// Sorted split parameters in [0,1] (0 and 1 included) for the given segments.
static std::vector<double> splitParameters(std::span<const Segment> segs, TexelGrid grid) {
    std::vector<double> g { 0.0, 1.0 };
    for (const Segment &s : segs) {
        addAxisCrossings(s.x0, s.x1, grid.width - 1, g);
        addAxisCrossings(s.y0, s.y1, grid.height - 1, g);
    }
    std::ranges::sort(g);
    return g;
}


/// Cell whose texels a lookup at @p gamma touches (clamp-to-edge, origin <= size-2).
static Cell cellAt(const Segment &seg, double gamma, TexelGrid grid) {
    const double x = std::clamp(seg.x0 + gamma * (seg.x1 - seg.x0), 0.0, grid.width - 1.0);
    const double y = std::clamp(seg.y0 + gamma * (seg.y1 - seg.y0), 0.0, grid.height - 1.0);
    return { std::min(static_cast<int>(std::floor(x)), grid.width - 2),
             std::min(static_cast<int>(std::floor(y)), grid.height - 2) };
}


/// One axis inside a piece: linear while inside [0,maxCoord], constant once clamped.
static AxisLine axisLine(double c0, double c1, double mid, int maxCoord, int origin) {
    // The piece is clamped as a whole, so its midpoint decides which regime it is in.
    const double cm = c0 + mid * (c1 - c0);
    if (cm < 0.0) {
        return { -static_cast<double>(origin), 0.0 };
    }
    if (cm > static_cast<double>(maxCoord)) {
        return { static_cast<double>(maxCoord - origin), 0.0 };
    }
    return { c0 - static_cast<double>(origin), c1 - c0 };
}


/// The segment's local cell coordinates (s,t) over the piece around @p mid.
static LocalLine localLine(const Segment &seg, double mid, Cell cell, TexelGrid grid) {
    return { axisLine(seg.x0, seg.x1, mid, grid.width - 1, cell.px),
             axisLine(seg.y0, seg.y1, mid, grid.height - 1, cell.py) };
}


/// Bilinear value (1-s)(1-t)p00 + s(1-t)p10 + (1-s)t p01 + st p11 as a polynomial in gamma.
static CellPolynomial valuePolynomial(const LocalLine &l) {
    const double   e = l.s.s0 * l.t.ds + l.t.s0 * l.s.ds; // d/dgamma of s*t at gamma = 0
    const double   st = l.s.s0 * l.t.s0;
    const double   q = l.s.ds * l.t.ds;
    CellPolynomial c;
    c.c2 = { q, -q, -q, q };
    c.c1 = { e - l.s.ds - l.t.ds, l.s.ds - e, l.t.ds - e, e };
    c.c0 = { 1.0 - l.s.s0 - l.t.s0 + st, l.s.s0 - st, l.t.s0 - st, st };
    return c;
}


/// Derivative of the bilinear along direction @p g (in cell axes) as a polynomial in gamma.
static CellPolynomial gradientPolynomial(const LocalLine &l, math::Vec2d g) {
    const double   slope = g.x * l.t.ds + g.y * l.s.ds;
    const double   base = g.x * l.t.s0 + g.y * l.s.s0;
    CellPolynomial c;
    c.c1 = { slope, -slope, -slope, slope };
    c.c0 = { base - g.x - g.y, g.x - base, g.y - base, base };
    return c;
}


/// UV of corner @p corner (taken mod 3) of the face @p he lies in.
static math::Vec2f cornerUv(const TexturedMesh &tm, HalfEdgeRef he, int corner) {
    const std::size_t slot =
        3 * static_cast<std::size_t>(he.face) + static_cast<std::size_t>(corner % 3);
    return tm.uvs[static_cast<std::size_t>(tm.faceU[slot])];
}


/// Unit UV-space direction perpendicular to @p he pointing into its own triangle.
static math::Vec2d inwardNormal(const TexturedMesh &tm, HalfEdgeRef he) {
    const math::Vec2f p0 = cornerUv(tm, he, he.corner);
    const math::Vec2f p1 = cornerUv(tm, he, he.corner + 1);
    const math::Vec2f p2 = cornerUv(tm, he, he.corner + 2);
    const double      vx = static_cast<double>(p1.x) - static_cast<double>(p0.x);
    const double      vy = static_cast<double>(p1.y) - static_cast<double>(p0.y);
    const double      area2 = (static_cast<double>(p1.x) - static_cast<double>(p0.x)) *
                                  (static_cast<double>(p2.y) - static_cast<double>(p0.y)) -
                              (static_cast<double>(p2.x) - static_cast<double>(p0.x)) *
                                  (static_cast<double>(p1.y) - static_cast<double>(p0.y));

    // The triangle's winding decides which perpendicular points inwards.
    const math::Vec2d perp = area2 > 0.0 ? math::Vec2d(-vy, vx) : math::Vec2d(vy, -vx);
    const double      len = std::hypot(perp.x, perp.y);
    return len > 1e-8 ? math::Vec2d(perp.x / len, perp.y / len) : perp;
}


/// Texel indices of a cell in the order (p00, p10, p01, p11).
static std::array<int, 4> cellTexels(Cell c, TexelGrid grid) {
    const int base = c.py * grid.width + c.px;
    return { base, base + 1, base + grid.width, base + grid.width + 1 };
}


/// T and N of @p he between its corner @p c0 (gamma 0) and @p c1 (gamma 1).
static EdgeFrame edgeFrame(const TexturedMesh &tm, const FrameField &fr, HalfEdgeRef he, int c0, int c1) {
    const std::size_t f = 3 * static_cast<std::size_t>(he.face);
    auto              tan = [&](int c) {
        return math::Vec3d(fr.tangents[f + c % 3].x, fr.tangents[f + c % 3].y, fr.tangents[f + c % 3].z);
    };
    auto nor = [&](int c) {
        return math::Vec3d(fr.normals[static_cast<std::size_t>(tm.faceV[f + c % 3])]);
    };
    return { tan(c0), tan(c1), nor(c0), nor(c1), fr.tangents[f].w < 0.0f ? -1.0 : 1.0 };
}


/// Both half-edges of @p pair as texel segments, B reversed so gamma addresses one 3D point.
static PairSegments pairSegments(const TexturedMesh &tm, const SeamPair &pair, TexelGrid grid,
                                 const FrameField *frames) {
    PairSegments ps;
    ps.a = halfEdgeSegment(tm, pair.a, grid);
    const Segment rb = halfEdgeSegment(tm, pair.b, grid);
    ps.b = { rb.x1, rb.y1, rb.x0, rb.y0 }; // B runs end-to-start: same 3D point as A
    ps.nA = inwardNormal(tm, pair.a);
    ps.nB = inwardNormal(tm, pair.b);
    if (frames != nullptr) {
        ps.frame = { edgeFrame(tm, *frames, pair.a, pair.a.corner, pair.a.corner + 1),
                     edgeFrame(tm, *frames, pair.b, pair.b.corner + 1, pair.b.corner) };
    }
    return ps;
}


/// The piece over gamma in [a,b]: both sides' cells and polynomial coefficients.
static SeamPiece makePiece(const PairSegments &ps, double a, double b, double length,
                           TexelGrid grid, PieceKind kind) {
    // Between two split parameters neither side changes cell, so the midpoint identifies them.
    const double         mid = 0.5 * (a + b);
    const Cell           cellA = cellAt(ps.a, mid, grid);
    const Cell           cellB = cellAt(ps.b, mid, grid);
    const LocalLine      la = localLine(ps.a, mid, cellA, grid);
    const LocalLine      lb = localLine(ps.b, mid, cellB, grid);
    const bool           value = kind == PieceKind::Value;
    const CellPolynomial pa = value ? valuePolynomial(la) : gradientPolynomial(la, ps.nA);
    const CellPolynomial pb = value ? valuePolynomial(lb) : gradientPolynomial(lb, ps.nB);
    const double         sign = value ? -1.0 : 1.0; // mismatch = A - B; C1 residual = dA + dB

    SeamPiece piece;
    piece.a = a;
    piece.b = b;
    piece.length3D = length;
    piece.frame = ps.frame;

    // Side A fills slots 0..3, side B slots 4..7 (pre-multiplied by the sign).
    const auto ta = cellTexels(cellA, grid);
    const auto tb = cellTexels(cellB, grid);
    for (std::size_t i = 0; i < 4; ++i) {
        piece.texel[i] = ta[i];
        piece.texel[i + 4] = tb[i];
        piece.c2[i] = pa.c2[i];
        piece.c1[i] = pa.c1[i];
        piece.c0[i] = pa.c0[i];
        piece.c2[i + 4] = sign * pb.c2[i];
        piece.c1[i + 4] = sign * pb.c1[i];
        piece.c0[i + 4] = sign * pb.c0[i];
    }
    return piece;
}


/// Split one seam pair at every texel-line crossing and append its pieces.
static void appendPairPieces(const TexturedMesh &tm, const SeamPair &pair, TexelGrid grid,
                             PieceKind kind, const FrameField *frames, std::vector<SeamPiece> &out) {
    const PairSegments           ps = pairSegments(tm, pair, grid, frames);
    const std::array<Segment, 2> both { ps.a, ps.b };
    const std::vector<double>    g = splitParameters(both, grid);
    for (std::size_t i = 0; i + 1 < g.size(); ++i) {
        if (g[i + 1] - g[i] > kMinPiece) {
            out.push_back(makePiece(ps, g[i], g[i + 1], pair.length3D, grid, kind));
        }
    }
}


Segment toTexelSegment(math::Vec2f uvStart, math::Vec2f uvEnd, TexelGrid grid) {
    const double w = grid.width;
    const double h = grid.height;
    return { static_cast<double>(uvStart.x) * w - 0.5, static_cast<double>(uvStart.y) * h - 0.5,
             static_cast<double>(uvEnd.x) * w - 0.5, static_cast<double>(uvEnd.y) * h - 0.5 };
}


Segment halfEdgeSegment(const TexturedMesh &tm, HalfEdgeRef he, TexelGrid grid) {
    return toTexelSegment(cornerUv(tm, he, he.corner), cornerUv(tm, he, he.corner + 1), grid);
}


std::vector<Cell> cellsAlong(const Segment &seg, TexelGrid grid) {
    const std::array<Segment, 1> one { seg };
    const std::vector<double>    g = splitParameters(one, grid);
    std::vector<Cell>            cells;
    for (std::size_t i = 0; i + 1 < g.size(); ++i) {
        if (g[i + 1] - g[i] > kMinPiece) {
            cells.push_back(cellAt(seg, 0.5 * (g[i] + g[i + 1]), grid));
        }
    }
    return cells;
}


std::vector<SeamPiece> seamPieces(const TexturedMesh &tm, const SeamTopology &topo, TexelGrid grid,
                                  PieceKind kind, const FrameField *frames) {
    std::vector<SeamPiece> pieces;
    for (const SeamPair &pair : topo.seams) {
        appendPairPieces(tm, pair, grid, kind, frames, pieces);
    }
    return pieces;
}


double totalSeamLength(const SeamTopology &topo) {
    double total = 0.0;
    for (const SeamPair &pair : topo.seams) {
        total += pair.length3D;
    }
    return total;
}


bool hasSeamLength(const SeamTopology &topo) {
    // `> 0` is false for a NaN total too.
    return totalSeamLength(topo) > 0.0;
}


std::array<double, 8> evaluatePiece(const SeamPiece &piece, double gamma) {
    std::array<double, 8> phi {};
    for (std::size_t i = 0; i < 8; ++i) {
        phi[i] = (piece.c2[i] * gamma + piece.c1[i]) * gamma + piece.c0[i]; // Horner
    }
    return phi;
}


GaussRule5 gauss5On(double a, double b) {
    constexpr std::array<double, 5> kNode { -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831,
                                            0.9061798459386640 };
    constexpr std::array<double, 5> kWeight { 0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                                              0.4786286704993665, 0.2369268850561891 };
    const double                    mid = 0.5 * (a + b);
    const double                    half = 0.5 * (b - a);
    GaussRule5                      rule;
    for (std::size_t i = 0; i < 5; ++i) {
        rule.gamma[i] = mid + half * kNode[i];
        rule.weight[i] = half * kWeight[i];
    }
    return rule;
}


/// Row-major F = [T | w N x T | N] of @p e at @p gamma.
static std::array<double, 9> frameMatrix(const EdgeFrame &e, double gamma) {
    const math::Vec3d t = e.t0 * (1.0 - gamma) + e.t1 * gamma;
    const math::Vec3d n = e.n0 * (1.0 - gamma) + e.n1 * gamma;
    const math::Vec3d b = n.cross(t) * e.w;
    return { t.x, b.x, n.x, t.y, b.y, n.y, t.z, b.z, n.z };
}


std::array<std::array<double, 9>, 8> frameBlocks(const SeamPiece &piece, double gamma) {
    const std::array<double, 8>                phi = evaluatePiece(piece, gamma);
    const std::array<std::array<double, 9>, 2> f { frameMatrix(piece.frame[0], gamma), frameMatrix(piece.frame[1], gamma) };
    std::array<std::array<double, 9>, 8>       blocks {};
    for (std::size_t i = 0; i < 8; ++i) {
        for (std::size_t k = 0; k < 9; ++k) {
            blocks[i][k] = phi[i] * f[i / 4][k];
        }
    }
    return blocks;
}


std::array<double, 3> integrateEnergy(std::span<const SeamPiece> pieces, double totalLength,
                                           std::span<const double> texels) {
    std::array<double, 3> acc {};
    for (const SeamPiece &piece : pieces) {
        const GaussRule5 rule = gauss5On(piece.a, piece.b);
        for (std::size_t q = 0; q < 5; ++q) {
            const auto            blocks = frameBlocks(piece, rule.gamma[q]);
            std::array<double, 3> r {};
            for (std::size_t i = 0; i < 8; ++i) {
                const std::size_t t = 3 * static_cast<std::size_t>(piece.texel[i]);
                for (std::size_t k = 0; k < 3; ++k) {
                    r[k] += blocks[i][3 * k] * texels[t] + blocks[i][3 * k + 1] * texels[t + 1] + blocks[i][3 * k + 2] * texels[t + 2];
                }
            }
            for (std::size_t k = 0; k < 3; ++k) {
                acc[k] += rule.weight[q] * piece.length3D * r[k] * r[k];
            }
        }
    }
    for (double &d : acc) {
        d /= totalLength;
    }
    return acc;
}

} // namespace decimate::seam::detail
