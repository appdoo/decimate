// This file is part of the decimate library; see decimate.h for version/license details.
/// @file erase_system.cpp
/// @brief Assembly of the erase normal equations (paper section 4.1).

#include "erase_system.h"

#include <array>
#include <cstddef>
#include <utility>

namespace decimate::seam::detail {

namespace {

constexpr std::size_t kChannels = 3;

using Matrix8 = std::array<std::array<double, 8>, 8>;


/// Accumulates symmetric terms into an EraseSystem.
template <int B>
class Assembler {
public:
    explicit Assembler(EraseSystem<B> &sys)
        : sys_(sys) {
    }

    /// A scalar term between unknowns @p row and @p col: once per channel of the block diagonal.
    void add(int row, int col, double value) {
        // The mirror block above the diagonal is not stored.
        double *block = sys_.blocks.slot(row, col);
        if (block == nullptr) {
            return;
        }
        for (int k = 0; k < B; ++k) {
            block[k * (B + 1)] += value;
        }
    }

    [[nodiscard]] BlockMatrix<B> &blocks() {
        return sys_.blocks;
    }

    void addRhs(int row, std::size_t channel, double value) {
        sys_.rhs[kChannels * static_cast<std::size_t>(row) + channel] += value;
    }

private:
    EraseSystem<B> &sys_;
};


/// Weight of one seam energy and the total seam length it is normalized by.
struct SeamWeight {
    double weight = 0.0;
    double total = 0.0;
};


/// Grid accessor for the Dirichlet loops.
template <int B>
struct Lattice {
    const EraseInput     &in;
    const EraseSystem<B> &sys;

    [[nodiscard]] bool isUnknown(int x, int y) const {
        if (x < 0 || y < 0 || x >= in.grid.width || y >= in.grid.height) {
            return false;
        }
        return in.domain.unknown[static_cast<std::size_t>(y) * in.grid.width + x] != 0;
    }

    [[nodiscard]] std::size_t texel(int x, int y) const {
        return static_cast<std::size_t>(y) * in.grid.width + x;
    }
};


/// A lattice point or step.
struct Point {
    int x = 0;
    int y = 0;
};

} // namespace


/// q += w * phi phi^T
static void addOuterProduct(Matrix8 &q, const std::array<double, 8> &phi, double w) {
    for (std::size_t i = 0; i < 8; ++i) {
        for (std::size_t j = 0; j < 8; ++j) {
            q[i][j] += w * phi[i] * phi[j];
        }
    }
}


/// scale * integral over the piece of phi phi^T, mapped to unknown ids.
template <int B>
static void addPieceQuadratic(Assembler<B> &as, const EraseSystem<B> &sys, const SeamPiece &piece,
                              double scale) {
    // Exact integration: the integrand is a polynomial of degree <= 4 in gamma.
    Matrix8          q {};
    const GaussRule5 rule = gauss5On(piece.a, piece.b);
    for (std::size_t g = 0; g < 5; ++g) {
        addOuterProduct(q, evaluatePiece(piece, rule.gamma[g]), scale * rule.weight[g]);
    }

    for (std::size_t i = 0; i < 8; ++i) {
        for (std::size_t j = 0; j < 8; ++j) {
            const int ui = sys.unknownOf[static_cast<std::size_t>(piece.texel[i])];
            const int uj = sys.unknownOf[static_cast<std::size_t>(piece.texel[j])];
            as.add(ui, uj, q[i][j]);
        }
    }
}


/// The 3x3 block of @p q at (@p r0, @p c0) added to the row-major @p block.
static void addBlock(double *block, const std::array<std::array<double, 24>, 24> &q, std::size_t r0, std::size_t c0) {
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            block[3 * r + c] += q[r0 + r][c0 + c];
        }
    }
}


/// scale * integral over the framed piece of Phi^T Phi: 24x24, mapped to the 3n numbering.
template <int B>
static void addFramedQuadratic(Assembler<B> &as, const EraseSystem<B> &sys, const SeamPiece &piece, double scale) {
    // Exact: the integrand has degree <= 8, so 5 Gauss nodes.
    std::array<std::array<double, 24>, 24> q {};
    const GaussRule5                       rule = gauss5On(piece.a, piece.b);
    for (std::size_t g = 0; g < 5; ++g) {
        const auto blocks = frameBlocks(piece, rule.gamma[g]);
        for (std::size_t i = 0; i < 24; ++i) {
            for (std::size_t j = 0; j < 24; ++j) {
                double dot = 0.0;
                for (std::size_t k = 0; k < 3; ++k) {
                    dot += blocks[i / 3][3 * k + i % 3] * blocks[j / 3][3 * k + j % 3];
                }
                q[i][j] += scale * rule.weight[g] * dot;
            }
        }
    }

    // Block (a, b) of the 24x24 goes to the lower triangle only; its mirror (b, a) is the same numbers.
    for (std::size_t a = 0; a < 8; ++a) {
        const int ua = sys.unknownOf[static_cast<std::size_t>(piece.texel[a])];
        for (std::size_t b = 0; b < 8; ++b) {
            double *block = as.blocks().slot(ua, sys.unknownOf[static_cast<std::size_t>(piece.texel[b])]);
            if (block != nullptr) {
                addBlock(block, q, 3 * a, 3 * b);
            }
        }
    }
}


/// Add the weighted quadratic of every piece, each scaled by its edge's share of the seam length.
template <int B>
static void addSeamTerms(Assembler<B> &as, const EraseSystem<B> &sys, std::span<const SeamPiece> pieces,
                         const SeamWeight &w, bool framed) {
    for (const SeamPiece &piece : pieces) {
        const double scale = w.weight * piece.length3D / w.total;
        if (framed) {
            addFramedQuadratic(as, sys, piece, scale);
        } else {
            addPieceQuadratic(as, sys, piece, scale);
        }
    }
}


/// w_change/|C| (x - y)^2 on trusted unknowns; |C| counts every non-free texel of the grid.
template <int B>
static void addAnchor(Assembler<B> &as, const EraseSystem<B> &sys, const EraseInput &in) {
    const std::size_t total = static_cast<std::size_t>(in.grid.width) * in.grid.height;
    const double      trusted = static_cast<double>(total) - in.domain.freeCount;
    const double      w = in.weights.change / trusted;
    for (int u = 0; u < sys.n; ++u) {
        const auto texel = static_cast<std::size_t>(sys.texelOf[static_cast<std::size_t>(u)]);

        // Free texels carry no trustworthy value, so they are not pulled towards the input.
        if (in.domain.free[texel] != 0) {
            continue;
        }

        as.add(u, u, w);
        for (std::size_t ch = 0; ch < kChannels; ++ch) {
            as.addRhs(u, ch, w * in.texels[kChannels * texel + ch]);
        }
    }
}


/// Mass of the grid edge from @p at by @p step: 1/8 per parallel neighbour edge fully in the mask.
template <int B>
static double edgeMass(const Lattice<B> &lat, Point at, Point step) {
    const auto [x, y] = at;
    const auto [dx, dy] = step;
    double mass = 0.0;
    for (const int side : { -1, 1 }) {
        const int sx = x + side * dy; // the parallel edge lies one cell sideways
        const int sy = y + side * dx;
        if (lat.isUnknown(sx, sy) && lat.isUnknown(sx + dx, sy + dy)) {
            mass += 0.125;
        }
    }
    return mass;
}


/// The quadratic of one edge between unknowns @p ua and @p ub.
template <int B>
static void addEdgeQuadratic(Assembler<B> &as, int ua, int ub, double mass) {
    as.add(ua, ua, mass);
    as.add(ub, ub, mass);
    as.add(ua, ub, -mass);
    as.add(ub, ua, -mass);
}


/// m (x_b - x_a)^2 - 2 s m (x_b - x_a)(y_b - y_a): the masked Dirichlet term of one edge.
template <int B>
static void addSmoothEdge(Assembler<B> &as, const Lattice<B> &lat, std::size_t ta, std::size_t tb,
                          double mass) {
    const int ua = lat.sys.unknownOf[ta];
    const int ub = lat.sys.unknownOf[tb];
    addEdgeQuadratic(as, ua, ub, mass);

    const bool trusted = lat.in.domain.free[ta] == 0 && lat.in.domain.free[tb] == 0;
    if (!trusted) {
        return; // target gradient is zero next to a free texel: pure smoothness
    }

    // Keep the input's own gradient across the edge, so smoothing does not flatten real detail.
    for (std::size_t ch = 0; ch < kChannels; ++ch) {
        const double dy = lat.in.texels[kChannels * tb + ch] - lat.in.texels[kChannels * ta + ch];
        as.addRhs(ua, ch, -mass * dy);
        as.addRhs(ub, ch, mass * dy);
    }
}


/// The smoothness edge from lattice point (x,y) by (dx,dy), if its far end is in the domain.
template <int B>
static void addSmoothEdgeFrom(Assembler<B> &as, const Lattice<B> &lat, Point at, Point step) {
    const auto [x, y] = at;
    const auto [dx, dy] = step;
    if (!lat.isUnknown(x + dx, y + dy)) {
        return;
    }
    const double mass = lat.in.weights.smooth * edgeMass(lat, at, step);
    if (mass > 0.0) {
        addSmoothEdge(as, lat, lat.texel(x, y), lat.texel(x + dx, y + dy), mass);
    }
}


/// The two smoothness edges leaving lattice point (x,y): towards +x and towards +y.
template <int B>
static void addSmoothEdgesFrom(Assembler<B> &as, const Lattice<B> &lat, int x, int y) {
    addSmoothEdgeFrom(as, lat, { x, y }, { 1, 0 });
    addSmoothEdgeFrom(as, lat, { x, y }, { 0, 1 });
}


/// The smoothness edges of every unknown texel in row @p y.
template <int B>
static void addSmoothRow(Assembler<B> &as, const Lattice<B> &lat, int y) {
    for (int x = 0; x < lat.in.grid.width; ++x) {
        if (lat.isUnknown(x, y)) {
            addSmoothEdgesFrom(as, lat, x, y);
        }
    }
}


/// The masked Dirichlet energy (w_grad) over the whole grid.
template <int B>
static void addSmoothness(Assembler<B> &as, const EraseSystem<B> &sys, const EraseInput &in) {
    const Lattice<B> lat { in, sys };
    for (int y = 0; y < in.grid.height; ++y) {
        addSmoothRow(as, lat, y);
    }
}


/// Number the unknown texels of the domain in row-major order.
template <int B>
static void numberUnknowns(EraseSystem<B> &sys, const EraseDomain &domain) {
    sys.unknownOf.assign(domain.unknown.size(), -1);
    for (std::size_t t = 0; t < domain.unknown.size(); ++t) {
        if (domain.unknown[t] != 0) {
            sys.unknownOf[t] = sys.n++;
            sys.texelOf.push_back(static_cast<int>(t));
        }
    }
    sys.rhs.assign(kChannels * static_cast<std::size_t>(sys.n), 0.0);
}


/// An unknown with an all-zero diagonal has no energy (zero-mass smoothness only): pin it.
template <int B>
static void pinDecoupled(EraseSystem<B> &sys, const EraseInput &in) {
    // The diagonal of channel 0 stands for the unknown: a framed unknown is pinned only if all its rows are empty.
    std::vector<double> diag(static_cast<std::size_t>(sys.n), 0.0);
    for (int u = 0; u < sys.n; ++u) {
        diag[static_cast<std::size_t>(u)] = sys.blocks.diag[B * B * static_cast<std::size_t>(u)];
    }

    // A zero row would make K singular; fixing the unknown to its input value is the neutral choice.
    Assembler<B> as(sys);
    for (int u = 0; u < sys.n; ++u) {
        if (diag[static_cast<std::size_t>(u)] != 0.0) {
            continue;
        }
        const auto texel = static_cast<std::size_t>(sys.texelOf[static_cast<std::size_t>(u)]);
        as.add(u, u, 1.0);
        for (std::size_t ch = 0; ch < kChannels; ++ch) {
            sys.rhs[kChannels * static_cast<std::size_t>(u) + ch] =
                in.texels[kChannels * texel + ch];
        }
        ++sys.pinned;
    }
}


/// Declare the unknowns a piece couples: all pairs of its 8 corner texels.
template <int B>
static void couplePiece(BlockPattern &pattern, const EraseSystem<B> &sys, const SeamPiece &piece) {
    for (std::size_t a = 0; a < 8; ++a) {
        const int ua = sys.unknownOf[static_cast<std::size_t>(piece.texel[a])];
        for (std::size_t b = a + 1; b < 8; ++b) {
            pattern.couple(ua, sys.unknownOf[static_cast<std::size_t>(piece.texel[b])]);
        }
    }
}


/// The block graph of the system before any value: piece footprints plus grid neighbours.
///
/// A superset of what the terms touch (a smoothness edge of zero mass still couples its ends here), which only
/// costs a zero block; a missing block would throw in `BlockMatrix::slot`.
template <int B>
static BlockMatrix<B> blockPattern(const EraseSystem<B> &sys, const EraseInput &in) {
    BlockPattern pattern(sys.n);
    for (const SeamPiece &piece : in.valuePieces) {
        couplePiece(pattern, sys, piece);
    }
    for (const SeamPiece &piece : in.gradPieces) {
        couplePiece(pattern, sys, piece);
    }

    const Lattice<B> lat { in, sys };
    for (int u = 0; u < sys.n; ++u) {
        const int t = sys.texelOf[static_cast<std::size_t>(u)];
        const int x = t % in.grid.width;
        const int y = t / in.grid.width;
        for (const Point step : { Point { 1, 0 }, Point { 0, 1 } }) {
            if (lat.isUnknown(x + step.x, y + step.y)) {
                pattern.couple(u, sys.unknownOf[lat.texel(x + step.x, y + step.y)]);
            }
        }
    }
    return pattern.template build<B>();
}


template <int B>
EraseSystem<B> buildEraseSystem(const EraseInput &in) {
    EraseSystem<B> sys;
    numberUnknowns(sys, in.domain);
    sys.blocks = blockPattern(sys, in);

    // Energy terms: seam mismatch (frame-coupled when B = 3), C1 continuity, anchor to input, smoothness.
    Assembler<B> as(sys);
    addSeamTerms(as, sys, in.valuePieces, { in.weights.seam, in.totalLength }, B == 3);
    addSeamTerms(as, sys, in.gradPieces, { in.weights.gradient, in.totalLength }, false);
    addAnchor(as, sys, in);
    addSmoothness(as, sys, in);

    // Last, because it needs the finished diagonal.
    pinDecoupled(sys, in);
    return sys;
}


template EraseSystem<1> buildEraseSystem<1>(const EraseInput &);
template EraseSystem<3> buildEraseSystem<3>(const EraseInput &);

} // namespace decimate::seam::detail
