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
class Assembler {
public:
    explicit Assembler(EraseSystem &sys)
        : sys_(sys) {
    }

    void add(int row, int col, double value) {
        sys_.entries.push_back({ row, col, value });
    }

    void addRhs(int row, std::size_t channel, double value) {
        sys_.rhs[kChannels * static_cast<std::size_t>(row) + channel] += value;
    }

private:
    EraseSystem &sys_;
};


/// Weight of one seam energy and the total seam length it is normalized by.
struct SeamWeight {
    double weight = 0.0;
    double total = 0.0;
};


/// Grid accessor for the Dirichlet loops.
struct Lattice {
    const EraseInput  &in;
    const EraseSystem &sys;

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
static void addPieceQuadratic(Assembler &as, const EraseSystem &sys, const SeamPiece &piece,
                              double scale) {
    // Exact integration: the integrand is a polynomial of degree <= 4 in gamma.
    Matrix8         q {};
    const GaussRule rule = gaussOn(piece.a, piece.b);
    for (std::size_t g = 0; g < 3; ++g) {
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


/// Add the weighted quadratic of every piece, each scaled by its edge's share of the seam length.
static void addSeamTerms(Assembler &as, const EraseSystem &sys, std::span<const SeamPiece> pieces,
                         const SeamWeight &w) {
    for (const SeamPiece &piece : pieces) {
        addPieceQuadratic(as, sys, piece, w.weight * piece.length3D / w.total);
    }
}


/// w_change/|C| (x - y)^2 on trusted unknowns; |C| counts every non-free texel of the grid.
static void addAnchor(Assembler &as, const EraseSystem &sys, const EraseInput &in) {
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
static double edgeMass(const Lattice &lat, Point at, Point step) {
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


/// m (x_b - x_a)^2 - 2 s m (x_b - x_a)(y_b - y_a): the masked Dirichlet term of one edge.
static void addSmoothEdge(Assembler &as, const Lattice &lat, std::size_t ta, std::size_t tb,
                          double mass) {
    const int ua = lat.sys.unknownOf[ta];
    const int ub = lat.sys.unknownOf[tb];
    as.add(ua, ua, mass);
    as.add(ub, ub, mass);
    as.add(ua, ub, -mass);
    as.add(ub, ua, -mass);

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


/// The two smoothness edges leaving lattice point (x,y): towards +x and towards +y.
static void addSmoothEdgesFrom(Assembler &as, const Lattice &lat, int x, int y) {
    for (const auto [dx, dy] : { std::pair { 1, 0 }, std::pair { 0, 1 } }) {
        if (!lat.isUnknown(x + dx, y + dy)) {
            continue;
        }
        const double mass = lat.in.weights.smooth * edgeMass(lat, { x, y }, { dx, dy });
        if (mass > 0.0) {
            addSmoothEdge(as, lat, lat.texel(x, y), lat.texel(x + dx, y + dy), mass);
        }
    }
}


/// The smoothness edges of every unknown texel in row @p y.
static void addSmoothRow(Assembler &as, const Lattice &lat, int y) {
    for (int x = 0; x < lat.in.grid.width; ++x) {
        if (lat.isUnknown(x, y)) {
            addSmoothEdgesFrom(as, lat, x, y);
        }
    }
}


/// The masked Dirichlet energy (w_grad) over the whole grid.
static void addSmoothness(Assembler &as, const EraseSystem &sys, const EraseInput &in) {
    const Lattice lat { in, sys };
    for (int y = 0; y < in.grid.height; ++y) {
        addSmoothRow(as, lat, y);
    }
}


/// Number the unknown texels in row-major order.
static void numberUnknowns(EraseSystem &sys, const EraseDomain &domain) {
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
static void pinDecoupled(EraseSystem &sys, const EraseInput &in) {
    std::vector<double> diag(static_cast<std::size_t>(sys.n), 0.0);
    for (const SparseEntry &e : sys.entries) {
        if (e.row == e.col) {
            diag[static_cast<std::size_t>(e.row)] += e.value;
        }
    }

    // A zero row would make K singular; fixing the unknown to its input value is the neutral choice.
    for (int u = 0; u < sys.n; ++u) {
        if (diag[static_cast<std::size_t>(u)] != 0.0) {
            continue;
        }
        const auto texel = static_cast<std::size_t>(sys.texelOf[static_cast<std::size_t>(u)]);
        sys.entries.push_back({ u, u, 1.0 });
        for (std::size_t ch = 0; ch < kChannels; ++ch) {
            sys.rhs[kChannels * static_cast<std::size_t>(u) + ch] =
                in.texels[kChannels * texel + ch];
        }
        ++sys.pinned;
    }
}


EraseSystem buildEraseSystem(const EraseInput &in) {
    EraseSystem sys;
    numberUnknowns(sys, in.domain);

    // Energy terms: seam mismatch, C1 continuity, anchor to input, smoothness.
    Assembler as(sys);
    addSeamTerms(as, sys, in.valuePieces, { in.weights.seam, in.totalLength });
    addSeamTerms(as, sys, in.gradPieces, { in.weights.gradient, in.totalLength });
    addAnchor(as, sys, in);
    addSmoothness(as, sys, in);

    // Last, because it needs the finished diagonal.
    pinDecoupled(sys, in);
    return sys;
}

} // namespace decimate::seam::detail
