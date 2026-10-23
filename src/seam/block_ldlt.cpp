// This file is part of the decimate library; see decimate.h for version/license details.
/// @file block_ldlt.cpp
/// @brief Block pattern, and the BxB-block LDL^T of the erase system (B = 1 object space, B = 3 framed).

#include "block_ldlt.h"

#include <Eigen/OrderingMethods>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace decimate::seam::detail {

namespace {

/// The block graph in the elimination order: the upper triangle by column, the elimination tree, the column counts of L.
struct Permuted {
    std::vector<int> iperm;  ///< block node -> elimination position.
    std::vector<int> ap;     ///< n + 1 column starts of the upper triangle.
    std::vector<int> ai;     ///< its rows (elimination positions), unsorted.
    std::vector<int> aslot;  ///< per input lower block: 2 * position in `ai` + 1 if the block must be transposed.
    std::vector<int> parent; ///< elimination tree, -1 at roots.
    std::vector<int> count;  ///< blocks of each column of L.
};


/// Hands out typed, 8-aligned spans of one arena; with a null base it only measures.
class Carver {
public:
    explicit Carver(std::byte *base)
        : base_(base) {
    }

    template <class T>
    std::span<T> take(std::size_t count) {
        const std::size_t at = used_;
        used_ = (used_ + count * sizeof(T) + 7) / 8 * 8;
        return { base_ == nullptr ? nullptr : reinterpret_cast<T *>(base_ + at), count };
    }

    [[nodiscard]] std::size_t used() const {
        return used_;
    }

private:
    std::byte  *base_;
    std::size_t used_ = 0;
};

} // namespace


/// Approximate-minimum-degree order of the block graph: result[k] = the node eliminated k-th.
template <int B>
static std::vector<int> minimumDegree(const BlockMatrix<B> &m) {
    // AMD needs the diagonal in the pattern it is handed, or it degenerates to a dense factor.
    std::vector<int> start(static_cast<std::size_t>(m.n) + 1, 0);
    std::vector<int> row;
    row.reserve(m.row.size() + static_cast<std::size_t>(m.n));
    for (int j = 0; j < m.n; ++j) {
        row.push_back(j);
        row.insert(row.end(), m.row.begin() + m.start[static_cast<std::size_t>(j)], m.row.begin() + m.start[static_cast<std::size_t>(j) + 1]);
        start[static_cast<std::size_t>(j) + 1] = static_cast<int>(row.size());
    }

    const std::vector<double>         ones(row.size(), 1.0);
    const Eigen::SparseMatrix<double> pattern = Eigen::Map<const Eigen::SparseMatrix<double>>(
        m.n, m.n, static_cast<Eigen::Index>(row.size()), start.data(), row.data(), ones.data());
    Eigen::AMDOrdering<int>                                       amd;
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> perm;
    amd(pattern, perm);
    return { perm.indices().data(), perm.indices().data() + m.n };
}


/// The graph, elimination tree and column counts for the elimination order @p perm (perm[k] = node).
template <int B>
static Permuted permute(const BlockMatrix<B> &m, const std::vector<int> &perm) {
    const auto n = static_cast<std::size_t>(m.n);
    Permuted   p;
    p.iperm.assign(n, 0);
    for (std::size_t k = 0; k < n; ++k) {
        p.iperm[static_cast<std::size_t>(perm[k])] = static_cast<int>(k);
    }

    // The input holds each block once, below the diagonal; the factorization wants the upper triangle by column.
    p.ap.assign(n + 1, 0);
    for (std::size_t j = 0; j < n; ++j) {
        for (int s = m.start[j]; s < m.start[j + 1]; ++s) {
            const int a = p.iperm[static_cast<std::size_t>(m.row[static_cast<std::size_t>(s)])];
            ++p.ap[static_cast<std::size_t>(std::max(a, p.iperm[j])) + 1];
        }
    }
    for (std::size_t k = 0; k < n; ++k) {
        p.ap[k + 1] += p.ap[k];
    }
    p.ai.assign(m.row.size(), 0);
    p.aslot.assign(m.row.size(), 0);
    std::vector<int> next(p.ap.begin(), p.ap.end() - 1);
    for (std::size_t j = 0; j < n; ++j) {
        for (int s = m.start[j]; s < m.start[j + 1]; ++s) {
            const int a = p.iperm[static_cast<std::size_t>(m.row[static_cast<std::size_t>(s)])];
            const int b = p.iperm[j];
            const int at = next[static_cast<std::size_t>(std::max(a, b))]++;
            p.ai[static_cast<std::size_t>(at)] = std::min(a, b);
            p.aslot[static_cast<std::size_t>(s)] = 2 * at + (a > b ? 1 : 0);
        }
    }

    // Elimination tree and column counts by the reach of each row (Davis, LDL).
    p.parent.assign(n, -1);
    p.count.assign(n, 0);
    std::vector<int> seen(n, -1);
    for (std::size_t k = 0; k < n; ++k) {
        seen[k] = static_cast<int>(k);
        for (int s = p.ap[k]; s < p.ap[k + 1]; ++s) {
            for (auto i = static_cast<std::size_t>(p.ai[static_cast<std::size_t>(s)]); seen[i] != static_cast<int>(k); i = static_cast<std::size_t>(p.parent[i])) {
                if (p.parent[i] < 0) {
                    p.parent[i] = static_cast<int>(k);
                }
                ++p.count[i];
                seen[i] = static_cast<int>(k);
            }
        }
    }
    return p;
}


/// Children before parents, siblings in index order: the postorder of a forest.
static std::vector<int> postorder(const std::vector<int> &parent) {
    const auto       n = static_cast<int>(parent.size());
    std::vector<int> firstChild(parent.size(), -1);
    std::vector<int> sibling(parent.size(), -1);
    for (int j = n - 1; j >= 0; --j) {
        const int up = parent[static_cast<std::size_t>(j)];
        if (up >= 0) {
            sibling[static_cast<std::size_t>(j)] = firstChild[static_cast<std::size_t>(up)];
            firstChild[static_cast<std::size_t>(up)] = j;
        }
    }

    std::vector<int> order;
    std::vector<int> pending;
    order.reserve(parent.size());
    for (int root = 0; root < n; ++root) {
        if (parent[static_cast<std::size_t>(root)] >= 0) {
            continue;
        }
        pending.push_back(root);
        while (!pending.empty()) {
            const auto v = static_cast<std::size_t>(pending.back());
            if (firstChild[v] < 0) {
                order.push_back(static_cast<int>(v));
                pending.pop_back();
                continue;
            }
            pending.push_back(firstChild[v]);
            firstChild[v] = sibling[static_cast<std::size_t>(firstChild[v])];
        }
    }
    return order;
}

template <int B>
double *BlockMatrix<B>::slot(int i, int j) {
    return const_cast<double *>(std::as_const(*this).slot(i, j));
}


template <int B>
const double *BlockMatrix<B>::slot(int i, int j) const {
    if (i == j) {
        return diag.data() + B * B * static_cast<std::size_t>(i);
    }
    if (i < j) {
        return nullptr;
    }

    const int *first = row.data() + start[static_cast<std::size_t>(j)];
    const int *last = row.data() + start[static_cast<std::size_t>(j) + 1];
    const int *hit = std::lower_bound(first, last, i);
    if (hit == last || *hit != i) {
        throw std::logic_error("BlockMatrix: block outside the declared pattern");
    }
    return off.data() + B * B * static_cast<std::size_t>(hit - row.data());
}


BlockPattern::BlockPattern(int n)
    : n_(n) {
}


void BlockPattern::couple(int a, int b) {
    if (a != b && a >= 0 && b >= 0) {
        pairs_.emplace_back(std::min(a, b), std::max(a, b));
    }
}


template <int B>
BlockMatrix<B> BlockPattern::build() {
    BlockMatrix<B> m;
    m.n = n_;
    m.start.assign(static_cast<std::size_t>(n_) + 1, 0);

    // Counting sort by column, then sort + unique inside each (short) column.
    for (const auto &p : pairs_) {
        ++m.start[static_cast<std::size_t>(p.first) + 1];
    }
    for (int j = 0; j < n_; ++j) {
        m.start[static_cast<std::size_t>(j) + 1] += m.start[static_cast<std::size_t>(j)];
    }
    std::vector<int> fill(m.start.begin(), m.start.end() - 1);
    m.row.resize(pairs_.size());
    for (const auto &p : pairs_) {
        m.row[static_cast<std::size_t>(fill[static_cast<std::size_t>(p.first)]++)] = p.second;
    }

    int kept = 0;
    for (int j = 0; j < n_; ++j) {
        int *first = m.row.data() + m.start[static_cast<std::size_t>(j)];
        int *last = m.row.data() + m.start[static_cast<std::size_t>(j) + 1];
        std::sort(first, last);
        last = std::unique(first, last);
        const int count = static_cast<int>(last - first);
        std::copy(first, last, m.row.begin() + kept);
        kept += count;
        m.start[static_cast<std::size_t>(j)] = kept - count;
    }
    m.start[static_cast<std::size_t>(n_)] = kept;
    m.row.resize(static_cast<std::size_t>(kept));
    m.diag.assign(B * B * static_cast<std::size_t>(n_), 0.0);
    m.off.assign(B * B * static_cast<std::size_t>(kept), 0.0);
    pairs_.clear();
    return m;
}


/// y -= l * yi, all BxB row-major.
template <int B>
static void subtractProduct(double *y, const double *l, const double *yi) {
    for (int r = 0; r < B; ++r) {
        for (int c = 0; c < B; ++c) {
            double dot = 0.0;
            for (int k = 0; k < B; ++k) {
                dot += l[B * r + k] * yi[B * k + c];
            }
            y[B * r + c] -= dot;
        }
    }
}


/// out = (C C^T)^-1 in for one B-vector read and written with @p stride; @p c is the factor of `factorDiagonal`.
template <int B>
static void solveColumn(const double *c, const double *in, double *out, int stride) {
    double z[B];
    double w[B];
    for (int i = 0; i < B; ++i) {
        double v = in[i * stride];
        for (int k = 0; k < i; ++k) {
            v -= c[B * i + k] * z[k];
        }
        z[i] = v * c[B * i + i];
    }
    for (int i = B - 1; i >= 0; --i) {
        double v = z[i];
        for (int k = i + 1; k < B; ++k) {
            v -= c[B * k + i] * w[k];
        }
        w[i] = v * c[B * i + i];
    }
    for (int i = 0; i < B; ++i) {
        out[i * stride] = w[i];
    }
}


/// The Cholesky factor C of the symmetric BxB @p d (lower triangle read), stored with the diagonal inverted so the
/// solves multiply; false unless @p d is positive definite.
template <int B>
static bool factorDiagonal(const double *d, double *c) {
    double root[B];
    for (int i = 0; i < B; ++i) {
        for (int j = 0; j <= i; ++j) {
            double v = d[B * i + j];
            for (int k = 0; k < j; ++k) {
                v -= c[B * i + k] * c[B * j + k];
            }
            if (j < i) {
                c[B * i + j] = v / root[j];
                continue;
            }
            if (!(v > 0.0)) {
                return false;
            }
            root[i] = std::sqrt(v);
            c[B * i + i] = 1.0 / root[i];
        }
    }
    return true;
}


template <int B>
BlockLdlt<B>::BlockLdlt(const BlockMatrix<B> &pattern)
    : n_(pattern.n) {
    if (n_ == 0) {
        return;
    }

    // Symbolic: AMD, then the postorder of its elimination tree so a subtree's columns sit together.
    std::vector<int> order = minimumDegree(pattern);
    const Permuted   first = permute(pattern, order);
    const auto       post = postorder(first.parent);
    std::vector<int> composed(order.size());
    for (std::size_t k = 0; k < order.size(); ++k) {
        composed[k] = order[static_cast<std::size_t>(post[k])];
    }
    order = std::move(composed);
    const Permuted p = permute(pattern, order);

    std::vector<int> colStart(static_cast<std::size_t>(n_) + 1, 0);
    for (std::size_t k = 0; k < p.count.size(); ++k) {
        colStart[k + 1] = colStart[k] + p.count[k];
    }

    // One arena for everything the numeric phase touches, sized once.
    const auto n = static_cast<std::size_t>(n_);
    const auto nnzL = static_cast<std::size_t>(colStart[n]);
    const auto nnzA = p.ai.size();
    const auto carve = [&](Carver &c) {
        perm_ = c.take<int>(n);
        parent_ = c.take<int>(n);
        fill_ = c.take<int>(n);
        flag_ = c.take<int>(n);
        stack_ = c.take<int>(n);
        lp_ = c.take<int>(n + 1);
        li_ = c.take<int>(nnzL);
        ap_ = c.take<int>(n + 1);
        ai_ = c.take<int>(nnzA);
        aslot_ = c.take<int>(nnzA);
        lx_ = c.take<double>(B * B * nnzL);
        dc_ = c.take<double>(B * B * n);
        y_ = c.take<double>(B * B * n);
        x_ = c.take<double>(3 * n);
        aval_ = c.take<double>(B * B * nnzA);
    };
    Carver measure(nullptr);
    carve(measure);
    arenaBytes_ = measure.used();
    arena_ = std::make_unique_for_overwrite<std::byte[]>(arenaBytes_);
    Carver real(arena_.get());
    carve(real);

    std::copy(order.begin(), order.end(), perm_.begin());
    std::copy(p.parent.begin(), p.parent.end(), parent_.begin());
    std::copy(colStart.begin(), colStart.end(), lp_.begin());
    std::copy(p.ap.begin(), p.ap.end(), ap_.begin());
    std::copy(p.ai.begin(), p.ai.end(), ai_.begin());
    std::copy(p.aslot.begin(), p.aslot.end(), aslot_.begin());
    std::fill(y_.begin(), y_.end(), 0.0);
}


template <int B>
void BlockLdlt<B>::scatter(const BlockMatrix<B> &k) {
    // Diagonal blocks and upper-A blocks, moved to elimination order (a block below the new diagonal is transposed).
    for (std::size_t e = 0; e < perm_.size(); ++e) {
        std::copy_n(k.diag.data() + B * B * static_cast<std::size_t>(perm_[e]), B * B, dc_.data() + B * B * e);
    }
    for (std::size_t s = 0; s < aslot_.size(); ++s) {
        const double *src = k.off.data() + B * B * s;
        double       *dst = aval_.data() + B * B * static_cast<std::size_t>(aslot_[s] >> 1);
        for (int r = 0; r < B; ++r) {
            for (int c = 0; c < B; ++c) {
                dst[B * r + c] = (aslot_[s] & 1) != 0 ? src[B * c + r] : src[B * r + c];
            }
        }
    }
}


template <int B>
int BlockLdlt<B>::gatherRow(int step) {
    // Seed y with row `step` of A, and walk each entry up the elimination tree to the first column already seen.
    int top = n_;
    flag_[static_cast<std::size_t>(step)] = step;
    for (int s = ap_[static_cast<std::size_t>(step)]; s < ap_[static_cast<std::size_t>(step) + 1]; ++s) {
        auto i = static_cast<std::size_t>(ai_[static_cast<std::size_t>(s)]);
        std::copy_n(aval_.data() + B * B * static_cast<std::size_t>(s), B * B, y_.data() + B * B * i);

        int len = 0;
        for (; flag_[i] != step; i = static_cast<std::size_t>(parent_[i])) {
            stack_[static_cast<std::size_t>(len++)] = static_cast<int>(i);
            flag_[i] = step;
        }
        while (len > 0) {
            stack_[static_cast<std::size_t>(--top)] = stack_[static_cast<std::size_t>(--len)];
        }
    }
    return top;
}


template <int B>
void BlockLdlt<B>::eliminate(int i, int step, double *d) {
    const auto col = static_cast<std::size_t>(i);
    double    *yi = y_.data() + B * B * col;
    double     yCopy[B * B];
    std::copy(yi, yi + B * B, yCopy);
    std::fill(yi, yi + B * B, 0.0);

    // Push column i's finished blocks into the rows they touch.
    const std::size_t first = static_cast<std::size_t>(lp_[col]);
    const std::size_t end = first + static_cast<std::size_t>(fill_[col]);
    for (std::size_t p = first; p < end; ++p) {
        subtractProduct<B>(y_.data() + B * B * static_cast<std::size_t>(li_[p]), lx_.data() + B * B * p, yCopy);
    }

    // w = D_i^-1 y_i is L_(step,i)^T; the Schur update of D_step uses it too.
    double w[B * B];
    for (int j = 0; j < B; ++j) {
        solveColumn<B>(dc_.data() + B * B * col, yCopy + j, w + j, B);
    }
    for (int r = 0; r < B; ++r) {
        for (int c = 0; c < B; ++c) {
            double dot = 0.0;
            for (int k = 0; k < B; ++k) {
                dot += yCopy[B * k + r] * w[B * k + c];
            }
            d[B * r + c] -= dot;
        }
    }

    double *l = lx_.data() + B * B * end;
    for (int r = 0; r < B; ++r) {
        for (int c = 0; c < B; ++c) {
            l[B * r + c] = w[B * c + r];
        }
    }
    li_[end] = step;
    ++fill_[col];
}


template <int B>
bool BlockLdlt<B>::factorize(const BlockMatrix<B> &k) {
    scatter(k);
    std::fill(fill_.begin(), fill_.end(), 0);
    std::fill(flag_.begin(), flag_.end(), -1);

    for (int step = 0; step < n_; ++step) {
        int top = gatherRow(step);

        double *dStep = dc_.data() + B * B * static_cast<std::size_t>(step);
        double  d[B * B];
        std::copy(dStep, dStep + B * B, d);
        for (; top < n_; ++top) {
            eliminate(stack_[static_cast<std::size_t>(top)], step, d);
        }
        if (!factorDiagonal<B>(d, dStep)) {
            return false;
        }
    }
    return true;
}


// The three values of a node are 3 / B groups of B: one group of 3 (framed), or 3 scalar channels.
template <int B>
void BlockLdlt<B>::forwardSubstitute() {
    for (std::size_t i = 0; i < x_.size() / 3; ++i) {
        for (auto p = static_cast<std::size_t>(lp_[i]); p < static_cast<std::size_t>(lp_[i]) + static_cast<std::size_t>(fill_[i]); ++p) {
            const double *l = lx_.data() + B * B * p;
            for (int g = 0; g < 3; g += B) {
                const double *xi = x_.data() + 3 * i + g;
                double       *target = x_.data() + 3 * static_cast<std::size_t>(li_[p]) + g;
                for (int r = 0; r < B; ++r) {
                    double dot = 0.0;
                    for (int k = 0; k < B; ++k) {
                        dot += l[B * r + k] * xi[k];
                    }
                    target[r] -= dot;
                }
            }
        }
    }
}


template <int B>
void BlockLdlt<B>::scaleByDiagonal() {
    for (std::size_t i = 0; i < x_.size() / 3; ++i) {
        for (int g = 0; g < 3; g += B) {
            solveColumn<B>(dc_.data() + B * B * i, x_.data() + 3 * i + g, x_.data() + 3 * i + g, 1);
        }
    }
}


template <int B>
void BlockLdlt<B>::backSubstitute() {
    for (std::size_t i = x_.size() / 3; i-- > 0;) {
        for (auto p = static_cast<std::size_t>(lp_[i]); p < static_cast<std::size_t>(lp_[i]) + static_cast<std::size_t>(fill_[i]); ++p) {
            const double *l = lx_.data() + B * B * p;
            for (int g = 0; g < 3; g += B) {
                double       *xi = x_.data() + 3 * i + g;
                const double *above = x_.data() + 3 * static_cast<std::size_t>(li_[p]) + g;
                for (int c = 0; c < B; ++c) {
                    double dot = 0.0;
                    for (int k = 0; k < B; ++k) {
                        dot += l[B * k + c] * above[k];
                    }
                    xi[c] -= dot;
                }
            }
        }
    }
}


template <int B>
void BlockLdlt<B>::solve(std::span<const double> b, std::span<double> x) {
    for (std::size_t e = 0; e < perm_.size(); ++e) {
        std::copy_n(b.data() + 3 * static_cast<std::size_t>(perm_[e]), 3, x_.data() + 3 * e);
    }
    forwardSubstitute();
    scaleByDiagonal();
    backSubstitute();
    for (std::size_t e = 0; e < perm_.size(); ++e) {
        std::copy_n(x_.data() + 3 * e, 3, x.data() + 3 * static_cast<std::size_t>(perm_[e]));
    }
}


template struct BlockMatrix<1>;
template struct BlockMatrix<3>;
template BlockMatrix<1> BlockPattern::build<1>();
template BlockMatrix<3> BlockPattern::build<3>();
template class BlockLdlt<1>;
template class BlockLdlt<3>;

} // namespace decimate::seam::detail
