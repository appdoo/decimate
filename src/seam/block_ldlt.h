// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file block_ldlt.h
/// @brief Internal: the erase system as a symmetric BxB-block matrix, and its block LDL^T.
///
/// Plain data and plain doubles (no Eigen in the interface).  B = 1: the object-space system, one scalar per
/// unknown texel, the three channels are three right-hand sides.  B = 3: the framed system couples the three
/// channels of a texel through its frame, so K is a graph of 3x3 blocks, one node per unknown texel.  Either
/// way a node carries 3 right-hand-side / solution values (`3*u + channel`).

#include <cstddef>
#include <memory>
#include <span>
#include <utility>
#include <vector>

namespace decimate::seam::detail {

/// @brief A symmetric block matrix, lower triangle only, every block BxB row-major (B*B doubles).
///
/// Column j lists its blocks (i, j), i > j, sorted by i; the diagonal blocks live apart.  The pattern is fixed at
/// construction (@ref BlockPattern), so assembling a term is a slot lookup and B*B additions, never an insertion.
template <int B>
struct BlockMatrix {
    int                 n = 0; ///< Block rows = block columns.
    std::vector<int>    start; ///< n + 1 offsets into `row`/`off`.
    std::vector<int>    row;   ///< Block row of each off-diagonal block.
    std::vector<double> diag;  ///< n * B*B: the diagonal blocks.
    std::vector<double> off;   ///< B*B per off-diagonal block.

    /// @brief The B*B values of block (i, j) for i >= j, or nullptr for i < j (the mirror holds it).
    /// @throws std::logic_error if the block is below the diagonal but not in the pattern.
    [[nodiscard]] double *slot(int i, int j);

    [[nodiscard]] const double *slot(int i, int j) const;
};


/// @brief Collects the block graph of a system before any value is known.
class BlockPattern {
public:
    explicit BlockPattern(int n);

    /// @brief Declare that blocks (a, b) and (b, a) may become non-zero.  a == b, or a negative index (a texel
    ///        fixed to its input, no unknown), is ignored.
    void couple(int a, int b);

    /// @brief The zero-valued matrix with exactly the declared pattern.
    template <int B>
    [[nodiscard]] BlockMatrix<B> build();

private:
    int                              n_;
    std::vector<std::pair<int, int>> pairs_; ///< (column, row), column < row.
};


/// @brief LDL^T of a @ref BlockMatrix: BxB blocks, up-looking, AMD + elimination-tree postorder.
///
/// The constructor is the symbolic phase (ordering, elimination tree, column counts); it then allocates ONE arena
/// that holds every number and index the numeric phase touches, so `factorize` and `solve` allocate nothing and a
/// block is B*B plain doubles at a computable offset.  Layout choice (measured, Spec 31): flat sections of
/// array-of-structs blocks beat per-column interleaving and padded 3x4 rows, because the step cost is the 3x3
/// products themselves, not the pointer chasing.
template <int B>
class BlockLdlt {
public:
    /// @param pattern the block graph of the system to factor (values are not read).
    explicit BlockLdlt(const BlockMatrix<B> &pattern);

    /// @brief Numeric factorization of @p k, which must have the pattern given to the constructor.
    /// @return false if a diagonal block is not positive definite; the caller retries once with a ridge.
    [[nodiscard]] bool factorize(const BlockMatrix<B> &k);

    /// @brief Solve K x = b after a successful `factorize`.
    /// @param b 3n values, row `3*u + channel`.
    /// @param x 3n values, same numbering; may not alias @p b.
    void solve(std::span<const double> b, std::span<double> x);

    /// @brief Bytes of the arena (receipt).
    [[nodiscard]] std::size_t arenaBytes() const {
        return arenaBytes_;
    }

private:
    void scatter(const BlockMatrix<B> &k);
    int  gatherRow(int step); ///< Row @p step's pattern into `stack_[top..n)`, topologically sorted; returns top.
    void eliminate(int i, int step, double *d);
    void forwardSubstitute();
    void scaleByDiagonal();
    void backSubstitute();

    int                          n_ = 0;
    std::size_t                  arenaBytes_ = 0;
    std::unique_ptr<std::byte[]> arena_;

    // Integer sections of the arena.
    std::span<int> perm_;   ///< n: elimination position -> block node.
    std::span<int> parent_; ///< n: elimination tree.
    std::span<int> fill_;   ///< n: blocks of column i written so far.
    std::span<int> flag_;   ///< n: last step that visited the column.
    std::span<int> stack_;  ///< n: the pattern of the current row.
    std::span<int> lp_;     ///< n + 1: column starts in `li_`/`lx_`.
    std::span<int> li_;     ///< L block rows.
    std::span<int> ap_;     ///< n + 1: column starts of the permuted upper A.
    std::span<int> ai_;     ///< rows of the permuted upper A.
    std::span<int> aslot_;  ///< per lower block of the input: 2 * position in `aval_` + transpose flag.

    // Double sections of the arena; a block is B*B row-major doubles.
    std::span<double> lx_;   ///< L blocks, column by column.
    std::span<double> dc_;   ///< n diagonal blocks: assembled D, then its Cholesky factor.
    std::span<double> y_;    ///< n scratch blocks of the current row.
    std::span<double> x_;    ///< 3n right-hand side / solution (3 per node), in elimination order.
    std::span<double> aval_; ///< permuted upper A blocks.
};

extern template struct BlockMatrix<1>;
extern template BlockMatrix<1> BlockPattern::build<1>();
extern template BlockMatrix<3> BlockPattern::build<3>();
extern template struct BlockMatrix<3>;
extern template class BlockLdlt<1>;
extern template class BlockLdlt<3>;

} // namespace decimate::seam::detail
