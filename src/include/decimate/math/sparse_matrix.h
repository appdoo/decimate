// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file sparse_matrix.h
/// @brief Sparse matrix with Gauss-Seidel iterative solver — RAII port.
///
/// ## Design
///
/// This is a behavior-preserving port of legacy `SMat` from
/// `legacy2003/math/sparsematrix.h` (Spec 05, Beat 2a).
///
/// ### Changes from legacy
/// - The Win32 platform header include is removed (was line 17; unused there).
/// - All `printf` debug calls removed.  The *numeric* result of
///   `solveIterative()` is byte-identical: the solver loop, Gauss-Seidel
///   formula, convergence criterion (ε=1e-4, maxiter=10000) are unchanged.
/// - `rows` / `cols` head-pointer arrays: `new Sparse*[m/n]` →
///   `std::vector<Sparse*>` (eliminates manual allocation).
/// - `Sparse` nodes are still allocated individually with `new`/`delete`
///   in the linked-list structure; `clear()` / `~SMat()` RAII-cleans them.
///
/// ### What is NOT changed (antipattern 11)
/// - The Gauss-Seidel iteration formula and convergence condition.
/// - `get()` / `insert()` / `mul()` / `mulrowbycolumn()` semantics.
/// - The doubly-linked sparse node layout (`nextRow` / `nextCol`).
/// - `lrTriangulate` and `inverse` (kept as-is; printf removed only).
///
/// ### Solver characterization (T5, T8, AC7)
/// `solveIterative()` is Gauss-Seidel, NOT conjugate gradient.
/// The spec explicitly forbids upgrading — characterize as-is.
///
/// ## Note on `real` typedef
/// `typedef double real` is preserved so that the LSCM code that will
/// eventually consume this type (Spec 08) compiles without changes.

#include <cmath>
#include <cstring>
#include <utility>
#include <vector>

namespace decimate::math {

/// @brief Sparse matrix node (doubly-linked by row and column).
///
/// `SMat` owns all nodes it creates; they are freed by `clear()` / `~SMat()`.
struct Sparse {
    double  value = 0.0;       ///< Stored entry value.
    int     row = 0;           ///< Row index.
    int     col = 0;           ///< Column index.
    Sparse *nextRow = nullptr; ///< Next non-zero in the same column (row-ascending).
    Sparse *nextCol = nullptr; ///< Next non-zero in the same row   (col-ascending).
};


/// @brief Element value type (double) — preserved from legacy typedef.
using real = double;


/// @brief Sparse matrix — doubly-linked list of non-zero `Sparse` nodes.
///
/// Rows are singly-linked along `nextCol`; columns along `nextRow`.
/// `SMat` owns all nodes it allocates.  Copying is a deep copy.
///
/// Invariant: every node is reachable from exactly one row list and one column list, both
/// sorted ascending. Methods taking (row, column) do NOT bounds-check; indices must be in
/// [0, m) x [0, n). The Gauss-elimination helpers (subrow, eliminate, inverse) only maintain the
/// row lists, so the column lists of their result are not valid afterwards.
///
/// ## Numeric correctness guarantee
/// `solveIterative()` runs Gauss-Seidel with ε=1e-4 / maxiter=10000.
/// The formula is:
///   x_k ← b̃_k − Σ_{j≠k} ã_{kj} x_j
/// where ã = A with diagonal divided out and b̃ = b / a_{ii}.
class SMat {
public:
    int m = 0; ///< Number of rows.
    int n = 0; ///< Number of columns.

    /// Head of each row's linked list (indexed by row; entries are sorted by col).
    std::vector<Sparse *> rows;
    /// Head of each column's linked list (indexed by col; entries are sorted by row).
    std::vector<Sparse *> cols;

    /// @brief Construct an m × n zero matrix.
    SMat(int m, int n) {
        init(m, n);
    }

    /// @brief Deep-copy constructor.
    SMat(const SMat &mat) {
        init(mat.m, mat.n);
        copy(mat);
    }

    /// @brief Move constructor — takes ownership.
    SMat(SMat &&other) noexcept
        : m(other.m), n(other.n),
          rows(std::move(other.rows)),
          cols(std::move(other.cols)) {
        other.m = 0;
        other.n = 0;
    }

    /// @brief Copy-assignment.
    SMat &operator=(const SMat &mat) {
        if (this != &mat) {
            clear();
            rows.assign(static_cast<std::size_t>(mat.m), nullptr);
            cols.assign(static_cast<std::size_t>(mat.n), nullptr);
            m = mat.m;
            n = mat.n;
            copy(mat);
        }
        return *this;
    }

    /// @brief Move-assignment.
    SMat &operator=(SMat &&other) noexcept {
        if (this != &other) {
            clear();
            m = other.m;
            n = other.n;
            rows = std::move(other.rows);
            cols = std::move(other.cols);
            other.m = 0;
            other.n = 0;
        }
        return *this;
    }

    /// @brief Destructor — frees all Sparse nodes.
    ~SMat() {
        clear();
    }

    /// @brief (Re-)initialise as an m × n zero matrix.
    ///
    /// Existing nodes are freed first.
    void init(int newM, int newN) {
        clear();
        m = newM;
        n = newN;
        rows.assign(static_cast<std::size_t>(m), nullptr);
        cols.assign(static_cast<std::size_t>(n), nullptr);
    }

    /// @brief Deep-copy all entries from `mat` into `this`.
    ///
    /// @pre `this` and `mat` have the same dimensions.
    void copy(const SMat &mat) {
        for (int i = 0; i < m; ++i) {
            Sparse *element = mat.rows[static_cast<std::size_t>(i)];
            while (element) {
                insert(i, element->col, element->value);
                element = element->nextCol;
            }
        }
    }

    /// @brief Free all Sparse nodes; leave row/col head arrays intact (zeroed).
    void clear() {
        // Every node is on exactly one row list, so freeing via the rows frees all of them.
        for (int i = 0; i < m; ++i) {
            Sparse *element = rows[static_cast<std::size_t>(i)];
            while (element) {
                Sparse *garbage = element;
                element = element->nextCol;
                delete garbage;
            }
        }

        // The nodes are gone; drop the dangling heads of both index arrays.
        for (Sparse *&head : rows) {
            head = nullptr;
        }
        for (Sparse *&head : cols) {
            head = nullptr;
        }
    }

    /// @brief Set diagonal to 1, all other entries to 0.
    void identity() {
        clear();
        for (int i = 0; i < m; ++i) {
            insert(i, i, 1.0);
        }
    }

    /// @brief Get or create the node at (row, column).
    ///
    /// If the node does not exist it is created with value 0.  The caller
    /// can then set `node->value`.
    ///
    /// @pre 0 <= row < m and 0 <= column < n (unchecked).
    /// @return Pointer to the Sparse node (never null).
    Sparse *get(int row, int column) {
        Sparse **left = &rows[static_cast<std::size_t>(row)];
        Sparse **top = &cols[static_cast<std::size_t>(column)];

        // Walk the row list to the first node at or past the column.
        while ((*left) && ((*left)->col < column)) {
            left = &((*left)->nextCol);
        }

        // Already stored: hand it out instead of creating a duplicate.
        if ((*left) && ((*left)->col == column)) {
            return *left;
        }

        // New node: it must also be linked into the column list at its row position.
        while ((*top) && ((*top)->row < row)) {
            top = &((*top)->nextRow);
        }

        Sparse *element = new Sparse;
        element->row = row;
        element->col = column;
        element->value = 0.0;
        element->nextRow = *top;
        element->nextCol = *left;

        *left = element;
        *top = element;
        return element;
    }

    /// @brief Return the value at (row, column), or 0.0 if not stored.
    [[nodiscard]] real getvalue(int row, int column) const {
        const Sparse *left = rows[static_cast<std::size_t>(row)];
        while (left && left->col < column) {
            left = left->nextCol;
        }
        if (left && left->col == column) {
            return left->value;
        }
        return 0.0;
    }

    /// @brief Insert value at (row, column).  Out-of-range indices and zero values are ignored
    ///        (so an existing entry cannot be overwritten with 0).
    void insert(int row, int column, real value) {
        if (row < 0 || column < 0 || row >= m || column >= n || value == 0.0) {
            return;
        }
        get(row, column)->value = value;
    }

    /// @brief Row-times-vector: return row `rowid` · `vec`.
    [[nodiscard]] real mul(int rowid, const real *vec) const {
        real          result = 0.0;
        const Sparse *row = rows[static_cast<std::size_t>(rowid)];
        while (row) {
            result += row->value * vec[row->col];
            row = row->nextCol;
        }
        return result;
    }

    /// @brief Compute `dst = A * vec` (allocates / replaces `*dst`; `vec` must have n entries).
    /// @note LEGACY-SHAPED API: the caller owns `*dst` and must `delete[]` it.
    ///       The raw owning pointer is preserved deliberately to keep the
    ///       signature the LSCM port (Spec 08) consumes; full RAII conversion
    ///       (e.g. `std::vector<real>` return) is deferred to that spec.
    void mul(real **dst, const real *vec) const {
        real *tmp = new real[static_cast<std::size_t>(m)];
        for (int i = 0; i < m; ++i) {
            tmp[i] = mul(i, vec);
        }

        // Build into tmp first so *dst may alias vec.
        if (*dst) {
            delete[] *dst;
        }
        *dst = tmp;
    }

    /// @brief Dot product of sparse row and sparse column.
    [[nodiscard]] static real mulrowbycolumn(const Sparse *row, const Sparse *col) {
        real result = 0.0;
        while (row && col) {
            if (row->col == col->row) {
                result += row->value * col->value;
                row = row->nextCol;
                col = col->nextRow;
            } else if (row->col < col->row) {
                row = row->nextCol;
            } else {
                col = col->nextRow;
            }
        }
        return result;
    }

    /// @brief Dot product of sparse row and sparse column, stopping at `rowbreak`.
    [[nodiscard]] static real mulrowbycolumn(const Sparse *row, const Sparse *col,
                                             int rowbreak) {
        real result = 0.0;
        while (row && col && row->col < rowbreak && col->row < rowbreak) {
            if (row->col == col->row) {
                result += row->value * col->value;
                row = row->nextCol;
                col = col->nextRow;
            } else if (row->col < col->row) {
                row = row->nextCol;
            } else {
                col = col->nextRow;
            }
        }
        return result;
    }

    /// @brief Dot product of dense row array and sparse column.
    [[nodiscard]] static real mulrowbycolumn(const real *row, const Sparse *col) {
        real result = 0.0;
        while (col) {
            result += row[col->row] * col->value;
            col = col->nextRow;
        }
        return result;
    }

    /// @brief Sparse matrix multiplication: return `this * m2`.
    /// @pre n == m2.m (unchecked).
    [[nodiscard]] SMat operator*(const SMat &m2) const {
        SMat              result(m, m2.n);
        std::vector<real> frows(static_cast<std::size_t>(n), 0.0);

        for (int i = 0; i < result.m; ++i) {
            // Scatter row i into a dense buffer so each column product is a plain lookup.
            std::fill(frows.begin(), frows.end(), 0.0);
            const Sparse *loop = rows[static_cast<std::size_t>(i)];
            while (loop) {
                frows[static_cast<std::size_t>(loop->col)] = loop->value;
                loop = loop->nextCol;
            }

            // Results arrive in column order, so append at the tail of the row list.
            Sparse **element = &result.rows[static_cast<std::size_t>(i)];
            for (int j = 0; j < result.n; ++j) {
                real value = mulrowbycolumn(frows.data(), m2.cols[static_cast<std::size_t>(j)]);
                if (value != 0.0) {
                    Sparse *sp = new Sparse;
                    sp->col = j;
                    sp->nextCol = nullptr;
                    sp->nextRow = nullptr;
                    sp->row = i;
                    sp->value = value;
                    *element = sp;
                    element = &(sp->nextCol);
                }
            }
        }
        return result;
    }

    /// @brief Ported product `rows[j] . cols[i]` per entry (legacy2003/math/sparsematrix.h:358).
    ///
    /// NOTE: that is sum_k A(j,k) * A(k,i), i.e. a square matrix product, which equals A^T A
    /// only when A is symmetric. Kept as-is (behavior-preserving, antipattern 11).
    [[nodiscard]] SMat ata() const {
        SMat result(n, n);

        for (int i = 0; i < result.m; ++i) {
            Sparse **element = &result.rows[static_cast<std::size_t>(i)];
            for (int j = 0; j < result.n; ++j) {
                real value = mulrowbycolumn(rows[static_cast<std::size_t>(j)],
                                            cols[static_cast<std::size_t>(i)]);
                if (value != 0.0) {
                    Sparse *sp = new Sparse;
                    sp->col = j;
                    sp->nextCol = nullptr;
                    sp->nextRow = nullptr;
                    sp->row = i;
                    sp->value = value;
                    *element = sp;
                    element = &(sp->nextCol);
                }
            }
        }
        return result;
    }

    /// @brief Transpose `this` in-place and return a reference.
    SMat &transpose() {
        // Swap the roles of the two links on every node; the loop then follows the old nextRow
        // link, which after the swap is the old nextCol, i.e. it still walks the row.
        for (int i = 0; i < m; ++i) {
            Sparse *element = rows[static_cast<std::size_t>(i)];
            while (element) {
                std::swap(element->col, element->row);
                std::swap(element->nextCol, element->nextRow);
                element = element->nextRow;
            }
        }

        std::swap(n, m);
        std::swap(rows, cols);
        return *this;
    }

    /// @brief Negate all entries in-place.
    void neg() {
        for (int i = 0; i < m; ++i) {
            Sparse *element = rows[static_cast<std::size_t>(i)];
            while (element) {
                element->value = -element->value;
                element = element->nextCol;
            }
        }
    }

    /// @brief Divide all entries in a sparse row by `divider`.
    static void divide(Sparse *row, real divider) {
        while (row) {
            row->value /= divider;
            row = row->nextCol;
        }
    }

    /// @brief Divide row `row` by `divider`.
    void divide(int row, real divider) {
        divide(rows[static_cast<std::size_t>(row)], divider);
    }

    /// @brief Subtract `modifier * src` from the row pointed to by `element`.
    ///
    /// Entries that become exactly 0 are unlinked and freed. Only the row links are maintained;
    /// new nodes get `row = 0` and no column links.
    static void subrow(Sparse **element, const Sparse *src, real modifier) {
        while (src) {
            // Advance in the destination row to src's column (both rows are col-ascending).
            while (*element && (*element)->col < src->col) {
                element = &((*element)->nextCol);
            }

            if (!(*element) || (*element)->col > src->col) {
                Sparse *newelem = new Sparse;
                newelem->col = src->col;
                newelem->nextCol = *element;
                newelem->row = 0; // row field unused here
                newelem->nextRow = nullptr;
                newelem->value = -modifier * src->value;
                *element = newelem;
            } else if (src->col == (*element)->col) {
                (*element)->value -= modifier * src->value;
            }

            // Cancelled to exactly 0: keep the sparse structure sparse.
            if (std::fabs((*element)->value) == 0.0) {
                Sparse *garbage = *element;
                *element = (*element)->nextCol;
                delete garbage;
            }
            src = src->nextCol;
        }
    }

    /// @brief Gauss elimination: remove the leading coefficient of `*element`
    ///        using `src` as the pivot row.
    void eliminate(Sparse **element, const Sparse *src) {
        if ((*element)->col != src->col) {
            return;
        }
        // Remove the leading entry (it becomes 0 by construction), then subtract the rest of the pivot row.
        real    modifier = (*element)->value / src->value;
        Sparse *garbage = *element;
        *element = (*element)->nextCol;
        delete garbage;
        src = src->nextCol;
        subrow(element, src, modifier);
    }

    /// @brief Find the pivot row in `column` between rows `start` and `end`.
    ///
    /// Despite the max-search bookkeeping this returns the FIRST row whose head is in `column`
    /// (early return, legacy2003/math/sparsematrix.h:512); `start` if none.
    static int pivotSearch(Sparse **rows, int column, int start, int end) {
        real maxVal = 0.0;
        int  pivot = start;
        for (int i = start; i < end; ++i) {
            if (rows[i] && rows[i]->col == column &&
                std::fabs(rows[i]->value) > maxVal) {
                maxVal = std::fabs(rows[i]->value);
                pivot = i;
                return pivot;
            }
        }
        return pivot;
    }

    /// @brief Gauss-Jordan back-substitution: clear above-diagonal entries and
    ///        normalise each row so the left block becomes the identity.
    /// @pre every localRows[i] is non-empty with a non-zero head (not guaranteed after forwardEliminate).
    void backSubstitute(Sparse **localRows) {
        for (int i = m - 1; i >= 0; --i) {
            const Sparse *src = localRows[i];

            for (int j = i - 1; j >= 0; --j) {
                Sparse **elem = &localRows[j];
                while (*elem && (*elem)->col < src->col) {
                    elem = &((*elem)->nextCol);
                }
                if (*elem && (*elem)->col == src->col) {
                    subrow(elem, src, (*elem)->value / src->value);
                }
            }
            divide(localRows[i], localRows[i]->value);
        }
    }

    /// @brief Compute A^{-1} via Gauss-Jordan elimination on [A | I].
    ///
    /// @param result  Optional output (may be null); existing contents are replaced.
    /// @return        The reduced matrix. For a singular A it is the PARTIALLY reduced [A | I],
    ///                not an error value; the caller cannot tell the two cases apart.
    SMat inverse(SMat *result) {
        SMat inv(m, n * 2);
        inv.copy(*this);

        // Augment with the identity in columns n..2n-1.
        for (int i = 0; i < m; ++i) {
            inv.insert(i, i + n, 1.0);
        }

        Sparse **localRows = inv.rows.data();
        if (forwardEliminate(inv, localRows)) {
            backSubstitute(localRows);
            dropIdentityOffset(inv);
        }

        // Singular case falls through with a partially-reduced matrix.
        if (result) {
            *result = inv;
        }
        return inv;
    }

    /// @brief True when `row` has a usable pivot in `column`: a stored, non-zero
    ///        entry whose column is exactly `column`.
    static bool hasPivot(const Sparse *row, int column) {
        return row && row->col == column && std::fabs(row->value) != 0.0;
    }

    /// @brief Gauss-Jordan forward elimination with partial pivoting.
    /// @return false if any pivot column is empty (matrix singular).  Every pivot is
    ///         checked, the last included: the elimination loop stops one row early,
    ///         yet backSubstitute() dereferences the last row.
    bool forwardEliminate(SMat &inv, Sparse **localRows) {
        for (int i = 0; i < m - 1; ++i) {
            const int pivot = pivotSearch(localRows, i, i, m);
            std::swap(localRows[i], localRows[pivot]);

            // No usable pivot in this column: the matrix is singular.
            if (!hasPivot(localRows[i], i)) {
                return false;
            }

            for (int j = i + 1; j < m; ++j) {
                if (localRows[j] && localRows[j]->col == i &&
                    localRows[j]->value != 0.0) {
                    inv.eliminate(&localRows[j], localRows[i]);
                }
            }
        }

        // The last row has no row below it to eliminate, but it is a pivot too.
        return m == 0 || hasPivot(localRows[m - 1], m - 1);
    }

    /// @brief Shift every entry's column left by n (undo the [A|I] augmentation).
    ///
    /// Only the nodes' `col` is rewritten; `inv.n`, `inv.cols` and the `row` fields are not.
    void dropIdentityOffset(SMat &inv) const {
        for (int i = 0; i < inv.m; ++i) {
            for (Sparse *loop = inv.rows[static_cast<std::size_t>(i)]; loop;
                 loop = loop->nextCol) {
                loop->col -= n;
            }
        }
    }

    /// @brief Solve A·x = b using Gauss-Seidel iteration.
    ///
    /// ## Algorithm (behavior-preserving from legacy)
    ///
    /// 1. Normalise b: b̃_i = b_i / a_{ii}.
    /// 2. Normalise A: ã_{ij} = a_{ij} / a_{ii}  (for j ≠ i).
    /// 3. Iterate until convergence or maxiter:
    ///      x_k ← b̃_k − Σ_{j≠k} ã_{kj} x_j
    ///
    /// Convergence criterion: ‖Δx‖₂ < 1e-4.
    /// Maximum iterations: 10000.
    ///
    /// @note All `printf` debug output has been removed (Spec 05 core-purity
    ///       constraint).  The numeric result is unchanged.
    ///
    /// @param srcb  Right-hand-side column vector (m × 1 SMat); must not be null.
    /// @return      Solution vector x as an (m × 1) SMat. Not-converged after the iteration cap is
    ///              returned silently; a zero diagonal entry produces inf/NaN.
    [[nodiscard]] SMat solveIterative(SMat *srcb) {
        SMat a(m, n);
        SMat b(*srcb);
        SMat x(b.m, 1);

        // Divide b by the diagonal and seed x with explicit zeros, so b and x have a node
        // in every row and the loop below can update values in place.
        for (int i = 0; i < b.m; ++i) {
            real aii = getvalue(i, i);
            if (b.rows[static_cast<std::size_t>(i)]) {
                b.rows[static_cast<std::size_t>(i)]->value /= aii;
            } else {
                b.insert(i, 0, 1.0);
                b.rows[static_cast<std::size_t>(i)]->value = 0.0;
            }
            x.insert(i, 0, 1.0);
            x.rows[static_cast<std::size_t>(i)]->value = 0.0;
        }


        // Off-diagonal part of A scaled by its diagonal; the diagonal itself is folded into b.
        for (int i = 0; i < m; ++i) {
            real aii = getvalue(i, i);
            for (int k = 0; k < n; ++k) {
                if (k != i) {
                    a.insert(i, k, getvalue(i, k) / aii);
                }
            }
        }

        constexpr int  kMaxIter = 10000;
        constexpr real kConv = 1.0e-4;
        int            iter = 0;
        real           diff = 0.0;

        // Gauss-Seidel sweeps update x in place, so later rows already see the new values.
        do {
            diff = 0.0;
            for (int k = 0; k < x.m; ++k) {
                real newvalue = b.getvalue(k, 0) -
                                mulrowbycolumn(a.rows[static_cast<std::size_t>(k)],
                                               x.cols[0]);
                real diff2 = x.rows[static_cast<std::size_t>(k)]->value - newvalue;
                diff += diff2 * diff2;
                x.rows[static_cast<std::size_t>(k)]->value = newvalue;
            }
            diff = std::sqrt(diff);
        } while (iter++ < kMaxIter && diff > kConv);

        return x;
    }

    /// @brief LR (LU) factorisation into lower and upper triangular matrices.
    void lrTriangulate(SMat &l, SMat &r) {
        l.identity();
        r.identity();

        // R's first row is A's first row.
        const Sparse *loop = rows[0];
        while (loop) {
            r.insert(0, loop->col, loop->value);
            loop = loop->nextCol;
        }

        for (int i = 0; i < m; ++i) {
            // L(i,k) for k < i: eliminate against the already known R rows.
            for (int k = 0; k < i; ++k) {
                real value = getvalue(i, k);
                value -= mulrowbycolumn(l.rows[static_cast<std::size_t>(i)],
                                        r.cols[static_cast<std::size_t>(k)], k);
                l.insert(i, k, value / r.getvalue(k, k));
            }

            // R(i,k) for k >= i.
            for (int k = i; k < m; ++k) {
                real value = getvalue(i, k);
                value -= mulrowbycolumn(l.rows[static_cast<std::size_t>(i)],
                                        r.cols[static_cast<std::size_t>(k)], i);
                r.insert(i, k, value);
            }
        }
    }
};

} // namespace decimate::math


namespace decimate { // transitional re-export (canonical: math::)
using math::SMat;
using math::Sparse;
} // namespace decimate