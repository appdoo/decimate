// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file sparse_bit_set.h
/// @brief Sparse bit set — behavior-preserving port of
///        legacy2003/math/csparsebitset.{h,cpp} (Spec 05, Beat 2a).
///
/// ## Design
///
/// `CSparseBitSet` stores bit-plane words sparsely: only the 32-bit words
/// that actually contain set bits are kept, identified by their word offset
/// (bit_index >> 5).  Each word is stored in a dynamic array that grows
/// by 5 entries on demand.
///
/// ## Port notes
///
/// Storage: `int *bits`, `int *offs` → `std::vector<int>` for both. The legacy
/// Init-5 / grow-by-5 raw-array steps were a representation detail (no observable
/// effect) and are dropped — the vector grows itself.
///
/// ## Known bug (characterized per antipattern 11 — do NOT fix)
///
/// `erase()` contains two defects from 2009:
///
/// 1. **`!`-vs-`~` bug (legacy2003/math/csparsebitset.cpp:103):**
///    `bits[offset] &= !(1 << bit & 0x1F)` uses logical-NOT instead of
///    bitwise-NOT (`~`).  For any nonzero mask `!(mask) == 0`, so the entire
///    word is cleared instead of just the target bit.
///
/// 2. **Off-by-one guard (`if (offset > 0)`):**
///    `Find()` returns the *array index* (not the word offset).  When the
///    word containing the target bit was the first word inserted (array
///    index 0), the condition `offset > 0` is false and `erase()` does
///    nothing at all — the bit is NOT cleared.
///
/// These bugs are preserved exactly.  `erase()` for a bit in array slot 0
/// is a no-op; `erase()` for a bit in array slot > 0 clears the whole word.
/// Tests T12 and T13 characterize both behaviours.
///
/// A correct erase would be:
///   `if (offset >= 0) bits[offset] &= ~(1 << (bit & 0x1F));`
/// That fix is deferred to Beat 2b / Spec 12.

#include <vector>

namespace decimate::math {

/// @brief Sparse bit set — stores only the populated 32-bit words.
///
/// ## Word layout
/// Bit `b` lives in word `b >> 5` at position `b & 0x1F`.
class CSparseBitSet {
public:
    /// @brief Construct an empty bitset. (Rule of Zero — vector owns storage.)
    CSparseBitSet() = default;

    /// @brief Set bit `bit`.
    ///
    /// If the word containing `bit` is not yet allocated, a new word entry
    /// is added (capacity grows by 5 if needed).
    void insert(int bit);

    /// @brief Clear bit `bit`.
    ///
    /// @warning This method contains the characterized 2009 bugs:
    ///   - If the word's array index is 0 → this is a no-op (bit not cleared).
    ///   - If the word's array index is > 0 → the *entire word* is cleared
    ///     (`!`-vs-`~` wipes the word to 0).
    ///
    /// Do NOT fix these — they are the source-of-truth behaviour (T13).
    void erase(int bit);

    /// @brief Test whether bit `bit` is set.
    ///
    /// @return `true` if the bit is set, `false` if the word is not stored
    ///         or the specific bit within the word is 0.
    [[nodiscard]] bool test(int bit) const;

private:
    std::vector<int> bits_; ///< 32-bit word values, indexed by array position.
    std::vector<int> offs_; ///< Word offsets (bit >> 5), one per entry.

    /// @brief Find the array index for word `wordOffset`, or -1 if absent.
    [[nodiscard]] int find(int wordOffset) const;

    /// @brief Append a new empty word entry for `wordOffset`; returns its index.
    int insertWord(int wordOffset);
};

} // namespace decimate::math


namespace decimate {
using math::CSparseBitSet;
} // namespace decimate