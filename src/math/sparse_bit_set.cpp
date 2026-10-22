// This file is part of the decimate library; see decimate.h for version/license details.
/// @file sparse_bit_set.cpp
/// @brief CSparseBitSet implementation — behavior-preserving port of
///        legacy2003/math/csparsebitset.cpp (Spec 05, Beat 2a).
///
/// All logic including the characterized bugs is preserved exactly.
/// See sparse_bit_set.h for the bug description.

#include <decimate/math/sparse_bit_set.h>

namespace decimate::math {

/// Linear scan for the array slot holding word @p wordOffset, or -1.
/// The sparse word count is small, so linear is fine.
int CSparseBitSet::find(int wordOffset) const {
    for (std::size_t i = 0; i < offs_.size(); ++i) {
        if (offs_[i] == wordOffset) {
            return static_cast<int>(i);
        }
    }

    return -1;
}


/// Appends a new empty word slot for @p wordOffset and returns its index.
/// The vector grows itself; the legacy grow-by-5 raw array had no observable effect.
int CSparseBitSet::insertWord(int wordOffset) {
    offs_.push_back(wordOffset);
    bits_.push_back(0);

    return static_cast<int>(offs_.size()) - 1;
}


void CSparseBitSet::insert(int bit) {
    int wordOffset = bit >> 5;
    int arrayIdx = find(wordOffset);

    // First bit in this word: allocate its slot.
    if (arrayIdx < 0) {
        arrayIdx = insertWord(wordOffset);
    }

    bits_[static_cast<std::size_t>(arrayIdx)] |= (1 << (bit & 0x1F));
}


/// Reproduces the 2009 erase() including both bugs (DO NOT FIX, characterized by tests):
///   Bug 1 (off-by-one guard): `if (offset > 0)` is strict, so array slot 0 is never erased.
///   Bug 2 (!-vs-~): `bits[offset] &= !(...)` — logical NOT of a nonzero mask is 0, so the whole word is zeroed.
void CSparseBitSet::erase(int bit) {
    int wordOffset = bit >> 5;
    int offset = find(wordOffset); // array index, not word offset

    // Bug 1: strict greater-than means array index 0 is skipped (no-op).
    if (offset > 0) {
        // Bug 2: logical-NOT wipes the entire word to 0.
        bits_[static_cast<std::size_t>(offset)] &=
            !(1 << (bit & 0x1F)); // NOLINT — intentional bug reproduction
    }
}


bool CSparseBitSet::test(int bit) const {
    int wordOffset = bit >> 5;
    int arrayIdx = find(wordOffset);

    // No word stored for this range: every bit in it is clear.
    if (arrayIdx < 0) {
        return false;
    }

    return (bits_[static_cast<std::size_t>(arrayIdx)] & (1 << (bit & 0x1F))) != 0;
}

} // namespace decimate::math
