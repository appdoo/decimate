// This file is part of the decimate library; see decimate.h for version/license details.
/// @file sparse_bit_set.cpp
/// @brief CSparseBitSet implementation — behavior-preserving port of
///        legacy/math/csparsebitset.cpp (Spec 05, Beat 2a).
///
/// All logic including the characterized bugs is preserved exactly.
/// See sparse_bit_set.h for the bug description.

#include <decimate/math/sparse_bit_set.h>

namespace decimate::math {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

// find — linear scan for the array slot holding word `wordOffset`, or -1.
// (Legacy Find(); the sparse word count is small, so linear is fine.)
int CSparseBitSet::find(int wordOffset) const {
    for (std::size_t i = 0; i < offs_.size(); ++i) {
        if (offs_[i] == wordOffset) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

// insertWord — append a new empty word slot for `wordOffset`, return its index.
// (Legacy grew a raw array by 5; std::vector handles growth — the legacy
// capacity steps were a representation detail with no behavioral effect.)
int CSparseBitSet::insertWord(int wordOffset) {
    offs_.push_back(wordOffset);
    bits_.push_back(0);
    return static_cast<int>(offs_.size()) - 1;
}

// ---------------------------------------------------------------------------
// insert — mirrors legacy insert()
// ---------------------------------------------------------------------------

void CSparseBitSet::insert(int bit) {
    int wordOffset = bit >> 5;
    int arrayIdx   = find(wordOffset);
    if (arrayIdx < 0) {
        arrayIdx = insertWord(wordOffset);
    }
    bits_[static_cast<std::size_t>(arrayIdx)] |= (1 << (bit & 0x1F));
}

// ---------------------------------------------------------------------------
// erase — mirrors legacy erase() including the bugs (DO NOT FIX)
// ---------------------------------------------------------------------------
//
// Bug 1 (off-by-one guard): `if (offset > 0)` uses strict-greater-than so
//   the entry at array index 0 is never erased.
// Bug 2 (!-vs-~): `bits[offset] &= !(1 << bit & 0x1F)` —
//   `!(nonzero) == 0` so the whole word is zeroed.
//
// Both bugs are reproduced faithfully below.

void CSparseBitSet::erase(int bit) {
    int wordOffset = bit >> 5;
    int offset     = find(wordOffset);  // array index, not word offset

    // Bug 1: strict greater-than means array index 0 is skipped (no-op).
    if (offset > 0) {
        // Bug 2: logical-NOT wipes the entire word to 0.
        bits_[static_cast<std::size_t>(offset)] &=
            !(1 << (bit & 0x1F));  // NOLINT — intentional bug reproduction
    }
}

// ---------------------------------------------------------------------------
// test
// ---------------------------------------------------------------------------

bool CSparseBitSet::test(int bit) const {
    int wordOffset = bit >> 5;
    int arrayIdx   = find(wordOffset);
    if (arrayIdx < 0) {
        return false;
    }
    return (bits_[static_cast<std::size_t>(arrayIdx)] & (1 << (bit & 0x1F))) != 0;
}

}  // namespace decimate::math