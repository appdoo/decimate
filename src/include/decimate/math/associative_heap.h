// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file associative_heap.h
/// @brief Intrusive associative min-heap for the QEM priority queue.
///
/// Elements store their own heap position (`heapPos`) and sort key
/// (`heapValue`) by inheriting from `Heapable`.  This intrusive design lets
/// the heap re-key an element in O(log n) without a separate handle map —
/// the decrease-key path the QEM edge-collapse loop depends on.
///
/// `AssociativeHeap<T>` is templated on the element type (which must derive
/// from `Heapable`) so `pop()` / `remove()` return `T*` directly — the QEM
/// loop gets a `VertexPair*` with no downcast.  Header-only + inline so the
/// collapse hot loop can inline the sift operations.
///
/// Storage is a `std::vector<Heapable*>` of NON-owning pointers: the heap
/// orders elements it does not own (callers manage element lifetime, which
/// antipattern 9 permits for observers).  `heap_.size()` is the sole source
/// of truth for the element count.

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace decimate::math {

/// @brief Intrusive base for heap elements.
///
/// A type stored in `AssociativeHeap` must publicly inherit `Heapable` and
/// set a valid `heapValue` before `push()` / `repair()`.  `heapPos` is owned
/// by the heap — callers must not write it.
class Heapable {
public:
    double heapValue = 0.0;  ///< Sort key — lower values bubble to the top.
    int    heapPos   = -1;   ///< Slot in the heap array (heap-managed).
};

/// @brief Intrusive binary min-heap with O(log n) decrease-key.
///
/// Push candidates, pop the minimum, and re-key (`repair`) an element whose
/// `heapValue` changed — the QEM collapse loop.  Holds non-owning pointers.
///
/// @tparam T  Element type; must publicly derive from `Heapable`.
///            Defaults to `Heapable` for backwards-compatible use where the
///            caller accesses elements via the `Heapable` base.
template <typename T = Heapable>
class AssociativeHeap {
    static_assert(std::is_base_of_v<Heapable, T>,
                  "AssociativeHeap<T>: T must derive from Heapable");

public:
    /// @param initialSize  Capacity to reserve up front (grows as needed).
    explicit AssociativeHeap(int initialSize = 100) {
        init(initialSize);
    }

    // Rule of Zero — the vector handles copy/move/destroy.

    /// @brief Discard all elements and reserve a fresh capacity.
    void init(int initialSize) {
        heap_.clear();
        heap_.reserve(static_cast<std::size_t>(initialSize));
    }

    /// @brief Insert @p element, set its `heapPos`, and sift it into place.
    void push(T* element) {
        element->heapPos = size();
        heap_.push_back(element);
        upHeap(element->heapPos);  // upHeap(0) is a no-op — no empty special case
    }

    /// @brief Remove and return the minimum element, or nullptr if empty.
    T* pop() { return remove(0); }

    /// @brief Remove and return the element at @p position, or nullptr if empty.
    T* remove(int position) {
        if (heap_.empty()) {
            return nullptr;
        }
        Heapable* removed = heap_[static_cast<std::size_t>(position)];
        Heapable* last    = heap_.back();
        heap_.pop_back();
        // If the removed slot still exists, move the former last element into
        // it and let it settle.  (When `position` was the last slot, it is gone.)
        if (position < size()) {
            at(position)  = last;
            last->heapPos = position;
            downHeap(position);
        }
        return static_cast<T*>(removed);
    }

    /// @brief Remove and return @p element (located via its `heapPos`).
    T* remove(T* element) { return remove(element->heapPos); }

    /// @brief Re-key the element at @p position after its `heapValue` changed.
    void repair(int pos) {
        // Smaller than parent → sift up; otherwise sift down.
        if (pos > 0 && at(pos)->heapValue < at((pos - 1) / 2)->heapValue) {
            upHeap(pos);
        } else {
            downHeap(pos);
        }
    }

    /// @brief Re-key @p element after its `heapValue` changed.
    void repair(T* element) { repair(element->heapPos); }

    /// @brief Number of elements currently in the heap.
    [[nodiscard]] int size() const noexcept {
        return static_cast<int>(heap_.size());
    }

    /// @brief True if the heap holds no elements.
    [[nodiscard]] bool empty() const noexcept { return heap_.empty(); }

private:
    /// Int-indexed accessor — heap arithmetic is in `int`; cast in one place.
    Heapable*& at(int i) { return heap_[static_cast<std::size_t>(i)]; }

    void upHeap(int k) {
        while (k > 0) {
            const int parent = (k - 1) / 2;
            if (at(k)->heapValue >= at(parent)->heapValue) {
                break;
            }
            std::swap(at(k), at(parent));
            at(k)->heapPos      = k;
            at(parent)->heapPos = parent;
            k = parent;
        }
    }

    void downHeap(int pos) {
        const int n = size();
        for (;;) {
            int child = (2 * pos) + 1;
            if (child >= n) {
                break;
            }
            if (child + 1 < n && at(child + 1)->heapValue < at(child)->heapValue) {
                ++child;  // pick the smaller child
            }
            if (at(pos)->heapValue < at(child)->heapValue) {
                break;  // already in order
            }
            std::swap(at(pos), at(child));
            at(pos)->heapPos   = pos;
            at(child)->heapPos = child;
            pos = child;
        }
    }

    std::vector<Heapable*> heap_;  ///< Non-owning storage; size() is the count.
};

}  // namespace decimate::math

namespace decimate {  // transitional re-export
using math::Heapable;
using math::AssociativeHeap;
}