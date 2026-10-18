// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file quadric_error_system.h
/// @brief Per-vertex quadric accumulator — port of legacy `CQuadricErrorSystem`.
///
/// `QuadricErrorSystem` owns a flat array of `QError` values (one per vertex),
/// accumulates face-plane quadrics during initialisation, merges on collapse,
/// and computes the optimal collapse position and error.
///
/// Modernisation from legacy:
///   - `new QError[]` → `std::vector<QError>` (no naked new/delete).
///   - No `Destroy()` method — destructor is the destructor.
///   - `Initialize` takes `const` pointers (callee does not mutate input data).

#include <decimate/reduce/qerror.h>

#include <decimate/math/vec3.h>
#include <vector>

#include <decimate/types/mesh_types.h>

namespace decimate {

/// @brief Quadric-error metric system for the QEM edge-collapse reducer.
///
/// Lifecycle: `initialize` → repeated `optimize` / `update` calls → destroyed.
///
/// The system holds no back-pointer to the mesh; vertex positions are passed
/// by reference at `initialize` time and stored as a non-owning span so that
/// `optimize` (which needs the midpoint fallback) can read current positions.
class QuadricErrorSystem {
public:
    QuadricErrorSystem() = default;

    // Rule of Zero: vector handles copy/move/destroy.

    /// @brief Build per-vertex quadrics from all face planes.
    ///
    /// For each face f, computes the face-plane quadric K_p and adds it
    /// to the three corner vertices.  O(F) time, O(V) space.
    ///
    /// @param numPoints  Number of vertices.
    /// @param points     Vertex positions array (non-owning pointer, lifetime
    ///                   must span the reducer's use of this system).
    /// @param numFaces   Number of faces.
    /// @param faces      Face index array.
    void initialize(int numPoints, const decimate::math::Vec3f* points,
                    int numFaces, const Face* faces);

    /// @brief Notify of an edge collapse (p1 absorbs p2).
    ///
    /// Merges p2's quadric into p1.  O(1) — 10 double additions.
    ///
    /// @param p1  Surviving vertex index.
    /// @param p2  Removed vertex index.
    void update(int p1, int p2);

    /// @brief Compute the optimal position and error for collapsing edge (p1, p2).
    ///
    /// Forms the combined quadric Q = Q[p1] + Q[p2], then:
    ///   - If `boundary` or the 4×4 inversion fails: position = midpoint of
    ///     the current vertex positions.  (Beat-2a behaviour — F5 deferred.)
    ///   - Otherwise: position = QEM-optimal via matrix inversion.
    ///
    /// @param p1        First vertex index.
    /// @param p2        Second vertex index.
    /// @param boundary  True if the edge is on the mesh boundary.
    /// @param[out] error  Quadric error at the chosen position.
    /// @param[out] p      Chosen position.
    void optimize(int p1, int p2, bool boundary, double& error, decimate::math::Vec3f& p) const;

private:
    std::vector<QError> errorMetrics_;
    const decimate::math::Vec3f*  pointRef_ = nullptr;  ///< Non-owning; set by initialize.
};

}  // namespace decimate
