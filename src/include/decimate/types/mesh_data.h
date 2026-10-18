// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file mesh_data.h
/// @brief Public API for the mesh-core I/O module (Spec 06).
///
/// MeshData is the value-type output of every loader in this module.
/// It owns its storage (std::vector) and follows normal value semantics —
/// copy, move, assign all work correctly without any manual resource management.
///
/// Representation contract (behavior-preserving from legacy):
///   positions  == legacy coord3, one decimate::math::Vec3f per vertex (the flat legacy
///                 coord3[0 .. 3*numCoord) packed XYZ, now typed: positions[i]
///                 holds vertex i. Vec3f is standard-layout and tightly packed,
///                 so positions.data() is byte-identical to the legacy float[].
///   indices    == legacy faces[0 .. 3*numFaces)    — 0-based triangle indices
///                 (legacy stored 1-based and decremented by -1 in read_face3)
///
/// FROZEN BUGS (Beat-2a — characterize, do NOT fix):
///   Bug 1 (ObjLdr.cpp:83-84): row[0]=='v' dispatch routes 'vn' lines to
///         read_vertex3, appending normal XYZ as a phantom position vertex.
///   Bug 2 (read_face3:24): atoi("N//M") == N — slash-separated face refs
///         silently discard vt/vn attributes; only the vertex index survives.
///
/// Locale-dependent parsing (atof/atoi) is also preserved from legacy —
/// TODO (Beat 2b): replace with locale-independent from_chars / fast_float.
///
/// Boundary contracts: no GL/GLUT/windows.h; no stdout/stderr I/O;
/// no global mutable state.

#include <string>
#include <string_view>
#include <vector>

#include <decimate/math/vec3.h>

namespace decimate {

/// @brief Mesh representation: vertex positions + triangle indices.
///
/// positions: one decimate::math::Vec3f per vertex.
///   positions[i] is vertex i's XYZ. positions.size() is the vertex count.
///   Vec3f is standard-layout and tightly packed (sizeof == 3*float), so
///   reinterpret_cast<const float*>(positions.data()) yields the legacy flat
///   coord3 buffer unchanged — safe for GL upload / raw-buffer interop.
///
/// indices: packed triangle indices (0-based), 3 per triangle.
///   indices[3*t], indices[3*t+1], indices[3*t+2] are the three vertex indices
///   of triangle t.
///
/// Both vectors are empty for an empty or failed load. The caller checks the
/// size (or the thrown exception) rather than a bool return.
struct MeshData {
    std::vector<decimate::math::Vec3f> positions;
    std::vector<int>         indices;
};

/// @brief Load an OBJ mesh from @p path.
///
/// Behavior-preserves ObjLdr::load exactly (Spec 06, Beat 2a):
///   - Dispatch on row[0] only → 'vn' lines become phantom positions (Bug 1).
///   - atoi on face tokens → slash refs truncated to vertex index (Bug 2).
///   - atof/atoi for numeric parsing (locale-dependent — deferred to Beat 2b).
///
/// @param path  File-system path to the .obj file.
/// @returns MeshData owning the vertex and index vectors.
/// @throws std::runtime_error if the file cannot be opened.
MeshData loadObj(std::string_view path);

/// @brief Load an Open Inventor ASCII (.iv) mesh from @p path.
///
/// Behavior-preserves InventorLdr::load exactly (Spec 06, Beat 2a):
///   - Recursive descent: Separator → Coordinate3 / IndexedFaceSet.
///   - Triangle-only: reads exactly three indices + the -1 terminator per face.
///   - Material and MaterialBinding blocks are silently skipped.
///
/// @param path  File-system path to the .iv file.
/// @returns MeshData owning the vertex and index vectors.
/// @throws std::runtime_error if the file cannot be opened or is malformed.
MeshData loadInventor(std::string_view path);

/// @brief Extension dispatcher: routes to loadObj or loadInventor by file suffix.
///
/// ".obj" (case-insensitive suffix) → loadObj
/// ".iv"  (case-insensitive suffix) → loadInventor
///
/// Mirrors the dispatch in the 2009 main.cpp:134-137.
///
/// @param path  File-system path; extension determines the parser.
/// @returns MeshData from the appropriate loader.
/// @throws std::runtime_error for unrecognised extensions or load failures.
MeshData load(std::string_view path);

}  // namespace decimate
