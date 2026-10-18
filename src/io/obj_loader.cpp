// This file is part of the decimate library; see decimate.h for version/license details.
/// @file obj_loader.cpp
/// @brief OBJ mesh loader — behavior-preserving port of legacy ObjLdr::load.
///
/// Ported from legacy/loader/ObjLdr.cpp (Spec 06, Beat 2a).
///
/// Representation changes (no numeric effect):
///   - Output is MeshData (owned std::vector) instead of raw new[]/out-params.
///   - The legacy vector → memcpy → new[] double-allocation is dropped: we
///     accumulate directly into the result vectors.
///   - std::string line buffer instead of fixed char[2048] — truncation at
///     2047 chars is not behavior worth preserving, and no test fixture or
///     model asset approaches that length.
///
/// Algorithm lines are UNCHANGED.  Both bugs are FROZEN (Beat 2a):
///
///   FROZEN BUG 1 (ObjLdr.cpp:83-84, row[0]=='v' only dispatch):
///     'vn' lines begin with 'v', so they route to read_vertex3.
///     Their XYZ normal components are appended as phantom vertex positions.
///     DO NOT add full-keyword matching here — that is the Beat-2b fix.
///
///   FROZEN BUG 2 (read_face3:24, atoi stops at '/'):
///     Face tokens like "251//255" → atoi("251//255") == 251.
///     Only the vertex index survives; vt/vn attributes are silently discarded.
///     DO NOT split on '/' to extract attribute indices — that is the Beat-2b fix.
///
/// NOTE (Beat 2b): atof/atoi are locale-dependent (LC_NUMERIC decimal separator).
/// This is the locale hazard from Finding B.  Replace with from_chars / fast_float.

#include <decimate/types/mesh_data.h>

#include <array>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace decimate {

// Legacy face-vertex cap: read_face3 stored at most 10 indices per face.
// Preserved (no model asset has a >10-gon; raising it is a Beat-2b decision).
constexpr int kMaxFaceVerts = 10;

// ---------------------------------------------------------------------------
// readVertex3 — port of ObjLdr::read_vertex3 (ObjLdr.cpp:50)
// ---------------------------------------------------------------------------
// Skips the first whitespace-separated token (the keyword, e.g. "v" or "vn"),
// then reads three whitespace-separated floats via atof.
//
// Called for every line whose first character is 'v' — including 'vn' and 'vt'
// lines (FROZEN BUG 1: row[0]-only dispatch in load()).
//
// NOTE (Beat 2b): atof is locale-dependent.

static void readVertex3(char* line, std::vector<decimate::math::Vec3f>& coords) {
    std::strtok(line, " \t");  // skip the keyword token
    decimate::math::Vec3f v;
    for (int i = 0; i < 3; ++i) {
        const char* token = std::strtok(nullptr, " \t");
        v[i] = token != nullptr ? static_cast<float>(std::atof(token)) : 0.0F;
    }
    coords.push_back(v);
}

// ---------------------------------------------------------------------------
// readFace3 — port of ObjLdr::read_face3 (ObjLdr.cpp:14)
// ---------------------------------------------------------------------------
// Reads up to kMaxFaceVerts whitespace-separated index tokens, converts each
// via atoi, decrements by 1 (1-based → 0-based), then fan-triangulates n-gons.
//
// FROZEN BUG 2: atoi("251//255") == 251 — slash-separated refs are silently
// truncated; only the vertex index survives.  DO NOT split on '/'.

static void readFace3(char* line, std::vector<int>& faces) {
    std::array<int, kMaxFaceVerts> indices{};
    std::strtok(line, " \t");  // skip the 'f' keyword
    int count = 0;
    for (const char* token = std::strtok(nullptr, " \t");
         token != nullptr && count < kMaxFaceVerts;
         token = std::strtok(nullptr, " \t")) {
        // FROZEN BUG 2: atoi stops at '/', discarding vt/vn attribute indices.
        indices[static_cast<std::size_t>(count++)] = std::atoi(token) - 1;
    }
    if (count < 3) {
        return;  // not a triangle
    }
    // First triangle, then fan-triangulate the rest (n-gon → n-2 triangles).
    const int anchor = indices[0];
    int       prev   = indices[1];
    for (int i = 2; i < count; ++i) {
        faces.push_back(anchor);
        faces.push_back(prev);
        faces.push_back(indices[static_cast<std::size_t>(i)]);
        prev = indices[static_cast<std::size_t>(i)];
    }
}

// ---------------------------------------------------------------------------
// loadObj — port of ObjLdr::load (ObjLdr.cpp:66)
// ---------------------------------------------------------------------------

MeshData loadObj(std::string_view path) {
    MeshData result;
    std::vector<decimate::math::Vec3f>& coords = result.positions;
    std::vector<int>&         faces  = result.indices;

    // RAII file handle.
    std::ifstream in{std::string(path)};
    if (!in) {
        throw std::runtime_error(
            std::string("loadObj: unable to open file '") + std::string(path) + "'");
    }

    // Read line by line, dispatch on the first character (FROZEN BUG 1:
    // 'v' only, not "v " — so 'vn'/'vt' also route to readVertex3).
    // strtok mutates its buffer, so we tokenize the line string's own storage
    // (std::string::data() is mutable and null-terminated since C++17) — no
    // fixed char[] copy, hence no 2047-char truncation.
    for (std::string row; std::getline(in, row);) {
        switch (row.empty() ? '\0' : row[0]) {
            case 'v':
                // FROZEN BUG 1: 'vn'/'vt' lines (row[0]=='v') append their XYZ
                // as phantom positions.  DO NOT add "v "/"vn " keyword matching.
                readVertex3(row.data(), coords);
                break;
            case 'f':
                readFace3(row.data(), faces);
                break;
            default:
                break;
        }
    }

    // Trust boundary: reject indices outside [0, positions.size()) so malformed
    // input (0, negative/relative, past-the-end) throws instead of reaching core.
    // Phantom vn/vt positions (FROZEN BUG 1) count toward the size, as before.
    for (const int idx : faces) {
        if (idx < 0 || static_cast<std::size_t>(idx) >= coords.size()) {
            throw std::runtime_error("loadObj: face index " + std::to_string(idx + 1) +
                                     " out of range (vertex count " +
                                     std::to_string(coords.size()) + ")");
        }
    }

    // No memcpy / new[] needed: we accumulated directly into result.positions
    // and result.indices (eliminating the legacy double-allocation, hot edge 4).
    // The numeric content is identical to legacy.
    return result;
}

}  // namespace decimate
