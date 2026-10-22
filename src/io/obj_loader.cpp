// This file is part of the decimate library; see decimate.h for version/license details.
/// @file obj_loader.cpp
/// @brief OBJ mesh loader — behavior-preserving port of legacy ObjLdr::load.
///
/// Ported from legacy2003/loader/ObjLdr.cpp (Spec 06, Beat 2a).
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
///   FROZEN BUG 1 (legacy2003/loader/ObjLdr.cpp:83-84, row[0]=='v' only dispatch):
///     'vn' lines begin with 'v', so they route to read_vertex3.
///     Their XYZ normal components are appended as phantom vertex positions.
///     DO NOT add full-keyword matching here — that is the Beat-2b fix.
///
///   FROZEN BUG 2 (legacy2003/loader/ObjLdr.cpp:24 in read_face3, atoi stops at '/'):
///     Face tokens like "251//255" → atoi("251//255") == 251.
///     Only the vertex index survives; vt/vn attributes are silently discarded.
///     DO NOT split on '/' to extract attribute indices — that is the Beat-2b fix.
///
/// Tokenizing is std::string_view and number conversion is parse_number.h, so the
/// loader keeps no hidden state (the 2003 strtok cursor) and ignores the process
/// locale (the 2003 atof/atoi).

#include <decimate/types/mesh_data.h>

#include "parse_number.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace decimate {

/// Face-vertex cap: legacy2003/loader/ObjLdr.cpp:14 (read_face3) stores `int indices[10]`.
/// Preserved (no model asset has a >10-gon; raising it is a Beat-2b decision).
constexpr int kMaxFaceVerts = 10;


/// Pops the next token (a run of non-space, non-tab characters) off the front of
/// @p rest.  Returns an empty view when none is left, which is how strtok reported
/// "no more tokens" (a token is never empty).  Stateless, unlike strtok.
static std::string_view nextToken(std::string_view &rest) {
    const std::size_t begin = rest.find_first_not_of(" \t");
    if (begin == std::string_view::npos) {
        rest = {};
        return {};
    }
    rest.remove_prefix(begin);

    const std::size_t      len = std::min(rest.find_first_of(" \t"), rest.size());
    const std::string_view token = rest.substr(0, len);
    rest.remove_prefix(len);
    return token;
}


/// Port of ObjLdr::read_vertex3 (legacy2003/loader/ObjLdr.cpp:50): skips the keyword token
/// ("v", "vn", ...) and reads three whitespace-separated floats; missing ones are 0.
///
/// Called for every line whose first character is 'v' — including 'vn' and 'vt' lines
/// (FROZEN BUG 1: row[0]-only dispatch in load()).
static void readVertex3(std::string_view line, std::vector<math::Vec3f> &coords) {
    nextToken(line); // skip the keyword token

    math::Vec3f v;
    for (int i = 0; i < 3; ++i) {
        // A missing token is empty and parses as 0.
        v[i] = io_detail::parseFloat(nextToken(line));
    }

    coords.push_back(v);
}


/// Port of ObjLdr::read_face3 (legacy2003/loader/ObjLdr.cpp:14): reads up to kMaxFaceVerts
/// whitespace-separated index tokens, converts each via atoi, decrements by 1 (1-based to
/// 0-based), then fan-triangulates n-gons. Lines with fewer than 3 indices are dropped.
///
/// FROZEN BUG 2: atoi("251//255") == 251 — slash-separated refs are silently
/// truncated; only the vertex index survives.  DO NOT split on '/'.
static void readFace3(std::string_view line, std::vector<int> &faces) {
    std::array<int, kMaxFaceVerts> indices {};
    nextToken(line); // skip the 'f' keyword

    int count = 0;
    for (std::string_view token = nextToken(line);
         !token.empty() && count < kMaxFaceVerts;
         token = nextToken(line)) {
        // FROZEN BUG 2: the integer parse stops at '/', discarding vt/vn attribute indices.
        indices[static_cast<std::size_t>(count++)] = io_detail::parseInt(token) - 1;
    }

    // Fewer than 3 corners is not a polygon.
    if (count < 3) {
        return;
    }

    // First triangle, then fan-triangulate the rest (n-gon → n-2 triangles).
    const int anchor = indices[0];
    int       prev = indices[1];
    for (int i = 2; i < count; ++i) {
        faces.push_back(anchor);
        faces.push_back(prev);
        faces.push_back(indices[static_cast<std::size_t>(i)]);
        prev = indices[static_cast<std::size_t>(i)];
    }
}


/// Port of ObjLdr::load (legacy2003/loader/ObjLdr.cpp:66); contract in mesh_data.h.
MeshData loadObj(std::string_view path) {
    MeshData                  result;
    std::vector<math::Vec3f> &coords = result.positions;
    std::vector<int>         &faces = result.indices;

    std::ifstream in { std::string(path) };
    if (!in) {
        throw std::runtime_error(
            std::string("loadObj: unable to open file '") + std::string(path) + "'");
    }

    // Read line by line, dispatch on the first character (FROZEN BUG 1:
    // 'v' only, not "v " — so 'vn'/'vt' also route to readVertex3).
    // Lines are tokenized as string_views of `row` — no copy, no fixed char[]
    // buffer, hence no 2047-char truncation.
    for (std::string row; std::getline(in, row);) {
        switch (row.empty() ? '\0' : row[0]) {
            case 'v':
                // FROZEN BUG 1: 'vn'/'vt' lines (row[0]=='v') append their XYZ
                // as phantom positions.  DO NOT add "v "/"vn " keyword matching.
                readVertex3(row, coords);
                break;
            case 'f':
                readFace3(row, faces);
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

} // namespace decimate
