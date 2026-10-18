// This file is part of the decimate library; see decimate.h for version/license details.
/// @file load.cpp
/// @brief Extension-based dispatcher: routes .obj → loadObj, .iv → loadInventor.
///
/// Mirrors the dispatch logic in the 2009 main.cpp:134-137.
/// Extension comparison is case-insensitive (e.g. ".OBJ" and ".IV" are accepted).

#include <decimate/types/mesh_data.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <string_view>

namespace decimate {

namespace {

/// Returns the file extension (including the dot), lower-cased.
std::string lowerExtension(std::string_view path) {
    auto pos = path.rfind('.');
    if (pos == std::string_view::npos) {
        return {};
    }
    std::string ext{path.substr(pos)};
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return ext;
}

}  // namespace

MeshData load(std::string_view path) {
    const std::string ext = lowerExtension(path);
    if (ext == ".obj") {
        return loadObj(path);
    }
    if (ext == ".iv") {
        return loadInventor(path);
    }
    throw std::runtime_error(
        std::string("load: unrecognised file extension '") + ext +
        "' for path '" + std::string(path) + "'");
}

}  // namespace decimate
