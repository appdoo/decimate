// This file is part of the decimate library; see decimate.h for version/license details.
// Standalone consumer smoke-test for an INSTALLED decimate package.
// Proves: find_package(decimate) works, the namespaced targets link, and the
// prefixed <decimate/...> public headers resolve from the install prefix.
#include <cstdio>

#include <decimate/types/mesh_data.h> // prefixed public header (install-tree path)
#include <decimate.h>                 // version smoke-test symbol (umbrella target)

/// Builds a one-triangle MeshData and prints the version: success means the install tree links.
int main() {
    decimate::MeshData mesh;    // value type from decimate::types
    mesh.indices = { 0, 1, 2 }; // one triangle
    std::printf("decimate consumer OK — version=%d, tris=%zu\n",
                decimate::version(), mesh.indices.size() / 3);
    return 0;
}
