// This file is part of the decimate library; see decimate.h for version/license details.
/// @file main.cpp
/// @brief mesh-viewer entry point — Interactive QEM collapse viewer.
///
/// Usage:
///   mesh-viewer [path.obj|path.iv]
///   mesh-viewer --frames N [path.obj]   headless smoke: render N frames, exit 0
///   mesh-viewer --shot <outdir> [--shot-size WxH] [path]   write each panel as a TGA, exit 0
///
/// Default mesh: spot.obj from the model/ directory (MESH_MODEL_DIR definition).
///
/// Keys:  - / down  decimate   + / up  grow back   r reset   w wireframe
///        k  cycle reducer (QEM -> octree clustering -> halfedge)   q quit
/// Mouse: left-drag orbit   scroll wheel zoom

#include "viewer.h"

#ifdef STACKTRACE_ENABLED
#  include "Stacktrace.h"
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#ifndef MESH_MODEL_DIR
#  define MESH_MODEL_DIR "."
#endif

// NOLINTNEXTLINE(bugprone-exception-escape): main may propagate unhandled exceptions
int main(int argc, char* argv[]) {
#ifdef STACKTRACE_ENABLED
    diag::Stacktrace::installCrashHandler();
#endif

    std::string meshPath = std::string(MESH_MODEL_DIR) + "/spot.obj";
    int smokeFrames = 0;
    ui::ShotOptions shot;
    bool doShot = false;

    // Parse arguments: [--frames N] [path]
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--frames") == 0 && i + 1 < argc) {
            smokeFrames = std::atoi(argv[++i]);  // NOLINT(cert-err34-c): simple parse
        } else if (std::strcmp(argv[i], "--shot") == 0 && i + 1 < argc) {
            doShot       = true;
            shot.outDir  = argv[++i];
        } else if (std::strcmp(argv[i], "--shot-size") == 0 && i + 1 < argc) {
            int w = 0, h = 0;
            if (std::sscanf(argv[++i], "%dx%d", &w, &h) != 2 || w < 16 || h < 16 || w > 8192 || h > 8192) {
                std::fprintf(stderr, "mesh-viewer: --shot-size expects WxH (16..8192)\n");
                return 2;
            }
            shot.width  = w;
            shot.height = h;
        } else {
            meshPath = argv[i];
        }
    }

    // Also accept MESH_VIEWER_FRAMES env var for headless CI
    if (smokeFrames == 0) {
        const char* envFrames = std::getenv("MESH_VIEWER_FRAMES");  // NOLINT(concurrency-mt-unsafe)
        if (envFrames != nullptr) {
            smokeFrames = std::atoi(envFrames);  // NOLINT(cert-err34-c)
        }
    }

    return ui::runViewer(meshPath, smokeFrames, doShot ? &shot : nullptr);
}
