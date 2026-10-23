// This file is part of the decimate library; see decimate.h for version/license details.
/// @file main.cpp
/// @brief mesh-viewer entry point — Interactive QEM collapse viewer.
///
/// Usage:
///   mesh-viewer [path.obj|path.iv]
///   mesh-viewer --frames N [path.obj]   headless smoke: render N frames, exit 0
///   mesh-viewer --shot <outdir> [--shot-size WxH] [--straighten] [--tangent] [path]   write each panel as a TGA, exit 0
///   mesh-viewer --max-charts N          segmenter chart cap (default 16; 1 while the stretch bound is on)
///   mesh-viewer --max-stretch L2[,Linf] per-chart stretch bound (default off; 0 = off): over-bound charts are cut open (Gu 2002)
///
/// Default mesh: spot.obj from the model/ directory (MESH_MODEL_DIR definition).
///
/// Keys:  - / down  decimate   + / up  grow back   r reset   w wireframe
///        k  cycle reducer (QEM -> octree clustering -> halfedge)   q quit
///        e  cycle gutter (dilate / seam-fill / erased)   x  toggle seam straightening
///        l  cycle seam-aware LOD (full / 50% / 10%) on the erased texture
/// Mouse: left-drag orbit   scroll wheel zoom

#include "viewer.h"

#ifdef STACKTRACE_ENABLED
#include "Stacktrace.h"
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#ifndef MESH_MODEL_DIR
#define MESH_MODEL_DIR "."
#endif


/// Parses the command line (see the file header for the flags) and hands off to ui::runViewer.
// NOLINTNEXTLINE(bugprone-exception-escape): main may propagate unhandled exceptions
int main(int argc, char *argv[]) {
#ifdef STACKTRACE_ENABLED
    diag::Stacktrace::installCrashHandler();
#endif

    std::string     meshPath = std::string(MESH_MODEL_DIR) + "/spot.obj";
    int             smokeFrames = 0;
    ui::ShotOptions shot;
    bool            doShot = false;
    bool            lscm = false;
    int             maxCharts = 16;
    double          maxStretchL2 = 0.0;
    double          maxStretchLinf = 0.0;

    // A bare argument is the mesh path; later ones override earlier ones.
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--frames") == 0 && i + 1 < argc) {
            smokeFrames = std::atoi(argv[++i]); // NOLINT(cert-err34-c): simple parse
        } else if (std::strcmp(argv[i], "--shot") == 0 && i + 1 < argc) {
            doShot = true;
            shot.outDir = argv[++i];
        } else if (std::strcmp(argv[i], "--lscm") == 0) {
            lscm = true;
        } else if (std::strcmp(argv[i], "--max-charts") == 0 && i + 1 < argc) {
            maxCharts = std::atoi(argv[++i]); // NOLINT(cert-err34-c): simple parse
            if (maxCharts < 1) {
                std::fprintf(stderr, "mesh-viewer: --max-charts expects N >= 1\n");
                return 2;
            }
        } else if (std::strcmp(argv[i], "--max-stretch") == 0 && i + 1 < argc) {
            // "L2" or "L2,Linf"; a missing Linf is 4, 0 switches a bound off ("0" alone: both).
            const int n = std::sscanf(argv[++i], "%lf,%lf", &maxStretchL2, &maxStretchLinf);
            if (n == 1) {
                maxStretchLinf = maxStretchL2 == 0.0 ? 0.0 : 4.0;
            }
            if (n < 1 || maxStretchL2 < 0.0 || maxStretchLinf < 0.0) {
                std::fprintf(stderr, "mesh-viewer: --max-stretch expects L2[,Linf] (>= 0, 0 = off)\n");
                return 2;
            }
        } else if (std::strcmp(argv[i], "--tangent") == 0) {
            shot.tangent = true;
        } else if (std::strcmp(argv[i], "--straighten") == 0) {
            shot.straighten = true;
        } else if (std::strcmp(argv[i], "--shot-size") == 0 && i + 1 < argc) {
            int w = 0, h = 0;
            if (std::sscanf(argv[++i], "%dx%d", &w, &h) != 2 || w < 16 || h < 16 || w > 8192 || h > 8192) {
                std::fprintf(stderr, "mesh-viewer: --shot-size expects WxH (16..8192)\n");
                return 2;
            }
            shot.width = w;
            shot.height = h;
        } else {
            meshPath = argv[i];
        }
    }

    // CI cannot always edit the command line, so the env var is the fallback for headless frames.
    if (smokeFrames == 0) {
        const char *envFrames = std::getenv("MESH_VIEWER_FRAMES"); // NOLINT(concurrency-mt-unsafe)
        if (envFrames != nullptr) {
            smokeFrames = std::atoi(envFrames); // NOLINT(cert-err34-c)
        }
    }

    return ui::runViewer(meshPath, smokeFrames, doShot ? &shot : nullptr, lscm, maxCharts, maxStretchL2, maxStretchLinf);
}
