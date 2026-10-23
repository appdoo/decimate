// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once
/// @file viewer.h
/// @brief Interactive GLFW/OpenGL collapse viewer — Milestone 2 UI layer.
///
/// Owns the GLFW window (RAII via GlfwWindow), the orbit camera, and the
/// two MeshData buffers (original + current).  All algorithm calls go through
/// decimate::load / decimate::reduce from mesh-core; no geometry math lives here.
///
/// Boundary contract: this file #includes GL/GLFW headers freely.
/// mesh-core headers MUST NOT be included from any core/ translation unit
/// via this path — the dependency is one-way: ui → core, never core → ui.

#include <string>

// Forward-declare GLFW types so callers don't need the GLFW header
// just to talk to Viewer.
struct GLFWwindow;

namespace ui {

/// @brief RAII owner of a GLFW window + OpenGL context.
///
/// Destructor calls glfwDestroyWindow then glfwTerminate.
/// Non-copyable; move is deleted too — one context per process.
struct GlfwWindow {
    /// Native window; null until `open` succeeded. Callers use it for the GLFW callbacks and title.
    GLFWwindow *handle = nullptr;

    GlfwWindow() noexcept = default;
    ~GlfwWindow();

    GlfwWindow(const GlfwWindow &) = delete;
    GlfwWindow &operator=(const GlfwWindow &) = delete;
    GlfwWindow(GlfwWindow &&) = delete;
    GlfwWindow &operator=(GlfwWindow &&) = delete;

    /// @brief Create a window of the given size with an OpenGL 3.3 core-profile context.
    /// @returns true on success; logs error to stderr on failure.
    [[nodiscard]] bool open(int width, int height, const char *title) noexcept;

    /// @returns true once the user asked to close the window (Esc/q or the close button). Requires a window.
    [[nodiscard]] bool shouldClose() const noexcept;

    /// Presents the frame and pumps GLFW events, which fires the key/mouse callbacks. Requires a window.
    void swapAndPoll() noexcept;
};

/// @brief `--shot` request: write each viewer panel as its own TGA into `outDir`, then exit.
struct ShotOptions {
    std::string outDir;
    int         width = 1600; ///< Pixel size of the 3D panel images (atlas images are square, min(w,h)).
    int         height = 1000;
    bool        tangent = false;    ///< Bake and show the shots in tangent space (MikkTSpace) instead of object space.
    bool        straighten = false; ///< Also write the seam-straightened `_straightened` atlas/bake panels.
};

/// @brief Run the interactive viewer loop.
///
/// @param meshPath   Path to the mesh file to load (.obj or .iv).
/// @param smokeFrames  If > 0 run that many frames offscreen then exit 0.
///                     Used by the headless smoke test (CI / this-env check).
/// @param shot       If non-null: run the smoke pipeline, write the panel images, exit
///                   (0 on success). Takes precedence over `smokeFrames`.
/// @param lscm       Flatten the charts with LSCM alone (`--lscm`) instead of the default bijective
///                   flattener.
/// @param maxCharts  Segmenter chart cap (`--max-charts N`, default 16); 1 while the stretch bound is on.
/// @param maxStretchL2    Per-chart L2 stretch bound (`--max-stretch L2[,Linf]`, default off); 0 = off.
/// @param maxStretchLinf  Per-chart worst-triangle stretch bound (default off); 0 = off.  With either bound
///                        on, over-bound charts are cut open (Spec 20).
/// @returns Exit code: 0 on success, non-zero on load failure.
int runViewer(const std::string &meshPath, int smokeFrames, const ShotOptions *shot = nullptr,
              bool lscm = false, int maxCharts = 16, double maxStretchL2 = 0.0, double maxStretchLinf = 0.0);

} // namespace ui
