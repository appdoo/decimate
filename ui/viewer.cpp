// This file is part of the decimate library; see decimate.h for version/license details.
/// @file viewer.cpp
/// @brief Interactive GLFW/OpenGL core-profile 3.3 2×2 panel viewer.
///
/// A VBO/VAO + GLSL 330 pipeline (the GL 2.1 immediate-mode path is gone).
/// Four shader programs (`Programs`) cover the four panels:
///   mesh   — panels ① and ②: lit mesh with a per-vertex colour attribute.
///   ortho  — panel ③: 2D ortho UV atlas (colour + line drawing).
///   quad   — panel ③ in 'n' mode: the baked normal map as a textured quad.
///   baked  — panel ④: textured mesh with tangent-space normal-map decode.
///
/// ## 2×2 Panel Layout (one window, four glViewport regions)
///
///   ┌────────────────┬────────────────┐
///   │ ① Decimated   │ ② Chart segs  │  top row
///   │   mesh 3D     │   3D by chart  │
///   ├────────────────┼────────────────┤
///   │ ③ UV atlas    │ ④ baked       │  bottom row
///   │   [0,1]² 2D   │   normal map   │
///   └────────────────┴────────────────┘
///
/// Panel ④ tangent decode:
///   The bake stores normals in the per-face TBN frame: X=snrm·T, Y=snrm·B, Z=snrm·N.
///   (T=Gram-Schmidt UV-gradient tangent, N=interpolated low normal, B=N×T.)
///   To reconstruct: decode texture → (dx,dy,dz), then world = dx*T + dy*B + dz*N.
///   For object-space (bakeTangentSpace=false) the stored normal is in world/object
///   space directly — display as-is (mapped to RGB).
///
/// NOLINT annotations at GL C-API sites where clang-tidy flags C-array usage or
/// C-style casts inside GL macros — the GL API gives us no choice there.

#define GL_SILENCE_DEPRECATION

// Core-profile GL 3.3 on macOS requires <OpenGL/gl3.h>, NOT <OpenGL/gl.h>.
// <OpenGL/gl.h> only exposes the legacy (compatibility) profile symbols.
#ifdef __APPLE__
#include <OpenGL/gl3.h> // NOLINT(modernize-deprecated-headers)
#else
#include <GL/gl.h> // NOLINT(modernize-deprecated-headers)
#endif

#include <GLFW/glfw3.h>

#include "viewer.h"

// mesh-core public API — the only algorithm surface this file touches.
#include <decimate/chart/chart.h>
#include <decimate/pipeline.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normal_map.h>
#include <decimate/bake/normals.h>
#include <decimate/chart/parametrize.h>
#include <decimate/seam/erase.h>
#include <decimate/seam/straighten.h>
#include <decimate/seam/seam_energy.h>
#include <decimate/seam/seam_topology.h>
#include <decimate/reduce/reducer.h>
#include <decimate/reduce/cluster_reduce.h>
#include <decimate/reduce/halfedge_reduce.h>
#include <decimate/reduce/seam_aware_reduce.h>

// UI-layer I/O — TGA writer (no GL dependency).
#include "shot_image.h"
#include "tga_writer.h"

#include <decimate/math/vec2.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

/// File-local shorthand (an alias, not a using-directive): this file is outside namespace decimate.
namespace math = decimate::math;

namespace ui {

namespace {

/// Plain float triple for flat normals and the bounding box (the viewer's own, not decimate::math's).
struct Vec3f {
    float x {}, y {}, z {};
};


/// Axis-aligned bounding box of the loaded mesh; sizes the initial orbit distance.
struct Aabb {
    Vec3f lo, hi;
};


/// Linear RGB tint in [0,1] for chart colours and overlay ink.
struct Rgb {
    float r, g, b;
};


/// Atlas UV point.
using Vec2f = math::Vec2f;


/// Column-major 4x4 matrix (OpenGL layout): just enough for perspective, ortho and the orbit camera.
struct GlMat4 {
    float m[16] {}; // NOLINT(modernize-avoid-c-arrays)

    static GlMat4 identity() noexcept {
        GlMat4 r;
        r.m[0] = r.m[5] = r.m[10] = r.m[15] = 1.0f;
        return r;
    }

    // column-major indexing: m[col*4 + row]
    float &at(int col, int row) noexcept {
        return m[col * 4 + row];
    }
    [[nodiscard]] float at(int col, int row) const noexcept {
        return m[col * 4 + row];
    }
};


/// The four linked GLSL programs; built once, after the GL context is current.
struct Programs {
    unsigned int mesh = 0;  // panels ① ②: kMeshVS + kMeshFS
    unsigned int ortho = 0; // panel ③ atlas: kOrthoVS + kOrthoFS
    unsigned int quad = 0;  // panel ③ normal-map preview: kQuadVS + kQuadFS
    unsigned int baked = 0; // panel ④: kBakedVS + kBakedFS
};


/// Everything the GLFW callbacks and the panel renderers share: meshes, camera, chartification, bake and LOD state.
struct ViewerState {
    // Meshes
    decimate::MeshData original;
    decimate::MeshData current;
    int                originalVertexCount {};
    int                currentTarget {};

    // Camera
    float azimuth {};
    float elevation { 30.0f };
    float distance {};
    float centreX {}, centreY {}, centreZ {};

    // Interaction
    bool   leftDragging { false };
    double lastMouseX {}, lastMouseY {};

    // Display
    bool wireframe { false };

    // Reducer strategy, cycled with 'K': 0 = QEM heap edge-collapse, 1 = octree
    // clustering (fast/sloppy), 2 = HALFEDGE corner-table edge-collapse (best).
    int reducer { 0 };

    // Chartification: charts, packed UVs and the charted mesh (see decimate::Atlas)
    decimate::Atlas atlas;
    bool            charted { false };

    // Baked normal map
    decimate::bake::NormalMap bakedMap;
    unsigned int              bakedTex = 0;
    bool                      baked = false;

    bool atlasShowNormalMap = false;
    bool bakeTangentSpace = false;
    bool panel4LitGray = false; // panel ④: lit grayscale vs abs-normal RGB
    /// Which single map `bakedMap` holds (see decimate::Gutter); cycled with 'E'.
    decimate::Gutter gutterMode = decimate::Gutter::Seamless;
    /// Seam straightening (Liu et al. 2017, section 6) between parametrize and bake; toggled with 'X'.
    bool straighten = false;
    /// Flattener of every chart (`--lscm` selects LSCM; default bijective).
    decimate::chart::Flattener flattener = decimate::chart::Flattener::Bijective;
    /// Segmenter chart cap (`--max-charts N`).
    int maxCharts = 16;
    /// Stretch bound per chart (`--max-stretch L2[,Linf]`, Spec 20); 0 = off.
    double maxStretchL2 = 0.0;
    double maxStretchLinf = 0.0;

    /// One seam-aware LOD (Liu et al. 2017, section 5), ready for `drawBakedMesh`.
    struct LodView {
        decimate::MeshData       mesh;         ///< positions + surface indices.
        std::vector<math::Vec2f> uv;           ///< Per-CORNER UV (3*faceCount).
        int                      vertices = 0; ///< Achieved vertex count (>= the target).
    };
    /// Panel-4 level of detail, cycled with 'L': 0 = full, 1 = 50%, 2 = 10% of the vertices.
    /// The LODs share the ONE uploaded texture: cycling never re-bakes.
    int                    lodLevel = 0;
    bool                   lodValid = false; ///< `lod` matches the current chartified mesh + UVs.
    std::array<LodView, 2> lod;
};


/// Projection, model-view and their product for one 3D panel.
struct CameraMatrices {
    GlMat4 proj;
    GlMat4 mv;
    GlMat4 mvp;
};


/// Offscreen colour+depth target. Rendering here instead of the (hidden) window
/// gives exactly the requested pixel size, independent of Retina scaling.
struct ShotTarget {
    unsigned int fbo {}, color {}, depth {};
    int          w {}, h {};
};


/// What the lit baked panel looks like for the current gutter mode, at frame size and at 3x
/// (the 3x render is where a seam close-up is cropped from).
struct LitShots {
    ui::Image frame;
    ui::Image big;
};


/// Triangle soup (x, y, r, g, b per vertex) for the 2D shot overlays; lines are quads whose
/// width is given in atlas-UV units (core-profile GL has no wide lines).
struct ShotBatch {
    std::vector<float> v;

    void vert(Vec2f p, Rgb c) {
        v.insert(v.end(), { p.x, p.y, c.r, c.g, c.b });
    }
    void tri(Vec2f a, Vec2f b, Vec2f c, Rgb col) {
        vert(a, col);
        vert(b, col);
        vert(c, col);
    }

    /// Segment a-b, `w` wide, square-capped so consecutive segments leave no gap at a joint.
    void seg(Vec2f a, Vec2f b, float w, Rgb c) {
        const float dx = b.x - a.x, dy = b.y - a.y;
        const float len = std::sqrt(dx * dx + dy * dy);
        if (len < 1e-9f) {
            return;
        }
        const float hx = dx / len * (w * 0.5f), hy = dy / len * (w * 0.5f);
        const Vec2f p0(a.x - hx, a.y - hy), p1(b.x + hx, b.y + hy);
        const Vec2f n(-hy, hx);
        const Vec2f q0(p0.x + n.x, p0.y + n.y), q1(p0.x - n.x, p0.y - n.y);
        const Vec2f q2(p1.x - n.x, p1.y - n.y), q3(p1.x + n.x, p1.y + n.y);
        tri(q0, q1, q2, c);
        tri(q0, q2, q3, c);
    }
    void loop(const std::vector<Vec2f> &pts, float w, Rgb c) {
        for (std::size_t i = 0; i < pts.size(); ++i) {
            seg(pts[i], pts[(i + 1) % pts.size()], w, c);
        }
    }
    void square(Vec2f p, float w, Rgb c) {
        const float h = w * 0.5f;
        const Vec2f a(p.x - h, p.y - h), b(p.x + h, p.y - h), cc(p.x + h, p.y + h), d(p.x - h, p.y + h);
        tri(a, b, cc, c);
        tri(a, cc, d, c);
    }
};


/// The sharpest convex gutterLoop corner of a chart, with what the close-up needs to frame it.
struct SharpCorner {
    int   chart = -1;
    int   vertex = -1; ///< index into that chart's gutterLoop
    float angleDeg = 180.0f;
    Vec2f nearest { 0.0f, 0.0f }; ///< the chart corner UV closest to that gutter vertex (its boundary vertex)
    float dist = 0.0f;            ///< distance to it (the realised miter offset)
};

} // namespace


/// Routes GLFW errors to stderr; GLFW calls this from inside glfwInit/glfwCreateWindow.
static void glfwErrorCallback(int code, const char *desc) noexcept {
    std::fprintf(stderr, "[GLFW error %d] %s\n", code, desc);
}


bool GlfwWindow::open(int width, int height, const char *title) noexcept {
    glfwSetErrorCallback(glfwErrorCallback);
    if (glfwInit() == GLFW_FALSE) {
        std::fprintf(stderr, "mesh-viewer: glfwInit failed\n");
        return false;
    }

    // Core-profile 3.3 — mandatory for GLSL 330 and VAO.
    // macOS also needs GLFW_OPENGL_FORWARD_COMPAT for the 3.3 core path.
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    handle = glfwCreateWindow(width, height, title, nullptr, nullptr);
    if (handle == nullptr) {
        std::fprintf(stderr, "mesh-viewer: glfwCreateWindow failed\n");
        glfwTerminate();
        return false;
    }

    glfwMakeContextCurrent(handle);
    glfwSwapInterval(1); // vsync
    return true;
}


GlfwWindow::~GlfwWindow() {
    if (handle != nullptr) {
        glfwDestroyWindow(handle);
    }
    glfwTerminate();
}


bool GlfwWindow::shouldClose() const noexcept {
    return glfwWindowShouldClose(handle) != 0;
}


void GlfwWindow::swapAndPoll() noexcept {
    glfwSwapBuffers(handle);
    glfwPollEvents();
}


/// Bounds of `positions`; an empty mesh gives the zero box.
[[nodiscard]] static Aabb computeAabb(const std::vector<math::Vec3f> &positions) noexcept {
    if (positions.empty()) {
        return { { 0, 0, 0 }, { 0, 0, 0 } };
    }
    Aabb box {
        { positions[0].x, positions[0].y, positions[0].z },
        { positions[0].x, positions[0].y, positions[0].z }
    };
    for (std::size_t i = 1; i < positions.size(); ++i) {
        const auto &v = positions[i];
        box.lo.x = std::min(box.lo.x, v.x);
        box.lo.y = std::min(box.lo.y, v.y);
        box.lo.z = std::min(box.lo.z, v.z);
        box.hi.x = std::max(box.hi.x, v.x);
        box.hi.y = std::max(box.hi.y, v.y);
        box.hi.z = std::max(box.hi.z, v.z);
    }
    return box;
}


/// Half the box diagonal: the radius of a sphere that encloses the mesh.
[[nodiscard]] static float aabbRadius(const Aabb &b) noexcept {
    const float dx = b.hi.x - b.lo.x;
    const float dy = b.hi.y - b.lo.y;
    const float dz = b.hi.z - b.lo.z;
    return 0.5f * std::sqrt(dx * dx + dy * dy + dz * dz);
}


/// Distinct tint per chart id; ids wrap around the 16-entry palette so any chart count gets a colour.
[[nodiscard]] static Rgb chartColor(int id) noexcept {
    static constexpr std::array<Rgb, 16> kPalette = { {
        { 0.90f, 0.20f, 0.20f }, // red
        { 0.20f, 0.65f, 0.90f }, // sky-blue
        { 0.20f, 0.80f, 0.35f }, // green
        { 0.95f, 0.80f, 0.10f }, // yellow
        { 0.75f, 0.30f, 0.90f }, // violet
        { 0.95f, 0.50f, 0.10f }, // orange
        { 0.15f, 0.75f, 0.75f }, // cyan
        { 0.90f, 0.20f, 0.70f }, // pink
        { 0.50f, 0.75f, 0.20f }, // lime
        { 0.20f, 0.35f, 0.85f }, // blue
        { 0.85f, 0.55f, 0.20f }, // amber
        { 0.55f, 0.20f, 0.65f }, // purple
        { 0.25f, 0.65f, 0.50f }, // teal
        { 0.80f, 0.80f, 0.20f }, // gold
        { 0.60f, 0.20f, 0.25f }, // maroon
        { 0.45f, 0.45f, 0.85f }, // periwinkle
    } };
    const int                            idx = std::abs(id) % static_cast<int>(kPalette.size());
    return kPalette[static_cast<std::size_t>(idx)];
}


/// a * b in column-major layout.
[[nodiscard]] static GlMat4 mat4Mul(const GlMat4 &a, const GlMat4 &b) noexcept {
    GlMat4 r;
    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 4; ++row) {
            float s = 0.0f;
            for (int k = 0; k < 4; ++k) {
                s += a.at(k, row) * b.at(col, k);
            }
            r.at(col, row) = s;
        }
    }
    return r;
}


/// Perspective matrix matching the old glFrustum(-right,right,-top,top,near,far).
[[nodiscard]] static GlMat4 makePerspective(float fovYDeg, float aspect, float nearP, float farP) noexcept {
    const float top = nearP * std::tan(fovYDeg * 0.5f * 3.14159265f / 180.0f);
    const float right = top * aspect;
    // glFrustum equivalent (column-major)
    GlMat4 m;
    m.at(0, 0) = nearP / right;
    m.at(1, 1) = nearP / top;
    m.at(2, 2) = -(farP + nearP) / (farP - nearP);
    m.at(2, 3) = -1.0f;
    m.at(3, 2) = -2.0f * farP * nearP / (farP - nearP);
    return m;
}


/// Ortho matching glOrtho(l,r,b,t,-1,1).
[[nodiscard]] static GlMat4 makeOrtho(float l, float r, float b, float t) noexcept {
    GlMat4 m = GlMat4::identity();
    m.at(0, 0) = 2.0f / (r - l);
    m.at(1, 1) = 2.0f / (t - b);
    m.at(2, 2) = -1.0f;
    m.at(3, 0) = -(r + l) / (r - l);
    m.at(3, 1) = -(t + b) / (t - b);
    return m;
}


/// Orbit model-view: translate(0,0,-dist) · rotX(elev) · rotY(azim) · translate(-cx,-cy,-cz).
[[nodiscard]] static GlMat4 makeOrbitMV(float azimDeg, float elevDeg, float dist,
                                        float cx, float cy, float cz) noexcept {
    const float az = azimDeg * 3.14159265f / 180.0f;
    const float el = elevDeg * 3.14159265f / 180.0f;
    const float saz = std::sin(az), caz = std::cos(az);
    const float sel = std::sin(el), cel = std::cos(el);

    // rotY(az): rotate azimuth around Y
    GlMat4 rY = GlMat4::identity();
    rY.at(0, 0) = caz;
    rY.at(2, 0) = saz;
    rY.at(0, 2) = -saz;
    rY.at(2, 2) = caz;

    // rotX(el): tilt around X
    GlMat4 rX = GlMat4::identity();
    rX.at(1, 1) = cel;
    rX.at(2, 1) = sel;
    rX.at(1, 2) = -sel;
    rX.at(2, 2) = cel;

    // translate(-cx, -cy, -cz)
    GlMat4 tC = GlMat4::identity();
    tC.at(3, 0) = -cx;
    tC.at(3, 1) = -cy;
    tC.at(3, 2) = -cz;

    // translate(0, 0, -dist)
    GlMat4 tD = GlMat4::identity();
    tD.at(3, 2) = -dist;

    // MV = tD · rX · rY · tC
    return mat4Mul(tD, mat4Mul(rX, mat4Mul(rY, tC)));
}


/// Compiles one shader stage; a compile error is logged but the (invalid) handle is still returned.
static unsigned int compileShader(GLenum type, const char *src) {
    const unsigned int sh = glCreateShader(type);
    glShaderSource(sh, 1, &src, nullptr);
    glCompileShader(sh);
    int ok = 0;
    glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
    if (ok == GL_FALSE) {
        char log[1024];
        glGetShaderInfoLog(sh, sizeof(log), nullptr, log);
        std::fprintf(stderr, "[GL] shader compile error:\n%s\n", log);
    }
    return sh;
}


/// Links a vertex + fragment shader into a program and frees the two shader objects.
static unsigned int linkProgram(unsigned int vs, unsigned int fs) {
    const unsigned int prog = glCreateProgram();
    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glLinkProgram(prog);
    int ok = 0;
    glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (ok == GL_FALSE) {
        char log[1024];
        glGetProgramInfoLog(prog, sizeof(log), nullptr, log);
        std::fprintf(stderr, "[GL] program link error:\n%s\n", log);
    }
    glDeleteShader(vs);
    glDeleteShader(fs);
    return prog;
}


/// Mesh program (panels ① ②). Attributes: aPos(0), aNrm(1), aCol(2).
/// Headlight from view direction (0,0.6,1) normalised — constant, no matrix needed.

static const char *kMeshVS = R"(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNrm;
layout(location=2) in vec3 aCol;
uniform mat4 uMVP;
uniform mat4 uMV;   // model-view for transforming normals
out vec3 vNrm;
out vec3 vCol;
void main() {
    gl_Position = uMVP * vec4(aPos, 1.0);
    // Transform normal to view space (no non-uniform scale, so MV is fine).
    vNrm = mat3(uMV) * aNrm;
    vCol = aCol;
}
)";


static const char *kMeshFS = R"(
#version 330 core
in vec3 vNrm;
in vec3 vCol;
out vec4 fragColor;
// Headlight direction in view space (matches old GL_LIGHT0 position (0,0.6,1,0)).
const vec3 kLightDir = normalize(vec3(0.0, 0.6, 1.0));
const vec3 kAmbient  = vec3(0.20, 0.20, 0.20);
const vec3 kDiffuse  = vec3(0.85, 0.80, 0.75);
void main() {
    vec3 n = normalize(vNrm);
    float diff = max(dot(n, kLightDir), 0.0);
    vec3 lit = vCol * (kAmbient + kDiffuse * diff);
    fragColor = vec4(lit, 1.0);
}
)";


/// 2D ortho program (panel ③). Attributes: aPos2(0) as vec2, aCol(1) as vec3.
/// Used for both filled chart triangles and line edges.

static const char *kOrthoVS = R"(
#version 330 core
layout(location=0) in vec2 aPos2;
layout(location=1) in vec3 aCol;
uniform mat4 uOrtho;
out vec3 vCol;
void main() {
    gl_Position = uOrtho * vec4(aPos2, 0.0, 1.0);
    vCol = aCol;
}
)";


static const char *kOrthoFS = R"(
#version 330 core
in vec3 vCol;
out vec4 fragColor;
void main() {
    fragColor = vec4(vCol, 1.0);
}
)";


/// Textured-quad program for the baked-map preview (panel ③ 'n' mode). Kept apart from the
/// ortho program so that one needs no sampler.

static const char *kQuadVS = R"(
#version 330 core
layout(location=0) in vec2 aPos2;
layout(location=1) in vec2 aUV;
uniform mat4 uOrtho;
out vec2 vUV;
void main() {
    gl_Position = uOrtho * vec4(aPos2, 0.0, 1.0);
    vUV = aUV;
}
)";


static const char *kQuadFS = R"(
#version 330 core
in vec2 vUV;
uniform sampler2D uTex;
out vec4 fragColor;
void main() {
    fragColor = vec4(texture(uTex, vUV).rgb, 1.0);
}
)";


/// Baked-mesh program (panel ④). Attributes: aPos(0), aNrm(1), aUV(2), aTan(3).
/// uniform uTangentSpace: 1 = decode TBN, 0 = display object-space normal as RGB.
///
/// TBN decode matches mapFaceRt convention exactly:
///   Stored in texture: x=dot(snrm,T), y=dot(snrm,B), z=dot(snrm,N)
///   T = Gram-Schmidt UV-gradient tangent, N = interpolated low normal, B = N×T
///   Reconstruct: world = dx*T + dy*B + dz*N
///
/// Handedness: B = N×T (no hand sign needed; the bake uses the same cross product).

static const char *kBakedVS = R"(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNrm;
layout(location=2) in vec2 aUV;
layout(location=3) in vec3 aTan;
uniform mat4 uMVP;
uniform mat4 uMV;
out vec3 vNrm;
out vec3 vTan;
out vec2 vUV;
void main() {
    gl_Position = uMVP * vec4(aPos, 1.0);
    mat3 mv3 = mat3(uMV);
    vNrm = mv3 * aNrm;
    vTan = mv3 * aTan;
    vUV  = aUV;
}
)";


static const char *kBakedFS = R"(
#version 330 core
in vec3 vNrm;
in vec3 vTan;
in vec2 vUV;
uniform sampler2D uTex;
uniform bool uTangentSpace;   // true = TBN decode; false = object-space show-as-RGB
uniform bool uLitGray;        // true = lit grayscale (read the bump); false = abs-normal RGB
uniform mat4 uMV;             // same uniform as in the VS; rotates object-space normals to view space
out vec4 fragColor;
void main() {
    // Decode: map [0,255] → [-1,1] (rangeCompress inverse)
    vec3 t = texture(uTex, vUV).rgb * 2.0 - 1.0;

    vec3 worldNrm;
    if (uTangentSpace) {
        // Reconstruct TBN in view space (matches mapFaceRt: T=tangent, N=normal, B=N×T).
        vec3 N = normalize(vNrm);
        vec3 T = normalize(vTan - dot(vTan, N) * N);  // Gram-Schmidt re-orthogonalise
        vec3 B = cross(N, T);                          // same as bake: B = N×T
        // t = (dot(snrm,T), dot(snrm,B), dot(snrm,N)) → reconstruct
        worldNrm = normalize(t.x * T + t.y * B + t.z * N);
    } else {
        // Object-space: texture already holds the world normal directly.
        worldNrm = normalize(t);
        // The light below is a camera-space headlight (as in the TBN branch, whose
        // normals are already view-space). Lit mode must rotate the object-space normal
        // into view space too, else the fixed object-space light leaves most of the model black.
        if (uLitGray) {
            worldNrm = normalize(mat3(uMV) * worldNrm);
        }
    }

    const vec3 kLightDir = normalize(vec3(0.0, 0.6, 1.0));
    float diff = max(dot(worldNrm, kLightDir), 0.0);
    vec3 col;
    if (uLitGray) {
        // Lit grayscale: a plain shaded surface lit THROUGH the normal map —
        // best for actually reading the baked bump detail.
        col = vec3(0.15 + 0.85 * diff);
    } else {
        // abs-normal RGB: colour encodes which way each texel faces; makes the
        // "seam is gone" result obvious. Chromatic by construction (abs → saturated).
        col = abs(worldNrm) * (0.4 + 0.6 * diff);
    }
    fragColor = vec4(col, 1.0);
}
)";


/// Compiles and links all four programs; needs a current GL context.
static Programs buildPrograms() {
    Programs p;
    p.mesh = linkProgram(compileShader(GL_VERTEX_SHADER, kMeshVS),
                         compileShader(GL_FRAGMENT_SHADER, kMeshFS));
    p.ortho = linkProgram(compileShader(GL_VERTEX_SHADER, kOrthoVS),
                          compileShader(GL_FRAGMENT_SHADER, kOrthoFS));
    p.quad = linkProgram(compileShader(GL_VERTEX_SHADER, kQuadVS),
                         compileShader(GL_FRAGMENT_SHADER, kQuadFS));
    p.baked = linkProgram(compileShader(GL_VERTEX_SHADER, kBakedVS),
                          compileShader(GL_FRAGMENT_SHADER, kBakedFS));
    return p;
}


/// The viewer state the GLFW callbacks reach (they carry no user pointer); set by runViewer, null outside it.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static ViewerState *gState = nullptr;


static void bakeAndUpload(ViewerState &s);


/// @brief Print the straightening stats of `st`.
static void printStraighten(const decimate::seam::StraightenStats &st) {
    std::fprintf(stderr,
                 "[straighten] chains seam=%d boundary=%d | protected=%d (rules %d/%d/%d/%d) | "
                 "ears=%d (same-neighbour %d) | reverted charts=%d\n",
                 st.seamChains, st.boundaryChains, st.protectedVertices, st.protectedByRule[0],
                 st.protectedByRule[1], st.protectedByRule[2], st.protectedByRule[3], st.earFaces,
                 st.sameNeighbourEars, st.revertedCharts);
}


/// @brief Injectivity-guard knobs of the viewer: the flattener (`--lscm`) and the stretch bound.
static decimate::chart::InjectivityOptions injectivityOptions(const ViewerState &s) {
    decimate::chart::InjectivityOptions o;
    o.flattener = s.flattener;
    o.maxStretchL2 = s.maxStretchL2;
    o.maxStretchLinf = s.maxStretchLinf;
    return o;
}


/// @brief Chart cap handed to the segmenter: `--max-charts`, or one chart when the bound cuts it open.
static int segmenterCap(const ViewerState &s) {
    const bool bounded = s.maxStretchL2 > 0.0 || s.maxStretchLinf > 0.0;
    return bounded ? 1 : s.maxCharts;
}


/// @brief Chartify `state.current` (decimate::chartify), log the stats, then bake.
static void computeChartification(ViewerState &s) {
    if (s.current.indices.empty() || s.current.positions.empty()) {
        std::fprintf(stderr, "[chart] mesh is empty — nothing to chartify\n");
        return;
    }

    decimate::ChartOptions opt;
    opt.maxCharts = segmenterCap(s);
    opt.injectivity = injectivityOptions(s);
    opt.straighten = s.straighten;
    s.atlas = decimate::chartify(s.current, opt);

    s.charted = true;
    s.lodValid = false; // the LODs were built from the previous atlas
    s.lodLevel = 0;

    const auto &a = s.atlas;
    const auto &inj = a.stats;
    if (s.straighten) {
        printStraighten(a.straightenStats);
    }
    const int numCharts = static_cast<int>(a.charts.charts.size());
    const auto ms = [](double sec) { return static_cast<long>(sec * 1000.0); };
    std::fprintf(stderr, "[chart] %zu faces → %d charts, %zu UVs | hem=%ldms seg=%ldms param=%ldms total=%ldms\n",
                 s.current.indices.size() / 3, numCharts, a.uv.size(), ms(a.hemSeconds), ms(a.segmentSeconds),
                 ms(a.flattenSeconds), ms(a.hemSeconds + a.segmentSeconds + a.flattenSeconds));
    std::fprintf(stderr,
                 "[chart] flattener = %s | selfOverlapping before/after = %d/%d, non-injective after = %d, flattened bijectively = %d, "
                 "repaired = %d (%.2fs), splits = %d\n",
                 s.flattener == decimate::chart::Flattener::Lscm ? "lscm" : "bijective",
                 inj.selfOverlappingBefore, inj.selfOverlappingAfter, inj.nonInjectiveAfter, inj.flattenedBijectively,
                 inj.repairedBijectively, inj.bijectiveSeconds, inj.splits);
    std::fprintf(stderr,
                 "[chart] stretch bound L2 <= %g, Linf <= %g (0 = off) | %d charts, %d cut rounds, %d closed charts opened, "
                 "max chart L2 = %.3f, Linf = %.3f\n",
                 s.maxStretchL2, s.maxStretchLinf, numCharts, inj.cuts, inj.closedOpened, inj.maxChartL2,
                 inj.maxChartLinf);

    bakeAndUpload(s);
}


/// Gutter name for log lines.
static const char *gutterName(decimate::Gutter g) {
    switch (g) {
        case decimate::Gutter::Seamless:
            return "seamless";
        case decimate::Gutter::Dilate:
            return "dilate";
        case decimate::Gutter::SeamFill:
            return "seam-fill";
        case decimate::Gutter::Erased:
            return "erased";
    }
    return "?";
}


/// seamless -> dilate -> seam-fill -> erased -> seamless.
static decimate::Gutter nextGutterMode(decimate::Gutter g) {
    using decimate::Gutter;
    return g == Gutter::Seamless   ? Gutter::Dilate
           : g == Gutter::Dilate   ? Gutter::SeamFill
           : g == Gutter::SeamFill ? Gutter::Erased
                                   : Gutter::Seamless;
}


/// Print the seam mismatch D_total (paper eq. 9) per channel, in texel units [0,1].
static void printSeamEnergy(const char *label, const std::array<double, 3> &d) {
    std::fprintf(stderr, "[seam] D_total %-9s R=%.4e G=%.4e B=%.4e\n", label, d[0], d[1], d[2]);
}


/// @brief Bake the normal map (decimate::bakeAtlas) and upload it as a cached GL texture.
static void bakeAndUpload(ViewerState &s) {
    // Drop the old texture first: if this bake comes back empty, a stale map must not stay on screen.
    if (s.bakedTex != 0) {
        glDeleteTextures(1, &s.bakedTex);
        s.bakedTex = 0;
        s.baked = false;
    }

    decimate::BakeOptions opt;
    opt.space = s.bakeTangentSpace ? decimate::BakeSpace::Tangent : decimate::BakeSpace::Object;
    opt.gutter = s.gutterMode;
    auto r = decimate::bakeAtlas(s.atlas, s.original, opt);

    if (r.gutterUsed == decimate::Gutter::Erased) {
        printSeamEnergy("dilate", r.dDilate);
        std::fprintf(stderr, "[seam] D_total erased    R=%.4e G=%.4e B=%.4e  (%d unknown texels)\n", r.dTotal[0],
                     r.dTotal[1], r.dTotal[2], r.unknownTexels);
    } else {
        printSeamEnergy(gutterName(r.gutterUsed), r.dTotal);
    }
    s.bakedMap = std::move(r.map);
    std::fprintf(stderr, "[bake] %s-space normals, %s gutter\n",
                 s.bakeTangentSpace ? "tangent" : "object", gutterName(r.gutterUsed));

    if (s.bakedMap.pixels.empty()) {
        std::fprintf(stderr, "[bake] WARNING: empty NormalMap — skipping GL upload\n");
        return;
    }

    // Linear filtering without mips; clamped so lookups at the atlas border do not wrap around.
    unsigned int tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // rows are 3*width bytes, not 4-byte aligned
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
                 s.bakedMap.width, s.bakedMap.height, 0,
                 GL_RGB, GL_UNSIGNED_BYTE,
                 s.bakedMap.pixels.data());
    glBindTexture(GL_TEXTURE_2D, 0);

    if (const GLenum glErr = glGetError(); glErr != GL_NO_ERROR) {
        std::fprintf(stderr, "[bake] GL error 0x%x after glTexImage2D\n",
                     static_cast<unsigned>(glErr));
    }

    s.bakedTex = tex;
    s.baked = true;
}


/// TARGET percent of the full vertex count of LOD level 1 and 2.  The reducer stops above the
/// target when no legal collapse remains, so the title shows the achieved count too.
static constexpr int kLodPercent[2] = { 50, 10 }; // NOLINT(modernize-avoid-c-arrays)

/// Decimate the chartified atlas (`s.atlas.mesh` + `s.atlas.uv`) seam-aware, once, for both LOD
/// levels, and print what each level costs.  Seams only collapse on a straightened atlas
/// ('x'); on a plain one the LODs still keep the seam, they just cannot simplify it.
static void computeLods(ViewerState &s) {
    const decimate::TexturedMesh tm = decimate::seam::weld(s.atlas.mesh, s.atlas.uv);
    const int                    full = static_cast<int>(tm.positions.size());

    // Both levels decimate the same welded atlas mesh, so cycling LODs never re-bakes.
    for (std::size_t i = 0; i < s.lod.size(); ++i) {
        const auto res = decimate::reduceSeamAware(tm, full * kLodPercent[i] / 100);

        s.lod[i].mesh.positions = res.mesh.positions;
        s.lod[i].mesh.indices = res.mesh.faceV;
        s.lod[i].uv = res.mesh.cornerUV();
        s.lod[i].vertices = static_cast<int>(res.mesh.positions.size());

        const auto &st = res.stats;
        std::fprintf(stderr,
                     "[lod] target %d%%: verts %d -> %zu (achieved), tris %zu | seams %d -> %d (collapsed %d, "
                     "uncollapsible %d)%s\n",
                     kLodPercent[i], full, res.mesh.positions.size(), res.mesh.faceCount(),
                     st.seamEdgesBefore, st.seamEdgesAfter, st.collapsedSeamEdges,
                     st.uncollapsibleSeamEdges,
                     s.straighten ? "" : "  (seams are not straightened: press 'x')");

        // The energy of the shown texture on the reduced mesh tells how well the LOD keeps the seam.
        if (s.baked) {
            printSeamEnergy("lod", decimate::seam::seamEnergy(res.mesh, s.bakedMap));
        }
    }
    s.lodValid = true;
}


/// Key 'L': full -> 50% -> 10% LOD of the model in panel 4, drawn with the ONE uploaded
/// texture (the erased map): the first press switches to the erased gutter and bakes once.
static void cycleLod(ViewerState &s) {
    if (!s.charted) {
        std::fprintf(stderr, "[lod] press 'c' first to chartify and bake\n");
        return;
    }

    // LODs are drawn with the erased texture, so switch to it (one bake) before building them.
    if (!s.lodValid) {
        if (s.gutterMode != decimate::Gutter::Erased) {
            s.gutterMode = decimate::Gutter::Erased;
            bakeAndUpload(s);
        }
        computeLods(s);
    }

    s.lodLevel = (s.lodLevel + 1) % 3;
    if (s.lodLevel == 0) {
        std::fprintf(stderr, "[lod] full\n");
        return;
    }

    const auto &shown = s.lod[static_cast<std::size_t>(s.lodLevel - 1)];
    std::fprintf(stderr, "[lod] target %d%%, achieved %d verts\n", kLodPercent[s.lodLevel - 1], shown.vertices);
}


/// Upload float data into a fresh VBO, bind it to attribute `attrib`. (Currently unused: every draw helper
/// below packs one interleaved buffer instead.)
/// Each element has `stride` floats total; this attribute starts at `offset`.
/// Returns the VBO handle so the caller can delete it.
static unsigned int makeVboAttrib(unsigned int vao, int attrib,
                                  const float *data, int dataCount,
                                  int stride, int offset) {
    (void)vao; // VAO already bound by caller
    unsigned int vbo = 0;
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(dataCount * sizeof(float)),
                 data, GL_STREAM_DRAW);
    glEnableVertexAttribArray(static_cast<unsigned>(attrib));
    glVertexAttribPointer(static_cast<unsigned>(attrib), stride, GL_FLOAT, GL_FALSE,
                          stride * static_cast<int>(sizeof(float)),
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void *>(static_cast<std::intptr_t>(offset * sizeof(float))));
    return vbo;
}


/// Unit normal of triangle (a,b,c); a degenerate triangle returns its (near-)zero cross product unnormalised.
static Vec3f flatNormal(const math::Vec3f &a, const math::Vec3f &b, const math::Vec3f &c) noexcept {
    const float ex = b.x - a.x, ey = b.y - a.y, ez = b.z - a.z;
    const float fx = c.x - a.x, fy = c.y - a.y, fz = c.z - a.z;
    float       nx = ey * fz - ez * fy;
    float       ny = ez * fx - ex * fz;
    float       nz = ex * fy - ey * fx;
    const float len = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (len > 1e-8f) {
        nx /= len;
        ny /= len;
        nz /= len;
    }
    return { nx, ny, nz };
}


/// Draw the mesh with flat per-face normals and a uniform colour (panel ①).
/// Each vertex in the interleaved buffer: pos(3) + nrm(3) + col(3) = 9 floats.
static void drawMeshFlat(const Programs &prog,
                         const GlMat4 &mvp, const GlMat4 &mv,
                         const decimate::MeshData &mesh,
                         float cr, float cg, float cb) {
    const auto       &pos = mesh.positions;
    const auto       &idx = mesh.indices;
    const std::size_t numTris = idx.size() / 3;
    const auto        nv = static_cast<int>(pos.size());

    // Flat shading needs one normal per face, so the buffer is un-indexed: 3 vertices per triangle,
    // interleaved pos(3) nrm(3) col(3).
    std::vector<float> buf;
    buf.reserve(numTris * 3 * 9);

    for (std::size_t t = 0; t < numTris; ++t) {
        const int ia = idx[t * 3 + 0];
        const int ib = idx[t * 3 + 1];
        const int ic = idx[t * 3 + 2];
        // A corrupt index would read outside `positions`: drop the triangle instead.
        if (ia >= nv || ib >= nv || ic >= nv || ia < 0 || ib < 0 || ic < 0) {
            continue;
        }
        const auto &a = pos[static_cast<std::size_t>(ia)];
        const auto &b = pos[static_cast<std::size_t>(ib)];
        const auto &c = pos[static_cast<std::size_t>(ic)];
        const auto  n = flatNormal(a, b, c);

        for (const auto &v : { a, b, c }) {
            buf.push_back(v.x);
            buf.push_back(v.y);
            buf.push_back(v.z);
            buf.push_back(n.x);
            buf.push_back(n.y);
            buf.push_back(n.z);
            buf.push_back(cr);
            buf.push_back(cg);
            buf.push_back(cb);
        }
    }

    if (buf.empty()) {
        return;
    }

    const int stride = 9; // 3 pos + 3 nrm + 3 col
    const int vertCount = static_cast<int>(buf.size()) / stride;

    // Core profile has no client-side arrays: a throwaway VAO + VBO per draw is fine for a debug viewer
    // and needs no invalidation tracking when the mesh changes.
    unsigned int vao = 0;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    unsigned int vbo = 0;
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(buf.size() * sizeof(float)),
                 buf.data(), GL_STREAM_DRAW);

    // pos @ attrib 0, nrm @ attrib 1, col @ attrib 2
    const int sBytes = stride * static_cast<int>(sizeof(float));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sBytes, nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void *>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void *>(6 * sizeof(float)));

    glUseProgram(prog.mesh);
    glUniformMatrix4fv(glGetUniformLocation(prog.mesh, "uMVP"), 1, GL_FALSE, mvp.m);
    glUniformMatrix4fv(glGetUniformLocation(prog.mesh, "uMV"), 1, GL_FALSE, mv.m);

    glDrawArrays(GL_TRIANGLES, 0, vertCount);

    glBindVertexArray(0);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}


/// Draw mesh with per-face chart colours (panel ②).
static void drawMeshChartColors(const Programs &prog,
                                const GlMat4 &mvp, const GlMat4 &mv,
                                const decimate::MeshData        &mesh,
                                const decimate::chart::ChartSet &cs) {
    const auto       &pos = mesh.positions;
    const auto       &idx = mesh.indices;
    const std::size_t numTris = idx.size() / 3;
    const auto        nv = static_cast<int>(pos.size());

    // Same un-indexed layout as drawMeshFlat, with the colour taken from the face's chart.
    std::vector<float> buf;
    buf.reserve(numTris * 3 * 9);

    for (std::size_t t = 0; t < numTris; ++t) {
        const int ia = idx[t * 3 + 0];
        const int ib = idx[t * 3 + 1];
        const int ic = idx[t * 3 + 2];
        // A corrupt index would read outside `positions`: drop the triangle instead.
        if (ia >= nv || ib >= nv || ic >= nv || ia < 0 || ib < 0 || ic < 0) {
            continue;
        }
        const auto &a = pos[static_cast<std::size_t>(ia)];
        const auto &b = pos[static_cast<std::size_t>(ib)];
        const auto &c = pos[static_cast<std::size_t>(ic)];
        const auto  n = flatNormal(a, b, c);

        // A face beyond the chart map (stale chartSet) falls back to chart 0's colour.
        int chartId = 0;
        if (t < cs.faceToChart.size()) {
            chartId = cs.faceToChart[t];
        }
        const Rgb col = chartColor(chartId);

        for (const auto &v : { a, b, c }) {
            buf.push_back(v.x);
            buf.push_back(v.y);
            buf.push_back(v.z);
            buf.push_back(n.x);
            buf.push_back(n.y);
            buf.push_back(n.z);
            buf.push_back(col.r);
            buf.push_back(col.g);
            buf.push_back(col.b);
        }
    }

    if (buf.empty()) {
        return;
    }

    const int stride = 9;
    const int vertCount = static_cast<int>(buf.size()) / stride;
    const int sBytes = stride * static_cast<int>(sizeof(float));

    // Throwaway VAO + VBO per draw (see drawMeshFlat).
    unsigned int vao = 0, vbo = 0;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(buf.size() * sizeof(float)),
                 buf.data(), GL_STREAM_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sBytes, nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void *>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void *>(6 * sizeof(float)));

    glUseProgram(prog.mesh);
    glUniformMatrix4fv(glGetUniformLocation(prog.mesh, "uMVP"), 1, GL_FALSE, mvp.m);
    glUniformMatrix4fv(glGetUniformLocation(prog.mesh, "uMV"), 1, GL_FALSE, mv.m);

    glDrawArrays(GL_TRIANGLES, 0, vertCount);

    glBindVertexArray(0);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}


/// Draw UV atlas in ortho [0,1]² with chart colours + wireframe edges (panel ③).
// NOLINTNEXTLINE(readability-function-size): two separate draw calls (fill + lines)
static void drawUVAtlas(const Programs                  &prog,
                        const GlMat4                    &ortho,
                        const decimate::MeshData        &mesh,
                        const decimate::chart::ChartSet &cs,
                        const std::vector<math::Vec2f>  &uv) {
    const auto       &idx = mesh.indices;
    const std::size_t numTris = idx.size() / 3;
    const std::size_t numCornerUV = uv.size();

    // Filled triangles, interleaved (u, v, r, g, b); the lines below reuse the same layout so one program draws both.
    std::vector<float> fillBuf;
    fillBuf.reserve(numTris * 3 * 5);

    for (std::size_t t = 0; t < numTris; ++t) {
        // Fewer corner UVs than faces (stale atlas): skip rather than index past `uv`.
        if (3 * t + 2 >= numCornerUV) {
            continue;
        }
        // A face beyond the chart map (stale chartSet) falls back to chart 0's colour.
        int chartId = 0;
        if (t < cs.faceToChart.size()) {
            chartId = cs.faceToChart[t];
        }
        const Rgb col = chartColor(chartId);

        for (int k = 0; k < 3; ++k) {
            const auto &uvc = uv[3 * t + static_cast<std::size_t>(k)];
            fillBuf.push_back(uvc.x);
            fillBuf.push_back(uvc.y);
            fillBuf.push_back(col.r);
            fillBuf.push_back(col.g);
            fillBuf.push_back(col.b);
        }
    }

    // Triangle edges as black lines.
    std::vector<float> lineBuf;
    lineBuf.reserve(numTris * 6 * 5);

    for (std::size_t t = 0; t < numTris; ++t) {
        // Fewer corner UVs than faces (stale atlas): skip rather than index past `uv`.
        if (3 * t + 2 >= numCornerUV) {
            continue;
        }
        const auto &a = uv[3 * t + 0];
        const auto &b = uv[3 * t + 1];
        const auto &c = uv[3 * t + 2];
        // 3 edges x 2 endpoints; each endpoint = pos2 + black.
        auto push2 = [&](const math::Vec2f &p0, const math::Vec2f &p1) {
            for (const auto &p : { p0, p1 }) {
                lineBuf.push_back(p.x);
                lineBuf.push_back(p.y);
                lineBuf.push_back(0.0f);
                lineBuf.push_back(0.0f);
                lineBuf.push_back(0.0f);
            }
        };
        push2(a, b);
        push2(b, c);
        push2(c, a);
    }

    // White unit-square outline, 4 segments, interleaved (u, v, r, g, b).
    const float bW[] = {
        0.0f,
        0.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        0.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        0.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        0.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        0.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        0.0f,
        0.0f,
        1.0f,
        1.0f,
        1.0f,
    }; // NOLINT(modernize-avoid-c-arrays)

    glUseProgram(prog.ortho);
    glUniformMatrix4fv(glGetUniformLocation(prog.ortho, "uOrtho"), 1, GL_FALSE, ortho.m);

    // Fills first so the wireframe and border land on top (depth test is off in this panel).
    if (!fillBuf.empty()) {
        unsigned int vao = 0, vbo = 0;
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     static_cast<GLsizeiptr>(fillBuf.size() * sizeof(float)),
                     fillBuf.data(), GL_STREAM_DRAW);
        const int sB = 5 * static_cast<int>(sizeof(float));
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sB, nullptr);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sB,
                              // NOLINTNEXTLINE(performance-no-int-to-ptr)
                              reinterpret_cast<const void *>(2 * sizeof(float)));
        glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(fillBuf.size()) / 5);
        glBindVertexArray(0);
        glDeleteBuffers(1, &vbo);
        glDeleteVertexArrays(1, &vao);
    }

    if (!lineBuf.empty()) {
        unsigned int vao = 0, vbo = 0;
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     static_cast<GLsizeiptr>(lineBuf.size() * sizeof(float)),
                     lineBuf.data(), GL_STREAM_DRAW);
        const int sB = 5 * static_cast<int>(sizeof(float));
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sB, nullptr);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sB,
                              // NOLINTNEXTLINE(performance-no-int-to-ptr)
                              reinterpret_cast<const void *>(2 * sizeof(float)));
        glDrawArrays(GL_LINES, 0, static_cast<int>(lineBuf.size()) / 5);
        glBindVertexArray(0);
        glDeleteBuffers(1, &vbo);
        glDeleteVertexArrays(1, &vao);
    }

    {
        unsigned int vao = 0, vbo = 0;
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(bW), bW, GL_STREAM_DRAW);
        const int sB = 5 * static_cast<int>(sizeof(float));
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sB, nullptr);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sB,
                              // NOLINTNEXTLINE(performance-no-int-to-ptr)
                              reinterpret_cast<const void *>(2 * sizeof(float)));
        glDrawArrays(GL_LINES, 0, 8);
        glBindVertexArray(0);
        glDeleteBuffers(1, &vbo);
        glDeleteVertexArrays(1, &vao);
    }
}


/// Draw the baked normal map as a textured quad over [0,1]² (panel ③ 'n' mode).
static void drawNormalMapQuad(const Programs &prog, const GlMat4 &ortho, unsigned int tex) {
    // pos2 + uv2 interleaved, 4 vertices as a triangle strip.
    static constexpr float kQuad[] = {
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        1.0f,
        0.0f,
        1.0f,
        0.0f,
        0.0f,
        1.0f,
        0.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
    }; // NOLINT(modernize-avoid-c-arrays)

    unsigned int vao = 0, vbo = 0;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(kQuad), kQuad, GL_STREAM_DRAW);
    const int sB = 4 * static_cast<int>(sizeof(float));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sB, nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sB,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void *>(2 * sizeof(float)));

    glUseProgram(prog.quad);
    glUniformMatrix4fv(glGetUniformLocation(prog.quad, "uOrtho"), 1, GL_FALSE, ortho.m);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
    glUniform1i(glGetUniformLocation(prog.quad, "uTex"), 0);

    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    glBindTexture(GL_TEXTURE_2D, 0);
    glBindVertexArray(0);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);

    // The unit-square outline as 4 line segments (pos2 only).
    static constexpr float kBorder[] = {
        0.0f,
        0.0f,
        1.0f,
        0.0f,
        1.0f,
        0.0f,
        1.0f,
        1.0f,
        1.0f,
        1.0f,
        0.0f,
        1.0f,
        0.0f,
        1.0f,
        0.0f,
        0.0f,
    }; // NOLINT(modernize-avoid-c-arrays)
    // The ortho program wants pos2 + col3, so widen the outline with white.
    float borderFull[8 * 5]; // NOLINT(modernize-avoid-c-arrays)
    for (int i = 0; i < 8; ++i) {
        borderFull[i * 5 + 0] = kBorder[i * 2 + 0];
        borderFull[i * 5 + 1] = kBorder[i * 2 + 1];
        borderFull[i * 5 + 2] = 1.0f;
        borderFull[i * 5 + 3] = 1.0f;
        borderFull[i * 5 + 4] = 1.0f;
    }

    unsigned int bvao = 0, bvbo = 0;
    glGenVertexArrays(1, &bvao);
    glBindVertexArray(bvao);
    glGenBuffers(1, &bvbo);
    glBindBuffer(GL_ARRAY_BUFFER, bvbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(borderFull), borderFull, GL_STREAM_DRAW);
    const int bsB = 5 * static_cast<int>(sizeof(float));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, bsB, nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, bsB,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void *>(2 * sizeof(float)));
    glUseProgram(prog.ortho);
    glUniformMatrix4fv(glGetUniformLocation(prog.ortho, "uOrtho"), 1, GL_FALSE, ortho.m);
    glDrawArrays(GL_LINES, 0, 8);
    glBindVertexArray(0);
    glDeleteBuffers(1, &bvbo);
    glDeleteVertexArrays(1, &bvao);
}


/// Smooth per-vertex normals: the normalised sum of the UNIT flat normals of a vertex's faces
/// (one equal vote per face; not angle- or area-weighted). Supplies the interpolated-N channel
/// for the TBN decode. The bake takes its low-poly normals from bake::computeVertexNormals,
/// so this is only a stand-in that keeps the viewer independent of that module.
static std::vector<Vec3f> computeSmoothNormals(const decimate::MeshData &mesh) {
    const auto       &pos = mesh.positions;
    const auto       &idx = mesh.indices;
    const std::size_t nv = pos.size();
    const std::size_t numTris = idx.size() / 3;

    std::vector<Vec3f> acc(nv, { 0.0f, 0.0f, 0.0f });

    // Accumulate each face's unit normal into its three corner vertices.
    for (std::size_t t = 0; t < numTris; ++t) {
        const int ia = idx[t * 3 + 0];
        const int ib = idx[t * 3 + 1];
        const int ic = idx[t * 3 + 2];
        if (ia < 0 || ib < 0 || ic < 0 ||
            static_cast<std::size_t>(ia) >= nv ||
            static_cast<std::size_t>(ib) >= nv ||
            static_cast<std::size_t>(ic) >= nv) {
            continue;
        }
        const auto fn = flatNormal(pos[static_cast<std::size_t>(ia)],
                                   pos[static_cast<std::size_t>(ib)],
                                   pos[static_cast<std::size_t>(ic)]);
        for (const int vi : { ia, ib, ic }) {
            auto &a = acc[static_cast<std::size_t>(vi)];
            a.x += fn.x;
            a.y += fn.y;
            a.z += fn.z;
        }
    }

    // An isolated vertex (zero sum) gets +Y so the shader never normalises a zero vector.
    std::vector<Vec3f> result(nv);
    for (std::size_t i = 0; i < nv; ++i) {
        const auto &a = acc[i];
        const float len = std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
        result[i] = (len > 1e-8f) ? Vec3f { a.x / len, a.y / len, a.z / len }
                                  : Vec3f { 0.0f, 1.0f, 0.0f };
    }
    return result;
}


/// Draw panel ④: chartedMesh textured with the baked normal map, with TBN decode.
///
/// VBO layout: un-indexed (3 distinct vertices per face), interleaved
///   pos(3) + nrm(3) + uv(2) + tan(3) = 11 floats per corner.
///
/// aTan carries the per-FACE tangent (exact mapFaceRt formula, constant over the
/// triangle), replicated to all 3 corners.  Per-vertex averaging is intentionally
/// NOT used: at a chart seam a shared vertex belongs to multiple charts with
/// different per-corner UVs, so averaging cancels/garbles the tangent → noisy TBN.
/// Carrying the per-face constant means T is C0 (not C1) across triangle boundaries,
/// but the shader's Gram-Schmidt re-orthogonalises against the interpolated N anyway,
/// so the decode is smooth everywhere the bake is smooth.
///
/// aNrm carries the smooth (accumulated) vertex normal so the shader interpolates N
/// continuously — matching the bake's lowNormals interpolation per texel.
// NOLINTNEXTLINE(readability-function-size): builds the un-indexed VBO with 11-component interleaved data
static void drawBakedMesh(const Programs &prog,
                          const GlMat4 &mvp, const GlMat4 &mv,
                          const ViewerState              &s,
                          const decimate::MeshData       &mesh,
                          const std::vector<math::Vec2f> &cornerUV) {
    if (!s.baked || s.bakedTex == 0) {
        return;
    }

    const auto       &pos = mesh.positions;
    const auto       &idx = mesh.indices;
    const std::size_t numTris = idx.size() / 3;
    const auto        nv = static_cast<int>(pos.size());

    // Smooth vertex normals — matches the bake's lowNormals per-texel interpolation.
    const std::vector<Vec3f> smoothN = computeSmoothNormals(mesh);

    // Build un-indexed VBO: pos(3) + nrm(3) + uv(2) + tan(3) = 11 floats
    // Per-face tangent (constant over the triangle, replicated to all 3 corners).
    std::vector<float> buf;
    buf.reserve(numTris * 3 * 11);

    for (std::size_t t = 0; t < numTris; ++t) {
        const int ia = idx[t * 3 + 0];
        const int ib = idx[t * 3 + 1];
        const int ic = idx[t * 3 + 2];
        // A corrupt index would read outside `positions`: drop the triangle instead.
        if (ia >= nv || ib >= nv || ic >= nv || ia < 0 || ib < 0 || ic < 0) {
            continue;
        }
        // Fewer corner UVs than faces (stale atlas): skip rather than index past `cornerUV`.
        if (3 * t + 2 >= cornerUV.size()) {
            continue;
        }

        const auto &pa = pos[static_cast<std::size_t>(ia)];
        const auto &pb = pos[static_cast<std::size_t>(ib)];
        const auto &pc = pos[static_cast<std::size_t>(ic)];

        // Per-face tangent: exact mapFaceRt formula using THIS face's per-corner UVs.
        // T = (pb-pa)*(dv2/det) - (pc-pa)*(dv1/det)
        // Degenerate UV (|det| < 1e-12 or |tFace| ≈ 0): fall back to an
        // N-orthogonal vector so the shader doesn't divide by zero.
        const auto &uA = cornerUV[3 * t + 0];
        const auto &uB = cornerUV[3 * t + 1];
        const auto &uC = cornerUV[3 * t + 2];
        const float du1 = uB.x - uA.x, dv1 = uB.y - uA.y;
        const float du2 = uC.x - uA.x, dv2 = uC.y - uA.y;
        const float det = du1 * dv2 - du2 * dv1;

        Vec3f tFace { 1.0f, 0.0f, 0.0f }; // fallback
        if (std::fabs(det) > 1e-12f) {
            const float r = 1.0f / det;
            const float tx = (pb.x - pa.x) * (dv2 * r) - (pc.x - pa.x) * (dv1 * r);
            const float ty = (pb.y - pa.y) * (dv2 * r) - (pc.y - pa.y) * (dv1 * r);
            const float tz = (pb.z - pa.z) * (dv2 * r) - (pc.z - pa.z) * (dv1 * r);
            const float tlen = std::sqrt(tx * tx + ty * ty + tz * tz);
            if (tlen > 1e-8f) {
                tFace = { tx / tlen, ty / tlen, tz / tlen };
            }
        }

        // Un-indexed corners: position and smooth normal per vertex, UV per face corner (seam vertices
        // carry different UVs per chart), and the face-constant tangent.
        struct Corner {
            int vi;
            int k;
        };
        const Corner corners[3] = { { ia, 0 }, { ib, 1 }, { ic, 2 } }; // NOLINT(modernize-avoid-c-arrays)
        for (const auto &[vi, k] : corners) {
            const auto &p = pos[static_cast<std::size_t>(vi)];
            const auto &n = smoothN[static_cast<std::size_t>(vi)];
            const auto &uvc = cornerUV[3 * t + static_cast<std::size_t>(k)];
            buf.push_back(p.x);
            buf.push_back(p.y);
            buf.push_back(p.z);
            buf.push_back(n.x);
            buf.push_back(n.y);
            buf.push_back(n.z);
            buf.push_back(uvc.x);
            buf.push_back(uvc.y);
            buf.push_back(tFace.x);
            buf.push_back(tFace.y);
            buf.push_back(tFace.z);
        }
    }

    if (buf.empty()) {
        return;
    }

    const int stride = 11;
    const int vertCount = static_cast<int>(buf.size()) / stride;
    const int sBytes = stride * static_cast<int>(sizeof(float));

    unsigned int vao = 0, vbo = 0;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(buf.size() * sizeof(float)),
                 buf.data(), GL_STREAM_DRAW);

    // aPos @ 0 (3 floats, offset 0)
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sBytes, nullptr);
    // aNrm @ 1 (3 floats, offset 3)
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void *>(3 * sizeof(float)));
    // aUV @ 2 (2 floats, offset 6)
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void *>(6 * sizeof(float)));
    // aTan @ 3 (3 floats, offset 8)
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void *>(8 * sizeof(float)));

    glUseProgram(prog.baked);
    glUniformMatrix4fv(glGetUniformLocation(prog.baked, "uMVP"), 1, GL_FALSE, mvp.m);
    glUniformMatrix4fv(glGetUniformLocation(prog.baked, "uMV"), 1, GL_FALSE, mv.m);
    glUniform1i(glGetUniformLocation(prog.baked, "uTangentSpace"),
                s.bakeTangentSpace ? GL_TRUE : GL_FALSE);
    glUniform1i(glGetUniformLocation(prog.baked, "uLitGray"),
                s.panel4LitGray ? GL_TRUE : GL_FALSE);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, s.bakedTex);
    glUniform1i(glGetUniformLocation(prog.baked, "uTex"), 0);

    glDrawArrays(GL_TRIANGLES, 0, vertCount);

    glBindTexture(GL_TEXTURE_2D, 0);
    glBindVertexArray(0);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}


/// Panel ④ fallback when not yet baked: draw abs(face_normal) as RGB colour.
static void drawFaceNormalColors(const Programs &prog,
                                 const GlMat4 &mvp, const GlMat4 &mv,
                                 const decimate::MeshData &mesh) {
    const auto       &pos = mesh.positions;
    const auto       &idx = mesh.indices;
    const std::size_t numTris = idx.size() / 3;
    const auto        nv = static_cast<int>(pos.size());

    std::vector<float> buf;
    buf.reserve(numTris * 3 * 9);

    for (std::size_t t = 0; t < numTris; ++t) {
        const int ia = idx[t * 3 + 0];
        const int ib = idx[t * 3 + 1];
        const int ic = idx[t * 3 + 2];
        // A corrupt index would read outside `positions`: drop the triangle instead.
        if (ia >= nv || ib >= nv || ic >= nv || ia < 0 || ib < 0 || ic < 0) {
            continue;
        }
        const auto &a = pos[static_cast<std::size_t>(ia)];
        const auto &b = pos[static_cast<std::size_t>(ib)];
        const auto &c = pos[static_cast<std::size_t>(ic)];
        const auto  fn = flatNormal(a, b, c);
        const float cr = std::abs(fn.x), cg = std::abs(fn.y), cb2 = std::abs(fn.z);

        for (const auto &v : { a, b, c }) {
            buf.push_back(v.x);
            buf.push_back(v.y);
            buf.push_back(v.z);
            buf.push_back(fn.x);
            buf.push_back(fn.y);
            buf.push_back(fn.z);
            buf.push_back(cr);
            buf.push_back(cg);
            buf.push_back(cb2);
        }
    }

    if (buf.empty()) {
        return;
    }

    const int stride = 9;
    const int vertCount = static_cast<int>(buf.size()) / stride;
    const int sBytes = stride * static_cast<int>(sizeof(float));

    // Throwaway VAO + VBO per draw (see drawMeshFlat).
    unsigned int vao = 0, vbo = 0;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(buf.size() * sizeof(float)),
                 buf.data(), GL_STREAM_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sBytes, nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void *>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void *>(6 * sizeof(float)));

    glUseProgram(prog.mesh);
    glUniformMatrix4fv(glGetUniformLocation(prog.mesh, "uMVP"), 1, GL_FALSE, mvp.m);
    glUniformMatrix4fv(glGetUniformLocation(prog.mesh, "uMV"), 1, GL_FALSE, mv.m);

    glDrawArrays(GL_TRIANGLES, 0, vertCount);

    glBindVertexArray(0);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}


/// Perspective + orbit matrices for a `vpW`x`vpH` viewport; the clip planes scale with the orbit distance.
[[nodiscard]] static CameraMatrices makeCamera(const ViewerState &s, int vpW, int vpH) noexcept {
    // A zero-height viewport (minimised window) must not divide by zero.
    const float aspect = (vpH > 0)
                             ? (static_cast<float>(vpW) / static_cast<float>(vpH))
                             : 1.0f;
    const float nearP = s.distance * 0.01f;
    const float farPlane = s.distance * 100.0f;

    CameraMatrices c;
    c.proj = makePerspective(45.0f, aspect, nearP, farPlane);
    c.mv = makeOrbitMV(s.azimuth, s.elevation, s.distance,
                       s.centreX, s.centreY, s.centreZ);
    c.mvp = mat4Mul(c.proj, c.mv);
    return c;
}


/// Render one 3D panel (panels ① and ②).
static void render3DPanel(const Programs &prog, const ViewerState &s,
                          int vpX, int vpY, int vpW, int vpH,
                          bool useChartColors) {
    glViewport(vpX, vpY, vpW, vpH);
    glEnable(GL_DEPTH_TEST);
    glClear(GL_DEPTH_BUFFER_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, s.wireframe ? GL_LINE : GL_FILL);

    const auto cam = makeCamera(s, vpW, vpH);

    if (!useChartColors) {
        drawMeshFlat(prog, cam.mvp, cam.mv, s.current, 0.7f, 0.75f, 0.82f);
    } else if (s.charted) {
        drawMeshChartColors(prog, cam.mvp, cam.mv, s.atlas.mesh, s.atlas.charts);
    } else {
        drawMeshFlat(prog, cam.mvp, cam.mv, s.current, 0.55f, 0.55f, 0.58f);
    }
}


/// Render the UV atlas panel (panel ③).
static void renderUVPanel(const Programs &prog, const ViewerState &s,
                          int vpX, int vpY, int vpW, int vpH) {
    glViewport(vpX, vpY, vpW, vpH);
    glClear(GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);

    constexpr float kMargin = 0.05f;
    const GlMat4    ortho = makeOrtho(-kMargin, 1.0f + kMargin, -kMargin, 1.0f + kMargin);

    if (s.atlasShowNormalMap && s.baked && s.bakedTex != 0) {
        drawNormalMapQuad(prog, ortho, s.bakedTex);
    } else if (s.charted) {
        drawUVAtlas(prog, ortho, s.atlas.mesh, s.atlas.charts, s.atlas.uv);
    }
    // Not charted: dark background only.

    glEnable(GL_DEPTH_TEST);
}


/// Render panel ④ — baked normal map on the low-poly mesh with TBN decode.
static void renderBakedPanel(const Programs &prog, const ViewerState &s,
                             int vpX, int vpY, int vpW, int vpH) {
    glViewport(vpX, vpY, vpW, vpH);
    glEnable(GL_DEPTH_TEST);
    glClear(GL_DEPTH_BUFFER_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    const auto cam = makeCamera(s, vpW, vpH);

    if (s.baked && s.bakedTex != 0) {
        const bool  lod = s.lodLevel > 0 && s.lodValid;
        const auto &view = s.lod[static_cast<std::size_t>(std::max(s.lodLevel, 1) - 1)];
        if (lod) {
            drawBakedMesh(prog, cam.mvp, cam.mv, s, view.mesh, view.uv);
        } else {
            drawBakedMesh(prog, cam.mvp, cam.mv, s, s.atlas.mesh, s.atlas.uv);
        }
    } else {
        drawFaceNormalColors(prog, cam.mvp, cam.mv, s.current);
    }
}


/// Rewrites the window title as a HUD: counts, atlas mode, active LOD and the key legend.
static void updateTitle(GLFWwindow *win, const ViewerState &s) {
    char      buf[384]; // NOLINT(modernize-avoid-c-arrays)
    const int vertexCount = static_cast<int>(s.current.positions.size());
    const int triangleCount = static_cast<int>(s.current.indices.size() / 3);

    char chartFrag[64]; // NOLINT(modernize-avoid-c-arrays)
    if (s.charted) {
        const int numCharts = static_cast<int>(s.atlas.charts.charts.size());
        std::snprintf(chartFrag, sizeof(chartFrag), "charts:%d", numCharts);
    } else {
        std::snprintf(chartFrag, sizeof(chartFrag), "(press c to chartify)");
    }

    char lodFrag[64] = ""; // NOLINT(modernize-avoid-c-arrays)
    if (s.lodLevel > 0 && s.lodValid) {
        std::snprintf(lodFrag, sizeof(lodFrag), "  [LOD target %d%%: %d verts]",
                      kLodPercent[s.lodLevel - 1], s.lod[static_cast<std::size_t>(s.lodLevel - 1)].vertices);
    }
    std::snprintf(buf, sizeof(buf),
                  "mesh-viewer  |  verts:%d  tris:%d  %s  [atlas:%s]%s"
                  "  | -/+ reduce  r reset  c chartify  t tangent  l LOD (erased map; e locked)  n atlas/normalmap  g shade  w wire  q quit",
                  vertexCount, triangleCount, chartFrag,
                  s.atlasShowNormalMap ? "normalmap" : "charts", lodFrag);
    glfwSetWindowTitle(win, buf);
}


/// Re-decimates `s.original` down to `s.currentTarget` vertices with the active reducer and logs the result.
static void applyReduce(ViewerState &s) {
    // At or above the original size there is nothing to reduce: show the untouched mesh.
    if (s.currentTarget >= s.originalVertexCount) {
        s.current = s.original;
        return;
    }
    switch (s.reducer) {
        case 1:
            s.current = decimate::reduceClustering(s.original, s.currentTarget);
            break;
        case 2:
            s.current = decimate::reduceHalfEdge(s.original, s.currentTarget);
            break;
        default:
            s.current = decimate::reduce(s.original, s.currentTarget);
            break;
    }

    static const char *const kNames[] = { "QEM", "CLUSTER", "HALFEDGE" };
    std::fprintf(stderr, "[reduce] %s  target=%d  actual verts=%d  tris=%d\n",
                 kNames[s.reducer],
                 s.currentTarget,
                 static_cast<int>(s.current.positions.size()),
                 static_cast<int>(s.current.indices.size() / 3));
}


/// Key '-': shrink to 70% of the displayed vertices, never below 4.
static void decimateStep(ViewerState &s) {
    const int displayed = static_cast<int>(s.current.positions.size());
    s.currentTarget = std::max(4, displayed * 7 / 10);
    applyReduce(s);
}


/// Key '+': grow by 40% (+1 so tiny meshes still advance), capped at the original size.
static void growStep(ViewerState &s) {
    const int displayed = static_cast<int>(s.current.positions.size());
    s.currentTarget = std::min(s.originalVertexCount, (displayed * 14 / 10) + 1);
    applyReduce(s);
}


/// Key 'X': flip seam straightening and, if already chartified, redo chartify + bake.
static void toggleStraighten(ViewerState &s) {
    s.straighten = !s.straighten;
    std::fprintf(stderr, "[straighten] %s\n", s.straighten ? "on" : "off");
    if (!s.charted) {
        std::fprintf(stderr, "[straighten] press 'c' to chartify and bake\n");
        return;
    }
    computeChartification(s);
}


/// GLFW key callback: dispatches the viewer keys (legend in main.cpp).
static void onKey(GLFWwindow *win, int key, int /*scancode*/, int action, int /*mods*/) {
    if (action != GLFW_PRESS && action != GLFW_REPEAT) {
        return;
    }
    ViewerState &s = *gState;

    if (key == GLFW_KEY_ESCAPE || key == GLFW_KEY_Q) {
        glfwSetWindowShouldClose(win, GLFW_TRUE);
        return;
    }

    // Reset drops the chartification too: it described the mesh that is being thrown away.
    if (key == GLFW_KEY_R) {
        s.currentTarget = s.originalVertexCount;
        s.current = s.original;
        s.charted = false;
        s.baked = false;
        s.lodLevel = 0;
        std::fprintf(stderr, "[reset] verts=%d  tris=%d\n",
                     static_cast<int>(s.current.positions.size()),
                     static_cast<int>(s.current.indices.size() / 3));
        return;
    }
    if (key == GLFW_KEY_W) {
        s.wireframe = !s.wireframe;
        return;
    }

    if (key == GLFW_KEY_K) {
        s.reducer = (s.reducer + 1) % 3; // QEM -> CLUSTER -> HALFEDGE -> QEM
        static const char *const kNames[] = { "QEM (heap edge-collapse)",
                                              "CLUSTER (octree, fast/sloppy)",
                                              "HALFEDGE (corner-table edge-collapse)" };
        std::fprintf(stderr, "[reducer] %s\n", kNames[s.reducer]);
        applyReduce(s); // re-decimate at the same target with the new strategy
        return;
    }

    if (key == GLFW_KEY_S) {
        if (s.baked) {
            const bool ok = ui::writeTGA("normalmap.tga", s.bakedMap);
            if (ok) { // NOLINT(bugprone-branch-clone)
                std::fprintf(stderr, "[save] normalmap.tga written (%dx%d)\n",
                             s.bakedMap.width, s.bakedMap.height);
            } else {
                std::fprintf(stderr, "[save] ERROR: could not write normalmap.tga\n");
            }
        } else {
            std::fprintf(stderr, "[save] press 'c' first to bake the normal map\n");
        }
        return;
    }

    // 'c' / 't' choose the normal space and (re)chartify + bake in it; the Seamless gutter picks the one
    // that hides its seams in that space.
    if (key == GLFW_KEY_C) {
        s.bakeTangentSpace = false;
        computeChartification(s);
        return;
    }
    if (key == GLFW_KEY_T) {
        s.bakeTangentSpace = true;
        computeChartification(s);
        return;
    }

    if (key == GLFW_KEY_MINUS || key == GLFW_KEY_DOWN) {
        decimateStep(s);
        return;
    }
    if (key == GLFW_KEY_EQUAL || key == GLFW_KEY_UP) {
        growStep(s);
        return;
    }

    if (key == GLFW_KEY_E) {
        if (!s.charted) {
            std::fprintf(stderr, "[gutter] press 'c' first to chartify and bake\n");
            return;
        }
        if (s.lodLevel > 0 && s.lodValid) { // the LODs are drawn with the erased map only
            std::fprintf(stderr, "[gutter] locked to 'erased' in the LOD view: press 'l' back to full\n");
            return;
        }
        s.gutterMode = nextGutterMode(s.gutterMode);
        std::fprintf(stderr, "[gutter] %s\n", gutterName(s.gutterMode));
        bakeAndUpload(s);
        return;
    }

    if (key == GLFW_KEY_X) {
        toggleStraighten(s);
        return;
    }
    if (key == GLFW_KEY_L) {
        cycleLod(s);
        return;
    }
    if (key == GLFW_KEY_N) {
        s.atlasShowNormalMap = !s.atlasShowNormalMap;
        return;
    }
    if (key == GLFW_KEY_G) {
        s.panel4LitGray = !s.panel4LitGray;
        std::fprintf(stderr, "[panel4] %s\n",
                     s.panel4LitGray ? "lit grayscale" : "abs-normal RGB");
        return;
    }
}


/// Left button starts/stops the orbit drag.
static void onMouseButton(GLFWwindow * /*win*/, int button, int action, int /*mods*/) {
    if (button != GLFW_MOUSE_BUTTON_LEFT) {
        return;
    }
    ViewerState &s = *gState;
    s.leftDragging = (action == GLFW_PRESS);
}


/// Orbit drag: 0.5 degree per pixel; elevation is clamped short of the poles so the camera never flips.
static void onCursorPos(GLFWwindow * /*win*/, double xpos, double ypos) {
    ViewerState &s = *gState;
    if (!s.leftDragging) {
        s.lastMouseX = xpos;
        s.lastMouseY = ypos;
        return;
    }

    const double dx = xpos - s.lastMouseX;
    const double dy = ypos - s.lastMouseY;
    s.azimuth += static_cast<float>(dx) * 0.5f;
    s.elevation += static_cast<float>(dy) * 0.5f;
    s.elevation = std::max(-89.0f, std::min(89.0f, s.elevation));
    s.lastMouseX = xpos;
    s.lastMouseY = ypos;
}


/// Wheel zoom: 10% per notch, never closer than 0.01.
static void onScroll(GLFWwindow * /*win*/, double /*xoff*/, double yoff) {
    ViewerState &s = *gState;
    s.distance *= (yoff > 0) ? 0.9f : 1.1f;
    s.distance = std::max(0.01f, s.distance);
}


/// Smoke-test line: did the bake produce a NormalMap with exactly width*height texels?
static void logBakeStatus(const ViewerState &s) {
    const int  nmW = s.bakedMap.width;
    const int  nmH = s.bakedMap.height;
    const auto nmPixels = static_cast<int>(s.bakedMap.pixels.size());
    const bool ok = (nmW > 0 && nmH > 0 && nmPixels == nmW * nmH);
    if (ok) { // NOLINT(bugprone-branch-clone)
        std::fprintf(stderr,
                     "[smoke] panel④ bake OK: %dx%d NormalMap (%d texels) baked=%s\n",
                     nmW, nmH, nmPixels, s.baked ? "true" : "false");
    } else {
        std::fprintf(stderr,
                     "[smoke] panel④ bake FAIL: w=%d h=%d pixels=%d baked=%s\n",
                     nmW, nmH, nmPixels, s.baked ? "true" : "false");
    }
}


/// Creates an offscreen RGBA8 + depth24 framebuffer (completeness is not checked).
[[nodiscard]] static ShotTarget makeShotTarget(int w, int h) {
    ShotTarget t;
    t.w = w;
    t.h = h;
    glGenRenderbuffers(1, &t.color);
    glBindRenderbuffer(GL_RENDERBUFFER, t.color);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, w, h);
    glGenRenderbuffers(1, &t.depth);
    glBindRenderbuffer(GL_RENDERBUFFER, t.depth);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
    glGenFramebuffers(1, &t.fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, t.fbo);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, t.color);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, t.depth);
    return t;
}


/// Releases the framebuffer and both renderbuffers of `t`.
static void destroyShotTarget(ShotTarget &t) {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glDeleteFramebuffers(1, &t.fbo);
    glDeleteRenderbuffers(1, &t.color);
    glDeleteRenderbuffers(1, &t.depth);
}


/// Clear the target to `bg`, run `draw` into the first `w`x`h` pixels and read them back
/// as a top-row-first Image (glReadPixels rows are bottom-up, so they are flipped here).
template <class Draw>
static ui::Image renderImage(int w, int h, float bg, Draw draw) {
    glDisable(GL_SCISSOR_TEST);
    glViewport(0, 0, w, h);
    glClearColor(bg, bg, bg, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    draw();

    ui::Image img;
    img.w = w;
    img.h = h;
    img.rgb.resize(static_cast<std::size_t>(w) * static_cast<std::size_t>(h) * 3);
    std::vector<unsigned char> raw(img.rgb.size());
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, raw.data());
    const std::size_t row = static_cast<std::size_t>(w) * 3;
    for (int y = 0; y < h; ++y) {
        std::copy_n(raw.begin() + static_cast<std::ptrdiff_t>(static_cast<std::size_t>(h - 1 - y) * row), row,
                    img.rgb.begin() + static_cast<std::ptrdiff_t>(static_cast<std::size_t>(y) * row));
    }
    return img;
}


/// Writes `img` as a TGA and logs the outcome.
static bool writeShot(const std::string &path, const ui::Image &img) {
    const bool ok = ui::writeTGA(path, img);
    std::fprintf(stderr, "[shot] %s %s (%dx%d)\n", ok ? "wrote" : "FAILED", path.c_str(), img.w, img.h);
    return ok;
}


/// Renders with `draw` and writes the result as a TGA.
template <class Draw>
static bool shotPanel(const std::string &path, int w, int h, float bg, Draw draw) {
    return writeShot(path, renderImage(w, h, bg, draw));
}


/// Frame the mesh for a README shot: iteratively re-centre the orbit pivot on the projected
/// extent and scale the distance until the model's projected bounds fill `fill` of the frame
/// (height- or width-limited, whichever is tighter). Uses the real vertices, not the AABB
/// sphere, which over-estimates by ~30-40% for compact models.
static void fitCameraToMesh(ViewerState &s, float aspect, float fill) {
    constexpr float kTanHalfFov = 0.41421356f; // tan(45 deg / 2), matches makeCamera
    for (int iter = 0; iter < 24; ++iter) {
        const GlMat4 mv = makeOrbitMV(s.azimuth, s.elevation, s.distance, s.centreX, s.centreY, s.centreZ);
        float        loX = 1e30f, hiX = -1e30f, loY = 1e30f, hiY = -1e30f;
        for (const auto &p : s.current.positions) {
            const float depth = -(mv.at(0, 2) * p.x + mv.at(1, 2) * p.y + mv.at(2, 2) * p.z + mv.at(3, 2));
            const float nx = (mv.at(0, 0) * p.x + mv.at(1, 0) * p.y + mv.at(2, 0) * p.z + mv.at(3, 0)) /
                             (depth * kTanHalfFov * aspect);
            const float ny = (mv.at(0, 1) * p.x + mv.at(1, 1) * p.y + mv.at(2, 1) * p.z + mv.at(3, 1)) /
                             (depth * kTanHalfFov);
            loX = std::min(loX, nx);
            hiX = std::max(hiX, nx);
            loY = std::min(loY, ny);
            hiY = std::max(hiY, ny);
        }
        // Pan the pivot along the camera's right/up axes (rows 0/1 of the view rotation)
        // by the offset of the extent's centre, in world units at the pivot depth.
        const float panX = 0.5f * (loX + hiX) * s.distance * kTanHalfFov * aspect;
        const float panY = 0.5f * (loY + hiY) * s.distance * kTanHalfFov;
        s.centreX += mv.at(0, 0) * panX + mv.at(0, 1) * panY;
        s.centreY += mv.at(1, 0) * panX + mv.at(1, 1) * panY;
        s.centreZ += mv.at(2, 0) * panX + mv.at(2, 1) * panY;
        const float half = std::max(0.5f * (hiX - loX), 0.5f * (hiY - loY));
        s.distance *= half / fill;
    }
}


/// Top-left of the `win`x`win` window (8-px grid) where the two bakes differ most:
/// that is where the gutter fill is visible.
[[nodiscard]] static std::pair<int, int> mostDifferentWindow(const decimate::bake::NormalMap &a,
                                                             const decimate::bake::NormalMap &b, int win) {
    std::pair<int, int> best { 0, 0 };
    long                bestScore = -1;
    for (int y0 = 0; y0 + win <= a.height; y0 += 8) {
        for (int x0 = 0; x0 + win <= a.width; x0 += 8) {
            long score = 0;
            for (int y = y0; y < y0 + win; ++y) {
                for (int x = x0; x < x0 + win; ++x) {
                    const auto i = static_cast<std::size_t>(y) * static_cast<std::size_t>(a.width) + static_cast<std::size_t>(x);
                    score += (a.pixels[i] != b.pixels[i]) ? 1 : 0;
                }
            }
            if (score > bestScore) {
                bestScore = score;
                best = { x0, y0 };
            }
        }
    }
    return best;
}


/// Nearest-neighbour `scale`x crop; rows flipped so v points up like the atlas panel.
[[nodiscard]] static decimate::bake::NormalMap cropZoom(const decimate::bake::NormalMap &src,
                                                        int x0, int y0, int win, int scale) {
    decimate::bake::NormalMap out;
    out.width = win * scale;
    out.height = win * scale;
    out.pixels.resize(static_cast<std::size_t>(out.width) * static_cast<std::size_t>(out.height));
    for (int y = 0; y < out.height; ++y) {
        const int sy = y0 + (win - 1 - y / scale);
        for (int x = 0; x < out.width; ++x) {
            const int sx = x0 + x / scale;
            out.pixels[static_cast<std::size_t>(y) * static_cast<std::size_t>(out.width) + static_cast<std::size_t>(x)] =
                src.pixels[static_cast<std::size_t>(sy) * static_cast<std::size_t>(src.width) + static_cast<std::size_t>(sx)];
        }
    }
    return out;
}


/// Rows of a NormalMap as an Image (row 0 = top), for the 2D composites.
[[nodiscard]] static ui::Image toImage(const decimate::bake::NormalMap &nm) {
    ui::Image im;
    im.w = nm.width;
    im.h = nm.height;
    im.rgb.reserve(nm.pixels.size() * 3);
    for (const auto &p : nm.pixels) {
        im.rgb.insert(im.rgb.end(), { p.r, p.g, p.b });
    }
    return im;
}


/// The baked-model panels (normal-map quad, abs-normal RGB, lit grayscale) for whatever
/// gutter mode `s.gutterMode` currently selects; `suffix` is "", "_dilate" or "_erased".
static LitShots shotBakePanels(const Programs &prog, ViewerState &s, const ShotTarget &t, const ShotTarget &big,
                               const std::string &base, const std::string &suffix, bool &ok) {
    const int sq = std::min(t.w, t.h);
    glBindFramebuffer(GL_FRAMEBUFFER, t.fbo);
    s.atlasShowNormalMap = true;
    ok &= shotPanel(base + "_normalmap" + suffix + ".tga", sq, sq, 0.82f,
                    [&] { renderUVPanel(prog, s, 0, 0, sq, sq); });
    s.atlasShowNormalMap = false;
    s.panel4LitGray = false;
    ok &= shotPanel(base + "_baked" + suffix + ".tga", t.w, t.h, 1.0f,
                    [&] { renderBakedPanel(prog, s, 0, 0, t.w, t.h); });

    s.panel4LitGray = true; // the G-key path: uLitGray in the baked fragment shader
    LitShots lit;
    lit.frame = renderImage(t.w, t.h, 1.0f, [&] { renderBakedPanel(prog, s, 0, 0, t.w, t.h); });
    ok &= writeShot(base + "_baked_lit" + suffix + ".tga", lit.frame);
    glBindFramebuffer(GL_FRAMEBUFFER, big.fbo);
    lit.big = renderImage(big.w, big.h, 1.0f, [&] { renderBakedPanel(prog, s, 0, 0, big.w, big.h); });
    glBindFramebuffer(GL_FRAMEBUFFER, t.fbo);
    s.panel4LitGray = false;
    return lit;
}


/// Draws a ShotBatch triangle soup with the ortho program (no depth test).
static void drawBatch(const Programs &prog, const GlMat4 &ortho, const ShotBatch &b) {
    if (b.v.empty()) {
        return;
    }
    glDisable(GL_DEPTH_TEST);
    glUseProgram(prog.ortho);
    glUniformMatrix4fv(glGetUniformLocation(prog.ortho, "uOrtho"), 1, GL_FALSE, ortho.m);
    unsigned int vao = 0, vbo = 0;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(b.v.size() * sizeof(float)), b.v.data(), GL_STREAM_DRAW);
    const int stride = 5 * static_cast<int>(sizeof(float));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, stride, nullptr);
    glEnableVertexAttribArray(1);
    // NOLINTNEXTLINE(performance-no-int-to-ptr)
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<const void *>(2 * sizeof(float)));
    glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(b.v.size() / 5));
    glBindVertexArray(0);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
    glEnable(GL_DEPTH_TEST);
}


/// The packer's output for the current chartification (`s.atlas.packed`, computed once by
/// computeChartification), kept so `gutterLoop` can be read.
/// Returns empty if its UVs disagree with the viewer's `s.atlas.uv` (straighten rewrites it; then the outlines would lie).
static std::vector<decimate::chart::HalfEdgeChart> packedCharts(const ViewerState &s) {
    std::vector<decimate::chart::HalfEdgeChart> charts = s.atlas.packed;
    float                                       worst = 0.0f;
    for (const auto &c : charts) {
        for (std::size_t i = 0; i < c.faces.size(); ++i) {
            const auto f = static_cast<std::size_t>(c.faces[i]);
            for (std::size_t k = 0; k < 3; ++k) {
                const Vec2f &a = c.uv[static_cast<std::size_t>(
                    decimate::chart::localCorner(c, static_cast<int>(i), static_cast<int>(k)))];
                const Vec2f &b = s.atlas.uv[3 * f + k];
                worst = std::max({ worst, std::fabs(a.x - b.x), std::fabs(a.y - b.y) });
            }
        }
    }
    std::fprintf(stderr, "[shot] packed %zu charts; max |packed uv - s.atlas.uv| = %.2e\n", charts.size(),
                 static_cast<double>(worst));
    if (worst > 1e-5f) {
        std::fprintf(stderr, "[shot] packed charts do not match s.atlas.uv; gutter images skipped\n");
        charts.clear();
    }
    return charts;
}


constexpr Rgb   kGutterInk { 0.03f, 0.03f, 0.08f };
constexpr Rgb   kWireInk { 0.30f, 0.30f, 0.30f };
constexpr Rgb   kBoundaryInk { 0.05f, 0.25f, 0.85f };
constexpr float kShotBg = 0.97f;


/// Add the filled triangles (and optionally edges) of faces whose chart passes `pick`.
template <class Pick>
static void addChartFaces(ShotBatch &fill, ShotBatch &wire, const ViewerState &s, float wireW, Pick pick) {
    const std::size_t nf = s.atlas.mesh.indices.size() / 3;
    for (std::size_t f = 0; f < nf && 3 * f + 2 < s.atlas.uv.size() && f < s.atlas.charts.faceToChart.size(); ++f) {
        const int ci = s.atlas.charts.faceToChart[f];
        if (!pick(ci)) {
            continue;
        }
        const Rgb col = chartColor(ci);
        fill.tri(s.atlas.uv[3 * f], s.atlas.uv[3 * f + 1], s.atlas.uv[3 * f + 2], col);
        wire.loop({ s.atlas.uv[3 * f], s.atlas.uv[3 * f + 1], s.atlas.uv[3 * f + 2] }, wireW, kWireInk);
    }
}


/// <model>_atlas_gutter: chart fills + triangle edges + every chart's gutterLoop, light background.
static bool shotAtlasGutter(const Programs &prog, const ViewerState &s,
                            const std::vector<decimate::chart::HalfEdgeChart> &charts, const ShotTarget &ss, int n,
                            const std::string &path) {
    constexpr float kMargin = 0.05f, kSS = 2.0f;
    const float     uvPerPx = (1.0f + 2.0f * kMargin) / (static_cast<float>(n) * kSS); // one supersampled px
    ShotBatch       fill, wire, ink;
    addChartFaces(fill, wire, s, 1.0f * kSS * uvPerPx, [](int) { return true; });
    ink.loop({ { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } }, 1.0f * kSS * uvPerPx, Rgb { 0.6f, 0.6f, 0.6f });
    for (const auto &c : charts) {
        ink.loop(c.gutterLoop, 1.5f * kSS * uvPerPx, kGutterInk);
    }
    glBindFramebuffer(GL_FRAMEBUFFER, ss.fbo);
    const int       big = n * static_cast<int>(kSS);
    const GlMat4    ortho = makeOrtho(-kMargin, 1.0f + kMargin, -kMargin, 1.0f + kMargin);
    const ui::Image img = renderImage(big, big, kShotBg, [&] {
        drawBatch(prog, ortho, fill);
        drawBatch(prog, ortho, wire);
        drawBatch(prog, ortho, ink);
    });
    return writeShot(path, ui::downscale(img, static_cast<int>(kSS)));
}


/// Smallest interior angle among CONVEX gutterLoop vertices of charts with >= `minFaces` faces.
/// Vertices that are miter-fold artefacts (the loop folds back on itself: edges shorter than the
/// gutter, offset well under the gutter or beyond the miter limit) are skipped: they are not corners.
static SharpCorner sharpestConvexCorner(const ViewerState &s, const std::vector<decimate::chart::HalfEdgeChart> &charts,
                                        const std::vector<int> &facesPerChart, int minFaces) {
    SharpCorner                     best;
    const float                     gutter = static_cast<float>(decimate::chart::kGutterPx) / 255.0f;
    std::vector<std::vector<Vec2f>> pts(charts.size());
    for (std::size_t f = 0; f < s.atlas.charts.faceToChart.size() && 3 * f + 2 < s.atlas.uv.size(); ++f) {
        const int ci = s.atlas.charts.faceToChart[f];
        if (ci >= 0 && static_cast<std::size_t>(ci) < charts.size()) {
            pts[static_cast<std::size_t>(ci)].insert(pts[static_cast<std::size_t>(ci)].end(), s.atlas.uv.begin() + 3 * static_cast<std::ptrdiff_t>(f),
                                                     s.atlas.uv.begin() + 3 * static_cast<std::ptrdiff_t>(f) + 3);
        }
    }
    for (std::size_t ci = 0; ci < charts.size(); ++ci) {
        const auto &g = charts[ci].gutterLoop;
        if (g.size() < 3 || facesPerChart[ci] < minFaces) {
            continue;
        }
        float area2 = 0.0f;
        for (std::size_t i = 0; i < g.size(); ++i) {
            const Vec2f &a = g[i];
            const Vec2f &b = g[(i + 1) % g.size()];
            area2 += a.x * b.y - b.x * a.y;
        }
        const float sign = area2 >= 0.0f ? 1.0f : -1.0f;
        for (std::size_t i = 0; i < g.size(); ++i) {
            const Vec2f &p = g[(i + g.size() - 1) % g.size()];
            const Vec2f &c = g[i];
            const Vec2f &n = g[(i + 1) % g.size()];
            const float  e1x = c.x - p.x, e1y = c.y - p.y, e2x = n.x - c.x, e2y = n.y - c.y;
            if (std::hypot(e1x, e1y) < gutter || std::hypot(e2x, e2y) < gutter) {
                continue; // edges shorter than the gutter: fold territory
            }
            const float cross = e1x * e2y - e1y * e2x;
            if (cross * sign <= 0.0f) {
                continue; // reflex or straight
            }
            const float turn = std::atan2(std::fabs(cross), e1x * e2x + e1y * e2y);
            const float ang = 180.0f - turn * 180.0f / 3.14159265f;
            if (ang >= best.angleDeg) {
                continue;
            }
            // Distance of a loop point to the chart's corner UVs (and which UV is closest).
            const auto nearestTo = [&](const Vec2f &p, Vec2f &at) {
                float d2 = 1e30f;
                for (const Vec2f &q : pts[ci]) {
                    const float dd = (q.x - p.x) * (q.x - p.x) + (q.y - p.y) * (q.y - p.y);
                    if (dd < d2) {
                        d2 = dd;
                        at = q;
                    }
                }
                return std::sqrt(d2);
            };
            // The whole stretch i-3..i+3 must be a clean offset: a fold anywhere nearby
            // (distance to the chart well under the gutter) makes a messy picture, not a miter.
            bool  clean = true;
            Vec2f near = c;
            float dist = 0.0f;
            for (int k = -3; k <= 3 && clean; ++k) {
                Vec2f       at = c;
                const float d = nearestTo(g[static_cast<std::size_t>((static_cast<long>(i) + k + static_cast<long>(g.size())) % static_cast<long>(g.size()))], at);
                clean = d >= 0.9f * gutter && d <= decimate::chart::kDefaultMiterLimit * gutter * 1.05f;
                if (k == 0) {
                    near = at;
                    dist = d;
                }
            }
            // ...and the loop must not cross itself there (a fold also crosses its own edges).
            const auto at = [&](int k) {
                return g[static_cast<std::size_t>((static_cast<long>(i) + k + static_cast<long>(g.size())) % static_cast<long>(g.size()))];
            };
            const auto side = [](const Vec2f &a, const Vec2f &b, const Vec2f &p) {
                return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
            };
            for (int a = -3; a < 3 && clean; ++a) {
                for (int b = a + 2; b < 3 && clean; ++b) {
                    const bool crosses = side(at(a), at(a + 1), at(b)) * side(at(a), at(a + 1), at(b + 1)) < 0.0f &&
                                         side(at(b), at(b + 1), at(a)) * side(at(b), at(b + 1), at(a + 1)) < 0.0f;
                    clean = !crosses;
                }
            }
            // A tip beyond the miter limit is clamped (a truncated spike, not a miter): skip it.
            if (!clean || dist > 0.97f * decimate::chart::kDefaultMiterLimit * gutter) {
                continue;
            }
            best = { static_cast<int>(ci), static_cast<int>(i), ang, near, dist };
        }
    }
    return best;
}


/// <model>_gutter_corner: vector close-up (rendered at high res, not a pixel crop) of the sharpest
/// convex gutterLoop corner: that chart's triangles, its boundary, and the gutter outline.
static bool shotGutterCorner(const Programs &prog, const ViewerState &s,
                             const std::vector<decimate::chart::HalfEdgeChart> &charts, const ShotTarget &ss, int n,
                             const std::string &path) {
    std::vector<int> facesPerChart(charts.size(), 0);
    for (const int ci : s.atlas.charts.faceToChart) {
        if (ci >= 0 && static_cast<std::size_t>(ci) < charts.size()) {
            ++facesPerChart[static_cast<std::size_t>(ci)];
        }
    }
    const int   maxFaces = *std::max_element(facesPerChart.begin(), facesPerChart.end());
    const int   minFaces = std::max(30, maxFaces / 4);
    SharpCorner pick = sharpestConvexCorner(s, charts, facesPerChart, minFaces);
    if (pick.chart < 0) {
        pick = sharpestConvexCorner(s, charts, facesPerChart, 0);
    }
    if (pick.chart < 0) {
        std::fprintf(stderr, "[shot] no convex gutterLoop corner found; no corner image\n");
        return true;
    }
    const auto &gl = charts[static_cast<std::size_t>(pick.chart)].gutterLoop;
    const Vec2f g = gl[static_cast<std::size_t>(pick.vertex)];

    // Centre between the gutter vertex and its chart boundary vertex so both are in frame.
    const Vec2f ctr((g.x + pick.nearest.x) * 0.5f, (g.y + pick.nearest.y) * 0.5f);
    const float half = std::max(3.0f * pick.dist, 12.0f / 255.0f);
    std::fprintf(stderr, "[shot] corner: chart %d (%d faces) gutter vertex %d, interior angle %.1f deg, view +-%.4f uv\n",
                 pick.chart, facesPerChart[static_cast<std::size_t>(pick.chart)], pick.vertex,
                 static_cast<double>(pick.angleDeg), static_cast<double>(half));

    constexpr float kSS = 2.0f;
    const float     uvPerPx = 2.0f * half / (static_cast<float>(n) * kSS);
    ShotBatch       fill, wire, ink;
    addChartFaces(fill, wire, s, 1.5f * kSS * uvPerPx, [&](int ci) { return ci == pick.chart; });

    // Chart boundary = edges owned by exactly one of the chart's faces (mesh vertex pair key).
    std::map<std::pair<int, int>, std::pair<int, std::pair<Vec2f, Vec2f>>> edges;
    const auto                                                            &idx = s.atlas.mesh.indices;
    for (std::size_t f = 0; f < s.atlas.charts.faceToChart.size() && 3 * f + 2 < s.atlas.uv.size(); ++f) {
        if (s.atlas.charts.faceToChart[f] != pick.chart) {
            continue;
        }
        for (std::size_t k = 0; k < 3; ++k) {
            const int a = static_cast<int>(idx[3 * f + k]), b = static_cast<int>(idx[3 * f + (k + 1) % 3]);
            auto     &e = edges[{ std::min(a, b), std::max(a, b) }];
            ++e.first;
            e.second = { s.atlas.uv[3 * f + k], s.atlas.uv[3 * f + (k + 1) % 3] };
        }
    }
    for (const auto &[key, e] : edges) {
        if (e.first == 1) {
            ink.seg(e.second.first, e.second.second, 3.0f * kSS * uvPerPx, kBoundaryInk);
        }
    }
    ink.loop(gl, 3.0f * kSS * uvPerPx, kGutterInk);
    ink.square(g, 9.0f * kSS * uvPerPx, kGutterInk);

    glBindFramebuffer(GL_FRAMEBUFFER, ss.fbo);
    const int       big = n * static_cast<int>(kSS);
    const GlMat4    ortho = makeOrtho(ctr.x - half, ctr.x + half, ctr.y - half, ctr.y + half);
    const ui::Image img = renderImage(big, big, kShotBg, [&] {
        drawBatch(prog, ortho, fill);
        drawBatch(prog, ortho, wire);
        drawBatch(prog, ortho, ink);
    });
    return writeShot(path, ui::downscale(img, static_cast<int>(kSS)));
}


/// Seam-straightened atlas + bake panels (`_straightened`): same cutter, UVs post-processed.
static bool shotStraightened(const Programs &prog, ViewerState &s, const ShotTarget &t, const ShotTarget &big,
                             int sq, const std::string &base) {
    bool ok = true;
    s.straighten = true;
    computeChartification(s); // prints [straighten] stats, bakes seam-fill
    glBindFramebuffer(GL_FRAMEBUFFER, t.fbo);
    s.atlasShowNormalMap = false;
    ok &= shotPanel(base + "_atlas_straightened.tga", sq, sq, 0.82f,
                    [&] { renderUVPanel(prog, s, 0, 0, sq, sq); });
    shotBakePanels(prog, s, t, big, base, "_straightened", ok);
    s.straighten = false;
    return ok;
}


/// One LOD level of `shotLods`: the abs-normal RGB panel and its lit grayscale twin.
static bool shotLodLevel(const Programs &prog, ViewerState &s, const ShotTarget &t,
                         const std::string &path) {
    bool ok = true;
    for (const bool lit : { false, true }) {
        s.panel4LitGray = lit;
        ok &= shotPanel(path + (lit ? "_lit" : "") + ".tga", t.w, t.h, 1.0f,
                        [&] { renderBakedPanel(prog, s, 0, 0, t.w, t.h); });
    }
    s.panel4LitGray = false;
    return ok;
}


/// Seam-aware LOD shots: the straightened atlas with the erased texture, drawn full / 50% /
/// 10% (`_lod100`, `_lod50`, `_lod10`; `_lit` = lit grayscale).  The ONE texture is baked
/// once; the LOD levels only swap the mesh and its UVs.
static bool shotLods(const Programs &prog, ViewerState &s, const ShotTarget &t, const std::string &base) {
    s.straighten = true;
    computeChartification(s);
    s.gutterMode = decimate::Gutter::Erased;
    bakeAndUpload(s);
    computeLods(s);
    glBindFramebuffer(GL_FRAMEBUFFER, t.fbo);
    s.atlasShowNormalMap = false;
    bool ok = true;
    for (const int level : { 0, 1, 2 }) {
        s.lodLevel = level;
        const std::string suffix = level == 0 ? "_lod100" : "_lod" + std::to_string(kLodPercent[level - 1]);
        ok &= shotLodLevel(prog, s, t, base + suffix);
    }
    s.lodLevel = 0;
    s.straighten = false;
    s.gutterMode = decimate::Gutter::SeamFill;
    return ok;
}


/// Reduce → chartify → bake (as in the smoke path), then write each panel to `opt.outDir`.
// NOLINTNEXTLINE(readability-function-size)
static int runShots(const Programs &prog, ViewerState &s, const ShotOptions &opt,
                    const std::string &meshPath) {
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(opt.outDir, ec);
    const std::string base = (fs::path(opt.outDir) / fs::path(meshPath).stem()).string();

    for (int i = 0; i < 4; ++i) {
        decimateStep(s);
    }
    computeChartification(s); // also bakes (seam-fill)
    if (!s.charted || !s.baked) {
        std::fprintf(stderr, "[shot] chartify/bake failed\n");
        return 1;
    }

    // 3/4 view, framed so the model's projected extent fills 85% of the frame.
    s.azimuth = 145.0f;
    s.elevation = 20.0f;
    fitCameraToMesh(s, static_cast<float>(opt.width) / static_cast<float>(opt.height), 0.85f);

    constexpr int kBig = 3; // seam close-ups are cropped from a 3x render
    ShotTarget    t = makeShotTarget(opt.width, opt.height);
    ShotTarget    big = makeShotTarget(opt.width * kBig, opt.height * kBig);
    glBindFramebuffer(GL_FRAMEBUFFER, t.fbo);
    const int sq = std::min(t.w, t.h);
    bool      ok = true;
    ok &= shotPanel(base + "_mesh.tga", t.w, t.h, 1.0f,
                    [&] { render3DPanel(prog, s, 0, 0, t.w, t.h, /*useChartColors=*/false); });
    ok &= shotPanel(base + "_charts.tga", t.w, t.h, 1.0f,
                    [&] { render3DPanel(prog, s, 0, 0, t.w, t.h, /*useChartColors=*/true); });
    ok &= shotPanel(base + "_atlas.tga", sq, sq, 0.82f,
                    [&] { renderUVPanel(prog, s, 0, 0, sq, sq); });
    const LitShots litSeam = shotBakePanels(prog, s, t, big, base, "", ok);

    // Before/after: same charts and UVs, plain blur dilation instead of seam-correct fill.
    const decimate::bake::NormalMap seamMap = s.bakedMap;
    s.gutterMode = decimate::Gutter::Dilate;
    bakeAndUpload(s);
    const LitShots                  litDilate = shotBakePanels(prog, s, t, big, base, "_dilate", ok);
    const decimate::bake::NormalMap dilateMap = s.bakedMap;
    s.gutterMode = decimate::Gutter::SeamFill;

    // Texture-space close-up where the two gutters differ most (crisp nearest-neighbour).
    constexpr int kWin = 64, kScale = 10;
    const auto [zx, zy] = mostDifferentWindow(seamMap, s.bakedMap, kWin);
    const auto zoomSeam = cropZoom(seamMap, zx, zy, kWin, kScale);
    const auto zoomDilate = cropZoom(s.bakedMap, zx, zy, kWin, kScale);
    ok &= ui::writeTGA(base + "_normalmap_zoom.tga", zoomSeam);
    ok &= ui::writeTGA(base + "_normalmap_dilate_zoom.tga", zoomDilate);
    ok &= writeShot(base + "_seamfill_vs_dilate.tga",
                    ui::sideBySide(toImage(zoomSeam), toImage(zoomDilate),
                                   "seam-correct fill: real neighbour", "plain dilation"));
    std::fprintf(stderr, "[shot] zoom window %dx%d at texel (%d,%d), %dx\n", kWin, kWin, zx, zy, kScale);

    // Model-space question: is the difference visible on the lit 3D model at all?
    constexpr int kThresh = 12; // summed |dR|+|dG|+|dB| out of 765
    const auto    st = ui::diffStats(litSeam.frame, litDilate.frame, kThresh);
    std::fprintf(stderr, "[shot] lit model, seam-fill vs dilate: %ld of %d px differ by >= %d/765 (max %d)\n",
                 st.pixelsOver, litSeam.frame.w * litSeam.frame.h, kThresh, st.maxDiff);
    const auto bst = ui::diffStats(litSeam.big, litDilate.big, kThresh);
    if (bst.pixelsOver == 0) {
        std::fprintf(stderr, "[shot] no visible seam difference at model level; no seam close-up written\n");
    } else {
        constexpr int kCrop = 640;
        const auto [cx, cy] = ui::mostDifferentWindow(litSeam.big, litDilate.big, kCrop, 16, kThresh);
        const auto cropSeam = ui::cropImage(litSeam.big, cx, cy, kCrop, kCrop);
        const auto cropDilate = ui::cropImage(litDilate.big, cx, cy, kCrop, kCrop);
        ok &= writeShot(base + "_baked_lit_seam_zoom.tga", cropSeam);
        ok &= writeShot(base + "_baked_lit_dilate_zoom.tga", cropDilate);
        ok &= writeShot(base + "_baked_lit_seamfill_vs_dilate.tga",
                        ui::sideBySide(cropSeam, cropDilate, "seam-correct fill: real neighbour", "plain dilation"));
        std::fprintf(stderr, "[shot] 3D seam crop %dx%d at (%d,%d) of the %dx%d render\n", kCrop, kCrop, cx, cy,
                     litSeam.big.w, litSeam.big.h);
    }

    // Erased: seam::erase on top of the same dilate bake, so the only difference is the erasure.
    s.gutterMode = decimate::Gutter::Erased;
    bakeAndUpload(s);
    const LitShots litErased = shotBakePanels(prog, s, t, big, base, "_erased", ok);
    const auto [ex, ey] = mostDifferentWindow(dilateMap, s.bakedMap, kWin);
    ok &= ui::writeTGA(base + "_normalmap_erased_zoom.tga", cropZoom(s.bakedMap, ex, ey, kWin, kScale));
    ok &= writeShot(base + "_dilate_vs_erased.tga",
                    ui::sideBySide(toImage(cropZoom(dilateMap, ex, ey, kWin, kScale)),
                                   toImage(cropZoom(s.bakedMap, ex, ey, kWin, kScale)),
                                   "plain dilation", "seam erased"));
    const auto est = ui::diffStats(litDilate.frame, litErased.frame, kThresh);
    std::fprintf(stderr, "[shot] lit model, dilate vs erased: %ld of %d px differ by >= %d/765 (max %d)\n",
                 est.pixelsOver, litDilate.frame.w * litDilate.frame.h, kThresh, est.maxDiff);
    if (ui::diffStats(litDilate.big, litErased.big, kThresh).pixelsOver > 0) {
        constexpr int kCrop = 640;
        const auto [cx, cy] = ui::mostDifferentWindow(litDilate.big, litErased.big, kCrop, 16, kThresh);
        ok &= writeShot(base + "_baked_lit_dilate_vs_erased.tga",
                        ui::sideBySide(ui::cropImage(litDilate.big, cx, cy, kCrop, kCrop),
                                       ui::cropImage(litErased.big, cx, cy, kCrop, kCrop),
                                       "plain dilation", "seam erased"));
    }
    s.gutterMode = decimate::Gutter::SeamFill;

    // Gutter views: the packer's real gutterLoop (core), drawn as vector geometry, 2x supersampled.
    const auto packed = packedCharts(s);
    if (!packed.empty()) {
        constexpr int kCornerPx = 1000;
        ShotTarget    ss = makeShotTarget(2 * std::max(sq, kCornerPx), 2 * std::max(sq, kCornerPx));
        ok &= shotAtlasGutter(prog, s, packed, ss, sq, base + "_atlas_gutter.tga");
        ok &= shotGutterCorner(prog, s, packed, ss, kCornerPx, base + "_gutter_corner.tga");
        destroyShotTarget(ss);
    } else {
        ok = false;
    }

    ok &= shotLods(prog, s, t, base); // re-chartifies (straightened) and bakes the erased map

    if (opt.straighten) { // last: it re-chartifies the viewer state
        ok &= shotStraightened(prog, s, t, big, sq, base);
    }

    destroyShotTarget(big);
    destroyShotTarget(t);
    return ok ? 0 : 1;
}


// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
int runViewer(const std::string &meshPath, int smokeFrames, const ShotOptions *shot, bool lscm, int maxCharts,
              double maxStretchL2, double maxStretchLinf) {
    decimate::MeshData original;
    try {
        original = decimate::load(meshPath);
    } catch (const std::exception &ex) {
        std::fprintf(stderr, "mesh-viewer: failed to load '%s': %s\n"
                             "  (pass a mesh path, or run tools/fetch_models.sh)\n",
                     meshPath.c_str(), ex.what());
        return 1;
    }

    if (original.positions.empty()) {
        std::fprintf(stderr, "mesh-viewer: '%s' loaded with no vertices.\n",
                     meshPath.c_str());
        return 1;
    }

    const Aabb  box = computeAabb(original.positions);
    const float radius = aabbRadius(box);

    ViewerState state;
    state.maxCharts = maxCharts;
    state.maxStretchL2 = maxStretchL2;
    state.maxStretchLinf = maxStretchLinf;
    state.flattener = lscm ? decimate::chart::Flattener::Lscm : decimate::chart::Flattener::Bijective;
    state.original = original;
    state.current = original;
    state.originalVertexCount = static_cast<int>(original.positions.size());
    state.currentTarget = state.originalVertexCount;
    state.centreX = (box.lo.x + box.hi.x) * 0.5f;
    state.centreY = (box.lo.y + box.hi.y) * 0.5f;
    state.centreZ = (box.lo.z + box.hi.z) * 0.5f;
    state.distance = radius * 2.5f;

    std::fprintf(stderr, "[load] '%s'  verts=%d  tris=%d  radius=%.3f\n",
                 meshPath.c_str(),
                 state.originalVertexCount,
                 static_cast<int>(original.indices.size() / 3),
                 static_cast<double>(radius));

    gState = &state;

    GlfwWindow win;

    // Smoke and --shot runs must work headless, so they use a hidden 1x1 window instead of the visible one.
    if (smokeFrames > 0 || shot != nullptr) {
        glfwSetErrorCallback(glfwErrorCallback);
        if (glfwInit() == GLFW_FALSE) {
            std::fprintf(stderr, "mesh-viewer: glfwInit failed in smoke path\n");
            gState = nullptr;
            return 1;
        }
        glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
        // Core-profile 3.3 required even in the smoke path (we compile shaders).
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
        win.handle = glfwCreateWindow(1, 1, "mesh-viewer-smoke", nullptr, nullptr);
        if (win.handle == nullptr) {
            std::fprintf(stderr, "mesh-viewer: smoke window creation failed\n");
            glfwTerminate();
            gState = nullptr;
            return 1;
        }
        glfwMakeContextCurrent(win.handle);
    } else if (!win.open(1200, 800, "mesh-viewer")) {
        gState = nullptr;
        return 1;
    }

    // The callbacks find the viewer state through gState.
    glfwSetKeyCallback(win.handle, onKey);
    glfwSetMouseButtonCallback(win.handle, onMouseButton);
    glfwSetCursorPosCallback(win.handle, onCursorPos);
    glfwSetScrollCallback(win.handle, onScroll);

    // Log the GL version so the user can confirm core-profile 3.3 landed.
    std::fprintf(stderr, "[GL] version: %s\n",
                 // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
                 reinterpret_cast<const char *>(glGetString(GL_VERSION)));

    // Build shader programs — must happen after the GL context is current.
    const Programs prog = buildPrograms();

    glEnable(GL_DEPTH_TEST);

    if (smokeFrames > 0) {
        std::fprintf(stderr, "[smoke] step 0: verts=%d  tris=%d\n",
                     static_cast<int>(state.current.positions.size()),
                     static_cast<int>(state.current.indices.size() / 3));
        for (int i = 0; i < 4; ++i) {
            decimateStep(state);
        }
        computeChartification(state);

        const int triCount = static_cast<int>(state.current.indices.size() / 3);
        const int numCharts = static_cast<int>(state.atlas.charts.charts.size());
        const int numUV = static_cast<int>(state.atlas.uv.size());

        std::fprintf(stderr,
                     "[smoke] panels: mesh=%dtris charts=%d uv=%d\n",
                     triCount, numCharts, numUV);

        logBakeStatus(state);
    }

    int exitCode = 0;
    if (shot != nullptr) {
        exitCode = runShots(prog, state, *shot, meshPath);
    }

    int frameCount = 0;

    while (shot == nullptr && !win.shouldClose()) {
        int width = 0, height = 0;
        glfwGetFramebufferSize(win.handle, &width, &height);

        // A minimised window reports 0: clamp so the half-size viewports stay valid.
        if (height == 0) {
            height = 1;
        }
        if (width == 0) {
            width = 1;
        }

        const int hw = width / 2;
        const int hh = height / 2;

        // Each panel clears only its own quarter (scissor) to its own tint, then draws unscissored.
        // Panel ① top-left — dark blue-grey
        glViewport(0, hh, hw, hh);
        glScissor(0, hh, hw, hh);
        glEnable(GL_SCISSOR_TEST);
        glClearColor(0.10f, 0.11f, 0.14f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glDisable(GL_SCISSOR_TEST);
        render3DPanel(prog, state, 0, hh, hw, hh, /*useChartColors=*/false);

        // Panel ② top-right — dark purple-grey
        glViewport(hw, hh, hw, hh);
        glScissor(hw, hh, hw, hh);
        glEnable(GL_SCISSOR_TEST);
        glClearColor(0.10f, 0.09f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glDisable(GL_SCISSOR_TEST);
        render3DPanel(prog, state, hw, hh, hw, hh, /*useChartColors=*/true);

        // Panel ③ bottom-left — very dark green
        glViewport(0, 0, hw, hh);
        glScissor(0, 0, hw, hh);
        glEnable(GL_SCISSOR_TEST);
        glClearColor(0.08f, 0.11f, 0.09f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glDisable(GL_SCISSOR_TEST);
        renderUVPanel(prog, state, 0, 0, hw, hh);

        // Panel ④ bottom-right — dark amber
        glViewport(hw, 0, hw, hh);
        glScissor(hw, 0, hw, hh);
        glEnable(GL_SCISSOR_TEST);
        glClearColor(0.13f, 0.10f, 0.06f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glDisable(GL_SCISSOR_TEST);
        renderBakedPanel(prog, state, hw, 0, hw, hh);

        updateTitle(win.handle, state);
        win.swapAndPoll();

        ++frameCount;
        if (smokeFrames > 0 && frameCount >= smokeFrames) {
            break;
        }
    }

    gState = nullptr;

    // The context is still current here, so the programs can be deleted.
    glDeleteProgram(prog.mesh);
    glDeleteProgram(prog.ortho);
    glDeleteProgram(prog.quad);
    glDeleteProgram(prog.baked);

    if (smokeFrames > 0) {
        std::fprintf(stderr, "[smoke] done — %d frames rendered\n", frameCount);
    }

    return exitCode;
}

} // namespace ui
