// This file is part of the decimate library; see decimate.h for version/license details.
/// @file viewer.cpp
/// @brief Interactive GLFW/OpenGL core-profile 3.3 2×2 panel viewer.
///
/// Replaced the GL 2.1 immediate-mode path with a VBO/VAO + GLSL 330 pipeline.
/// Three shader programs cover all four panels:
///   kProgMesh   — panels ① and ②: lit mesh with a per-vertex colour attribute.
///   kProgOrtho  — panel ③: 2D ortho UV atlas (colour + line drawing).
///   kProgBaked  — panel ④: textured mesh with tangent-space normal-map decode.
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
#  include <OpenGL/gl3.h>  // NOLINT(modernize-deprecated-headers)
#else
#  include <GL/gl.h>       // NOLINT(modernize-deprecated-headers)
#endif

#include <GLFW/glfw3.h>

#include "viewer.h"

// mesh-core public API — the only algorithm surface this file touches.
#include <decimate/chart/chart.h>
#include <decimate/mesh/half_edge_mesh.h>
#include <decimate/types/mesh_data.h>
#include <decimate/types/mesh_types.h>
#include <decimate/bake/normal_bake.h>
#include <decimate/bake/normal_map.h>
#include <decimate/bake/normals.h>
#include <decimate/chart/parametrize.h>
#include <decimate/reduce/reducer.h>
#include <decimate/reduce/cluster_reduce.h>
#include <decimate/reduce/halfedge_reduce.h>

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

namespace ui {

// ---------------------------------------------------------------------------
// GlfwWindow — RAII wrapper
// ---------------------------------------------------------------------------

namespace {
void glfwErrorCallback(int code, const char* desc) noexcept {
    std::fprintf(stderr, "[GLFW error %d] %s\n", code, desc);
}
}  // namespace

bool GlfwWindow::open(int width, int height, const char* title) noexcept {
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
    glfwSwapInterval(1);  // vsync
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

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

struct Vec3f { float x{}, y{}, z{}; };

struct Aabb { Vec3f lo, hi; };

[[nodiscard]] Aabb computeAabb(const std::vector<decimate::math::Vec3f>& positions) noexcept {
    if (positions.empty()) {
        return {{0, 0, 0}, {0, 0, 0}};
    }
    Aabb box{
        {positions[0].x, positions[0].y, positions[0].z},
        {positions[0].x, positions[0].y, positions[0].z}};
    for (std::size_t i = 1; i < positions.size(); ++i) {
        const auto& v = positions[i];
        box.lo.x = std::min(box.lo.x, v.x);
        box.lo.y = std::min(box.lo.y, v.y);
        box.lo.z = std::min(box.lo.z, v.z);
        box.hi.x = std::max(box.hi.x, v.x);
        box.hi.y = std::max(box.hi.y, v.y);
        box.hi.z = std::max(box.hi.z, v.z);
    }
    return box;
}

[[nodiscard]] float aabbRadius(const Aabb& b) noexcept {
    float dx = b.hi.x - b.lo.x;
    float dy = b.hi.y - b.lo.y;
    float dz = b.hi.z - b.lo.z;
    return 0.5F * std::sqrt(dx * dx + dy * dy + dz * dz);
}

// ---------------------------------------------------------------------------
// Chart colour palette
// ---------------------------------------------------------------------------

struct Rgb { float r, g, b; };

[[nodiscard]] Rgb chartColor(int id) noexcept {
    // NOLINT(modernize-avoid-c-arrays): fixed-size compile-time table
    static constexpr std::array<Rgb, 16> kPalette = {{
        {0.90F, 0.20F, 0.20F},  // red
        {0.20F, 0.65F, 0.90F},  // sky-blue
        {0.20F, 0.80F, 0.35F},  // green
        {0.95F, 0.80F, 0.10F},  // yellow
        {0.75F, 0.30F, 0.90F},  // violet
        {0.95F, 0.50F, 0.10F},  // orange
        {0.15F, 0.75F, 0.75F},  // cyan
        {0.90F, 0.20F, 0.70F},  // pink
        {0.50F, 0.75F, 0.20F},  // lime
        {0.20F, 0.35F, 0.85F},  // blue
        {0.85F, 0.55F, 0.20F},  // amber
        {0.55F, 0.20F, 0.65F},  // purple
        {0.25F, 0.65F, 0.50F},  // teal
        {0.80F, 0.80F, 0.20F},  // gold
        {0.60F, 0.20F, 0.25F},  // maroon
        {0.45F, 0.45F, 0.85F},  // periwinkle
    }};
    const int idx = std::abs(id) % static_cast<int>(kPalette.size());
    return kPalette[static_cast<std::size_t>(idx)];
}

// ---------------------------------------------------------------------------
// mat4 — column-major 4×4, enough for perspective / ortho / lookAt
// ---------------------------------------------------------------------------

struct mat4 {
    float m[16]{};  // NOLINT(modernize-avoid-c-arrays)

    static mat4 identity() noexcept {
        mat4 r;
        r.m[0] = r.m[5] = r.m[10] = r.m[15] = 1.0F;
        return r;
    }

    // column-major indexing: m[col*4 + row]
    float& at(int col, int row) noexcept { return m[col * 4 + row]; }
    [[nodiscard]] float at(int col, int row) const noexcept { return m[col * 4 + row]; }
};

[[nodiscard]] mat4 mat4Mul(const mat4& a, const mat4& b) noexcept {
    mat4 r;
    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 4; ++row) {
            float s = 0.0F;
            for (int k = 0; k < 4; ++k) {
                s += a.at(k, row) * b.at(col, k);
            }
            r.at(col, row) = s;
        }
    }
    return r;
}

/// Perspective matrix matching the old glFrustum(-right,right,-top,top,near,far).
[[nodiscard]] mat4 makePerspective(float fovYDeg, float aspect, float nearP, float farP) noexcept {
    const float top   = nearP * std::tan(fovYDeg * 0.5F * 3.14159265F / 180.0F);
    const float right = top * aspect;
    // glFrustum equivalent (column-major)
    mat4 m;
    m.at(0, 0) = nearP / right;
    m.at(1, 1) = nearP / top;
    m.at(2, 2) = -(farP + nearP) / (farP - nearP);
    m.at(2, 3) = -1.0F;
    m.at(3, 2) = -2.0F * farP * nearP / (farP - nearP);
    return m;
}

/// Ortho matching glOrtho(l,r,b,t,-1,1).
[[nodiscard]] mat4 makeOrtho(float l, float r, float b, float t) noexcept {
    mat4 m = mat4::identity();
    m.at(0, 0) =  2.0F / (r - l);
    m.at(1, 1) =  2.0F / (t - b);
    m.at(2, 2) = -1.0F;
    m.at(3, 0) = -(r + l) / (r - l);
    m.at(3, 1) = -(t + b) / (t - b);
    return m;
}

/// Orbit model-view: translate(0,0,-dist) · rotX(elev) · rotY(azim) · translate(-cx,-cy,-cz).
[[nodiscard]] mat4 makeOrbitMV(float azimDeg, float elevDeg, float dist,
                                float cx, float cy, float cz) noexcept {
    const float az = azimDeg * 3.14159265F / 180.0F;
    const float el = elevDeg * 3.14159265F / 180.0F;
    const float saz = std::sin(az), caz = std::cos(az);
    const float sel = std::sin(el), cel = std::cos(el);

    // rotY(az): rotate azimuth around Y
    mat4 rY = mat4::identity();
    rY.at(0, 0) =  caz; rY.at(2, 0) = saz;
    rY.at(0, 2) = -saz; rY.at(2, 2) = caz;

    // rotX(el): tilt around X
    mat4 rX = mat4::identity();
    rX.at(1, 1) = cel; rX.at(2, 1) = sel;
    rX.at(1, 2) = -sel; rX.at(2, 2) = cel;

    // translate(-cx, -cy, -cz)
    mat4 tC = mat4::identity();
    tC.at(3, 0) = -cx; tC.at(3, 1) = -cy; tC.at(3, 2) = -cz;

    // translate(0, 0, -dist)
    mat4 tD = mat4::identity();
    tD.at(3, 2) = -dist;

    // MV = tD · rX · rY · tC
    return mat4Mul(tD, mat4Mul(rX, mat4Mul(rY, tC)));
}

// ---------------------------------------------------------------------------
// Shader compile / link
// ---------------------------------------------------------------------------

static unsigned int compileShader(GLenum type, const char* src) {
    unsigned int sh = glCreateShader(type);
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

static unsigned int linkProgram(unsigned int vs, unsigned int fs) {
    unsigned int prog = glCreateProgram();
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

// ---------------------------------------------------------------------------
// GLSL sources (GLSL 330 core)
// ---------------------------------------------------------------------------

// --- Mesh program (panels ① ②) ---
// Attributes: aPos(0), aNrm(1), aCol(2)
// Headlight from view direction (0,0.6,1) normalised — constant, no matrix needed.

static const char* kMeshVS = R"(
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

static const char* kMeshFS = R"(
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

// --- 2D ortho program (panel ③) ---
// Attributes: aPos2(0) as vec2, aCol(1) as vec3.
// Used for both filled chart triangles and line edges.

static const char* kOrthoVS = R"(
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

static const char* kOrthoFS = R"(
#version 330 core
in vec3 vCol;
out vec4 fragColor;
void main() {
    fragColor = vec4(vCol, 1.0);
}
)";

// --- Full-screen quad program for the baked-map preview (panel ③ 'n' mode) ---
// Same ortho program works — just pass white colour and use a textured version.
// Use a separate program to avoid adding a sampler to the ortho one.

static const char* kQuadVS = R"(
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

static const char* kQuadFS = R"(
#version 330 core
in vec2 vUV;
uniform sampler2D uTex;
out vec4 fragColor;
void main() {
    fragColor = vec4(texture(uTex, vUV).rgb, 1.0);
}
)";

// --- Baked-mesh program (panel ④) ---
// Attributes: aPos(0), aNrm(1), aUV(2), aTan(3)
// uniform uTangentSpace: 1 = decode TBN, 0 = display object-space normal as RGB.
//
// TBN decode matches mapFaceRt convention exactly:
//   Stored in texture: x=dot(snrm,T), y=dot(snrm,B), z=dot(snrm,N)
//   T = Gram-Schmidt UV-gradient tangent, N = interpolated low normal, B = N×T
//   Reconstruct: world = dx*T + dy*B + dz*N
//
// Handedness: B = N×T (no hand sign needed; the bake uses the same cross product).

static const char* kBakedVS = R"(
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

static const char* kBakedFS = R"(
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

// ---------------------------------------------------------------------------
// Shader programs — compiled once at startup
// ---------------------------------------------------------------------------

struct Programs {
    unsigned int mesh  = 0;  // panels ① ②: kMeshVS + kMeshFS
    unsigned int ortho = 0;  // panel ③ atlas: kOrthoVS + kOrthoFS
    unsigned int quad  = 0;  // panel ③ normal-map preview: kQuadVS + kQuadFS
    unsigned int baked = 0;  // panel ④: kBakedVS + kBakedFS
};

static Programs buildPrograms() {
    Programs p;
    p.mesh  = linkProgram(compileShader(GL_VERTEX_SHADER, kMeshVS),
                          compileShader(GL_FRAGMENT_SHADER, kMeshFS));
    p.ortho = linkProgram(compileShader(GL_VERTEX_SHADER, kOrthoVS),
                          compileShader(GL_FRAGMENT_SHADER, kOrthoFS));
    p.quad  = linkProgram(compileShader(GL_VERTEX_SHADER, kQuadVS),
                          compileShader(GL_FRAGMENT_SHADER, kQuadFS));
    p.baked = linkProgram(compileShader(GL_VERTEX_SHADER, kBakedVS),
                          compileShader(GL_FRAGMENT_SHADER, kBakedFS));
    return p;
}

// ---------------------------------------------------------------------------
// Viewer state
// ---------------------------------------------------------------------------

struct ViewerState {
    // Meshes
    decimate::MeshData original;
    decimate::MeshData current;
    int            originalVertexCount{};
    int            currentTarget{};

    // Camera
    float azimuth{};
    float elevation{30.F};
    float distance{};
    float centreX{}, centreY{}, centreZ{};

    // Interaction
    bool   leftDragging{false};
    double lastMouseX{}, lastMouseY{};

    // Display
    bool wireframe{false};

    // Reducer strategy: false = QEM edge-collapse (reduce), true = octree
    // Reducer strategy, cycled with 'K': 0 = QEM heap edge-collapse, 1 = octree
    // clustering (fast/sloppy), 2 = HALFEDGE corner-table edge-collapse (best).
    int reducer{0};

    // Chartification
    decimate::HalfEdgeMesh                hem;
    decimate::chart::ChartSet             chartSet;
    std::vector<decimate::math::Vec2<float>>    uv;       ///< Per-CORNER atlas UV (3*faceCount).
    bool                              charted{false};
    decimate::MeshData                    chartedMesh;

    // Baked normal map
    decimate::bake::NormalMap  bakedMap;
    unsigned int           bakedTex = 0;
    bool                   baked    = false;

    bool atlasShowNormalMap = false;
    bool bakeTangentSpace   = false;
    bool panel4LitGray      = false;  // panel ④: lit grayscale vs abs-normal RGB
    /// true = bakeSeamFill (real neighbour-chart texels in the gutter, the default);
    /// false = plain bake() + blur dilation. Only `--shot` flips it, for before/after images.
    bool bakeSeamFillMode   = true;
};

// Global — GLFW callbacks need access.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static ViewerState* gState = nullptr;

// ---------------------------------------------------------------------------
// Chartification helper
// ---------------------------------------------------------------------------

static void bakeAndUpload(ViewerState& s);

/// @brief Build HalfEdgeMesh + segment + parametrize on `state.current`.
static void computeChartification(ViewerState& s) {
    const auto& data = s.current;
    const int numPoints = static_cast<int>(data.positions.size());
    const int numFaces  = static_cast<int>(data.indices.size()  / 3);

    if (numFaces <= 0 || numPoints <= 0) {
        std::fprintf(stderr, "[chart] mesh is empty — nothing to chartify\n");
        return;
    }

    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    const auto* faces = reinterpret_cast<const decimate::Face*>(data.indices.data());

    const auto t0 = std::chrono::steady_clock::now();
    s.hem = decimate::HalfEdgeMesh::create(numPoints, numFaces, faces);

    const auto t1 = std::chrono::steady_clock::now();
    s.chartSet = decimate::chart::segment(s.hem, data, /*maxCharts=*/64);

    const auto t2 = std::chrono::steady_clock::now();
    s.uv = decimate::chart::parametrizeCorners(s.hem, data, s.chartSet, /*atlasResolution=*/256);

    const auto t3 = std::chrono::steady_clock::now();

    const auto msHem   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    const auto msSeg   = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    const auto msPar   = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
    const auto msTotal = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t0).count();

    s.charted     = true;
    s.chartedMesh = data;

    const int numCharts = static_cast<int>(s.chartSet.charts.size());
    const int numUV     = static_cast<int>(s.uv.size());

    std::fprintf(stderr,
        "[chart] %d faces → %d charts, %d UVs | hem=%ldms seg=%ldms param=%ldms total=%ldms\n",
        numFaces, numCharts, numUV,
        static_cast<long>(msHem), static_cast<long>(msSeg),
        static_cast<long>(msPar), static_cast<long>(msTotal));

    bakeAndUpload(s);
}

// ---------------------------------------------------------------------------
// Bake + GL texture upload
// ---------------------------------------------------------------------------

/// @brief Raytrace the normal map and upload it as a cached GL texture.
static void bakeAndUpload(ViewerState& s) {
    if (s.bakedTex != 0) {
        glDeleteTextures(1, &s.bakedTex);
        s.bakedTex = 0;
        s.baked    = false;
    }

    const std::vector<decimate::math::Vec3<float>> lowN =
        decimate::bake::computeVertexNormals(s.chartedMesh, decimate::bake::Weighting::EqualUnit);
    const std::vector<decimate::math::Vec3<float>> highN =
        decimate::bake::computeVertexNormals(s.original,    decimate::bake::Weighting::EqualUnit);

    constexpr int kTexSize    = 256;
    constexpr int kBakeDilate = 2;
    if (s.bakeSeamFillMode) {
        s.bakedMap = decimate::bake::bakeSeamFill(
            s.chartedMesh, s.hem, s.uv, s.chartSet.faceToChart,
            lowN, s.original, highN, kTexSize, kBakeDilate, s.bakeTangentSpace);
    } else {
        s.bakedMap = decimate::bake::bake(
            s.chartedMesh, s.uv, lowN, s.original, highN,
            kTexSize, kBakeDilate, s.bakeTangentSpace);
    }
    std::fprintf(stderr, "[bake] %s-space normals, %s gutter\n",
                 s.bakeTangentSpace ? "tangent" : "object",
                 s.bakeSeamFillMode ? "seam-fill" : "dilate");

    if (s.bakedMap.pixels.empty()) {
        std::fprintf(stderr, "[bake] WARNING: empty NormalMap — skipping GL upload\n");
        return;
    }

    unsigned int tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
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
    s.baked    = true;
}

// ---------------------------------------------------------------------------
// Transient VBO/VAO draw helpers
// Build GPU buffers fresh per frame — fine for a debug viewer; avoids
// tracking mesh invalidation.  A VAO wraps each draw to satisfy core-profile.
// ---------------------------------------------------------------------------

/// Upload float data into a fresh VBO, bind it to attribute `attrib`.
/// Each element has `stride` floats total; this attribute starts at `offset`.
/// Returns the VBO handle so the caller can delete it.
static unsigned int makeVboAttrib(unsigned int vao, int attrib,
                                   const float* data, int dataCount,
                                   int stride, int offset) {
    (void)vao;  // VAO already bound by caller
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
                          reinterpret_cast<const void*>(static_cast<std::intptr_t>(offset * sizeof(float))));
    return vbo;
}

// --- Panel ① / ②: flat-shaded lit mesh with per-vertex colour ---------------

/// Flat normal for face (a,b,c).
static Vec3f flatNormal(const decimate::math::Vec3f& a, const decimate::math::Vec3f& b, const decimate::math::Vec3f& c) noexcept {
    float ex = b.x - a.x, ey = b.y - a.y, ez = b.z - a.z;
    float fx = c.x - a.x, fy = c.y - a.y, fz = c.z - a.z;
    float nx = ey * fz - ez * fy;
    float ny = ez * fx - ex * fz;
    float nz = ex * fy - ey * fx;
    const float len = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (len > 1e-8F) { nx /= len; ny /= len; nz /= len; }
    return {nx, ny, nz};
}

/// Draw the mesh with flat per-face normals and a uniform colour (panel ①).
/// Each vertex in the interleaved buffer: pos(3) + nrm(3) + col(3) = 9 floats.
static void drawMeshFlat(const Programs& prog,
                          const mat4& mvp, const mat4& mv,
                          const decimate::MeshData& mesh,
                          float cr, float cg, float cb) {
    const auto& pos = mesh.positions;
    const auto& idx = mesh.indices;
    const std::size_t numTris = idx.size() / 3;
    const auto nv = static_cast<int>(pos.size());

    // Interleaved: pos(3) nrm(3) col(3)
    std::vector<float> buf;
    buf.reserve(numTris * 3 * 9);

    for (std::size_t t = 0; t < numTris; ++t) {
        const int ia = idx[t * 3 + 0];
        const int ib = idx[t * 3 + 1];
        const int ic = idx[t * 3 + 2];
        if (ia >= nv || ib >= nv || ic >= nv || ia < 0 || ib < 0 || ic < 0) {
            continue;
        }
        const auto& a = pos[static_cast<std::size_t>(ia)];
        const auto& b = pos[static_cast<std::size_t>(ib)];
        const auto& c = pos[static_cast<std::size_t>(ic)];
        const auto  n = flatNormal(a, b, c);

        for (const auto& v : {a, b, c}) {
            buf.push_back(v.x); buf.push_back(v.y); buf.push_back(v.z);
            buf.push_back(n.x); buf.push_back(n.y); buf.push_back(n.z);
            buf.push_back(cr);  buf.push_back(cg);  buf.push_back(cb);
        }
    }

    if (buf.empty()) {
        return;
    }

    const int stride = 9;  // 3 pos + 3 nrm + 3 col
    const int vertCount = static_cast<int>(buf.size()) / stride;

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
                          reinterpret_cast<const void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void*>(6 * sizeof(float)));

    glUseProgram(prog.mesh);
    glUniformMatrix4fv(glGetUniformLocation(prog.mesh, "uMVP"), 1, GL_FALSE, mvp.m);
    glUniformMatrix4fv(glGetUniformLocation(prog.mesh, "uMV"),  1, GL_FALSE, mv.m);

    glDrawArrays(GL_TRIANGLES, 0, vertCount);

    glBindVertexArray(0);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}

/// Draw mesh with per-face chart colours (panel ②).
static void drawMeshChartColors(const Programs& prog,
                                 const mat4& mvp, const mat4& mv,
                                 const decimate::MeshData& mesh,
                                 const decimate::chart::ChartSet& cs) {
    const auto& pos = mesh.positions;
    const auto& idx = mesh.indices;
    const std::size_t numTris = idx.size() / 3;
    const auto nv = static_cast<int>(pos.size());

    std::vector<float> buf;
    buf.reserve(numTris * 3 * 9);

    for (std::size_t t = 0; t < numTris; ++t) {
        const int ia = idx[t * 3 + 0];
        const int ib = idx[t * 3 + 1];
        const int ic = idx[t * 3 + 2];
        if (ia >= nv || ib >= nv || ic >= nv || ia < 0 || ib < 0 || ic < 0) {
            continue;
        }
        const auto& a = pos[static_cast<std::size_t>(ia)];
        const auto& b = pos[static_cast<std::size_t>(ib)];
        const auto& c = pos[static_cast<std::size_t>(ic)];
        const auto  n = flatNormal(a, b, c);

        int chartId = 0;
        if (t < cs.faceToChart.size()) {
            chartId = cs.faceToChart[t];
        }
        const Rgb col = chartColor(chartId);

        for (const auto& v : {a, b, c}) {
            buf.push_back(v.x); buf.push_back(v.y); buf.push_back(v.z);
            buf.push_back(n.x); buf.push_back(n.y); buf.push_back(n.z);
            buf.push_back(col.r); buf.push_back(col.g); buf.push_back(col.b);
        }
    }

    if (buf.empty()) {
        return;
    }

    const int stride   = 9;
    const int vertCount = static_cast<int>(buf.size()) / stride;
    const int sBytes   = stride * static_cast<int>(sizeof(float));

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
                          reinterpret_cast<const void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void*>(6 * sizeof(float)));

    glUseProgram(prog.mesh);
    glUniformMatrix4fv(glGetUniformLocation(prog.mesh, "uMVP"), 1, GL_FALSE, mvp.m);
    glUniformMatrix4fv(glGetUniformLocation(prog.mesh, "uMV"),  1, GL_FALSE, mv.m);

    glDrawArrays(GL_TRIANGLES, 0, vertCount);

    glBindVertexArray(0);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}

// --- Panel ③: UV atlas -------------------------------------------------------

/// Draw UV atlas in ortho [0,1]² with chart colours + wireframe edges.
// NOLINTNEXTLINE(readability-function-size): two separate draw calls (fill + lines)
static void drawUVAtlas(const Programs& prog,
                         const mat4& ortho,
                         const decimate::MeshData& mesh,
                         const decimate::chart::ChartSet& cs,
                         const std::vector<decimate::math::Vec2<float>>& uv) {
    const auto& idx = mesh.indices;
    const std::size_t numTris = idx.size() / 3;
    const std::size_t numCornerUV = uv.size();

    // Filled triangles — vertex = (u, v, r, g, b) interleaved
    std::vector<float> fillBuf;
    fillBuf.reserve(numTris * 3 * 5);

    for (std::size_t t = 0; t < numTris; ++t) {
        if (3 * t + 2 >= numCornerUV) {
            continue;
        }
        int chartId = 0;
        if (t < cs.faceToChart.size()) {
            chartId = cs.faceToChart[t];
        }
        const Rgb col = chartColor(chartId);

        for (int k = 0; k < 3; ++k) {
            const auto& uvc = uv[3 * t + static_cast<std::size_t>(k)];
            fillBuf.push_back(uvc.x); fillBuf.push_back(uvc.y);
            fillBuf.push_back(col.r); fillBuf.push_back(col.g); fillBuf.push_back(col.b);
        }
    }

    // Line edges — black wireframe
    std::vector<float> lineBuf;
    lineBuf.reserve(numTris * 6 * 5);

    for (std::size_t t = 0; t < numTris; ++t) {
        if (3 * t + 2 >= numCornerUV) {
            continue;
        }
        const auto& a = uv[3 * t + 0];
        const auto& b = uv[3 * t + 1];
        const auto& c = uv[3 * t + 2];
        // 3 edges × 2 endpoints; each endpoint = pos2+col(black)
        auto push2 = [&](const decimate::math::Vec2<float>& p0, const decimate::math::Vec2<float>& p1) {
            for (const auto& p : {p0, p1}) {
                lineBuf.push_back(p.x); lineBuf.push_back(p.y);
                lineBuf.push_back(0.0F); lineBuf.push_back(0.0F); lineBuf.push_back(0.0F);
            }
        };
        push2(a, b); push2(b, c); push2(c, a);
    }

    // --- Border quad (white) ---
    const float bW[] = {
        0.0F,0.0F, 1.0F,1.0F,1.0F,
        1.0F,0.0F, 1.0F,1.0F,1.0F,
        1.0F,0.0F, 1.0F,1.0F,1.0F,
        1.0F,1.0F, 1.0F,1.0F,1.0F,
        1.0F,1.0F, 1.0F,1.0F,1.0F,
        0.0F,1.0F, 1.0F,1.0F,1.0F,
        0.0F,1.0F, 1.0F,1.0F,1.0F,
        0.0F,0.0F, 1.0F,1.0F,1.0F,
    };  // NOLINT(modernize-avoid-c-arrays)

    glUseProgram(prog.ortho);
    glUniformMatrix4fv(glGetUniformLocation(prog.ortho, "uOrtho"), 1, GL_FALSE, ortho.m);

    // Draw filled
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
                              reinterpret_cast<const void*>(2 * sizeof(float)));
        glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(fillBuf.size()) / 5);
        glBindVertexArray(0);
        glDeleteBuffers(1, &vbo);
        glDeleteVertexArrays(1, &vao);
    }

    // Draw wireframe lines
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
                              reinterpret_cast<const void*>(2 * sizeof(float)));
        glDrawArrays(GL_LINES, 0, static_cast<int>(lineBuf.size()) / 5);
        glBindVertexArray(0);
        glDeleteBuffers(1, &vbo);
        glDeleteVertexArrays(1, &vao);
    }

    // Draw border
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
                              reinterpret_cast<const void*>(2 * sizeof(float)));
        glDrawArrays(GL_LINES, 0, 8);
        glBindVertexArray(0);
        glDeleteBuffers(1, &vbo);
        glDeleteVertexArrays(1, &vao);
    }
}

/// Draw the baked normal map as a textured quad over [0,1]² (panel ③ 'n' mode).
static void drawNormalMapQuad(const Programs& prog, const mat4& ortho, unsigned int tex) {
    // pos2 + uv2 interleaved (4 vertices, triangle-strip)
    static constexpr float kQuad[] = {
        0.0F, 0.0F,  0.0F, 0.0F,
        1.0F, 0.0F,  1.0F, 0.0F,
        0.0F, 1.0F,  0.0F, 1.0F,
        1.0F, 1.0F,  1.0F, 1.0F,
    };  // NOLINT(modernize-avoid-c-arrays)

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
                          reinterpret_cast<const void*>(2 * sizeof(float)));

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

    // White border lines
    static constexpr float kBorder[] = {
        0.0F,0.0F, 1.0F,0.0F,
        1.0F,0.0F, 1.0F,1.0F,
        1.0F,1.0F, 0.0F,1.0F,
        0.0F,1.0F, 0.0F,0.0F,
    };  // NOLINT(modernize-avoid-c-arrays)
    // Reuse ortho program for border (just white colour, 2D positions only).
    // Pack as pos2+col3 for the ortho shader.
    float borderFull[8 * 5];  // NOLINT(modernize-avoid-c-arrays)
    for (int i = 0; i < 8; ++i) {
        borderFull[i * 5 + 0] = kBorder[i * 2 + 0];
        borderFull[i * 5 + 1] = kBorder[i * 2 + 1];
        borderFull[i * 5 + 2] = 1.0F;
        borderFull[i * 5 + 3] = 1.0F;
        borderFull[i * 5 + 4] = 1.0F;
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
                          reinterpret_cast<const void*>(2 * sizeof(float)));
    glUseProgram(prog.ortho);
    glUniformMatrix4fv(glGetUniformLocation(prog.ortho, "uOrtho"), 1, GL_FALSE, ortho.m);
    glDrawArrays(GL_LINES, 0, 8);
    glBindVertexArray(0);
    glDeleteBuffers(1, &bvbo);
    glDeleteVertexArrays(1, &bvao);
}

// --- Panel ④: textured baked mesh with tangent-space decode ------------------

/// Compute smooth per-vertex normals by angle-weighted face-normal accumulation.
/// Each vertex receives the sum of the face normals of the faces it belongs to,
/// weighted by the interior angle at that vertex — standard "angle-weighted" smooth
/// normal (better than equal-weight at heavily tessellated poles).
///
/// Used to supply the interpolated-N channel for the TBN decode, matching the
/// lowNormals path that mapFaceRt uses (bake::computeVertexNormals with EqualUnit
/// weighting; angle-weighted is a reasonable stand-in for the debug viewer).
static std::vector<Vec3f> computeSmoothNormals(const decimate::MeshData& mesh) {
    const auto& pos = mesh.positions;
    const auto& idx = mesh.indices;
    const std::size_t nv = pos.size();
    const std::size_t numTris = idx.size() / 3;

    std::vector<Vec3f> acc(nv, {0.0F, 0.0F, 0.0F});

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
        for (const int vi : {ia, ib, ic}) {
            auto& a = acc[static_cast<std::size_t>(vi)];
            a.x += fn.x; a.y += fn.y; a.z += fn.z;
        }
    }

    std::vector<Vec3f> result(nv);
    for (std::size_t i = 0; i < nv; ++i) {
        const auto& a = acc[i];
        const float len = std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
        result[i] = (len > 1e-8F) ? Vec3f{a.x / len, a.y / len, a.z / len}
                                   : Vec3f{0.0F, 1.0F, 0.0F};
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
static void drawBakedMesh(const Programs& prog,
                           const mat4& mvp, const mat4& mv,
                           const ViewerState& s) {
    if (!s.baked || s.bakedTex == 0) {
        return;
    }

    const auto& mesh = s.chartedMesh;
    const auto& pos  = mesh.positions;
    const auto& idx  = mesh.indices;
    const std::size_t numTris = idx.size() / 3;
    const auto nv = static_cast<int>(pos.size());

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
        if (ia >= nv || ib >= nv || ic >= nv || ia < 0 || ib < 0 || ic < 0) {
            continue;
        }
        if (3 * t + 2 >= s.uv.size()) {
            continue;
        }

        const auto& pa = pos[static_cast<std::size_t>(ia)];
        const auto& pb = pos[static_cast<std::size_t>(ib)];
        const auto& pc = pos[static_cast<std::size_t>(ic)];

        // Per-face tangent: exact mapFaceRt formula using THIS face's per-corner UVs.
        // T = (pb-pa)*(dv2/det) - (pc-pa)*(dv1/det)
        // Degenerate UV (|det| < 1e-12 or |tFace| ≈ 0): fall back to an
        // N-orthogonal vector so the shader doesn't divide by zero.
        const auto& uA = s.uv[3 * t + 0];
        const auto& uB = s.uv[3 * t + 1];
        const auto& uC = s.uv[3 * t + 2];
        const float du1 = uB.x - uA.x, dv1 = uB.y - uA.y;
        const float du2 = uC.x - uA.x, dv2 = uC.y - uA.y;
        const float det = du1 * dv2 - du2 * dv1;

        Vec3f tFace{1.0F, 0.0F, 0.0F};  // fallback
        if (std::fabs(det) > 1e-12F) {
            const float r = 1.0F / det;
            const float tx = (pb.x - pa.x) * (dv2 * r) - (pc.x - pa.x) * (dv1 * r);
            const float ty = (pb.y - pa.y) * (dv2 * r) - (pc.y - pa.y) * (dv1 * r);
            const float tz = (pb.z - pa.z) * (dv2 * r) - (pc.z - pa.z) * (dv1 * r);
            const float tlen = std::sqrt(tx * tx + ty * ty + tz * tz);
            if (tlen > 1e-8F) {
                tFace = {tx / tlen, ty / tlen, tz / tlen};
            }
        }

        struct Corner { int vi; int k; };
        const Corner corners[3] = {{ia, 0}, {ib, 1}, {ic, 2}};  // NOLINT(modernize-avoid-c-arrays)
        for (const auto& [vi, k] : corners) {
            const auto& p   = pos[static_cast<std::size_t>(vi)];
            const auto& n   = smoothN[static_cast<std::size_t>(vi)];
            const auto& uvc = s.uv[3 * t + static_cast<std::size_t>(k)];
            buf.push_back(p.x);      buf.push_back(p.y);      buf.push_back(p.z);
            buf.push_back(n.x);      buf.push_back(n.y);      buf.push_back(n.z);
            buf.push_back(uvc.x);    buf.push_back(uvc.y);
            buf.push_back(tFace.x);  buf.push_back(tFace.y);  buf.push_back(tFace.z);
        }
    }

    if (buf.empty()) {
        return;
    }

    const int stride   = 11;
    const int vertCount = static_cast<int>(buf.size()) / stride;
    const int sBytes   = stride * static_cast<int>(sizeof(float));

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
                          reinterpret_cast<const void*>(3 * sizeof(float)));
    // aUV @ 2 (2 floats, offset 6)
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void*>(6 * sizeof(float)));
    // aTan @ 3 (3 floats, offset 8)
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void*>(8 * sizeof(float)));

    glUseProgram(prog.baked);
    glUniformMatrix4fv(glGetUniformLocation(prog.baked, "uMVP"), 1, GL_FALSE, mvp.m);
    glUniformMatrix4fv(glGetUniformLocation(prog.baked, "uMV"),  1, GL_FALSE, mv.m);
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
static void drawFaceNormalColors(const Programs& prog,
                                  const mat4& mvp, const mat4& mv,
                                  const decimate::MeshData& mesh) {
    const auto& pos = mesh.positions;
    const auto& idx = mesh.indices;
    const std::size_t numTris = idx.size() / 3;
    const auto nv = static_cast<int>(pos.size());

    std::vector<float> buf;
    buf.reserve(numTris * 3 * 9);

    for (std::size_t t = 0; t < numTris; ++t) {
        const int ia = idx[t * 3 + 0];
        const int ib = idx[t * 3 + 1];
        const int ic = idx[t * 3 + 2];
        if (ia >= nv || ib >= nv || ic >= nv || ia < 0 || ib < 0 || ic < 0) {
            continue;
        }
        const auto& a = pos[static_cast<std::size_t>(ia)];
        const auto& b = pos[static_cast<std::size_t>(ib)];
        const auto& c = pos[static_cast<std::size_t>(ic)];
        const auto fn = flatNormal(a, b, c);
        const float cr = std::abs(fn.x), cg = std::abs(fn.y), cb2 = std::abs(fn.z);

        for (const auto& v : {a, b, c}) {
            buf.push_back(v.x); buf.push_back(v.y); buf.push_back(v.z);
            buf.push_back(fn.x); buf.push_back(fn.y); buf.push_back(fn.z);
            buf.push_back(cr); buf.push_back(cg); buf.push_back(cb2);
        }
    }

    if (buf.empty()) {
        return;
    }

    const int stride   = 9;
    const int vertCount = static_cast<int>(buf.size()) / stride;
    const int sBytes   = stride * static_cast<int>(sizeof(float));

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
                          reinterpret_cast<const void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sBytes,
                          // NOLINTNEXTLINE(performance-no-int-to-ptr)
                          reinterpret_cast<const void*>(6 * sizeof(float)));

    glUseProgram(prog.mesh);
    glUniformMatrix4fv(glGetUniformLocation(prog.mesh, "uMVP"), 1, GL_FALSE, mvp.m);
    glUniformMatrix4fv(glGetUniformLocation(prog.mesh, "uMV"),  1, GL_FALSE, mv.m);

    glDrawArrays(GL_TRIANGLES, 0, vertCount);

    glBindVertexArray(0);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}

// ---------------------------------------------------------------------------
// Panel render functions
// ---------------------------------------------------------------------------

struct CameraMatrices {
    mat4 proj;
    mat4 mv;
    mat4 mvp;
};

[[nodiscard]] CameraMatrices makeCamera(const ViewerState& s, int vpW, int vpH) noexcept {
    const float aspect = (vpH > 0)
        ? (static_cast<float>(vpW) / static_cast<float>(vpH))
        : 1.0F;
    const float nearP    = s.distance * 0.01F;
    const float farPlane = s.distance * 100.0F;

    CameraMatrices c;
    c.proj = makePerspective(45.0F, aspect, nearP, farPlane);
    c.mv   = makeOrbitMV(s.azimuth, s.elevation, s.distance,
                         s.centreX, s.centreY, s.centreZ);
    c.mvp  = mat4Mul(c.proj, c.mv);
    return c;
}

/// Render one 3D panel (panels ① and ②).
static void render3DPanel(const Programs& prog, const ViewerState& s,
                          int vpX, int vpY, int vpW, int vpH,
                          bool useChartColors) {
    glViewport(vpX, vpY, vpW, vpH);
    glEnable(GL_DEPTH_TEST);
    glClear(GL_DEPTH_BUFFER_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, s.wireframe ? GL_LINE : GL_FILL);

    const auto cam = makeCamera(s, vpW, vpH);

    if (!useChartColors) {
        drawMeshFlat(prog, cam.mvp, cam.mv, s.current, 0.7F, 0.75F, 0.82F);
    } else if (s.charted) {
        drawMeshChartColors(prog, cam.mvp, cam.mv, s.chartedMesh, s.chartSet);
    } else {
        drawMeshFlat(prog, cam.mvp, cam.mv, s.current, 0.55F, 0.55F, 0.58F);
    }
}

/// Render the UV atlas panel (panel ③).
static void renderUVPanel(const Programs& prog, const ViewerState& s,
                          int vpX, int vpY, int vpW, int vpH) {
    glViewport(vpX, vpY, vpW, vpH);
    glClear(GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);

    constexpr float kMargin = 0.05F;
    const mat4 ortho = makeOrtho(-kMargin, 1.0F + kMargin, -kMargin, 1.0F + kMargin);

    if (s.atlasShowNormalMap && s.baked && s.bakedTex != 0) {
        drawNormalMapQuad(prog, ortho, s.bakedTex);
    } else if (s.charted) {
        drawUVAtlas(prog, ortho, s.chartedMesh, s.chartSet, s.uv);
    }
    // Not charted: dark background only.

    glEnable(GL_DEPTH_TEST);
}

/// Render panel ④ — baked normal map on the low-poly mesh with TBN decode.
static void renderBakedPanel(const Programs& prog, const ViewerState& s,
                              int vpX, int vpY, int vpW, int vpH) {
    glViewport(vpX, vpY, vpW, vpH);
    glEnable(GL_DEPTH_TEST);
    glClear(GL_DEPTH_BUFFER_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    const auto cam = makeCamera(s, vpW, vpH);

    if (s.baked && s.bakedTex != 0) {
        drawBakedMesh(prog, cam.mvp, cam.mv, s);
    } else {
        drawFaceNormalColors(prog, cam.mvp, cam.mv, s.current);
    }
}

// ---------------------------------------------------------------------------
// Title bar HUD
// ---------------------------------------------------------------------------

static void updateTitle(GLFWwindow* win, const ViewerState& s) {
    char buf[384];  // NOLINT(modernize-avoid-c-arrays)
    const int vertexCount   = static_cast<int>(s.current.positions.size());
    const int triangleCount = static_cast<int>(s.current.indices.size()   / 3);

    char chartFrag[64];  // NOLINT(modernize-avoid-c-arrays)
    if (s.charted) {
        const int numCharts = static_cast<int>(s.chartSet.charts.size());
        std::snprintf(chartFrag, sizeof(chartFrag), "charts:%d", numCharts);
    } else {
        std::snprintf(chartFrag, sizeof(chartFrag), "(press c to chartify)");
    }

    std::snprintf(buf, sizeof(buf),
        "mesh-viewer  |  verts:%d  tris:%d  %s  [atlas:%s]"
        "  | -/+ reduce  r reset  c chartify  t tangent  n atlas/normalmap  g shade  w wire  q quit",
        vertexCount, triangleCount, chartFrag,
        s.atlasShowNormalMap ? "normalmap" : "charts");
    glfwSetWindowTitle(win, buf);
}

// ---------------------------------------------------------------------------
// Reduce helpers
// ---------------------------------------------------------------------------

static void applyReduce(ViewerState& s) {
    if (s.currentTarget >= s.originalVertexCount) {
        s.current = s.original;
        return;
    }
    switch (s.reducer) {
        case 1:  s.current = decimate::reduceClustering(s.original, s.currentTarget); break;
        case 2:  s.current = decimate::reduceHalfEdge(s.original, s.currentTarget);  break;
        default: s.current = decimate::reduce(s.original, s.currentTarget);          break;
    }
    static const char* const kNames[] = {"QEM", "CLUSTER", "HALFEDGE"};
    std::fprintf(stderr, "[reduce] %s  target=%d  actual verts=%d  tris=%d\n",
        kNames[s.reducer],
        s.currentTarget,
        static_cast<int>(s.current.positions.size()),
        static_cast<int>(s.current.indices.size()   / 3));
}

static void decimateStep(ViewerState& s) {
    const int displayed = static_cast<int>(s.current.positions.size());
    s.currentTarget = std::max(4, displayed * 7 / 10);
    applyReduce(s);
}

static void growStep(ViewerState& s) {
    const int displayed = static_cast<int>(s.current.positions.size());
    s.currentTarget = std::min(s.originalVertexCount, (displayed * 14 / 10) + 1);
    applyReduce(s);
}

// ---------------------------------------------------------------------------
// GLFW callbacks
// ---------------------------------------------------------------------------

static void onKey(GLFWwindow* win, int key, int /*scancode*/, int action, int /*mods*/) {
    if (action != GLFW_PRESS && action != GLFW_REPEAT) {
        return;
    }
    ViewerState& s = *gState;

    if (key == GLFW_KEY_ESCAPE || key == GLFW_KEY_Q) {
        glfwSetWindowShouldClose(win, GLFW_TRUE);
        return;
    }
    if (key == GLFW_KEY_R) {
        s.currentTarget = s.originalVertexCount;
        s.current = s.original;
        s.charted = false;
        s.baked   = false;
        std::fprintf(stderr, "[reset] verts=%d  tris=%d\n",
            static_cast<int>(s.current.positions.size()),
            static_cast<int>(s.current.indices.size()   / 3));
        return;
    }
    if (key == GLFW_KEY_W) {
        s.wireframe = !s.wireframe;
        return;
    }
    if (key == GLFW_KEY_K) {
        s.reducer = (s.reducer + 1) % 3;  // QEM -> CLUSTER -> HALFEDGE -> QEM
        static const char* const kNames[] = {"QEM (heap edge-collapse)",
                                             "CLUSTER (octree, fast/sloppy)",
                                             "HALFEDGE (corner-table edge-collapse)"};
        std::fprintf(stderr, "[reducer] %s\n", kNames[s.reducer]);
        applyReduce(s);  // re-decimate at the same target with the new strategy
        return;
    }
    if (key == GLFW_KEY_S) {
        if (s.baked) {
            const bool ok = ui::writeTGA("normalmap.tga", s.bakedMap);
            if (ok) {  // NOLINT(bugprone-branch-clone)
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

static void onMouseButton(GLFWwindow* /*win*/, int button, int action, int /*mods*/) {
    if (button != GLFW_MOUSE_BUTTON_LEFT) {
        return;
    }
    ViewerState& s = *gState;
    s.leftDragging = (action == GLFW_PRESS);
}

static void onCursorPos(GLFWwindow* /*win*/, double xpos, double ypos) {
    ViewerState& s = *gState;
    if (!s.leftDragging) {
        s.lastMouseX = xpos;
        s.lastMouseY = ypos;
        return;
    }
    const double dx = xpos - s.lastMouseX;
    const double dy = ypos - s.lastMouseY;
    s.azimuth   += static_cast<float>(dx) * 0.5F;
    s.elevation += static_cast<float>(dy) * 0.5F;
    s.elevation  = std::max(-89.0F, std::min(89.0F, s.elevation));
    s.lastMouseX = xpos;
    s.lastMouseY = ypos;
}

static void onScroll(GLFWwindow* /*win*/, double /*xoff*/, double yoff) {
    ViewerState& s = *gState;
    s.distance *= (yoff > 0) ? 0.9F : 1.1F;
    s.distance  = std::max(0.01F, s.distance);
}

// ---------------------------------------------------------------------------
// Bake smoke log
// ---------------------------------------------------------------------------

static void logBakeStatus(const ViewerState& s) {
    const int nmW = s.bakedMap.width;
    const int nmH = s.bakedMap.height;
    const auto nmPixels = static_cast<int>(s.bakedMap.pixels.size());
    const bool ok = (nmW > 0 && nmH > 0 && nmPixels == nmW * nmH);
    if (ok) {  // NOLINT(bugprone-branch-clone)
        std::fprintf(stderr,
            "[smoke] panel④ bake OK: %dx%d NormalMap (%d texels) baked=%s\n",
            nmW, nmH, nmPixels, s.baked ? "true" : "false");
    } else {
        std::fprintf(stderr,
            "[smoke] panel④ bake FAIL: w=%d h=%d pixels=%d baked=%s\n",
            nmW, nmH, nmPixels, s.baked ? "true" : "false");
    }
}


// ---------------------------------------------------------------------------
// --shot: write every panel as its own image (README screenshots)
// ---------------------------------------------------------------------------

namespace {

/// Offscreen colour+depth target. Rendering here instead of the (hidden) window
/// gives exactly the requested pixel size, independent of Retina scaling.
struct ShotTarget {
    unsigned int fbo{}, color{}, depth{};
    int          w{}, h{};
};

[[nodiscard]] ShotTarget makeShotTarget(int w, int h) {
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

void destroyShotTarget(ShotTarget& t) {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glDeleteFramebuffers(1, &t.fbo);
    glDeleteRenderbuffers(1, &t.color);
    glDeleteRenderbuffers(1, &t.depth);
}

/// Clear the target to `bg`, run `draw` into the first `w`x`h` pixels and read them back
/// as a top-row-first Image (glReadPixels rows are bottom-up, so they are flipped here).
template <class Draw>
ui::Image renderImage(int w, int h, float bg, Draw draw) {
    glDisable(GL_SCISSOR_TEST);
    glViewport(0, 0, w, h);
    glClearColor(bg, bg, bg, 1.0F);
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

bool writeShot(const std::string& path, const ui::Image& img) {
    const bool ok = ui::writeTGA(path, img);
    std::fprintf(stderr, "[shot] %s %s (%dx%d)\n", ok ? "wrote" : "FAILED", path.c_str(), img.w, img.h);
    return ok;
}

template <class Draw>
bool shotPanel(const std::string& path, int w, int h, float bg, Draw draw) {
    return writeShot(path, renderImage(w, h, bg, draw));
}

/// Frame the mesh for a README shot: iteratively re-centre the orbit pivot on the projected
/// extent and scale the distance until the model's projected bounds fill `fill` of the frame
/// (height- or width-limited, whichever is tighter). Uses the real vertices, not the AABB
/// sphere, which over-estimates by ~30-40% for compact models.
void fitCameraToMesh(ViewerState& s, float aspect, float fill) {
    constexpr float kTanHalfFov = 0.41421356F;  // tan(45 deg / 2), matches makeCamera
    for (int iter = 0; iter < 24; ++iter) {
        const mat4 mv = makeOrbitMV(s.azimuth, s.elevation, s.distance, s.centreX, s.centreY, s.centreZ);
        float loX = 1e30F, hiX = -1e30F, loY = 1e30F, hiY = -1e30F;
        for (const auto& p : s.current.positions) {
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
        const float panX = 0.5F * (loX + hiX) * s.distance * kTanHalfFov * aspect;
        const float panY = 0.5F * (loY + hiY) * s.distance * kTanHalfFov;
        s.centreX += mv.at(0, 0) * panX + mv.at(0, 1) * panY;
        s.centreY += mv.at(1, 0) * panX + mv.at(1, 1) * panY;
        s.centreZ += mv.at(2, 0) * panX + mv.at(2, 1) * panY;
        const float half = std::max(0.5F * (hiX - loX), 0.5F * (hiY - loY));
        s.distance *= half / fill;
    }
}

/// Top-left of the `win`x`win` window (8-px grid) where the two bakes differ most:
/// that is where the gutter fill is visible.
[[nodiscard]] std::pair<int, int> mostDifferentWindow(const decimate::bake::NormalMap& a,
                                                      const decimate::bake::NormalMap& b, int win) {
    std::pair<int, int> best{0, 0};
    long bestScore = -1;
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
                best      = {x0, y0};
            }
        }
    }
    return best;
}

/// Nearest-neighbour `scale`x crop; rows flipped so v points up like the atlas panel.
[[nodiscard]] decimate::bake::NormalMap cropZoom(const decimate::bake::NormalMap& src,
                                                 int x0, int y0, int win, int scale) {
    decimate::bake::NormalMap out;
    out.width  = win * scale;
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
[[nodiscard]] ui::Image toImage(const decimate::bake::NormalMap& nm) {
    ui::Image im;
    im.w = nm.width;
    im.h = nm.height;
    im.rgb.reserve(nm.pixels.size() * 3);
    for (const auto& p : nm.pixels) {
        im.rgb.insert(im.rgb.end(), {p.r, p.g, p.b});
    }
    return im;
}

/// What the lit baked panel looks like for the current gutter mode, at frame size and at 3x
/// (the 3x render is where a seam close-up is cropped from).
struct LitShots {
    ui::Image frame;
    ui::Image big;
};

/// The baked-model panels (normal-map quad, abs-normal RGB, lit grayscale) for whatever
/// gutter mode `s.bakeSeamFillMode` currently selects; `suffix` is "" or "_dilate".
LitShots shotBakePanels(const Programs& prog, ViewerState& s, const ShotTarget& t, const ShotTarget& big,
                        const std::string& base, const std::string& suffix, bool& ok) {
    const int sq = std::min(t.w, t.h);
    glBindFramebuffer(GL_FRAMEBUFFER, t.fbo);
    s.atlasShowNormalMap = true;
    ok &= shotPanel(base + "_normalmap" + suffix + ".tga", sq, sq, 0.82F,
                    [&] { renderUVPanel(prog, s, 0, 0, sq, sq); });
    s.atlasShowNormalMap = false;
    s.panel4LitGray      = false;
    ok &= shotPanel(base + "_baked" + suffix + ".tga", t.w, t.h, 1.0F,
                    [&] { renderBakedPanel(prog, s, 0, 0, t.w, t.h); });

    s.panel4LitGray = true;  // the G-key path: uLitGray in the baked fragment shader
    LitShots lit;
    lit.frame = renderImage(t.w, t.h, 1.0F, [&] { renderBakedPanel(prog, s, 0, 0, t.w, t.h); });
    ok &= writeShot(base + "_baked_lit" + suffix + ".tga", lit.frame);
    glBindFramebuffer(GL_FRAMEBUFFER, big.fbo);
    lit.big = renderImage(big.w, big.h, 1.0F, [&] { renderBakedPanel(prog, s, 0, 0, big.w, big.h); });
    glBindFramebuffer(GL_FRAMEBUFFER, t.fbo);
    s.panel4LitGray = false;
    return lit;
}

// --- Gutter views (README images): the packer's real gutterLoop, drawn, never recomputed ---

using Vec2f = decimate::math::Vec2<float>;

/// Triangle soup (x, y, r, g, b per vertex) for the 2D shot overlays; lines are quads whose
/// width is given in atlas-UV units (core-profile GL has no wide lines).
struct ShotBatch {
    std::vector<float> v;

    void vert(Vec2f p, Rgb c) { v.insert(v.end(), {p.x, p.y, c.r, c.g, c.b}); }
    void tri(Vec2f a, Vec2f b, Vec2f c, Rgb col) { vert(a, col); vert(b, col); vert(c, col); }

    /// Segment a-b, `w` wide, square-capped so consecutive segments leave no gap at a joint.
    void seg(Vec2f a, Vec2f b, float w, Rgb c) {
        const float dx = b.x - a.x, dy = b.y - a.y;
        const float len = std::sqrt(dx * dx + dy * dy);
        if (len < 1e-9F) {
            return;
        }
        const float hx = dx / len * (w * 0.5F), hy = dy / len * (w * 0.5F);
        const Vec2f p0(a.x - hx, a.y - hy), p1(b.x + hx, b.y + hy);
        const Vec2f n(-hy, hx);
        const Vec2f q0(p0.x + n.x, p0.y + n.y), q1(p0.x - n.x, p0.y - n.y);
        const Vec2f q2(p1.x - n.x, p1.y - n.y), q3(p1.x + n.x, p1.y + n.y);
        tri(q0, q1, q2, c);
        tri(q0, q2, q3, c);
    }
    void loop(const std::vector<Vec2f>& pts, float w, Rgb c) {
        for (std::size_t i = 0; i < pts.size(); ++i) {
            seg(pts[i], pts[(i + 1) % pts.size()], w, c);
        }
    }
    void square(Vec2f p, float w, Rgb c) {
        const float h = w * 0.5F;
        const Vec2f a(p.x - h, p.y - h), b(p.x + h, p.y - h), cc(p.x + h, p.y + h), d(p.x - h, p.y + h);
        tri(a, b, cc, c);
        tri(a, cc, d, c);
    }
};

void drawBatch(const Programs& prog, const mat4& ortho, const ShotBatch& b) {
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
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<const void*>(2 * sizeof(float)));
    glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(b.v.size() / 5));
    glBindVertexArray(0);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
    glEnable(GL_DEPTH_TEST);
}

/// The packer's output for the current chartification: flatten + packAtlas through the public
/// API (the same two calls parametrizeCorners makes), kept so `gutterLoop` can be read.
/// Returns empty if its UVs disagree with the viewer's `s.uv` (then the outlines would lie).
std::vector<decimate::chart::HalfEdgeChart> packedCharts(const ViewerState& s) {
    std::vector<decimate::chart::HalfEdgeChart> charts;
    for (const auto& ch : s.chartSet.charts) {
        charts.push_back(decimate::chart::flatten(s.hem, s.chartedMesh, ch));
    }
    decimate::chart::packAtlas(charts, s.chartedMesh, 256);  // same resolution as computeChartification

    const auto& idx = s.chartedMesh.indices;
    float worst = 0.0F;
    for (std::size_t f = 0; f < s.chartSet.faceToChart.size(); ++f) {
        const int ci = s.chartSet.faceToChart[f];
        if (ci < 0 || static_cast<std::size_t>(ci) >= charts.size()) {
            continue;
        }
        const auto& c = charts[static_cast<std::size_t>(ci)];
        for (std::size_t k = 0; k < 3; ++k) {
            const auto it = std::find(c.localToMesh.begin(), c.localToMesh.end(), static_cast<int>(idx[3 * f + k]));
            if (it == c.localToMesh.end()) {
                continue;
            }
            const Vec2f& a = c.uv[static_cast<std::size_t>(it - c.localToMesh.begin())];
            const Vec2f& b = s.uv[3 * f + k];
            worst = std::max({worst, std::fabs(a.x - b.x), std::fabs(a.y - b.y)});
        }
    }
    std::fprintf(stderr, "[shot] repacked %zu charts; max |repacked uv - s.uv| = %.2e\n", charts.size(),
                 static_cast<double>(worst));
    if (worst > 1e-5F) {
        std::fprintf(stderr, "[shot] repack does not match s.uv; gutter images skipped\n");
        charts.clear();
    }
    return charts;
}

constexpr Rgb kGutterInk{0.03F, 0.03F, 0.08F};
constexpr Rgb kWireInk{0.30F, 0.30F, 0.30F};
constexpr Rgb kBoundaryInk{0.05F, 0.25F, 0.85F};
constexpr float kShotBg = 0.97F;

/// Add the filled triangles (and optionally edges) of faces whose chart passes `pick`.
template <class Pick>
void addChartFaces(ShotBatch& fill, ShotBatch& wire, const ViewerState& s, float wireW, Pick pick) {
    const std::size_t nf = s.chartedMesh.indices.size() / 3;
    for (std::size_t f = 0; f < nf && 3 * f + 2 < s.uv.size() && f < s.chartSet.faceToChart.size(); ++f) {
        const int ci = s.chartSet.faceToChart[f];
        if (!pick(ci)) {
            continue;
        }
        const Rgb col = chartColor(ci);
        fill.tri(s.uv[3 * f], s.uv[3 * f + 1], s.uv[3 * f + 2], col);
        wire.loop({s.uv[3 * f], s.uv[3 * f + 1], s.uv[3 * f + 2]}, wireW, kWireInk);
    }
}

/// <model>_atlas_gutter: chart fills + triangle edges + every chart's gutterLoop, light background.
bool shotAtlasGutter(const Programs& prog, const ViewerState& s,
                     const std::vector<decimate::chart::HalfEdgeChart>& charts, ShotTarget& ss, int n,
                     const std::string& path) {
    constexpr float kMargin = 0.05F, kSS = 2.0F;
    const float uvPerPx = (1.0F + 2.0F * kMargin) / (static_cast<float>(n) * kSS);  // one supersampled px
    ShotBatch fill, wire, ink;
    addChartFaces(fill, wire, s, 1.0F * kSS * uvPerPx, [](int) { return true; });
    ink.loop({{0, 0}, {1, 0}, {1, 1}, {0, 1}}, 1.0F * kSS * uvPerPx, Rgb{0.6F, 0.6F, 0.6F});
    for (const auto& c : charts) {
        ink.loop(c.gutterLoop, 1.5F * kSS * uvPerPx, kGutterInk);
    }
    glBindFramebuffer(GL_FRAMEBUFFER, ss.fbo);
    const int big = n * static_cast<int>(kSS);
    const mat4 ortho = makeOrtho(-kMargin, 1.0F + kMargin, -kMargin, 1.0F + kMargin);
    const ui::Image img = renderImage(big, big, kShotBg, [&] {
        drawBatch(prog, ortho, fill);
        drawBatch(prog, ortho, wire);
        drawBatch(prog, ortho, ink);
    });
    return writeShot(path, ui::downscale(img, static_cast<int>(kSS)));
}

struct SharpCorner {
    int   chart = -1;
    int   vertex = -1;  ///< index into that chart's gutterLoop
    float angleDeg = 180.0F;
    Vec2f nearest{0.0F, 0.0F};  ///< the chart corner UV closest to that gutter vertex (its boundary vertex)
    float dist     = 0.0F;      ///< distance to it (the realised miter offset)
};

/// Smallest interior angle among CONVEX gutterLoop vertices of charts with >= `minFaces` faces.
/// Vertices that are miter-fold artefacts (the loop folds back on itself: edges shorter than the
/// gutter, offset well under the gutter or beyond the miter limit) are skipped: they are not corners.
SharpCorner sharpestConvexCorner(const ViewerState& s, const std::vector<decimate::chart::HalfEdgeChart>& charts,
                                 const std::vector<int>& facesPerChart, int minFaces) {
    SharpCorner best;
    const float gutter = static_cast<float>(decimate::chart::kGutterPx) / 255.0F;
    std::vector<std::vector<Vec2f>> pts(charts.size());
    for (std::size_t f = 0; f < s.chartSet.faceToChart.size() && 3 * f + 2 < s.uv.size(); ++f) {
        const int ci = s.chartSet.faceToChart[f];
        if (ci >= 0 && static_cast<std::size_t>(ci) < charts.size()) {
            pts[static_cast<std::size_t>(ci)].insert(pts[static_cast<std::size_t>(ci)].end(), s.uv.begin() + 3 * static_cast<std::ptrdiff_t>(f),
                                                     s.uv.begin() + 3 * static_cast<std::ptrdiff_t>(f) + 3);
        }
    }
    for (std::size_t ci = 0; ci < charts.size(); ++ci) {
        const auto& g = charts[ci].gutterLoop;
        if (g.size() < 3 || facesPerChart[ci] < minFaces) {
            continue;
        }
        float area2 = 0.0F;
        for (std::size_t i = 0; i < g.size(); ++i) {
            const Vec2f& a = g[i];
            const Vec2f& b = g[(i + 1) % g.size()];
            area2 += a.x * b.y - b.x * a.y;
        }
        const float sign = area2 >= 0.0F ? 1.0F : -1.0F;
        for (std::size_t i = 0; i < g.size(); ++i) {
            const Vec2f& p = g[(i + g.size() - 1) % g.size()];
            const Vec2f& c = g[i];
            const Vec2f& n = g[(i + 1) % g.size()];
            const float e1x = c.x - p.x, e1y = c.y - p.y, e2x = n.x - c.x, e2y = n.y - c.y;
            if (std::hypot(e1x, e1y) < gutter || std::hypot(e2x, e2y) < gutter) {
                continue;  // edges shorter than the gutter: fold territory
            }
            const float cross = e1x * e2y - e1y * e2x;
            if (cross * sign <= 0.0F) {
                continue;  // reflex or straight
            }
            const float turn = std::atan2(std::fabs(cross), e1x * e2x + e1y * e2y);
            const float ang  = 180.0F - turn * 180.0F / 3.14159265F;
            if (ang >= best.angleDeg) {
                continue;
            }
            // Distance of a loop point to the chart's corner UVs (and which UV is closest).
            const auto nearestTo = [&](const Vec2f& p, Vec2f& at) {
                float d2 = 1e30F;
                for (const Vec2f& q : pts[ci]) {
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
            Vec2f near  = c;
            float dist  = 0.0F;
            for (int k = -3; k <= 3 && clean; ++k) {
                Vec2f at = c;
                const float d = nearestTo(g[static_cast<std::size_t>((static_cast<long>(i) + k + static_cast<long>(g.size())) % static_cast<long>(g.size()))], at);
                clean = d >= 0.9F * gutter && d <= decimate::chart::kDefaultMiterLimit * gutter * 1.05F;
                if (k == 0) {
                    near = at;
                    dist = d;
                }
            }
            // ...and the loop must not cross itself there (a fold also crosses its own edges).
            const auto at = [&](int k) {
                return g[static_cast<std::size_t>((static_cast<long>(i) + k + static_cast<long>(g.size())) % static_cast<long>(g.size()))];
            };
            const auto side = [](const Vec2f& a, const Vec2f& b, const Vec2f& p) {
                return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
            };
            for (int a = -3; a < 3 && clean; ++a) {
                for (int b = a + 2; b < 3 && clean; ++b) {
                    const bool crosses = side(at(a), at(a + 1), at(b)) * side(at(a), at(a + 1), at(b + 1)) < 0.0F &&
                                         side(at(b), at(b + 1), at(a)) * side(at(b), at(b + 1), at(a + 1)) < 0.0F;
                    clean = !crosses;
                }
            }
            // A tip beyond the miter limit is clamped (a truncated spike, not a miter): skip it.
            if (!clean || dist > 0.97F * decimate::chart::kDefaultMiterLimit * gutter) {
                continue;
            }
            best = {static_cast<int>(ci), static_cast<int>(i), ang, near, dist};
        }
    }
    return best;
}

/// <model>_gutter_corner: vector close-up (rendered at high res, not a pixel crop) of the sharpest
/// convex gutterLoop corner: that chart's triangles, its boundary, and the gutter outline.
bool shotGutterCorner(const Programs& prog, const ViewerState& s,
                      const std::vector<decimate::chart::HalfEdgeChart>& charts, ShotTarget& ss, int n,
                      const std::string& path) {
    std::vector<int> facesPerChart(charts.size(), 0);
    for (const int ci : s.chartSet.faceToChart) {
        if (ci >= 0 && static_cast<std::size_t>(ci) < charts.size()) {
            ++facesPerChart[static_cast<std::size_t>(ci)];
        }
    }
    const int maxFaces = *std::max_element(facesPerChart.begin(), facesPerChart.end());
    const int minFaces = std::max(30, maxFaces / 4);
    SharpCorner pick   = sharpestConvexCorner(s, charts, facesPerChart, minFaces);
    if (pick.chart < 0) {
        pick = sharpestConvexCorner(s, charts, facesPerChart, 0);
    }
    if (pick.chart < 0) {
        std::fprintf(stderr, "[shot] no convex gutterLoop corner found; no corner image\n");
        return true;
    }
    const auto& gl = charts[static_cast<std::size_t>(pick.chart)].gutterLoop;
    const Vec2f g  = gl[static_cast<std::size_t>(pick.vertex)];

    // Centre between the gutter vertex and its chart boundary vertex so both are in frame.
    const Vec2f ctr((g.x + pick.nearest.x) * 0.5F, (g.y + pick.nearest.y) * 0.5F);
    const float half = std::max(3.0F * pick.dist, 12.0F / 255.0F);
    std::fprintf(stderr, "[shot] corner: chart %d (%d faces) gutter vertex %d, interior angle %.1f deg, view +-%.4f uv\n",
                 pick.chart, facesPerChart[static_cast<std::size_t>(pick.chart)], pick.vertex,
                 static_cast<double>(pick.angleDeg), static_cast<double>(half));

    constexpr float kSS = 2.0F;
    const float uvPerPx = 2.0F * half / (static_cast<float>(n) * kSS);
    ShotBatch fill, wire, ink;
    addChartFaces(fill, wire, s, 1.5F * kSS * uvPerPx, [&](int ci) { return ci == pick.chart; });

    // Chart boundary = edges owned by exactly one of the chart's faces (mesh vertex pair key).
    std::map<std::pair<int, int>, std::pair<int, std::pair<Vec2f, Vec2f>>> edges;
    const auto& idx = s.chartedMesh.indices;
    for (std::size_t f = 0; f < s.chartSet.faceToChart.size() && 3 * f + 2 < s.uv.size(); ++f) {
        if (s.chartSet.faceToChart[f] != pick.chart) {
            continue;
        }
        for (std::size_t k = 0; k < 3; ++k) {
            const int a = static_cast<int>(idx[3 * f + k]), b = static_cast<int>(idx[3 * f + (k + 1) % 3]);
            auto& e = edges[{std::min(a, b), std::max(a, b)}];
            ++e.first;
            e.second = {s.uv[3 * f + k], s.uv[3 * f + (k + 1) % 3]};
        }
    }
    for (const auto& [key, e] : edges) {
        if (e.first == 1) {
            ink.seg(e.second.first, e.second.second, 3.0F * kSS * uvPerPx, kBoundaryInk);
        }
    }
    ink.loop(gl, 3.0F * kSS * uvPerPx, kGutterInk);
    ink.square(g, 9.0F * kSS * uvPerPx, kGutterInk);

    glBindFramebuffer(GL_FRAMEBUFFER, ss.fbo);
    const int big    = n * static_cast<int>(kSS);
    const mat4 ortho = makeOrtho(ctr.x - half, ctr.x + half, ctr.y - half, ctr.y + half);
    const ui::Image img = renderImage(big, big, kShotBg, [&] {
        drawBatch(prog, ortho, fill);
        drawBatch(prog, ortho, wire);
        drawBatch(prog, ortho, ink);
    });
    return writeShot(path, ui::downscale(img, static_cast<int>(kSS)));
}

}  // namespace

/// Reduce → chartify → bake (as in the smoke path), then write each panel to `opt.outDir`.
// NOLINTNEXTLINE(readability-function-size)
static int runShots(const Programs& prog, ViewerState& s, const ShotOptions& opt,
                    const std::string& meshPath) {
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(opt.outDir, ec);
    const std::string base = (fs::path(opt.outDir) / fs::path(meshPath).stem()).string();

    for (int i = 0; i < 4; ++i) {
        decimateStep(s);
    }
    computeChartification(s);  // also bakes (seam-fill)
    if (!s.charted || !s.baked) {
        std::fprintf(stderr, "[shot] chartify/bake failed\n");
        return 1;
    }

    // 3/4 view, framed so the model's projected extent fills 85% of the frame.
    s.azimuth   = 145.0F;
    s.elevation = 20.0F;
    fitCameraToMesh(s, static_cast<float>(opt.width) / static_cast<float>(opt.height), 0.85F);

    constexpr int kBig = 3;  // seam close-ups are cropped from a 3x render
    ShotTarget t   = makeShotTarget(opt.width, opt.height);
    ShotTarget big = makeShotTarget(opt.width * kBig, opt.height * kBig);
    glBindFramebuffer(GL_FRAMEBUFFER, t.fbo);
    const int sq = std::min(t.w, t.h);
    bool ok = true;
    ok &= shotPanel(base + "_mesh.tga", t.w, t.h, 1.0F,
                    [&] { render3DPanel(prog, s, 0, 0, t.w, t.h, /*useChartColors=*/false); });
    ok &= shotPanel(base + "_charts.tga", t.w, t.h, 1.0F,
                    [&] { render3DPanel(prog, s, 0, 0, t.w, t.h, /*useChartColors=*/true); });
    ok &= shotPanel(base + "_atlas.tga", sq, sq, 0.82F,
                    [&] { renderUVPanel(prog, s, 0, 0, sq, sq); });
    const LitShots litSeam = shotBakePanels(prog, s, t, big, base, "", ok);

    // Before/after: same charts and UVs, plain blur dilation instead of seam-correct fill.
    const decimate::bake::NormalMap seamMap = s.bakedMap;
    s.bakeSeamFillMode = false;
    bakeAndUpload(s);
    const LitShots litDilate = shotBakePanels(prog, s, t, big, base, "_dilate", ok);
    s.bakeSeamFillMode = true;

    // Texture-space close-up where the two gutters differ most (crisp nearest-neighbour).
    constexpr int kWin = 64, kScale = 10;
    const auto [zx, zy] = mostDifferentWindow(seamMap, s.bakedMap, kWin);
    const auto zoomSeam   = cropZoom(seamMap, zx, zy, kWin, kScale);
    const auto zoomDilate = cropZoom(s.bakedMap, zx, zy, kWin, kScale);
    ok &= ui::writeTGA(base + "_normalmap_zoom.tga", zoomSeam);
    ok &= ui::writeTGA(base + "_normalmap_dilate_zoom.tga", zoomDilate);
    ok &= writeShot(base + "_seamfill_vs_dilate.tga",
                    ui::sideBySide(toImage(zoomSeam), toImage(zoomDilate),
                                   "seam-correct fill: real neighbour", "plain dilation"));
    std::fprintf(stderr, "[shot] zoom window %dx%d at texel (%d,%d), %dx\n", kWin, kWin, zx, zy, kScale);

    // Model-space question: is the difference visible on the lit 3D model at all?
    constexpr int kThresh = 12;  // summed |dR|+|dG|+|dB| out of 765
    const auto st = ui::diffStats(litSeam.frame, litDilate.frame, kThresh);
    std::fprintf(stderr, "[shot] lit model, seam-fill vs dilate: %ld of %d px differ by >= %d/765 (max %d)\n",
                 st.pixelsOver, litSeam.frame.w * litSeam.frame.h, kThresh, st.maxDiff);
    const auto bst = ui::diffStats(litSeam.big, litDilate.big, kThresh);
    if (bst.pixelsOver == 0) {
        std::fprintf(stderr, "[shot] no visible seam difference at model level; no seam close-up written\n");
    } else {
        constexpr int kCrop = 640;
        const auto [cx, cy] = ui::mostDifferentWindow(litSeam.big, litDilate.big, kCrop, 16, kThresh);
        const auto cropSeam   = ui::cropImage(litSeam.big, cx, cy, kCrop, kCrop);
        const auto cropDilate = ui::cropImage(litDilate.big, cx, cy, kCrop, kCrop);
        ok &= writeShot(base + "_baked_lit_seam_zoom.tga", cropSeam);
        ok &= writeShot(base + "_baked_lit_dilate_zoom.tga", cropDilate);
        ok &= writeShot(base + "_baked_lit_seamfill_vs_dilate.tga",
                        ui::sideBySide(cropSeam, cropDilate, "seam-correct fill: real neighbour", "plain dilation"));
        std::fprintf(stderr, "[shot] 3D seam crop %dx%d at (%d,%d) of the %dx%d render\n", kCrop, kCrop, cx, cy,
                     litSeam.big.w, litSeam.big.h);
    }

    // Gutter views: the packer's real gutterLoop (core), drawn as vector geometry, 2x supersampled.
    const auto packed = packedCharts(s);
    if (!packed.empty()) {
        constexpr int kCornerPx = 1000;
        ShotTarget ss = makeShotTarget(2 * std::max(sq, kCornerPx), 2 * std::max(sq, kCornerPx));
        ok &= shotAtlasGutter(prog, s, packed, ss, sq, base + "_atlas_gutter.tga");
        ok &= shotGutterCorner(prog, s, packed, ss, kCornerPx, base + "_gutter_corner.tga");
        destroyShotTarget(ss);
    } else {
        ok = false;
    }

    destroyShotTarget(big);
    destroyShotTarget(t);
    return ok ? 0 : 1;
}

// ---------------------------------------------------------------------------
// runViewer — public entry point
// ---------------------------------------------------------------------------

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
int runViewer(const std::string& meshPath, int smokeFrames, const ShotOptions* shot) {
    decimate::MeshData original;
    try {
        original = decimate::load(meshPath);
    } catch (const std::exception& ex) {
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

    const Aabb box    = computeAabb(original.positions);
    const float radius = aabbRadius(box);

    ViewerState state;
    state.original            = original;
    state.current             = original;
    state.originalVertexCount = static_cast<int>(original.positions.size());
    state.currentTarget       = state.originalVertexCount;
    state.centreX             = (box.lo.x + box.hi.x) * 0.5F;
    state.centreY             = (box.lo.y + box.hi.y) * 0.5F;
    state.centreZ             = (box.lo.z + box.hi.z) * 0.5F;
    state.distance            = radius * 2.5F;

    std::fprintf(stderr, "[load] '%s'  verts=%d  tris=%d  radius=%.3f\n",
        meshPath.c_str(),
        state.originalVertexCount,
        static_cast<int>(original.indices.size() / 3),
        static_cast<double>(radius));

    gState = &state;

    GlfwWindow win;

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

    glfwSetKeyCallback(win.handle,         onKey);
    glfwSetMouseButtonCallback(win.handle, onMouseButton);
    glfwSetCursorPosCallback(win.handle,   onCursorPos);
    glfwSetScrollCallback(win.handle,      onScroll);

    // Log the GL version so the user can confirm core-profile 3.3 landed.
    std::fprintf(stderr, "[GL] version: %s\n",
                 // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
                 reinterpret_cast<const char*>(glGetString(GL_VERSION)));

    // Build shader programs — must happen after the GL context is current.
    const Programs prog = buildPrograms();

    glEnable(GL_DEPTH_TEST);

    if (smokeFrames > 0) {
        std::fprintf(stderr, "[smoke] step 0: verts=%d  tris=%d\n",
            static_cast<int>(state.current.positions.size()),
            static_cast<int>(state.current.indices.size()   / 3));
        for (int i = 0; i < 4; ++i) {
            decimateStep(state);
        }
        computeChartification(state);

        const int triCount  = static_cast<int>(state.current.indices.size()   / 3);
        const int numCharts = static_cast<int>(state.chartSet.charts.size());
        const int numUV     = static_cast<int>(state.uv.size());

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
        if (height == 0) { height = 1; }
        if (width  == 0) { width  = 1; }

        const int hw = width  / 2;
        const int hh = height / 2;

        // Panel ① top-left — dark blue-grey
        glViewport(0, hh, hw, hh);
        glScissor(0, hh, hw, hh);
        glEnable(GL_SCISSOR_TEST);
        glClearColor(0.10F, 0.11F, 0.14F, 1.0F);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glDisable(GL_SCISSOR_TEST);
        render3DPanel(prog, state, 0, hh, hw, hh, /*useChartColors=*/false);

        // Panel ② top-right — dark purple-grey
        glViewport(hw, hh, hw, hh);
        glScissor(hw, hh, hw, hh);
        glEnable(GL_SCISSOR_TEST);
        glClearColor(0.10F, 0.09F, 0.15F, 1.0F);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glDisable(GL_SCISSOR_TEST);
        render3DPanel(prog, state, hw, hh, hw, hh, /*useChartColors=*/true);

        // Panel ③ bottom-left — very dark green
        glViewport(0, 0, hw, hh);
        glScissor(0, 0, hw, hh);
        glEnable(GL_SCISSOR_TEST);
        glClearColor(0.08F, 0.11F, 0.09F, 1.0F);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glDisable(GL_SCISSOR_TEST);
        renderUVPanel(prog, state, 0, 0, hw, hh);

        // Panel ④ bottom-right — dark amber
        glViewport(hw, 0, hw, hh);
        glScissor(hw, 0, hw, hh);
        glEnable(GL_SCISSOR_TEST);
        glClearColor(0.13F, 0.10F, 0.06F, 1.0F);
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

    // Clean up programs — context is still current here.
    glDeleteProgram(prog.mesh);
    glDeleteProgram(prog.ortho);
    glDeleteProgram(prog.quad);
    glDeleteProgram(prog.baked);

    if (smokeFrames > 0) {
        std::fprintf(stderr, "[smoke] done — %d frames rendered\n", frameCount);
    }

    return exitCode;
}

}  // namespace ui
