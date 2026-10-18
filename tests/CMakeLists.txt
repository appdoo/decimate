# Model directory for tests/tools. Overridable so a clone without the large models can be
# exercised against a directory holding only spot.obj (model tests must then SKIP, not fail).
set(MESH_MODEL_DIR "${CMAKE_SOURCE_DIR}/model" CACHE PATH
    "Directory holding test meshes (spot.obj always; cow/horse/bunny optional)")

# Absolute per-component tolerance for the spot.obj UV golden (B1 + uv_golden_reproducible).
# Measured: -O0 vs -O2 (Release) max |dUV| = 1.19e-7 on spot.obj (clang/arm64). eps leaves
# headroom for other compilers' FMA contraction while still catching a 1e-3 drift 100x over.
set(MESH_UV_GOLDEN_EPS 1e-5)

add_executable(mesh_tests
    main.cpp
    puremath_vec3_test.cpp
    diag_stacktrace_test.cpp
    # Spec 05: half-edge boundary-walk API (firstBoundaryEdge/nextBoundaryEdge/boundaryPointIndices)
    half_edge_boundary_test.cpp
    raytracing_characterization.cpp
    core_math_characterization.cpp
    # Spec 06: mesh loader characterization (OBJ + Open Inventor)
    loader_characterization.cpp
    # Spec 07: QEM edge-collapse reducer characterization
    reduce_characterization.cpp
    # Octree vertex-clustering decimator (heap-free strategy) — C1-C6
    cluster_reduce_test.cpp
    # Half-edge edge-collapse decimator (modernized 2003 CHalfedgeCollapse) — HE1-HE2
    halfedge_reduce_test.cpp
    # Spec 08-01: Region-growing chart segmentation
    chart_segmentation_test.cpp
    # Spec 08-02: LSCM conformal mapper + atlas packer + CG solver
    cg_solver_test.cpp
    parametrize_pack_test.cpp
    # Spec 08-02: UV-atlas quality on the real cow (seam correctness + fold-overs)
    uv_quality_test.cpp
    # Spec 08-03: Chartification-fidelity integration tests on .iv models (D1/D2/D3)
    chart_iv_integration_test.cpp
    # Spec 09-01: HalfEdgeChart + boundary-driven footprint (B1-B2, S1-S5)
    halfedge_chart_test.cpp
    # Spec 09-02 (+amendment): miterOffsetLoop + placed gutterLoop (S1-S3, S5', S6, B1, B2, B3')
    miter_gutter_test.cpp
    # Spec 10: raytraced normal bake characterization (B1-B4, S1, S2)
    bake_characterization.cpp
    # Seam-correct gutter fill (bakeSeamFill) — 2-pass cross-chart sampling
    seam_bake_test.cpp
    # Spec 11: TGA writer round-trip (B1) — tga_writer.cpp has no GL dep, safe here
    tga_writer_test.cpp
    "${CMAKE_SOURCE_DIR}/ui/tga_writer.cpp"
)

# Link Catch2 (we supply our own main() in main.cpp so the executable can
# install the diag crash handler before any test runs — Spec 03-03).
# diag is linked here (executable boundary) — never inside decimate itself.
target_link_libraries(mesh_tests PRIVATE
    Catch2::Catch2
    decimate
    decimate::math
    diag
)

# Spec 11: tga_writer.h lives in ui/ — expose it for the TGA test TU.
# Only the header is needed; the .cpp source is compiled directly above.
target_include_directories(mesh_tests PRIVATE
    "${CMAKE_SOURCE_DIR}/ui"
)

# Spec 06: provide fixture + model paths to loader_characterization.cpp.
# MESH_FIXTURES_DIR — absolute path to tests/fixtures/
# MESH_MODEL_DIR    — absolute path to model/ (shared fixtures at repo root)
target_compile_definitions(mesh_tests PRIVATE
    MESH_FIXTURES_DIR="${CMAKE_SOURCE_DIR}/tests/fixtures"
    MESH_MODEL_DIR="${MESH_MODEL_DIR}"
    MESH_UV_GOLDEN_EPS=${MESH_UV_GOLDEN_EPS}
)

# Register tests with CTest via Catch2's CMake integration
catch_discover_tests(mesh_tests)

# ---------------------------------------------------------------------------
# atlas_dump — diagnostic: packed-atlas silhouette PPM + density (NOT a test).
# Used to compare current packing against the 2003 thesis atlas figures.
# ---------------------------------------------------------------------------
add_executable(atlas_dump atlas_dump.cpp)
target_link_libraries(atlas_dump PRIVATE decimate decimate::math)
target_compile_definitions(atlas_dump PRIVATE
    MESH_MODEL_DIR="${MESH_MODEL_DIR}")
target_include_directories(atlas_dump PRIVATE "${CMAKE_SOURCE_DIR}/src/reduce")

# raster_bench — diagnostic: real packAtlas() timing on cow charts (calls rasterLine).
add_executable(raster_bench raster_bench.cpp)
target_link_libraries(raster_bench PRIVATE decimate decimate::math)
target_compile_definitions(raster_bench PRIVATE
    MESH_MODEL_DIR="${MESH_MODEL_DIR}"
    MESH_FIXTURES_DIR="${CMAKE_SOURCE_DIR}/tests/fixtures")
target_include_directories(raster_bench PRIVATE "${CMAKE_SOURCE_DIR}/src/reduce")

# reduce_chartify_timing — diagnostic: reduce→chartify timing vs target res.
add_executable(reduce_chartify_timing reduce_chartify_timing.cpp)
target_link_libraries(reduce_chartify_timing PRIVATE decimate decimate::math)
target_compile_definitions(reduce_chartify_timing PRIVATE
    MESH_MODEL_DIR="${MESH_MODEL_DIR}")

# reduce_bench — speed benchmark for decimate::reduce (QEM). best-of-N throughput,
# sampler-friendly for the per-phase profile (heap vs quadric vs adjacency-cache).
add_executable(reduce_bench reduce_bench.cpp)
target_link_libraries(reduce_bench PRIVATE decimate decimate::math)
target_compile_definitions(reduce_bench PRIVATE
    MESH_MODEL_DIR="${MESH_MODEL_DIR}")

# mesh_defects — topology defect detector (non-manifold/pinch/spike). Diagnostic.
add_executable(mesh_defects mesh_defects.cpp)
target_link_libraries(mesh_defects PRIVATE decimate decimate::math)

# ---------------------------------------------------------------------------
# capture_golden — legacy rtrace compiled on clang, emits the golden fixture.
# NOT part of decimate. Links only the legacy raytracing cone (pure functions).
# ---------------------------------------------------------------------------
add_executable(capture_golden
    capture_golden.cpp
    "${CMAKE_SOURCE_DIR}/legacy2003/raytracing/raytracing.cpp"
    "${CMAKE_SOURCE_DIR}/legacy2003/math/vec.cpp"
)

# Legacy headers require include paths for the legacy tree.
# This target deliberately reaches into legacy2003/; it is the ONLY target allowed to.
target_include_directories(capture_golden PRIVATE
    "${CMAKE_SOURCE_DIR}/legacy2003/math"
    "${CMAKE_SOURCE_DIR}/legacy2003"
)

# C++14 is sufficient for the legacy code; avoid forcing C++20 features onto
# the unmodified legacy sources (they use register, pre-C++17 struct init, etc.)
set_target_properties(capture_golden PROPERTIES CXX_STANDARD 14)

# Suppress warnings from legacy code — we only care about compile success.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(capture_golden PRIVATE
        -Wno-deprecated-register
        -Wno-unused-variable
        -Wno-sign-compare
        -Wno-unused-parameter
    )
endif()

# CTest fixture: re-run capture to verify reproducibility (T5).
# The primary golden is the checked-in file; this test re-generates and diffs.
add_test(
    NAME golden_reproducible
    COMMAND ${CMAKE_COMMAND}
        -DCAPTURE_EXE=$<TARGET_FILE:capture_golden>
        -DGOLDEN_FILE=${CMAKE_SOURCE_DIR}/tests/fixtures/raytracing_golden.txt
        -DBUILD_DIR=${CMAKE_BINARY_DIR}
        -P ${CMAKE_SOURCE_DIR}/tests/check_golden_reproducible.cmake
)
set_tests_properties(golden_reproducible PROPERTIES
    DEPENDS capture_golden
)

# ---------------------------------------------------------------------------
# capture_bake_golden — Spec 10 Phase A golden capture.
# Links: unmodified legacy texturesample.cpp + raytracing.cpp + vec.cpp
#        PLUS decimate (for reduce/segment/parametrizeCorners/loadObj).
# The shim/ dir (utility.h with no windows.h) is prepended to the include
# path so it intercepts legacy2003/utility.h before the legacy tree is searched.
# ---------------------------------------------------------------------------
add_executable(capture_bake_golden
    capture_bake_golden.cpp
    "${CMAKE_SOURCE_DIR}/legacy2003/raytracing/raytracing.cpp"
    "${CMAKE_SOURCE_DIR}/legacy2003/math/vec.cpp"
)

# Shim dir FIRST so tests/shim/utility.h overrides legacy2003/utility.h.
# Then legacy include dirs so the legacy code compiles unmodified.
target_include_directories(capture_bake_golden PRIVATE
    "${CMAKE_SOURCE_DIR}/tests/shim"          # bake_capture_pre.h; force-included
    "${CMAKE_SOURCE_DIR}/legacy2003/math"         # math3D.h, vec.h
    "${CMAKE_SOURCE_DIR}/legacy2003/raytracing"   # raytracing.h
    "${CMAKE_SOURCE_DIR}/legacy2003"              # def.h
    # NOTE: legacy2003/split is intentionally OMITTED.
    # texturesample.h and chart.h in legacy2003/split/ are found via the relative
    # #include in texturesample.cpp (GCC/Clang search the including file's dir
    # first for quoted includes), so they don't need to be on the -I path.
    # Adding legacy2003/split would cause the legacy chart.h (which pulls sparsematrix.h
    # → windows.h) to shadow the core chart.h when capture_bake_golden.cpp
    # includes "chart.h".
)

# Core public headers (reduce, chart, mesh, io, etc.)
target_link_libraries(capture_bake_golden PRIVATE decimate decimate::math)

# C++14 — legacy code uses register, pre-C++17 struct init, etc.
set_target_properties(capture_bake_golden PROPERTIES CXX_STANDARD 14)

# Suppress legacy warnings (we only care about compile success + correctness).
# -include forces the shim pre-header before any TU source — this is the only
# portable way to block windows.h and the ill-formed CSet::operator[] without
# editing legacy2003/ files.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(capture_bake_golden PRIVATE
        -Wno-deprecated-register
        -Wno-register            # clang 10+: 'register' in C++17+
        -Wno-unused-variable
        -Wno-sign-compare
        -Wno-unused-parameter
        -Wno-reorder-ctor
        -include ${CMAKE_SOURCE_DIR}/tests/shim/bake_capture_pre.h
    )
endif()

# Supply model + fixture path macros
target_compile_definitions(capture_bake_golden PRIVATE
    MESH_MODEL_DIR="${MESH_MODEL_DIR}"
    MESH_FIXTURES_DIR="${CMAKE_SOURCE_DIR}/tests/fixtures"
)

# NOTE: the golden_bake_reproducible CTest was REMOVED — it pinned a SHA-256 of
# the whole reduce→chartify→pack→bake chain, so any upstream improvement (the
# fused Tetris packer changed the atlas UVs) re-broke it despite the bake KERNEL
# being unchanged — the same brittle whole-pipeline pin we dropped for the two
# bake buffer-hash Catch2 cases.  The `capture_bake_golden` target is KEPT: it
# still writes `bake_golden.txt`, whose `legacy_vertex_normal` anchors feed the
# surviving normals-equivalence test (computeVertexNormals vs legacy createNormals,
# which is UV-independent and thus stable).  Bake correctness is covered by the
# kernel tests (normals-equivalence, miss-fallback, coverage), not a chain hash.

# ---------------------------------------------------------------------------
# capture_uv_golden — Spec 09-01 Task 1: spot.obj per-corner UV golden.
# Loads spot.obj, runs segment -> parametrizeCorners, and either writes
# tests/fixtures/spot_uv_golden.txt (regenerate; use a Debug build) or, with
# --check <golden> <eps>, compares a fresh run against it by tolerance.
# NOT part of decimate. Links only decimate (no legacy reach).
# ---------------------------------------------------------------------------
add_executable(capture_uv_golden
    capture_uv_golden.cpp
)

target_link_libraries(capture_uv_golden PRIVATE decimate decimate::math)

target_compile_definitions(capture_uv_golden PRIVATE
    MESH_MODEL_DIR="${MESH_MODEL_DIR}"
    MESH_FIXTURES_DIR="${CMAKE_SOURCE_DIR}/tests/fixtures"
)

# CTest: a fresh capture must match the checked-in golden within MESH_UV_GOLDEN_EPS
# (S6 / AC5). Tolerance, not bytes: Release (-O2) reorders float math vs Debug.
add_test(
    NAME uv_golden_reproducible
    COMMAND capture_uv_golden --check
        ${CMAKE_SOURCE_DIR}/tests/fixtures/spot_uv_golden.txt
        ${MESH_UV_GOLDEN_EPS}
)
