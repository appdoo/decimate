# Model directory for tests/tools. Overridable so a clone without the large models can be
# exercised against a directory holding only spot.obj (model tests must then SKIP, not fail).
set(MESH_MODEL_DIR "${CMAKE_SOURCE_DIR}/model" CACHE PATH
    "Directory holding test meshes (spot.obj always; cow/horse/bunny optional)")

# Absolute per-component tolerance for the spot.obj UV golden (B1 + uv_golden_reproducible).
# Measured: -O0 vs -O2 (Release) max |dUV| = 1.19e-7 on spot.obj (clang/arm64). eps leaves
# headroom for other compilers' FMA contraction while still catching a 1e-3 drift 100x over.
set(MESH_UV_GOLDEN_EPS 1e-5)

# RELEASE-PINNED tests: cases tagged [release-pinned] (bijective B2 energy sums, chart_injectivity I2
# table, seam_aware_reduce I1/I2 LOD counts) and the ctests uv_golden_reproducible[_lscm] pin exact
# numbers that are float-order sensitive across optimisation levels.  They are pinned in build-rel
# (Release); in a Debug build they may legitimately differ.  Run the rest with
# `mesh_tests "~[release-pinned]"`.
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
    # Spec 24-03: hidden state and math (M1 OBJ strtok/locale, M2 lexer, M3 SMat pivots, M4 Quat)
    correctness_io_math_test.cpp
    # Spec 07: QEM edge-collapse reducer characterization
    reduce_characterization.cpp
    # Octree vertex-clustering decimator (heap-free strategy) — C1-C6
    cluster_reduce_test.cpp
    # Half-edge edge-collapse decimator (modernized 2003 CHalfedgeCollapse) — HE1-HE2
    halfedge_reduce_test.cpp
    # Spec 08-01: Region-growing chart segmentation
    chart_segmentation_test.cpp
    # Spec 08-02: LSCM conformal mapper + atlas packer + CG solver
    parametrize_pack_test.cpp
    # Spec 08-02: UV-atlas quality on the real cow (seam correctness + fold-overs)
    uv_quality_test.cpp
    # Spec 08-03: Chartification-fidelity integration tests on .iv models (D1/D2/D3)
    chart_iv_integration_test.cpp
    # Spec 09-01: HalfEdgeChart + boundary-driven footprint (B1-B2, S1-S5)
    halfedge_chart_test.cpp
    # Spec 09-02 (+amendment): miterOffsetLoop + placed gutterLoop (S1-S3, S5', S6, B1, B2, B3')
    miter_gutter_test.cpp
    # Spec 19 Amendment 4: conservative silhouette primitive coverSegment (S1-S4)
    cover_segment_test.cpp
    # Spec 19-03: chart injectivity guard (selfOverlapPairs, splitChart, flattenCharts)
    chart_injectivity_test.cpp
    # Spec 20: error-bounded chart count (stretch splits) + [.chartcost] measurement
    chart_stretch_test.cpp
    # Spec 25: chartify + bakeAtlas facade == the viewer's hand-assembled pipeline
    pipeline_facade_test.cpp
    # Spec 19-02: self-seam cut (seamHalfEdges export, corner-built flatten) - B1-B3, S1-S4
    self_seam_cut_test.cpp
    # Bijective free-boundary flattening (Smith & Schaefer 2015) as the guard's repair - S1-S5, B1, I1-I2
    bijective_flatten_test.cpp
    # Spec 10: raytraced normal bake characterization (B1-B4, S1, S2)
    bake_characterization.cpp
    # Spec 22-01: exact BVH nearest-hit + deterministic parallel bake (I1, I3 cow, I4, B1, S1-S7, [.bench])
    bake_bvh_test.cpp
    # Seam-correct gutter fill (bakeSeamFill) — 2-pass cross-chart sampling
    seam_bake_test.cpp
    # Spec 17-01: seam energy metric + seam erasure (goldens from the MIT reference)
    seam_erase_test.cpp
    # Spec 17-02: seam straightening (section 6) + cutter threshold parameter
    seam_straighten_test.cpp
    # Spec 17-03: seam-aware decimation (section 5)
    seam_aware_reduce_test.cpp
    # Spec 24-02: NaN / non-finite / out-of-range input paths (N1-N5)
    correctness_nan_test.cpp
    # Spec 11: TGA writer round-trip (B1) — tga_writer.cpp has no GL dep, safe here
    tga_writer_test.cpp
    # Spec 24-01: chart hangs and stack depth (H1-H5)
    correctness_chart_test.cpp
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
    "${CMAKE_SOURCE_DIR}/src/seam"   # Spec 17-01: erase_domain.h / seam_pieces.h (white-box checks)
    "${CMAKE_SOURCE_DIR}/src/bake"   # Spec 22-01: bake_detail.h (bakeWith test seam)
    "${CMAKE_SOURCE_DIR}/src/chart"  # bijective_detail.h (gradient / flip-bound white-box checks)
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
# chart_legacy_golden — Spec 19-01: live A/B of segment() against the UNMODIFIED 2003
# CChartification, compiled at build time (legacy2003/ is never edited).
#   legacy_chart_cut  C++14 OBJECT lib: legacy2003/chartification/*.cpp + math/{halfedge,vec}.cpp
#                     + tests/legacy_cut/legacy_side.cpp (timer / mapper / packer stubs).
#   -fms-extensions resolves the backslash includes, -fdelayed-template-parsing the lazy template
#   in cset.h:65; tests/shim/legacy_cut/windows.h is an empty <windows.h>.
# Like capture_golden / capture_bake_golden, these are the only targets that reach legacy2003/.
# ---------------------------------------------------------------------------
add_library(legacy_chart_cut OBJECT
    "${CMAKE_SOURCE_DIR}/legacy2003/chartification/cchartification.cpp"
    "${CMAKE_SOURCE_DIR}/legacy2003/chartification/cchart.cpp"
    "${CMAKE_SOURCE_DIR}/legacy2003/chartification/cchartset.cpp"
    "${CMAKE_SOURCE_DIR}/legacy2003/chartification/cchartboundary.cpp"
    "${CMAKE_SOURCE_DIR}/legacy2003/math/halfedge.cpp"
    "${CMAKE_SOURCE_DIR}/legacy2003/math/vec.cpp"
    legacy_cut/legacy_side.cpp
)
target_include_directories(legacy_chart_cut PRIVATE
    "${CMAKE_SOURCE_DIR}/tests/shim/legacy_cut"   # empty windows.h FIRST
    "${CMAKE_SOURCE_DIR}/legacy2003"
    "${CMAKE_SOURCE_DIR}/legacy2003/math"
)
set_target_properties(legacy_chart_cut PROPERTIES CXX_STANDARD 14 CXX_EXTENSIONS OFF)
target_compile_options(legacy_chart_cut PRIVATE -fms-extensions -fdelayed-template-parsing -w)

add_executable(chart_legacy_golden
    legacy_cut/chart_legacy_golden.cpp
    $<TARGET_OBJECTS:legacy_chart_cut>
)
target_link_libraries(chart_legacy_golden PRIVATE decimate decimate::math)
target_compile_definitions(chart_legacy_golden PRIVATE MESH_MODEL_DIR="${MESH_MODEL_DIR}")
add_test(NAME chart_legacy_golden COMMAND chart_legacy_golden)
set_tests_properties(chart_legacy_golden PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 600)

# ---------------------------------------------------------------------------
# capture_uv_golden — Spec 09-01 Task 1: spot.obj per-corner UV golden.
# Loads spot.obj, runs segment -> parametrizeCorners, and either writes
# tests/fixtures/spot_uv_golden.txt (regenerate; use the Release build, build-rel/ --
# the bijective flattener is not bit-stable across optimisation levels) or, with
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

# Spec 19 Amendment 3: Flattener::Lscm reproduces the pre-amendment spot atlas.
add_test(
    NAME uv_golden_reproducible_lscm
    COMMAND capture_uv_golden --lscm --check
        ${CMAKE_SOURCE_DIR}/tests/fixtures/spot_uv_golden_lscm.txt
        ${MESH_UV_GOLDEN_EPS}
)

# ---------------------------------------------------------------------------
# atlas_invariant_test — Spec 19-02 Amendment 2: no texel is covered by two charts' UV triangles
# on cow/spot/bunny/horse x caps 16..128 (full pipeline, 1024^2 texel grid).  Heavy in the
# unoptimised build, hence the `slow` label and the long timeout; it runs in both build dirs.
# ---------------------------------------------------------------------------
add_executable(atlas_invariant_test atlas_invariant_test.cpp)
target_link_libraries(atlas_invariant_test PRIVATE Catch2::Catch2WithMain decimate decimate::math)
target_compile_definitions(atlas_invariant_test PRIVATE MESH_MODEL_DIR="${MESH_MODEL_DIR}")
add_test(NAME atlas_invariant_test COMMAND atlas_invariant_test)
set_tests_properties(atlas_invariant_test PROPERTIES LABELS slow TIMEOUT 3600)

# ---------------------------------------------------------------------------
# bake_bvh_equality_slow - Spec 22-01 I2 + I3 on bunny/horse (brute-force oracle: minutes).
# ---------------------------------------------------------------------------
add_executable(bake_bvh_equality_slow bake_bvh_equality_slow.cpp)
target_link_libraries(bake_bvh_equality_slow PRIVATE Catch2::Catch2WithMain decimate decimate::math)
target_include_directories(bake_bvh_equality_slow PRIVATE "${CMAKE_SOURCE_DIR}/src/bake")
target_compile_definitions(bake_bvh_equality_slow PRIVATE MESH_MODEL_DIR="${MESH_MODEL_DIR}")
add_test(NAME bake_bvh_equality_slow COMMAND bake_bvh_equality_slow)
set_tests_properties(bake_bvh_equality_slow PROPERTIES LABELS slow TIMEOUT 3600)

# ---------------------------------------------------------------------------
# seamless_sweep — Spec 17-03: chartification-vs-seam-cost sweep (NOT a test).  The
# 4 x 3 setting run is for the Release build (build-rel); ctest runs one setting on spot (I3).
# ---------------------------------------------------------------------------
add_executable(seamless_sweep seamless_sweep.cpp)
target_link_libraries(seamless_sweep PRIVATE decimate decimate::math)
target_compile_definitions(seamless_sweep PRIVATE MESH_MODEL_DIR="${MESH_MODEL_DIR}")

add_test(
    NAME seamless_sweep_csv
    COMMAND ${CMAKE_COMMAND}
        -DSWEEP_EXE=$<TARGET_FILE:seamless_sweep>
        -DMESH=${CMAKE_SOURCE_DIR}/model/spot.obj
        -DOUT=${CMAKE_BINARY_DIR}/sweep.csv
        -P ${CMAKE_SOURCE_DIR}/tests/check_sweep_csv.cmake
)
set_tests_properties(seamless_sweep_csv PROPERTIES DEPENDS seamless_sweep)
