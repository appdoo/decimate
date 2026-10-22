# ---------------------------------------------------------------------------
# decimate — modular component libraries.
#
# Each pipeline module is its own target with a namespaced ALIAS
# `decimate::<module>`, so a consumer can pull exactly what it needs:
#   need a half-edge mesh?  -> decimate::mesh
#   need the atlas packer?  -> decimate::chart
#   need normal baking?     -> decimate::bake
#
# Dependency DAG (acyclic — verified):
#   math (leaf: vec math + numeric)   types (-> math)   raytracing (leaf)
#   mesh   -> types, math
#   io     -> types, math
#   chart  -> mesh, math, types     (Eigen: private, conformal_mapper.cpp + bijective_mapper.cpp only)
#   reduce -> mesh, math, types     (Eigen: private, seam_aware_reduce.cpp only)
#   bake   -> mesh, raytracing, types, math
#   seam   -> bake, types, math     (Eigen: private, .cpp files only)
#   pipeline -> chart, bake, seam   (chartify + bakeAtlas)
#
# `decimate` (alias `decimate::all`) is an umbrella that re-exports every
# component, so a consumer wanting the whole pipeline links a single target.
# It carries the version smoke-test TU.
# ---------------------------------------------------------------------------

# --- decimate::math — foundation: vector math + numeric utilities -----------
# Vec2/3/4, Mat4, Quat, AABB, Plane (absorbed from the former pure-math) +
# sparse_matrix / sparse_bit_set / associative_heap / mat4-invert.
# All namespace decimate::math. The single numeric foundation.
add_library(decimate_math STATIC
    math/sparse_bit_set.cpp
)
add_library(decimate::math ALIAS decimate_math)
target_include_directories(decimate_math PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_compile_features(decimate_math PUBLIC cxx_std_20)

# --- decimate::types — value types MeshData / Face (header-only) ------------
add_library(decimate_types INTERFACE)
add_library(decimate::types ALIAS decimate_types)
target_include_directories(decimate_types INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(decimate_types INTERFACE decimate::math)

# --- decimate::mesh — half-edge mesh ---------------------------------------
add_library(decimate_mesh STATIC mesh/half_edge_mesh.cpp)
add_library(decimate::mesh ALIAS decimate_mesh)
target_include_directories(decimate_mesh PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(decimate_mesh PUBLIC decimate::types decimate::math)

# --- decimate::io — OBJ + Inventor loaders + lexer -------------------------
add_library(decimate_io STATIC
    io/obj_loader.cpp
    io/inventor_loader.cpp
    io/lexer.cpp
    io/load.cpp
)
add_library(decimate::io ALIAS decimate_io)

# fast_float v8.3.2 parses the loaders' floating-point text.  std::from_chars for
# float/double is missing from Apple libc++, and the C replacements are unusable here:
# atof/strtod read the process-global LC_NUMERIC locale (a "de_DE" comma-decimal
# process misreads "1.5") and strtok keeps hidden static state (not thread-safe).
# fast_float is locale-independent, thread-safe and correctly rounded.  Header-only and
# populated only (its own CMake is never run), added as a SYSTEM PRIVATE include like
# Eigen below: no public header names it and the export set needs no extra target.
# Override the download with -DFETCHCONTENT_SOURCE_DIR_FAST_FLOAT=<path>.
include(FetchContent)
FetchContent_Declare(
    fast_float
    GIT_REPOSITORY https://github.com/fastfloat/fast_float.git
    GIT_TAG        v8.3.2
    GIT_SHALLOW    TRUE
    SOURCE_SUBDIR  decimate-populate-only   # header-only: populate, don't build
)
FetchContent_MakeAvailable(fast_float)

target_include_directories(decimate_io PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_include_directories(decimate_io SYSTEM PRIVATE $<BUILD_INTERFACE:${fast_float_SOURCE_DIR}/include>)
target_link_libraries(decimate_io PUBLIC decimate::types decimate::math)

# --- decimate::raytracing — pure ray functions -----------------------------
add_library(decimate_raytracing STATIC raytracing/raytracing.cpp raytracing/bvh.cpp)
add_library(decimate::raytracing ALIAS decimate_raytracing)
target_include_directories(decimate_raytracing PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(decimate_raytracing PUBLIC decimate::types decimate::math)

# --- decimate::chart — segmentation + LSCM + atlas packing -----------------
add_library(decimate_chart STATIC
    chart/segmenter.cpp
    chart/conformal_mapper.cpp
    chart/atlas_packer.cpp
    chart/injectivity.cpp
    chart/chart_cut.cpp
    chart/bijective_mapper.cpp
)
add_library(decimate::chart ALIAS decimate_chart)
target_include_directories(decimate_chart PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(decimate_chart PUBLIC decimate::mesh decimate::math decimate::types decimate::math)

# --- decimate::reduce — QEM edge-collapse reducer --------------------------
add_library(decimate_reduce STATIC
    reduce/qerror.cpp
    reduce/quadric_error_system.cpp
    reduce/reducer.cpp
    reduce/cluster_reduce.cpp
    reduce/corner_mesh.cpp
    reduce/halfedge_reduce.cpp
    reduce/seam_aware_reduce.cpp
)
add_library(decimate::reduce ALIAS decimate_reduce)
target_include_directories(decimate_reduce PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(decimate_reduce PUBLIC decimate::mesh decimate::math decimate::types decimate::math)

# --- decimate::bake — raytraced normal bake --------------------------------
add_library(decimate_bake STATIC
    bake/normals.cpp
    bake/normal_bake.cpp
)
add_library(decimate::bake ALIAS decimate_bake)
target_include_directories(decimate_bake PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(decimate_bake PUBLIC decimate::mesh decimate::raytracing decimate::types decimate::math)
# Parallel bake (Spec 22-01): std::thread needs the platform thread library on Linux; OFF compiles
# the same code with one worker.
if(DECIMATE_THREADS)
    find_package(Threads REQUIRED)
    target_link_libraries(decimate_bake PUBLIC Threads::Threads)
else()
    target_compile_definitions(decimate_bake PRIVATE DECIMATE_NO_THREADS)
endif()

# --- decimate::seam — seam metric + seam erasure (Liu et al. 2017) ---------
# Eigen 3.4.0 is populated only (its own CMake is never run) and added as a SYSTEM
# PRIVATE include of the .cpp files: no public header names Eigen, and the install /
# export set needs no Eigen target. This is the first FetchContent dependency of the
# core, so the first core configure needs the network once.
include(FetchContent)
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE
    SOURCE_SUBDIR  decimate-populate-only   # no CMakeLists.txt there: populate, don't build
)
FetchContent_MakeAvailable(eigen)

add_library(decimate_seam STATIC
    seam/weld.cpp
    seam/find_seams.cpp
    seam/seam_pieces.cpp
    seam/seam_energy.cpp
    seam/erase_domain.cpp
    seam/erase_system.cpp
    seam/erase.cpp
    seam/straighten.cpp
    seam/straighten_atlas.cpp
    seam/straighten_chains.cpp
    seam/straighten_warp.cpp
)
add_library(decimate::seam ALIAS decimate_seam)
# The bijective flattener's Tutte solve (chart/bijective_mapper.cpp) uses the same private Eigen.
target_include_directories(decimate_chart SYSTEM PRIVATE $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>)
target_compile_definitions(decimate_chart PRIVATE EIGEN_MPL2_ONLY)
# Spec 17-03: the seam-aware reducer's closed-form placement uses the same private Eigen.
target_include_directories(decimate_reduce SYSTEM PRIVATE $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>)
target_compile_definitions(decimate_reduce PRIVATE EIGEN_MPL2_ONLY)
target_include_directories(decimate_seam PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_include_directories(decimate_seam SYSTEM PRIVATE $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>)
target_compile_definitions(decimate_seam PRIVATE EIGEN_MPL2_ONLY)
target_link_libraries(decimate_seam PUBLIC decimate::bake decimate::types decimate::math)

# --- decimate::pipeline — chartify + bakeAtlas, the two-call facade ---------
add_library(decimate_pipeline STATIC pipeline/pipeline.cpp)
add_library(decimate::pipeline ALIAS decimate_pipeline)
target_include_directories(decimate_pipeline PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(decimate_pipeline PUBLIC decimate::chart decimate::bake decimate::seam decimate::mesh decimate::types decimate::math)

# ---------------------------------------------------------------------------
# decimate / decimate::all — umbrella that re-exports every component, so a
# consumer that wants the whole pipeline links one target. Carries the version
# smoke-test TU. Pull individual modules via decimate::<module> instead.
# ---------------------------------------------------------------------------
add_library(decimate STATIC version.cpp)
add_library(decimate::all ALIAS decimate)
target_include_directories(decimate PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(decimate PUBLIC
    decimate::chart
    decimate::reduce
    decimate::bake
    decimate::seam
    decimate::pipeline
    decimate::io
    decimate::types decimate::math
)

# ---------------------------------------------------------------------------
# Install / export — makes the components consumable via find_package(decimate).
# decimate::math + decimate::types are exported too (PUBLIC deps of every component).
# ---------------------------------------------------------------------------
include(GNUInstallDirs)

# Export under the short, namespaced names decimate::<module> (NAMESPACE prepends
# to the EXPORT_NAME, not the internal target name — which stays decimate_<module>
# to avoid global-target collisions when consumed via add_subdirectory).
set_target_properties(decimate_types       PROPERTIES EXPORT_NAME types)
set_target_properties(decimate_math       PROPERTIES EXPORT_NAME math)
set_target_properties(decimate_mesh       PROPERTIES EXPORT_NAME mesh)
set_target_properties(decimate_io         PROPERTIES EXPORT_NAME io)
set_target_properties(decimate_raytracing PROPERTIES EXPORT_NAME raytracing)
set_target_properties(decimate_chart      PROPERTIES EXPORT_NAME chart)
set_target_properties(decimate_reduce     PROPERTIES EXPORT_NAME reduce)
set_target_properties(decimate_bake       PROPERTIES EXPORT_NAME bake)
set_target_properties(decimate_seam       PROPERTIES EXPORT_NAME seam)
set_target_properties(decimate_pipeline  PROPERTIES EXPORT_NAME pipeline)
set_target_properties(decimate            PROPERTIES EXPORT_NAME all)

install(TARGETS
        decimate_types decimate_math decimate_mesh decimate_io
        decimate_raytracing decimate_chart decimate_reduce decimate_bake decimate_seam decimate_pipeline decimate
    EXPORT decimateTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

# Public header tree: core/include/decimate/<module>/*.h  ->  <prefix>/include/decimate/...
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        FILES_MATCHING PATTERN "*.h")
