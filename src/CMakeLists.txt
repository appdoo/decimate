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
#   chart  -> mesh, math, types
#   reduce -> mesh, math, types
#   bake   -> mesh, raytracing, types, math
#
# `decimate` (alias `decimate::all`) is an umbrella that re-exports every
# component, so a consumer wanting the whole pipeline links a single target.
# It carries the version smoke-test TU.
# ---------------------------------------------------------------------------

# --- decimate::math — foundation: vector math + numeric utilities -----------
# Vec2/3/4, Mat4, Quat, AABB, Plane (absorbed from the former pure-math) +
# cg_solver / sparse_matrix / sparse_bit_set / associative_heap / mat4-invert.
# All namespace decimate::math. The single numeric foundation.
add_library(decimate_math STATIC
    math/cg_solver.cpp
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
target_include_directories(decimate_io PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(decimate_io PUBLIC decimate::types decimate::math)

# --- decimate::raytracing — pure ray functions -----------------------------
add_library(decimate_raytracing STATIC raytracing/raytracing.cpp)
add_library(decimate::raytracing ALIAS decimate_raytracing)
target_include_directories(decimate_raytracing PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(decimate_raytracing PUBLIC decimate::types decimate::math)

# --- decimate::chart — segmentation + LSCM + atlas packing -----------------
add_library(decimate_chart STATIC
    chart/segmenter.cpp
    chart/conformal_mapper.cpp
    chart/atlas_packer.cpp
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
set_target_properties(decimate            PROPERTIES EXPORT_NAME all)

install(TARGETS
        decimate_types decimate_math decimate_mesh decimate_io
        decimate_raytracing decimate_chart decimate_reduce decimate_bake decimate
    EXPORT decimateTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

# Public header tree: core/include/decimate/<module>/*.h  ->  <prefix>/include/decimate/...
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        FILES_MATCHING PATTERN "*.h")
