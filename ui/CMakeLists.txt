cmake_minimum_required(VERSION 3.21)

# GLFW via FetchContent — only fetched when MESH_BUILD_UI=ON.
# Set doc/test/example build OFF for a lean dependency.
include(FetchContent)
FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        3.4
    GIT_SHALLOW    TRUE
)
set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)  # keep glfw out of the decimate package install
FetchContent_MakeAvailable(glfw)

find_package(OpenGL REQUIRED)

add_executable(mesh-viewer
    main.cpp
    viewer.cpp
    tga_writer.cpp
    shot_image.cpp
)

target_include_directories(mesh-viewer PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(mesh-viewer PRIVATE
    decimate
    decimate::math
    glfw
    OpenGL::GL
    diag
)

target_compile_definitions(mesh-viewer PRIVATE
    MESH_MODEL_DIR="${CMAKE_SOURCE_DIR}/model"
)

# Propagate STACKTRACE_ENABLED from the diag target.
# diag is STATIC (or INTERFACE no-op) depending on host — get_target_property
# gracefully returns NOTFOUND when the definition isn't set.
get_target_property(_diag_defs diag INTERFACE_COMPILE_DEFINITIONS)
if(_diag_defs MATCHES "STACKTRACE_ENABLED")
    target_compile_definitions(mesh-viewer PRIVATE STACKTRACE_ENABLED)
endif()
