cmake_minimum_required(VERSION 3.21)

project(decimate VERSION 0.7.0 LANGUAGES CXX)

# NOTE: the canonical build is left at the default optimization level, because the
# exact-hash characterization goldens (LSCM UV hashes, stacktrace) drift under -O2
# float reordering. Build BENCHMARKS in a dedicated Release tree instead:
#     cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release && cmake --build build-rel
# (Making Release the default would need tolerance-based float goldens.)

# C++20, portable, no compiler extensions
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Emit compile_commands.json — required by clang-tidy and clangd (editor diagnostics).
# Spec 03-04 owns this flag; see also scripts/clang-tidy-diff.sh.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(GNUInstallDirs)

# Consumers (add_subdirectory / FetchContent of decimate) default these OFF, so they
# don't build our tests/viewer or fetch Catch2. They are ON only for a top-level build.
option(DECIMATE_BUILD_TESTS "Build decimate tests" ${PROJECT_IS_TOP_LEVEL})
# Worker threads for the normal bake (std::thread only). OFF = identical code, one worker on the
# caller's thread. Forced OFF for Emscripten without pthreads (-pthread in the flags) - Spec 22-01.
option(DECIMATE_THREADS     "Parallel normal bake on std::thread (OFF: one worker, same code)" ON)
if(EMSCRIPTEN AND NOT CMAKE_CXX_FLAGS MATCHES "pthread")
    set(DECIMATE_THREADS OFF CACHE BOOL "Parallel normal bake on std::thread" FORCE)
endif()
option(MESH_BUILD_UI        "Build the GLFW interactive viewer (mesh-viewer)" ${PROJECT_IS_TOP_LEVEL})

# Library subdirectories. diag is a host-only leaf; the decimate library itself
# (all decimate::<module> components) lives in src/.
add_subdirectory(libs/diag)
add_subdirectory(src)

if(DECIMATE_BUILD_TESTS)
    # Catch2 v3 via FetchContent — only fetched for a tests build.
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.15.0
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(Catch2)
    enable_testing()
    include(CTest)
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
    include(Catch)
    add_subdirectory(tests)
endif()

if(MESH_BUILD_UI)
    add_subdirectory(ui)
endif()

# ---------------------------------------------------------------------------
# Package config — find_package(decimate [COMPONENTS ...]) for external consumers.
# The component install(TARGETS ... EXPORT decimateTargets) calls live in
# core/; here we install the export set + generated config.
# ---------------------------------------------------------------------------
include(CMakePackageConfigHelpers)
set(DECIMATE_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/decimate")

install(EXPORT decimateTargets
        NAMESPACE decimate::
        FILE decimateTargets.cmake
        DESTINATION "${DECIMATE_CMAKE_DIR}")

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/decimateConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/decimateConfig.cmake"
    INSTALL_DESTINATION "${DECIMATE_CMAKE_DIR}")

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/decimateConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/decimateConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/decimateConfigVersion.cmake"
    DESTINATION "${DECIMATE_CMAKE_DIR}")
