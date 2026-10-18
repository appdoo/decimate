# diag — host-only stacktrace / crash-handler library.
#
# Platform gate mirrors SetupStacktrace.cmake from the author's unreal_mmpg project.
# On hosts where execinfo.h / backtrace() are available (Linux, macOS/Darwin) this
# produces a real STATIC library compiled with STACKTRACE_ENABLED.
# On other hosts (Windows, Emscripten, …) it produces a no-op INTERFACE target so
# the rest of the build graph resolves without change.
#
# Boundary Contract: decimate must NOT link diag.
# Only executables (test runner, future UI app) link diag.

cmake_minimum_required(VERSION 3.21)

include(CheckIncludeFileCXX)
include(CheckFunctionExists)

check_include_file_cxx(execinfo.h HAVE_EXECINFO_H)
check_function_exists(backtrace HAVE_BACKTRACE)

if(HAVE_EXECINFO_H AND HAVE_BACKTRACE)
    message(STATUS "[diag] backtrace support detected — building real stacktrace library")

    add_library(diag STATIC
        Stacktrace.cpp
    )

    target_compile_definitions(diag PUBLIC STACKTRACE_ENABLED)

    # Consumers can #include "Stacktrace.h" via the diag include dir.
    target_include_directories(diag PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

    target_compile_features(diag PUBLIC cxx_std_20)

    # Export all symbols from the executable so dladdr() can resolve demangled
    # function names in the test runner and other executables that link diag.
    # -rdynamic on Linux; -Wl,-export_dynamic on macOS/clang (AppleClang ignores
    # -rdynamic silently, but -Wl,-export_dynamic is the correct spelling).
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_options(diag INTERFACE -rdynamic)
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        target_link_options(diag INTERFACE -Wl,-export_dynamic)
    endif()

    # Always include debug symbols so demangled caller symbols resolve in T3.
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(diag PRIVATE -g)
    elseif(MSVC)
        target_compile_options(diag PRIVATE /Zi)
    endif()

else()
    message(WARNING "[diag] backtrace not available — producing no-op diag target")

    # No-op INTERFACE target so the CMake graph resolves on every host.
    # Stacktrace.cpp still compiles (the #else branch provides stubs) but
    # rather than compiling a stub object on non-backtrace hosts we just
    # leave the INTERFACE empty — no sources, no defines, no link flags.
    add_library(diag INTERFACE)
endif()
