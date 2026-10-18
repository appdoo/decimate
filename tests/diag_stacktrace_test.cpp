// This file is part of the decimate library; see decimate.h for version/license details.
// Spec 03-03 — T3 integration test for diag::Stacktrace::capture().
//
// On a backtrace-capable host (Linux, macOS) this test:
//   1. Calls diag::Stacktrace::capture() from a known named function.
//   2. Asserts at least one frame is returned.
//   3. Asserts that the calling function's demangled symbol appears in some
//      frame's symbol string.
//
// On hosts without backtrace support (STACKTRACE_ENABLED not defined) the
// test is compiled as a no-op so the build and test suite stay green.

#include <catch2/catch_test_macros.hpp>

#ifdef STACKTRACE_ENABLED

#include "Stacktrace.h"

#include <algorithm>
#include <iostream>
#include <string>

// A named function so its demangled symbol is predictable in the stacktrace.
// The name must be long enough to survive demangling unambiguously.
//
// Optimized builds would inline this into the TEST_CASE (dropping its frame) or turn the
// final call into a tail call (replacing its frame), so it is pinned noinline and the call
// is made non-tail by a side effect that follows it.
#if defined(_MSC_VER)
#  define MESH_TEST_NOINLINE __declspec(noinline)
#else
#  define MESH_TEST_NOINLINE [[gnu::noinline]]
#endif

static volatile int g_afterCapture = 0;

MESH_TEST_NOINLINE static std::vector<diag::Stacktrace::Frame> captureFromKnownCallee() {
    // skip=0 keeps captureFromKnownCallee itself in the frames.
    auto frames = diag::Stacktrace::capture(0);
    g_afterCapture = g_afterCapture + 1;  // volatile write: forbids tail-call / frame elision
    return frames;
}

TEST_CASE("diag::Stacktrace::capture returns frames with caller symbol", "[diag][integration]") {
    auto frames = captureFromKnownCallee();

    // Print the captured frames for evidence in the test output / transcript.
    std::cout << diag::Stacktrace::format(frames, /*detailed=*/true) << std::flush;

    REQUIRE_FALSE(frames.empty());

    // At least one frame must contain the known callee's name.
    const std::string needle = "captureFromKnownCallee";
    bool found = std::any_of(frames.begin(), frames.end(), [&](const diag::Stacktrace::Frame& f) {
        return f.symbol.find(needle) != std::string::npos;
    });

    INFO("Captured " << frames.size() << " frame(s). Looking for '" << needle << "' in symbols:");
    for (const auto& f : frames) {
        INFO("  " << f.symbol);
    }

    REQUIRE(found);
}

#else // STACKTRACE_ENABLED not defined

TEST_CASE("diag::Stacktrace::capture skipped — backtrace not available on this host", "[diag][integration]") {
    // No-op: platform gate excluded the real implementation.
    // The test reports as passing so the suite stays green on non-backtrace hosts.
    SUCCEED("backtrace not available — T3 skipped");
}

#endif // STACKTRACE_ENABLED
