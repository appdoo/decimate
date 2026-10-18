// This file is part of the decimate library; see decimate.h for version/license details.
// Custom Catch2 test runner entry point.
//
// We provide our own main() (linking Catch2::Catch2, not Catch2WithMain) so
// we can install the diag crash handler before any test runs.  The crash
// handler prints a demangled stacktrace on SIGSEGV/SIGABRT/SIGFPE/SIGILL/SIGBUS,
// which is invaluable when a test triggers a fatal signal.
//
// Spec 03-03 — Boundary Contract: installCrashHandler() is called here, at the
// executable entry point.  mesh-core never links diag and never calls into it.

#include <catch2/catch_session.hpp>

#include <decimate.h>

#ifdef STACKTRACE_ENABLED
#include "Stacktrace.h"
#endif

int main(int argc, char* argv[]) {
#ifdef STACKTRACE_ENABLED
    diag::Stacktrace::installCrashHandler();
#endif

    return Catch::Session().run(argc, argv);
}

// ---------------------------------------------------------------------------
// mesh-core smoke tests
// ---------------------------------------------------------------------------
// These tests verify that mesh-core is a real, linkable STATIC library and
// that its public API is reachable from the test runner.

#include <catch2/catch_test_macros.hpp>

TEST_CASE("decimate::version returns expected version", "[core][smoke]") {
    REQUIRE(decimate::version() == 6);
}
