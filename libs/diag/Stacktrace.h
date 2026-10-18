// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <string>
#include <vector>

namespace diag {

/**
 * Stacktrace utility for debugging and error reporting.
 *
 * Provides POSIX backtrace()-based stack capture with demangled symbol
 * resolution and a crash handler for common fatal signals.
 *
 * Only available on hosts where execinfo.h / backtrace() are present
 * (Linux, macOS). On other platforms the CMake platform gate produces a
 * no-op stub. The crash handler must be installed at the executable entry
 * point (never inside mesh-core) by calling installCrashHandler().
 */
class Stacktrace {
public:
    struct Frame {
        void*       address;        ///< Raw instruction pointer
        std::string symbol;         ///< Demangled function name
        std::string module;         ///< Library / executable basename
        std::string sourceLocation; ///< File:line if debug info is available
        std::size_t offset;         ///< Byte offset within the function
    };

    /**
     * Capture the current stacktrace.
     * @param skip       Number of innermost frames to drop (use 1 to omit
     *                   capture() itself from the result).
     * @param maxFrames  Upper bound on frames returned (default 64).
     * @return Vector of resolved stack frames; empty when unavailable.
     */
    static std::vector<Frame> capture(std::size_t skip = 0, std::size_t maxFrames = 64);

    /**
     * Format a captured stacktrace as a human-readable multi-line string.
     * @param frames   Frames from capture().
     * @param detailed Include raw addresses and byte offsets.
     * @return Formatted string terminated by a newline.
     */
    static std::string format(const std::vector<Frame>& frames, bool detailed = true);

    /**
     * Print the current stacktrace to stderr.
     * @param skip       Extra frames to skip on top of print() itself.
     * @param maxFrames  Upper bound on frames printed.
     */
    static void print(std::size_t skip = 1, std::size_t maxFrames = 64);

    /**
     * Install signal handlers that print a stacktrace on fatal crash.
     *
     * Handles: SIGSEGV, SIGABRT, SIGFPE, SIGILL, SIGBUS.
     * Call once from the executable entry point. Must NOT be called from
     * library code — this is an executable-boundary concern.
     */
    static void installCrashHandler();

private:
    static std::string demangle(const char* mangled);
    static void        signalHandler(int signal);
};

} // namespace diag
