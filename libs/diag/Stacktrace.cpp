// This file is part of the decimate library; see decimate.h for version/license details.
#include "Stacktrace.h"

#ifdef STACKTRACE_ENABLED

#include <execinfo.h>
#include <cxxabi.h>
#include <dlfcn.h>
#include <signal.h>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <iomanip>

#ifdef __APPLE__
#include <unistd.h>
#endif

namespace diag {

std::string Stacktrace::demangle(const char* mangled) {
    if (!mangled) {
        return "??";
    }

    int   status    = 0;
    char* demangled = abi::__cxa_demangle(mangled, nullptr, nullptr, &status);

    if (status == 0 && demangled) {
        std::string result(demangled);
        free(demangled);
        return result;
    }

    return mangled;
}

std::vector<Stacktrace::Frame> Stacktrace::capture(std::size_t skip, std::size_t maxFrames) {
    std::vector<Frame> frames;

    std::vector<void*> buffer(maxFrames + skip);
    int count = backtrace(buffer.data(), static_cast<int>(buffer.size()));

    if (count <= 0) {
        return frames;
    }

    char** symbols = backtrace_symbols(buffer.data(), count);
    if (!symbols) {
        return frames;
    }

    for (int i = static_cast<int>(skip); i < count; ++i) {
        Frame frame;
        frame.address = buffer[i];
        frame.offset  = 0;

        Dl_info info;
        if (dladdr(buffer[i], &info)) {
            if (info.dli_fname) {
                frame.module = info.dli_fname;
                const char* lastSlash = strrchr(info.dli_fname, '/');
                if (lastSlash) {
                    frame.module = lastSlash + 1;
                }
            }

            if (info.dli_sname) {
                frame.symbol = demangle(info.dli_sname);

                if (info.dli_saddr) {
                    frame.offset = static_cast<std::size_t>(
                        reinterpret_cast<char*>(buffer[i]) -
                        reinterpret_cast<char*>(info.dli_saddr)
                    );
                }
            } else {
                // Fall back to parsing backtrace_symbols output.
                // Format: "<module> <address> <symbol> + <offset>"
                std::string symbolStr(symbols[i]);

                std::size_t symbolStart = symbolStr.find(' ');
                if (symbolStart != std::string::npos) {
                    symbolStart = symbolStr.find(' ', symbolStart + 1);
                    if (symbolStart != std::string::npos) {
                        std::size_t symbolEnd = symbolStr.find(" + ", symbolStart);
                        if (symbolEnd != std::string::npos) {
                            std::string mangledSymbol =
                                symbolStr.substr(symbolStart + 1, symbolEnd - symbolStart - 1);
                            frame.symbol = demangle(mangledSymbol.c_str());

                            std::size_t offsetStart = symbolEnd + 3;
                            std::size_t offsetEnd   = symbolStr.find(' ', offsetStart);
                            if (offsetEnd == std::string::npos) {
                                offsetEnd = symbolStr.length();
                            }
                            std::string offsetStr =
                                symbolStr.substr(offsetStart, offsetEnd - offsetStart);
                            frame.offset = std::stoull(offsetStr, nullptr, 0);
                        }
                    }
                }
            }
        }

        if (frame.symbol.empty()) {
            frame.symbol = symbols[i];
        }

        frames.push_back(frame);
    }

    free(symbols);
    return frames;
}

std::string Stacktrace::format(const std::vector<Frame>& frames, bool detailed) {
    std::ostringstream oss;
    oss << "Stacktrace (" << frames.size() << " frames):\n";

    for (std::size_t i = 0; i < frames.size(); ++i) {
        const auto& frame = frames[i];

        oss << "  #" << std::setw(2) << i << " ";

        if (detailed) {
            oss << "0x" << std::hex << std::setw(16) << std::setfill('0')
                << reinterpret_cast<uintptr_t>(frame.address) << std::dec << " ";
        }

        if (!frame.module.empty()) {
            oss << "[" << frame.module << "] ";
        }

        oss << frame.symbol;

        if (detailed && frame.offset > 0) {
            oss << " + 0x" << std::hex << frame.offset << std::dec;
        }

        if (!frame.sourceLocation.empty()) {
            oss << " at " << frame.sourceLocation;
        }

        oss << "\n";
    }

    return oss.str();
}

void Stacktrace::print(std::size_t skip, std::size_t maxFrames) {
    auto        frames    = capture(skip + 1, maxFrames); // +1 to skip print() itself
    std::string formatted = format(frames, true);
    std::cerr << formatted << std::flush;
}

void Stacktrace::signalHandler(int sig) {
    const char* signalName = "UNKNOWN";
    switch (sig) {
        case SIGSEGV: signalName = "SIGSEGV (Segmentation Fault)"; break;
        case SIGABRT: signalName = "SIGABRT (Abort)";               break;
        case SIGFPE:  signalName = "SIGFPE (Floating Point Exception)"; break;
        case SIGILL:  signalName = "SIGILL (Illegal Instruction)";   break;
        case SIGBUS:  signalName = "SIGBUS (Bus Error)";             break;
    }

    std::cerr << "\n========================================\n";
    std::cerr << "FATAL ERROR: Caught signal " << sig << " - " << signalName << "\n";
    std::cerr << "========================================\n\n";

    print(1, 64); // skip the signal handler frame itself

    std::cerr << "\n========================================\n";
    std::cerr << "Terminating...\n";
    std::cerr << "========================================\n" << std::flush;

    ::signal(sig, SIG_DFL);
    ::raise(sig);
}

void Stacktrace::installCrashHandler() {
    signal(SIGSEGV, signalHandler);
    signal(SIGABRT, signalHandler);
    signal(SIGFPE,  signalHandler);
    signal(SIGILL,  signalHandler);
    signal(SIGBUS,  signalHandler);
}

} // namespace diag

#else // STACKTRACE_ENABLED not defined — no-op stubs

namespace diag {

std::vector<Stacktrace::Frame> Stacktrace::capture(std::size_t, std::size_t) {
    return {};
}

std::string Stacktrace::format(const std::vector<Frame>&, bool) {
    return "Stacktrace not available on this platform\n";
}

void Stacktrace::print(std::size_t, std::size_t) {
    std::cerr << "Stacktrace not available on this platform\n";
}

void Stacktrace::installCrashHandler() {
    // no-op on platforms without backtrace support
}

std::string Stacktrace::demangle(const char* mangled) {
    return mangled ? mangled : "??";
}

void Stacktrace::signalHandler(int sig) {
    (void)sig;
}

} // namespace diag

#endif // STACKTRACE_ENABLED
