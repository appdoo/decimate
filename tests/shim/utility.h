// This file is part of the decimate library; see decimate.h for version/license details.
// Shim utility.h for the capture_bake_golden target.
//
// Placed on the include path BEFORE the legacy tree, this header replaces
// legacy2003/utility.h (which includes <windows.h>) with three no-op timer
// prototypes so the legacy texturesample.cpp compiles on macOS/clang without
// any modification to legacy2003/.
//
// The actual stubs (setup_timer, start_timer, stop_timer) are defined
// directly in capture_bake_golden.cpp, matching the in-TU shim pattern
// used by capture_golden.cpp for barycentricBadouelf.

#ifndef UTILITY_H
#define UTILITY_H

// No windows.h — the three functions are enough for compile-time satisfaction.
void setup_timer();
void start_timer();
double stop_timer();

#endif  // UTILITY_H
