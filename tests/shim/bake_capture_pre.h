// This file is part of the decimate library; see decimate.h for version/license details.
// bake_capture_pre.h — force-included before compiling legacy texturesample.cpp.
//
// Defines the include guards for headers that cannot compile on macOS/clang:
//   - utility.h: includes <windows.h>
//   - algorithm.h → halfedge.h → cset.h: CSet has an ill-formed template body
//     (uses undeclared 'array' in CSet::operator[]) that causes a hard error
//     in strict C++ mode even without instantiation.
//
// Also provides the timer prototypes (so texturesample.cpp sees them) and
// suppresses the 'register' keyword warning from vec.h.
//
// This header is injected via -include, so it is processed before any #include
// in the compiled TU.  Legacy source files are NOT modified.

// -- Block windows.h / utility.h
#ifndef UTILITY_H
#define UTILITY_H
void setup_timer();
void start_timer();
double stop_timer();
#endif

// -- Block algorithm.h (pulls halfedge.h → cset.h which has broken CSet::operator[])
#ifndef ALGORITHM_H
#define ALGORITHM_H
#endif

// -- Block halfedge.h (transitively pulled by algorithm.h)
#ifndef HALFEDGE_H
#define HALFEDGE_H
#endif

// -- Block cset.h (transitively pulled by halfedge.h)
#ifndef CSET_H
#define CSET_H
#endif
