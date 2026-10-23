/**
 * decimate — appearance-preserving mesh simplification
 *
 * version 0.7
 *
 * Copyright (C) 2003-2026 Andor Greißl (Greissl)
 * Report bugs and download new versions at https://github.com/appdoo/decimate
 *
 * This library is distributed under the MIT License. See notice at the end of
 * this file.
 */
#pragma once

/// @file decimate.h
/// @brief Umbrella header + canonical version/license notice for the decimate
///        library.  Every other decimate source carries a one-line reference
///        back to this file; see the banner above for version and the MIT
///        notice below for license terms.
///
/// decimate is a modern C++20 port of a 2003 diploma thesis: an
/// appearance-preserving mesh-simplification + texture-atlas toolkit.  The
/// pipeline is:
///
///   load → reduce (QEM) → chartify → parametrize (LSCM) → pack (atlas)
///        → bake (raytraced normals) → render
///
/// The core library carries zero rendering/OS dependencies (no OpenGL, no
/// windowing, no globals) so it is testable and WebAssembly-compilable; all
/// rendering lives in the separate UI layer.
///
/// This header is the canonical version/license anchor and the public entry
/// point for the library version symbol (decimate::version), included as
/// `#include <decimate.h>`.

#define DECIMATE_VERSION_MAJOR  0
#define DECIMATE_VERSION_MINOR  8
#define DECIMATE_VERSION_STRING "0.8"

namespace decimate {
/// @brief Library version as MAJOR*100 + MINOR. A real (linkable) symbol, so it
///        also serves as the smoke-test that the decimate archive linked.
[[nodiscard]] int version() noexcept;
} // namespace decimate

/*
 * MIT License
 *
 * Copyright (C) 2003-2026 Andor Greißl (Greissl)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
