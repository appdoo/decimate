// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file tga_writer.h
/// @brief UI-layer 24-bit uncompressed TARGA writer.
///
/// Writes a `decimate::bake::NormalMap` to disk as a standard uncompressed true-colour
/// TGA (image type 2, 24-bit BGR, top-left origin).  This is intentionally in the
/// UI layer — core has no file I/O (anti-pattern 14 / boundary contract Spec 11).

#include <string>
#include <vector>

#include <decimate/bake/normal_map.h>  // decimate::bake::NormalMap — core value type (no GL, no OS dep)

namespace ui {

/// @brief Write `nm` as an uncompressed 24-bit BGR TARGA file.
///
/// TGA header: image type 2 (uncompressed true-colour), 24 bpp, top-left origin
/// (descriptor byte 0x20).  Pixel data written row-major B,G,R so a standard TGA
/// viewer decodes it correctly.  Row 0 of `nm.pixels` maps to the top row of the
/// image (top-left-origin flag makes this explicit to the reader).
///
/// @param path  Destination path (created/overwritten).
/// @param nm    Source map; `nm.pixels.size()` must equal `nm.width * nm.height`.
/// @returns     `true` on success, `false` if the file could not be opened.
[[nodiscard]] bool writeTGA(const std::string& path, const decimate::bake::NormalMap& nm);

/// @brief Write tightly packed 8-bit RGB pixels (3 bytes/pixel, `w*h*3` bytes) as a 24-bit TGA.
///
/// @param bottomUp  true if row 0 of `rgb` is the BOTTOM row (glReadPixels order); the
///                  header then declares a bottom-left origin so viewers show it upright.
/// @returns         `true` on success; `false` on bad size or unwritable path.
[[nodiscard]] bool writeTGA(const std::string& path, int w, int h,
                            const std::vector<unsigned char>& rgb, bool bottomUp);

}  // namespace ui
