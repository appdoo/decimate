// This file is part of the decimate library; see decimate.h for version/license details.
/// @file tga_writer.cpp
/// @brief Uncompressed 24-bit TARGA writer — UI layer, no GL, no core dependency
///        beyond the `NormalMap` value type.

#include "tga_writer.h"

#include <array>
#include <cstdio>
#include <cstring>

namespace ui {

/// Writes a NormalMap as a top-left-origin 24-bit TGA (contract: tga_writer.h).
bool writeTGA(const std::string &path, const decimate::bake::NormalMap &nm) {
    if (nm.width <= 0 || nm.height <= 0) {
        return false;
    }

    // NOLINTNEXTLINE(cppcoreguidelines-owning-memory): FILE* — no C++ RAII wrapper
    std::FILE *fp = std::fopen(path.c_str(), "wb");
    if (fp == nullptr) {
        return false;
    }

    // TGA 18-byte header (uncompressed true-colour, image type 2, 24 bpp).
    // Reference: Paul Bourke, "TGA File Format Summary" (1996).
    //
    // Byte layout (all little-endian multi-byte fields):
    //   [0]   ID field length — 0 (no image ID)
    //   [1]   colour map type — 0 (no colour map)
    //   [2]   image type — 2 (uncompressed RGB)
    //   [3-7] colour map spec — all zeros
    //   [8-9] X-origin — 0
    //  [10-11] Y-origin — 0
    //  [12-13] width  (little-endian uint16)
    //  [14-15] height (little-endian uint16)
    //   [16]  pixel depth — 24
    //   [17]  image descriptor — 0x20 (bit 5 = top-left origin, bits[3:0]=0 attr bits)
    std::array<unsigned char, 18> hdr {};
    hdr[2] = 2;                                             // image type: uncompressed RGB
    hdr[12] = static_cast<unsigned char>(nm.width & 0xFF);  // width  low byte
    hdr[13] = static_cast<unsigned char>(nm.width >> 8);    // width  high byte
    hdr[14] = static_cast<unsigned char>(nm.height & 0xFF); // height low byte
    hdr[15] = static_cast<unsigned char>(nm.height >> 8);   // height high byte
    hdr[16] = 24;                                           // bits per pixel
    hdr[17] = 0x20;                                         // top-left origin

    std::fwrite(hdr.data(), 1, hdr.size(), fp);

    // Pixel data goes out row-major in B,G,R order (what TGA readers expect); row 0 is the
    // top of the image because header byte [17] declares a top-left origin.
    const std::size_t nPixels = static_cast<std::size_t>(nm.width) * static_cast<std::size_t>(nm.height);
    for (std::size_t i = 0; i < nPixels; ++i) {
        const auto                        &px = nm.pixels[i];
        const std::array<unsigned char, 3> bgr = { px.b, px.g, px.r };
        std::fwrite(bgr.data(), 1, bgr.size(), fp);
    }

    std::fclose(fp);
    return true;
}


/// Writes raw RGB bytes as a 24-bit TGA (contract: tga_writer.h).
bool writeTGA(const std::string &path, int w, int h,
              const std::vector<unsigned char> &rgb, bool bottomUp) {
    if (w <= 0 || h <= 0 || w > 0xFFFF || h > 0xFFFF ||
        rgb.size() != static_cast<std::size_t>(w) * static_cast<std::size_t>(h) * 3) {
        return false;
    }

    // NOLINTNEXTLINE(cppcoreguidelines-owning-memory): FILE* — no C++ RAII wrapper
    std::FILE *fp = std::fopen(path.c_str(), "wb");
    if (fp == nullptr) {
        return false;
    }

    std::array<unsigned char, 18> hdr {};
    hdr[2] = 2; // uncompressed true-colour
    hdr[12] = static_cast<unsigned char>(w & 0xFF);
    hdr[13] = static_cast<unsigned char>(w >> 8);
    hdr[14] = static_cast<unsigned char>(h & 0xFF);
    hdr[15] = static_cast<unsigned char>(h >> 8);
    hdr[16] = 24;
    hdr[17] = bottomUp ? 0x00 : 0x20; // bit 5: top-left origin
    std::fwrite(hdr.data(), 1, hdr.size(), fp);

    // TGA stores B,G,R; swap in one buffer so the pixel data goes out with a single fwrite.
    std::vector<unsigned char> bgr(rgb.size());
    for (std::size_t i = 0; i < rgb.size(); i += 3) {
        bgr[i] = rgb[i + 2];
        bgr[i + 1] = rgb[i + 1];
        bgr[i + 2] = rgb[i];
    }

    // fclose can fail on a buffered flush (disk full), so its result counts as part of success.
    const bool ok = std::fwrite(bgr.data(), 1, bgr.size(), fp) == bgr.size();
    return (std::fclose(fp) == 0) && ok;
}

} // namespace ui
