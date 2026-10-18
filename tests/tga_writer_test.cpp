// This file is part of the decimate library; see decimate.h for version/license details.
/// @file tga_writer_test.cpp
/// @brief Spec 11 B1 — TGA writer round-trip: 18-byte header validation and
///        4-pixel BGR encoding check.  Fully off-GPU; no GL or GLFW required.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include <decimate/bake/normal_map.h>  // decimate::bake::NormalMap, Rgb
#include "tga_writer.h"  // ui::writeTGA

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

namespace {

/// Read a complete file into a byte buffer.
std::vector<uint8_t> readFile(const std::string& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) {
        return {};
    }
    const auto size = static_cast<std::size_t>(f.tellg());
    f.seekg(0, std::ios::beg);
    std::vector<uint8_t> buf(size);
    f.read(reinterpret_cast<char*>(buf.data()),  // NOLINT: reading raw bytes
           static_cast<std::streamsize>(size));
    return buf;
}

/// Little-endian uint16 from two bytes.
uint16_t le16(const std::vector<uint8_t>& buf, std::size_t off) {
    return static_cast<uint16_t>(buf[off]) |
           (static_cast<uint16_t>(buf[off + 1]) << 8);
}

}  // namespace

// ---------------------------------------------------------------------------
// B1 — TGA header + 4-pixel BGR round-trip
// ---------------------------------------------------------------------------

TEST_CASE("writeTGA round-trips a 2x2 NormalMap correctly", "[tga][B1]") {
    // Build a 2×2 NormalMap with four known, distinct RGB values.
    // pixel (x,y) lives at pixels[y * width + x].
    decimate::bake::NormalMap nm;
    nm.width  = 2;
    nm.height = 2;
    nm.pixels = {
        {  10,  20,  30 },   // (0,0) top-left
        {  40,  50,  60 },   // (1,0) top-right
        {  70,  80,  90 },   // (0,1) bottom-left
        { 100, 110, 120 },   // (1,1) bottom-right
    };

    // Write to a temp file.
    const std::string tmpPath = "/tmp/tga_writer_test_b1.tga";
    REQUIRE(ui::writeTGA(tmpPath, nm));

    // Re-read the raw bytes.
    const std::vector<uint8_t> raw = readFile(tmpPath);

    // -----------------------------------------------------------------------
    // Header validation (18 bytes).
    // -----------------------------------------------------------------------
    REQUIRE(raw.size() >= 18 + 4 * 3);  // header + 4 pixels × 3 bytes

    // Byte [0] — ID field length (we write 0).
    CHECK(raw[0] == 0);

    // Byte [1] — colour map type (0 = no colour map).
    CHECK(raw[1] == 0);

    // Byte [2] — image type 2 = uncompressed true-colour RGB.
    CHECK(raw[2] == 2);

    // Bytes [3-7] — colour map spec, all zero for type-2 images.
    CHECK(raw[3] == 0);
    CHECK(raw[4] == 0);
    CHECK(raw[5] == 0);
    CHECK(raw[6] == 0);
    CHECK(raw[7] == 0);

    // Bytes [8-9] — X-origin (0).
    CHECK(le16(raw, 8) == 0);

    // Bytes [10-11] — Y-origin (0).
    CHECK(le16(raw, 10) == 0);

    // Bytes [12-13] — width (little-endian uint16).
    CHECK(le16(raw, 12) == 2);

    // Bytes [14-15] — height (little-endian uint16).
    CHECK(le16(raw, 14) == 2);

    // Byte [16] — bits per pixel: 24.
    CHECK(raw[16] == 24);

    // Byte [17] — image descriptor: 0x20 = top-left origin.
    CHECK(raw[17] == 0x20);

    // -----------------------------------------------------------------------
    // Pixel data: 4 pixels, each 3 bytes, B G R order.
    // Row 0 = top of image (top-left-origin set in header).
    // pixel[y*w+x] → raw byte at 18 + (y*w+x)*3.
    // -----------------------------------------------------------------------
    // Known RGB → expected BGR bytes on disk.
    struct PixelCheck {
        uint8_t expectedB;
        uint8_t expectedG;
        uint8_t expectedR;
        uint8_t srcR;
        uint8_t srcG;
        uint8_t srcB;
    };

    const std::array<PixelCheck, 4> checks = {{
        {  30,  20,  10,  10,  20,  30 },  // (0,0)
        {  60,  50,  40,  40,  50,  60 },  // (1,0)
        {  90,  80,  70,  70,  80,  90 },  // (0,1)
        { 120, 110, 100, 100, 110, 120 },  // (1,1)
    }};

    for (std::size_t i = 0; i < checks.size(); ++i) {
        const std::size_t byteOff = 18 + i * 3;
        CAPTURE(i);
        // File byte order: B, G, R.
        CHECK(raw[byteOff + 0] == checks[i].expectedB);
        CHECK(raw[byteOff + 1] == checks[i].expectedG);
        CHECK(raw[byteOff + 2] == checks[i].expectedR);
    }

    // Cleanup.
    std::remove(tmpPath.c_str());
}

TEST_CASE("writeTGA returns false for empty NormalMap", "[tga][B1]") {
    decimate::bake::NormalMap empty;  // width=0, height=0
    CHECK_FALSE(ui::writeTGA("/tmp/should_not_exist_11b1.tga", empty));
}
