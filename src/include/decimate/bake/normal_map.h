// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

/// @file normal_map.h
/// @brief Pure in-memory normal map value type (Spec 10).
///
/// `NormalMap` is a flat RGB pixel buffer representing the baked normal map
/// produced by `bake()`. It is a pure value type — no file I/O, no GL, no OS
/// dependencies. Writing to disk (TGA/PNG) is the UI layer's job (Spec 11).
///
/// Index convention: pixel (x, y) lives at `pixels[y * width + x]`.

#include <cstdint>
#include <vector>

namespace decimate::bake {

/// @brief 8-bit RGB colour with components in [0, 255].
///
/// Range-encoded normals use the legacy formula:
///   r = (unsigned char)((normal.x + 1.0f) * 127.0f)
/// so the neutral (flat surface facing the ray) encodes as (127, 127, 254)
/// and the background / no-data sentinel is (128, 128, 128).
struct Rgb {
    std::uint8_t r{};
    std::uint8_t g{};
    std::uint8_t b{};

    [[nodiscard]] constexpr bool operator==(const Rgb& o) const noexcept {
        return r == o.r && g == o.g && b == o.b;
    }
    [[nodiscard]] constexpr bool operator!=(const Rgb& o) const noexcept {
        return !(*this == o);
    }
};

/// @brief Range-encoded normal map: width × height pixels in row-major order.
///
/// The background sentinel (128, 128, 128) marks texels not covered by any
/// low-poly face triangle.  `dilate()` fills those from their neighbours.
///
/// ## Invariants
///
/// - `pixels.size() == static_cast<size_t>(width) * height`.
/// - All pixel values produced by `bake()` are in [0, 255] per channel.
struct NormalMap {
    int              width  = 0;
    int              height = 0;
    std::vector<Rgb> pixels;  ///< Row-major; pixel (x,y) at pixels[y*width+x].
};

}  // namespace decimate::bake
