#pragma once
/// @file shot_image.h
/// @brief GL-free RGB image helpers for the `--shot` README images: crop, per-pixel diff,
/// captioned side-by-side composition (tiny built-in bitmap font). UI layer only.

#include <string>
#include <utility>
#include <vector>

namespace ui {

/// @brief Tightly packed 8-bit RGB image, row 0 = TOP row.
struct Image {
    int                        w = 0;
    int                        h = 0;
    std::vector<unsigned char> rgb;  ///< w*h*3 bytes
};

/// @brief Copy the `cw`x`ch` rectangle whose top-left is (x0,y0). The rectangle must lie inside `src`.
[[nodiscard]] Image cropImage(const Image& src, int x0, int y0, int cw, int ch);

/// @brief How much two equally sized images differ.
struct DiffStats {
    long pixelsOver = 0;  ///< pixels whose summed |dR|+|dG|+|dB| is >= the threshold
    int  maxDiff    = 0;  ///< largest summed |dR|+|dG|+|dB| over all pixels
};

/// @brief Compare two equally sized images; `threshold` is on the summed channel difference.
[[nodiscard]] DiffStats diffStats(const Image& a, const Image& b, int threshold);

/// @brief Top-left of the `win`x`win` window (on a `stride` grid) holding the most pixels that
/// differ by >= `threshold`. Equal sizes required; `win` must fit inside the image.
[[nodiscard]] std::pair<int, int> mostDifferentWindow(const Image& a, const Image& b,
                                                      int win, int stride, int threshold);

/// @brief Box-filter downscale by an integer `factor` (anti-aliases a supersampled render).
/// Width/height must be multiples of `factor`.
[[nodiscard]] Image downscale(const Image& src, int factor);

/// @brief Two images side by side on a white canvas, each under a burned-in caption.
/// Captions use lowercase a-z, space, '-', ':' (the built-in font); anything else renders blank.
[[nodiscard]] Image sideBySide(const Image& left, const Image& right,
                               const std::string& leftCaption, const std::string& rightCaption);

/// @brief Write `img` as a top-left-origin 24-bit TGA.
[[nodiscard]] bool writeTGA(const std::string& path, const Image& img);

}  // namespace ui
