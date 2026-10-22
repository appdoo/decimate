// This file is part of the decimate library; see decimate.h for version/license details.
/// @file shot_image.cpp
/// @brief GL-free RGB image helpers for the `--shot` images (contracts: shot_image.h).

#include "shot_image.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <string_view>

#include "tga_writer.h"

namespace ui {

namespace {

/// Text-rendering metrics of the built-in bitmap font, in pixels.
constexpr int kTextScale = 3;
constexpr int kGlyphW = 5;
constexpr int kGlyphH = 9; // 2 ascender rows + 5 x-height rows + 2 descender rows
constexpr int kAdvance = (kGlyphW + 1) * kTextScale;
constexpr int kBarH = (kGlyphH + 4) * kTextScale;
constexpr int kGap = 24;

/// One glyph: optional ascender rows (2), x-height rows (5), optional descender rows (2),
/// each row a 5-char string of '#' / '.', comma separated within a part.
struct Glyph {
    char             c;
    std::string_view asc;
    std::string_view xh;
    std::string_view desc;
};

constexpr std::array<Glyph, 23> kFont { {
    { ' ', "", ".....,.....,.....,.....,.....", "" },
    { '-', "", ".....,.....,#####,.....,.....", "" },
    { ':', "", ".....,..#..,.....,..#..,.....", "" },
    { 'a', "", ".###.,....#,.####,#...#,.####", "" },
    { 'b', "#....,#....", "#.##.,##..#,#...#,##..#,#.##.", "" },
    { 'c', "", ".###.,#...#,#....,#...#,.###.", "" },
    { 'd', "....#,....#", ".##.#,#..##,#...#,#..##,.##.#", "" },
    { 'e', "", ".###.,#...#,#####,#....,.###.", "" },
    { 'f', "..##.,.#...", "####.,.#...,.#...,.#...,.#...", "" },
    { 'g', "", ".####,#...#,#...#,#...#,.####", "....#,.###." },
    { 'h', "#....,#....", "#.##.,##..#,#...#,#...#,#...#", "" },
    { 'i', ".....,..#..", ".##..,..#..,..#..,..#..,.###.", "" },
    { 'l', ".##..,..#..", "..#..,..#..,..#..,..#..,.###.", "" },
    { 'm', "", "##.#.,#.#.#,#.#.#,#...#,#...#", "" },
    { 'n', "", "#.##.,##..#,#...#,#...#,#...#", "" },
    { 'o', "", ".###.,#...#,#...#,#...#,.###.", "" },
    { 'p', "", "#.##.,##..#,#...#,#...#,####.", "#....,#...." },
    { 'r', "", "#.##.,##..#,#....,#....,#....", "" },
    { 's', "", ".####,#....,.###.,....#,####.", "" },
    { 't', ".#...,.#...", "####.,.#...,.#...,.#..#,..##.", "" },
    { 'u', "", "#...#,#...#,#...#,#..##,.##.#", "" },
    { 'y', "", "#...#,#...#,#...#,#...#,.####", "....#,.###." },
    { 'v', "", "#...#,#...#,#...#,.#.#.,..#..", "" },
} };

} // namespace


/// Expands a glyph into 9 rows of 5 chars (blank where a part is absent).
static std::array<std::string_view, kGlyphH> glyphRows(const Glyph &g) {
    std::array<std::string_view, kGlyphH> rows;
    rows.fill(".....");
    // Parts have fixed row slots (ascender 0, x-height 2, descender 7) so letters share a baseline.
    const auto split = [&rows](std::string_view part, int first) {
        for (int i = 0; !part.empty(); ++i) {
            const auto comma = part.find(',');
            rows[static_cast<std::size_t>(first + i)] = part.substr(0, 5);
            part = (comma == std::string_view::npos) ? std::string_view {} : part.substr(comma + 1);
        }
    };
    split(g.asc, 0);
    split(g.xh, 2);
    split(g.desc, 7);
    return rows;
}


/// Writes a grey pixel; out-of-range coordinates are clipped so callers can draw past the edge.
static void putPixel(Image &im, int x, int y, unsigned char v) {
    if (x < 0 || y < 0 || x >= im.w || y >= im.h) {
        return;
    }
    const auto i = (static_cast<std::size_t>(y) * static_cast<std::size_t>(im.w) + static_cast<std::size_t>(x)) * 3;
    im.rgb[i] = im.rgb[i + 1] = im.rgb[i + 2] = v;
}


/// Fills a grey rectangle (clipped by putPixel).
static void fillRect(Image &im, int x0, int y0, int rw, int rh, unsigned char v) {
    for (int y = y0; y < y0 + rh; ++y) {
        for (int x = x0; x < x0 + rw; ++x) {
            putPixel(im, x, y, v);
        }
    }
}


/// Copies `src` into `dst` at (x0,y0) row by row; the rectangle must lie inside `dst`.
static void blit(Image &dst, const Image &src, int x0, int y0) {
    for (int y = 0; y < src.h; ++y) {
        const auto s = static_cast<std::size_t>(y) * static_cast<std::size_t>(src.w) * 3;
        const auto d = (static_cast<std::size_t>(y0 + y) * static_cast<std::size_t>(dst.w) + static_cast<std::size_t>(x0)) * 3;
        std::copy_n(src.rgb.begin() + static_cast<std::ptrdiff_t>(s), static_cast<std::size_t>(src.w) * 3,
                    dst.rgb.begin() + static_cast<std::ptrdiff_t>(d));
    }
}


/// Burns `text` into `im` with the built-in font; characters outside the font advance the pen but draw nothing.
static void drawText(Image &im, int x0, int y0, const std::string &text, unsigned char ink) {
    int pen = x0;
    for (const char ch : text) {
        const auto it = std::find_if(kFont.begin(), kFont.end(), [ch](const Glyph &g) { return g.c == ch; });
        if (it != kFont.end()) {
            const auto rows = glyphRows(*it);
            for (int gy = 0; gy < kGlyphH; ++gy) {
                for (int gx = 0; gx < kGlyphW; ++gx) {
                    if (rows[static_cast<std::size_t>(gy)][static_cast<std::size_t>(gx)] == '#') {
                        fillRect(im, pen + gx * kTextScale, y0 + gy * kTextScale, kTextScale, kTextScale, ink);
                    }
                }
            }
        }
        pen += kAdvance;
    }
}


Image cropImage(const Image &src, int x0, int y0, int cw, int ch) {
    Image out;
    out.w = cw;
    out.h = ch;
    out.rgb.resize(static_cast<std::size_t>(cw) * static_cast<std::size_t>(ch) * 3);
    for (int y = 0; y < ch; ++y) {
        const auto s = (static_cast<std::size_t>(y0 + y) * static_cast<std::size_t>(src.w) + static_cast<std::size_t>(x0)) * 3;
        std::copy_n(src.rgb.begin() + static_cast<std::ptrdiff_t>(s), static_cast<std::size_t>(cw) * 3,
                    out.rgb.begin() + static_cast<std::ptrdiff_t>(static_cast<std::size_t>(y) * static_cast<std::size_t>(cw) * 3));
    }
    return out;
}


/// Summed |dR|+|dG|+|dB| of pixel number `px` (row-major index).
[[nodiscard]] static int pixelDiff(const Image &a, const Image &b, std::size_t px) {
    const std::size_t i = px * 3;
    return std::abs(a.rgb[i] - b.rgb[i]) + std::abs(a.rgb[i + 1] - b.rgb[i + 1]) +
           std::abs(a.rgb[i + 2] - b.rgb[i + 2]);
}


DiffStats diffStats(const Image &a, const Image &b, int threshold) {
    DiffStats st;

    const auto n = static_cast<std::size_t>(a.w) * static_cast<std::size_t>(a.h);
    for (std::size_t p = 0; p < n; ++p) {
        const int d = pixelDiff(a, b, p);
        st.maxDiff = std::max(st.maxDiff, d);
        st.pixelsOver += (d >= threshold) ? 1 : 0;
    }
    return st;
}


std::pair<int, int> mostDifferentWindow(const Image &a, const Image &b, int win, int stride, int threshold) {
    // Summed-area table of "differs" so every candidate window costs O(1).
    const auto       satWidth = static_cast<std::size_t>(a.w) + 1;
    std::vector<int> sat(satWidth * (static_cast<std::size_t>(a.h) + 1), 0);
    for (int y = 0; y < a.h; ++y) {
        for (int x = 0; x < a.w; ++x) {
            const auto hit = pixelDiff(a, b, static_cast<std::size_t>(y) * static_cast<std::size_t>(a.w) + static_cast<std::size_t>(x)) >= threshold ? 1 : 0;
            const auto at = [&](int xx, int yy) -> int & {
                return sat[static_cast<std::size_t>(yy) * satWidth + static_cast<std::size_t>(xx)];
            };
            at(x + 1, y + 1) = hit + at(x, y + 1) + at(x + 1, y) - at(x, y);
        }
    }
    std::pair<int, int> best { 0, 0 };
    int                 bestScore = -1;

    // Strict '>' keeps the first (top-left-most) window on ties.
    for (int y0 = 0; y0 + win <= a.h; y0 += stride) {
        for (int x0 = 0; x0 + win <= a.w; x0 += stride) {
            const auto at = [&](int xx, int yy) {
                return sat[static_cast<std::size_t>(yy) * satWidth + static_cast<std::size_t>(xx)];
            };
            const int score = at(x0 + win, y0 + win) - at(x0, y0 + win) - at(x0 + win, y0) + at(x0, y0);
            if (score > bestScore) {
                bestScore = score;
                best = { x0, y0 };
            }
        }
    }
    return best;
}


Image sideBySide(const Image &left, const Image &right, const std::string &leftCaption,
                 const std::string &rightCaption) {
    Image out;
    out.w = left.w + right.w + 3 * kGap;
    out.h = std::max(left.h, right.h) + kBarH + 2 * kGap;
    out.rgb.assign(static_cast<std::size_t>(out.w) * static_cast<std::size_t>(out.h) * 3, 255);

    // Pictures sit below a caption bar; kGap pads left, between, right, top and bottom.
    const int xl = kGap;
    const int xr = 2 * kGap + left.w;
    const int ty = kGap + kBarH;
    blit(out, left, xl, ty);
    blit(out, right, xr, ty);

    drawText(out, xl, kGap + 2 * kTextScale, leftCaption, 20);
    drawText(out, xr, kGap + 2 * kTextScale, rightCaption, 20);
    return out;
}


bool writeTGA(const std::string &path, const Image &img) {
    return writeTGA(path, img.w, img.h, img.rgb, /*bottomUp=*/false);
}


Image downscale(const Image &src, int factor) {
    Image out;
    out.w = src.w / factor;
    out.h = src.h / factor;
    out.rgb.resize(static_cast<std::size_t>(out.w) * static_cast<std::size_t>(out.h) * 3);
    const int area = factor * factor;

    // Box filter with round-to-nearest: (sum + area/2) / area.
    for (int y = 0; y < out.h; ++y) {
        for (int x = 0; x < out.w; ++x) {
            for (int c = 0; c < 3; ++c) {
                int sum = 0;
                for (int dy = 0; dy < factor; ++dy) {
                    const std::size_t row = static_cast<std::size_t>(y * factor + dy) * static_cast<std::size_t>(src.w);
                    for (int dx = 0; dx < factor; ++dx) {
                        sum += src.rgb[(row + static_cast<std::size_t>(x * factor + dx)) * 3 + static_cast<std::size_t>(c)];
                    }
                }
                out.rgb[(static_cast<std::size_t>(y) * static_cast<std::size_t>(out.w) + static_cast<std::size_t>(x)) * 3 +
                        static_cast<std::size_t>(c)] = static_cast<unsigned char>((sum + area / 2) / area);
            }
        }
    }
    return out;
}

} // namespace ui
