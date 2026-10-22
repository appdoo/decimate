// This file is part of the decimate library; see decimate.h for version/license details.
/// @file cover_segment_test.cpp
/// @brief Spec 19 Amendment 4: unit tests of the packer's conservative silhouette primitive.
///
/// S1 [S]  steep edge inside one column: every row it passes is covered (2003 took ONE sample)
/// S2 [S]  steep edge across a column boundary: both columns cover the shared boundary row
/// S3 [S]  end point exactly on a column boundary: the next column is covered too
/// S4 [S]  rows clamp to [0, res-1]; untouched columns keep their initial values

#include <catch2/catch_test_macros.hpp>

#include <vector>

#include <decimate/chart/parametrize.h>

using decimate::chart::coverSegment;

namespace {

/// Per-column row range, as the packer keeps it: `lower` starts at res (empty), `upper` at -1 (empty).
struct Silhouette {
    std::vector<int> lower;
    std::vector<int> upper;
    explicit Silhouette(int res)
        : lower(static_cast<std::size_t>(res), res), upper(static_cast<std::size_t>(res), -1) {
    }
};

} // namespace

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S1: a steep edge inside one column covers every row it passes", "[chart][packer][cover]") {
    Silhouette s(16);
    coverSegment(1.2f, 0.5f, 1.8f, 9.5f, 16, s.lower, s.upper);
    CHECK(s.lower[1] == 0);
    CHECK(s.upper[1] == 9);

    // Columns the edge never enters keep their empty sentinels.
    CHECK(s.upper[0] == -1); // neighbours untouched
    CHECK(s.upper[2] == -1);
    CHECK(s.lower[0] == 16);
    CHECK(s.lower[2] == 16);
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S2: a steep edge across a column boundary covers the boundary row in both columns",
          "[chart][packer][cover]") {
    Silhouette s(16);
    coverSegment(1.5f, 0.0f, 2.5f, 10.0f, 16, s.lower, s.upper); // slope 10, crosses x = 2 at v = 5
    CHECK(s.lower[1] == 0);
    CHECK(s.upper[1] == 5);
    CHECK(s.lower[2] == 5);
    CHECK(s.upper[2] == 10);

    CHECK(s.upper[0] == -1);
    CHECK(s.upper[3] == -1);
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S3: an end point on a column boundary covers the column that starts there",
          "[chart][packer][cover]") {
    Silhouette s(16);
    coverSegment(5.0f, 4.0f, 3.0f, 2.0f, 16, s.lower, s.upper); // given right-to-left on purpose
    CHECK(s.lower[3] == 2);
    CHECK(s.upper[3] == 3);
    CHECK(s.lower[4] == 3);
    CHECK(s.upper[4] == 4);

    // The end point on x = 5.0 touches column 5 only in a single point.
    CHECK(s.lower[5] == 4); // x = 5.0 is the left border of column 5: a single point, row 4
    CHECK(s.upper[5] == 4);
    CHECK(s.upper[2] == -1);
    CHECK(s.upper[6] == -1);
}

// NOLINTNEXTLINE(readability-function-size,readability-function-cognitive-complexity)
TEST_CASE("S4: rows clamp to the grid and accumulate as min/max over several segments",
          "[chart][packer][cover]") {
    Silhouette s(8);
    coverSegment(2.0f, 0.0f, 2.9f, 7.0f, 8, s.lower, s.upper);
    coverSegment(2.1f, 3.0f, 2.8f, 3.5f, 8, s.lower, s.upper); // inside the first: no change
    CHECK(s.lower[2] == 0);
    CHECK(s.upper[2] == 7);

    // A second segment at the grid edge: row 7 is the last valid row, no clamp overshoot.
    coverSegment(6.0f, 7.0f, 6.5f, 7.0f, 8, s.lower, s.upper); // top row, at the grid edge
    CHECK(s.lower[6] == 7);
    CHECK(s.upper[6] == 7);
}
