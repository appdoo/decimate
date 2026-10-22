// This file is part of the decimate library; see decimate.h for version/license details.
/// @file puremath_vec3_test.cpp
/// @brief Exact-value tests of decimate::math::Vec2/3/4: cross/dot/length, union aliasing, converting ctor.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <type_traits>

#include <decimate/math/vec2.h>
#include <decimate/math/vec3.h>
#include <decimate/math/vec4.h>

TEST_CASE("decimate::math::Vec3 cross product", "[puremath][vec3]") {
    // T4a: (1,0,0) x (0,1,0) == (0,0,1)
    decimate::math::Vec3<double> result = decimate::math::Vec3<double> { 1, 0, 0 }.cross(decimate::math::Vec3<double> { 0, 1, 0 });
    CHECK(result.x == 0.0);
    CHECK(result.y == 0.0);
    CHECK(result.z == 1.0);
}


TEST_CASE("decimate::math::Vec3 dot product", "[puremath][vec3]") {
    // T4b: (1,0,0) . (0,1,0) == 0
    double result = decimate::math::Vec3<double> { 1, 0, 0 }.dot(decimate::math::Vec3<double> { 0, 1, 0 });
    CHECK(result == 0.0);
}


TEST_CASE("decimate::math::Vec3 length", "[puremath][vec3]") {
    // T4c: (3,4,0).length() == 5 (3-4-5 right triangle)
    float result = decimate::math::Vec3f { 3, 4, 0 }.length();
    CHECK_THAT(result, Catch::Matchers::WithinAbs(5.0f, 1e-5f));
}


TEST_CASE("decimate::math::Vec3 union: .x/.y/.z and data[i]/operator[] alias same storage",
          "[puremath][vec3][union]") {
    decimate::math::Vec3f v { 1.0f, 2.0f, 3.0f };

    // Named access and indexed access read the same bytes.
    CHECK(v[0] == 1.0f);
    CHECK(v[1] == 2.0f);
    CHECK(v[2] == 3.0f);
    CHECK(v.data[0] == 1.0f);
    CHECK(v.data[2] == 3.0f);

    // Writing through one view is visible through the other.
    v[1] = 9.0f;
    CHECK(v.y == 9.0f);
    v.z = 7.0f;
    CHECK(v[2] == 7.0f);

    // Layout contract the reinterpret_cast in MeshData/bake relies on.
    static_assert(std::is_standard_layout<decimate::math::Vec3f>::value);
    static_assert(sizeof(decimate::math::Vec3f) == 3 * sizeof(float));

    // An array of Vec3 is contiguous floats: arr[1].x follows arr[0].z.
    decimate::math::Vec3f arr[2] = { { 1, 2, 3 }, { 4, 5, 6 } };
    const float          *flat = reinterpret_cast<const float *>(arr);
    CHECK(flat[3] == 4.0f); // arr[1].x
    CHECK(flat[5] == 6.0f); // arr[1].z
}


TEST_CASE("decimate::math::Vec2/3/4 explicit cross-precision converting ctor",
          "[puremath][vec][convert]") {
    // Vec3<double> -> Vec3f: components static_cast'd, value preserved.
    decimate::math::Vec3<double> d3 { 1.5, -2.25, 3.0 };
    decimate::math::Vec3f        f3 { d3 };
    CHECK(f3.x == 1.5f);
    CHECK(f3.y == -2.25f);
    CHECK(f3.z == 3.0f);

    // Vec2 and Vec4 carry the same converting ctor.
    decimate::math::Vec2f f2 { decimate::math::Vec2<double> { 0.5, 4.0 } };
    CHECK(f2.x == 0.5f);
    CHECK(f2.y == 4.0f);

    decimate::math::Vec4f f4 { decimate::math::Vec4<double> { 1.0, 2.0, 3.0, 4.0 } };
    CHECK(f4.w == 4.0f);

    // The conversion is explicit — no silent Vec3<double> -> Vec3f.
    static_assert(!std::is_convertible_v<decimate::math::Vec3<double>, decimate::math::Vec3f>);
    static_assert(std::is_constructible_v<decimate::math::Vec3f, decimate::math::Vec3<double>>);
}
