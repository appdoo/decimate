// This file is part of the decimate library; see decimate.h for version/license details.
/// @file model_fixture.h
/// @brief Test helper: resolve a model under MESH_MODEL_DIR, or SKIP the test if absent.
///
/// The large models are not part of a fresh clone: `horse.iv`/`bunny.iv` are fetched by
/// `tools/fetch_models.sh`, and `cow.obj` is a local-only fixture that is never
/// redistributed.  Tests that need them call `modelOrSkip("cow.obj")` instead of building
/// the path by hand, so a clone without the models reports "skipped", not "failed".
/// `spot.obj` ships in the repo and is loaded directly (its absence must FAIL).
#pragma once

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <string>

#ifndef MESH_MODEL_DIR
#error "MESH_MODEL_DIR must be defined by CMake"
#endif

/// @brief Absolute path of `MESH_MODEL_DIR/<name>`; SKIPs the running test case if the
///        file does not exist.  SKIP unwinds by exception, so this is safe to call from
///        helper functions reached from a TEST_CASE.
/// @param name  File name relative to the model directory (e.g. "cow.obj").
inline std::string modelOrSkip(const std::string &name) {
    const std::string path = std::string(MESH_MODEL_DIR) + "/" + name;
    if (std::filesystem::exists(path)) {
        return path;
    }
    if (name == "cow.obj") {
        SKIP("model/cow.obj is a local-only fixture (not redistributable)");
    }
    SKIP("model/" << name << " not present \xE2\x80\x94 run tools/fetch_models.sh");
    return path; // not reached: SKIP throws
}
