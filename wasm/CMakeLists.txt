# Wasm build of the unchanged core (spec 34). Included from the root only under Emscripten.
# Output: decimate.js (ES6 module factory) + decimate.wasm.
add_executable(decimate_wasm bindings.cpp)
target_link_libraries(decimate_wasm PRIVATE decimate::pipeline decimate::io decimate::reduce)
set_target_properties(decimate_wasm PROPERTIES OUTPUT_NAME decimate SUFFIX ".js")
target_link_options(decimate_wasm PRIVATE
    -lembind
    -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createDecimate
    -sALLOW_MEMORY_GROWTH=1 -sMAXIMUM_MEMORY=4GB -sENVIRONMENT=web,worker,node
    -sFORCE_FILESYSTEM=1
    # Native-like stack: the 64 KB default overflows on Armadillo and, without assertions, silently
    # corrupts the heap. 4 GB: the tangent erase at 1024 needs ~2.8 GB.
    -sSTACK_SIZE=8MB)

# Copy the artefacts next to the static site (web/), if it exists.
if(EXISTS "${CMAKE_SOURCE_DIR}/web")
    add_custom_target(web
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE_DIR:decimate_wasm>/decimate.js
                $<TARGET_FILE_DIR:decimate_wasm>/decimate.wasm "${CMAKE_SOURCE_DIR}/web/"
        DEPENDS decimate_wasm)
endif()

# W1 + W2: one ctest entry, skipped (not failed) when node or the native reference is missing.
find_program(NODE_EXECUTABLE node NO_CMAKE_FIND_ROOT_PATH)
# W2 needs a native reference that rounds like Wasm: wasm has no FMA, arm64 clang contracts a*b+c by
# default, and that alone makes the QEM collapse order (hence the whole atlas) diverge. Build it with
#   cmake -S . -B build-w2ref -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS=-ffp-contract=off \
#         -DCMAKE_CXX_FLAGS=-ffp-contract=off -DMESH_BUILD_UI=OFF && cmake --build build-w2ref --target wasm_ref
# With that, wasm and native are bit-identical (max map diff 0). Missing -> W2 is skipped (exit 77).
set(DECIMATE_NATIVE_REF "${CMAKE_SOURCE_DIR}/build-w2ref/wasm/ref/wasm_ref" CACHE FILEPATH "Host wasm_ref binary for W2")
if(NODE_EXECUTABLE)
    add_test(NAME wasm_w1_w2
        COMMAND ${NODE_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/w1_w2.mjs
                $<TARGET_FILE_DIR:decimate_wasm>/decimate.js ${CMAKE_SOURCE_DIR}/model/spot.obj ${DECIMATE_NATIVE_REF})
    set_tests_properties(wasm_w1_w2 PROPERTIES SKIP_RETURN_CODE 77)
endif()
