# W3 + W4 (spec 34): export.js against the glTF validator and three.js OBJLoader. Included under Emscripten only,
# next to wasm_w1_w2. Needs `npm ci` in web/test (node_modules is not committed); without it the test exits 77 = skipped.
find_program(NODE_EXECUTABLE node NO_CMAKE_FIND_ROOT_PATH)
if(NODE_EXECUTABLE)
    add_test(NAME web_w3_w4
        COMMAND ${NODE_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/w3_w4.mjs
                $<TARGET_FILE_DIR:decimate_wasm>/decimate.js ${CMAKE_SOURCE_DIR}/model/spot.obj)
    set_tests_properties(web_w3_w4 PROPERTIES SKIP_RETURN_CODE 77)
endif()
