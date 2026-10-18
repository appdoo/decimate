# check_golden_reproducible.cmake
# Called by CTest to verify the golden fixture is reproducible (T5).
# Runs capture_golden into a temp file and diffs against the checked-in golden.

# BUILD_DIR is passed in from CMakeLists; strip any surrounding quotes.
string(REPLACE "\"" "" BUILD_DIR_CLEAN "${BUILD_DIR}")
set(TEMP_FILE "${BUILD_DIR_CLEAN}/raytracing_golden_regen.txt")

execute_process(
    COMMAND "${CAPTURE_EXE}" "${TEMP_FILE}"
    RESULT_VARIABLE CAPTURE_RESULT
)
if(NOT CAPTURE_RESULT EQUAL 0)
    message(FATAL_ERROR "capture_golden exited with code ${CAPTURE_RESULT}")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files "${GOLDEN_FILE}" "${TEMP_FILE}"
    RESULT_VARIABLE DIFF_RESULT
)
if(NOT DIFF_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Golden fixture is NOT reproducible!\n"
        "  Expected: ${GOLDEN_FILE}\n"
        "  Got:      ${TEMP_FILE}\n"
        "The checked-in golden does not match a fresh capture. "
        "Re-run capture_golden and commit the new fixture."
    )
endif()

message(STATUS "Golden fixture is reproducible (byte-identical).")
