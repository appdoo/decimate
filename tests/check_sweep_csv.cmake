# check_sweep_csv.cmake
# Called by CTest (Spec 17-03, I3): runs seamless_sweep on one setting and checks the CSV
# artifact -- header byte-for-byte, exactly one row per stage, every numeric field finite.
#
#   -DSWEEP_EXE=<tool> -DMESH=<obj> -DOUT=<csv>

set(EXPECTED_HEADER "maxCharts,cutThreshold,texSize,space,stage,targetFraction,achievedVertices,earFaces,sameNeighbourEars,protectedVertices,seamEdges,seamLength3D,seamTexels,dR,dG,dB,collapsibleSeamEdges,uncollapsibleSeamEdges,collapsedSeamEdges,revertedCharts")
set(EXPECTED_STAGES chartify straighten dilate seamfill erased lod50 lod10)

file(REMOVE "${OUT}")
execute_process(
    COMMAND "${SWEEP_EXE}" --mesh "${MESH}" --maxcharts 64 --thresholds 0.3 --out "${OUT}"
    RESULT_VARIABLE SWEEP_RESULT
)
if(NOT SWEEP_RESULT EQUAL 0)
    message(FATAL_ERROR "seamless_sweep exited with code ${SWEEP_RESULT}")
endif()

file(STRINGS "${OUT}" LINES)
list(LENGTH LINES NUM_LINES)
list(GET LINES 0 HEADER)
if(NOT HEADER STREQUAL EXPECTED_HEADER)
    message(FATAL_ERROR "CSV header differs:\n  got:      ${HEADER}\n  expected: ${EXPECTED_HEADER}")
endif()
math(EXPR NUM_ROWS "${NUM_LINES} - 1")
if(NOT NUM_ROWS EQUAL 7)
    message(FATAL_ERROR "expected 7 rows (one per stage), got ${NUM_ROWS}")
endif()

set(NUMBER_RE "^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$")
set(ROW_INDEX 0)
foreach(LINE IN LISTS LINES)
    if(ROW_INDEX GREATER 0)
        string(REPLACE "," ";" FIELDS "${LINE}")
        list(LENGTH FIELDS NUM_FIELDS)
        if(NOT NUM_FIELDS EQUAL 20)
            message(FATAL_ERROR "row ${ROW_INDEX} has ${NUM_FIELDS} fields: ${LINE}")
        endif()
        math(EXPR STAGE_INDEX "${ROW_INDEX} - 1")
        list(GET EXPECTED_STAGES ${STAGE_INDEX} WANT_STAGE)
        list(GET FIELDS 4 GOT_STAGE)
        if(NOT GOT_STAGE STREQUAL WANT_STAGE)
            message(FATAL_ERROR "row ${ROW_INDEX}: stage '${GOT_STAGE}', expected '${WANT_STAGE}'")
        endif()
        foreach(COLUMN 0 1 2 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19)
            list(GET FIELDS ${COLUMN} VALUE)
            if(NOT VALUE STREQUAL "" AND NOT VALUE MATCHES "${NUMBER_RE}")
                message(FATAL_ERROR "row ${ROW_INDEX} column ${COLUMN}: '${VALUE}' is not a finite number")
            endif()
        endforeach()
    endif()
    math(EXPR ROW_INDEX "${ROW_INDEX} + 1")
endforeach()

message(STATUS "sweep CSV ok: ${NUM_ROWS} rows, header exact, all numeric fields finite.")
