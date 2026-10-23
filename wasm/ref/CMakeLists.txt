# Host-side reference for W2: the same pipeline_run.h the Wasm module uses, as a native binary.
add_executable(wasm_ref ref.cpp)
target_include_directories(wasm_ref PRIVATE ..)
target_link_libraries(wasm_ref PRIVATE decimate::pipeline decimate::io decimate::reduce)
