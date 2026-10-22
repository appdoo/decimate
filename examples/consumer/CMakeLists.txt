# Standalone project — built AGAINST an installed decimate, NOT part of the main
# build. Proves external consumption:
#   cmake --install <build> --prefix <P>
#   cmake -S examples/consumer -B <cb> -DCMAKE_PREFIX_PATH=<P> && cmake --build <cb>
cmake_minimum_required(VERSION 3.21)
project(decimate_consumer LANGUAGES CXX)

find_package(decimate 0.7 REQUIRED COMPONENTS chart bake mesh)

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE decimate::all)
