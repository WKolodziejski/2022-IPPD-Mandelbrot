cmake_minimum_required(VERSION 3.9)
project(2022_IPPD_Mandelbrot)

set(CMAKE_CXX_STANDARD 17)

add_executable(2022_IPPD_Mandelbrot main.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(2022_IPPD_Mandelbrot PUBLIC OpenMP::OpenMP_CXX)
endif()