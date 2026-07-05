cmake_minimum_required(VERSION 3.0...3.5)

project(MFHMT LANGUAGES CXX Fortran)

add_subdirectory(fortran)
add_subdirectory(cpp)
add_subdirectory(openmp)