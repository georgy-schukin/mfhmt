cmake_minimum_required(VERSION 3.10)

project(MFHMT LANGUAGES CXX Fortran)

add_subdirectory(fortran)
add_subdirectory(cpp)
add_subdirectory(openmp)
add_subdirectory(mpi)