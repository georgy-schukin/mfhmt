cmake_minimum_required(VERSION 3.10)

project(MFHMT_OMP CXX)

find_package(OpenMP REQUIRED)

add_executable(mfhmt_omp
    mfhmt_omp.cpp
    ../common/output.cpp
)

target_link_libraries(mfhmt_omp PRIVATE OpenMP::OpenMP_CXX)

install(TARGETS mfhmt_omp DESTINATION openmp)