cmake_minimum_required(VERSION 3.10)

project(MFHMT_F Fortran)

add_executable(mfhmt_f
    smt45.f
)

install(TARGETS mfhmt_f DESTINATION fortran)
