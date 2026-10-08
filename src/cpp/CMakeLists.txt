cmake_minimum_required(VERSION 3.10)

project(MFHMT_CPP CXX)

add_executable(mfhmt
    mfhmt.cpp
    ../common/output.cpp
)

install(TARGETS mfhmt DESTINATION cpp)
