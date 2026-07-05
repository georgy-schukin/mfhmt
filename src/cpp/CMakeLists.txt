cmake_minimum_required(VERSION 3.0...3.5)

project(MFHMT_CPP CXX)

add_executable(mfhmt
    mfhmt.cpp
    output.cpp
)

install(TARGETS mfhmt DESTINATION cpp)
