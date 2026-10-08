cmake_minimum_required(VERSION 3.10)

project(MFHMT_MPI CXX)

find_package(MPI REQUIRED COMPONENTS CXX)

add_executable(mfhmt_mpi
    block_decomp.cpp
    distributed_array3d.cpp
    comm_util.cpp
    mfhmt_mpi.cpp
    ../common/output.cpp
)

target_link_libraries(mfhmt_mpi PRIVATE MPI::MPI_CXX)

install(TARGETS mfhmt_mpi DESTINATION mpi)
