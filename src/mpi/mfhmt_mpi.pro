TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX = mpicxx
QMAKE_CC = mpicc
QMAKE_LINK = mpicxx

TARGET = mfhmt_mpi

exists(local.pri) {
    include(local.pri)
}

SOURCES +=   \
    ../common/output.cpp \
    async_op.cpp \
    block_decomp.cpp \
    comm_util.cpp \
    distributed_array3d.cpp \
    mfhmt_mpi.cpp \
    reducer.cpp

HEADERS += \
    ../common/output.h \
    ../common/timer.h \
    async_op.h \
    block_decomp.h \
    cart_topology.h \
    comm_util.h \
    defs.h \
    distributed_array3d.h \
    reducer.h \
    shadowed_array3d.h \
    trace_event.h