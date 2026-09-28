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
    block_decomp.cpp \
    comm_util.cpp \
    distributed_array3d.cpp \
    mfhmt_mpi.cpp

HEADERS += \
    ../common/output.h \
    ../common/timer.h \
    block_decomp.h \
    comm_util.h \
    defs.h \
    distributed_array3d.h \
    shadowed_array3d.h