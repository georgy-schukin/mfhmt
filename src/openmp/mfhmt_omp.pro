TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

TARGET = mfhmt_omp

SOURCES +=   \
    ../common/output.cpp \
    mfhmt_omp.cpp

HEADERS +=   \
    ../common/output.h \
    ../common/timer.h \
    array3d.h \
    defs.h

win32-msvc* {
    QMAKE_CXXFLAGS += /openmp
}

win32-g++*|linux* {
    QMAKE_CXXFLAGS += -fopenmp
    LIBS += -fopenmp
}
