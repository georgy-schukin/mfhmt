TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

TARGET = mfhmt_omp

SOURCES +=   \
    mfhmt_omp.cpp \
    output.cpp \
    timer.cpp

HEADERS +=   \
    array3d.h \
    common.h \
    output.h \
    timer.h

win32-msvc* {
    QMAKE_CXXFLAGS += /openmp
}

win32-g++*|linux* {
    QMAKE_CXXFLAGS += -fopenmp
    LIBS += -fopenmp
}
