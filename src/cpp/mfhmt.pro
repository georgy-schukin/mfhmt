TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

TARGET = mfhmt

SOURCES +=   \
    ../common/output.cpp \
    mfhmt.cpp

HEADERS +=   \
    ../common/output.h \
    array3d.h \
    defs.h
