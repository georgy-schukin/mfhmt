TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

TARGET = mfhmt

SOURCES +=   \
    mfhmt.cpp \
    output.cpp

HEADERS +=   \
    array3d.h \
    common.h \
    output.h
