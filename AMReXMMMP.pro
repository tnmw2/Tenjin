TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    timestep.cpp \
    cellarray.cpp \
    initialise.cpp

HEADERS += \
    simulationheader.h \
    timestep.h \
    cellarray.h \
    amrexheader.h

