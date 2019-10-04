TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    timestep.cpp \
    cellarray.cpp \
    initialise.cpp \
    fluxarray.cpp \
    boxaccesscellarray.cpp \
    hllc.cpp \
    print.cpp

HEADERS += \
    simulationheader.h \
    timestep.h \
    cellarray.h \
    amrexheader.h \
    fluxarray.h \
    structdefinitions.h

