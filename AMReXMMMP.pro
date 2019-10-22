TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    timestep.cpp \
    cellarray.cpp \
    initialise.cpp \
    boxaccesscellarray.cpp \
    hllc.cpp \
    print.cpp \
    cell.cpp \
    accesspattern.cpp \
    reactive.cpp \
    equationofstate.cpp \
    tensor.cpp \
    thinc.cpp \
    radial.cpp \
    flux.cpp \
    boundary.cpp

HEADERS += \
    simulationheader.h \
    timestep.h \
    cellarray.h \
    amrexheader.h \
    structdefinitions.h \
    cell.h \
    accesspattern.h \
    equationofstate.h \
    tensor.h \
    thinc.h \
    flux.h

