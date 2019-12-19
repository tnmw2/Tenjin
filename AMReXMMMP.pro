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
    boundary.cpp \
    flux.cpp \
    radial.cpp \
    plastic.cpp \
    AmrLevelAdv.cpp \
    newAMRFunctions.cpp \
    LevelBldAdv.cpp \
    pathconservativehllc.cpp

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
    flux.h \
    plastic.h \
    Adv_F.H \
    AmrLevelAdv.H

DISTFILES += \
    Make.package \
    ../build-AMReXMMMP-Desktop-Debug/GNUmakefile \
    bc_fill_nd.f90

