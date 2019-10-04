#ifndef SIMULATIONHEADER
#define SIMULATIONHEADER


#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
#include <iostream>
#include <ctime>

#include "amrexheader.h"

using namespace amrex;


#include "structdefinitions.h"
#include "cellarray.h"

void HLLCadvance(CellArray& U,CellArray& U1,CellArray& UStar,Array<MultiFab, AMREX_SPACEDIM>& flux_arr,Geometry const& geom, ParameterStruct& parameters);

void PrintAllVarsTo1DGnuplotFile(MultiFab& phi, int picture, std::__cxx11::string filename);



void initialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial);
void setInitialConditions(CellArray& U, ParameterStruct& parameters);


#endif // SIMULATIONHEADER

