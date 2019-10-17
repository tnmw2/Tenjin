#ifndef SIMULATIONHEADER
#define SIMULATIONHEADER


#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
#include <iostream>
#include <ctime>
#include <utility>
#include <map>
#include <vector>

#include "tensor.h"


#include "amrexheader.h"

using namespace amrex;


#include "structdefinitions.h"
#include "accesspattern.h"
#include "equationofstate.h"
#include "cellarray.h"
#include "cell.h"



void advance(CellArray& U, CellArray& U1, CellArray &U2, CellArray &MUSCLgrad, CellArray &UL, CellArray &UR, CellArray& UStar, Array<MultiFab, AMREX_SPACEDIM>& flux_arr, Geometry const& geom, ParameterStruct& parameters, Vector<BCRec> &bc);

void PrintAllVarsTo1DGnuplotFile(CellArray &U, int picture, std::__cxx11::string filename);

void initialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial);
void setInitialConditions(CellArray& U, ParameterStruct& parameters, InitialStruct &initial);

void reactiveUpdate(CellArray& U, CellArray& U1, CellArray& U2, ParameterStruct& parameters);
#endif // SIMULATIONHEADER

