#ifndef SIMULATIONHEADER
#define SIMULATIONHEADER


#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <utility>
#include <map>
#include <vector>
#include <libconfig.h++>

#include "tensor.h"
#include "amrexheader.h"

using namespace amrex;


#include "structdefinitions.h"
#include "accesspattern.h"
#include "equationofstate.h"
#include "cellarray.h"
#include "cell.h"
#include "thinc.h"
#include "flux.h"
#include "plastic.h"




void advance(CellArray& U, CellArray& U1, CellArray &U2, CellArray &MUSCLgrad, CellArray &UL, CellArray &UR, CellArray& ULStar, CellArray &URStar, CellArray& UStarStar, Array<MultiFab, AMREX_SPACEDIM>& flux_arr, Geometry const& geom, ParameterStruct& parameters, Vector<BCRec> &bc, THINCArray& THINC);

void PrintAllVarsTo1DGnuplotFile(CellArray &U, int picture, std::__cxx11::string filename);

void libConfigInitialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial, PlasticEOS &plastic);
//void initialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial);
void setInitialConditions(CellArray& U, ParameterStruct& parameters, InitialStruct &initial);
void setBoundaryConditions(Vector<BCRec>& bc, ParameterStruct& parameters, InitialStruct& initial, AccessPattern &accessPattern);

void geometricSourceTerm(CellArray& U, ParameterStruct& parameters);

void reactiveUpdate(CellArray& U, CellArray& U1, CellArray& U2, ParameterStruct& parameters);
#endif // SIMULATIONHEADER

