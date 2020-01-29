#ifndef SIMULATIONHEADER
#define SIMULATIONHEADER

#include <new>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
#include <iostream>
#include <sstream>
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
#include "Adv_F.H"


void advance(CellArray& U, CellArray& U1, CellArray &U2, CellArray &MUSCLgrad, CellArray &UL, CellArray &UR, CellArray& ULStar, CellArray &URStar, CellArray& UStarStar, Array<MultiFab, AMREX_SPACEDIM>& flux_arr, Geometry const& geom, ParameterStruct& parameters, Vector<BCRec> &bc, THINCArray& THINC);

void PrintAllVarsTo1DGnuplotFile(CellArray &U, int picture, std::__cxx11::string filename, const Real *dx, const Real *prob_lo);

void libConfigInitialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial, PlasticEOS& plastic);
void setInitialConditions(CellArray& U, ParameterStruct& parameters, InitialStruct &initial, const Real* dx, const Real* prob_lo);
void setBoundaryConditions(Vector<BCRec>& bc, ParameterStruct& parameters, InitialStruct& initial, AccessPattern& accessPattern);
void geometricSourceTerm(CellArray& U, ParameterStruct& parameters, const Real *dx, Real dt, const Real *prob_lo, MultiFab &S_new);

void reactiveUpdate(CellArray& U, CellArray& U1, CellArray& U2, ParameterStruct& parameters, Real dt, MultiFab &S_new);
void reactiveUpdateInHLLC(BoxAccessCellArray& U, ParameterStruct& parameters, Real dt);
void reactiveUpdateOutHLLC(CellArray& U, CellArray &U1, CellArray &U2, ParameterStruct& parameters, Real dt);


void customAbort(Vector<Real>& values, std::string& Message);
void getStarStateAlone(BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, ParameterStruct& parameters, Direction_enum d, const Real *dx);
void calculatePathConservativeFluxes(BoxAccessCellArray& fluxboxL, BoxAccessCellArray& fluxboxR, BoxAccessCellArray& Ubox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, ParameterStruct& parameters, Direction_enum d, const Real *dx);
void PCupdate(BoxAccessCellArray& fluxboxL, BoxAccessCellArray& fluxboxR, BoxAccessCellArray& Ubox, BoxAccessCellArray& U1box, ParameterStruct& parameters, Direction_enum d, Real dt, const Real* dx);

#endif // SIMULATIONHEADER

