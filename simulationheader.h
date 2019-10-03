#ifndef SIMULATIONHEADER
#define SIMULATIONHEADER


#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
#include <iostream>

#include "amrexheader.h"
#include "timestep.h"
#include "cellarray.h"



using namespace amrex;



struct InitialStruct
{
    Real startT;
    Real finalT;
};

struct ParameterStruct
{
    Vector<Real> dimL;
    Vector<Real> dx;
    Vector<int>  n_cells;

    int Ncomp;
    int Nghost;

    int max_grid_size;
    int plot_int;

    Real phiL;
    Real phiR;
    Real CFL;
    Real x0;
    Real a;
    Real adiabaticIndex;
    Real dt;

    ParameterStruct()
    {
        dimL 	= Vector<Real>(AMREX_SPACEDIM);
        dx	 	= Vector<Real>(AMREX_SPACEDIM);
        n_cells = Vector<int> (AMREX_SPACEDIM);
    }
};


void initialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial);


#endif // SIMULATIONHEADER

