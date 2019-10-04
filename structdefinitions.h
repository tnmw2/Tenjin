#ifndef STRUCTDEFINITIONS
#define STRUCTDEFINITIONS

#include "amrexheader.h"

enum Variable
{
    RHO,
    RHOU,
    TOTAL_E,
    VELOCITY,
    P,
    SOUNDSPEED
};

struct InitialStruct
{
    Real startT;
    Real finalT;

    std::string filename;
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

#endif // STRUCTDEFINITIONS

