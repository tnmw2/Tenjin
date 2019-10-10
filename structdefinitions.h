#ifndef STRUCTDEFINITIONS
#define STRUCTDEFINITIONS

#include "amrexheader.h"

enum Variable
{
    ALPHA,
    ALPHARHO,
    RHO_K,
    RHO,
    RHOU,
    TOTAL_E,
    VELOCITY,
    P,
    SOUNDSPEED
};

enum Material_type
{
    solid,
    fluid
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
    int numberOfMaterials;

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

struct MaterialSpecifier
{
    Variable var;
    int      mat;
    int      row;
    int      col;

    MaterialSpecifier(Variable v=RHO, int m=0, int r=0, int c=0)
    {
        var=v;
        mat=m;
        row=r;
        col=c;
    }
};

#endif // STRUCTDEFINITIONS

