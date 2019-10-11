#ifndef STRUCTDEFINITIONS
#define STRUCTDEFINITIONS

#include "amrexheader.h"

/** An enumeration of all the Variables used. This is passed
 * to things like AccessPattern to return variables from a
 * BoxAccessCellArray
 */
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
    SOUNDSPEED,
    USTAR
};

/** Not yet used...
 */
enum Material_type
{
    solid,
    fluid
};

/** Holds data about initial conditions etc.
 */
struct InitialStruct
{
    Real startT;
    Real finalT;

    Real rhoL;
    Real rhoR;
    Real uL;
    Real uR;
    Real pL;
    Real pR;

    std::string filename;
};

/** Holds simulation data that needs to be passed around.
 *  Contains the current timestep dt and the CFL number
 *  for example.
 */
struct ParameterStruct
{
    Vector<Real> dimL;
    Vector<Real> dx;
    Vector<int>  n_cells;

    int Ncomp;
    int Nghost;
    int numberOfMaterials;

    int max_grid_size;

    Real CFL;
    Real x0;
    Real dt;

    Vector<Real> adiabaticIndex;

    ParameterStruct()
    {
        dimL 	= Vector<Real>(AMREX_SPACEDIM);
        dx	 	= Vector<Real>(AMREX_SPACEDIM);
        n_cells = Vector<int> (AMREX_SPACEDIM);
    }
};

/** A struct that can uniquely specify a thermodynamic varible.
 *  mat - the material the variable relates to
 *  row, col - used to get components of tensor variables like velocity
 *             and (later) deformation tensor.
 */
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

