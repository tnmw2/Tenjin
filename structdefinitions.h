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
    USTAR,
    RHO_MIX,
    ALPHARHOLAMBDA,
    LAMBDA,
    SIGMA
};

/** Not yet used...
 */
enum Material_type
{
    solid,
    fluid
};

enum Direction_enum
{
    x,
    y,
    z
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

enum Var_type
{
    PRIMITIVE,
    CONSERVATIVE,
    BOTH,
    NEITHER
};

class MieGruneisenEOS;

struct MaterialDescriptor
{
    Material_type   phase;
    bool            mixture;
    int             mixtureIndex;

    MieGruneisenEOS* EOS;
};


/** Holds data about initial conditions etc.
 */
struct InitialStruct
{
    int numberOfStates;

    Real startT;
    Real finalT;

    Vector<Real> rho;
    Vector<Real> u;
    Vector<Real> v;
    Vector<Real> w;
    Vector<Real> p;
    Vector<int> alpha;
    Vector<int> lambda;
    Vector<Real> interfaces;

    std::string filename;

    InitialStruct(){}
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
    int numberOfMixtures;

    int max_grid_size;

    Real CFL;
    Real x0;
    Real dt;

    Vector<Real> adiabaticIndex;
    Vector<Real> CV;

    Vector<Real> mixtureAdiabaticIndex;
    Vector<Real> mixtureCV;

    Vector<MaterialDescriptor>  materialInfo;

    ParameterStruct()
    {
        dimL 	= Vector<Real>(AMREX_SPACEDIM);
        dx	 	= Vector<Real>(AMREX_SPACEDIM);
        n_cells = Vector<int> (AMREX_SPACEDIM);
    }
};




#endif // STRUCTDEFINITIONS

