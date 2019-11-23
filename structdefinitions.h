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
    SIGMA,
    V_TENSOR,
    VSTAR,
    DEVH,
    HJ2,
    EPSILON,
    ALPHARHOEPSILON
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
    z,
    LEFT,
    RIGHT
};

enum Boundary_type
{
    TRANMISSIVE,
    REFLECTIVE
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
    NEITHER,
    CELL,
    NOTCELL,
    REFINE
};

class MieGruneisenEOS;

struct MaterialDescriptor
{
    Material_type   phase;
    bool            mixture = 0;
    int             mixtureIndex;
    int             plastic = false;

    MieGruneisenEOS* EOS;
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

    int SOLID;
    int THINC;
    int RADIAL;
    int PLASTIC;
    int REACTIVE;
    int MUSCL;

    int Ncomp;
    int Nghost;
    int numberOfMaterials;
    int numberOfMixtures;
    int NLevelSets;

    int max_grid_size;

    Real CFL;
    Real x0;
    Real dt;
    Real THINCbeta;

    /*Vector<Real> adiabaticIndex;
    Vector<Real> CV;

    Vector<Real> mixtureAdiabaticIndex;
    Vector<Real> mixtureCV;*/

    Vector<MaterialDescriptor>  materialInfo;

    ParameterStruct()
    {
        dimL 	= Vector<Real>(AMREX_SPACEDIM);
        dx	 	= Vector<Real>(AMREX_SPACEDIM);
        n_cells = Vector<int> (AMREX_SPACEDIM);
    }
};

/** Holds data about initial conditions etc.
 */
struct InitialStruct
{
    int numberOfStates;
    int numberOfPictures;

    Real startT;
    Real finalT;

    Vector<Real> u;
    Vector<Real> v;
    Vector<Real> w;
    Vector<Real> p;
    Vector<Vector<Real> > F;
    Vector<Vector<Real> > rho;
    Vector<Vector<Real> > alpha;
    Vector<Vector<Real> > lambda;
    Real interface;

    Vector<int> lowBoundary;
    Vector<int> highBoundary;

    std::string filename;

    InitialStruct()
    {
        lowBoundary.resize(AMREX_SPACEDIM);
        highBoundary.resize(AMREX_SPACEDIM);
    }

    void resize(ParameterStruct& parameters)
    {
        u.resize(numberOfStates);
        v.resize(numberOfStates);
        w.resize(numberOfStates);
        p.resize(numberOfStates);
        F.resize(numberOfStates);
        rho.resize(numberOfStates);
        alpha.resize(numberOfStates);
        lambda.resize(numberOfStates);

        for(int i=0;i<numberOfStates;i++)
        {
            F[i].resize(9);
            rho[i].resize(parameters.numberOfMaterials);
            alpha[i].resize(parameters.numberOfMaterials);
            lambda[i].resize(parameters.numberOfMaterials);
        }



    }
};



#endif // STRUCTDEFINITIONS

