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

enum Interface_type
{
    DIFFUSE,
    SHARP
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
    int      com;
    int      row;
    int      col;

    MaterialSpecifier(Variable v=RHO, int m=0, int p = 0, int r=0, int c=0)
    {
        var=v;
        mat=m;
        com=p;
        row=r;
        col=c;
    }

    bool operator <  (const MaterialSpecifier& b) const;
    bool operator == (const MaterialSpecifier& b) const;
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
    Interface_type  interface;
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
    int SOLID;
    int THINC;
    int RADIAL;
    int PLASTIC;
    int REACTIVE;
    int MUSCL;

    int numberOfMaterials;

    int Ncomp;

    int         numberOfSharpMaterials;
    Vector<int> numberOfDiffuseMaterials;
    Vector<int> diffuseMaterialContainsMixture;
    Vector<int> interfaceType;

    Real THINCbeta;

    Vector< Vector<MaterialDescriptor> >  materialInfo;

    ParameterStruct(){}
};

/** Holds data about initial conditions etc.
 */
struct InitialStruct
{
    int numberOfStates;
    int numberOfPictures;

    Real startT;
    Real finalT;

    Vector<Vector<Vector<Real> > > u;
    Vector<Vector<Vector<Real> > > v;
    Vector<Vector<Vector<Real> > > w;
    Vector<Vector<Vector<Real> > > p;
    Vector<Vector<Vector<Real> > > F;
    Vector<Vector<Vector<Real> > > rho;
    Vector<Vector<Vector<Real> > > rhoa;
    Vector<Vector<Vector<Real> > > rhob;
    Vector<Vector<Vector<Real> > > alpha;
    Vector<Vector<Vector<Real> > > lambda;

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
        rhoa.resize(numberOfStates);
        rhob.resize(numberOfStates);
        alpha.resize(numberOfStates);
        lambda.resize(numberOfStates);


        for(int i=0;i<numberOfStates;i++)
        {
            F     [i].resize(parameters.numberOfSharpMaterials);
            u     [i].resize(parameters.numberOfSharpMaterials);
            v     [i].resize(parameters.numberOfSharpMaterials);
            w     [i].resize(parameters.numberOfSharpMaterials);
            p     [i].resize(parameters.numberOfSharpMaterials);
            rho   [i].resize(parameters.numberOfSharpMaterials);
            rhoa  [i].resize(parameters.numberOfSharpMaterials);
            rhob  [i].resize(parameters.numberOfSharpMaterials);
            alpha [i].resize(parameters.numberOfSharpMaterials);
            lambda[i].resize(parameters.numberOfSharpMaterials);

            for(int s = 0; s < parameters.numberOfSharpMaterials; s++)
            {
                F     [i][s].resize(9);
                u     [i][s].resize(1);
                v     [i][s].resize(1);
                w     [i][s].resize(1);
                p     [i][s].resize(1);
                rho   [i][s].resize(parameters.numberOfDiffuseMaterials[s]);
                rhoa  [i][s].resize(parameters.numberOfDiffuseMaterials[s]);
                rhob  [i][s].resize(parameters.numberOfDiffuseMaterials[s]);
                alpha [i][s].resize(parameters.numberOfDiffuseMaterials[s]);
                lambda[i][s].resize(parameters.numberOfDiffuseMaterials[s]);
            }
        }



    }
};



#endif // STRUCTDEFINITIONS

