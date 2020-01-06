#ifndef CELL_H
#define CELL_H

#include "amrexheader.h"
#include "structdefinitions.h"
#include "cellarray.h"

/** \class Cell
 * Because AMReX stores things as n,i,j,k, rather than i,j,k,n
 * if we want to act on variables with the same index, then we
 * first group them into a Cell class for ease.
 *
 * \class Material
 * A sub-set of a cell which only contains information relating
 * to one material
 */

class Material
{
    public:


    Material(){}


    Material_type		phase;


    //Thermodynamic Data ----------------

    Real*	alpha;
    Real* 	rho;
    Real*	alphaRho;
    Real*   lambda;
    Real*   alphaRhoLambda;
    Real*   epsilon;
    Real*   alphaRhoEpsilon;


};

class Cell
{
public:

    Material_type phase;

    static const int numberOfComponents = 3;

    AccessPattern& accessPattern;

    Cell(BoxAccessCellArray& U, int i, int j, int k,Material_type _phase);
    Cell(const BoxAccessCellArray& U, int i, int j, int k,Material_type _phase);


    Real&  operator()(Variable var, int mat=0, int row=0, int col=0);
    Real&  operator()(MaterialSpecifier m);

    void  operator= (Cell& U);

    bool contains_nan();
    bool check(MaterialSpecifier &n);

    BoxAccessCellArray* parent;

    int parent_i;
    int parent_j;
    int parent_k;
};

#endif // CELL_H
