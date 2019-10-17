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


};

class Cell
{
public:

    Real* rho;
    std::vector<Real*> rhoU;
    Real* E;
    std::vector<Real*> u;
    Real* p;
    Vector<Real*> sigma;

    Real* a;
    Real* uStar;

    static const int numberOfComponents = 3;

    std::vector<Material> materials;

    Cell(BoxAccessCellArray& U, int i, int j, int k);
    Cell(){}

    Real& operator()(Variable var, int mat=0, int row=0, int col=0);
    Real& operator()(MaterialSpecifier m);
};

#endif // CELL_H
