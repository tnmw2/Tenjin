#ifndef CELL_H
#define CELL_H

#include "amrexheader.h"
#include "structdefinitions.h"
#include "cellarray.h"

class Material
{
    public:


    Material(){}


    Material_type		phase;


    //Thermodynamic Data ----------------

    Real*	alpha;
    Real* 	rho;
    Real*	alphaRho;


};

class Cell
{
public:

    Real* rho;
    std::vector<Real*> rhoU;
    Real* E;
    std::vector<Real*> u;
    Real* p;
    Real* a;
    Real* uStar;

    static const int numberOfComponents = 3;

    std::vector<Material> materials;

    Cell(BoxAccessCellArray& U, int i, int j, int k);

    //Cell(Array4<Real> const& arr, int i, int j, int k);

    Real& operator()(Variable var, int mat=0, int row=0, int col=0);
    Real& operator()(MaterialSpecifier m);
};

#endif // CELL_H
