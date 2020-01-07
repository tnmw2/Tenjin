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

class Cell;

class Material
{
    public:


    Material(){}

    void allocateSpace(Cell *ptr, int n);
    void operator= (Material& M);

    Material_type		phase;


    //Thermodynamic Data ----------------

    Real*           rho;
    Vector<Real*>   rhoU;
    Vector<Real*>   u;
    Real*           E;
    Real*           p;
    Vector<Real*>   sigma;
    Vector<Real*>   V;
    Vector<Real*>   VStar;

    Real*           epsilon;
    Real*           rhoEpsilon;

    Real*           a;
    Real*           uStar;

    Cell*           parent;

    int             materialNumber;
};

class Cell
{
public:

    Material_type phase;

    static const int numberOfComponents = 3;

    //std::vector<Material> materials;

    AccessPattern& accessPattern;

    Cell(BoxAccessCellArray& U, int i, int j, int k,Material_type _phase);

    Real&  operator()(Variable var, int mat, int row=0, int col=0);
    Real&  operator()(MaterialSpecifier m);

    void  assignPointer(BoxAccessCellArray& U, int i, int j, int k, MaterialSpecifier m);

    void  operator= (Cell& U);

    bool contains_nan();
    bool check(MaterialSpecifier &n);

    void setMaterial(Cell& U, int m);

    BoxAccessCellArray* parent;

    int parent_i;
    int parent_j;
    int parent_k;
};

#endif // CELL_H
