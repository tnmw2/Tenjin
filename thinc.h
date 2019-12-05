#ifndef THINC_H
#define THINC_H

#include "amrexheader.h"
#include "structdefinitions.h"


class THINCArray
{
public:

    THINCArray(BoxArray& ba, DistributionMapping& dm, const int Nghost, ParameterStruct &parameters);

    void addTHINCvariable(Variable var, int materialNumber=1, int rowNumber=1, int colNumber=1);

    iMultiFab data;

    Vector<MaterialSpecifier> THINCvariables;
};

class BoxAccessTHINCArray
{
public:

    BoxAccessTHINCArray(MFIter& mfi, const Box &bx, THINCArray &U);

    const Box& box;

    IArrayBox& iab;

    Vector<MaterialSpecifier>& THINCvariables;

    int& mixedCellFlag  (int i, int j, int k);
    int& TBVFlag        (int i, int j, int k, int m);

    void THINCreconstruction(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& UTHINC_L, BoxAccessCellArray& UTHINC_R, ParameterStruct& parameters, const Real *dx, Direction_enum d);
    void cautiousTHINCreconstruction(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& UTHINC_L, BoxAccessCellArray& UTHINC_R, ParameterStruct& parameters, const Real *dx, Direction_enum d);
    void variableTHINCreconstruction(MaterialSpecifier& n, BoxAccessCellArray& U, BoxAccessCellArray& UTHINC_L, BoxAccessCellArray& UTHINC_R, Vector<Real>& beta, Vector<Real>& coshBeta, Vector<Real>& tanhBeta, Real epsilon, int i, int j ,int k, Direction_enum d);


};

#endif // THINC_H

