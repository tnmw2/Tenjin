#ifndef THINC_H
#define THINC_H

#include "amrexheader.h"
#include "structdefinitions.h"


class THINCArray
{
public:

    THINCArray(BoxArray& ba, DistributionMapping& dm, const int Nghost, ParameterStruct &parameters);

    iMultiFab data;
};

class BoxAccessTHINCArray
{
public:

    BoxAccessTHINCArray(MFIter& mfi, const Box &bx, THINCArray &U);

    const Box& box;

    IArrayBox& iab;

    int& mixedCellFlag  (int i, int j, int k);
    int& TBVFlag        (int i, int j, int k, int m);

    void THINCreconstruction(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& UTHINC_L, BoxAccessCellArray& UTHINC_R, ParameterStruct& parameters, const Real *dx, Direction_enum d);

};

#endif // THINC_H

