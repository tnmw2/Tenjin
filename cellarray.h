#ifndef CELLARRAY_H
#define CELLARRAY_H

#include "amrexheader.h"
#include "structdefinitions.h"

class CellArray;

class BoxAccessCellArray
{
public:

    BoxAccessCellArray(const Box& bx, FArrayBox& fb, Array4<Real> const& prop_arr);

    const Box& box;

    FArrayBox& fab;

    Array4<Real> const& arr;
};

class CellArray
{
public:

    CellArray(BoxArray& ba, DistributionMapping& dm, const int Ncomp, const int Nghost);

    void primitiveToConservative(ParameterStruct& parameters);
    void primitiveToConservative(BoxAccessCellArray& U, ParameterStruct& parameters);

    void conservativeToPrimitive(ParameterStruct& parameters);
    void conservativeToPrimitive(BoxAccessCellArray& U, ParameterStruct& parameters);

    void getSoundSpeed(ParameterStruct& parameters);
    void getSoundSpeed(BoxAccessCellArray& U, ParameterStruct& parameters);


    void operator=(CellArray& U);
    CellArray& operator*(Real d);
    CellArray& operator+(CellArray& U);

    MultiFab data;


};

#endif // CELLARRAY_H
