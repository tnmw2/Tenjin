#ifndef CELLARRAY_H
#define CELLARRAY_H

#include "amrexheader.h"
#include "structdefinitions.h"

class CellArray;

class BoxAccessCellArray
{
public:

    BoxAccessCellArray(const Box& bx, FArrayBox& fb, Array4<Real> const& prop_arr, CellArray &U);

    const Box& box;

    FArrayBox& fab;

    Array4<Real> const& arr;

    std::map<Variable,int>& accessPattern;

    Real& operator()(Variable var, int mat, int row, int col, int i, int j=0, int k=0);
    Real& operator()(MaterialSpecifier& m, int i, int j=0, int k=0);
    Real& operator()(MaterialSpecifier  m, int i, int j=0, int k=0);
    Real& operator()(int i, int j, int k, Variable var, int mat=0, int row=0, int col=0);
    Real& operator()(int i, int j, int k, int var, int mat=0, int row=0, int col=0);



};

class CellArray
{
public:

    CellArray(BoxArray& ba, DistributionMapping& dm, const int Ncomp, const int Nghost, std::map<Variable,int>& _accessPattern);

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

    std::map<Variable,int>& accessPattern;



};

#endif // CELLARRAY_H
