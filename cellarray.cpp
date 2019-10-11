#include "cellarray.h"

CellArray::CellArray(BoxArray& ba, DistributionMapping& dm, const int Ncomp, const int Nghost, AccessPattern &_accessPattern, ParameterStruct& parameters) : data(ba,dm,Ncomp,Nghost), accessPattern(_accessPattern), numberOfMaterials{parameters.numberOfMaterials}{}

void CellArray::conservativeToPrimitive(ParameterStruct& parameters)
{
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,(*this));

        baca.conservativeToPrimitive(parameters);
    }
}

void CellArray::primitiveToConservative(ParameterStruct& parameters)
{
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,(*this));

        baca.primitiveToConservative(parameters);
    }
}

void CellArray::getSoundSpeed(ParameterStruct& parameters)
{
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,(*this));

        baca.getSoundSpeed(parameters);
    }
}

void CellArray::operator=(CellArray& U)
{
    MultiFab::Copy(data, U.data, 0, 0, data.nComp(), data.nGrow());

    return;
}

CellArray& CellArray::operator*(Real d)
{
    data.mult(d, 0, data.nComp(), data.nGrow());

    return *this;
}

CellArray& CellArray::operator+(CellArray& U)
{
    MultiFab::Add(data, U.data, 0, 0, data.nComp(), data.nGrow());

    return *this;
}



