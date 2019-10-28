#include "cellarray.h"

CellArray::CellArray(BoxArray& ba, DistributionMapping& dm, const int Ncomp, const int Nghost, AccessPattern &_accessPattern, ParameterStruct& parameters) : data(ba,dm,Ncomp,Nghost), accessPattern(_accessPattern), numberOfMaterials{parameters.numberOfMaterials}{}

void CellArray::conservativeToPrimitive()
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,(*this));

        baca.conservativeToPrimitive();
    }
}

void CellArray::primitiveToConservative()
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,(*this));

        baca.primitiveToConservative();
    }
}

void CellArray::getSoundSpeed()
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,(*this));

        baca.getSoundSpeed();
    }
}

void CellArray::cleanUpV()
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,(*this));

        baca.cleanUpV();
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

bool CellArray::contains_nan()
{
    bool checker = false;
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,(*this));

        checker = baca.contains_nan();
    }

    return checker;
}



