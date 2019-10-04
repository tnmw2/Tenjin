#include "cellarray.h"

CellArray::CellArray(BoxArray& ba, DistributionMapping& dm, const int Ncomp, const int Nghost) : data(ba,dm,Ncomp,Nghost){}

void CellArray::conservativeToPrimitive(ParameterStruct& parameters)
{
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        FArrayBox& fab = data[mfi];

        Array4<Real> const& prop_arr = fab.array();

        BoxAccessCellArray baca(bx,fab,prop_arr);

        conservativeToPrimitive(baca,parameters);
    }
}

void CellArray::primitiveToConservative(ParameterStruct& parameters)
{
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        FArrayBox& fab = data[mfi];

        Array4<Real> const& prop_arr = fab.array();

        BoxAccessCellArray baca(bx,fab,prop_arr);

        primitiveToConservative(baca,parameters);
    }
}

void CellArray::conservativeToPrimitive(BoxAccessCellArray& U, ParameterStruct& parameters)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                U.arr(i,j,k,VELOCITY) = U.arr(i,j,k,RHOU)/U.arr(i,j,k,RHO);
                U.arr(i,j,k,P) = (U.arr(i,j,k,TOTAL_E)-0.5*U.arr(i,j,k,RHO)*U.arr(i,j,k,VELOCITY)*U.arr(i,j,k,VELOCITY))*(parameters.adiabaticIndex-1.0);
            }
        }
    }
}

void CellArray::primitiveToConservative(BoxAccessCellArray& U, ParameterStruct& parameters)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                U.arr(i,j,k,RHOU) 		= U.arr(i,j,k,RHO)*U.arr(i,j,k,VELOCITY);
                U.arr(i,j,k,TOTAL_E)   	= U.arr(i,j,k,P)/(parameters.adiabaticIndex-1.0) + 0.5*U.arr(i,j,k,RHO)*U.arr(i,j,k,VELOCITY)*U.arr(i,j,k,VELOCITY);
            }
        }
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


void CellArray::getSoundSpeed(ParameterStruct& parameters)
{
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        FArrayBox& fab = data[mfi];

        Array4<Real> const& prop_arr = fab.array();

        BoxAccessCellArray baca(bx,fab,prop_arr);

        getSoundSpeed(baca,parameters);
    }
}

void CellArray::getSoundSpeed(BoxAccessCellArray& U, ParameterStruct& parameters)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                U.arr(i,j,k,SOUNDSPEED) = sqrt(U.arr(i,j,k,P)*(parameters.adiabaticIndex)/U.arr(i,j,k,RHO));
            }
        }
    }
}
