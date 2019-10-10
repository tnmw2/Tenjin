#include "cellarray.h"


CellArray::CellArray(BoxArray& ba, DistributionMapping& dm, const int Ncomp, const int Nghost, std::map<Variable,int>& _accessPattern, ParameterStruct& parameters) : data(ba,dm,Ncomp,Nghost), accessPattern(_accessPattern), numberOfMaterials{parameters.numberOfMaterials}{}

void CellArray::conservativeToPrimitive(ParameterStruct& parameters)
{
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        FArrayBox& fab = data[mfi];

        Array4<Real> const& prop_arr = fab.array();

        BoxAccessCellArray baca(bx,fab,prop_arr,(*this));

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

        BoxAccessCellArray baca(bx,fab,prop_arr,(*this));

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
                U(i,j,k,VELOCITY) = U(i,j,k,RHOU)/U(i,j,k,RHO);
                U(i,j,k,P) = (U(i,j,k,TOTAL_E)-0.5*U(i,j,k,RHO)*U(i,j,k,VELOCITY)*U(i,j,k,VELOCITY))*(parameters.adiabaticIndex-1.0);
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
                U(i,j,k,RHOU) 		= U(i,j,k,RHO)*U(i,j,k,VELOCITY);
                U(i,j,k,TOTAL_E)   	= U(i,j,k,P)/(parameters.adiabaticIndex-1.0) + 0.5*U(i,j,k,RHO)*U(i,j,k,VELOCITY)*U(i,j,k,VELOCITY);
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

        BoxAccessCellArray baca(bx,fab,prop_arr,(*this));

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
                U(i,j,k,SOUNDSPEED) = sqrt(U(i,j,k,P)*(parameters.adiabaticIndex)/U(i,j,k,RHO));
            }
        }
    }
}
