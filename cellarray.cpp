#include "cellarray.h"
#include "simulationheader.h"

CellArray::CellArray(MultiFab& S, AccessPattern &_accessPattern, ParameterStruct& parameters) : data(S), accessPattern(_accessPattern), numberOfMaterials{parameters.numberOfMaterials}{}

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

void CellArray::cleanUpAlpha()
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,(*this));

        baca.cleanUpAlpha();
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

int CellArray::getArrayPosition(MaterialSpecifier& m)
{
    switch(m.var)
    {
        case RHO:               return (accessPattern[m.var]+m.mat);
            break;
        case RHOU:              return (accessPattern[m.var]+m.mat*numberOfComponents+m.row);
            break;
        case TOTAL_E:           return (accessPattern[m.var]+m.mat);
            break;
        case VELOCITY:          return (accessPattern[m.var]+m.mat*numberOfComponents+m.row);
            break;
        case P:                 return (accessPattern[m.var]+m.mat);
            break;
        case SOUNDSPEED:        return (accessPattern[m.var]+m.mat);
            break;
        case USTAR:             return (accessPattern[m.var]+m.mat);
            break;
        case SIGMA:             return (accessPattern[m.var]+m.mat*numberOfComponents*numberOfComponents+m.row*numberOfComponents+m.col);
            break;
        case V_TENSOR:          return (accessPattern[m.var]+m.mat*numberOfComponents*numberOfComponents+m.row*numberOfComponents+m.col);
            break;
        case DEVH:              return (accessPattern[m.var]+m.mat*numberOfComponents*numberOfComponents+m.row*numberOfComponents+m.col);
            break;
        case VSTAR:             return (accessPattern[m.var]+m.mat*numberOfComponents*numberOfComponents+m.row*numberOfComponents+m.col);
            break;
        case HJ2:               return (accessPattern[m.var]+m.mat);
            break;
        case EPSILON:           return (accessPattern[m.var]+m.mat);
            break;
        case ALPHARHOEPSILON:   return (accessPattern[m.var]+m.mat);
            break;
    default: Print() << "Incorrect Access variable " << m.var << " in boxaccesscellarray: " << accessPattern.variableNames[getArrayPosition(m)] << std::endl;
        exit(1);
    }
}


