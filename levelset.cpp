#include "levelset.h"

LevelSet::LevelSet(MultiFab& S, ParameterStruct& parameters) : data(S), NLevelSets(parameters.NLevelSets){}

void LevelSet::initialise(const Real* dx, const Real* prob_lo)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessLevelSet bals(mfi,bx,(*this));

        bals.initialise(dx,prob_lo);
    }


}

void  LevelSet::resetLevelSet(MultiFab& S_new)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessLevelSet  bals(mfi,bx,(*this));

        bals.resetLevelSet();
    }
}

void  LevelSet::sweep(MultiFab& S_new, const Real* dx, Geometry& geom, Vector<BCRec>& levelSet_bc, int dir, int sense, int sign)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.validbox();

        BoxAccessLevelSet bals(mfi,box,(*this));

        bals.fastSweep(dx,dir,sense,sign);
    }

    data.FillBoundary(geom.periodicity());
    FillDomainBoundary(data, geom, levelSet_bc);
}

