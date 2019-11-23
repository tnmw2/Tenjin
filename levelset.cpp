#include "simulationheader.h"

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

void  LevelSet::advanceLevelSet(MultiFab& S, CellArray& U, LevelSet& LS, Real dt, const Real* dx, Vector<BCRec>& LevelSet_bc, Geometry& geom)
{
    LS.data.FillBoundary(geom.periodicity());
    FillDomainBoundary(LS.data, geom, LevelSet_bc);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,U);
        BoxAccessLevelSet  bals_new(mfi,bx,(*this));
        BoxAccessLevelSet  bals_old(mfi,bx,LS);

        bals_new.advanceLevelSet(baca,bals_old,dt,dx);
    }

    data.FillBoundary(geom.periodicity());
    FillDomainBoundary(data, geom, LevelSet_bc);
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

