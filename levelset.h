#ifndef LEVELSET_H
#define LEVELSET_H

#include "amrexheader.h"
#include "structdefinitions.h"
#include "simulationheader.h"


class LevelSet
{
public:

    LevelSet(MultiFab& S, ParameterStruct& parameters);

    MultiFab& data;

    int NLevelSets;

    void initialise         (const Real* dx, const Real* prob_lo);
    void advanceLevelSet    (MultiFab& S, CellArray& U, LevelSet& LS, Real dt, const Real* dx, Vector<BCRec> &LevelSet_bc, Geometry &geom);
    void resetLevelSet      (MultiFab &S_new);
    void sweep              (MultiFab &S_new, const Real* dx, Geometry& geom, Vector<BCRec>& levelSet_bc, int dir, int sense, int sign);

};

class BoxAccessLevelSet
{
public:

    BoxAccessLevelSet(MFIter& mfi, const Box &bx, LevelSet &U);

    const Box& box;

    FArrayBox& fab;

    int NLevelSets;

    void initialise(const Real* dx, const Real* prob_lo);

    Real& operator()(int i, int j, int k, int mat);

    void  advanceLevelSet        (BoxAccessCellArray& U, BoxAccessLevelSet& LS, Real dt, const Real* dx);
    Real  levelSetDerivative     (BoxAccessLevelSet& LS, Real v, const Real* dx, int dir, int i, int j, int k, int n);
    Real  D1                     (BoxAccessLevelSet& LS, int dir, int i, int j, int k, int n, int sign, const Real *dx);
    void  resetLevelSet          ();
    void  fastSweep              (const Real *dx, int dir, int sense, int sign);
    bool  cellIsNextToAnInterface(int i, int j, int k, int n);


};


#endif // LEVELSET_H
