#ifndef LEVELSET_H
#define LEVELSET_H

#include "amrexheader.h"
#include "structdefinitions.h"

class LevelSet
{
public:

    LevelSet(MultiFab& S, ParameterStruct& parameters);

    MultiFab& data;

    int NLevelSets;

    void initialise         (InitialStruct& initial, const Real* dx, const Real* prob_lo);
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

    void initialise(InitialStruct& initial, const Real* dx, const Real* prob_lo);

    Real& operator()(int i, int j, int k, int mat=0);

    Real  levelSetDerivative        (Real v, const Real* dx, int dir, int i, int j, int k, int n);
    Real  D1                        (int dir, int i, int j, int k, int n, int sign, const Real *dx);
    Real  D2                        (int dir, int i, int j, int k, int n, int sign, const Real* dx);
    Real  D3                        (int dir, int i, int j, int k, int n, int sign, const Real* dx);
    void  resetLevelSet             ();
    void  fastSweep                 (const Real *dx, int xsense, int ysense, int sign);
    bool  cellIsNextToAnInterface   (int i , int j , int k, int n, int limiter=-1);
    void  calculateNormal           (int i , int j , int k, int n, const Real* dx, Real& nx, Real& ny);
    void  calculateInterpolationPoint(int i , int j , int k, int n, const Real* dx, Real& nx, Real& ny, Real& cx, Real& cy, Real& ix, Real& iy);
    void  calculateProbes            (int i , int j , int k, int n, const Real* dx, Real& nx, Real& ny, Real& ix, Real& iy, Vector<Real>& px, Vector<Real>& py);

    bool customComparator(int i, int lim, int sense);
    void customChanger(int& i, int sense);

    bool cellIsValid(int i, int j, int k, int m);
    int whatMaterialIsValid(int i, int j, int k);
    bool cellIsNearInterface(int i, int j, int k, const Real* dx);




};


#endif // LEVELSET_H
