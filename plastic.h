#ifndef PLASTICEOS_H
#define PLASTICEOS_H

#include "cellarray.h"


class PlasticEOS
{
public:
    PlasticEOS();

    Real plasticStrainRate(double Jnew, double J, BoxAccessCellArray& U, int i, int j, int k, ParameterStruct& parameters, int m, Real Tstar);
    Real bisection(BoxAccessCellArray& U, int i, int j, int k, Real J, ParameterStruct& parameters, int m, Real dt);
    Real bisectionFunction(Real Jnew, Real J, BoxAccessCellArray& U, int i, int j, int k, ParameterStruct& parameters, int m, Real dt, Real Tstar);
    Real epsilonFunction(double J, double Jnew, BoxAccessCellArray& U, int i, int j, int k, int m);

    bool overYieldStress    (Real& J, BoxAccessCellArray& U, int i, int j, int k, int m);
    void boxPlasticUpdate   (BoxAccessCellArray& U, ParameterStruct& parameters, Real dt);
    void SVDboxPlasticUpdate(BoxAccessCellArray& U, ParameterStruct& parameters, Real dt);
    void plasticUpdate      (CellArray& U,          ParameterStruct& parameters, Real dt, MultiFab &S_new);


    Vector<Real> yieldStress;

    /*Real c1 = 0.324E9; //AL
    Real c2 = 0.114E9;
    Real c3 = 0.002;
    Real n  = 0.42;
    */

    Real c1 = 0.4E9;  //Cu
    Real c2 = 0.177E9;
    Real c3 = 0.025;
    Real n  = 1.0;

    Real mt = 1.0;

    Real udaykumarConstant = 100E6;

};


#endif // PLASTICEOS_H
