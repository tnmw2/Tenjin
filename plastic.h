#ifndef PLASTICEOS_H
#define PLASTICEOS_H

#include "cellarray.h"


class PlasticEOS
{
public:
    PlasticEOS();

    Real plasticStrainRate(double Jnew, double J, BoxAccessCellArray& U, int i, int j, int k, ParameterStruct& parameters,int m);
    Real bisection(BoxAccessCellArray& U, int i, int j, int k, Real J, ParameterStruct& parameters, int m, Real dt);
    Real bisectionFunction(Real Jnew, Real J, BoxAccessCellArray& U, int i, int j, int k, ParameterStruct& parameters, int m, Real dt);
    Real epsilonFunction(double J, double Jnew, BoxAccessCellArray& U, int i, int j, int k, int m);

    bool    overYieldStress(Real& J, BoxAccessCellArray& U, int i, int j, int k, int m);
    void  	boxPlasticUpdate(BoxAccessCellArray& U, ParameterStruct& parameters, Real dt);
    void  	plasticUpdate   (CellArray& U,          ParameterStruct& parameters, Real dt, MultiFab &S_new);

    Vector<Real> yieldStress;

    Real c1 = 0.324E9; //60E6;  //
    Real c2 = 0.114E9; //500E6; //
    Real c3 = 0.002;   //0.02;  //
    Real n  = 0.42;    //0.3;   //

};


#endif // PLASTICEOS_H
