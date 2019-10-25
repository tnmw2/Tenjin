#ifndef PLASTICEOS_H
#define PLASTICEOS_H

#include "cellarray.h"


class PlasticEOS
{
public:
    PlasticEOS();

    Real plasticStrainRate(double Jnew, double J, BoxAccessCellArray& U, int i, int j, int k, ParameterStruct& parameters,int m);
    Real bisection(BoxAccessCellArray& U, int i, int j, int k, Real J, ParameterStruct& parameters, int m);
    Real bisectionFunction(Real Jnew, Real J, BoxAccessCellArray& U, int i, int j, int k, ParameterStruct& parameters, int m);

    bool    overYieldStress(Real& J, BoxAccessCellArray& U, int i, int j, int k, int m);
    void  	boxPlasticUpdate(BoxAccessCellArray& U, ParameterStruct& parameters);
    void  	plasticUpdate   (CellArray& U,          ParameterStruct& parameters);

    Vector<Real> yieldStress;

};


#endif // PLASTICEOS_H
