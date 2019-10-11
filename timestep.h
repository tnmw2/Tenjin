#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "amrexheader.h"
#include "structdefinitions.h"
#include "cellarray.h"

/** \class Calculates the timestep.
 */
class TimeStep
{
public:

    TimeStep(BoxArray& ba, DistributionMapping& dm, const int Ncomp, const int Nghost);
    void boxTimeStepFunction(BoxAccessCellArray& U, Array4<Real> const& prop_arr_time, ParameterStruct& parameters);
    Real getTimeStep(CellArray& U, ParameterStruct& parameters);

    MultiFab data;
};

#endif // TIMESTEP_H
