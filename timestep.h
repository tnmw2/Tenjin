#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "amrexheader.h"

class TimeStep
{
public:

    TimeStep(BoxArray& ba, DistributionMapping& dm, const int Ncomp, const int Nghost);

    MultiFab data;
};

#endif // TIMESTEP_H
