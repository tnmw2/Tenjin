#include "timestep.h"

TimeStep::TimeStep(BoxArray &ba, DistributionMapping &dm, const int Ncomp, const int Nghost) : data(ba,dm,Ncomp,Nghost)
{

}

