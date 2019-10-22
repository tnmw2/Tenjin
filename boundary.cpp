#include "simulationheader.h"

void setBoundaryConditions(Vector<BCRec>& bc, ParameterStruct& parameters, InitialStruct& initial, AccessPattern& accessPattern)
{
    for(int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        for(int n = 0; n < parameters.Ncomp; ++n)
        {
            if(initial.lowBoundary[dir] == REFLECTIVE)
            {
                if(n == accessPattern[VELOCITY]+1 || n == accessPattern[RHOU]+1)
                {
                    bc[n].setLo(dir, BCType::reflect_odd);
                }
                else
                {
                    bc[n].setLo(dir, BCType::reflect_even);
                }
            }
            else
            {
                bc[n].setLo(dir, BCType::foextrap);
            }

            if(initial.highBoundary[dir] == REFLECTIVE)
            {
                if(n == accessPattern[VELOCITY]+1 || n == accessPattern[RHOU]+1)
                {
                    bc[n].setHi(dir, BCType::reflect_odd);
                }
                else
                {
                    bc[n].setHi(dir, BCType::reflect_even);
                }
            }
            else
            {
                bc[n].setHi(dir, BCType::foextrap);
            }
        }
    }
}
