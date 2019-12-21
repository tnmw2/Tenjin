#include "simulationheader.h"

void setBoundaryConditions(Vector<BCRec>& bc, ParameterStruct& parameters, InitialStruct& initial, AccessPattern& accessPattern)
{
    int row, col;

    for(int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        int counter = 0;

        for(auto n : accessPattern.allVariables)
        {
            if(initial.lowBoundary[dir] == REFLECTIVE)
            {
                if( (n.var == VELOCITY && n.row == dir) || (n.var == RHOU && n.row == dir))
                {
                    bc[counter].setLo(dir, BCType::reflect_odd);
                }
                else
                {
                    bc[counter].setLo(dir, BCType::reflect_even);

                }
            }
            else
            {
                bc[counter].setLo(dir, BCType::foextrap);
            }

            if(initial.highBoundary[dir] == REFLECTIVE)
            {
                if( (n.var == VELOCITY && n.row == dir) || (n.var == RHOU && n.row == dir))
                {
                    bc[counter].setHi(dir, BCType::reflect_odd);
                }
                else
                {
                    bc[counter].setHi(dir, BCType::reflect_even);

                }
            }
            else
            {
                bc[counter].setHi(dir, BCType::foextrap);
            }

            counter++;
        }
    }
}
