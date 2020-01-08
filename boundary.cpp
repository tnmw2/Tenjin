#include "simulationheader.h"

void setBoundaryConditions(CellArray& U, Vector<BCRec>& bc, ParameterStruct& parameters, InitialStruct& initial, AccessPattern& accessPattern)
{

    for(int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        for(auto n : accessPattern.allVariables)
        {

            if(initial.lowBoundary[dir] == REFLECTIVE)
            {
                if((n.var == VELOCITY || n.var == RHOU) && n.row == dir)
                {
                    bc[U.getArrayPosition(n)].setLo(dir, BCType::reflect_odd);
                }
                else if((n.var == V_TENSOR || n.var == DEVH || n.var == SIGMA) && (((n.row == dir) && (n.col != dir)) || ((n.row != dir) && (n.col == dir))))
                {
                    bc[U.getArrayPosition(n)].setLo(dir, BCType::reflect_odd);
                }
                else
                {
                    bc[U.getArrayPosition(n)].setLo(dir, BCType::reflect_even);
                }
            }
            else
            {
                bc[U.getArrayPosition(n)].setLo(dir, BCType::foextrap);
            }

            if(initial.highBoundary[dir] == REFLECTIVE)
            {
                if((n.var == VELOCITY || n.var == RHOU) && n.row == dir)
                {
                    bc[U.getArrayPosition(n)].setHi(dir, BCType::reflect_odd);
                }
                else if((n.var == V_TENSOR || n.var == DEVH || n.var == SIGMA) && (((n.row == dir) && (n.col != dir)) || ((n.row != dir) && (n.col == dir))))
                {
                    bc[U.getArrayPosition(n)].setHi(dir, BCType::reflect_odd);
                }
                else
                {
                    bc[U.getArrayPosition(n)].setHi(dir, BCType::reflect_even);
                }
            }
            else
            {
                bc[U.getArrayPosition(n)].setHi(dir, BCType::foextrap);
            }
        }
    }
}
