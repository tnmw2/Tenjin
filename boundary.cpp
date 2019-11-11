#include "simulationheader.h"

void setBoundaryConditions(Vector<BCRec>& bc, ParameterStruct& parameters, InitialStruct& initial, AccessPattern& accessPattern)
{
    int row, col;

    for(int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        for(int n = 0; n < parameters.Ncomp; ++n)
        {
            if(initial.lowBoundary[dir] == REFLECTIVE)
            {
                if(n == accessPattern[VELOCITY]+dir || n == accessPattern[RHOU]+dir)
                {
                    bc[n].setLo(dir, BCType::reflect_odd);
                }
                else if(parameters.SOLID)
                {
                    if( n >= accessPattern[V_TENSOR] && n < accessPattern[V_TENSOR] + 9)
                    {
                        row = (n - accessPattern[V_TENSOR])/3;
                        col = (n - accessPattern[V_TENSOR])%3;

                        if( ((row == dir) && (col != dir)) || ((row != dir) && (col == dir)))
                        {
                            bc[n].setLo(dir, BCType::reflect_odd);
                        }
                        else
                        {
                            bc[n].setLo(dir, BCType::reflect_even);
                        }
                    }
                    else if( n >= accessPattern[DEVH] && n < accessPattern[DEVH] + 9)
                    {
                        row = (n - accessPattern[DEVH])/3;
                        col = (n - accessPattern[DEVH])%3;

                        if( ((row == dir) && (col != dir)) || ((row != dir) && (col == dir)))
                        {
                            bc[n].setLo(dir, BCType::reflect_odd);
                        }
                        else
                        {
                            bc[n].setLo(dir, BCType::reflect_even);
                        }
                    }
                    else if( n >= accessPattern[SIGMA] && n < accessPattern[SIGMA] + 9)
                    {
                        row = (n - accessPattern[SIGMA])/3;
                        col = (n - accessPattern[SIGMA])%3;

                        if( ((row == dir) && (col != dir)) || ((row != dir) && (col == dir)))
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
                        bc[n].setLo(dir, BCType::reflect_even);
                    }
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
                if(n == accessPattern[VELOCITY]+dir || n == accessPattern[RHOU]+dir)
                {
                    bc[n].setHi(dir, BCType::reflect_odd);
                }
                else if(parameters.SOLID)
                {
                    if( n >= accessPattern[V_TENSOR] && n < accessPattern[V_TENSOR] + 9)
                    {
                        row = (n - accessPattern[V_TENSOR])/3;
                        col = (n - accessPattern[V_TENSOR])%3;

                        if( ((row == dir) && (col != dir)) || ((row != dir) && (col == dir)))
                        {
                            bc[n].setHi(dir, BCType::reflect_odd);
                        }
                        else
                        {
                            bc[n].setHi(dir, BCType::reflect_even);
                        }
                    }
                    else if( n >= accessPattern[DEVH] && n < accessPattern[DEVH] + 9)
                    {
                        row = (n - accessPattern[DEVH])/3;
                        col = (n - accessPattern[DEVH])%3;

                        if( ((row == dir) && (col != dir)) || ((row != dir) && (col == dir)))
                        {
                            bc[n].setHi(dir, BCType::reflect_odd);
                        }
                        else
                        {
                            bc[n].setHi(dir, BCType::reflect_even);
                        }
                    }
                    else if( n >= accessPattern[SIGMA] && n < accessPattern[SIGMA] + 9)
                    {
                        row = (n - accessPattern[SIGMA])/3;
                        col = (n - accessPattern[SIGMA])%3;

                        if( ((row == dir) && (col != dir)) || ((row != dir) && (col == dir)))
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
                        bc[n].setHi(dir, BCType::reflect_even);
                    }
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
