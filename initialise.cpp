#include "simulationheader.h"

void initial_conditions(BoxAccessCellArray& U, ParameterStruct& parameters)
{

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);


    int int_x0 = (parameters.x0/parameters.dimL[0])*parameters.n_cells[0];

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                if(i <= int_x0)
                {
                    U(i,j,k,RHO)         = 1.0;
                    U(i,j,k,VELOCITY)    = 0.0;
                    U(i,j,k,P)           = 1.0;
                }
                else
                {
                    U(i,j,k,RHO)         = 0.125;
                    U(i,j,k,VELOCITY)    = 0.0;
                    U(i,j,k,P)           = 0.1;
                }
                /*if(i <= int_x0)
                {
                    U(i,j,k,RHO)         = 1.0;
                    U(i,j,k,VELOCITY)    = -2.0;
                    U(i,j,k,P)           = 0.4;
                }
                else
                {
                    U(i,j,k,RHO)         = 1.0;
                    U(i,j,k,VELOCITY)    = 2.0;
                    U(i,j,k,P)           = 0.4;
                }*/
            }
        }
    }

}

void initialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial)
{
    ParmParse pp;

    pp.get("max_grid_size",parameters.max_grid_size);

    parameters.plot_int = 1;
    pp.query("plot_int",parameters.plot_int);

    //pp.queryarr("is_periodic",is_periodic);

    pp.get("startT", initial.startT);
    pp.get("finalT", initial.finalT);

    //Test Case specific variables:

    pp.get("phiL", parameters.phiL);
    pp.get("phiR", parameters.phiR);

    //initial discontinuity
    pp.get("x0",  parameters.x0);
    pp.get("CFL", parameters.CFL);

    pp.get("gamma", parameters.adiabaticIndex);
    pp.getarr("dimL", parameters.dimL);
    pp.getarr("n_cells" , parameters.n_cells);
    pp.get("Ncomp", parameters.Ncomp);
    pp.get("Nghost", parameters.Nghost);

    pp.get("plotDirectory",initial.filename);
}

void setInitialConditions(CellArray& U, ParameterStruct& parameters)
{
    for(MFIter mfi(U.data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        FArrayBox& fab = U.data[mfi];

        Array4<Real> const& prop_arr = fab.array();

        BoxAccessCellArray baca(bx,fab,prop_arr,U);

        initial_conditions(baca, parameters);

        U.primitiveToConservative(baca, parameters);
    }
}
