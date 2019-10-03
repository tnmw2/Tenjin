#include "simulationheader.h"

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
    pp.get("a", parameters.a);
    pp.getarr("dimL", parameters.dimL);
    pp.getarr("n_cells" , parameters.n_cells);
    pp.get("Ncomp", parameters.Ncomp);
    pp.get("Nghost", parameters.Nghost);
}

void setInitialConditions(MultiFab& phi, ParameterStruct& parameters, InitialStruct& initial)
{
    for(MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        FArrayBox& fab = phi[mfi];

        Array4<Real> const& prop_arr = fab.array();

        initial_conditions(bx, prop_arr, parameters,initial);

        primtiveToConservative(bx, prop_arr ,parameters);

    }
}
