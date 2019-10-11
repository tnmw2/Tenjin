#include "simulationheader.h"

void initial_conditions(BoxAccessCellArray& U, ParameterStruct& parameters, InitialStruct& initial)
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
                if(parameters.numberOfMaterials == 1)
                {
                    if(i <= int_x0)
                    {
                        U(i,j,k,ALPHA)     = 1.0;
                        U(i,j,k,RHO_K)     = initial.rhoL;
                        U(i,j,k,VELOCITY)  = initial.uL;
                        U(i,j,k,P)         = initial.pL;
                    }
                    else
                    {
                        U(i,j,k,ALPHA)     = 1.0;
                        U(i,j,k,RHO_K)     = initial.rhoR;
                        U(i,j,k,VELOCITY)  = initial.uR;
                        U(i,j,k,P)         = initial.pR;
                    }
                }
                else if(parameters.numberOfMaterials == 2)
                {
                    if(i <= int_x0)
                    {
                        U(i,j,k,ALPHA,0)     = 0.999999;
                        U(i,j,k,ALPHA,1)     = 0.000001;
                        U(i,j,k,RHO_K,0)     = initial.rhoL;
                        U(i,j,k,RHO_K,1)     = initial.rhoR;
                        U(i,j,k,VELOCITY)    = initial.uL;
                        U(i,j,k,P)           = initial.pL;
                    }
                    else
                    {
                        U(i,j,k,ALPHA,1)     = 0.999999;
                        U(i,j,k,ALPHA,0)     = 0.000001;
                        U(i,j,k,RHO_K,0)     = initial.rhoL;
                        U(i,j,k,RHO_K,1)     = initial.rhoR;
                        U(i,j,k,VELOCITY)    = initial.uR;
                        U(i,j,k,P)           = initial.pR;
                    }
                }
            }
        }
    }

}

void initialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial)
{
    ParmParse pp;

    pp.get("max_grid_size",parameters.max_grid_size);

    pp.get("startT", initial.startT);
    pp.get("finalT", initial.finalT);

    pp.get("rhoL",  initial.rhoL);
    pp.get("rhoR",  initial.rhoR);
    pp.get("uL",    initial.uL);
    pp.get("uR",    initial.uR);
    pp.get("pL",    initial.pL);
    pp.get("pR",    initial.pR);

    pp.getarr("gamma",parameters.adiabaticIndex);

    pp.get("x0",  parameters.x0);
    pp.get("CFL", parameters.CFL);


    pp.get("Nghost", parameters.Nghost);
    pp.get("Nmat", parameters.numberOfMaterials);

    pp.getarr("dimL", parameters.dimL);
    pp.getarr("n_cells" , parameters.n_cells);


    int m = parameters.numberOfMaterials;

    parameters.Ncomp = m+m+m+3+3+1+1+1+1+1; //alpha,alpharho,rho_k,u,rhou,E,p,soundspeed,ustar


    pp.get("plotDirectory",initial.filename);
}

void setInitialConditions(CellArray& U, ParameterStruct& parameters, InitialStruct& initial)
{
    for(MFIter mfi(U.data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,U);

        initial_conditions(baca, parameters, initial);

        baca.primitiveToConservative(parameters);
    }
}


