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
                    if(i<int_x0)
                    {
                        U(i,j,k,ALPHA)     = 1.0;
                        U(i,j,k,RHO_K)     = initial.rhoL;
                        U(i,j,k,P)         = initial.pL;

                        for(int dir = 0; dir< AMREX_SPACEDIM ; dir++)
                        {
                            U(i,j,k,VELOCITY,0,dir)  = initial.uL[dir];
                        }

                        if(parameters.numberOfMixtures == 1)
                        {
                            U(i,j,k,RHO_MIX,0,0)  = initial.rhoL;
                            U(i,j,k,RHO_MIX,0,1)  = initial.rhoR;
                            U(i,j,k,LAMBDA ,0)    = 0.999999;
                        }


                    }
                    else
                    {
                        U(i,j,k,ALPHA)     = 1.0;
                        U(i,j,k,RHO_K)     = initial.rhoR;
                        U(i,j,k,P)         = initial.pR;

                        for(int dir=0;dir<AMREX_SPACEDIM;dir++)
                        {
                            U(i,j,k,VELOCITY,0,dir)  = initial.uR[dir];
                        }

                        if(parameters.numberOfMixtures == 1)
                        {
                            U(i,j,k,RHO_MIX,0,0)  = initial.rhoL;
                            U(i,j,k,RHO_MIX,0,1)  = initial.rhoR;
                            U(i,j,k,LAMBDA ,0)    = 0.000001;
                        }

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
                        U(i,j,k,P)           = initial.pL;

                        for(int dir=0; dir< AMREX_SPACEDIM ; dir++)
                        {
                            U(i,j,k,VELOCITY,0,dir)  = initial.uL[dir];
                        }
                    }
                    else
                    {
                        U(i,j,k,ALPHA,1)     = 0.999999;
                        U(i,j,k,ALPHA,0)     = 0.000001;
                        U(i,j,k,RHO_K,0)     = initial.rhoL;
                        U(i,j,k,RHO_K,1)     = initial.rhoR;
                        U(i,j,k,P)           = initial.pR;

                        for(int dir=0;dir<AMREX_SPACEDIM;dir++)
                        {
                            U(i,j,k,VELOCITY,0,dir)  = initial.uR[dir];
                        }
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
    pp.getarr("uL", initial.uL);
    pp.getarr("uR", initial.uR);
    pp.get("pL",    initial.pL);
    pp.get("pR",    initial.pR);

    pp.getarr("gamma",parameters.adiabaticIndex);
    pp.getarr("CV",parameters.CV);

    Vector<Real> pref;
    Vector<Real> eref;
    Vector<Real> mixpref;
    Vector<Real> mixeref;

    pp.getarr("pref",pref);
    pp.getarr("eref",eref);


    pp.get("x0",  parameters.x0);
    pp.get("CFL", parameters.CFL);

    pp.get("Nghost", parameters.Nghost);
    pp.get("Nmat", parameters.numberOfMaterials);
    pp.get("Nmix", parameters.numberOfMixtures);

    pp.getarr("mixGamma",   parameters.mixtureAdiabaticIndex);
    pp.getarr("mixCV",      parameters.mixtureCV);
    pp.getarr("mixpref",    mixpref);
    pp.getarr("mixeref",    mixeref);


    pp.getarr("dimL", parameters.dimL);
    pp.getarr("n_cells" , parameters.n_cells);

    pp.get("plotDirectory",initial.filename);

    int m = parameters.numberOfMaterials;
    int mix = parameters.numberOfMixtures;

    parameters.Ncomp = ((2+2)*mix)+m+m+m+AMREX_SPACEDIM+AMREX_SPACEDIM+1+1+1+1+1; //(rho_mix,alpharholambda,lambda),alpha,alpharho,rho_k,u,rhou,E,p,soundspeed,ustar

    parameters.materialInfo.resize(parameters.numberOfMaterials);

    Vector<int> temp(parameters.numberOfMaterials);

    pp.getarr("mixture",temp);

    int counter = 0;

    for(int m = 0; m<parameters.numberOfMaterials; m++)
    {
        parameters.materialInfo[m].mixture = (bool)temp[m];

        if(parameters.materialInfo[m].mixture)
        {
            parameters.materialInfo[m].EOS = new MixtureEOS();

            Vector<Real> temp{parameters.adiabaticIndex[m],pref[m],eref[m],parameters.CV[m],parameters.mixtureAdiabaticIndex[counter],mixpref[counter],mixeref[counter],parameters.mixtureCV[counter]};

            parameters.materialInfo[m].EOS->define(temp);

            parameters.materialInfo[m].mixtureIndex = counter;
            counter++;

         }
        else
        {
            parameters.materialInfo[m].EOS = new MieGruneisenEOS(parameters.adiabaticIndex[m],pref[m],eref[m],parameters.CV[m]);
        }
    }


}

void setInitialConditions(CellArray& U, ParameterStruct& parameters, InitialStruct& initial)
{
    for(MFIter mfi(U.data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,U);

        initial_conditions(baca, parameters, initial);

        baca.primitiveToConservative();
    }
}


