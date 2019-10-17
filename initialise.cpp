#include "simulationheader.h"

void initial_conditions(BoxAccessCellArray& U, ParameterStruct& parameters, InitialStruct& initial)
{

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    Vector<int> boundary;

    boundary.push_back(0);

    for(int s = 0; s < initial.numberOfStates-1; s++)
    {
        boundary.push_back((int)((initial.interfaces[s]/parameters.dimL[0])*parameters.n_cells[0]));
    }

    boundary.push_back(parameters.n_cells[0]);

    for                 (int k = lo.z; k <= hi.z; ++k)
    {
        for             (int j = lo.y; j <= hi.y; ++j)
        {
            for         (int i = lo.x; i <= hi.x; ++i)
            {
                for     (int s = 0; s < initial.numberOfStates; s++)
                {
                    if(boundary[s] <= i && i < boundary[s+1])
                    {
                        for (int m = 0; m < parameters.numberOfMaterials; m++)
                        {
                            if(initial.alpha[s] == m)
                            {
                                U(i,j,k,ALPHA,m)     = 1.0 - std::max(1E-6,((Real)(parameters.numberOfMaterials-1))*1E-6);
                            }
                            else
                            {
                                U(i,j,k,ALPHA,m)     = 1E-6;
                            }

                            U(i,j,k,RHO_K,m)    = initial.rho[s];
                            U(i,j,k,P)          = initial.p[s];

                            U(i,j,k,VELOCITY,0,0)  = initial.u[s];
                            U(i,j,k,VELOCITY,0,1)  = initial.v[s];
                            U(i,j,k,VELOCITY,0,2)  = initial.w[s];

                            if(U.accessPattern.materialInfo[m].mixture)
                            {
                                U(i,j,k,RHO_MIX,m,0)    = initial.rho[s];
                                U(i,j,k,RHO_MIX,m,1)    = initial.rho[s];
                                if(initial.lambda[s])
                                {
                                    U(i,j,k,LAMBDA,m)   = 0.999999;
                                }
                                else
                                {
                                    U(i,j,k,LAMBDA,m)   = 0.0;
                                }
                            }
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

    pp.get("Nstates",initial.numberOfStates);

    pp.getarr("interfaces",initial.interfaces);

    pp.getarr("rho",    initial.rho);
    pp.getarr("u",      initial.u);
    pp.getarr("v",      initial.v);
    pp.getarr("w",      initial.w);
    pp.getarr("p",      initial.p);
    pp.getarr("alpha",  initial.alpha);
    pp.getarr("lambda", initial.lambda);
    pp.getarr("gamma",parameters.adiabaticIndex);
    pp.getarr("CV",parameters.CV);

    Vector<Real> pref;
    Vector<Real> eref;
    Vector<Real> mixpref;
    Vector<Real> mixeref;

    pp.getarr("pref",pref);
    pp.getarr("eref",eref);


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

    parameters.Ncomp = ((2+2)*mix)+m+m+m+3+3+1+1+1+1+1+9; //(rho_mix,alpharholambda,lambda),alpha,alpharho,rho_k,u,rhou,E,p,soundspeed,ustar, sigma

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


