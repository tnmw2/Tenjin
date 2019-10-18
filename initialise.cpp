#include "simulationheader.h"

void initial_conditions(BoxAccessCellArray& U, ParameterStruct& parameters, InitialStruct& initial)
{

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);


    std::vector<double> F1{1.0,0.0,0.0,-0.01,0.95,0.02,-0.015,0.0,0.9};
    std::vector<double> F2{1.0,0.0,0.0,0.015,0.95,0.0,-0.01,0.0,0.9};

    std::vector<double> V1;
    std::vector<double> V2;

    V1.resize(9);
    V2.resize(9);

    getStretchTensor(&V1[0],&F1[0]);
    getStretchTensor(&V2[0],&F2[0]);

    double norm1 = std::pow(det(&V1[0],U.numberOfComponents),-1.0/3.0);
    double norm2 = std::pow(det(&V2[0],U.numberOfComponents),-1.0/3.0);

    for(int row=0;row<U.numberOfComponents;row++)
    {
        for(int col=0;col<U.numberOfComponents;col++)
        {
            V1[row*U.numberOfComponents+col] *= norm1;
            V2[row*U.numberOfComponents+col] *= norm2;
        }
    }

    std::vector< std::vector<double> > Vs{V1,V2};
    std::vector< std::vector<double> > Fs{F1,F2};


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

                            U.accessPattern.materialInfo[m].EOS->setRhoFromDeformationTensor(U,i,j,k,m,&(Fs[s][0]));
                        }


                        for(int row = 0; row<U.numberOfComponents; row++)
                        {
                            for(int col = 0; col<U.numberOfComponents; col++)
                            {
                                //U(i,j,k,V_TENSOR,0,row,col) = delta<Real>(row,col);

                                U(i,j,k,V_TENSOR,0,row,col) = Vs[s][row*U.numberOfComponents+col];

                            }
                        }

                    }
                }
            }
        }
    }

    //U.normaliseV();
}

void initialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial)
{
    ParmParse pp;

    pp.get("max_grid_size",parameters.max_grid_size);

    pp.get("startT", initial.startT);
    pp.get("finalT", initial.finalT);

    pp.get("SOLID",  parameters.SOLID);

    pp.get("Nstates",initial.numberOfStates);

    pp.getarr("interfaces",initial.interfaces);

    pp.getarr("rho",    initial.rho);
    pp.getarr("u",      initial.u);
    pp.getarr("v",      initial.v);
    pp.getarr("w",      initial.w);
    pp.getarr("p",      initial.p);
    pp.getarr("alpha",  initial.alpha);
    pp.getarr("lambda", initial.lambda);
    pp.getarr("gamma",  parameters.adiabaticIndex);
    pp.getarr("CV",     parameters.CV);

    Vector<Real> pref;
    Vector<Real> eref;
    Vector<Real> mixpref;
    Vector<Real> mixeref;

    pp.getarr("pref",pref);
    pp.getarr("eref",eref);


    pp.get("CFL", parameters.CFL);

    pp.get("Nghost",    parameters.Nghost);
    pp.get("Nmat",      parameters.numberOfMaterials);
    pp.get("Nmix",      parameters.numberOfMixtures);

    pp.getarr("mixGamma",   parameters.mixtureAdiabaticIndex);
    pp.getarr("mixCV",      parameters.mixtureCV);
    pp.getarr("mixpref",    mixpref);
    pp.getarr("mixeref",    mixeref);


    pp.getarr("dimL",       parameters.dimL);
    pp.getarr("n_cells" ,   parameters.n_cells);

    pp.get("plotDirectory",initial.filename);

    int m = parameters.numberOfMaterials;
    int mix = parameters.numberOfMixtures;

    parameters.Ncomp = ((2+2)*mix)+m+m+m+3+3+1+1+1+1+1+9+(9+9+9+1)*parameters.SOLID; //(rho_mix,alpharholambda,lambda),alpha,alpharho,rho_k,u,rhou,E,p,soundspeed,ustar,sigma,V,VStar,devH,HJ2

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
            //parameters.materialInfo[m].EOS = new MieGruneisenEOS(parameters.adiabaticIndex[m],pref[m],eref[m],parameters.CV[m]);


            parameters.materialInfo[m].EOS = new RomenskiiSolidEOS();

            if(m==0)
            {
                Vector<Real> temp{parameters.adiabaticIndex[m],pref[m],eref[m],parameters.CV[m],2700.0,76.3E9,0.63,2.29,26.1E9};

                parameters.materialInfo[m].EOS->define(temp);
            }
            else
            {
                Vector<Real> temp{parameters.adiabaticIndex[m],pref[m],eref[m],parameters.CV[m],8930.0,136.5E9,1.0,3.0,39.4E9};

                parameters.materialInfo[m].EOS->define(temp);
            }
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


