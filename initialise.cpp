#include "simulationheader.h"

/** Convert a c++ string to a char array for file I/O operations
 */
void stringtochar(std::string string,char* file)
{
    int i=0;
    for(; string[i] ; i++)
    {
        file[i]=string[i];
    }
    file[i]='\0';

    return;
}

void BANKS_initial_conditions(BoxAccessCellArray& U, ParameterStruct& parameters, InitialStruct& initial)
{

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    int s = 0;

    int radiusInt       = (int)((initial.interfaces[0]/parameters.dimL[1])*parameters.n_cells[1]);
    int middleInt       = (int)(parameters.n_cells[1]);
    int startOfTubeInt  = (int)((initial.interfaces[0]/parameters.dimL[0])*parameters.n_cells[0]);
    int endOfBoosterInt = (int)((2.0*initial.interfaces[0]/parameters.dimL[0])*parameters.n_cells[0]);

    for                 (int k = lo.z; k <= hi.z; ++k)
    {
        for             (int j = lo.y; j <= hi.y; ++j)
        {
            for         (int i = lo.x; i <= hi.x; ++i)
            {
                if(i<startOfTubeInt ||  j< middleInt-radiusInt) //j>middleInt+radiusInt ||
                {
                    s=2;
                }
                else if(i<endOfBoosterInt)
                {
                    s=0;
                }
                else
                {
                    s=1;
                }

                for (int m = 0; m < parameters.numberOfMaterials; m++)
                {
                    U(i,j,k,ALPHA,m)    = initial.alpha[s][m];
                    U(i,j,k,RHO_K,m)    = initial.rho[s][m];

                    if(U.accessPattern.materialInfo[m].mixture)
                    {
                        U(i,j,k,RHO_MIX,m,0)    = initial.rho[s][m];
                        U(i,j,k,RHO_MIX,m,1)    = initial.rho[s][m];

                        U(i,j,k,LAMBDA,m)       = initial.lambda[s][m];
                    }
                }

                U(i,j,k,P)             = initial.p[s];

                U(i,j,k,VELOCITY,0,0)  = initial.u[s];
                U(i,j,k,VELOCITY,0,1)  = initial.v[s];
                U(i,j,k,VELOCITY,0,2)  = initial.w[s];

            }
        }
    }
}

void BANKS_SOLID_initial_conditions(BoxAccessCellArray& U, ParameterStruct& parameters, InitialStruct& initial)
{

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    int s = 0;

    int radiusInt       = (int)((initial.interfaces[0]/parameters.dimL[1])*parameters.n_cells[1]);
    int middleInt       = (int)(parameters.n_cells[1]);
    int startOfTubeInt  = (int)((initial.interfaces[0]/parameters.dimL[0])*parameters.n_cells[0]);
    int endOfBoosterInt = (int)((2.0*initial.interfaces[0]/parameters.dimL[0])*parameters.n_cells[0]);
    int chamfer         = (int)((0.2*initial.interfaces[0]/parameters.dimL[0])*parameters.n_cells[0]);
    std::vector< std::vector<double> > V(initial.numberOfStates);

    if(parameters.SOLID)
    {
        for(int s = 0; s < initial.numberOfStates; s++)
        {
            V[s].resize(9);

            getStretchTensor(&V[s][0],&initial.F[s][0]);

            double norm1 = std::pow(det(&V[s][0]),-1.0/3.0);

            for(int row=0;row<U.numberOfComponents;row++)
            {
                for(int col=0;col<U.numberOfComponents;col++)
                {
                    V[s][row*U.numberOfComponents+col] *= norm1;

                }
            }
        }
    }

    for                 (int k = lo.z; k <= hi.z; ++k)
    {
        for             (int j = lo.y; j <= hi.y; ++j)
        {
            for         (int i = lo.x; i <= hi.x; ++i)
            {
                if(i<startOfTubeInt ||  j< middleInt-radiusInt) //j>middleInt+radiusInt ||
                {
                    s=2;
                }
                else if(i>startOfTubeInt+chamfer && i < endOfBoosterInt)
                {
                    s=0;
                }
                else if(i<=startOfTubeInt+chamfer)
                {
                    if(j>(middleInt-radiusInt+chamfer))
                    {
                        s=0;
                    }
                    else if((j-(middleInt-radiusInt+chamfer))*(j-(middleInt-radiusInt+chamfer))+(i-(startOfTubeInt+chamfer))*(i-(startOfTubeInt+chamfer)) < chamfer*chamfer)
                    {
                        s=0;
                    }
                    else
                    {
                        s=2;
                    }
                }
                else
                {
                    s=1;
                }

                for (int m = 0; m < parameters.numberOfMaterials; m++)
                {
                    U(i,j,k,ALPHA,m)    = initial.alpha[s][m];
                    U(i,j,k,RHO_K,m)    = initial.rho[s][m];

                    if(U.accessPattern.materialInfo[m].mixture)
                    {
                        U(i,j,k,RHO_MIX,m,0)    = initial.rho[s][m];
                        U(i,j,k,RHO_MIX,m,1)    = initial.rho[s][m];

                        U(i,j,k,LAMBDA,m)       = initial.lambda[s][m];
                    }

                    if(parameters.materialInfo[m].phase == solid)
                    {
                        U.accessPattern.materialInfo[m].EOS->setRhoFromDeformationTensor(U,i,j,k,m,&initial.F[s][0]);
                    }
                }

                U(i,j,k,P)             = initial.p[s];

                U(i,j,k,VELOCITY,0,0)  = initial.u[s];
                U(i,j,k,VELOCITY,0,1)  = initial.v[s];
                U(i,j,k,VELOCITY,0,2)  = initial.w[s];

                if(parameters.SOLID)
                {
                    for(int row = 0; row<U.numberOfComponents; row++)
                    {
                        for(int col = 0; col<U.numberOfComponents; col++)
                        {
                            U(i,j,k,V_TENSOR,0,row,col) = V[s][row*U.numberOfComponents+col];
                        }
                    }
                }
            }
        }
    }
}

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

    for(int s = 0; s < initial.numberOfStates; s++)
    {
        std::vector<double> V;

        if(parameters.SOLID)
        {
            V.resize(9);

            getStretchTensor(&V[0],&initial.F[s][0]);

            double norm1 = std::pow(det(&V[0]),-1.0/3.0);

            for(int row=0;row<U.numberOfComponents;row++)
            {
                for(int col=0;col<U.numberOfComponents;col++)
                {
                    V[row*U.numberOfComponents+col] *= norm1;
                }
            }
        }

        for                 (int k = lo.z; k <= hi.z; ++k)
        {
            for             (int j = lo.y; j <= hi.y; ++j)
            {
                for         (int i = lo.x; i <= hi.x; ++i)
                {
                    if(boundary[s] <= i && i < boundary[s+1])
                    {
                        for (int m = 0; m < parameters.numberOfMaterials; m++)
                        {
                            U(i,j,k,ALPHA,m)    = initial.alpha[s][m];
                            U(i,j,k,RHO_K,m)    = initial.rho[s][m];

                            if(U.accessPattern.materialInfo[m].mixture)
                            {
                                U(i,j,k,RHO_MIX,m,0)    = initial.rho[s][m];
                                U(i,j,k,RHO_MIX,m,1)    = initial.rho[s][m];

                                U(i,j,k,LAMBDA,m)       = initial.lambda[s][m];
                            }

                            if(parameters.materialInfo[m].phase == solid)
                            {
                                U.accessPattern.materialInfo[m].EOS->setRhoFromDeformationTensor(U,i,j,k,m,&initial.F[s][0]);
                            }
                        }

                        U(i,j,k,P)             = initial.p[s];

                        U(i,j,k,VELOCITY,0,0)  = initial.u[s];
                        U(i,j,k,VELOCITY,0,1)  = initial.v[s];
                        U(i,j,k,VELOCITY,0,2)  = initial.w[s];

                        if(parameters.SOLID)
                        {
                            for(int row = 0; row<U.numberOfComponents; row++)
                            {
                                for(int col = 0; col<U.numberOfComponents; col++)
                                {
                                    U(i,j,k,V_TENSOR,0,row,col) = V[row*U.numberOfComponents+col];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //U.normaliseV();
}

/*void initialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial)
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
}*/

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

void getMaterialParameters(libconfig::Setting& materialname, ParameterStruct& parameters, int m)
{
    using namespace libconfig;

    std::string EOSstring;
    std::string material1string;

    int tempa;

    materialname[m].lookupValue("material",material1string);

    materialname[m].lookupValue("mixture" ,tempa);

    parameters.materialInfo[m].mixture = tempa;

    materialname[m].lookupValue("EOS",EOSstring);

    Vector<Real> temp;

    temp.push_back(materialname[m]["adiabaticIndex"]);
    temp.push_back(materialname[m]["pref"]);
    temp.push_back(materialname[m]["eref"]);
    temp.push_back(materialname[m]["CV"]);

    if(EOSstring == "RomenskiiSolid")
    {
        parameters.materialInfo[m].EOS = new RomenskiiSolidEOS();

        temp.push_back(materialname[m]["rho0"]);
        temp.push_back(materialname[m]["K0"]);
        temp.push_back(materialname[m]["EOSalpha"]);
        temp.push_back(materialname[m]["EOSbeta"]);

        if(parameters.materialInfo[m].mixture)
        {
            std::cout << "Error: solid can't be a mixture (atm)" << std::endl;
            exit(1);
        }
    }
    /*else if(EOSstring == "MieGruneisenSolid")
    {
        parameters.materialInfo[m].EOS = new MieGruneisenSolidEOS();

        temp.push_back(materialname[m]["rho0"]);

        if(parameters.materialInfo[m].mixture)
        {
            std::cout << "Error: solid can't be a mixture (atm)" << std::endl;
            exit(1);
        }
    }*/
    else if(EOSstring == "MieGruneisen")
    {
        if(parameters.materialInfo[m].mixture)
        {
            parameters.materialInfo[m].EOS = new MixtureEOS();

            parameters.materialInfo[m].mixtureIndex = 0;

            temp.push_back(materialname[m]["adiabaticIndexmix"]);
            temp.push_back(materialname[m]["prefmix"]);
            temp.push_back(materialname[m]["erefmix"]);
            temp.push_back(materialname[m]["CVmix"]);
        }
        else
        {
            parameters.materialInfo[m].EOS = new MieGruneisenEOS();
        }
    }
    else
    {
        std::cout << "No valid EOS selected for material " << m << std::endl;
        exit(1);
    }

    if(material1string == "solid")
    {
        parameters.materialInfo[m].phase = solid;

        temp.push_back(materialname[m]["G0"]);

    }
    else if(material1string == "fluid")
    {
        parameters.materialInfo[m].phase = fluid;
    }

    parameters.materialInfo[m].EOS->define(temp);

    return;
}

void getState(libconfig::Setting& state, ParameterStruct& parameters, InitialStruct& initial, int s)
{
    using namespace libconfig;

    for(int m=0;m<parameters.numberOfMaterials;m++)
    {
        initial.alpha[s][m] = state[s]["material"][m]["alpha"];

        if(parameters.materialInfo[m].mixture)
        {
            initial.lambda[s][m] = state[s]["material"][m]["lambda"];
            initial.rho[s][m] = state[s]["material"][m]["rhoa"];
        }
        else if(parameters.materialInfo[m].phase == fluid)
        {
            initial.rho[s][m] = state[s]["material"][m]["rho"];
        }
        else
        {
            for(int i=0;i<9;i++)
            {
                initial.F[s][i] = state[s]["material"][m]["F"][i];

            }
        }
    }

    initial.p[s] = state[s]["p"];
    initial.u[s] = state[s]["u"];
    initial.v[s] = state[s]["v"];
    initial.w[s] = state[s]["w"];

    return;
}

void libConfigInitialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial)
{
    using namespace libconfig;

    try
    {
        ParmParse pp;

        std::string settingsFileString;

        pp.get("SettingsFile",settingsFileString);

        char settingsFile[settingsFileString.length()+1];

        stringtochar(settingsFileString,settingsFile);

        Config cfg;

        std::cout << "File: " << settingsFile << std::endl;

        try
        {
            cfg.readFile(settingsFile);
        }
        catch(ParseException except)
        {
            std::cout << "Incorrect Settings file" << std::endl;
            exit(1);
        }

        /*****************************************************
         * General Simulation Parameters
         ****************************************************/

        pp.get("max_grid_size", parameters.max_grid_size);
        pp.get("Nghost",        parameters.Nghost);

        pp.get("plotDirectory", initial.filename);

        cfg.lookupValue("finalTime",initial.finalT);

        cfg.lookupValue("numberOfxCells",parameters.n_cells[0]);
        cfg.lookupValue("numberOfyCells",parameters.n_cells[1]);
        cfg.lookupValue("numberOfzCells",parameters.n_cells[2]);

        cfg.lookupValue("xdomainLength",parameters.dimL[0]);
        cfg.lookupValue("ydomainLength",parameters.dimL[1]);
        cfg.lookupValue("zdomainLength",parameters.dimL[2]);

        cfg.lookupValue("CFL",parameters.CFL);

        cfg.lookupValue("materials",    parameters.numberOfMaterials);
        cfg.lookupValue("mixtures",     parameters.numberOfMixtures);
        cfg.lookupValue("states",       initial.numberOfStates);

        initial.interfaces.resize(1);

        cfg.lookupValue("interface",initial.interfaces[0]);

        cfg.lookupValue("numberOfPictures",initial.numberOfPictures);

        cfg.lookupValue("SOLID",    parameters.SOLID);
        cfg.lookupValue("THINC",    parameters.THINC);
        cfg.lookupValue("REACTIVE", parameters.REACTIVE);
        cfg.lookupValue("RADIAL",   parameters.RADIAL);
        cfg.lookupValue("PLASTIC",  parameters.PLASTIC);
        cfg.lookupValue("MUSCL",    parameters.MUSCL);

        cfg.lookupValue("THINCbeta",parameters.THINCbeta);

        Setting& root = cfg.getRoot();

        for(int dir; dir<AMREX_SPACEDIM;dir++)
        {
            initial.lowBoundary[dir]  = root["lowBoundary"][dir];
            initial.highBoundary[dir] = root["highBoundary"][dir];
        }


        int m = parameters.numberOfMaterials;
        int mix = parameters.numberOfMixtures;

        parameters.Ncomp = ((2+2)*mix)+m+m+m+3+3+1+1+1+1+1+9+(9+9+9+1)*parameters.SOLID;

        /**********************************************************
         * Material Parameters
         **********************************************************/


        parameters.materialInfo.resize(parameters.numberOfMaterials);
        initial.resize(parameters);

        Setting* materials = &root["listOfMaterials"];

        for(int m=0;m<parameters.numberOfMaterials;m++)
        {
            getMaterialParameters(*materials,parameters,m);
        }

        Setting* states = &root["listOfStates"];

        for(int s=0;s<initial.numberOfStates;s++)
        {
            getState(*states,parameters,initial,s);
        }


    }
    catch( SettingException except)
    {
        std::cout << except.getPath() << std::endl;
    }

    return;
}
