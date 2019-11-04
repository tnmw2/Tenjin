#include "simulationheader.h"

void chooseStateBasedOnInitialCondition(int& s, int i, int j, int k, InitialStruct& initial, ParameterStruct& parameters);

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


void initial_conditions(BoxAccessCellArray& U, ParameterStruct& parameters, InitialStruct& initial)
{

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    std::vector< std::vector<Real> > V(initial.numberOfStates);

    for(int s = 0; s < initial.numberOfStates; s++)
    {
        if(parameters.SOLID)
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

    int s=0;

    for                 (int k = lo.z; k <= hi.z; ++k)
    {
        for             (int j = lo.y; j <= hi.y; ++j)
        {
            for         (int i = lo.x; i <= hi.x; ++i)
            {
                chooseStateBasedOnInitialCondition(s,i,j,k,initial,parameters);

                //Print()<< s << " ";

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

                    if(parameters.PLASTIC)
                    {
                        U(i,j,k,EPSILON,m) = 0.0;
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

                U(i,j,k,P) += U.getEffectiveNonThermalPressure(i,j,k)/U.getEffectiveInverseGruneisen(i,j,k);

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

void getMaterialParameters(libconfig::Setting& materialname, ParameterStruct& parameters, int m, PlasticEOS& plastic)
{
    using namespace libconfig;

    std::string EOSstring;
    std::string material1string;

    materialname[m].lookupValue("material",material1string);

    int mix;

    materialname[m].lookupValue("mixture" ,mix);

    parameters.materialInfo[m].mixture = mix;

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
    else if(EOSstring == "WilkinsSolid")
    {
        parameters.materialInfo[m].EOS = new WilkinsSolidEOS();

        temp.push_back(materialname[m]["rho0"]);
        temp.push_back(materialname[m]["e1"]);
        temp.push_back(materialname[m]["e2"]);
        temp.push_back(materialname[m]["e3"]);
        temp.push_back(materialname[m]["e4"]);
        temp.push_back(materialname[m]["e5"]);

        if(parameters.materialInfo[m].mixture)
        {
            std::cout << "Error: solid can't be a mixture (atm)" << std::endl;
            exit(1);
        }
    }
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

    if(parameters.PLASTIC)
    {
        plastic.yieldStress[m] = 0.0;
    }

    if(material1string == "solid")
    {
        parameters.materialInfo[m].phase = solid;

        temp.push_back(materialname[m]["G0"]);

        if(parameters.PLASTIC)
        {
            plastic.yieldStress[m] = materialname[m]["yieldStress"];
            parameters.materialInfo[m].plastic = true;
        }

    }
    else if(material1string == "fluid")
    {
        parameters.materialInfo[m].phase = fluid;
        parameters.materialInfo[m].plastic = false;
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

void libConfigInitialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial, PlasticEOS &plastic)
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

        cfg.lookupValue("interface",initial.interface);

        cfg.lookupValue("numberOfPictures",initial.numberOfPictures);

        cfg.lookupValue("SOLID",    parameters.SOLID);
        cfg.lookupValue("THINC",    parameters.THINC);
        cfg.lookupValue("RADIAL",   parameters.RADIAL);
        cfg.lookupValue("PLASTIC",  parameters.PLASTIC);
        cfg.lookupValue("MUSCL",    parameters.MUSCL);
        cfg.lookupValue("REACTIVE", parameters.REACTIVE);

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
        plastic.yieldStress.resize(parameters.numberOfMaterials);


        Setting* materials = &root["listOfMaterials"];

        for(int m=0;m<parameters.numberOfMaterials;m++)
        {
            getMaterialParameters(*materials,parameters,m,plastic);
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

void chooseStateBasedOnInitialCondition(int& s, int i, int j, int k, InitialStruct& initial, ParameterStruct& parameters)
{

    /******************************************
     * 1D RP
     *****************************************/
    /*{
        if(i < (int)((initial.interface/parameters.dimL[0])*parameters.n_cells[0]))
        {
            s=0;
        }
        else
        {
            s=1;
        }
    }*/

    /******************************************
     * Wilkins
     *****************************************/
    /*{
        if(i < (int)((initial.interface/parameters.dimL[0])*parameters.n_cells[0]))
        {
            s=0;
        }
        else if(i < (int)((2.0*initial.interface/parameters.dimL[0])*parameters.n_cells[0]))
        {
            s=1;
        }
        else
        {
            s=2;
        }
    }*/


    /******************************************
     * Solid RMI - sine
     *****************************************/
    /*{
        Real amplitude = 2E-4;

        Real pi = 3.14;

        if(i < (int)((0.5*initial.interface/parameters.dimL[0])*parameters.n_cells[0]))
        {
            s=1;
        }
        else if(i < (int)(((0.1)*parameters.n_cells[0]))+ (int)((amplitude*sin(((double)j/(double)(2*parameters.n_cells[1]))*pi*2.0+pi/2.0)/parameters.dimL[0])*parameters.n_cells[0]))
        {
            s=2;
        }
        else
        {
            s=3;
        }
    }*/

    /******************************************
     * Solid RMI - saw
     *****************************************/
    {
        Real gradient = -0.2E-3;

        if(i < (int)((0.5*initial.interface/parameters.dimL[0])*parameters.n_cells[0]))
        {
            s=1;
        }
        else if(i < (int)(((0.1)*parameters.n_cells[0]))+ (int)((gradient*((double)j/(double)(parameters.n_cells[1]))/parameters.dimL[0])*parameters.n_cells[0]))
        {
            s=2;
        }
        else
        {
            s=3;
        }
    }

    /******************************************
     * Udaykunar Groove
     *****************************************/
    {
        Real x = i*parameters.dx[0];
        Real y = j*parameters.dx[1];

        Real shock = 1E-3;
        Real interface = initial.interface;
        Real radius = 15E-3;

        /*Real shock = std::max((int)((1E-3/parameters.dimL[1])*parameters.n_cells[1]),5);
        Real interface =      (int)((initial.interface/parameters.dimL[1])*parameters.n_cells[1]);
        Real radius =         (int)((15E-3/parameters.dimL[0])*parameters.n_cells[0]);*/

        if(y < shock)
        {
            s=1;
        }
        else if(y < interface)
        {
            if( (y-interface)*(y-interface) + (x)*(x) < radius*radius  )
            {
                s=3;
            }
            else
            {
                s=2;
            }
        }
        else
        {
            s=3;
        }
    }



    /******************************************
     * RateStick
     *****************************************/

    /*{
        int radiusInt       = (int)((initial.interface/parameters.dimL[0])*parameters.n_cells[0]);
        int startOfTubeInt  = (int)((initial.interface/parameters.dimL[1])*parameters.n_cells[1]);
        int endOfBoosterInt = (int)((2.0*initial.interface/parameters.dimL[1])*parameters.n_cells[1]);
        int chamfer         = (int)((0.2*initial.interface/parameters.dimL[1])*parameters.n_cells[1]);

        if(j<startOfTubeInt ||  i> radiusInt)
        {
            s=2;
        }
        else if(j>startOfTubeInt+chamfer && j <= endOfBoosterInt)
        {
            s=0;
        }
        else if(j<=startOfTubeInt+chamfer)
        {
            if(i<radiusInt-chamfer)
            {
                s=0;
            }
            else if((i-(radiusInt-chamfer))*(i-(radiusInt-chamfer))+(j-(startOfTubeInt+chamfer))*(j-(startOfTubeInt+chamfer)) < chamfer*chamfer)
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
    }*/


    /******************************************
     * Rod Impact
     *****************************************/

    /*{
        int chamfer = (int)((0.06E-2/parameters.dimL[0])*parameters.n_cells[0]);
        int length  = (int)((2.347E-2/parameters.dimL[1])*parameters.n_cells[1]);
        int radius  = (int)((initial.interface/parameters.dimL[0])*parameters.n_cells[0]);

        if(i<radius)
        {
            if(j< length-chamfer)
            {
                s=2;
            }
            else if(j<length)
            {
                if( (i< radius-chamfer) || (i - (radius-chamfer))*(i - (radius-chamfer))+(j - (length-chamfer))*(j - (length-chamfer)) < chamfer*chamfer     )
                {
                    s=2;
                }
                else
                {
                    s=1;
                }
            }
            else
            {
                s=1;
            }
        }
        else
        {
            s=0;
        }
    }*/


}
