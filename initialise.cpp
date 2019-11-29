#include "simulationheader.h"

void chooseStateBasedOnInitialCondition(int& s, int i, int j, int k, InitialStruct& initial, ParameterStruct& parameters);
void AMR_chooseStateBasedOnInitialCondition(int& s, Real x, Real y, Real z, InitialStruct& initial, ParameterStruct& parameters);
Real solidVolumeFractionWeight(int& s, Real x, Real y, Real z, InitialStruct& initial, ParameterStruct& parameters, const Real* dx);
Real densityWeight(int& s, Real x, Real y, Real z, InitialStruct& initial, ParameterStruct& parameters, const Real* dx, Real top, Real bottom);

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

void initial_conditions(BoxAccessCellArray& U, ParameterStruct& parameters, InitialStruct& initial,const Real* dx, const Real* prob_lo)
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

    Real volfrac;

    Real x,y,z;

    for 		(int k = lo.z; k <= hi.z; ++k)
    {
                z = prob_lo[2] + (Real(k)+0.5)*dx[2];

        for 	(int j = lo.y; j <= hi.y; ++j)
        {
                y = prob_lo[1] + (Real(j)+0.5)*dx[1];

            for (int i = lo.x; i <= hi.x; ++i)
            {
                x = prob_lo[0] + (Real(i)+0.5)*dx[0];


                AMR_chooseStateBasedOnInitialCondition(s,x,y,z,initial,parameters);

                for (int m = 0; m < parameters.numberOfMaterials; m++)
                {
                    U(i,j,k,RHO,m)    = initial.rho[s][m];

                    U(i,j,k,P,m)      = initial.p[s][m];

                    U(i,j,k,VELOCITY,m,0)  = initial.u[s][m];
                    U(i,j,k,VELOCITY,m,1)  = initial.v[s][m];
                    U(i,j,k,VELOCITY,m,2)  = initial.w[s][m];

                    U(i,j,k,P,m) += U.getEffectiveNonThermalPressure(i,j,k,m)/U.getEffectiveInverseGruneisen(i,j,k,m);
                }

                //U(i,j,k,RHO_K,0) = densityWeight(s,x,y,z,initial,parameters,dx,1.0,0.125);
                //U(i,j,k,P) = densityWeight(s,x,y,z,initial,parameters,dx,1.0,0.1);

                // Need for UdaykumarGroove Test
                //U(i,j,k,ALPHA,0) = solidVolumeFractionWeight(s,x,y,z,initial,parameters,dx);
                //U(i,j,k,ALPHA,1) = 1.0-U(i,j,k,ALPHA,0);

                //U(i,j,k,P)             = initial.p[s];
            }
        }
    }
}

void setInitialConditions(CellArray& U, ParameterStruct& parameters, InitialStruct& initial,const Real* dx, const Real* prob_lo)
{
    U.data.setVal(0.0);

    for(MFIter mfi(U.data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,U);

        initial_conditions(baca, parameters, initial,dx,prob_lo);

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

    /*if(EOSstring == "RomenskiiSolid")
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
    else */if(EOSstring == "SINGLEMATERIAL_MieGruneisen")
    {
        parameters.materialInfo[m].EOS = new SINGLEMATERIAL_MieGruneisenEOS();
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
        initial.rho[s][m] = state[s]["material"][m]["rho"];
        initial.p  [s][m] = state[s]["material"][m]["p"];
        initial.u  [s][m] = state[s]["material"][m]["u"];
        initial.v  [s][m] = state[s]["material"][m]["v"];
        initial.w  [s][m] = state[s]["material"][m]["w"];

    }
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

        //std::cout << "File: " << settingsFile << std::endl;

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

        //pp.get("max_grid_size", parameters.max_grid_size);
        //pp.get("Nghost",        parameters.Nghost);

        //pp.get("plotDirectory", initial.filename);


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
        Print() << except.getPath() << std::endl;
    }

    return;
}

void chooseStateBasedOnInitialCondition(int& s, int i, int j, int k, InitialStruct& initial, ParameterStruct& parameters)
{

    /******************************************
     * 1D RP
     *****************************************/
    {
        Print() << i << " " << (int)((initial.interface/parameters.dimL[0])*parameters.n_cells[0]) << std::endl;
        if(i < (int)((initial.interface/parameters.dimL[0])*parameters.n_cells[0]))
        {
            s=0;
        }
        else
        {
            Print() << "HERE" << std::endl;
            s=1;
        }
    }

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
    /*{
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
    }*/

    /******************************************
     * Udaykunar Groove
     *****************************************/
    /*{
        Real x = i*parameters.dx[0];
        Real y = j*parameters.dx[1];

        Real shock = 1E-3;
        Real interface = initial.interface;
        Real radius = 15E-3;

        //Real shock = std::max((int)((1E-3/parameters.dimL[1])*parameters.n_cells[1]),5);
        //Real interface =      (int)((initial.interface/parameters.dimL[1])*parameters.n_cells[1]);
        //Real radius =         (int)((15E-3/parameters.dimL[0])*parameters.n_cells[0]);

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
    }*/



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

void AMR_chooseStateBasedOnInitialCondition(int& s, Real x, Real y, Real z, InitialStruct& initial, ParameterStruct& parameters)
{

    /******************************************
     * 1D RP
     *****************************************/
    {
        if(x + y < initial.interface)
        {
            s=0;
        }
        else
        {
             s=1;
        }
    }

    /******************************************
     * 2D Sod
     *****************************************/
    /*{
        if((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) < 0.2*0.2)
        {
            s=0;
        }
        else
        {
             s=1;
        }
    }*/

    /******************************************
     * RateStick
     *****************************************/

    /*{
        Real radius       = initial.interface;
        Real startOfTube  = initial.interface;
        Real endOfBooster = 2.0*initial.interface;
        Real chamfer      = 0.2*initial.interface;

        if(y<startOfTube ||  x > radius)
        {
            s=2;
        }
        else if(y>startOfTube+chamfer && y <= endOfBooster)
        {
            s=0;
        }
        else if(y<=startOfTube+chamfer)
        {
            if(x<radius-chamfer)
            {
                s=0;
            }
            else if((x-(radius-chamfer))*(x-(radius-chamfer))+(y-(startOfTube+chamfer))*(y-(startOfTube+chamfer)) < chamfer*chamfer)
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
     * Udaykunar Groove
     *****************************************/
    /*{

        Real shock = 0.5E-3;
        Real interface = initial.interface;
        Real radius = 4E-3; //15E-3;

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
    }*/

    /******************************************
     * Udaykunar Groove 2D
     *****************************************/
    /*{

        Real shock = initial.interface-16E-3;
        Real interface = initial.interface;
        Real radius = 15E-3;
        Real middle = 25E-3;


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
    }*/


}

Real solidVolumeFractionWeight(int& s, Real x, Real y, Real z, InitialStruct& initial, ParameterStruct& parameters, const Real* dx)
{

        Real shock = 0.5E-3;
        Real interface = initial.interface;
        Real radius = 4E-3; //15E-3;

        int sub = 11;
        int counter = 0.0;

        if(y < shock)
        {
            return 0.999999;
        }
        else if(y < interface)
        {
            for(int row = -(sub-1)/2; row< (sub+1)/2; row++)
            {
                for(int col = -(sub-1)/2; col < (sub+1)/2; col++)
                {
                    if( (y+((Real)col)/((Real)sub)*dx[1]-interface)*(y+((Real)col)/((Real)sub)*dx[1]-interface) + (x+((Real)row/((Real)sub))*dx[0])*(x+((Real)row/((Real)sub))*dx[0]) < radius*radius  )
                    {
                        counter++;
                    }
                }
            }

            return 0.999998*(1.0 - ((Real) counter)/((Real) sub*sub))+0.000001;
        }
        else
        {
            return 0.000001;
        }
}

Real densityWeight(int& s, Real x, Real y, Real z, InitialStruct& initial, ParameterStruct& parameters, const Real* dx, Real top, Real bottom)
{

        Real interface = initial.interface;
        Real radius = 0.2;

        int sub = 11;
        int counter = 0.0;


        for(int row = -(sub-1)/2; row< (sub+1)/2; row++)
        {
            for(int col = -(sub-1)/2; col < (sub+1)/2; col++)
            {
                if( (y+((Real)col)/((Real)sub)*dx[1]-interface)*(y+((Real)col)/((Real)sub)*dx[1]-interface) + (x+((Real)row/((Real)sub))*dx[0]-interface)*(x+((Real)row/((Real)sub))*dx[0]-interface) < radius*radius  )
                {
                    counter++;
                }
            }
        }

        return (top-bottom)*(((Real) counter)/((Real) sub*sub))+bottom;
}
