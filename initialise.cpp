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
/*
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

    int radiusInt       = (int)((initial.interfaces[0]/parameters.dimL[0])*parameters.n_cells[0]);
    int startOfTubeInt  = (int)((initial.interfaces[0]/parameters.dimL[1])*parameters.n_cells[1]);
    int endOfBoosterInt = (int)((2.0*initial.interfaces[0]/parameters.dimL[1])*parameters.n_cells[1]);
    int chamfer         = (int)((0.2*initial.interfaces[0]/parameters.dimL[1])*parameters.n_cells[1]);
    int solidThickness  = (int)((0.5*initial.interfaces[0]/parameters.dimL[1])*parameters.n_cells[1]);

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
                if(j<startOfTubeInt-solidThickness ||  i> radiusInt+solidThickness)
                {
                    s=3;
                }
                else if(j<startOfTubeInt || i > radiusInt )
                {
                    if(i<radiusInt+solidThickness-chamfer)
                    {
                        s=2;
                    }
                    else if(j>startOfTubeInt-solidThickness+chamfer)
                    {
                        s=2;
                    }
                    else if( (i-(radiusInt+solidThickness-chamfer))*(i-(radiusInt+solidThickness-chamfer))+(j-(startOfTubeInt-solidThickness+chamfer))*(j-(startOfTubeInt-solidThickness+chamfer))< chamfer*chamfer   )
                    {
                        s=2;
                    }
                    else
                    {
                        s=3;
                    }
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
            }
        }
    }
}

void Schoch_Can_initial_conditions(BoxAccessCellArray& U, ParameterStruct& parameters, InitialStruct& initial)
{

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    int s = 0;

    int extRadius       = (int)((14.0E-2/parameters.dimL[0])*parameters.n_cells[0]);
    int extLength       = (int)((25.0E-2/parameters.dimL[1])*parameters.n_cells[1]);
    int boosterThickness= (int)(( 1.0E-2/parameters.dimL[1])*parameters.n_cells[1]);
    int chamfer         = (int)(( 1.2E-2/parameters.dimL[1])*parameters.n_cells[1]);
    int wallThickness   = (int)((   2E-2/parameters.dimL[1])*parameters.n_cells[1]);
    int sideSeparation  = (int)((((parameters.dimL[1] - 25.0E-2)/2.0)/parameters.dimL[1])*parameters.n_cells[1]);


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

                if(j<sideSeparation || j>parameters.n_cells[1]-sideSeparation || i>extRadius)
                {
                    s=0; //air outside
                }
                else if( j>=(sideSeparation+wallThickness+chamfer) && j<parameters.n_cells[1]-(sideSeparation+wallThickness+chamfer) && i<extRadius-(wallThickness) )
                {
                    if(j<sideSeparation+wallThickness+boosterThickness)
                    {
                       s=3; //booster beyond chamfer
                    }
                    else
                    {
                        s=2; //explosive between chamfers
                    }
                }
                else if( j>sideSeparation+wallThickness && j<=sideSeparation+wallThickness+chamfer && i<extRadius-(wallThickness+chamfer))
                {
                    s=3; //booster below chamfer
                }
                else if( j<parameters.n_cells[1]-(sideSeparation+wallThickness) && j>=parameters.n_cells[1]-(sideSeparation+wallThickness+chamfer) && i<extRadius-(wallThickness+chamfer))
                {
                    s=2; //right explosive below chamfer
                }
                else if( (j-(sideSeparation+wallThickness+chamfer))*(j-(sideSeparation+wallThickness+chamfer))+(i-(extRadius-(wallThickness+chamfer)))*(i-(extRadius-(wallThickness+chamfer)))<chamfer*chamfer )
                {
                    s=3; //left chamfer
                }
                else if( (j-(parameters.n_cells[1]-(sideSeparation+wallThickness+chamfer)))*(j-(parameters.n_cells[1]-(sideSeparation+wallThickness+chamfer)))+(i-(extRadius-(wallThickness+chamfer)))*(i-(extRadius-(wallThickness+chamfer)))<chamfer*chamfer )
                {
                    s=2; //right chamfer
                }
                else if( i<(extRadius-chamfer) || (j >sideSeparation+chamfer && j < parameters.n_cells[1]-(sideSeparation+chamfer)))
                {
                    s=1; //Copper
                }
                else if( (j-(sideSeparation+chamfer))*(j-(sideSeparation+chamfer))+ (i-(extRadius-chamfer))*(i-(extRadius-chamfer)) < chamfer*chamfer || (j-(parameters.n_cells[1]-(sideSeparation+chamfer)))*(j-(parameters.n_cells[1]-(sideSeparation+chamfer)))+ (i-(extRadius-chamfer))*(i-(extRadius-chamfer)) < chamfer*chamfer      )
                {
                    s=1;
                }
                else
                {
                    s=0;
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
            }
        }
    }
}
*/

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
                    std::cout << std::setw(10) << V[s][row*U.numberOfComponents+col] << " ";
                }

                std::cout << std::endl ;
            }

            std::cout << std::endl ;
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

                for (int m = 0; m < parameters.numberOfMaterials; m++)
                {
                    U(i,j,k,ALPHA,m)    = initial.alpha[s][m];
                    U(i,j,k,RHO_K,m)    = initial.rho[s][m];

                    if(parameters.numberOfMixtures>0)
                    {
                        if(U.accessPattern.materialInfo[m].mixture)
                        {
                            U(i,j,k,RHO_MIX,m,0)    = initial.rho[s][m];
                            U(i,j,k,RHO_MIX,m,1)    = initial.rho[s][m];
                            U(i,j,k,LAMBDA,m)       = initial.lambda[s][m];

                        }
                        else
                        {
                            U(i,j,k,RHO_MIX,m,0)    = 0.0;
                            U(i,j,k,RHO_MIX,m,1)    = 0.0;
                            U(i,j,k,LAMBDA,m)       = 0.0;
                        }
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
/*#ifdef _OPENMP
#pragma omp parallel
#endif*/
    for(MFIter mfi(U.data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray baca(mfi,bx,U);

        initial_conditions(baca, parameters, initial);

        baca.primitiveToConservative();

        std::cout << std::setprecision(5);

        std::cout << std::setw(10) << baca(0,0,0,HJ2) << std::endl;
        std::cout << std::setw(10) << baca(499,0,0,HJ2) << std::endl;


        std::cout << "F" << std::endl;

        for(int row=0;row<baca.numberOfComponents;row++)
        {
            for(int col=0;col<baca.numberOfComponents;col++)
            {

                std::cout << std::setw(10) << initial.F[0][row*3+col] << " ";
            }

            std::cout << std::endl ;
        }

        std::cout << std::endl ;

        std::cout << "F" << std::endl;

        for(int row=0;row<baca.numberOfComponents;row++)
        {
            for(int col=0;col<baca.numberOfComponents;col++)
            {

               std::cout << std::setw(10) << initial.F[1][row*3+col] << " ";
            }

            std::cout << std::endl ;
        }

        std::cout << std::endl ;




        std::cout << "V" << std::endl;

        for(int row=0;row<baca.numberOfComponents;row++)
        {
            for(int col=0;col<baca.numberOfComponents;col++)
            {

                std::cout << std::setw(10) << baca(0,0,0,V_TENSOR,0,row,col) << " ";
            }

            std::cout << std::endl ;
        }

        std::cout << std::endl ;

        std::cout << "V" << std::endl;

        for(int row=0;row<baca.numberOfComponents;row++)
        {
            for(int col=0;col<baca.numberOfComponents;col++)
            {

                std::cout << std::setw(10) << baca(499,0,0,V_TENSOR,0,row,col) << " ";
            }

            std::cout << std::endl ;
        }

        std::cout << std::endl ;

        std::cout << "dev" << std::endl;

        for(int row=0;row<baca.numberOfComponents;row++)
        {
            for(int col=0;col<baca.numberOfComponents;col++)
            {

                std::cout << std::setw(10) << baca(0,0,0,DEVH,0,row,col) << " ";
            }

            std::cout << std::endl ;
        }

        std::cout << std::endl ;

        for(int row=0;row<baca.numberOfComponents;row++)
        {
            for(int col=0;col<baca.numberOfComponents;col++)
            {

                std::cout << std::setw(10) << baca(499,0,0,DEVH,0,row,col) << " ";
            }

            std::cout << std::endl ;
        }

        std::cout << std::endl ;


        std::cout << "sigma" << std::endl;

        for(int row=0;row<baca.numberOfComponents;row++)
        {
            for(int col=0;col<baca.numberOfComponents;col++)
            {

                std::cout << std::setw(10) << baca(0,0,0,SIGMA,0,row,col) << " ";
            }

            std::cout << std::endl ;
        }

        std::cout << std::endl ;

        for(int row=0;row<baca.numberOfComponents;row++)
        {
            for(int col=0;col<baca.numberOfComponents;col++)
            {

                std::cout << std::setw(10) << baca(499,0,0,SIGMA,0,row,col) << " ";
            }

            std::cout << std::endl ;
        }

        std::cout << std::endl ;




    }
}

void getMaterialParameters(libconfig::Setting& materialname, ParameterStruct& parameters, PlasticEOS& plastic, int m)
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

        //Print() << "HEre" << std::endl;

        if(parameters.materialInfo[m].mixture)
        {
            std::cout << "Error: solid can't be a mixture (atm)" << std::endl;
            exit(1);
        }
    }
    /*else if(EOSstring == "WilkinsSolid")
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
    }*/
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
        /*if(parameters.materialInfo[m].mixture)
        {
            parameters.materialInfo[m].EOS = new MixtureEOS();

            temp.push_back(materialname[m]["adiabaticIndexmix"]);
            temp.push_back(materialname[m]["prefmix"]);
            temp.push_back(materialname[m]["erefmix"]);
            temp.push_back(materialname[m]["CVmix"]);
        }
        else*/
        {
            parameters.materialInfo[m].EOS = new MieGruneisenEOS();

            //Print() << "HEre" << std::endl;
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

        if(parameters.PLASTIC)
        {
            parameters.materialInfo[m].plastic = true;
        }
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

void libConfigInitialiseDataStructs(ParameterStruct& parameters, InitialStruct& initial, PlasticEOS& plastic)
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

        //parameters.Ncomp = ((2+2)*mix)+m+m+m+3+3+1+1+1+1+1+9+(9+9+9+1+(m+m)*parameters.PLASTIC)*parameters.SOLID+(m+m)*parameters.PLASTIC;

        /**********************************************************
         * Material Parameters
         **********************************************************/


        parameters.materialInfo.resize(parameters.numberOfMaterials);
        initial.resize(parameters);
        plastic.yieldStress.resize(parameters.numberOfMaterials);

        Setting* materials = &root["listOfMaterials"];

        for(int m=0;m<parameters.numberOfMaterials;m++)
        {
            getMaterialParameters(*materials,parameters,plastic,m);
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
        exit(1);
    }

    return;
}

void chooseStateBasedOnInitialCondition(int& s, int i, int j, int k, InitialStruct& initial, ParameterStruct& parameters)
{

    /******************************************
     * 1D RP
     *****************************************/
    {
        if(i < (int)((initial.interface/parameters.dimL[0])*parameters.n_cells[0]))
        {
            s=0;
        }
        else
        {
            s=1;
        }
    }

    /******************************************
     * Wilkins
     *****************************************/
    /*{
        if(j < (int)((initial.interface/parameters.dimL[1])*parameters.n_cells[1]))
        {
            s=0;
        }
        else if(j < (int)((2.0*initial.interface/parameters.dimL[1])*parameters.n_cells[1]))
        {
            s=1;
        }
        else
        {
            s=2;
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

