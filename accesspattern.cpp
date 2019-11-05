#include "simulationheader.h"

void AccessPattern::addVariable(int& position, std::string nameBase, Var_type type, Var_type INCELL, Variable var, int materialNumber=1, int rowNumber=1, int colNumber=1)
{
    if(materialNumber*rowNumber*colNumber > 0)
    {
        data.insert(std::pair<Variable,int>(var,position));

        position += materialNumber*rowNumber*colNumber;

        for(int m = 0; m < materialNumber ;m++)
        {
            for(int row = 0; row < rowNumber ; row++)
            {
                for(int col = 0; col < colNumber ; col++)
                {
                    variableNames.push_back((nameBase + std::to_string(m)+ "_" + std::to_string(row) + std::to_string(col)));

                    switch(type)
                    {
                    case CONSERVATIVE:  conservativeVariables.push_back(MaterialSpecifier(var,m,row,col));
                        break;
                    case PRIMITIVE:     primitiveVariables.push_back(MaterialSpecifier(var,m,row,col));
                        break;
                    case BOTH:          conservativeVariables.push_back(MaterialSpecifier(var,m,row,col));
                                        primitiveVariables.push_back(MaterialSpecifier(var,m,row,col));
                        break;
                    case NEITHER:
                        break;
                    default: Print() << "Didn't specifiy variable type in access pattern" << std::endl; exit(1);
                    }

                    if(INCELL == CELL)
                    {
                        cellVariables.push_back(MaterialSpecifier(var,m,row,col));
                    }

                    allVariables.push_back(MaterialSpecifier(var,m,row,col));


                }
            }
        }
    }
}

AccessPattern::AccessPattern(ParameterStruct& parameters) : materialInfo(parameters.materialInfo)
{
    define(parameters);
}

void AccessPattern::define(ParameterStruct& parameters)
{
    int n=0;

    data.clear();
    conservativeVariables.clear();
    primitiveVariables.clear();
    allVariables.clear();
    variableNames.clear();
    cellVariables.clear();

    addVariable(n,"alpha"   ,       BOTH,           CELL,    ALPHA,          parameters.numberOfMaterials);
    addVariable(n,"alphaRho",       CONSERVATIVE,   CELL,    ALPHARHO,       parameters.numberOfMaterials);
    addVariable(n,"rho_k",          PRIMITIVE,      CELL,    RHO_K,          parameters.numberOfMaterials);
    addVariable(n,"rho",            NEITHER,        CELL,    RHO);
    addVariable(n,"rhoU",           CONSERVATIVE,   CELL,    RHOU,           1,3);
    addVariable(n,"u",              PRIMITIVE,      CELL,    VELOCITY,       1,3);
    addVariable(n,"E",              CONSERVATIVE,   CELL,    TOTAL_E);
    addVariable(n,"p",              PRIMITIVE,      CELL,    P);
    addVariable(n,"a",              NEITHER,        CELL,    SOUNDSPEED);
    addVariable(n,"uStar",          NEITHER,        CELL,    USTAR);
    addVariable(n,"sigma",          NEITHER,        CELL,    SIGMA,          1,3,3);


    addVariable(n,"rhomix",         NEITHER,        NOTCELL,    RHO_MIX,        parameters.numberOfMixtures,2);
    addVariable(n,"lambda",         PRIMITIVE,      CELL,       LAMBDA,         parameters.numberOfMixtures);
    addVariable(n,"alpharholambda", CONSERVATIVE,   CELL,       ALPHARHOLAMBDA, parameters.numberOfMixtures);

    if(parameters.SOLID)
    {
        addVariable(n,"V",          BOTH,           CELL,    V_TENSOR,       1,3,3);
        addVariable(n,"VStar",      NEITHER,        CELL,    VSTAR,          1,3,3);
        addVariable(n,"devH",       NEITHER,        NOTCELL, DEVH,           1,3,3);
        addVariable(n,"HenckyJ2",   NEITHER,        NOTCELL, HJ2);

        if(parameters.PLASTIC)
        {
            addVariable(n,"epsilon",            PRIMITIVE,    CELL,  EPSILON,         parameters.numberOfMaterials);
            addVariable(n,"alphaRhoEpsilon",    CONSERVATIVE, CELL,  ALPHARHOEPSILON, parameters.numberOfMaterials);
        }
    }

}


int& AccessPattern::operator[](Variable var)
{
    return data[var];
}
