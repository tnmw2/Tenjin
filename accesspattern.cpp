#include "simulationheader.h"

void AccessPattern::addVariable(int& position, std::string nameBase, Var_type type, Var_type INCELL, Var_type INREFINE, Variable var, int materialNumber=1, int rowNumber=1, int colNumber=1)
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

                    if(INREFINE == REFINE)
                    {
                        refineVariables.push_back(MaterialSpecifier(var,m,row,col));
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
    refineVariables.clear();
    variableNames.clear();
    cellVariables.clear();

    addVariable(n,"alpha"   ,       BOTH,           CELL, REFINE,   ALPHA,          parameters.numberOfMaterials);
    addVariable(n,"alphaRho",       CONSERVATIVE,   CELL, REFINE,   ALPHARHO,       parameters.numberOfMaterials);
    addVariable(n,"rho_k",          PRIMITIVE,      CELL, NEITHER,  RHO_K,          parameters.numberOfMaterials);
    addVariable(n,"rho",            NEITHER,        CELL, NEITHER,  RHO);
    addVariable(n,"rhoU",           CONSERVATIVE,   CELL, NEITHER,  RHOU,           1,3);
    addVariable(n,"u",              PRIMITIVE,      CELL, NEITHER,  VELOCITY,       1,3);
    addVariable(n,"E",              CONSERVATIVE,   CELL, NEITHER,  TOTAL_E);
    addVariable(n,"p",              PRIMITIVE,      CELL, NEITHER,  P);
    addVariable(n,"a",              NEITHER,        CELL, NEITHER,  SOUNDSPEED);
    addVariable(n,"uStar",          NEITHER,        CELL, NEITHER,  USTAR);
    addVariable(n,"sigma",          NEITHER,        CELL, NEITHER,  SIGMA,          1,3,3);


    addVariable(n,"rhomix",         NEITHER,        NOTCELL,NEITHER,    RHO_MIX,        parameters.numberOfMixtures,2);
    addVariable(n,"lambda",         PRIMITIVE,      CELL,   NEITHER,    LAMBDA,         parameters.numberOfMixtures);
    addVariable(n,"alpharholambda", CONSERVATIVE,   CELL,   NEITHER,    ALPHARHOLAMBDA, parameters.numberOfMixtures);

    if(parameters.SOLID)
    {
        addVariable(n,"V",          BOTH,           CELL,    NEITHER, V_TENSOR,       1,3,3);
        addVariable(n,"VStar",      NEITHER,        CELL,    NEITHER, VSTAR,          1,3,3);
        addVariable(n,"devH",       NEITHER,        NOTCELL, NEITHER, DEVH,           1,3,3);
        addVariable(n,"HenckyJ2",   NEITHER,        NOTCELL, REFINE , HJ2);

        if(parameters.PLASTIC)
        {
            addVariable(n,"epsilon",            PRIMITIVE,    CELL, NEITHER,  EPSILON,         parameters.numberOfMaterials);
            addVariable(n,"alphaRhoEpsilon",    CONSERVATIVE, CELL, NEITHER,  ALPHARHOEPSILON, parameters.numberOfMaterials);
        }
    }

}

int& AccessPattern::operator[](Variable var)
{
    return data[var];
}
