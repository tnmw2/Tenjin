#include "simulationheader.h"

void AccessPattern::addVariable(int& position, std::string nameBase, Var_type type, Variable var, int materialNumber=1, int rowNumber=1, int colNumber=1)
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
                }
            }
        }
    }
}

AccessPattern::AccessPattern(ParameterStruct& parameters) : materialInfo(parameters.materialInfo)
{
    int n=0;

    addVariable(n,"alpha"   ,       BOTH,           ALPHA,          parameters.numberOfMaterials);
    addVariable(n,"alphaRho",       CONSERVATIVE,   ALPHARHO,       parameters.numberOfMaterials);
    addVariable(n,"rho_k",          PRIMITIVE,      RHO_K,          parameters.numberOfMaterials);
    addVariable(n,"rho",            NEITHER,        RHO);
    addVariable(n,"rhoU",           CONSERVATIVE,   RHOU,           1,3);
    addVariable(n,"u",              PRIMITIVE,      VELOCITY,       1,3);
    addVariable(n,"E",              CONSERVATIVE,   TOTAL_E);
    addVariable(n,"p",              PRIMITIVE,      P);
    addVariable(n,"a",              NEITHER,        SOUNDSPEED);
    addVariable(n,"uStar",          NEITHER,        USTAR);
    addVariable(n,"rhomix",         NEITHER,        RHO_MIX,        parameters.numberOfMixtures,2);
    addVariable(n,"lambda",         PRIMITIVE,      LAMBDA,         parameters.numberOfMixtures);
    addVariable(n,"alpharholambda", CONSERVATIVE,   ALPHARHOLAMBDA, parameters.numberOfMixtures);
    addVariable(n,"sigma",          NEITHER,        SIGMA,          1,3,3);
    addVariable(n,"V",              BOTH,           V_TENSOR,       1,3,3);
    addVariable(n,"VStar",          NEITHER,        VSTAR,          1,3,3);
    addVariable(n,"devH",           NEITHER,        DEVH,           1,3,3);
    addVariable(n,"HenckyJ2",       NEITHER,        HJ2);

}

int& AccessPattern::operator[](Variable var)
{
    return data[var];
}
