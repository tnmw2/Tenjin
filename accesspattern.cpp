#include "simulationheader.h"

void AccessPattern::addVariable(int& position, std::string nameBase, Var_type type, Variable var, int materialNumber=1, int rowNumber=1, int colNumber=1)
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

AccessPattern::AccessPattern(ParameterStruct& parameters)
{
    int n=0;

    addVariable(n,"alpha"   ,   BOTH,           ALPHA,      parameters.numberOfMaterials);
    addVariable(n,"alphaRho",   CONSERVATIVE,   ALPHARHO,   parameters.numberOfMaterials);
    addVariable(n,"rho_k",      PRIMITIVE,      RHO_K,      parameters.numberOfMaterials);
    addVariable(n,"rho",        NEITHER,        RHO);
    addVariable(n,"rhoU",       CONSERVATIVE,   RHOU,       1,AMREX_SPACEDIM);
    addVariable(n,"u",          PRIMITIVE,      VELOCITY,   1,AMREX_SPACEDIM);
    addVariable(n,"E",          CONSERVATIVE,   TOTAL_E);
    addVariable(n,"p",          PRIMITIVE,      P);
    addVariable(n,"a",          NEITHER,        SOUNDSPEED);
    addVariable(n,"uStar",      NEITHER,        USTAR);

    /*data.insert(std::pair<Variable,int>(ALPHA,n));

    n+=parameters.numberOfMaterials;

    data.insert(std::pair<Variable,int>(ALPHARHO,n));

    n+=parameters.numberOfMaterials;

    data.insert(std::pair<Variable,int>(RHO_K,n));

    n+=parameters.numberOfMaterials;

    data.insert(std::pair<Variable,int>(RHO,n));

    n+=1;

    data.insert(std::pair<Variable,int>(RHOU,n));

    n+=AMREX_SPACEDIM;

    data.insert(std::pair<Variable,int>(VELOCITY,n));

    n+=AMREX_SPACEDIM;

    data.insert(std::pair<Variable,int>(TOTAL_E,n));

    n+=1;

    data.insert(std::pair<Variable,int>(P,n));

    n+=1;

    data.insert(std::pair<Variable,int>(SOUNDSPEED,n));

    n+=1;

    data.insert(std::pair<Variable,int>(USTAR,n));

    for(int m=0;m<parameters.numberOfMaterials;m++)
    {
        conservativeVariables.push_back(MaterialSpecifier(ALPHA,m));
        primitiveVariables.push_back(MaterialSpecifier(ALPHA,m));

        conservativeVariables.push_back(MaterialSpecifier(ALPHARHO,m));
        primitiveVariables.push_back(MaterialSpecifier(RHO_K,m));
    }

    conservativeVariables.push_back(MaterialSpecifier(TOTAL_E));
    primitiveVariables.push_back(MaterialSpecifier(P));

    for(int row=0;row < AMREX_SPACEDIM ; row++)
    {
        conservativeVariables.push_back(MaterialSpecifier(RHOU,0,row));
        primitiveVariables.push_back(MaterialSpecifier(VELOCITY,0,row));
    }*/

}

int& AccessPattern::operator[](Variable var)
{
    return data[var];
}
