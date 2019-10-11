#include "simulationheader.h"

AccessPattern::AccessPattern(ParameterStruct& parameters)
{
    int n=0;

    data.insert(std::pair<Variable,int>(ALPHA,n));

    n+=parameters.numberOfMaterials;

    data.insert(std::pair<Variable,int>(ALPHARHO,n));

    n+=parameters.numberOfMaterials;

    data.insert(std::pair<Variable,int>(RHO_K,n));

    n+=parameters.numberOfMaterials;

    data.insert(std::pair<Variable,int>(RHO,n));

    n+=1;

    data.insert(std::pair<Variable,int>(RHOU,n));

    n+=3;

    data.insert(std::pair<Variable,int>(VELOCITY,n));

    n+=3;

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
        conservativeVariables.push_back(MaterialSpecifier(RHOU,row));
        primitiveVariables.push_back(MaterialSpecifier(VELOCITY,row));
    }

}

int& AccessPattern::operator[](Variable var)
{
    return data[var];
}
