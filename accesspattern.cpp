#include "simulationheader.h"

void AccessPattern::addVariable(int& position, std::string nameBase, Var_type type, Var_type INREFINE, Real lo, Real hi, Variable var, int m, int componentNumber=1, int rowNumber=1, int colNumber=1)
{
    if(componentNumber*rowNumber*colNumber > 0)
    {

        data[m][var] = position;

        numberOfComponentsForVariable[m][var]= ( (componentNumber > 1 && rowNumber*colNumber > 1) ? rowNumber*colNumber : 1);
        numberOfRowsForVariable[m][var]      = ( (rowNumber       > 1 &&           colNumber > 1) ?           colNumber : 1);

        limits.insert(std::pair<Variable,std::pair<Real,Real> >(var,std::make_pair(lo,hi)));



        position += componentNumber*rowNumber*colNumber;

        for(int c = 0; c < componentNumber;c++)
        {
            for(int row = 0; row < rowNumber ; row++)
            {
                for(int col = 0; col < colNumber ; col++)
                {
                    variableNames.push_back((nameBase + std::to_string(m)+ "_" + std::to_string(c)+ "_" + std::to_string(row) + std::to_string(col)));

                    switch(type)
                    {
                    case CONSERVATIVE:  conservativeVariables[m].push_back(MaterialSpecifier(var,m,c,row,col));
                        break;
                    case PRIMITIVE:     primitiveVariables[m].push_back(MaterialSpecifier(var,m,c,row,col));
                        break;
                    case BOTH:          conservativeVariables[m].push_back(MaterialSpecifier(var,m,c,row,col));
                                        primitiveVariables[m].push_back(MaterialSpecifier(var,m,c,row,col));
                        break;
                    case NEITHER:
                        break;
                    default: Abort("Didn't specifiy variable type in access pattern");
                    }

                    if(INREFINE == REFINE)
                    {
                        refineVariables.push_back(MaterialSpecifier(var,m,c,row,col));
                    }

                    allVariables.push_back(MaterialSpecifier(var,m,c,row,col));
                }
            }
        }
    }
}

AccessPattern::AccessPattern(ParameterStruct& parameters) : materialInfo(parameters.materialInfo), interface(parameters.interfaceType)
{
    //define(parameters);
}

void AccessPattern::define(ParameterStruct& parameters)
{
    int n=0;

    data.clear();
    limits.clear();
    conservativeVariables.clear();
    primitiveVariables.clear();
    allVariables.clear();
    refineVariables.clear();
    variableNames.clear();
    numberOfComponentsForVariable.clear();
    numberOfRowsForVariable.clear();


    Real c   = 3E8;
    Real min = 1E-20;
    Real max = 1E20;

    data.                         resize(parameters.numberOfSharpMaterials);
    numberOfComponentsForVariable.resize(parameters.numberOfSharpMaterials);
    numberOfRowsForVariable.      resize(parameters.numberOfSharpMaterials);
    conservativeVariables.        resize(parameters.numberOfSharpMaterials);
    primitiveVariables.           resize(parameters.numberOfSharpMaterials);


    for(int m = 0; m < parameters.numberOfSharpMaterials; m++)
    {
        data[m].resize(50);

        numberOfComponentsForVariable[m].resize(50);
        numberOfRowsForVariable[m].resize(50);

        if(parameters.interfaceType[m] == SHARP)
        {
            Print() << "Here " << m << std::endl;

           addVariable(n,"rho",   BOTH,         REFINE,  min,  max, RHO,        m);
           addVariable(n,"u",     PRIMITIVE,    NEITHER, -c,   c,   VELOCITY,   m, 1, 3);
           addVariable(n,"rhoU",  CONSERVATIVE, NEITHER, -c,   c,   RHOU,       m, 1, 3);
           addVariable(n,"p",     PRIMITIVE,    NEITHER, -max, max, P,          m);
           addVariable(n,"E",     CONSERVATIVE, NEITHER, -max, max, TOTAL_E,    m);
           addVariable(n,"a",     NEITHER,      NEITHER,  min, c,   SOUNDSPEED, m);
           addVariable(n,"sigma", NEITHER,      NEITHER, -max, max, SIGMA,      m, 1, 3, 3);

        }
        else
        {
            Print() << "Here " << m << std::endl;

            addVariable(n,"alpha"   , BOTH,         REFINE,   min,    1.0,  ALPHA,      m,     parameters.numberOfDiffuseMaterials[m]);
            addVariable(n,"alphaRho", CONSERVATIVE, REFINE,   min,    max,  ALPHARHO,   m,     parameters.numberOfDiffuseMaterials[m]);
            addVariable(n,"rho_k",    PRIMITIVE,    NEITHER,  min,    max,  RHO_K,      m,     parameters.numberOfDiffuseMaterials[m]);
            addVariable(n,"rho",      NEITHER,      NEITHER,  min,    max,  RHO,        m);
            addVariable(n,"rhoU",     CONSERVATIVE, NEITHER, -max,    max,  RHOU,       m,     1,3);
            addVariable(n,"u",        PRIMITIVE,    NEITHER ,-c ,     c,    VELOCITY,   m,     1,3);
            addVariable(n,"E",        CONSERVATIVE, NEITHER, -max,    max,  TOTAL_E,    m);
            addVariable(n,"p",        PRIMITIVE,    REFINE,  -max,    max,  P,          m);
            addVariable(n,"a",        NEITHER,      NEITHER,  min,    c,    SOUNDSPEED, m);
            addVariable(n,"sigma",    NEITHER,      NEITHER, -max,    max,  SIGMA,      m,     1,3,3);

        }
    }
}

int& AccessPattern::operator[](MaterialSpecifier& m)
{
    return data[m.mat][m.var];
}
