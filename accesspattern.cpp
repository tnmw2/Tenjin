#include "simulationheader.h"

void AccessPattern::addVariable(int& position, std::string nameBase, Var_type type, Var_type INCELL, Var_type INREFINE, Real lo, Real hi, Variable var, int materialNumber=1, int rowNumber=1, int colNumber=1)
{
    if(materialNumber*rowNumber*colNumber > 0)
    {
        data.insert(std::pair<Variable,int>(var,position));

        limits.insert(std::pair<Variable,std::pair<Real,Real> >(var,std::make_pair(lo,hi)));

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
    limits.clear();
    conservativeVariables.clear();
    primitiveVariables.clear();
    allVariables.clear();
    refineVariables.clear();
    variableNames.clear();
    cellVariables.clear();

    Real c   = 3E8;
    Real min = 1E-20;
    Real max = 1E20;




    addVariable(n,"alpha"   ,       BOTH,           CELL, REFINE,   min,    1.0,  ALPHA,          parameters.numberOfMaterials);
    addVariable(n,"alphaRho",       CONSERVATIVE,   CELL, REFINE,   min,    max,  ALPHARHO,       parameters.numberOfMaterials);
    addVariable(n,"rho_k",          PRIMITIVE,      CELL, NEITHER,  min,    max,  RHO_K,          parameters.numberOfMaterials);
    addVariable(n,"rho",            NEITHER,        CELL, NEITHER,  min,    max,  RHO);
    addVariable(n,"rhoU",           CONSERVATIVE,   CELL, NEITHER,  -max,   max,  RHOU,           1,3);
    addVariable(n,"u",              PRIMITIVE,      CELL, NEITHER,   -c,     c,    VELOCITY,       1,3);
    addVariable(n,"E",              CONSERVATIVE,   CELL, NEITHER,  min,    max,  TOTAL_E);
    addVariable(n,"p",              PRIMITIVE,      CELL, NEITHER,   -max,   max,  P);
    addVariable(n,"a",              NEITHER,        CELL, NEITHER,  min,    c,    SOUNDSPEED);
    addVariable(n,"uStar",          NEITHER,        CELL, NEITHER,  -c,     c,    USTAR);
    addVariable(n,"sigma",          NEITHER,        CELL, NEITHER,  -max,   max,  SIGMA,          1,3,3);


    addVariable(n,"rhomix",         NEITHER,        NOTCELL,NEITHER,   min, max,  RHO_MIX,        parameters.numberOfMixtures,2);
    addVariable(n,"lambda",         PRIMITIVE,      CELL,   REFINE,    0.0, 1.0,  LAMBDA,         parameters.numberOfMixtures);
    addVariable(n,"alpharholambda", CONSERVATIVE,   CELL,   NEITHER,   min, max,  ALPHARHOLAMBDA, parameters.numberOfMixtures);

    if(parameters.SOLID)
    {
        addVariable(n,"V",          BOTH,           CELL,    NEITHER,  -max, max,  V_TENSOR,       1,3,3);
        addVariable(n,"VStar",      NEITHER,        CELL,    NEITHER, -max, max,  VSTAR,          1,3,3);
        addVariable(n,"devH",       NEITHER,        NOTCELL, NEITHER, -max, max,  DEVH,           1,3,3);
        addVariable(n,"HenckyJ2",   NEITHER,        NOTCELL, NEITHER , -max, max,  HJ2);

        if(parameters.PLASTIC)
        {
            addVariable(n,"epsilon",            PRIMITIVE,    CELL, REFINE,   0.0, max,  EPSILON,         parameters.numberOfMaterials);
            addVariable(n,"alphaRhoEpsilon",    CONSERVATIVE, CELL, NEITHER,  0.0, max,  ALPHARHOEPSILON, parameters.numberOfMaterials);
        }
    }

}

int& AccessPattern::operator[](Variable var)
{
    return data[var];
}
