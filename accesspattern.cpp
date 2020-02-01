#include "simulationheader.h"

void AccessPattern::addVariable(int& position, std::string nameBase, Var_type type, Var_type INCELL, Var_type INREFINE, Real refineVal, Real lo, Real hi, Variable var, int materialNumber=1, int rowNumber=1, int colNumber=1)
{
    if(materialNumber*rowNumber*colNumber > 0)
    {
        //data.insert(std::pair<Variable,int>(var,position));

        data[var] = position;

        numberOfMaterialsForVariable[var] = ( (materialNumber > 1 && rowNumber*colNumber > 1) ? rowNumber*colNumber : 1  );
        numberOfRowsForVariable[var]      = ( (rowNumber > 1 && colNumber > 1) ? colNumber : 1);

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
                        refineValue.push_back(refineVal);
                    }

                    allVariables.push_back(MaterialSpecifier(var,m,row,col));


                }
            }
        }
    }
}

AccessPattern::AccessPattern(ParameterStruct& _parameters) : materialInfo(_parameters.materialInfo), parameters(_parameters)
{
    define(_parameters);
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
    refineValue.clear();
    variableNames.clear();
    cellVariables.clear();
    numberOfMaterialsForVariable.clear();
    numberOfRowsForVariable.clear();

    Real c   = 3E8;
    Real min = 1E-20;
    Real max = 1E20;

    data.resize(50);

    numberOfMaterialsForVariable.resize(50);
    numberOfRowsForVariable.resize(50);


    std::fill(data.begin(),data.end(),0);
    std::fill(numberOfMaterialsForVariable.begin(),numberOfMaterialsForVariable.end(),0);
    std::fill(numberOfRowsForVariable.begin(),numberOfRowsForVariable.end(),0);



    addVariable(n,"alpha"   ,       BOTH,           CELL, REFINE,  0.1,  min,    1.0,  ALPHA,          parameters.numberOfMaterials);
    addVariable(n,"alphaRho",       CONSERVATIVE,   CELL, REFINE, 0.0,  min,    max,  ALPHARHO,       parameters.numberOfMaterials);
    addVariable(n,"rho_k",          PRIMITIVE,      CELL, NEITHER, 0.0,  min,    max,  RHO_K,          parameters.numberOfMaterials);
    addVariable(n,"rho",            NEITHER,        CELL, NEITHER,  200.0,  min,    max,  RHO);
    addVariable(n,"rhoU",           CONSERVATIVE,   CELL, NEITHER, 0.0, -max,    max,  RHOU,           1,3);
    addVariable(n,"u",              PRIMITIVE,      CELL, NEITHER, 0.0,-c ,     c,    VELOCITY,       1,3);
    addVariable(n,"E",              CONSERVATIVE,   CELL, NEITHER, 0.0,-max,    max,  TOTAL_E);
    addVariable(n,"p",              PRIMITIVE,      CELL, REFINE,  1E10, -max,    max,  P);
    addVariable(n,"a",              NEITHER,        CELL, NEITHER, 0.0,  min,    c,    SOUNDSPEED);
    addVariable(n,"uStar",          NEITHER,        CELL, NEITHER, 0.0,-c,      c,    USTAR);
    addVariable(n,"sigma",          NEITHER,        CELL, NEITHER, 0.0,-max,    max,  SIGMA,          1,3,3);

    if(parameters.numberOfMixtures > 0)
    {
        addVariable(n,"rhomix",         NEITHER,        NOTCELL,NEITHER, 0.0,  min, max,  RHO_MIX,        parameters.numberOfMaterials,2);
        addVariable(n,"lambda",         PRIMITIVE,      CELL,   NEITHER, 0.0,  0.0, 1.0,  LAMBDA,         parameters.numberOfMaterials);
        addVariable(n,"alpharholambda", CONSERVATIVE,   CELL,   NEITHER, 0.0,  min, max,  ALPHARHOLAMBDA, parameters.numberOfMaterials);
    }

    if(parameters.SOLID)
    {
        addVariable(n,"V",          BOTH,           CELL,    NEITHER, 0.0, -max, max,  V_TENSOR,       1,3,3);
        addVariable(n,"VStar",      NEITHER,        CELL,    NEITHER, 0.0, -max, max,  VSTAR,          1,3,3);
        addVariable(n,"devH",       NEITHER,        NOTCELL, NEITHER, 0.0, -max, max,  DEVH,           1,3,3);
        addVariable(n,"HenckyJ2",   NEITHER,        NOTCELL, NEITHER, 0.0, -max, max,  HJ2);

        if(parameters.PLASTIC)
        {
            addVariable(n,"epsilon",            PRIMITIVE,    CELL, NEITHER, 0.0, 0.0, max,  EPSILON,         parameters.numberOfMaterials);
            addVariable(n,"alphaRhoEpsilon",    CONSERVATIVE, CELL, NEITHER, 0.0, 0.0, max,  ALPHARHOEPSILON, parameters.numberOfMaterials);
        }
    }
}

int& AccessPattern::operator[](Variable var)
{
    return data[var];
}


Real AccessPattern::refineCriterion(int m)
{
    return refineValue[m];
}
