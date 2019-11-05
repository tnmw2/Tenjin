#ifndef ACCESSPATTERN_H
#define ACCESSPATTERN_H

#include "structdefinitions.h"

/** \class AccessPattern
 * Because the number of materials in a simulation can change,
 * the location of each thermodynamic variable in a Multifab
 * will also change. This class keeps track of this my providing
 * a map from the variable name to the an int giving start of the
 * block where that variable is stored.
 *
 * Variables are stored contiguously. For example rho_0 is followed
 * by rho_1 and u_x is followed by u_y, hence we only need the int
 * giving us the start of the contiguous block.
 *
 * We also store a list of the conservative and primitve variabels
 * for easy looping.
 */

/** A struct that contains the particular information of a material.
 */

class AccessPattern
{
public:

    AccessPattern(ParameterStruct &parameters);

    void define(ParameterStruct &parameters);
    void addVariable(int& position, std::string nameBase, Var_type type, Var_type INCELL, Variable var, int materialNumber, int rowNumber, int colNumber);

    int& operator[](Variable var);

    std::map<Variable,int> data;

    Vector<MaterialSpecifier>   conservativeVariables;
    Vector<MaterialSpecifier>   primitiveVariables;
    Vector<MaterialSpecifier>   allVariables;


    Vector<std::string>         variableNames;
    Vector<MaterialDescriptor>&  materialInfo;

    Vector<MaterialSpecifier>   cellVariables;

};

#endif // ACCESSPATTERN_H
