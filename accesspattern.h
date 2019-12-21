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
    void addVariable(int& position, std::string nameBase, Var_type type, Var_type INREFINE, Real lo, Real hi, Variable var, int m, int componentNumber, int rowNumber, int colNumber);
    int& operator[](MaterialSpecifier& m);

    std::vector<std::vector<int> > data;
    std::map<Variable,std::pair<Real,Real> > limits;

    int numberOfSharpMaterials;

    Vector<int>                 numberOfDiffuseMaterials;

    Vector<Vector<MaterialSpecifier> >   conservativeVariables;
    Vector<Vector<MaterialSpecifier> >   primitiveVariables;
    Vector<MaterialSpecifier>            allVariables;
    Vector<MaterialSpecifier>            refineVariables;

    Vector<Vector<int> >                 numberOfComponentsForVariable;
    Vector<Vector<int> >                 numberOfRowsForVariable;

    Vector<std::string>         variableNames;
    Vector< Vector<MaterialDescriptor> >&  materialInfo;

    Vector<Interface_type>& interface;


    Vector< Vector<MaterialSpecifier> > material_conservativeVariables;
    Vector< Vector<MaterialSpecifier> > material_primitiveVariables;



};

#endif // ACCESSPATTERN_H
