#ifndef CELLARRAY_H
#define CELLARRAY_H

#include "amrexheader.h"
#include "structdefinitions.h"
#include "accesspattern.h"


/** \class CellArray
 * A wrapper for a Multifab which allows for easy addition and multiplication.
 * Also stores AccessPatterns and numbers of materials.
 *
 * \class BoxAccessCellArray
 * A wrapper for a Box and an FArrayBox that can return data using an AccessPattern.conservativeVariables
 * Contains functions to do things like convert primitive variables to conservative variables etc.
 */

class CellArray
{
public:

    CellArray(MultiFab& S, AccessPattern &_accessPattern, ParameterStruct& parameters);

    void primitiveToConservative();
    void conservativeToPrimitive();
    void getSoundSpeed          ();
    void cleanUpV               ();
    void cleanUpAlpha           ();

    void        operator=(CellArray& U);
    CellArray&  operator*(Real d);
    CellArray&  operator+(CellArray& U);
    bool        contains_nan();


    MultiFab& data;

    AccessPattern& accessPattern;

    int numberOfMaterials;



};

class BoxAccessCellArray
{
public:

    BoxAccessCellArray(const Box& bx, FArrayBox& fb, CellArray &U);

    BoxAccessCellArray(MFIter& mfi, const Box &bx, CellArray &U);

    const Box& box;

    FArrayBox& fab;

    AccessPattern& accessPattern;

    int numberOfMaterials;

    static const int numberOfComponents = 3;

    Real& left (Direction_enum d, int i, int j, int k, MaterialSpecifier& m);
    Real& right(Direction_enum d, int i, int j, int k, MaterialSpecifier& m);

    Real& left (Direction_enum d, int i, int j, int k, Variable var, int mat=0, int row=0, int col=0);
    Real& right(Direction_enum d, int i, int j, int k, Variable var, int mat=0, int row=0, int col=0);

    Real& neighbour(int di, int dj, int dk, int i, int j, int k, MaterialSpecifier& m);
    Real& neighbour(int di, int dj, int dk, int i, int j, int k, Variable var, int mat=0, int row=0, int col=0);

    Real& operator()(int i, int j, int k, MaterialSpecifier& m);
    Real& operator()(int i, int j, int k, Variable var, int mat=0, int row=0, int col=0);
    Real& operator()(int i, int j, int k, int var, int mat=0, int row=0, int col=0);


    void  conservativeToPrimitive();
    void  primitiveToConservative();

    void stressTensor                           (int i, int j, int k);

    Real getEffectiveInverseGruneisen           (int i, int j, int k);
    Real getEffectiveNonThermalPressure         (int i, int j, int k);
    Real getEffectiveNonThermalInternalEnergy   (int i, int j, int k);

    void getSoundSpeed();
    Real transverseWaveSpeed(int i, int j, int k);

    void amrexToArray(int i, int j, int k, Variable var, int m, double* copy, int nx=3, int ny=3);

    void getHenckyJ2(int i, int j, int k);
    void getDeviatoricHenckyStrain(int i, int j, int k);
    void normaliseV();
    void cleanUpV();
    void cleanUpAlpha();




    bool check(MaterialSpecifier& m);
    bool contains_nan();

    bool cellIsMostlyFluid(int i, int j, int k);


};



#endif // CELLARRAY_H
