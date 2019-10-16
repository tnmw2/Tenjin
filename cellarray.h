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

    CellArray(BoxArray& ba, DistributionMapping& dm, const int Ncomp, const int Nghost, AccessPattern& _accessPattern, ParameterStruct &parameters);

    void primitiveToConservative();
    void conservativeToPrimitive();
    void getSoundSpeed          ();


    void        operator=(CellArray& U);
    CellArray&  operator*(Real d);
    CellArray&  operator+(CellArray& U);

    MultiFab data;

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

    Real& left (Direction_enum d, int i, int j, int k, MaterialSpecifier& m);
    Real& right(Direction_enum d, int i, int j, int k, MaterialSpecifier& m);

    Real& left (Direction_enum d, int i, int j, int k, Variable var, int mat=0, int row=0, int col=0);
    Real& right(Direction_enum d, int i, int j, int k, Variable var, int mat=0, int row=0, int col=0);

    Real& operator()(int i, int j, int k, MaterialSpecifier& m);
    Real& operator()(int i, int j, int k, Variable var, int mat=0, int row=0, int col=0);
    Real& operator()(int i, int j, int k, int var, int mat=0, int row=0, int col=0);


    void  conservativeToPrimitive();
    void  primitiveToConservative();

    Real getEffectiveInverseGruneisen           (int i, int j, int k);
    Real getEffectiveNonThermalPressure         (int i, int j, int k);
    Real getEffectiveNonThermalInternalEnergy   (int i, int j, int k);

    void getSoundSpeed();


};



#endif // CELLARRAY_H
