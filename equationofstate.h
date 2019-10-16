#ifndef EQUATIONOFSTATE_H
#define EQUATIONOFSTATE_H

#include "amrexheader.h"
#include "cellarray.h"

/*
class EquationOfState
{
public:
    EquationOfState();
    virtual ~EquationOfState(){}

    Real adiabaticIndex;
    Real GruneisenGamma;
    Real CV;

    bool mixture = false;
    int  mixtureIndex;

    virtual void define(Vector<Real>& params) = 0;
    virtual Real coldCompressionInternalEnergy  (BoxAccessCellArray& U, int i, int j, int k, int m) = 0;
    virtual Real coldCompressionPressure        (BoxAccessCellArray& U, int i, int j, int k, int m) = 0;
    virtual Real inverseGruneisen               (BoxAccessCellArray& U, int i, int j, int k, int m) = 0;
    virtual void rootFind                       (BoxAccessCellArray& U, int i, int j, int k, int m) = 0;
    virtual Real getSoundSpeedContribution      (BoxAccessCellArray& U, int i, int j, int k, int m) = 0;

};
*/

class MieGruneisenEOS // : public EquationOfState
{

public:

    ~MieGruneisenEOS(){}

    MieGruneisenEOS(){}

    MieGruneisenEOS(Real gamma, Real _pref, Real _eref, Real _CV);

    virtual void define(Vector<Real>& params);
    virtual Real coldCompressionInternalEnergy  (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real coldCompressionPressure        (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real inverseGruneisen               (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual void rootFind                       (BoxAccessCellArray& U, int i, int j, int k, int m, Real kineticEnergy);
    virtual Real getSoundSpeedContribution      (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real getTemp                        (BoxAccessCellArray& U, int i, int j, int k, int m, int mixidx);
    virtual Real xi                             (BoxAccessCellArray& U, int i, int j, int k, int m);

    void copy(MieGruneisenEOS& C);


    Real adiabaticIndex;
    Real GruneisenGamma;
    Real CV;
    Real pref;
    Real eref;

    bool mixture = false;
    int  mixtureIndex;

};


class MixtureEOS : public MieGruneisenEOS
{

public:

    ~MixtureEOS(){}

    MixtureEOS();

    virtual void define(Vector<Real>& params);


    virtual Real coldCompressionInternalEnergy  (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real coldCompressionPressure        (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real inverseGruneisen               (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual void rootFind                       (BoxAccessCellArray& U, int i, int j, int k, int m, Real kineticEnergy);
    virtual Real getSoundSpeedContribution      (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real getTemp                        (BoxAccessCellArray& U, int i, int j, int k, int m, int mixidx);
    virtual Real xi                             (BoxAccessCellArray& U, int i, int j, int k, int m);

    Real bisectionFunction(BoxAccessCellArray& U, int i, int j, int k, int m, Real rhoTry, Real kineticEnergy, Real& p);
    Real pressureFunction (BoxAccessCellArray& U, int i, int j, int k, int m, Real& p);
    Real mixtureSoundSpeed(BoxAccessCellArray& U, int i, int j, int k, int m);


    MieGruneisenEOS first;
    MieGruneisenEOS second;

};

#endif // EQUATIONOFSTATE_H
