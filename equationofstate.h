#ifndef EQUATIONOFSTATE_H
#define EQUATIONOFSTATE_H

#include "amrexheader.h"
#include "cellarray.h"


class MieGruneisenEOS
{

public:

    ~MieGruneisenEOS(){}

    MieGruneisenEOS(){}

    MieGruneisenEOS(Real gamma, Real _pref, Real _eref, Real _CV);

    virtual void define(Vector<Real>& params);
    virtual Real coldCompressionInternalEnergy  (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real coldCompressionPressure        (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real shearInternalEnergy            (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real shearPressure                  (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real inverseGruneisen               (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual void rootFind                       (BoxAccessCellArray& U, int i, int j, int k, int m, Real kineticEnergy);
    virtual Real getSoundSpeedContribution      (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real getTemp                        (BoxAccessCellArray& U, int i, int j, int k, int m, int mixidx);
    virtual Real xi                             (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real componentShearModulus          (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real dGdrho                         (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real dG2drho2                       (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real transverseWaveSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual void setRhoFromDeformationTensor    (BoxAccessCellArray& U, int i, int j, int k, int m, double* F);
    virtual void defineMixtureDensities         (BoxAccessCellArray& U, int i, int j, int k, int m);

    void copy(MieGruneisenEOS& C);


    Real adiabaticIndex;
    Real GruneisenGamma;
    Real CV;
    Real pref;
    Real eref;

    bool mixture = false;
    int  mixtureIndex;

};

class SINGLEMATERIAL_MieGruneisenEOS : public MieGruneisenEOS
{

public:

    ~SINGLEMATERIAL_MieGruneisenEOS(){}

    SINGLEMATERIAL_MieGruneisenEOS(){}

    SINGLEMATERIAL_MieGruneisenEOS(Real gamma, Real _pref, Real _eref, Real _CV);

    virtual void define(Vector<Real>& params);
    virtual Real coldCompressionInternalEnergy  (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real coldCompressionPressure        (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real shearInternalEnergy            (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real shearPressure                  (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real inverseGruneisen               (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual void rootFind                       (BoxAccessCellArray& U, int i, int j, int k, int m, Real kineticEnergy);
    virtual Real getSoundSpeedContribution      (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real getTemp                        (BoxAccessCellArray& U, int i, int j, int k, int m, int mixidx);
    virtual Real xi                             (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real componentShearModulus          (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real dGdrho                         (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real dG2drho2                       (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real transverseWaveSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual void setRhoFromDeformationTensor    (BoxAccessCellArray& U, int i, int j, int k, int m, double* F);
    virtual void defineMixtureDensities         (BoxAccessCellArray& U, int i, int j, int k, int m);
};


class SINGLEMATERIAL_RomenskiiSolidEOS : public MieGruneisenEOS
{

public:

    ~SINGLEMATERIAL_RomenskiiSolidEOS(){}

    SINGLEMATERIAL_RomenskiiSolidEOS(){}

    virtual void define(Vector<Real>& params);
    virtual Real coldCompressionInternalEnergy  (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real coldCompressionPressure        (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real shearInternalEnergy            (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real shearPressure                  (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real getSoundSpeedContribution      (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real componentShearModulus          (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real dGdrho                         (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real dG2drho2                       (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual Real transverseWaveSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual void setRhoFromDeformationTensor    (BoxAccessCellArray& U, int i, int j, int k, int m, double* F);
    virtual Real inverseGruneisen               (BoxAccessCellArray& U, int i, int j, int k, int m);
    virtual void defineMixtureDensities         (BoxAccessCellArray& U, int i, int j, int k, int m);



    Real dpcdrho                                (BoxAccessCellArray& U, int i, int j, int k, int m);
    Real dpsdrho                                (BoxAccessCellArray& U, int i, int j, int k, int m);


    Real rho0;
    Real K0;
    Real EOSalpha;
    Real EOSbeta;
    Real G0;

};



#endif // EQUATIONOFSTATE_H
