#include "simulationheader.h"

MieGruneisenEOS::MieGruneisenEOS(Real gamma, Real _pref=0.0, Real _eref=0.0, Real _CV=0.0)
{
    adiabaticIndex 	= gamma;
    GruneisenGamma 	= gamma-1.0;
    pref 			= _pref;
    eref	 		= _eref;
    CV 				= _CV;
}

void MieGruneisenEOS::define(Vector<Real> &params)
{


    adiabaticIndex 	= params[0];
    GruneisenGamma 	= adiabaticIndex-1.0;
    pref 			= params[1];
    eref	 		= params[2];
    CV 				= params[3];



}

void MieGruneisenEOS::copy(MieGruneisenEOS& C)
{
    adiabaticIndex = C.adiabaticIndex;
    GruneisenGamma = C.GruneisenGamma;
    pref = C.pref;
    eref = C.eref;
    CV = C.CV;
}

Real MieGruneisenEOS::coldCompressionInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return eref*U(i,j,k,ALPHARHO,m);
}

Real MieGruneisenEOS::coldCompressionPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return pref*U(i,j,k,ALPHA,m)/GruneisenGamma;
}

Real MieGruneisenEOS::inverseGruneisen(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,ALPHA,m)/GruneisenGamma;
}

void MieGruneisenEOS::rootFind(BoxAccessCellArray& U, int i, int j, int k, int m, Real kineticEnergy)
{
    U(i,j,k,RHO_MIX,m,0) *= 1.0;

    Print() << "Trying to root find on a non-mixture EOS" << std::endl;

    exit(1);
}

Real MieGruneisenEOS::getSoundSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,P)*(GruneisenGamma+1.0)/U(i,j,k,RHO_K,m) - pref/U(i,j,k,RHO_K,m);
}

Real MieGruneisenEOS::getTemp(BoxAccessCellArray& U, int i, int j, int k, int m, int mixidx)
{
    mixidx *= 1;
    return U(i,j,k,P)/(U(i,j,k,RHO_K,m)*GruneisenGamma*CV);
}

Real MieGruneisenEOS::xi(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 1.0/GruneisenGamma;
}

Real MieGruneisenEOS::componentShearModulus(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0*U(i,j,k,0,m);
}

Real MieGruneisenEOS::dGdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0*U(i,j,k,0,m);
}

Real MieGruneisenEOS::dG2drho2(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0*U(i,j,k,0,m);
}

Real MieGruneisenEOS::transverseWaveSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0*U(i,j,k,0,m);
}

Real MieGruneisenEOS::shearInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0*U(i,j,k,0,m);
}

Real MieGruneisenEOS::shearPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0*U(i,j,k,0,m);
}

void MieGruneisenEOS::setRhoFromDeformationTensor(BoxAccessCellArray& U, int i, int j, int k, int m, double* F)
{
    return;
}

void MieGruneisenEOS::defineMixtureDensities(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return;
}


/*****************************************************
 * Single Material EOS
 ****************************************************/

SINGLEMATERIAL_MieGruneisenEOS::SINGLEMATERIAL_MieGruneisenEOS(Real gamma, Real _pref=0.0, Real _eref=0.0, Real _CV=0.0)
{
    adiabaticIndex 	= gamma;
    GruneisenGamma 	= gamma-1.0;
    pref 			= _pref;
    eref	 		= _eref;
    CV 				= _CV;
}

void SINGLEMATERIAL_MieGruneisenEOS::define(Vector<Real> &params)
{


    adiabaticIndex 	= params[0];
    GruneisenGamma 	= adiabaticIndex-1.0;
    pref 			= params[1];
    eref	 		= params[2];
    CV 				= params[3];



}

Real SINGLEMATERIAL_MieGruneisenEOS::coldCompressionInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return eref*U(i,j,k,RHO,m);
}

Real SINGLEMATERIAL_MieGruneisenEOS::coldCompressionPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return pref/GruneisenGamma;
}

Real SINGLEMATERIAL_MieGruneisenEOS::inverseGruneisen(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 1.0/GruneisenGamma;
}

void SINGLEMATERIAL_MieGruneisenEOS::rootFind(BoxAccessCellArray& U, int i, int j, int k, int m, Real kineticEnergy)
{
    amrex::Abort("Trying to root find on a non-mixture EOS");
}

Real SINGLEMATERIAL_MieGruneisenEOS::getSoundSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,P,m)*(GruneisenGamma+1.0)/U(i,j,k,RHO,m) - pref/U(i,j,k,RHO,m);
}

Real SINGLEMATERIAL_MieGruneisenEOS::getTemp(BoxAccessCellArray& U, int i, int j, int k, int m, int mixidx)
{
    return U(i,j,k,P,m)/(U(i,j,k,RHO,m)*GruneisenGamma*CV);
}

Real SINGLEMATERIAL_MieGruneisenEOS::xi(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 1.0/GruneisenGamma;
}

Real SINGLEMATERIAL_MieGruneisenEOS::componentShearModulus(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real SINGLEMATERIAL_MieGruneisenEOS::dGdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real SINGLEMATERIAL_MieGruneisenEOS::dG2drho2(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real SINGLEMATERIAL_MieGruneisenEOS::transverseWaveSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real SINGLEMATERIAL_MieGruneisenEOS::shearInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real SINGLEMATERIAL_MieGruneisenEOS::shearPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

void SINGLEMATERIAL_MieGruneisenEOS::setRhoFromDeformationTensor(BoxAccessCellArray& U, int i, int j, int k, int m, double* F)
{
    return;
}

void SINGLEMATERIAL_MieGruneisenEOS::defineMixtureDensities(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return;
}


