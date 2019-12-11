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



void SINGLEMATERIAL_RomenskiiSolidEOS::define(Vector<Real> &params)
{
    adiabaticIndex 	= params[0];
    GruneisenGamma 	= adiabaticIndex-1.0;
    pref 			= params[1];
    eref	 		= params[2];
    CV 				= params[3];

    rho0            = params[4];
    K0              = params[5];
    EOSalpha        = params[6];
    EOSbeta         = params[7];
    G0              = params[8];
}

Real SINGLEMATERIAL_RomenskiiSolidEOS::coldCompressionInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    Real x = U(i,j,k,RHO,m)/rho0;

    return U(i,j,k,RHO,m)*((K0/(2.0*rho0*EOSalpha*EOSalpha))*(std::pow(x,EOSalpha)-1.0)*(std::pow(x,EOSalpha)-1.0));
}

Real SINGLEMATERIAL_RomenskiiSolidEOS::coldCompressionPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    Real x = U(i,j,k,RHO,m)/rho0;

    return ((K0/EOSalpha)*std::pow(x,EOSalpha+1.0)*(std::pow(x,EOSalpha)-1.0))/GruneisenGamma;
}

Real SINGLEMATERIAL_RomenskiiSolidEOS::getSoundSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return (U(i,j,k,P)*(GruneisenGamma+1.0)-(coldCompressionPressure(U,i,j,k,m)+shearPressure(U,i,j,k,m))*(GruneisenGamma+1.0))/U(i,j,k,RHO,m) + dpcdrho(U,i,j,k,m) + dpsdrho(U,i,j,k,m) + (4.0/3.0)*componentShearModulus(U,i,j,k,m)/U(i,j,k,RHO,m);
}

Real SINGLEMATERIAL_RomenskiiSolidEOS::componentShearModulus(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return G0*std::pow((U(i,j,k,RHO,m)/rho0),EOSbeta+1.0);
}

Real SINGLEMATERIAL_RomenskiiSolidEOS::dGdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return G0*(EOSbeta+1.0)*std::pow((U(i,j,k,RHO,m)/rho0),EOSbeta)/rho0;
}

Real SINGLEMATERIAL_RomenskiiSolidEOS::dG2drho2(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return G0*(EOSbeta+1.0)*(EOSbeta)*std::pow((U(i,j,k,RHO,m)/rho0),EOSbeta-1.0)/(rho0*rho0);
}

Real SINGLEMATERIAL_RomenskiiSolidEOS::transverseWaveSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
     return componentShearModulus(U,i,j,k,m)/(U(i,j,k,RHO,m));
}

Real SINGLEMATERIAL_RomenskiiSolidEOS::shearInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return (componentShearModulus(U,i,j,k,m)*U(i,j,k,HJ2,m));
}

Real SINGLEMATERIAL_RomenskiiSolidEOS::shearPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return ((U(i,j,k,RHO,m)*dGdrho(U,i,j,k,m)-componentShearModulus(U,i,j,k,m))*U(i,j,k,HJ2,m))/GruneisenGamma;
}

Real SINGLEMATERIAL_RomenskiiSolidEOS::dpcdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    Real x = U(i,j,k,RHO,m)/rho0;

    return (K0/(EOSalpha*rho0))*std::pow(x,EOSalpha)*( (2.0*EOSalpha+1.0)*std::pow(x,EOSalpha) -  (EOSalpha+1.0) );
}

Real SINGLEMATERIAL_RomenskiiSolidEOS::dpsdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,RHO,m)*dG2drho2(U,i,j,k,m)*U(i,j,k,HJ2,m);
}

void SINGLEMATERIAL_RomenskiiSolidEOS::setRhoFromDeformationTensor(BoxAccessCellArray& U, int i, int j, int k, int m, double* F)
{
    double determinant = det(F);

    U(i,j,k,RHO,m) = rho0/determinant;

}

Real SINGLEMATERIAL_RomenskiiSolidEOS::inverseGruneisen(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 1.0/GruneisenGamma;
}

void SINGLEMATERIAL_RomenskiiSolidEOS::defineMixtureDensities(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return;
}
