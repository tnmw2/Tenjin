#include "simulationheader.h"

void RomenskiiSolidEOS::define(Vector<Real> &params)
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

Real RomenskiiSolidEOS::coldCompressionInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    Real x = U(i,j,k,RHO_K,m)/rho0;

    return U(i,j,k,ALPHARHO,m)*((K0/(2.0*rho0*EOSalpha*EOSalpha))*(std::pow(x,EOSalpha)-1.0)*(std::pow(x,EOSalpha)-1.0));
}

Real RomenskiiSolidEOS::coldCompressionPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    Real x = U(i,j,k,RHO_K,m)/rho0;

    return (U(i,j,k,ALPHA,m)/GruneisenGamma)*( (K0/EOSalpha)*std::pow(x,EOSalpha+1.0)*(std::pow(x,EOSalpha)-1.0));
}

Real RomenskiiSolidEOS::getSoundSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return (U(i,j,k,P)*(GruneisenGamma+1.0)-(coldCompressionPressure(U,i,j,k,m)+shearPressure(U,i,j,k,m)))/U(i,j,k,RHO_K,m) + dpcdrho(U,i,j,k,m) + dpsdrho(U,i,j,k,m) + (4.0/3.0)*componentShearModulus(U,i,j,k,m)/U(i,j,k,RHO_K,m);
}

Real RomenskiiSolidEOS::componentShearModulus(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return G0*std::pow((U(i,j,k,RHO_K,m)/rho0),EOSbeta+1.0);
}

Real RomenskiiSolidEOS::dGdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return G0*(EOSbeta+1.0)*std::pow((U(i,j,k,RHO_K,m)/rho0),EOSbeta)/rho0;
}

Real RomenskiiSolidEOS::dG2drho2(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return G0*(EOSbeta+1.0)*(EOSbeta)*std::pow((U(i,j,k,RHO_K,m)/rho0),EOSbeta-1.0)/(rho0*rho0);
}

Real RomenskiiSolidEOS::transverseWaveSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
     return U(i,j,k,ALPHA,m)*componentShearModulus(U,i,j,k,m)/(U(i,j,k,RHO_K,m)*GruneisenGamma);
}

Real RomenskiiSolidEOS::shearInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,ALPHA,m)*(componentShearModulus(U,i,j,k,m)*U(i,j,k,HJ2));
}

Real RomenskiiSolidEOS::shearPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return (U(i,j,k,ALPHA,m)/GruneisenGamma)*((U(i,j,k,RHO_K,m)*dGdrho(U,i,j,k,m)-componentShearModulus(U,i,j,k,m))*U(i,j,k,HJ2));
}

Real RomenskiiSolidEOS::dpcdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    Real x = U(i,j,k,RHO_K,m)/rho0;

    return (K0/(EOSalpha*rho0))*std::pow(x,EOSalpha)*( (2.0*EOSalpha+1.0)*std::pow(x,EOSalpha) -  (EOSalpha+1.0) );
}

Real RomenskiiSolidEOS::dpsdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,RHO_K,m)*dG2drho2(U,i,j,k,m)*U(i,j,k,HJ2);
}

void RomenskiiSolidEOS::setRhoFromDeformationTensor(BoxAccessCellArray& U, int i, int j, int k, int m, double* F)
{
    double determinant = det(F);

    U(i,j,k,RHO_K,m) = rho0/determinant;

}

Real RomenskiiSolidEOS::inverseGruneisen(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,ALPHA,m)/GruneisenGamma;
}

void RomenskiiSolidEOS::rootFind(BoxAccessCellArray& U, int i, int j, int k, int m, Real kineticEnergy)
{
    Print() << "Trying to root find on a non-mixture EOS" << std::endl;

    exit(1);
}

Real RomenskiiSolidEOS::getTemp(BoxAccessCellArray& U, int i, int j, int k, int m, int mixidx)
{
    return U(i,j,k,P)/(U(i,j,k,RHO_K,m)*GruneisenGamma*CV);
}

Real RomenskiiSolidEOS::xi(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 1.0/GruneisenGamma;
}

