#include "simulationheader.h"

//EquationOfState::EquationOfState(){}


/*****************************************************
 * Mixture EOS
 ****************************************************/

MixtureEOS::MixtureEOS()
{
    first.mixture  = true;
    second.mixture = true;

    //first.mixtureIndex = 0;
    //second.mixtureIndex = 1;
}

Real MixtureEOS::coldCompressionInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,LAMBDA,m)*first.coldCompressionInternalEnergy(U,i,j,k,m) + (1.0-U(i,j,k,LAMBDA,m))*second.coldCompressionInternalEnergy(U,i,j,k,m);
    //return U(i,j,k,ALPHARHOLAMBDA,m)*first.coldCompressionInternalEnergy(U,i,j,k,m) + (U(i,j,k,ALPHARHO,m)-U(i,j,k,ALPHARHOLAMBDA,m))*second.coldCompressionInternalEnergy(U,i,j,k,m);
}

Real MixtureEOS::coldCompressionPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return (first.coldCompressionPressure(U,i,j,k,m))*U(i,j,k,LAMBDA,m)*U(i,j,k,RHO_K,m)/(U(i,j,k,RHO_MIX,m,0))+(second.coldCompressionPressure(U,i,j,k,m))*(1.0-U(i,j,k,LAMBDA,m))*U(i,j,k,RHO_K,m)/(U(i,j,k,RHO_MIX,m,1));
    //return (first.coldCompressionPressure(U,i,j,k,m))*U(i,j,k,ALPHARHOLAMBDA,m)/(U(i,j,k,RHO_MIX,m,0)*first.GruneisenGamma)+(second.coldCompressionPressure(U,i,j,k,m))*(U(i,j,k,ALPHARHO,m)-U(i,j,k,ALPHARHOLAMBDA,m))/(U(i,j,k,RHO_MIX,m,1)*second.GruneisenGamma);
}

Real MixtureEOS::inverseGruneisen(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,ALPHARHOLAMBDA,m)/(first.GruneisenGamma*U(i,j,k,RHO_MIX,m,0))+(U(i,j,k,ALPHARHO,m)-U(i,j,k,ALPHARHOLAMBDA,m))/(second.GruneisenGamma*U(i,j,k,RHO_MIX,m,1));
}

void MixtureEOS::define(Vector<Real> &params)
{
    Vector<Real> _first;
    Vector<Real> _second;

    for(unsigned int i=0;i<params.size();i++)
    {
        if(i<4)
        {
            _first.push_back(params[i]);
        }
        else
        {
            _second.push_back(params[i]);
        }
    }

    first.define(_first);
    second.define(_second);

    return;
}

Real MixtureEOS::getSoundSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return mixtureSoundSpeed(U,i,j,k,m);
}

Real MixtureEOS::mixtureSoundSpeed(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    /******************************
     * Assumed constant pref terms
     ********************************/

    double dedrho = (U(i,j,k,LAMBDA,m))*((first.pref -U(i,j,k,P))/(first.GruneisenGamma*U(i,j,k,RHO_K,m)*U(i,j,k,RHO_K,m)));

    dedrho +=   (1.0-U(i,j,k,LAMBDA,m))*((second.pref -U(i,j,k,P))/(second.GruneisenGamma*U(i,j,k,RHO_K,m)*U(i,j,k,RHO_K,m)));

    return (  U(i,j,k,P)/(U(i,j,k,RHO_K,m)*U(i,j,k,RHO_K,m))  - dedrho )/(U(i,j,k,LAMBDA,m)/(U(i,j,k,RHO_MIX,m,0)*first.GruneisenGamma)+(1.0-U(i,j,k,LAMBDA,m))/(U(i,j,k,RHO_MIX,m,1)*second.GruneisenGamma));

}

Real MixtureEOS::xi(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 1.0/((U(i,j,k,LAMBDA,m)*first.adiabaticIndex*first.CV+(1.0-U(i,j,k,LAMBDA,m))*second.adiabaticIndex*second.CV)/(U(i,j,k,LAMBDA,m)*first.CV+(1.0-U(i,j,k,LAMBDA,m))*second.CV)-1.0);
}


/******************************************************
 * Root Finding stuff
 *****************************************************/

/**Sign function.
 */
template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}


Real rhoaFunc(BoxAccessCellArray& U, int i, int j, int k, int m, Real rhob)
{
    return U(i,j,k,LAMBDA,m)/(1.0/U(i,j,k,RHO_K,m)-(1.0-U(i,j,k,LAMBDA,m))/rhob);
}

Real rhobFunc(BoxAccessCellArray& U, int i, int j, int k, int m, Real rhoa)
{
    return (1.0-U(i,j,k,LAMBDA,m))/(1.0/U(i,j,k,RHO_K,m)-U(i,j,k,LAMBDA,m)/rhoa);
}

Real MixtureEOS::pressureFunction(BoxAccessCellArray& U, int i, int j, int k, int m, Real& p)
{

    //double prefa = getColdCompressionPressure(parameters,m,1,0,m);
    //double prefb = getColdCompressionPressure(parameters,m,1,1,mixtureIndex())+parameters.pref[mixtureIndex()];

    return ((p-first.pref)/(U(i,j,k,RHO_MIX,m,0)*first.GruneisenGamma*first.CV))-((p-second.pref)/(U(i,j,k,RHO_MIX,m,1)*second.GruneisenGamma*second.CV));
    //return ((p-prefa)/(rho_mix[0]*parameters.GruneisenGamma[0]*parameters.CV[m]))-((p-prefb)/(rho_mix[1]*parameters.GruneisenGamma[mixtureIndex()]*parameters.CV[mixtureIndex()]));

}

Real MixtureEOS::bisectionFunction(BoxAccessCellArray& U, int i, int j, int k, int m, Real rhoTry, Real kineticEnergy, Real& p)
{
    U(i,j,k,RHO_MIX,m,0) = rhoTry;
    U(i,j,k,RHO_MIX,m,1) = rhobFunc(U,i,j,k,m,rhoTry);

    p = (U(i,j,k,TOTAL_E) - kineticEnergy - U.getEffectiveNonThermalInternalEnergy(i,j,k)+ U.getEffectiveNonThermalPressure(i,j,k))/(U.getEffectiveInverseGruneisen(i,j,k));
    //p = (U(i,j,k,TOTAL_E)-kineticEnergy)/(U.getEffectiveInverseGruneisen(i,j,k));
    //p = (E- kineticEnergy(conservative) - getEffectiveNonThermalInternalEnergy(parameters) + getEffectiveNonThermalPressure(parameters))/getEffectiveInverseGruneisen(parameters) ;

    return pressureFunction(U,i,j,k,m,p);
}

void MixtureEOS::rootFind(BoxAccessCellArray& U, int i, int j, int k, int m, Real kineticEnergy)
{
    static const Real toleranceMixtureNotPresent = 1E-1;
    static const Real toleranceForSinglePhaseTreatment = 1E-2;
    static const Real toleranceForConvergence = 1E-3;
    static const Real toleranceForBeingNearRoot = 1E2;

    /****************************************************
     * NB: Changed these if clauses as they were giving NaNs
     * **************************************************/

    if(first.pref == 0.0 && second.pref == 0.0)
    {
        U(i,j,k,RHO_MIX,m,0) = U(i,j,k,RHO_K,m)*(U(i,j,k,LAMBDA,m)+(1.0-U(i,j,k,LAMBDA,m))*(second.GruneisenGamma*second.CV)/(first.GruneisenGamma*first.CV));
        U(i,j,k,RHO_MIX,m,1) = rhobFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,0));

        return;
    }

    if(std::isnan(U(i,j,k,LAMBDA,m)))
    {
        Print() << "Nan in Lambda in root finding" << std::endl;
    }

    if(U(i,j,k,ALPHA,m)<toleranceMixtureNotPresent)
    {
        U(i,j,k,RHO_MIX,m,0) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,1) = U(i,j,k,RHO_K,m);

        return;
    }

    if(U(i,j,k,LAMBDA,m)> 1.0-toleranceForSinglePhaseTreatment)
    {
        U(i,j,k,RHO_MIX,m,0) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,1) = U(i,j,k,RHO_K,m);

        return;
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        U(i,j,k,RHO_MIX,m,1) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,0) = U(i,j,k,RHO_K,m);

        return;
    }


    Real A = U(i,j,k,LAMBDA,m)*U(i,j,k,RHO_K,m)+1E-10;
    Real B = 200000.0; //100 for banks nondim

    Real p;

    Real mid = 0.5*(A+B);

    if(std::abs(bisectionFunction(U,i,j,k,m,A,kineticEnergy,p))<= toleranceForBeingNearRoot)
    {
        U(i,j,k,RHO_MIX,m,0) = A;
        U(i,j,k,RHO_MIX,m,1) = rhobFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,0));

        return;
    }
    else if(std::abs(bisectionFunction(U,i,j,k,m,B,kineticEnergy,p))<= toleranceForBeingNearRoot)
    {
        U(i,j,k,RHO_MIX,m,0) = B;
        U(i,j,k,RHO_MIX,m,1) = rhobFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,0));

        return;
    }


    if(sgn(bisectionFunction(U,i,j,k,m,A,kineticEnergy,p)) == sgn(bisectionFunction(U,i,j,k,m,B,kineticEnergy,p)) )
    {
        Print() << "Error in root Bisection " << std::endl;

        Print() << " A: " << A << " B: " << B << " at " << i << " " <<j << " "<< k << std::endl;
        Print() << bisectionFunction(U,i,j,k,m,A,kineticEnergy,p) << " " << bisectionFunction(U,i,j,k,m,B,kineticEnergy,p)<< std::endl;
        Print() << " alpha: " << U(i,j,k,ALPHA,m) << " lambda: " << U(i,j,k,LAMBDA,m) << std::endl;
        exit(1);
    }



    for(int n=0;;n++)
    {
        //Print() << A << " " << B << std::endl;

        if(sgn(bisectionFunction(U,i,j,k,m,mid,kineticEnergy,p)) == sgn(bisectionFunction(U,i,j,k,m,A,kineticEnergy,p)))
        {
            A = mid;
        }
        else
        {
            B = mid;
        }

        mid = 0.5*(A+B);


        if(std::abs(A-B)/(A+B)<toleranceForConvergence)
        {
            U(i,j,k,RHO_MIX,m,0) = mid;
            U(i,j,k,RHO_MIX,m,1) = rhobFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,0));

            return;
        }

        if(n>1000)
        {
            Print() << "Root finding is taking a while" << std::endl;
        }
    }
}

Real MixtureEOS::getTemp(BoxAccessCellArray& U, int i, int j, int k, int m, int mixidx)
{
    if(mixidx == 0)
    {
        return U(i,j,k,P)/(U(i,j,k,RHO_MIX,m,mixidx)*first.GruneisenGamma*first.CV);
    }
    else
    {
        return U(i,j,k,P)/(U(i,j,k,RHO_MIX,m,mixidx)*second.GruneisenGamma*second.CV);
    }
}


/*******************************************
 * Solid EOS
 ******************************************/

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







void WilkinsSolidEOS::define(Vector<Real> &params)
{
    adiabaticIndex 	= params[0];
    GruneisenGamma 	= adiabaticIndex-1.0;
    pref 			= params[1];
    eref	 		= params[2];
    CV 				= params[3];

    rho0            = params[4];
    e1              = params[5];
    e2              = params[6];
    e3              = params[7];
    e4              = params[8];
    e5              = params[9];
    G0              = params[10];
}

Real WilkinsSolidEOS::coldCompressionInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    Real x = U(i,j,k,RHO_K,m)/rho0;

    return U(i,j,k,ALPHARHO,m)*((e1+e2*x*x+e3*x+e4/x-e5*std::log(x))/rho0);
}

Real WilkinsSolidEOS::coldCompressionPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    Real x = U(i,j,k,RHO_K,m)/rho0;

    return (U(i,j,k,ALPHA,m)/GruneisenGamma)*( (2.0*e2*x*x*x+e3*x*x-e4-e5*x));
}

Real WilkinsSolidEOS::getSoundSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return (U(i,j,k,P)*(GruneisenGamma+1.0)-(coldCompressionPressure(U,i,j,k,m)+shearPressure(U,i,j,k,m)))/U(i,j,k,RHO_K,m) + dpcdrho(U,i,j,k,m) + dpsdrho(U,i,j,k,m) + (4.0/3.0)*componentShearModulus(U,i,j,k,m)/U(i,j,k,RHO_K,m);
}

Real WilkinsSolidEOS::componentShearModulus(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return G0;

    U(i,j,k,0,m) *=1.0;
}

Real WilkinsSolidEOS::dGdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;

    U(i,j,k,0,m) *=1.0;
}

Real WilkinsSolidEOS::dG2drho2(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;

    U(i,j,k,0,m) *=1.0;
}

Real WilkinsSolidEOS::transverseWaveSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
     return U(i,j,k,ALPHA,m)*componentShearModulus(U,i,j,k,m)/(U(i,j,k,RHO_K,m)*GruneisenGamma);
}

Real WilkinsSolidEOS::shearInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,ALPHA,m)*(componentShearModulus(U,i,j,k,m)*U(i,j,k,HJ2));
}

Real WilkinsSolidEOS::shearPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return (U(i,j,k,ALPHA,m)/GruneisenGamma)*((U(i,j,k,RHO_K,m)*dGdrho(U,i,j,k,m)-componentShearModulus(U,i,j,k,m))*U(i,j,k,HJ2));
}

Real WilkinsSolidEOS::dpcdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    Real x = U(i,j,k,RHO_K,m)/rho0;

    return (6.0*e2*x*x+2.0*e3*x-e5)/rho0;
}

Real WilkinsSolidEOS::dpsdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,RHO_K,m)*dG2drho2(U,i,j,k,m)*U(i,j,k,HJ2);
}

void WilkinsSolidEOS::setRhoFromDeformationTensor(BoxAccessCellArray& U, int i, int j, int k, int m, double* F)
{
    double determinant = det(F);

    U(i,j,k,RHO_K,m) = rho0/determinant;

}





/*void MieGruneisenSolidEOS::define(Vector<Real> &params)
{
    adiabaticIndex 	= params[0];
    GruneisenGamma 	= adiabaticIndex-1.0;
    pref 			= params[1];
    eref	 		= params[2];
    CV 				= params[3];

    rho0            = params[4];
    G0              = params[5];
}

Real MieGruneisenSolidEOS::getSoundSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return (U(i,j,k,P)*(GruneisenGamma+1.0)-(coldCompressionPressure(U,i,j,k,m)+shearPressure(U,i,j,k,m)))/U(i,j,k,RHO_K,m) + dpcdrho(U,i,j,k,m) + dpsdrho(U,i,j,k,m) + (4.0/3.0)*componentShearModulus(U,i,j,k,m)/U(i,j,k,RHO_K,m);
}

Real MieGruneisenSolidEOS::componentShearModulus(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return G0;

    U(i,j,k,RHO_K,m)*=1.0;
}

Real MieGruneisenSolidEOS::dGdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;

    U(i,j,k,RHO_K,m)*=1.0;
}

Real MieGruneisenSolidEOS::dG2drho2(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;

    U(i,j,k,RHO_K,m)*=1.0;
}

Real MieGruneisenSolidEOS::transverseWaveSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
     return U(i,j,k,ALPHA,m)*componentShearModulus(U,i,j,k,m)/(U(i,j,k,RHO_K,m)*GruneisenGamma);
}

Real MieGruneisenSolidEOS::shearInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,ALPHA,m)*(componentShearModulus(U,i,j,k,m)*U(i,j,k,HJ2));
}

Real MieGruneisenSolidEOS::shearPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return (U(i,j,k,ALPHA,m)/GruneisenGamma)*((U(i,j,k,RHO_K,m)*dGdrho(U,i,j,k,m)-componentShearModulus(U,i,j,k,m))*U(i,j,k,HJ2));
}

Real MieGruneisenSolidEOS::dpcdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;

    U(i,j,k,RHO_K,m)*=1.0;
}

Real MieGruneisenSolidEOS::dpsdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,RHO_K,m)*dG2drho2(U,i,j,k,m)*U(i,j,k,HJ2);
}

void MieGruneisenSolidEOS::setRhoFromDeformationTensor(BoxAccessCellArray& U, int i, int j, int k, int m, double* F)
{
    double determinant = det(F);

    U(i,j,k,RHO_K,m) = rho0/determinant;

}*/
