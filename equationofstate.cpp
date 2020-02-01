#include "simulationheader.h"

//EquationOfState::EquationOfState(){}

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
    //return U(i,j,k,P)*(GruneisenGamma+1.0)/U(i,j,k,RHO_K,m) - pref/U(i,j,k,RHO_K,m);

    return ((U(i,j,k,P))*(GruneisenGamma+1.0)-pref)/U(i,j,k,RHO_K,m);
}

Real MieGruneisenEOS::getTemp(BoxAccessCellArray& U, int i, int j, int k, int m, int mixidx)
{
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
 * JWL EOS
 ****************************************************/


void JWLEOS::define(Vector<Real> &params)
{

    adiabaticIndex 	= params[0];
    GruneisenGamma 	= adiabaticIndex-1.0;
    pref 			= params[1];
    eref	 		= params[2];
    CV 				= params[3];

    A               = params[4];
    B               = params[5];
    R1              = params[6];
    R2              = params[7];
    rho0            = params[8];

}

void JWLEOS::copy(JWLEOS& C)
{
    adiabaticIndex = C.adiabaticIndex;
    GruneisenGamma = C.GruneisenGamma;
    pref = C.pref;
    eref = C.eref;
    CV = C.CV;

    A = C.A;
    B = C.B;
    R1 = C.R1;
    R2 = C.R2;
    rho0 = C.rho0;
}

Real JWLEOS::coldCompressionInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,ALPHARHO,m)*(eref + (A/(rho0*R1))*std::exp(-R1*rho0/U(i,j,k,RHO_K,m)) + (B/(rho0*R2))*std::exp(-R2*rho0/U(i,j,k,RHO_K,m))  );
}

Real JWLEOS::coldCompressionPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return (U(i,j,k,ALPHA,m)/GruneisenGamma)*referencePressure(U,i,j,k,m);
}

Real JWLEOS::inverseGruneisen(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,ALPHA,m)/GruneisenGamma;
}

void JWLEOS::rootFind(BoxAccessCellArray& U, int i, int j, int k, int m, Real kineticEnergy)
{
    Abort("Trying to root find on a non-mixture EOS");
}

Real JWLEOS::getSoundSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    //return (U(i,j,k,P)-referencePressure(U,i,j,k,m) )*(GruneisenGamma+1.0)/U(i,j,k,RHO_K,m) + dpcdrho(U,i,j,k,m);
    return (U(i,j,k,P)*(GruneisenGamma+1.0) -referencePressure(U,i,j,k,m) )/U(i,j,k,RHO_K,m) + dprefdrho(U,i,j,k,m) - U(i,j,k,RHO_K,m)*GruneisenGamma*derefdrho(U,i,j,k,m);
}

Real JWLEOS::dprefdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    Real rho = U(i,j,k,RHO_K,m);

    return A*(R1*rho0/(rho*rho))*std::exp(-R1*rho0/rho) + B*(R2*rho0/(rho*rho))*std::exp(-R2*rho0/rho);
}

Real JWLEOS::derefdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    Real rho = U(i,j,k,RHO_K,m);

    return (referencePressure(U,i,j,k,m)-pref)/(rho*rho);
}

Real JWLEOS::getTemp(BoxAccessCellArray& U, int i, int j, int k, int m, int mixidx)
{
    return (U(i,j,k,P)- referencePressure(U,i,j,k,m))/(U(i,j,k,RHO_K,m)*GruneisenGamma*CV);
}

Real JWLEOS::xi(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 1.0/GruneisenGamma;
}

Real JWLEOS::referencePressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return (pref + A*std::exp(-R1*rho0/U(i,j,k,RHO_K,m)) + B*std::exp(-R2*rho0/U(i,j,k,RHO_K,m)));
}

Real JWLEOS::componentShearModulus(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real JWLEOS::dGdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real JWLEOS::dG2drho2(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real JWLEOS::transverseWaveSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real JWLEOS::shearInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real JWLEOS::shearPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

void JWLEOS::setRhoFromDeformationTensor(BoxAccessCellArray& U, int i, int j, int k, int m, double* F)
{
    return;
}

void JWLEOS::defineMixtureDensities(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return;
}


/*****************************************************
 * Mixture EOS
 ****************************************************/

Real  MixtureEOS::toleranceForSinglePhaseTreatment = 1E-3;
Real  MixtureEOS::toleranceForConvergence          = 1E-5;
Real  MixtureEOS::toleranceForBeingNearRoot        = 1E-3;


Real MixtureEOS::dGdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real MixtureEOS::dG2drho2(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real MixtureEOS::transverseWaveSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real MixtureEOS::shearInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real MixtureEOS::shearPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real MixtureEOS::componentShearModulus(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}


MixtureEOS::MixtureEOS()
{
    first.mixture  = true;
    second.mixture = true;

    first.mixtureIndex = 0;
    second.mixtureIndex = 1;

}

Real MixtureEOS::coldCompressionInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)
    {
        return first.coldCompressionInternalEnergy(U,i,j,k,m);
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        return second.coldCompressionInternalEnergy(U,i,j,k,m);
    }
    else
    {
        return U(i,j,k,LAMBDA,m)*first.coldCompressionInternalEnergy(U,i,j,k,m) + (1.0-U(i,j,k,LAMBDA,m))*second.coldCompressionInternalEnergy(U,i,j,k,m);
    }
}

Real MixtureEOS::coldCompressionPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)
    {
        return first.coldCompressionPressure(U,i,j,k,m);
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        return second.coldCompressionPressure(U,i,j,k,m);
    }
    else
    {
        return (first.coldCompressionPressure(U,i,j,k,m))*U(i,j,k,LAMBDA,m)*U(i,j,k,RHO_K,m)/(U(i,j,k,RHO_MIX,m,0))+(second.coldCompressionPressure(U,i,j,k,m))*(1.0-U(i,j,k,LAMBDA,m))*U(i,j,k,RHO_K,m)/(U(i,j,k,RHO_MIX,m,1));
    }
}

Real MixtureEOS::inverseGruneisen(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)
    {
        return first.inverseGruneisen(U,i,j,k,m);
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        return second.inverseGruneisen(U,i,j,k,m);
    }
    else
    {
        return U(i,j,k,ALPHARHOLAMBDA,m)/(first.GruneisenGamma*U(i,j,k,RHO_MIX,m,0))+(U(i,j,k,ALPHARHO,m)-U(i,j,k,ALPHARHOLAMBDA,m))/(second.GruneisenGamma*U(i,j,k,RHO_MIX,m,1));
    }

    //return U(i,j,k,ALPHA,m)/((U(i,j,k,LAMBDA,m)*first.adiabaticIndex*first.CV+(1.0-U(i,j,k,LAMBDA,m))*second.adiabaticIndex*second.CV)/(U(i,j,k,LAMBDA,m)*first.CV+(1.0-U(i,j,k,LAMBDA,m))*second.CV)-1.0);
    //return U(i,j,k,ALPHA,m)/first.GruneisenGamma;
    //return U(i,j,k,ALPHARHOLAMBDA,m)/(first.GruneisenGamma*U(i,j,k,RHO_MIX,m,0))+(U(i,j,k,ALPHARHO,m)-U(i,j,k,ALPHARHOLAMBDA,m))/(second.GruneisenGamma*U(i,j,k,RHO_MIX,m,1));
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

    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)
    {
        return first.getSoundSpeedContribution(U,i,j,k,m) ;
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        return second.getSoundSpeedContribution(U,i,j,k,m);
    }
    else
    {
        Real dedrho = (U(i,j,k,LAMBDA,m))*((first.pref -U(i,j,k,P))/(first.GruneisenGamma*U(i,j,k,RHO_K,m)*U(i,j,k,RHO_K,m)));

        dedrho +=   (1.0-U(i,j,k,LAMBDA,m))*((second.pref -U(i,j,k,P))/(second.GruneisenGamma*U(i,j,k,RHO_K,m)*U(i,j,k,RHO_K,m)));

        return (  U(i,j,k,P)/(U(i,j,k,RHO_K,m)*U(i,j,k,RHO_K,m))  - dedrho )/(U(i,j,k,LAMBDA,m)/(U(i,j,k,RHO_MIX,m,0)*first.GruneisenGamma)+(1.0-U(i,j,k,LAMBDA,m))/(U(i,j,k,RHO_MIX,m,1)*second.GruneisenGamma));
    }

}

Real MixtureEOS::xi(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)
    {
        return 1.0/first.GruneisenGamma;
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        return 1.0/second.GruneisenGamma;
    }
    else
    {
        return 1.0/((U(i,j,k,LAMBDA,m)*first.adiabaticIndex*first.CV+(1.0-U(i,j,k,LAMBDA,m))*second.adiabaticIndex*second.CV)/(U(i,j,k,LAMBDA,m)*first.CV+(1.0-U(i,j,k,LAMBDA,m))*second.CV)-1.0);
    }
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
    return std::min(std::max(U(i,j,k,LAMBDA,m)/(1.0/U(i,j,k,RHO_K,m)-(1.0-U(i,j,k,LAMBDA,m))/rhob),U(i,j,k,LAMBDA,m)*U(i,j,k,RHO_K,m)),1E4);
}

Real rhobFunc(BoxAccessCellArray& U, int i, int j, int k, int m, Real rhoa)
{
    return std::min(std::max((1.0-U(i,j,k,LAMBDA,m))/(1.0/U(i,j,k,RHO_K,m)-U(i,j,k,LAMBDA,m)/rhoa),(1.0-U(i,j,k,LAMBDA,m))*U(i,j,k,RHO_K,m)),1E4);
}

Real MixtureEOS::pressureFunction(BoxAccessCellArray& U, int i, int j, int k, int m, Real& p)
{
    return ((p-first.pref)/(U(i,j,k,RHO_MIX,m,0)*first.GruneisenGamma*first.CV))-((p-second.pref)/(U(i,j,k,RHO_MIX,m,1)*second.GruneisenGamma*second.CV));
}

Real MixtureEOS::bisectionFunction(BoxAccessCellArray& U, int i, int j, int k, int m, Real rhoTry, Real kineticEnergy, Real& p)
{
    U(i,j,k,RHO_MIX,m,0) = rhoTry;
    U(i,j,k,RHO_MIX,m,1) = rhobFunc(U,i,j,k,m,rhoTry);

    p = (U(i,j,k,TOTAL_E) - kineticEnergy - U.getEffectiveNonThermalInternalEnergy(i,j,k)+ U.getEffectiveNonThermalPressure(i,j,k))/(U.getEffectiveInverseGruneisen(i,j,k));

    return pressureFunction(U,i,j,k,m,p);
}

void MixtureEOS::rootFind(BoxAccessCellArray& U, int i, int j, int k, int m, Real kineticEnergy)
{

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

    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)
    {
        U(i,j,k,RHO_MIX,m,0) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,1) = rhobFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,0)); //U(i,j,k,RHO_K,m);

        return;
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        U(i,j,k,RHO_MIX,m,1) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,0) = rhoaFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,1));//U(i,j,k,RHO_K,m);

        return;
    }

    if(U(i,j,k,ALPHA,m) < 0.5)
    {
        U(i,j,k,RHO_MIX,m,0) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,1) = rhobFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,0)); //U(i,j,k,RHO_K,m);

        return;
    }



    Real A = U(i,j,k,LAMBDA,m)*U(i,j,k,RHO_K,m)+1E-10; //parameters.initialMixtureGuesses[0];
    Real B = 20000.0;

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
        std::string err = "Error in root Bisection";

        Vector<Real> vec;

        vec.push_back(A);
        vec.push_back(B);
        vec.push_back(bisectionFunction(U,i,j,k,m,A,kineticEnergy,p));
        vec.push_back(bisectionFunction(U,i,j,k,m,B,kineticEnergy,p));
        vec.push_back(U(i,j,k,ALPHA,m));


        customAbort(vec,err);

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
        return (U(i,j,k,P)-first.pref)/(U(i,j,k,RHO_MIX,m,mixidx)*first.GruneisenGamma*first.CV);
    }
    else
    {
        return (U(i,j,k,P)-second.pref)/(U(i,j,k,RHO_MIX,m,mixidx)*second.GruneisenGamma*second.CV);
    }
}

void MixtureEOS::defineMixtureDensities(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)
    {
        U(i,j,k,RHO_MIX,m,0) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,1) = U(i,j,k,RHO_K,m); //rhobFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,0));

    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        U(i,j,k,RHO_MIX,m,1) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,0) = U(i,j,k,RHO_K,m); //rhoaFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,1));

    }
    else
    {
        U(i,j,k,RHO_MIX,m,0) = U(i,j,k,RHO_K,m)*(U(i,j,k,LAMBDA,m)+(1.0-U(i,j,k,LAMBDA,m))*((U(i,j,k,P)-first.pref)/(U(i,j,k,P)-second.pref))*((second.CV*second.GruneisenGamma)/(first.CV*first.GruneisenGamma)));
        U(i,j,k,RHO_MIX,m,1) = rhobFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,0));
    }


    return;
}


/*****************************************************
 * JWL Mixture EOS
 ****************************************************/

Real  JWLMixtureEOS::toleranceForSinglePhaseTreatment = 1E-2;
Real  JWLMixtureEOS::toleranceForConvergence          = 1E-2;
Real  JWLMixtureEOS::toleranceForBeingNearRoot        = 1E-3;
Real  JWLMixtureEOS::maximumAllowableDensity          = 1E4;



Real JWLMixtureEOS::dGdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real JWLMixtureEOS::dG2drho2(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real JWLMixtureEOS::transverseWaveSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real JWLMixtureEOS::shearInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real JWLMixtureEOS::shearPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real JWLMixtureEOS::componentShearModulus(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}


JWLMixtureEOS::JWLMixtureEOS()
{
    first.mixture  = true;
    second.mixture = true;

    first.mixtureIndex = 0;
    second.mixtureIndex = 1;

}

Real JWLMixtureEOS::mixture_coldCompressionInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m, int idx)
{
    if(idx == 0)
    {
        return (first.eref + (first.A/(first.rho0*first.R1))*std::exp(-first.R1*first.rho0/U(i,j,k,RHO_MIX,m,0)) + (first.B/(first.rho0*first.R2))*std::exp(-first.R2*first.rho0/U(i,j,k,RHO_MIX,m,0))  );
    }
    else
    {
        return (second.eref + (second.A/(second.rho0*second.R1))*std::exp(-second.R1*second.rho0/U(i,j,k,RHO_MIX,m,1)) + (second.B/(second.rho0*second.R2))*std::exp(-second.R2*second.rho0/U(i,j,k,RHO_MIX,m,1))  );
    }
}

Real JWLMixtureEOS::mixture_coldCompressionPressure(BoxAccessCellArray& U, int i, int j, int k, int m, int idx)
{
    if(idx == 0)
    {
        return (1.0/(first.GruneisenGamma*U(i,j,k,RHO_MIX,m,0)))*mixture_referencePressure(U,i,j,k,m,0);
    }
    else
    {
        return (1.0/(second.GruneisenGamma*U(i,j,k,RHO_MIX,m,1)))*mixture_referencePressure(U,i,j,k,m,1);
    }
}

Real JWLMixtureEOS::mixture_referencePressure(BoxAccessCellArray& U, int i, int j, int k, int m, int idx)
{
    if(idx == 0)
    {
        return (first.pref + first.A*std::exp(-first.R1*first.rho0/U(i,j,k,RHO_MIX,m,0)) + first.B*std::exp(-first.R2*first.rho0/U(i,j,k,RHO_MIX,m,0)));
    }
    else
    {
        return (second.pref + second.A*std::exp(-second.R1*second.rho0/U(i,j,k,RHO_MIX,m,1)) + second.B*std::exp(-second.R2*second.rho0/U(i,j,k,RHO_MIX,m,1)));
    }
}

Real JWLMixtureEOS::coldCompressionInternalEnergy(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)// || U(i,j,k,ALPHA,m) < 0.5 )
    {
        return first.coldCompressionInternalEnergy(U,i,j,k,m);
        //return U(i,j,k,ALPHARHO,m)*mixture_coldCompressionInternalEnergy(U,i,j,k,m,0);
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        return second.coldCompressionInternalEnergy(U,i,j,k,m);
        //return U(i,j,k,ALPHARHO,m)*mixture_coldCompressionInternalEnergy(U,i,j,k,m,1);
    }
    else
    {
        //return U(i,j,k,ALPHA,m)*(U(i,j,k,LAMBDA,m)*mixture_coldCompressionInternalEnergy(U,i,j,k,m,0) +(1.0-U(i,j,k,LAMBDA,m))*mixture_coldCompressionInternalEnergy(U,i,j,k,m,1));
        return U(i,j,k,ALPHARHOLAMBDA,m)*mixture_coldCompressionInternalEnergy(U,i,j,k,m,0) + (U(i,j,k,ALPHARHO,m)-U(i,j,k,ALPHARHOLAMBDA,m))*mixture_coldCompressionInternalEnergy(U,i,j,k,m,1);
        //return U(i,j,k,LAMBDA,m)*first.coldCompressionInternalEnergy(U,i,j,k,m) + (1.0-U(i,j,k,LAMBDA,m))*second.coldCompressionInternalEnergy(U,i,j,k,m);
    }

}

Real JWLMixtureEOS::coldCompressionPressure(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)//  || U(i,j,k,ALPHA,m) < 0.5 )
    {
        //return U(i,j,k,ALPHARHO,m)*mixture_coldCompressionPressure(U,i,j,k,m,0);
        return first.coldCompressionPressure(U,i,j,k,m);
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        //return U(i,j,k,ALPHARHO,m)*mixture_coldCompressionPressure(U,i,j,k,m,1);
        return second.coldCompressionPressure(U,i,j,k,m);
    }
    else
    {
        return U(i,j,k,ALPHARHOLAMBDA,m)*mixture_coldCompressionPressure(U,i,j,k,m,0) + (U(i,j,k,ALPHARHO,m)-U(i,j,k,ALPHARHOLAMBDA,m))*mixture_coldCompressionPressure(U,i,j,k,m,1);
        //return (first.coldCompressionPressure(U,i,j,k,m))*U(i,j,k,LAMBDA,m)*U(i,j,k,RHO_K,m)/(U(i,j,k,RHO_MIX,m,0))+(second.coldCompressionPressure(U,i,j,k,m))*(1.0-U(i,j,k,LAMBDA,m))*U(i,j,k,RHO_K,m)/(U(i,j,k,RHO_MIX,m,1));
    }
}

Real JWLMixtureEOS::inverseGruneisen(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)//  || U(i,j,k,ALPHA,m) < 0.5 )
    {
        return first.inverseGruneisen(U,i,j,k,m);
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        return second.inverseGruneisen(U,i,j,k,m);
    }
    else
    {
        return U(i,j,k,ALPHARHOLAMBDA,m)/(first.GruneisenGamma*U(i,j,k,RHO_MIX,m,0))+(U(i,j,k,ALPHARHO,m)-U(i,j,k,ALPHARHOLAMBDA,m))/(second.GruneisenGamma*U(i,j,k,RHO_MIX,m,1));
    }
}

void JWLMixtureEOS::define(Vector<Real> &params)
{
    Vector<Real> _first;
    Vector<Real> _second;

    for(unsigned int i=0;i<params.size();i++)
    {
        if(i<9)
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

Real JWLMixtureEOS::getSoundSpeedContribution(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return mixtureSoundSpeed(U,i,j,k,m);
}

Real JWLMixtureEOS::mixtureSoundSpeed(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)//  || U(i,j,k,ALPHA,m) < 0.5 )
    {
        return first.getSoundSpeedContribution(U,i,j,k,m) ;
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        return second.getSoundSpeedContribution(U,i,j,k,m);
    }
    else
    {
        Real rho = U(i,j,k,RHO_MIX,m,0);
        Real prefa = (first.pref + first.A*std::exp(-first.R1*first.rho0/rho) + first.B*std::exp(-first.R2*first.rho0/rho)) ;
        Real dpca = first.A*(first.R1*first.rho0/(rho*rho))*std::exp(-first.R1*first.rho0/rho) + first.B*(first.R2*first.rho0/(rho*rho))*std::exp(-first.R2*first.rho0/rho);
        Real decdrhoa = (prefa-first.pref)/(rho*rho);

        rho = U(i,j,k,RHO_MIX,m,1);
        Real prefb = mixture_referencePressure(U,i,j,k,m,1);
        Real dpcb = second.A*(second.R1*second.rho0/(rho*rho))*std::exp(-second.R1*second.rho0/rho) + second.B*(second.R2*second.rho0/(rho*rho))*std::exp(-second.R2*second.rho0/rho);
        Real decdrhob = (prefb-second.pref)/(rho*rho);

        Real dedrho = (U(i,j,k,LAMBDA,m)*( prefa - U(i,j,k,P) - U(i,j,k,RHO_MIX,m,0)*dpca  + U(i,j,k,RHO_MIX,m,0)*U(i,j,k,RHO_MIX,m,0)*first.GruneisenGamma*decdrhoa )/(first.GruneisenGamma) + (1.0-U(i,j,k,LAMBDA,m))*( prefb - U(i,j,k,P) - U(i,j,k,RHO_MIX,m,1)*dpcb + U(i,j,k,RHO_MIX,m,1)*U(i,j,k,RHO_MIX,m,1)*second.GruneisenGamma*decdrhob )/(second.GruneisenGamma))/(U(i,j,k,RHO_K,m)*U(i,j,k,RHO_K,m));
        //Real dedrho = (( prefa - U(i,j,k,P) - U(i,j,k,RHO_MIX,m,0)*dpca  + U(i,j,k,RHO_MIX,m,0)*U(i,j,k,RHO_MIX,m,0)*first.GruneisenGamma*decdrhoa )/(first.GruneisenGamma) + ( prefb - U(i,j,k,P) - U(i,j,k,RHO_MIX,m,1)*dpcb + U(i,j,k,RHO_MIX,m,1)*U(i,j,k,RHO_MIX,m,1)*second.GruneisenGamma*decdrhob )/(second.GruneisenGamma))/(U(i,j,k,RHO_K,m)*U(i,j,k,RHO_K,m));

        //Real dedrho = (U(i,j,k,LAMBDA,m)*( (first.GruneisenGamma+1.0)*prefa - U(i,j,k,P) - U(i,j,k,RHO_MIX,m,0)*dpca)/(first.GruneisenGamma) + (1.0-U(i,j,k,LAMBDA,m))*( (second.GruneisenGamma+1.0)*prefb - U(i,j,k,P) - U(i,j,k,RHO_MIX,m,1)*dpcb)/(second.GruneisenGamma))/(U(i,j,k,RHO_K,m)*U(i,j,k,RHO_K,m));

        return (  U(i,j,k,P)/(U(i,j,k,RHO_K,m)*U(i,j,k,RHO_K,m))  - dedrho )/(U(i,j,k,LAMBDA,m)/(U(i,j,k,RHO_MIX,m,0)*first.GruneisenGamma)+(1.0-U(i,j,k,LAMBDA,m))/(U(i,j,k,RHO_MIX,m,1)*second.GruneisenGamma));
    }

}

Real JWLMixtureEOS::xi(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)//  || U(i,j,k,ALPHA,m) < 0.5 )
    {
        return 1.0/first.GruneisenGamma;
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        return 1.0/second.GruneisenGamma;
    }
    else
    {
        return 1.0/((U(i,j,k,LAMBDA,m)*first.adiabaticIndex*first.CV+(1.0-U(i,j,k,LAMBDA,m))*second.adiabaticIndex*second.CV)/(U(i,j,k,LAMBDA,m)*first.CV+(1.0-U(i,j,k,LAMBDA,m))*second.CV)-1.0);
    }
}


/******************************************************
 * Root Finding stuff
 *****************************************************/


Real JWLMixtureEOS::pressureFunction(BoxAccessCellArray& U, int i, int j, int k, int m, Real& p)
{
    Real prefa = mixture_referencePressure(U,i,j,k,m,0);

    Real prefb = mixture_referencePressure(U,i,j,k,m,1);


    return ((p-prefa)/(U(i,j,k,RHO_MIX,m,0)*first.GruneisenGamma*first.CV))-((p-prefb)/(U(i,j,k,RHO_MIX,m,1)*second.GruneisenGamma*second.CV));
}

Real JWLMixtureEOS::bisectionFunction(BoxAccessCellArray& U, int i, int j, int k, int m, Real rhoTry, Real kineticEnergy, Real& p)
{
    U(i,j,k,RHO_MIX,m,0) = rhoTry;
    U(i,j,k,RHO_MIX,m,1) = rhobFunc(U,i,j,k,m,rhoTry);

    p = (U(i,j,k,TOTAL_E) - kineticEnergy - U.getEffectiveNonThermalInternalEnergy(i,j,k)+ U.getEffectiveNonThermalPressure(i,j,k))/(U.getEffectiveInverseGruneisen(i,j,k));

    return pressureFunction(U,i,j,k,m,p);
}

void JWLMixtureEOS::rootFind(BoxAccessCellArray& U, int i, int j, int k, int m, Real kineticEnergy)
{

    if(std::isnan(U(i,j,k,LAMBDA,m)))
    {
        Print() << "Nan in Lambda in root finding" << std::endl;
    }

    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)
    {
        U(i,j,k,RHO_MIX,m,0) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,1) = rhobFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,0));

        return;
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        U(i,j,k,RHO_MIX,m,1) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,0) = rhoaFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,1));
        return;
    }

    if(U(i,j,k,ALPHA,m) < toleranceForSinglePhaseTreatment)
    {
        U(i,j,k,RHO_MIX,m,1) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,0) = rhoaFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,1));
        return;
    }



    Real A = U(i,j,k,LAMBDA,m)*U(i,j,k,RHO_K,m)+1E-5; //parameters.initialMixtureGuesses[0];
    Real B = maximumAllowableDensity;

    Real p;

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



    for(int it = 0; it < 50; it++)
    {
        int stop = 0;

        if(sgn(bisectionFunction(U,i,j,k,m,A,kineticEnergy,p)) == sgn(bisectionFunction(U,i,j,k,m,B,kineticEnergy,p)) )
        {
            for(int s = 0; s < 1000; s++)
            {
                Real temp = (B+(A-B)*((Real)s/1000.0));

                if(sgn(bisectionFunction(U,i,j,k,m,A,kineticEnergy,p)) != sgn(bisectionFunction(U,i,j,k,m,temp,kineticEnergy,p)))
                {
                    B = temp;
                    stop = 1;
                    break;
                }
            }

            if(stop)
            {
                break;
            }
            else
            {
                B = A + (B-A)/2.0;

                //B *= 2.0;
            }
        }
        else
        {
            break;
        }
    }

    Real mid = 0.5*(A+B);


    if(sgn(bisectionFunction(U,i,j,k,m,A,kineticEnergy,p)) == sgn(bisectionFunction(U,i,j,k,m,B,kineticEnergy,p)) )
    {
        std::string err = "Error in root Bisection, even after checking";

        Vector<Real> vec;

        vec.push_back(A);
        vec.push_back(B);
        vec.push_back(bisectionFunction(U,i,j,k,m,A,kineticEnergy,p));
        vec.push_back(bisectionFunction(U,i,j,k,m,B,kineticEnergy,p));
        vec.push_back(U(i,j,k,ALPHA,m));
        vec.push_back(U(i,j,k,LAMBDA,m));
        vec.push_back(U(i,j,k,P));
        customAbort(vec,err);

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
            //Print() << n << std::endl;
            return;
        }

        if(n>1000)
        {
            Print() << "Root finding is taking a while" << std::endl;
        }
    }
}

Real JWLMixtureEOS::getTemp(BoxAccessCellArray& U, int i, int j, int k, int m, int mixidx)
{
    if(mixidx == 0)
    {
        Real prefa = mixture_referencePressure(U,i,j,k,m,0);

        return (U(i,j,k,P)-prefa)/(U(i,j,k,RHO_MIX,m,mixidx)*first.GruneisenGamma*first.CV);
    }
    else
    {
        Real prefb = mixture_referencePressure(U,i,j,k,m,1);

        return (U(i,j,k,P)-prefb)/(U(i,j,k,RHO_MIX,m,mixidx)*second.GruneisenGamma*second.CV);
    }
}

Real JWLMixtureEOS::primitiveBisectionFunction(BoxAccessCellArray& U, int i, int j, int k, int m, Real& rhoaGuess)
{
    Real rhoa = rhoaGuess;
    Real prefa = (first.pref + first.A*std::exp(-first.R1*first.rho0/rhoa) + first.B*std::exp(-first.R2*first.rho0/rhoa)) ;

    Real rhob = rhobFunc(U,i,j,k,m,rhoaGuess);
    Real prefb = (second.pref + second.A*std::exp(-second.R1*second.rho0/rhob) + second.B*std::exp(-second.R2*second.rho0/rhob)) ;


    return ((U(i,j,k,P)-prefa)/(rhoa*first.GruneisenGamma*first.CV))-((U(i,j,k,P)-prefb)/(rhob*second.GruneisenGamma*second.CV));
}

void JWLMixtureEOS::defineMixtureDensities(BoxAccessCellArray& U, int i, int j, int k, int m)
{

    if(U(i,j,k,LAMBDA,m) > 1.0-toleranceForSinglePhaseTreatment)
    {
        U(i,j,k,RHO_MIX,m,0) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,1) = rhobFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,0));

        return;
    }
    else if(U(i,j,k,LAMBDA,m) < toleranceForSinglePhaseTreatment)
    {
        U(i,j,k,RHO_MIX,m,1) = U(i,j,k,RHO_K,m);
        U(i,j,k,RHO_MIX,m,0) = rhoaFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,1));
        return;
    }
    else
    {
        if(U(i,j,k,ALPHA,m) < toleranceForSinglePhaseTreatment)
        {
            U(i,j,k,RHO_MIX,m,1) = U(i,j,k,RHO_K,m);
            U(i,j,k,RHO_MIX,m,0) = rhoaFunc(U,i,j,k,m,U(i,j,k,RHO_MIX,m,1));
            return;
        }


        Real A = U(i,j,k,LAMBDA,m)*U(i,j,k,RHO_K,m)+1E-5; //parameters.initialMixtureGuesses[0];
        Real B = maximumAllowableDensity;



        for(int it = 0; it < 50; it++)
        {
            int stop = 0;

            if(sgn(primitiveBisectionFunction(U,i,j,k,m,A)) == sgn(primitiveBisectionFunction(U,i,j,k,m,B)) )
            {
                for(int s = 0; s < 1000; s++)
                {
                    Real temp = (B+(A-B)*((Real)s/1000.0));

                    if(sgn(primitiveBisectionFunction(U,i,j,k,m,A)) != sgn(primitiveBisectionFunction(U,i,j,k,m,temp)))
                    {
                        B = temp;
                        stop = 1;
                        break;
                    }
                }

                if(stop)
                {
                    break;
                }
                else
                {
                    B = A + (B-A)/2.0;
                    //B *= 2.0;
                }
            }
            else
            {
                break;
            }
        }

        Real mid = 0.5*(A+B);

        if(sgn(primitiveBisectionFunction(U,i,j,k,m,A)) == sgn(primitiveBisectionFunction(U,i,j,k,m,B)) )
        {

            std::string err = "Error in primitive root Bisection";

            Vector<Real> vec;

            vec.push_back(A);
            vec.push_back(B);
            vec.push_back(primitiveBisectionFunction(U,i,j,k,m,A));
            vec.push_back(primitiveBisectionFunction(U,i,j,k,m,B));
            vec.push_back(U(i,j,k,ALPHA,m));
            vec.push_back(U(i,j,k,LAMBDA,m));
            vec.push_back(U(i,j,k,P));

            customAbort(vec,err);

        }



        for(int n=0;;n++)
        {

            if(sgn(primitiveBisectionFunction(U,i,j,k,m,mid)) == sgn(primitiveBisectionFunction(U,i,j,k,m,A)))
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


    return;
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
    return (U(i,j,k,P)*(GruneisenGamma+1.0)-(coldCompressionPressure(U,i,j,k,m)+shearPressure(U,i,j,k,m))*(GruneisenGamma+1.0))/U(i,j,k,RHO_K,m) + dpcdrho(U,i,j,k,m) + dpsdrho(U,i,j,k,m) + (4.0/3.0)*componentShearModulus(U,i,j,k,m)/U(i,j,k,RHO_K,m);
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

void RomenskiiSolidEOS::defineMixtureDensities(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return;
}


/*******************************************
 * Wilkins EOS
 ******************************************/

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
    return (U(i,j,k,P)*(GruneisenGamma+1.0)-(coldCompressionPressure(U,i,j,k,m)+shearPressure(U,i,j,k,m))*(GruneisenGamma+1.0))/U(i,j,k,RHO_K,m) + dpcdrho(U,i,j,k,m) + dpsdrho(U,i,j,k,m) + (4.0/3.0)*componentShearModulus(U,i,j,k,m)/U(i,j,k,RHO_K,m);
}

Real WilkinsSolidEOS::componentShearModulus(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return G0;
}

Real WilkinsSolidEOS::dGdrho(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
}

Real WilkinsSolidEOS::dG2drho2(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return 0.0;
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

Real WilkinsSolidEOS::inverseGruneisen(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,ALPHA,m)/GruneisenGamma;
}

void WilkinsSolidEOS::defineMixtureDensities(BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return;
}

Real WilkinsSolidEOS::getTemp(BoxAccessCellArray& U, int i, int j, int k, int m, int mixidx)
{
    Real x = U(i,j,k,RHO_K,m)/rho0;

    return (U(i,j,k,P)-(2.0*e2*x*x*x+e3*x*x-e4-e5*x))/(U(i,j,k,RHO_K,m)*GruneisenGamma*CV);
}
