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

MixtureEOS::MixtureEOS()
{
    first.mixture  = true;
    second.mixture = true;

    first.mixtureIndex = 0;
    second.mixtureIndex = 1;
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
    //return U(i,j,k,ALPHA,m)/((U(i,j,k,LAMBDA,m)*first.adiabaticIndex*first.CV+(1.0-U(i,j,k,LAMBDA,m))*second.adiabaticIndex*second.CV)/(U(i,j,k,LAMBDA,m)*first.CV+(1.0-U(i,j,k,LAMBDA,m))*second.CV)-1.0);
    //return U(i,j,k,ALPHA,m)/first.GruneisenGamma;
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
    static const Real toleranceForSinglePhaseTreatment = 1E-3;
    static const Real toleranceForConvergence = 1E-7;
    static const Real toleranceForBeingNearRoot = 1E-4;

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


    Real A = U(i,j,k,LAMBDA,m)*U(i,j,k,RHO_K,m)+1E-10; //parameters.initialMixtureGuesses[0];
    Real B = 1000000.0;

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

        Print() << " A: " << A << " B: " << B << std::endl;
        Print() << bisectionFunction(U,i,j,k,m,A,kineticEnergy,p) << " " << bisectionFunction(U,i,j,k,m,B,kineticEnergy,p)<< std::endl;
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
