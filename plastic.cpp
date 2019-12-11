#include "simulationheader.h"

PlasticEOS::PlasticEOS(){}

Real PlasticEOS::epsilonFunction(double J, double Jnew, BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,EPSILON,m) - sqrt(2.0/3.0)*(Jnew-J);
}

/** Calculate the plastic flow rate, \f$ \chi \f$.
 */
Real PlasticEOS::plasticStrainRate(double Jnew, double J, BoxAccessCellArray& U, int i, int j, int k, ParameterStruct& parameters,int m, Real Tstar)
{
    if(parameters.PLASTIC==1)
    {
        return 1E20*Heaviside<Real,Real>(U.accessPattern.materialInfo[m].EOS->componentShearModulus(U,i,j,k,m)*Jnew*sqrt(3.0) - (*this).yieldStress[m]);
    }
    else if(parameters.PLASTIC==2)                                                                                          //((1.0-std::pow(Tstar,mt))*
    {
        return 1.0*std::exp((sqrt(3.0/2.0)*2.0*U.accessPattern.materialInfo[m].EOS->componentShearModulus(U,i,j,k,m)*Jnew/(c1+c2*std::pow(epsilonFunction(J,Jnew,U,i,j,k,m),n))-1.0)/c3);
    }
    else if(parameters.PLASTIC==3)                                                                                          //((1.0-std::pow(Tstar,mt))*
    {
        return 1E20*Heaviside<Real,Real>(U.accessPattern.materialInfo[m].EOS->componentShearModulus(U,i,j,k,m)*Jnew*sqrt(3.0) - ((*this).yieldStress[m]+udaykumarConstant*epsilonFunction(J,Jnew,U,i,j,k,m)));
    }
    else
    {
        std::cout << "Incorrect plasticity type" << std::endl;

        exit(1);
    }
}

/** Function used in bisection of which we find the root.
 */
Real PlasticEOS::bisectionFunction(Real Jnew, Real J, BoxAccessCellArray& U, int i, int j, int k, ParameterStruct& parameters, int m, Real dt, Real Tstar)
{
    return Jnew-J+   sqrt(3.0/2.0)*  dt*plasticStrainRate(Jnew,J,U,i,j,k,parameters,m,Tstar);
}

Real PlasticEOS::bisection(BoxAccessCellArray& U, int i, int j, int k, Real J, ParameterStruct& parameters, int m, Real dt)
{
    Real JA = 0.0;
    Real JB = J;
    Real Jmid = J*0.5;

    Real tolerance = 1E-10; //1E-5;

    if(J < tolerance)
    {
        return 0.0;
    }

    Real Tstar = 0.0;

    if(parameters.PLASTIC == 2)
    {
       Tstar = std::max((293.0+U.accessPattern.materialInfo[m].EOS->getTemp(U,i,j,k,m,0) - 293.0)/(273.0+600.0-293.0),tolerance);

       //Print() << Tstar << std::endl;
    }

    if(sgn<Real,int>(bisectionFunction(JA,J,U,i,j,k,parameters,m,dt,Tstar)) == sgn<Real,int>(bisectionFunction(JB,J,U,i,j,k,parameters,m,dt,Tstar)) )
    {
        Print() << "Error in plastic Bisection " << std::endl;
        exit(1);
    }

    while(true)
    {
        if(sgn<Real,int>(bisectionFunction(Jmid,J,U,i,j,k,parameters,m,dt,Tstar)) == sgn<Real,int>(bisectionFunction(JA,J,U,i,j,k,parameters,m,dt,Tstar)))
        {
            JA = Jmid;
        }
        else
        {
            JB = Jmid;
        }

        Jmid = 0.5*(JA+JB);

        if(std::abs(JA-JB)/(JA+JB)<tolerance)
        {
            break;
        }
    }

    return Jmid;
}

/** Calculates the von-Mises yield stress criterion and returns whether the current state exceeds the yield surface.
 */
bool PlasticEOS::overYieldStress(Real& J, BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U.accessPattern.materialInfo[m].EOS->componentShearModulus(U,i,j,k,m)*J*sqrt(3.0) > (*this).yieldStress[m];
}

/** Loops over all cells performing the plastic update.
 */
void PlasticEOS::boxPlasticUpdate(BoxAccessCellArray& U,ParameterStruct& parameters, Real dt)
{

    Real  J;
    Real  Jnew;

    Real  temp1[9];
    Real  temp2[9];
    Real  temp2inv[9];

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    for(int m = 0; m < U.numberOfMaterials; m++)
    {
        if(U.accessPattern.materialInfo[m].phase == solid)
        {
            for 		(int k = lo.z; k <= hi.z; ++k)
            {
                for 	(int j = lo.y; j <= hi.y; ++j)
                {
                    for (int i = lo.x; i <= hi.x; ++i)
                    {

                        U.getHenckyJ2(i,j,k,m);

                        J = sqrt(U(i,j,k,HJ2,m));

                        Jnew = bisection(U,i,j,k,J,parameters,m,dt);


                        U(i,j,k,EPSILON,m)          += sqrt(2.0/3.0)*(J-Jnew);
                        U(i,j,k,RHOEPSILON,m)   = U(i,j,k,EPSILON,m)*U(i,j,k,RHO,m);

                        if(J == 0.0 && Jnew == 0.0 )
                        {
                            J    = 1.0;
                            Jnew = 1.0;
                        }


                        for(int row =0;row<U.numberOfComponents;row++)
                        {
                            for(int col =0;col<U.numberOfComponents;col++)
                            {
                                if(U(i,j,k,DEVH,m,row,col) != 0.0 && J == 0.0 && Jnew != 0.0)
                                {
                                    std::cout << "Divide by zero error" << std::endl;
                                    exit(1);
                                }
                                else
                                {
                                    U(i,j,k,DEVH,m,row,col)	*= Jnew/J;
                                }


                                temp1[row*U.numberOfComponents+col]				 = delta<Real>(row,col)+0.5*U(i,j,k,DEVH,m,row,col);
                                temp2[row*U.numberOfComponents+col]				 = delta<Real>(row,col)-0.5*U(i,j,k,DEVH,m,row,col);

                            }
                        }

                        invert(temp2,temp2inv);

                        squareMatrixMultiply(temp2inv,temp1,temp2);

                        for(int row =0;row<U.numberOfComponents;row++)
                        {
                            for(int col =0;col<U.numberOfComponents;col++)
                            {
                                U(i,j,k,V_TENSOR,m,row,col) = temp2[row*U.numberOfComponents+col];
                            }
                        }
                    }
                }
            }
        }
    }


    return;
}


/** Loops over all cells performing the plastic update.
 */
void PlasticEOS::plasticUpdate(CellArray& U, ParameterStruct& parameters, Real dt, MultiFab& S_new)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(S_new); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray  Ubox(mfi,bx,U);

        boxPlasticUpdate(Ubox,parameters,dt);

        Ubox.conservativeToPrimitive();
    }

    return;
}
