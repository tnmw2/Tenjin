#include "simulationheader.h"

PlasticEOS::PlasticEOS(){}

Real PlasticEOS::epsilonFunction(double J, double Jnew, BoxAccessCellArray& U, int i, int j, int k, int m)
{
    return U(i,j,k,EPSILON,m) - sqrt(2.0/3.0)*(Jnew-J);
}

/** Calculate the plastic flow rate, \f$ \chi \f$.
 */
Real PlasticEOS::plasticStrainRate(double Jnew, double J, BoxAccessCellArray& U, int i, int j, int k, ParameterStruct& parameters,int m)
{
    if(parameters.PLASTIC==1)
    {
        return 1E20*Heaviside<Real,Real>(U.accessPattern.materialInfo[m].EOS->componentShearModulus(U,i,j,k,m)*Jnew*sqrt(3.0) - (*this).yieldStress[m]);
    }
    else if(parameters.PLASTIC==2)
    {
        return std::exp((1.0/c3)*(sqrt(3.0/2.0)*2.0*U.accessPattern.materialInfo[m].EOS->componentShearModulus(U,i,j,k,m)*Jnew/(c1+c2*std::pow(epsilonFunction(J,Jnew,U,i,j,k,m),n))-1.0));
    }
    else
    {
        std::cout << "Incorrect plasticity type" << std::endl;

        exit(1);
    }
}

/** Function used in bisection of which we find the root.
 */
Real PlasticEOS::bisectionFunction(Real Jnew, Real J, BoxAccessCellArray& U, int i, int j, int k, ParameterStruct& parameters, int m)
{
    return Jnew-J+parameters.dt*plasticStrainRate(Jnew,J,U,i,j,k,parameters,m);
}

Real PlasticEOS::bisection(BoxAccessCellArray& U, int i, int j, int k, Real J, ParameterStruct& parameters, int m)
{
    Real JA = 0.0;
    Real JB = J;
    Real Jmid = J*0.5;

    Real tolerance = 1E-5;

    if(sgn<Real,int>(bisectionFunction(JA,J,U,i,j,k,parameters,m)) == sgn<Real,int>(bisectionFunction(JB,J,U,i,j,k,parameters,m)) )
    {
        Print() << "Error in plastic Bisection " << std::endl;
        exit(1);
    }


    while(true)
    {
        if(sgn<Real,int>(bisectionFunction(Jmid,J,U,i,j,k,parameters,m)) == sgn<Real,int>(bisectionFunction(JA,J,U,i,j,k,parameters,m)))
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
void PlasticEOS::boxPlasticUpdate(BoxAccessCellArray& U,ParameterStruct& parameters)
{

    Real  JBefore;
    Real  JTotal;
    Real  norm;

    Real  J;
    Real  Jnew;

    Real  temp1[9];
    Real  temp2[9];
    Real  temp2inv[9];

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    for 		(int k = lo.z; k <= hi.z; ++k)
    {
        for 	(int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                if(U.cellIsMostlyFluid(i,j,k))
                {
                    for(int m=0;m<U.numberOfMaterials;m++)
                    {
                        U(i,j,k,EPSILON,m)          = 0.0;
                        U(i,j,k,ALPHARHOEPSILON,m)  = 0.0;
                    }

                    continue;
                }

                U.getHenckyJ2(i,j,k);

                JTotal  = 0.0;
                JBefore = 0.0;
                norm    = 0.0;


                for(int m=0;m<U.numberOfMaterials;m++)
                {
                    if(U.accessPattern.materialInfo[m].phase == solid)
                    {

                        J = sqrt(U(i,j,k,HJ2));


                        if(overYieldStress(J,U,i,j,k,m))
                        {
                            Jnew = yieldStress[m]/(sqrt(3.0)*U.accessPattern.materialInfo[m].EOS->componentShearModulus(U,i,j,k,m));

                            //Jnew = bisection(U,i,j,k,J,parameters,m);

                            /*if(parameters.PLASTIC==1)
                            {
                                U(i,j,k).Jnew[m] = parameters.yieldStress[m]/(sqrt3*U(i,j,k).componentShearModulus(parameters,m));
                            }
                            else
                            {
                                U(i,j,k).Jnew[m] = bisection(U(i,j,k),U(i,j,k).J[m],parameters,m);
                            }*/
                        }
                        else
                        {
                            Jnew = J;
                        }

                        U(i,j,k,EPSILON,m)          += sqrt(2.0/3.0)*(J-Jnew);
                        U(i,j,k,ALPHARHOEPSILON,m)   = U(i,j,k,EPSILON,m)*U(i,j,k,ALPHARHO,m);
                    }
                    else
                    {
                        J       = 0.0;
                        Jnew    = 0.0;

                        U(i,j,k,EPSILON,m)          = 0.0;
                        U(i,j,k,ALPHARHOEPSILON,m)  = 0.0;
                    }



                    /*if(sqrt(2.0/3.0)*(J-Jnew)  > 0.0)
                    {
                        Print() << sqrt(2.0/3.0)*(J-Jnew) << std::endl;
                    }*/

                    JBefore += (U.accessPattern.materialInfo[m].EOS->inverseGruneisen(U,i,j,k,m))*J;
                    JTotal  += (U.accessPattern.materialInfo[m].EOS->inverseGruneisen(U,i,j,k,m))*Jnew;
                    norm    += (U.accessPattern.materialInfo[m].EOS->inverseGruneisen(U,i,j,k,m));
                }

                JBefore *= 1.0/norm;
                JTotal  *= 1.0/norm;

                if(JBefore == 0.0 && JTotal == 0.0 )
                {
                    JBefore = 1.0;
                    JTotal = 1.0;
                }

                for(int row =0;row<U.numberOfComponents;row++)
                {
                    for(int col =0;col<U.numberOfComponents;col++)
                    {
                        if(U(i,j,k,DEVH,0,row,col) != 0.0 && JBefore == 0.0 && JTotal != 0.0)
                        {
                            std::cout << "Divide by zero error" << std::endl;
                            exit(1);
                        }
                        else
                        {
                            U(i,j,k,DEVH,0,row,col)	*= JTotal/JBefore;
                        }


                        temp1[row*U.numberOfComponents+col]				 = delta<Real>(row,col)+0.5*U(i,j,k,DEVH,0,row,col);
                        temp2[row*U.numberOfComponents+col]				 = delta<Real>(row,col)-0.5*U(i,j,k,DEVH,0,row,col);

                    }
                }

                invert(temp2,temp2inv);

                squareMatrixMultiply(temp2inv,temp1,temp2);

                for(int row =0;row<U.numberOfComponents;row++)
                {
                    for(int col =0;col<U.numberOfComponents;col++)
                    {
                        U(i,j,k,V_TENSOR,0,row,col) = temp2[row*U.numberOfComponents+col];
                    }
                }
            }
        }
    }


    return;
}

/** Loops over all cells performing the plastic update.
 */
void PlasticEOS::plasticUpdate(CellArray& U, ParameterStruct& parameters)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(U.data); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray  Ubox(mfi,bx,U);

        boxPlasticUpdate(Ubox,parameters);

        Ubox.conservativeToPrimitive();
    }

    return;
}
