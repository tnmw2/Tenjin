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
        return 1E20*Heaviside<Real,Real>(U.accessPattern.materialInfo[m].EOS->componentShearModulus(U,i,j,k,m)*Jnew*sqrt(3.0/2.0)*2.0 - (*this).yieldStress[m]); //Changed!!! from sqrt(3) -> sqrt(3/2)*2
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
                    /*for(int m=0;m<U.numberOfMaterials;m++)
                    {
                        U(i,j,k,EPSILON,m)          = 0.0;
                        U(i,j,k,ALPHARHOEPSILON,m)  = 0.0;
                    }*/

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
                        //Print() << (293.0+U.accessPattern.materialInfo[m].EOS->getTemp(U,i,j,k,m,0) - 293.0)/(273.0+600.0-293.0) << std::endl;

                        J = sqrt(U(i,j,k,HJ2));

                        Jnew = bisection(U,i,j,k,J,parameters,m,dt);

                        /*if(overYieldStress(J,U,i,j,k,m))
                        {
                            //Jnew = yieldStress[m]/(sqrt(3.0)*U.accessPattern.materialInfo[m].EOS->componentShearModulus(U,i,j,k,m));

                            Jnew = bisection(U,i,j,k,J,parameters,m,dt);
                            */
                            /*if(parameters.PLASTIC==1)
                            {
                                Jnew = yieldStress[m]/(sqrt(3.0)*U.accessPattern.materialInfo[m].EOS->componentShearModulus(U,i,j,k,m));
                            }
                            else
                            {
                                Jnew = bisection(U,i,j,k,J,parameters,m,dt);
                            }*/
                        /*}
                        else
                        {
                            Jnew = J;
                        }*/

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

void matrixPrinter(double* U, int N =3)
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            Print() <<U[i*N+j] << "  ";
        }

        Print()<<std::endl;
    }

    Print()<<std::endl;

    return;
}

/** Loops over all cells performing the plastic update.
 */
void PlasticEOS::SVDboxPlasticUpdate(BoxAccessCellArray& U,ParameterStruct& parameters, Real dt)
{

    Real  JBefore;
    Real  JTotal;
    Real  norm;

    Real  J;
    Real  Jnew;

    Real  temp1[9];
    Real  temp2[9];
    Real  temp2inv[9];

    double devH[U.numberOfComponents*U.numberOfComponents];
    double Vmat[U.numberOfComponents*U.numberOfComponents];

    double h[U.numberOfComponents];
    double kmat[U.numberOfComponents];

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

                        Jnew = bisection(U,i,j,k,J,parameters,m,dt);


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



                U.amrexToArray(i,j,k,DEVH,0,devH);

                singularValueDecomposition(h,devH,temp1,temp2);

                for(int row = 0; row < U.numberOfComponents; row++)
                {
                    h[row]     = devH[row*3+row];
                    h[row]    *= JTotal/JBefore;
                    kmat[row]  = std::exp(h[row]);

                    /*if(kmat[row] < 0.0 || h[row] < 0.0 )
                    {
                        Abort("Negative stretch singular value");
                    }*/
                }

                for(int row = 0; row < U.numberOfComponents; row++)
                {
                    for(int col = 0; col < U.numberOfComponents; col++)
                    {
                        devH[row*3+col] = 0.0;
                        Vmat[row*3+col] = 0.0;


                        for(int a = 0; a < U.numberOfComponents; a++)
                        {
                            devH[row*3+col] += h[a];//*temp2[row*3+a]*temp1[a*3+col];
                            Vmat[row*3+col] += kmat[a]*delta<Real>(row,a)*delta<Real>(col,a);//*temp2[row*3+a]*temp1[a*3+col];
                        }

                        /*if(Vmat[row*3+col]<0.0 ||  Vmat[row*3+col] > 2.0)
                        {
                            Vector<Real> a(1);

                            a[0] = Vmat[row*3+col];

                            std::string err = "Stupid V";

                            customAbort(a,err);
                        }*/
                    }
                }



                for(int row =0;row<U.numberOfComponents;row++)
                {
                    for(int col =0;col<U.numberOfComponents;col++)
                    {

                        temp1[row*U.numberOfComponents+col]				 = delta<Real>(row,col)+0.5*devH[row*3+col];
                        temp2[row*U.numberOfComponents+col]				 = delta<Real>(row,col)-0.5*devH[row*3+col];

                    }
                }

                invert(temp2,temp2inv);

                squareMatrixMultiply(temp2inv,temp1,temp2);

                for(int row =0;row<U.numberOfComponents;row++)
                {
                    for(int col =0;col<U.numberOfComponents;col++)
                    {
                        //U(i,j,k,V_TENSOR,0,row,col) = temp2[row*U.numberOfComponents+col];

                        U(i,j,k,V_TENSOR,0,row,col) = Vmat[row*U.numberOfComponents+col];


                        //Print() << U(i,j,k,V_TENSOR,0,row,col) - Vmat[row*3+col] << " " << U(i,j,k,V_TENSOR,0,row,col) + Vmat[row*3+col] <<  " " << U(i,j,k,V_TENSOR,0,row,col) << " " << Vmat[row*3+col] <<  std::endl;

                    }
                }

                for(int row = 0; row < U.numberOfComponents; row++)
                {
                    for(int col = 0; col < U.numberOfComponents; col++)
                    {
                        if(Vmat[row*3+col]<0.0 ||  Vmat[row*3+col] > 2.0)
                        {
                            Print() << JTotal/JBefore << std::endl;

                            Print() << "What the svd way gives:" << std::endl;

                            U.amrexToArray(i,j,k,DEVH,0,devH);

                            singularValueDecomposition(h,devH,temp1,temp2);


                            for(int a = 0; a < U.numberOfComponents; a++)
                            {
                               Print() << kmat[a]  << std::endl;
                            }

                            matrixPrinter(temp1);
                            Print() << " " << std::endl;
                            matrixPrinter(temp2);

                            Print() << "What the old way gives:" << std::endl;

                            U.amrexToArray(i,j,k,V_TENSOR,0,Vmat);

                            singularValueDecomposition(kmat,Vmat,temp1,temp2);

                            for(int a = 0; a < U.numberOfComponents; a++)
                            {
                               Print() << kmat[a]  << std::endl;
                            }

                            matrixPrinter(temp1);
                            Print() << " " << std::endl;
                            matrixPrinter(temp2);

                            Print() << "They should agree on devH:" << std::endl;

                            matrixPrinter(devH);

                            singularValueDecomposition(h,devH,temp1,temp2);

                            for(int a = 0; a < U.numberOfComponents; a++)
                            {
                               Print() << h[a] << " " << std::exp(h[a]) << "  " << std::log(kmat[a])  << std::endl;
                            }

                            for(int row = 0; row < U.numberOfComponents; row++)
                            {
                                for(int col = 0; col < U.numberOfComponents; col++)
                                {
                                    devH[row*3+col] = 0.0;
                                    Vmat[row*3+col] = 0.0;


                                    for(int a = 0; a < U.numberOfComponents; a++)
                                    {
                                        devH[row*3+col] += h[a]*temp2[row*3+a]*temp1[a*3+col];
                                        Vmat[row*3+col] += kmat[a]*temp2[row*3+a]*temp1[a*3+col];
                                    }
                                }
                            }

                            matrixPrinter(devH);
                            Print() << " " << std::endl;
                            matrixPrinter(temp1);
                            Print() << " " << std::endl;
                            matrixPrinter(temp2);





                            Abort(" ");
                        }
                    }
                }

                /*for(int row = 0; row < U.numberOfComponents; row++)
                {
                    for(int col = 0; col < U.numberOfComponents; col++)
                    {
                        U(i,j,k,V_TENSOR,0,row,col) = 0.0;

                        for(int a = 0; a < U.numberOfComponents; a++)
                        {
                            U(i,j,k,V_TENSOR,0,row,col) += kmat[a]*temp2[row*3+a]*temp1[a*3+col];
                        }
                    }
                }*/
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
