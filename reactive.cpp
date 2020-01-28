#include "simulationheader.h"

void pressureBased_UpdateMassFraction(BoxAccessCellArray& U, BoxAccessCellArray& U1, ParameterStruct& parameters, Real dt, int i, int j, int k, int m)
{

    Real sigma = 20.0;//      //20.0   ;//10.0;
    Real nu    = 0.5;
    Real n     = 2.5;

    if(U(i,j,k,ALPHA,m) < 0.5)
    {
        U1(i,j,k,ALPHARHOLAMBDA,m) = U(i,j,k,ALPHARHOLAMBDA,m);
        return;
    }


    U1(i,j,k,ALPHARHOLAMBDA,m) = U(i,j,k,ALPHARHOLAMBDA,m)-dt*(sigma*std::pow(U(i,j,k,ALPHARHO,m),1.0-nu)*std::pow(U(i,j,k,ALPHARHOLAMBDA,m),nu)*std::pow((U(i,j,k,P)/1E6),n));

    if(U1(i,j,k,ALPHARHOLAMBDA,m) < 0.0)
    {
        U1(i,j,k,ALPHARHOLAMBDA,m) = 0.0;
    }

}

void Arrhenius_UpdateMassFraction(BoxAccessCellArray& U, BoxAccessCellArray& U1, ParameterStruct& parameters, Real dt, int i, int j, int k, int m)
{

    Real Tc = 150.0; //Reference temp;

    Real A = 2E6;    //Pre-exponential factor

    Real T = U.accessPattern.materialInfo[m].EOS->getTemp(U,i,j,k,m,0);

    if(T<=0.0)
    {
        T = 1E-20;
    }

    Real temp = -A*U(i,j,k,ALPHARHOLAMBDA,m)*std::exp(-Tc/T);

    if(std::isnan(temp))
    {
        temp = 0.0;
    }

    U1(i,j,k,ALPHARHOLAMBDA,m) = U(i,j,k,ALPHARHOLAMBDA,m)+dt*temp;

    if(U1(i,j,k,ALPHARHOLAMBDA,m) < 0.0)
    {
        U1(i,j,k,ALPHARHOLAMBDA,m) = 0.0;
    }


}

void pressureBased_UpdateMassFraction_single(BoxAccessCellArray& U, ParameterStruct& parameters, Real dt, int i, int j, int k, int m)
{

    /*Real sigma = 20.0;
    Real nu    = 0.5;
    Real n     = 1.5;*/

    Real sigma = 20.0;
    Real nu    = 0.5;
    Real n     = 1.5;




    if(U(i,j,k,ALPHARHOLAMBDA,m) < 0.0)
    {
        U(i,j,k,ALPHARHOLAMBDA,m) = 1E-6;
    }

    U(i,j,k,ALPHARHOLAMBDA,m) -=  dt*(sigma*std::pow(U(i,j,k,ALPHARHO,m),1.0-nu)*std::pow(U(i,j,k,ALPHARHOLAMBDA,m),nu)*std::pow( ((U(i,j,k,P) > 1E6 ? U(i,j,k,P) : 0.0)/1E6),n));



    if(U(i,j,k,ALPHARHOLAMBDA,m) < 0.0 || std::isnan(U(i,j,k,ALPHARHOLAMBDA,m)))
    {
        U(i,j,k,ALPHARHOLAMBDA,m) = 1E-6;
    }
}

void IandG_UpdateMassFraction_single(BoxAccessCellArray& U, ParameterStruct& parameters, Real dt, int i, int j, int k, int m)
{
    if(U(i,j,k,ALPHARHOLAMBDA,m) < 0.0)
    {
        U(i,j,k,ALPHARHOLAMBDA,m) = 1E-6;
    }

    if(U(i,j,k,LAMBDA,m) < 1E-6)
    {
        return;
    }

    /*static Real rho0 = 1905.0;
    static Real I  = 4E12;
    static Real G1 = 4500E-27;
    static Real G2 = 30E-9;
    static Real a = 0.22;
    static Real b = 0.667;
    static Real c = 0.667;
    static Real d = 1.0;
    static Real e = 0.667;
    static Real g = 0.667;
    static int  X = 7;
    static int  Y = 3;
    static int  Z = 1;
    static Real phi_ig = 0.02;
    static Real phi_G1 = 0.8;
    static Real phi_G2 = 0.8;*/

    Real rho0 = 1905.0;
    Real I  = 4E12;
    Real G1 = 4500.0;//E-27;
    Real G2 = 0.3E6;
    Real a = 0.22;
    Real b = 0.667;
    Real c = 0.667;
    Real d = 1.0;
    Real e = 0.667;
    Real g = 0.667;
    int  X = 7;
    int  Y = 3;
    int  Z = 1;
    Real phi_ig = 0.02;
    Real phi_G1 = 0.8;
    Real phi_G2 = 0.8;

    Real phi = 1.0 - U(i,j,k,LAMBDA,m);

    Real r = 0.0;

    if((phi_ig-phi) > 0.0 && U(i,j,k,RHO_K,m)/rho0-(1.0+a) > 0.0)
    {
        r += I*std::pow(1.0-phi,b)*std::pow(U(i,j,k,RHO_K,m)/rho0-(1.0+a),X);
    }

    if( phi_G1-phi > 1E-6)
    {
        r += G1*std::pow(1.0-phi,c)*std::pow(phi,d)*std::pow(U(i,j,k,P)/1E9,Y);
    }

    if(phi-phi_G2 > 0.0)
    {
        r += G2*std::pow(1.0-phi,e)*std::pow(phi,g)*std::pow(U(i,j,k,P)/1E9,Z);
    }

    r = std::max(0.0,r);

    U(i,j,k,ALPHARHOLAMBDA,m) -=  dt*U(i,j,k,ALPHARHO,m)*r;

    if(U(i,j,k,ALPHARHOLAMBDA,m) < 0.0 || std::isnan(U(i,j,k,ALPHARHOLAMBDA,m)))
    {
        U(i,j,k,ALPHARHOLAMBDA,m) = 1E-6;
    }
}

void Schoch_ProgrammedBurn_UpdateMassFraction_single(BoxAccessCellArray& U, ParameterStruct& parameters, Real dt, int i, int j, int k, int m)
{

    Real A = 2.0E6;

    if(U(i,j,k,ALPHARHOLAMBDA,m) < 0.0)
    {
        U(i,j,k,ALPHARHOLAMBDA,m) = 1E-6;
    }

    //U(i,j,k,ALPHARHOLAMBDA,m) -=  dt*A*std::pow(U(i,j,k,ALPHARHOLAMBDA,m)*U(i,j,k,ALPHARHOLAMBDA,m)/U(i,j,k,ALPHARHO,m),0.75)*(U(i,j,k,P) > 1E8 ? 1.0 : 0.0);
    U(i,j,k,ALPHARHOLAMBDA,m) -=  dt*A*U(i,j,k,ALPHARHOLAMBDA,m)*std::pow(U(i,j,k,ALPHARHOLAMBDA,m)/U(i,j,k,ALPHARHO,m),0.5)*(U(i,j,k,P) > 9.9E8 ? 1.0 : 0.0);

    if(U(i,j,k,ALPHARHOLAMBDA,m) < 0.0 || std::isnan(U(i,j,k,ALPHARHOLAMBDA,m)))
    {
        U(i,j,k,ALPHARHOLAMBDA,m) = 1E-6;
    }
}

void Arrhenius_UpdateMassFraction_single(BoxAccessCellArray& U, ParameterStruct& parameters, Real dt, int i, int j, int k, int m)
{

    Real Tc = 210.0;   //Reference temp;
    Real A = 2.6E6;    //Pre-exponential factor
    Real T = U.accessPattern.materialInfo[m].EOS->getTemp(U,i,j,k,m,0);

    if(U(i,j,k,ALPHARHOLAMBDA,m) < 1E-6)
    {
        U(i,j,k,ALPHARHOLAMBDA,m) = 1E-6;

        return;
    }

    if(T<=0.0)
    {
        T = 1E-20;
    }

    Real temp = -A*U(i,j,k,ALPHARHOLAMBDA,m)*( (T > Tc) ? std::exp(-Tc/T) : 0.0);

    if(std::isnan(temp))
    {
        temp = 0.0;
    }

    U(i,j,k,ALPHARHOLAMBDA,m) +=  dt*temp;

    if(U(i,j,k,ALPHARHOLAMBDA,m) < 0.0 || std::isnan(U(i,j,k,ALPHARHOLAMBDA,m)))
    {
        U(i,j,k,ALPHARHOLAMBDA,m) = 1E-6;
    }
}


/** Adds the reactive source term with RK2.
 */
void reactiveUpdate(CellArray& U, CellArray& U1, CellArray& U2, ParameterStruct& parameters, Real dt, MultiFab& S_new)
{
    if(parameters.numberOfMixtures > 0)
    {
        U = U1;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(S_new); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            BoxAccessCellArray Ubox(mfi,bx,U);
            BoxAccessCellArray U1box(mfi,bx,U1);
            BoxAccessCellArray U2box(mfi,bx,U2);

            const auto lo = lbound(bx);
            const auto hi = ubound(bx);

            for                 (int m = 0;    m < parameters.numberOfMaterials; m++)
            {
                if(U1box.accessPattern.materialInfo[m].mixture)
                {
                    for    		(int k = lo.z; k <= hi.z; ++k)
                    {
                        for     (int j = lo.y; j <= hi.y; ++j)
                        {
                            for (int i = lo.x; i <= hi.x; ++i)
                            {
                                if(U1box(i,j,k,ALPHA,m) < 0.5 )// && U1box(i,j,k,ALPHARHOLAMBDA,m) < 0.9)
                                {
                                    //Abort("Exploding outside");
                                    U1box(i,j,k,ALPHARHOLAMBDA,m) = Ubox(i,j,k,ALPHARHOLAMBDA,m);

                                    continue;
                                }

                                pressureBased_UpdateMassFraction(Ubox,U1box,parameters,dt,i,j,k,m);

                                pressureBased_UpdateMassFraction(U1box,U2box,parameters,dt,i,j,k,m);

                                U1box(i,j,k,ALPHARHOLAMBDA,m) = ((Ubox(i,j,k,ALPHARHOLAMBDA,m)*(1.0/2.0))+(U2box(i,j,k,ALPHARHOLAMBDA,m)*(1.0/2.0)));


                            }
                        }
                    }
                }
            }

            U1box.conservativeToPrimitive();
        }
    }
}

void reactiveUpdateInHLLC(BoxAccessCellArray& U, ParameterStruct& parameters, Real dt)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    for                 (int m = 0;    m < parameters.numberOfMaterials; m++)
    {
        if(U.accessPattern.materialInfo[m].mixture)
        {
            for    		(int k = lo.z; k <= hi.z; ++k)
            {
                for     (int j = lo.y; j <= hi.y; ++j)
                {
                    for (int i = lo.x; i <= hi.x; ++i)
                    {
                        if(U(i,j,k,ALPHA,m) < 1E-3 )
                        {
                            continue;
                        }

                        IandG_UpdateMassFraction_single(U,parameters,dt,i,j,k,m);

                        //pressureBased_UpdateMassFraction_single(U,parameters,dt,i,j,k,m);

                        //Arrhenius_UpdateMassFraction_single(U,parameters,dt,i,j,k,m);

                        //Schoch_ProgrammedBurn_UpdateMassFraction_single(U,parameters,dt,i,j,k,m);
                    }
                }
            }
        }
    }
}
