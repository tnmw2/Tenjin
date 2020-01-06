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

    Real sigma = 20.0;
    Real nu    = 0.5;
    Real n     = 1.5;


    if(U(i,j,k,ALPHARHOLAMBDA,m) < 0.0)
    {
        U(i,j,k,ALPHARHOLAMBDA,m) = 0.0;
    }

    U(i,j,k,ALPHARHOLAMBDA,m) -=  dt*(sigma*std::pow(U(i,j,k,ALPHARHO,m),1.0-nu)*std::pow(U(i,j,k,ALPHARHOLAMBDA,m),nu)*std::pow( ((U(i,j,k,P) > 1E6 ? U(i,j,k,P) : 0.0)/1E6),n));

    if(U(i,j,k,ALPHARHOLAMBDA,m) < 0.0 || std::isnan(U(i,j,k,ALPHARHOLAMBDA,m)))
    {
        U(i,j,k,ALPHARHOLAMBDA,m) = 0.0;
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
                        /*if(U(i,j,k,ALPHA,m) < 0.5 )
                        {
                            continue;
                        }*/

                        pressureBased_UpdateMassFraction_single(U,parameters,dt,i,j,k,m);
                    }
                }
            }
        }
    }
}
