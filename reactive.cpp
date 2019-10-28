#include "simulationheader.h"

void pressureBased_UpdateMassFraction(BoxAccessCellArray& U, BoxAccessCellArray& U1, ParameterStruct& parameters, int i, int j, int k, int m)
{

    Real sigma = 20.0   ;//10.0;
    Real nu    = 0.5;
    Real n     = 1.2;

    if(U(i,j,k,ALPHA,m) < 0.5)
    {
        U1(i,j,k,ALPHARHOLAMBDA,m) = U(i,j,k,ALPHARHOLAMBDA,m);
        return;
    }


    U1(i,j,k,ALPHARHOLAMBDA,m) = U(i,j,k,ALPHARHOLAMBDA,m)-(parameters.dt)*(sigma*std::pow(U(i,j,k,ALPHARHO,m),1.0-nu)*std::pow(U(i,j,k,ALPHARHOLAMBDA,m),nu)*std::pow((U(i,j,k,P)/1E6),n));

    if(U1(i,j,k,ALPHARHOLAMBDA,m) < 0.0)
    {
        U1(i,j,k,ALPHARHOLAMBDA,m) = 0.0;
    }

}

void Arrhenius_UpdateMassFraction(BoxAccessCellArray& U, BoxAccessCellArray& U1, ParameterStruct& parameters, int i, int j, int k, int m)
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

    U1(i,j,k,ALPHARHOLAMBDA,m) = U(i,j,k,ALPHARHOLAMBDA,m)+(parameters.dt)*temp;

    if(U1(i,j,k,ALPHARHOLAMBDA,m) < 0.0)
    {
        U1(i,j,k,ALPHARHOLAMBDA,m) = 0.0;
    }


}


/** Adds the reactive source term with RK2.
 */
void reactiveUpdate(CellArray& U, CellArray& U1, CellArray& U2, ParameterStruct& parameters)
{
    if(parameters.numberOfMixtures > 0)
    {
        U = U1;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(U1.data); mfi.isValid(); ++mfi )
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
                                pressureBased_UpdateMassFraction(Ubox,U1box,parameters,i,j,k,m);
                                pressureBased_UpdateMassFraction(U1box,U2box,parameters,i,j,k,m);

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
