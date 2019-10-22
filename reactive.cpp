#include "simulationheader.h"

void ignition_updateMassFraction(BoxAccessCellArray& U, BoxAccessCellArray& U1, ParameterStruct& parameters)
{

    Real K      = 2.5E5;
    Real p_ig   = 0.01E6;

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

                        U1(i,j,k,ALPHARHOLAMBDA,m) = U(i,j,k,ALPHARHOLAMBDA,m)-(parameters.dt)*(K*U(i,j,k,ALPHARHO,m)*sqrt(U(i,j,k,LAMBDA,m))*Heaviside<Real,Real>(U(i,j,k,P)-p_ig));

                        if(U1(i,j,k,ALPHARHOLAMBDA,m) < 0.0)
                        {
                            U1(i,j,k,ALPHARHOLAMBDA,m) = 0.0;
                        }
                    }
                }
            }
        }
    }
}

void pressureBased_updateMassFraction(BoxAccessCellArray& U, BoxAccessCellArray& U1, ParameterStruct& parameters)
{

    Real sigma = 20.0   ;//10.0;
    Real nu    = 0.5;
    Real n     = 1.2;

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
                        U1(i,j,k,ALPHARHOLAMBDA,m) = U(i,j,k,ALPHARHOLAMBDA,m)-(parameters.dt)*(sigma*std::pow(U(i,j,k,ALPHARHO,m),1.0-nu)*std::pow(U(i,j,k,ALPHARHOLAMBDA,m),nu)*std::pow((U(i,j,k,P)/1E6),n));
                        //U1(i,j,k,ALPHARHOLAMBDA,m) = U(i,j,k,ALPHARHOLAMBDA,m)-(parameters.dt)*(sigma*U(i,j,k,ALPHARHO,m)*std::pow(U(i,j,k,LAMBDA,m),nu)*std::pow(U(i,j,k,P),n));

                        if(U1(i,j,k,ALPHARHOLAMBDA,m) < 0.0)
                        {
                            U1(i,j,k,ALPHARHOLAMBDA,m) = 0.0;
                        }
                    }
                }
            }
        }
    }
}

void ArrheniusRateLaw_updateMassFraction(BoxAccessCellArray& U, BoxAccessCellArray& U1, ParameterStruct& parameters)
{

    Real Tc = 150.0; //Reference temp;

    Real A = 2E6;    //Pre-exponential factor

    Real T;

    Real temp;

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
                        T = U.accessPattern.materialInfo[m].EOS->getTemp(U,i,j,k,m,0);

                        if(T<=0.0)
                        {
                            T = 1E-20;
                        }

                        temp = -A*U(i,j,k,ALPHARHOLAMBDA,m)*std::exp(-Tc/T);

                        if(std::isnan(temp))
                        {
                            temp = 0.0;
                        }

                        U1(i,j,k,ALPHARHOLAMBDA,m) = U(i,j,k,ALPHARHOLAMBDA,m)+(parameters.dt)*temp;

                    }
                }
            }
        }
    }
}

/** Adds the reactive source term.
 */
void RKreactiveUpdate(CellArray& U, CellArray& U1, ParameterStruct& parameters)
{
    U1 = U;

    for(MFIter mfi(U.data); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray Ubox(mfi,bx,U);
        BoxAccessCellArray U1box(mfi,bx,U1);

        pressureBased_updateMassFraction(Ubox,U1box,parameters);

        U1box.conservativeToPrimitive();

    }
}

/** Adds the reactive source term with RK2.
 */
void reactiveUpdate(CellArray& U, CellArray& U1, CellArray& U2, ParameterStruct& parameters)
{
    if(parameters.numberOfMixtures > 0)
    {
        U = U1;

        RKreactiveUpdate(U,U1,parameters);
        RKreactiveUpdate(U1,U2,parameters);

        //U1 = ((U1*(1.0/2.0))+(U2*(1.0/2.0)));

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
                                U1box(i,j,k,ALPHARHO,m) = ((Ubox(i,j,k,ALPHARHO,m)*(1.0/2.0))+(U2box(i,j,k,ALPHARHO,m)*(1.0/2.0)));
                            }
                        }
                    }
                }
            }

            U1box.conservativeToPrimitive();
        }
    }
}
