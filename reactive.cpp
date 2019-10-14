#include "simulationheader.h"

void updateMassFraction(BoxAccessCellArray& U, BoxAccessCellArray& U1, ParameterStruct& parameters)
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
                        T = U(i,j,k,P)/(U(i,j,k,RHO_MIX,m,0)*(parameters.adiabaticIndex[m]-1.0)*parameters.CV[m]);

                        if(T<=0.0)
                        {
                            T = 1E-20;
                        }

                        temp = -A*U(i,j,k,ALPHARHOLAMBDA,m)*std::exp(-Tc/T);

                        if(std::isnan(temp))
                        {
                            temp = 0.0;
                        }

                        //Print() << temp << std::endl;

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

        updateMassFraction(Ubox,U1box,parameters);

        U1box.conservativeToPrimitive();

    }
}


/** Adds the reactive source term with RK2.
 */
void reactiveUpdate(CellArray& U, CellArray& U1, CellArray& U2, ParameterStruct& parameters)
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
