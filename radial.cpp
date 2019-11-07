#include "simulationheader.h"

/** Adds the Geometric source term using the Geometric flux function to a box.
 */
void addGeometricSourceTerm(BoxAccessCellArray& U, ParameterStruct& parameters, const Real* dx, Real dt, const Real* prob_lo)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    Real inverse_r;

    const static Real multiplier = (AMREX_SPACEDIM == 1 ? -2.0 : -1.0) ;

    Real x,y,z;

    for    			(auto n : U.accessPattern.conservativeVariables)
    {  
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
                    z = prob_lo[2] + (Real(k)+0.5)*dx[2];

            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                    y = prob_lo[1] + (Real(j)+0.5)*dx[1];

                for (int i = lo.x; i <= hi.x; ++i)
                {
                    x = prob_lo[0] + (Real(i)+0.5)*dx[0];

                    inverse_r = multiplier/(x);

                    U(i,j,k,n) += dt*inverse_r*geometricFlux(U,i,j,k,n);

                }
            }
        }
    }

    return;
}


/** Creates boxes to add the Geometric source term using the Geometric flux function.
 */
void geometricSourceTerm(CellArray& U, ParameterStruct& parameters, const Real* dx, Real dt, const Real* prob_lo)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(U.data); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray  Ubox(mfi,bx,U);

        addGeometricSourceTerm(Ubox,parameters,dx,dt,prob_lo);

        Ubox.conservativeToPrimitive();
    }

    return;
}
