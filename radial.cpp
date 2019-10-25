#include "simulationheader.h"

/** Adds the Geometric source term using the Geometric flux function to a box.
 */
void addGeometricSourceTerm(BoxAccessCellArray& U, ParameterStruct& parameters)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    Real inverse_r;

    const static Real multiplier = (AMREX_SPACEDIM == 1 ? -2.0 : -1.0) ;


    for    			(auto n : U.accessPattern.conservativeVariables)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    inverse_r = multiplier/(((Real)i+0.5)*parameters.dx[0]);

                    U(i,j,k,n) += parameters.dt*inverse_r*geometricFlux(U,i,j,k,n);
                }
            }
        }
    }

    return;
}


/** Creates boxes to add the Geometric source term using the Geometric flux function.
 */
void geometricSourceTerm(CellArray& U, ParameterStruct& parameters)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(U.data); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray  Ubox(mfi,bx,U);

        addGeometricSourceTerm(Ubox,parameters);

        Ubox.conservativeToPrimitive();
    }
    return;
}
