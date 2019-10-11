#include "timestep.h"

TimeStep::TimeStep(BoxArray &ba, DistributionMapping &dm, const int Ncomp, const int Nghost) : data(ba,dm,Ncomp,Nghost){}

/** Calulates the Timestep in a given box.
 */
void TimeStep::boxTimeStepFunction(BoxAccessCellArray& U, Array4<Real> const& prop_arr_time, ParameterStruct& parameters)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);


    for 		(int k = lo.z; k <= hi.z; ++k)
    {
        for 	(int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                for(int n=0;n<AMREX_SPACEDIM;n++)
                {
                    prop_arr_time(i,j,k,n) = parameters.dx[n]/(U(i,j,k,SOUNDSPEED) + fabs(U(i,j,k,VELOCITY,0,n)));
                }

            }
        }
    }

    return;
}

/** Calulates the Timestep over the entire domain.
 * Calls boxTimestepFunction to calculate the steps
 * in each box before finding the minimumover the whole
 * multifab.
 */
Real TimeStep::getTimeStep(CellArray& U, ParameterStruct& parameters)
{
    for(MFIter mfi(U.data); mfi.isValid(); ++mfi)
    {
        const Box& bx   = mfi.validbox();

        FArrayBox& time_fab = data[mfi];

        Array4<Real> const& prop_arr_time = time_fab.array();

        BoxAccessCellArray baca(mfi,bx,U);

        baca.getSoundSpeed(parameters);

        boxTimeStepFunction(baca,prop_arr_time,parameters);
    }

    Real dt = std::numeric_limits<Real>::max();

    for(int n=0;n<AMREX_SPACEDIM;n++)
    {
        Real temp = data.min(n);

        dt = (temp < dt ? temp : dt);
    }


    return (parameters.CFL*dt);
}

