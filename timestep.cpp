#include "timestep.h"

TimeStep::TimeStep(BoxArray &ba, DistributionMapping &dm, const int Ncomp, const int Nghost) : data(ba,dm,Ncomp,Nghost){}

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

Real TimeStep::getTimeStep(CellArray& U, ParameterStruct& parameters)
{
    for(MFIter mfi(U.data); mfi.isValid(); ++mfi)
    {
        const Box& bx   = mfi.validbox();

        FArrayBox& fab  = U.data[mfi];
        FArrayBox& time_fab = data[mfi];

        Array4<Real> const& prop_arr_phi  = fab.array();
        Array4<Real> const& prop_arr_time = time_fab.array();


        BoxAccessCellArray baca(bx,fab,prop_arr_phi,U);

        U.getSoundSpeed(baca,parameters);

        boxTimeStepFunction(baca,prop_arr_time,parameters);
    }

    Real dt = std::numeric_limits<Real>::max();

    for(int n=0;n<AMREX_SPACEDIM;n++)
    {
        dt = (data.min(n) < dt ? data.min(n) : dt);
    }


    return (parameters.CFL*dt);
}

