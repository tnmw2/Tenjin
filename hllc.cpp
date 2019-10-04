#include "simulationheader.h"

class Cell
{
    public:

    Real& rho;
    Real& rhoU;
    Real& E;
    Real& u;
    Real& p;

    Cell(Array4<Real> const& arr, int i, int j, int k) : rho(arr(i,j,k,RHO)), rhoU(arr(i,j,k,RHOU)), E(arr(i,j,k,TOTAL_E)), u(arr(i,j,k,VELOCITY)), p(arr(i,j,k,P)){}
};

double getSstar(Array4<Real> const& prop_arr, Real SL, Real SR, int i, int j, int k)
{
    return (prop_arr(i,j,k,P)-prop_arr(i-1,j,k,P)+prop_arr(i-1,j,k,RHOU)*(SL-prop_arr(i-1,j,k,VELOCITY))-prop_arr(i,j,k,RHOU)*(SR-prop_arr(i,j,k,VELOCITY)))/(prop_arr(i-1,j,k,RHO)*(SL-prop_arr(i-1,j,k,VELOCITY)) -  prop_arr(i,j,k,RHO)*(SR-prop_arr(i,j,k,VELOCITY)));
}

void getStarState(Cell& U, Cell& UStar, Real SK, Real Sstar, ParameterStruct& parameters)
{
    static const Real tolerance  = 1E-20;

    Real multiplier = ( std::abs((Sstar-U.u)/Sstar) < tolerance ? 1.0 : (SK-U.u)/(SK-Sstar)   );

    UStar.rho  = multiplier*U.rho;
    UStar.rhoU = multiplier*U.rho*Sstar;
    UStar.E    = multiplier*U.rho*(U.E/U.rho + (Sstar-U.u)*(Sstar + (U.p)/(U.rho*(SK-U.u)) ));


    UStar.u    = UStar.rhoU/UStar.rho;
    UStar.p	   = (UStar.E-0.5*UStar.rho*UStar.u*UStar.u)*(parameters.adiabaticIndex-1.0);

    return;
}

Real flux(int n, Cell& U)
{
    switch(n)
    {
        case RHO:  			return U.rhoU;
        case RHOU: 			return U.rhoU*U.u+U.p;
        case TOTAL_E:	   	return U.u*(U.E+U.p);
        default:   amrex::Print() << "Bad flux variable" << std::endl; exit(1);
    }

}

void calc_fluxes(Box const& box, Array4<Real> const& flux_prop_arr, Array4<Real> const& prop_arr, Array4<Real> const& prop_arr_Star, ParameterStruct& parameters)
{

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Real SR;
    Real SL;
    Real Sstar;


    for 		(int k = lo.z; k <= hi.z; ++k)
    {
        for 	(int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x +1; ++i)
            {
                Cell UL(prop_arr,i-1,j,k);
                Cell UR(prop_arr,i,j,k);
                Cell UStar(prop_arr_Star,i,j,k);

                SR = std::max(std::abs(prop_arr(i-1,j,k,VELOCITY))+prop_arr(i-1,j,k,SOUNDSPEED),std::abs(prop_arr(i,j,k,VELOCITY))+prop_arr(i,j,k,SOUNDSPEED));
                SL = -SR;


                Sstar = getSstar(prop_arr,SL,SR,i,j,k);

                //amrex::Print() << SL << " " << Sstar << " " << SR << std::endl;


                if(SL>=0.0)
                {
                    for(int n = 0 ; n <= TOTAL_E; ++n)
                    {
                        flux_prop_arr(i,j,k,n)	= flux(n,UL);
                    }
                }
                else if(Sstar>=0.0)
                {
                    getStarState(UL,UStar,SL,Sstar,parameters);

                    for(int n = 0 ; n <= TOTAL_E; ++n)
                    {
                        flux_prop_arr(i,j,k,n)	= flux(n,UL)+SL*(prop_arr_Star(i,j,k,n)-prop_arr(i-1,j,k,n));
                    }
                }
                else if(SR>=0.0)
                {
                    getStarState(UR,UStar,SR,Sstar,parameters);

                    for(int n = 0 ; n <= TOTAL_E; ++n)
                    {
                        flux_prop_arr(i,j,k,n)	= flux(n,UR)+SR*(prop_arr_Star(i,j,k,n)-prop_arr(i,j,k,n));
                    }
                }
                else
                {
                    for(int n = 0 ; n <= TOTAL_E; ++n)
                    {
                        flux_prop_arr(i,j,k,n)	= flux(n,UR);
                    }
                }
            }
        }
    }

    return;

}

void update(Box const& box, Array4<Real> const& flux_prop_arr, Array4<Real> const& prop_arr, Array4<Real> const& prop_arr_new, ParameterStruct& parameters)
{

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for    			(int n = 0   ; n <= TOTAL_E; ++n)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    prop_arr_new(i,j,k,n) = prop_arr(i,j,k,n) + (parameters.dt/parameters.dx[0])*(flux_prop_arr(i,j,k,n) - flux_prop_arr(i+1,j,k,n));

                }
            }
        }
    }


}

void HLLCadvance(CellArray& U,CellArray& U1,CellArray& UStar,Array<MultiFab, AMREX_SPACEDIM>& flux_arr,Geometry const& geom, ParameterStruct& parameters)
{

    for(MFIter mfi(U.data); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        //access current array box
        FArrayBox& fab_old = U.data[mfi];
        FArrayBox& fab_new = U1.data[mfi];
        FArrayBox& fab_Star= UStar.data[mfi];
        FArrayBox& flux_fab_x = flux_arr[0][mfi];

        //access the array from the array box:
        Array4<Real> const& prop_arr 		= fab_old.array();
        Array4<Real> const& prop_arr_new	= fab_new.array();
        Array4<Real> const& prop_arr_Star	= fab_Star.array();
        Array4<Real> const& flux_prop_arr_x = flux_fab_x.array();

        BoxAccessCellArray baca(bx,fab_new,prop_arr_new);

        calc_fluxes	(bx, flux_prop_arr_x, prop_arr, prop_arr_Star, parameters);
        update		(bx, flux_prop_arr_x, prop_arr, prop_arr_new, parameters);

        U1.conservativeToPrimitive(baca,parameters);
    }

}

