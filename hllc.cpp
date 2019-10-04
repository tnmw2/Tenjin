#include "simulationheader.h"

class Cell
{
    public:

    Real& rho;
    Real& rhoU;
    Real& E;
    Real& u;
    Real& p;
    Real& a;

    Cell(Array4<Real> const& arr, int i, int j, int k) : rho(arr(i,j,k,RHO)), rhoU(arr(i,j,k,RHOU)), E(arr(i,j,k,TOTAL_E)), u(arr(i,j,k,VELOCITY)), p(arr(i,j,k,P)), a(arr(i,j,k,SOUNDSPEED)){}
};

double getSstar(Cell& UL, Cell& UR, Real SL, Real SR, int i, int j, int k)
{
    return (UR.p-UL.p+UL.rhoU*(SL-UL.u)-UR.rhoU*(SR-UR.u))/(UL.rho*(SL-UL.u) -  UR.rho*(SR-UR.u));
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

void calc_fluxes(Box const& box, Array4<Real> const& flux_prop_arr, Array4<Real> const& prop_arr_L, Array4<Real> const& prop_arr_R, Array4<Real> const& prop_arr_Star, ParameterStruct& parameters)
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

                Cell UL(prop_arr_R,i-1,j,k);
                Cell UR(prop_arr_L,i,j,k);
                Cell UStar(prop_arr_Star,i,j,k);


                SR = std::max(std::abs(UL.u)+UL.a,std::abs(UR.u)+UR.a);
                SL = -SR;


                Sstar = getSstar(UL,UR,SL,SR,i,j,k);

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
                        flux_prop_arr(i,j,k,n)	= flux(n,UL)+SL*(prop_arr_Star(i,j,k,n)-prop_arr_R(i-1,j,k,n));
                    }
                }
                else if(SR>=0.0)
                {
                    getStarState(UR,UStar,SR,Sstar,parameters);

                    for(int n = 0 ; n <= TOTAL_E; ++n)
                    {
                        flux_prop_arr(i,j,k,n)	= flux(n,UR)+SR*(prop_arr_Star(i,j,k,n)-prop_arr_L(i,j,k,n));
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

Real vanLeerlimiter(Real r)
{
    if(r<=0.0)
    {
        return 0.0;
    }
    else
    {
        return std::min(2.0/(1.0+r),2.0*r/(1.0+r));
    }
}

void MUSCLextrapolate(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& grad)
{
    Real r;

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    for    			(int n = 0   ; n < U.arr.nComp(); ++n)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    r = (U.arr(i,j,k,n)-U.arr(i-1,j,k,n))/(U.arr(i+1,j,k,n)-U.arr(i,j,k,n));

                    if(std::isinf(r) || std::isnan(r))
                    {
                        r = 0.0;
                    }

                    grad.arr(i,j,k,n) = vanLeerlimiter(r)*0.5*(U.arr(i+1,j,k,n)-U.arr(i-1,j,k,n));

                    UL.arr(i,j,k,n) = U.arr(i,j,k,n) - 0.5*grad.arr(i,j,k,n);
                    UR.arr(i,j,k,n) = U.arr(i,j,k,n) + 0.5*grad.arr(i,j,k,n);

                    //UL.arr(i,j,k,n) = U.arr(i,j,k,n);
                    //UR.arr(i,j,k,n) = U.arr(i,j,k,n);


                }
            }
        }
    }

    return;
}

void HLLCadvance(CellArray& U,CellArray& U1, CellArray& UL, CellArray& UR, CellArray& MUSCLgrad, CellArray& UStar,Array<MultiFab, AMREX_SPACEDIM>& flux_arr,Geometry const& geom, ParameterStruct& parameters,Vector<BCRec>& bc)
{
    for(MFIter mfi(U.data); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        //access current array box
        FArrayBox& fab_old      = U.data[mfi];
        FArrayBox& fab_L        = UL.data[mfi];
        FArrayBox& fab_R        = UR.data[mfi];
        FArrayBox& fab_grad     = MUSCLgrad.data[mfi];


        //access the array from the array box:
        Array4<Real> const& prop_arr 		= fab_old.array();
        Array4<Real> const& prop_arr_L      = fab_L.array();
        Array4<Real> const& prop_arr_R      = fab_R.array();
        Array4<Real> const& prop_arr_grad	= fab_grad.array();

        BoxAccessCellArray Ubox(bx,fab_old,prop_arr);
        BoxAccessCellArray ULbox(bx,fab_L,prop_arr_L);
        BoxAccessCellArray URbox(bx,fab_R,prop_arr_R);
        BoxAccessCellArray gradbox(bx,fab_grad,prop_arr_grad);

        MUSCLextrapolate(Ubox,ULbox,URbox,gradbox);

        //UL.conservativeToPrimitive(ULbox,parameters);
        //UR.conservativeToPrimitive(URbox,parameters);

        UL.primitiveToConservative(ULbox,parameters);
        UR.primitiveToConservative(URbox,parameters);


        UL.getSoundSpeed(ULbox,parameters);
        UR.getSoundSpeed(URbox,parameters);
    }

    FillDomainBoundary(UL.data, geom, bc);
    FillDomainBoundary(UR.data, geom, bc);

    UL.data.FillBoundary(geom.periodicity());
    UR.data.FillBoundary(geom.periodicity());

    for(MFIter mfi(U.data); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        //access current array box
        FArrayBox& fab_old      = U.data[mfi];
        FArrayBox& fab_new      = U1.data[mfi];
        FArrayBox& fab_L        = UL.data[mfi];
        FArrayBox& fab_R        = UR.data[mfi];
        FArrayBox& fab_Star     = UStar.data[mfi];
        FArrayBox& flux_fab_x   = flux_arr[0][mfi];

        //access the array from the array box:
        Array4<Real> const& prop_arr 		= fab_old.array();
        Array4<Real> const& prop_arr_new	= fab_new.array();
        Array4<Real> const& prop_arr_L      = fab_L.array();
        Array4<Real> const& prop_arr_R      = fab_R.array();
        Array4<Real> const& prop_arr_Star	= fab_Star.array();
        Array4<Real> const& flux_prop_arr_x = flux_fab_x.array();

        BoxAccessCellArray Ubox(bx,fab_old,prop_arr);
        BoxAccessCellArray U1box(bx,fab_new,prop_arr_new);
        BoxAccessCellArray ULbox(bx,fab_L,prop_arr_L);
        BoxAccessCellArray URbox(bx,fab_R,prop_arr_R);

        calc_fluxes	(bx, flux_prop_arr_x, prop_arr_L, prop_arr_R, prop_arr_Star, parameters);
        update		(bx, flux_prop_arr_x, prop_arr, prop_arr_new, parameters);

        U1.conservativeToPrimitive(U1box,parameters);
    }

    FillDomainBoundary(U1.data, geom, bc);
    U1.data.FillBoundary(geom.periodicity());

}

void RKadvance(CellArray& U,CellArray& U1,CellArray& U2, CellArray& MUSCLgrad, CellArray& UL, CellArray& UR, CellArray& UStar,Array<MultiFab, AMREX_SPACEDIM>& flux_arr,Geometry const& geom, ParameterStruct& parameters, Vector<BCRec>& bc)
{

    HLLCadvance(U, U1, UL, UR, MUSCLgrad, UStar, flux_arr, geom, parameters, bc);
    HLLCadvance(U1, U2, UL, UR, MUSCLgrad, UStar, flux_arr, geom, parameters, bc);

    U1 = ((U*(1.0/2.0))+(U2*(1.0/2.0)));

    //HLLCadvance(U, U1, UL, UR, MUSCLgrad, UStar, flux_arr, geom, parameters, bc);

}


