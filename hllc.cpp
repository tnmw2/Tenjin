#include "simulationheader.h"


double getSstar(Cell& UL, Cell& UR, Real SL, Real SR, int i, int j, int k)
{
    return (UR(P)-UL(P)+UL(RHOU)*(SL-UL(VELOCITY))-UR(RHOU)*(SR-UR(VELOCITY)))/(UL(RHO)*(SL-UL(VELOCITY)) -  UR(RHO)*(SR-UR(VELOCITY)));
}

void getStarState(Cell& U, Cell& UStar, Real SK, Real Sstar, ParameterStruct& parameters)
{
    static const Real tolerance  = 1E-20;

    Real multiplier = ( std::abs((Sstar-U(VELOCITY))/Sstar) < tolerance ? 1.0 : (SK-U(VELOCITY))/(SK-Sstar)   );

    UStar(RHO)          = multiplier*U(RHO);
    UStar(RHOU)         = multiplier*U(RHO)*Sstar;
    UStar(TOTAL_E)      = multiplier*U(RHO)*(U(TOTAL_E)/U(RHO) + (Sstar-U(VELOCITY))*(Sstar + (U(P))/(U(RHO)*(SK-U(VELOCITY))) ));


    UStar(VELOCITY)     = UStar(RHOU)/UStar(RHO);
    UStar(P)            = (UStar(TOTAL_E)-0.5*UStar(RHO)*UStar(VELOCITY)*UStar(VELOCITY))*(parameters.adiabaticIndex-1.0);

    return;
}

Real flux(MaterialSpecifier n, Cell& U)
{
    switch(n.var)
    {
        case RHO:  			return U(RHOU);
        case RHOU: 			return U(RHOU,n.mat,n.row)*U(VELOCITY,n.mat,n.row)+ U(P);
        case TOTAL_E:	   	return U(VELOCITY,n.mat,n.row)*(U(TOTAL_E)+U(P));
        default:   amrex::Print() << "Bad flux variable" << std::endl; exit(1);
    }

}

void calc_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, ParameterStruct& parameters)
{
    const auto lo = lbound(ULbox.box);
    const auto hi = ubound(ULbox.box);

    Real SR;
    Real SL;
    Real Sstar;

    std::vector<MaterialSpecifier> varsThatNeedUpdating = {MaterialSpecifier(RHO),MaterialSpecifier(RHOU),MaterialSpecifier(TOTAL_E)};

    for 		(int k = lo.z; k <= hi.z; ++k)
    {
        for 	(int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x +1; ++i)
            {

                Cell UL(URbox,i-1,j,k);
                Cell UR(ULbox,i,j,k);
                Cell UStar(UStarbox,i,j,k);


                SR = std::max(std::abs(UL(VELOCITY))+UL(SOUNDSPEED),std::abs(UR(VELOCITY))+UR(SOUNDSPEED));
                SL = -SR;


                Sstar = getSstar(UL,UR,SL,SR,i,j,k);

                if(SL>=0.0)
                {
                    for(auto n : varsThatNeedUpdating)
                    {
                        fluxbox(n,i,j,k)	= flux(n,UL);
                    }
                }
                else if(Sstar>=0.0)
                {
                    getStarState(UL,UStar,SL,Sstar,parameters);

                    for(auto n : varsThatNeedUpdating)
                    {
                        fluxbox(n,i,j,k)	= flux(n,UL)+SL*(UStarbox(n,i,j,k)-URbox(n,i-1,j,k));
                    }
                }
                else if(SR>=0.0)
                {
                    getStarState(UR,UStar,SR,Sstar,parameters);

                    for(auto n : varsThatNeedUpdating)
                    {
                        fluxbox(n,i,j,k)	= flux(n,UR)+SR*(UStarbox(n,i,j,k)-ULbox(n,i-1,j,k));
                    }
                }
                else
                {
                    for(auto n : varsThatNeedUpdating)
                    {
                        fluxbox(n,i,j,k)	= flux(n,UR);
                    }
                }
            }
        }
    }

    return;
}

/*void calc_fluxes(Box const& box, Array4<Real> const& flux_prop_arr, Array4<Real> const& prop_arr_L, Array4<Real> const& prop_arr_R, Array4<Real> const& prop_arr_Star, ParameterStruct& parameters)
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

}*/

void update(BoxAccessCellArray& fluxbox, BoxAccessCellArray& Ubox, BoxAccessCellArray& U1box, ParameterStruct& parameters) //(Box const& box, Array4<Real> const& flux_prop_arr, Array4<Real> const& prop_arr, Array4<Real> const& prop_arr_new, ParameterStruct& parameters)
{

    const auto lo = lbound(Ubox.box);
    const auto hi = ubound(Ubox.box);

    std::vector<MaterialSpecifier> varsThatNeedUpdating = {MaterialSpecifier(RHO),MaterialSpecifier(RHOU),MaterialSpecifier(TOTAL_E)};

    for    			(auto n : varsThatNeedUpdating)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    U1box(n,i,j,k) = Ubox(n,i,j,k) + (parameters.dt/parameters.dx[0])*(fluxbox(n,i,j,k) - fluxbox(n,i+1,j,k));
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

    std::vector<MaterialSpecifier> varsThatNeedUpdating = {MaterialSpecifier(RHO),MaterialSpecifier(VELOCITY),MaterialSpecifier(P)};

    for    			(auto n : varsThatNeedUpdating)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    r = (U(n,i,j,k)-U(n,i-1,j,k))/(U(n,i+1,j,k)-U(n,i,j,k));

                    if(std::isinf(r) || std::isnan(r))
                    {
                        r = 0.0;
                    }

                    grad(n,i,j,k) = vanLeerlimiter(r)*0.5*(U(n,i+1,j,k)-U(n,i-1,j,k));

                    UL(n,i,j,k) = U(n,i,j,k) - 0.5*grad(n,i,j,k);
                    UR(n,i,j,k) = U(n,i,j,k) + 0.5*grad(n,i,j,k);

                    //UL(n,i,j,k) = U(n,i,j,k);
                    //UR(n,i,j,k) = U(n,i,j,k);


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

        BoxAccessCellArray Ubox(bx,fab_old,prop_arr,U);
        BoxAccessCellArray ULbox(bx,fab_L,prop_arr_L,UL);
        BoxAccessCellArray URbox(bx,fab_R,prop_arr_R,UR);
        BoxAccessCellArray gradbox(bx,fab_grad,prop_arr_grad,MUSCLgrad);

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

        BoxAccessCellArray Ubox(bx,fab_old,prop_arr,U);
        BoxAccessCellArray U1box(bx,fab_new,prop_arr_new,U1);
        BoxAccessCellArray ULbox(bx,fab_L,prop_arr_L,UL);
        BoxAccessCellArray URbox(bx,fab_R,prop_arr_R,UR);
        BoxAccessCellArray UStarbox(bx,fab_Star,prop_arr_Star,UStar);
        BoxAccessCellArray fluxbox(bx,flux_fab_x,flux_prop_arr_x,U);

        //calc_fluxes (fluxbox, ULbox, URbox, UStarbox, parameters);
        calc_fluxes (fluxbox, ULbox, URbox, UStarbox, parameters);
        update      (fluxbox, Ubox, U1box, parameters);


        //calc_fluxes	(bx, flux_prop_arr_x, prop_arr_L, prop_arr_R, prop_arr_Star, parameters);
        //update		(bx, flux_prop_arr_x, prop_arr, prop_arr_new, parameters);

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


