#include "simulationheader.h"

/** Calculates the contact wave speed estimate.
 */
double getSstar(Cell& UL, Cell& UR, Real SL, Real SR, int i, int j, int k)
{
    return (UR(P)-UL(P)+UL(RHOU)*(SL-UL(VELOCITY))-UR(RHOU)*(SR-UR(VELOCITY)))/(UL(RHO)*(SL-UL(VELOCITY)) -  UR(RHO)*(SR-UR(VELOCITY)));
}

/** Calculates the HLLC intermediate state.
 */
void getStarState(Cell& U, Cell& UStar, Real SK, Real Sstar, ParameterStruct& parameters)
{
    static const Real tolerance  = 1E-20;

    Real multiplier = ( std::abs((Sstar-U(VELOCITY))/Sstar) < tolerance ? 1.0 : (SK-U(VELOCITY))/(SK-Sstar)   );

    UStar(RHO) = 0.0;

    for(int m=0;m<parameters.numberOfMaterials;m++)
    {
        UStar(ALPHA,m)          = U(ALPHA,m);
        UStar(ALPHARHO,m)       = multiplier*U(ALPHARHO,m);

        UStar(RHO_K,m)          = UStar(ALPHARHO,m)/UStar(ALPHA,m);

        UStar(RHO)             += UStar(ALPHARHO,m);

    }

    UStar(RHOU)         = multiplier*U(RHO)*Sstar;
    UStar(TOTAL_E)      = multiplier*U(RHO)*(U(TOTAL_E)/U(RHO) + (Sstar-U(VELOCITY))*(Sstar + (U(P))/(U(RHO)*(SK-U(VELOCITY))) ));


    UStar(VELOCITY)     = UStar(RHOU)/UStar(RHO);

    return;
}

/** Returns the equation flux of the conservative variables.
 */
Real flux(MaterialSpecifier n, Cell& U)
{
    switch(n.var)
    {
        case ALPHA:         return U(ALPHA,n.mat)*U(VELOCITY);
        case ALPHARHO:  	return U(ALPHARHO,n.mat)*U(VELOCITY);
        case RHOU: 			return U(RHOU)*U(VELOCITY)+ U(P);
        case TOTAL_E:	   	return U(VELOCITY)*(U(TOTAL_E)+U(P));
        default:   amrex::Print() << "Bad flux variable" << std::endl; exit(1);
    }

}

/** Calculates the HLLC fluxes.
 */
void calc_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, ParameterStruct& parameters)
{
    const auto lo = lbound(ULbox.box);
    const auto hi = ubound(ULbox.box);

    Real SR;
    Real SL;
    Real Sstar;

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

                fluxbox(i,j,k,USTAR) = Sstar;

                if(SL>=0.0)
                {
                    getStarState(UL,UStar,SL,Sstar,parameters);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,UStar);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL);
                        }
                    }
                }
                else if(Sstar>=0.0)
                {
                    getStarState(UL,UStar,SL,Sstar,parameters);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,UStar);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL)+SL*(UStar(n)-UL(n));
                        }
                    }

                }
                else if(SR>=0.0)
                {
                    getStarState(UR,UStar,SR,Sstar,parameters);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,UStar);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UR)+SR*(UStar(n)-UR(n));
                        }
                    }
                }
                else
                {
                    getStarState(UR,UStar,SR,Sstar,parameters);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,UStar);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UR);
                        }
                    }
                }
            }
        }
    }

    return;
}

/** Update the conservative variables using the calculated fluxes.
 *  The volume fraction alpha also needs a special
 *  non-conservative update.
 */
void update(BoxAccessCellArray& fluxbox, BoxAccessCellArray& Ubox, BoxAccessCellArray& U1box, ParameterStruct& parameters) //(Box const& box, Array4<Real> const& flux_prop_arr, Array4<Real> const& prop_arr, Array4<Real> const& prop_arr_new, ParameterStruct& parameters)
{

    const auto lo = lbound(Ubox.box);
    const auto hi = ubound(Ubox.box);

    for    			(auto n : Ubox.accessPattern.conservativeVariables)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    if(n.var == ALPHA)
                    {
                        U1box(i,j,k,n) = Ubox(i,j,k,n) + (parameters.dt/parameters.dx[0])*(fluxbox(i,j,k,n) - fluxbox(i+1,j,k,n) -Ubox(i,j,k,n)*(fluxbox(i,j,k,USTAR) - fluxbox(i+1,j,k,USTAR)));
                    }
                    else
                    {
                        U1box(i,j,k,n) = Ubox(i,j,k,n) + (parameters.dt/parameters.dx[0])*(fluxbox(i,j,k,n) - fluxbox(i+1,j,k,n));
                    }
                }
            }
        }
    }


}

/** van Leer slope limiter used in MUSCL extrapolation.
 */
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

/** Performs MUSCL extrapolation on the primitive variables.
 */
void MUSCLextrapolate(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& grad)
{
    Real r;

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    for    			(auto n : U.accessPattern.primitiveVariables)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    r = (U(i,j,k,n)-U(i-1,j,k,n))/(U(i+1,j,k,n)-U(i,j,k,n));

                    if(std::isinf(r) || std::isnan(r))
                    {
                        r = 0.0;
                    }

                    grad(i,j,k,n) = vanLeerlimiter(r)*0.5*(U(i+1,j,k,n)-U(i-1,j,k,n));

                    UL(i,j,k,n) = U(i,j,k,n) - 0.5*grad(i,j,k,n);
                    UR(i,j,k,n) = U(i,j,k,n) + 0.5*grad(i,j,k,n);

                }
            }
        }
    }

    return;
}

/** Calculate the HLLC flux and update new array
 */
void HLLCadvance(CellArray& U,CellArray& U1, CellArray& UL, CellArray& UR, CellArray& MUSCLgrad, CellArray& UStar,Array<MultiFab, AMREX_SPACEDIM>& flux_arr,Geometry const& geom, ParameterStruct& parameters,Vector<BCRec>& bc)
{
    /*-------------------------------------------------------------
     * Perform MUSCL extrapolation.
     * -----------------------------------------------------------*/

    for(MFIter mfi(U.data); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        /*-------------------------------------------------------------
         * Data can't be accessed straight from a Multifab so we make
         * some wrappers to hold the FArrayBoxes that can access the
         * data called BoxAccessCellArray.
         * -----------------------------------------------------------*/

        BoxAccessCellArray Ubox(mfi,bx,U);
        BoxAccessCellArray ULbox(mfi,bx,UL);
        BoxAccessCellArray URbox(mfi,bx,UR);
        BoxAccessCellArray gradbox(mfi,bx,MUSCLgrad);

        MUSCLextrapolate(Ubox,ULbox,URbox,gradbox);

        ULbox.primitiveToConservative(parameters);
        URbox.primitiveToConservative(parameters);

        ULbox.getSoundSpeed(parameters);
        URbox.getSoundSpeed(parameters);

    }

    FillDomainBoundary(UL.data, geom, bc);
    FillDomainBoundary(UR.data, geom, bc);

    UL.data.FillBoundary(geom.periodicity());
    UR.data.FillBoundary(geom.periodicity());

    /*-------------------------------------------------------------
     * Calulate HLLC flux and update the new array.
     * -----------------------------------------------------------*/

    for(MFIter mfi(U.data); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        FArrayBox& flux_fab_x   = flux_arr[0][mfi];

        BoxAccessCellArray Ubox(mfi,bx,U);
        BoxAccessCellArray U1box(mfi,bx,U1);
        BoxAccessCellArray ULbox(mfi,bx,UL);
        BoxAccessCellArray URbox(mfi,bx,UR);
        BoxAccessCellArray UStarbox(mfi,bx,UStar);
        BoxAccessCellArray fluxbox(bx,flux_fab_x,U);

        calc_fluxes (fluxbox, ULbox, URbox, UStarbox, parameters);
        update      (fluxbox, Ubox, U1box, parameters);

        U1box.conservativeToPrimitive(parameters);
    }

    FillDomainBoundary(U1.data, geom, bc);
    U1.data.FillBoundary(geom.periodicity());



}

/** 2nd Order Runge-Kutta time integration
 */
void advance(CellArray& U,CellArray& U1,CellArray& U2, CellArray& MUSCLgrad, CellArray& UL, CellArray& UR, CellArray& UStar,Array<MultiFab, AMREX_SPACEDIM>& flux_arr,Geometry const& geom, ParameterStruct& parameters, Vector<BCRec>& bc)
{

    HLLCadvance(U, U1, UL, UR, MUSCLgrad, UStar, flux_arr, geom, parameters, bc);

    HLLCadvance(U1, U2, UL, UR, MUSCLgrad, UStar, flux_arr, geom, parameters, bc);

    U1 = ((U*(1.0/2.0))+(U2*(1.0/2.0)));

}


