#include "simulationheader.h"

/** Calculates the contact wave speed estimate.
 */
double getSstar(Cell& UL, Cell& UR, Real SL, Real SR, int i, int j, int k, Direction_enum d)
{
    return (UR(P)-UL(P)+UL(RHOU,0,d)*(SL-UL(VELOCITY,0,d))-UR(RHOU,0,d)*(SR-UR(VELOCITY,0,d)))/(UL(RHO)*(SL-UL(VELOCITY,0,d)) -  UR(RHO)*(SR-UR(VELOCITY,0,d)));
}

/** Calculates the HLLC intermediate state.
 */
void getStarState(Cell& U, Cell& UStar, Real SK, Real Sstar, ParameterStruct& parameters, Direction_enum d)
{
    static const Real tolerance  = 1E-20;

    Real multiplier = ( std::abs((Sstar-U(VELOCITY,0,d))/Sstar) < tolerance ? 1.0 : (SK-U(VELOCITY,0,d))/(SK-Sstar)   );

    UStar(RHO) = 0.0;

    for(int m=0;m<parameters.numberOfMaterials;m++)
    {
        UStar(ALPHA,m)          = U(ALPHA,m);
        UStar(ALPHARHO,m)       = multiplier*U(ALPHARHO,m);

        UStar(RHO_K,m)          = UStar(ALPHARHO,m)/UStar(ALPHA,m);

        UStar(RHO)             += UStar(ALPHARHO,m);

    }

    for(int row = 0; row < AMREX_SPACEDIM; row++)
    {
        if(row == d)
        {
            UStar(RHOU,0,row)         = multiplier*U(RHO)*Sstar;
        }
        else
        {
            UStar(RHOU,0,row)         = multiplier*U(RHOU,0,row);
        }

        UStar(VELOCITY,0,row)     = UStar(RHOU,0,row)/UStar(RHO);
    }

    UStar(TOTAL_E)      = multiplier*U(RHO)*(U(TOTAL_E)/U(RHO) + (Sstar-U(VELOCITY,0,d))*(Sstar + (U(P))/(U(RHO)*(SK-U(VELOCITY,0,d))) ));

    return;
}

Real delta(int i, int j)
{
    return (Real)(i==j);
}

/** Returns the equation flux of the conservative variables.
 */
Real flux(MaterialSpecifier n, Cell& U, Direction_enum d)
{
    switch(n.var)
    {
        case ALPHA:         return U(ALPHA,n.mat)*      U(VELOCITY,0,d);
        case ALPHARHO:  	return U(ALPHARHO,n.mat)*   U(VELOCITY,0,d);
        case RHOU: 			return U(RHOU,0,n.row)*     U(VELOCITY,0,d)+U(P)*delta(n.row,(int)d);
        case TOTAL_E:	   	return U(VELOCITY,0,d)*(U(TOTAL_E)+U(P));
        default:   amrex::Print() << "Bad flux variable" << std::endl; exit(1);
    }

}

/** Calculates the HLLC fluxes.
 */
void calc_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, ParameterStruct& parameters, Direction_enum d)
{
    const auto lo = lbound(ULbox.box);
    const auto hi = ubound(ULbox.box);

    Real SR;
    Real SL;
    Real Sstar;

    IntVect extra(AMREX_D_DECL(0,0,0));

    extra[d]=1;

    for 		(int k = lo.z; k <= hi.z+extra[z]; ++k)
    {
        for 	(int j = lo.y; j <= hi.y+extra[y]; ++j)
        {
            for (int i = lo.x; i <= hi.x+extra[x]; ++i)
            {
                Cell UL;

                if(d==x)
                {
                    UL = Cell(URbox,i-1,j,k);
                }
                else if(d==y)
                {
                    UL = Cell(URbox,i,j-1,k);
                }
                else if(d==z)
                {
                    UL = Cell(URbox,i,j,k-1);
                }
                else
                {
                    Print() << "Incorrect direction in calc_fluxes" << std::endl; exit(1);
                }

                //Cell UL(URbox,i-1,j,k);
                Cell UR(ULbox,i,j,k);
                Cell UStar(UStarbox,i,j,k);


                SR = std::max(std::abs(UL(VELOCITY,0,d))+UL(SOUNDSPEED),std::abs(UR(VELOCITY,0,d))+UR(SOUNDSPEED));
                SL = -SR;

                Sstar = getSstar(UL,UR,SL,SR,i,j,k,d);

                fluxbox(i,j,k,USTAR) = Sstar;

                if(SL>=0.0)
                {
                    getStarState(UL,UStar,SL,Sstar,parameters,d);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,UStar,d);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL,d);
                        }
                    }
                }
                else if(Sstar>=0.0)
                {
                    getStarState(UL,UStar,SL,Sstar,parameters,d);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,UStar,d);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL,d)+SL*(UStar(n)-UL(n));
                        }
                    }

                }
                else if(SR>=0.0)
                {
                    getStarState(UR,UStar,SR,Sstar,parameters,d);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,UStar,d);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UR,d)+SR*(UStar(n)-UR(n));
                        }
                    }
                }
                else
                {
                    getStarState(UR,UStar,SR,Sstar,parameters,d);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,UStar,d);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UR,d);
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
void update(BoxAccessCellArray& fluxbox, BoxAccessCellArray& Ubox, BoxAccessCellArray& U1box, ParameterStruct& parameters, Direction_enum d)
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
                        U1box(i,j,k,n) += (parameters.dt/parameters.dx[d])*(fluxbox(i,j,k,n) - fluxbox.right(d,i,j,k,n) -Ubox(i,j,k,n)*(fluxbox(i,j,k,USTAR) - fluxbox.right(d,i,j,k,USTAR)));
                    }
                    else
                    {
                        U1box(i,j,k,n) += (parameters.dt/parameters.dx[d])*(fluxbox(i,j,k,n) - fluxbox.right(d,i,j,k,n));
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
void MUSCLextrapolate(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& grad, Direction_enum d)
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
                    r = (U(i,j,k,n)-U.left(d,i,j,k,n))/(U.right(d,i,j,k,n)-U(i,j,k,n));

                    if(std::isinf(r) || std::isnan(r))
                    {
                        r = 0.0;
                    }

                    grad(i,j,k,n) = vanLeerlimiter(r)*0.5*(U.right(d,i,j,k,n)-U.left(d,i,j,k,n));

                    UL(i,j,k,n) = U(i,j,k,n) - 0.5*grad(i,j,k,n);
                    UR(i,j,k,n) = U(i,j,k,n) + 0.5*grad(i,j,k,n);

                    //UL(i,j,k,n) = U(i,j,k,n);
                    //UR(i,j,k,n) = U(i,j,k,n);

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
    Direction_enum d;

    U1 = U;

    for(int dir = 0; dir < AMREX_SPACEDIM ; dir++)
    {
        d = (Direction_enum)dir;

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

            MUSCLextrapolate(Ubox,ULbox,URbox,gradbox,d);

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

            FArrayBox& flux_fab   = flux_arr[d][mfi];

            BoxAccessCellArray Ubox(mfi,bx,U);
            BoxAccessCellArray U1box(mfi,bx,U1);
            BoxAccessCellArray ULbox(mfi,bx,UL);
            BoxAccessCellArray URbox(mfi,bx,UR);
            BoxAccessCellArray UStarbox(mfi,bx,UStar);
            BoxAccessCellArray fluxbox(bx,flux_fab,U);

            calc_fluxes (fluxbox, ULbox, URbox, UStarbox, parameters,d);
            update      (fluxbox, Ubox, U1box, parameters,d);


        }


    }

    U1.conservativeToPrimitive(parameters);
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


