#include "simulationheader.h"

/** Calculates the contact wave speed estimate.
 */
double getSstar(Cell& UL, Cell& UR, Real SL, Real SR, int i, int j, int k, Direction_enum d)
{
    return (-UR(SIGMA,0,d,d)+UL(SIGMA,0,d,d)+UL(RHOU,0,d)*(SL-UL(VELOCITY,0,d))-UR(RHOU,0,d)*(SR-UR(VELOCITY,0,d)))/(UL(RHO)*(SL-UL(VELOCITY,0,d)) -  UR(RHO)*(SR-UR(VELOCITY,0,d)));
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

        if(parameters.materialInfo[m].mixture)
        {
            UStar(ALPHARHOLAMBDA,m) = multiplier*U(ALPHARHOLAMBDA,m);
        }

        if(parameters.materialInfo[m].plastic)
        {
            UStar(ALPHARHOEPSILON,m) = multiplier*U(ALPHARHOEPSILON,m);
        }

        UStar(RHO)             += UStar(ALPHARHO,m);

    }

    for(int row = 0; row < U.numberOfComponents; row++)
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

        if(parameters.SOLID)
        {
            for(int col=0;col<U.numberOfComponents; col++)
            {
                if(row == d)
                {
                    UStar(V_TENSOR,0,row,col) = U(V_TENSOR,0,row,col);
                }
                else
                {
                    UStar(V_TENSOR,0,row,col) = multiplier*U(V_TENSOR,0,row,col);
                }
            }
        }
    }


    UStar(TOTAL_E)      = multiplier*U(RHO)*(U(TOTAL_E)/U(RHO) + (Sstar-U(VELOCITY,0,d))*(Sstar - (U(SIGMA,0,d,d))/(U(RHO)*(SK-U(VELOCITY,0,d))) ));

    return;
}

void getSigmaStar(Cell& UKStar, Real Sstar, Real SLT, Real SRT, Cell& UL, Cell& UR, Cell& ULStar, Cell& URStar, Direction_enum d, ParameterStruct& parameters)
{
    for(int i=0;i<UL.numberOfComponents;i++)
    {
        UKStar(SIGMA,0,i,d) = (ULStar(RHO)*(Sstar-SLT)*URStar(RHO)*(Sstar-SRT)*(UL(VELOCITY,0,i) -UR(VELOCITY,0,i)) + ULStar(RHO)*(Sstar-SLT)*UR(SIGMA,0,i,d) - URStar(RHO)*(Sstar-SRT)*UL(SIGMA,0,i,d))/(ULStar(RHO)*(Sstar-SLT)-URStar(RHO)*(Sstar-SRT));
    }

    return;
}

/** Finds the (inner) intermediate states for the HLLD solver.
 */
void getStarStarState(Cell& UL, Cell& UR, Cell& ULStar, Cell& URStar, Cell& UStarStar, Real SL, Real SR, Real SLT, Real SRT, Real Sstar, ParameterStruct& parameters, Direction_enum d, Cell& UK, Cell& UKStar, Real SKT)
{
    int N = UL.numberOfComponents;

    UStarStar = UKStar;

    getSigmaStar(UKStar,Sstar,SLT,SRT,UL,UR,ULStar,URStar,d,parameters);

    for(int row =0;row<N;row++)
    {
        if(row == d)
        {
            continue;
        }
        else
        {
            UStarStar(RHOU,0,row) += (UKStar(SIGMA,0,row,d)-UK(SIGMA,0,row,d))/(Sstar-SKT);

        }

        UStarStar(VELOCITY,0,row) = UStarStar(RHOU,0,row)/UStarStar(RHO);

    }

    for(int row =0;row<N;row++)
    {
        if(row == d)
        {
            continue;
        }
        else
        {
            UStarStar(TOTAL_E) += (UStarStar(VELOCITY,0,row)*UKStar(SIGMA,0,row,d)-UK(VELOCITY,0,row)*UK(SIGMA,0,row,d))/(Sstar-SKT);

            for(int col =0;col<N;col++)
            {
                UStarStar(V_TENSOR,0,row,col) += (UStarStar(V_TENSOR,0,d,col)*(UStarStar(VELOCITY,0,row)-UK(VELOCITY,0,row)))/(Sstar-SKT);
            }
        }
    }

    return;
}

/** Calculates the HLLC fluxes.
 */
void calc_5Wave_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& ULStarbox, BoxAccessCellArray& URStarbox, BoxAccessCellArray& UStarStarbox, ParameterStruct& parameters, Direction_enum d)
{
    const auto lo = lbound(ULbox.box);
    const auto hi = ubound(ULbox.box);

    Real SR     = 0.0;
    Real SL     = 0.0;
    Real SLT    = 0.0;
    Real SRT    = 0.0;
    Real Sstar  = 0.0;

    IntVect extra(AMREX_D_DECL(0,0,0));

    extra[d]=1;

    for 		(int k = lo.z; k <= hi.z+extra[z]; ++k)
    {
        for 	(int j = lo.y; j <= hi.y+extra[y]; ++j)
        {
            for (int i = lo.x; i <= hi.x+extra[x]; ++i)
            {

                Cell UL(URbox,i-extra[x],j-extra[y],k-extra[z],solid);
                Cell UR(ULbox,i,j,k,solid);

                Cell ULStar(ULStarbox,i,j,k,solid);
                Cell URStar(URStarbox,i,j,k,solid);

                Cell UStarStar(UStarStarbox,i,j,k,solid);


                SR = std::max(std::abs(UL(VELOCITY,0,d))+UL(SOUNDSPEED),std::abs(UR(VELOCITY,0,d))+UR(SOUNDSPEED));
                SL = -SR;

                Sstar = getSstar(UL,UR,SL,SR,i,j,k,d);

                if(std::isnan(SL) || std::isnan(SR) || std::isnan(Sstar))
                {
                    UL.parent->checkLimits(UL.accessPattern.allVariables);
                    UR.parent->checkLimits(UL.accessPattern.allVariables);

                    SR = std::max(std::abs(UL(VELOCITY,0,d))+UL(SOUNDSPEED),std::abs(UR(VELOCITY,0,d))+UR(SOUNDSPEED));
                    SL = -SR;

                    Sstar = getSstar(UL,UR,SL,SR,i,j,k,d);

                    if(std::isnan(SL) || std::isnan(SR) || std::isnan(Sstar))
                    {
                        amrex::Abort("Nan in wavespeeds before flux");
                    }
                }

                fluxbox(i,j,k,USTAR) = Sstar;

                if(SL>=0.0)
                {
                    getStarState(UL,ULStar,SL,Sstar,parameters,d);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,ULStar,d);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL,d);
                        }
                    }

                    for(int row=0;row<ULbox.numberOfComponents;row++)
                    {
                        for(int col=0;col<ULbox.numberOfComponents;col++)
                        {
                            fluxbox(i,j,k,VSTAR,0,row,col) = UL(V_TENSOR,0,row,col);
                        }
                    }
                }
                else if(Sstar>=0.0)
                {
                    getStarState(UL,ULStar,SL,Sstar,parameters,d);

                    SLT = Sstar - ULStar.parent->transverseWaveSpeed(ULStar.parent_i,ULStar.parent_j,ULStar.parent_k);

                    if(SLT>=0.0)
                    {
                        for(auto n : ULbox.accessPattern.conservativeVariables)
                        {
                            if(n.var == ALPHA)
                            {
                                 fluxbox(i,j,k,n)	= flux(n,ULStar,d);
                            }
                            else
                            {
                                fluxbox(i,j,k,n)	= flux(n,UL,d)+SL*(ULStar(n)-UL(n));
                            }
                        }

                        for(int row=0;row<ULbox.numberOfComponents;row++)
                        {
                            for(int col=0;col<ULbox.numberOfComponents;col++)
                            {
                                fluxbox(i,j,k,VSTAR,0,row,col) = ULStar(V_TENSOR,0,row,col);
                            }
                        }
                    }
                    else
                    {
                        getStarState(UR,URStar,SR,Sstar,parameters,d);

                        SRT = Sstar + URStar.parent->transverseWaveSpeed(URStar.parent_i,URStar.parent_j,URStar.parent_k);

                        getStarStarState(UL,UR,ULStar,URStar,UStarStar,SL,SR,SLT,SRT,Sstar,parameters,d,UL,ULStar,SLT);

                        for(auto n : ULbox.accessPattern.conservativeVariables)
                        {
                            if(n.var == ALPHA)
                            {
                                 fluxbox(i,j,k,n)	= flux(n,ULStar,d);
                            }
                            else
                            {
                                fluxbox(i,j,k,n)	= flux(n,UL,d)+SL*(ULStar(n)-UL(n))+SLT*(UStarStar(n)-ULStar(n));
                            }
                        }

                        for(int row=0;row<ULbox.numberOfComponents;row++)
                        {
                            for(int col=0;col<ULbox.numberOfComponents;col++)
                            {
                                fluxbox(i,j,k,VSTAR,0,row,col) = UStarStar(V_TENSOR,0,row,col);
                            }
                        }
                    }
                }
                else if(SR>=0.0)
                {
                    getStarState(UR,URStar,SR,Sstar,parameters,d);

                    SRT = Sstar + URStar.parent->transverseWaveSpeed(URStar.parent_i,URStar.parent_j,URStar.parent_k); 

                    if(SRT>=0.0)
                    {
                        getStarState(UL,ULStar,SL,Sstar,parameters,d);

                        SLT = Sstar - ULStar.parent->transverseWaveSpeed(ULStar.parent_i,ULStar.parent_j,ULStar.parent_k);

                        getStarStarState(UL,UR,ULStar,URStar,UStarStar,SL,SR,SLT,SRT,Sstar,parameters,d,UR,URStar,SRT);

                        for(auto n : ULbox.accessPattern.conservativeVariables)
                        {
                            if(n.var == ALPHA)
                            {
                                 fluxbox(i,j,k,n)	= flux(n,URStar,d);
                            }
                            else
                            {
                                fluxbox(i,j,k,n)	= flux(n,UR,d)+SR*(URStar(n)-UR(n))+SRT*(UStarStar(n)-URStar(n));
                            }
                        }

                        for(int row=0;row<ULbox.numberOfComponents;row++)
                        {
                            for(int col=0;col<ULbox.numberOfComponents;col++)
                            {
                                fluxbox(i,j,k,VSTAR,0,row,col) = UStarStar(V_TENSOR,0,row,col);
                            }
                        }
                    }
                    else
                    {
                        for(auto n : ULbox.accessPattern.conservativeVariables)
                        {
                            if(n.var == ALPHA)
                            {
                                 fluxbox(i,j,k,n)	= flux(n,URStar,d);
                            }
                            else
                            {
                                fluxbox(i,j,k,n)	= flux(n,UR,d)+SR*(URStar(n)-UR(n));
                            }
                        }

                        for(int row=0;row<ULbox.numberOfComponents;row++)
                        {
                            for(int col=0;col<ULbox.numberOfComponents;col++)
                            {
                                fluxbox(i,j,k,VSTAR,0,row,col) = URStar(V_TENSOR,0,row,col);
                            }
                        }
                    }
                }
                else
                {
                    getStarState(UR,URStar,SR,Sstar,parameters,d);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,URStar,d);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UR,d);
                        }
                    }

                    for(int row=0;row<ULbox.numberOfComponents;row++)
                    {
                        for(int col=0;col<ULbox.numberOfComponents;col++)
                        {
                            fluxbox(i,j,k,VSTAR,0,row,col) = UR(V_TENSOR,0,row,col);
                        }
                    }
                }

                if(std::isnan(SL) || std::isnan(SR) || std::isnan(Sstar) || std::isnan(SLT) || std::isnan(SRT))
                {
                    Print() << "Nan in Wavespeeds " << std::endl;
                    Print() << SL << " " <<  SR << " " << Sstar << " " << SLT  << " "  << SRT << std::endl;
                    Print() << "At " <<  i << " " << j << " " <<k  << std::endl;
                    amrex::Abort("Nan in waveSpeeds");
                }

            }
        }
    }

    return;
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

                Cell UL(URbox,i-extra[x],j-extra[y],k-extra[z],fluid);
                Cell UR(ULbox,i,j,k,fluid);
                Cell UStar(UStarbox,i,j,k,fluid);



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
void update(BoxAccessCellArray& fluxbox, BoxAccessCellArray& Ubox, BoxAccessCellArray& U1box, ParameterStruct& parameters, Direction_enum d, Real dt, const Real* dx)
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
                        U1box(i,j,k,n) += (dt/dx[d])*(fluxbox(i,j,k,n) - fluxbox.right(d,i,j,k,n) -Ubox(i,j,k,n)*(fluxbox(i,j,k,USTAR) - fluxbox.right(d,i,j,k,USTAR)));
                    }
                    else if(n.var == V_TENSOR)
                    {
                        U1box(i,j,k,n) += (dt/dx[d])*(fluxbox(i,j,k,n) - fluxbox.right(d,i,j,k,n) -(2.0/3.0)*Ubox(i,j,k,n)*(fluxbox(i,j,k,USTAR) - fluxbox.right(d,i,j,k,USTAR)) + Ubox(i,j,k,VELOCITY,0,n.row)*(fluxbox(i,j,k,VSTAR,0,d,n.col) - fluxbox.right(d,i,j,k,VSTAR,0,d,n.col)));
                    }
                    else
                    {
                        U1box(i,j,k,n) += (dt/dx[d])*(fluxbox(i,j,k,n) - fluxbox.right(d,i,j,k,n));
                    }
                }
            }
        }
    }

    U1box.checkLimits(Ubox.accessPattern.conservativeVariables);
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

                }
            }
        }
    }

    UL.checkLimits(U.accessPattern.primitiveVariables);
    UR.checkLimits(U.accessPattern.primitiveVariables);

    return;
}

/** Calculate the HLLC flux and update new array
 */
void HLLCadvance(CellArray& U,CellArray& U1, CellArray& UL, CellArray& UR, CellArray& MUSCLgrad, CellArray& ULStar, CellArray& URStar, CellArray& UStarStar, Array<MultiFab, AMREX_SPACEDIM>& flux_arr,Geometry const& geom, ParameterStruct& parameters,Vector<BCRec>& bc, THINCArray& THINC)
{
    Direction_enum d;

    U1 = U;

    for(int dir = 0; dir < AMREX_SPACEDIM ; dir++)
    {
        d = (Direction_enum)dir;

        UL = U;
        UR = U;

        /*-------------------------------------------------------------
         * Perform MUSCL extrapolation.
         * -----------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(U.data); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            /*-------------------------------------------------------------
             * Data can't be accessed straight from a Multifab so we make
             * some wrappers to hold the FArrayBoxes that can access the
             * data called BoxAccessCellArray.
             * -----------------------------------------------------------*/

            BoxAccessCellArray  Ubox(mfi,bx,U);
            BoxAccessCellArray  ULbox(mfi,bx,UL);
            BoxAccessCellArray  URbox(mfi,bx,UR);


            if(parameters.MUSCL)
            {
                BoxAccessCellArray  gradbox(mfi,bx,MUSCLgrad);

                MUSCLextrapolate(Ubox,ULbox,URbox,gradbox,d);
            }

            if(parameters.THINC) // && d == x
            {
                BoxAccessCellArray ULTHINC(mfi,bx,ULStar);
                BoxAccessCellArray URTHINC(mfi,bx,URStar);
                BoxAccessTHINCArray THINCbox(mfi,bx,THINC);

                THINCbox.THINCreconstruction(Ubox,ULbox,URbox,ULTHINC,URTHINC,parameters, parameters.dx.dataPtr(),d);

                UL.cleanUpAlpha();
                UR.cleanUpAlpha();

            }


            ULbox.primitiveToConservative();
            URbox.primitiveToConservative();

            ULbox.getSoundSpeed();
            URbox.getSoundSpeed();

        }

        FillDomainBoundary(UL.data, geom, bc);
        FillDomainBoundary(UR.data, geom, bc);

        UL.data.FillBoundary(geom.periodicity());
        UR.data.FillBoundary(geom.periodicity());


        /*-------------------------------------------------------------
         * Calulate HLLC flux and update the new array.
         * -----------------------------------------------------------*/

#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(U.data); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            FArrayBox& flux_fab   = flux_arr[d][mfi];

            BoxAccessCellArray Ubox(mfi,bx,U);
            BoxAccessCellArray U1box(mfi,bx,U1);
            BoxAccessCellArray ULbox(mfi,bx,UL);
            BoxAccessCellArray URbox(mfi,bx,UR);
            BoxAccessCellArray ULStarbox(mfi,bx,ULStar);
            BoxAccessCellArray URStarbox(mfi,bx,URStar);
            BoxAccessCellArray UStarStarbox(mfi,bx,UStarStar);
            BoxAccessCellArray fluxbox(bx,flux_fab,U); 

            if(parameters.SOLID)
            {
                calc_5Wave_fluxes(fluxbox, ULbox, URbox, ULStarbox, URStarbox, UStarStarbox, parameters,d);
            }
            else
            {
                calc_fluxes(fluxbox, ULbox, URbox, ULStarbox, parameters,d);
            }

            update(fluxbox, Ubox, U1box, parameters,d,parameters.dt,parameters.dx.dataPtr());

        }



    }


    U1.conservativeToPrimitive();
    FillDomainBoundary(U1.data, geom, bc);
    U1.data.FillBoundary(geom.periodicity());



}

/** 2nd Order Runge-Kutta time integration
 */
void advance(CellArray& U, CellArray& U1, CellArray& U2, CellArray& MUSCLgrad, CellArray& UL, CellArray& UR, CellArray& ULStar, CellArray& URStar, CellArray& UStarStar, Array<MultiFab, AMREX_SPACEDIM>& flux_arr, Geometry const& geom, ParameterStruct& parameters, Vector<BCRec>& bc, THINCArray &THINC)
{
    HLLCadvance(U,  U1, UL, UR, MUSCLgrad, ULStar, URStar, UStarStar, flux_arr, geom, parameters, bc, THINC);

    HLLCadvance(U1, U2, UL, UR, MUSCLgrad, ULStar, URStar, UStarStar, flux_arr, geom, parameters, bc, THINC);

    U1 = ((U*(1.0/2.0))+(U2*(1.0/2.0)));

}


