#include "simulationheader.h"

/** Calculates the contact wave speed estimate.
 */
double getSstar(Cell& UL, Cell& UR, Real SL, Real SR, int i, int j, int k, Direction_enum d, int m)
{
    if( std::abs(UL(RHO,m)*(SL-UL(VELOCITY,m,d)) -  UR(RHO,m)*(SR-UR(VELOCITY,m,d))) < 1E-10 )
    {
        return 0.0;
    }
    else
    {
        return (-UR(SIGMA,m,d,d)+UL(SIGMA,m,d,d)+UL(RHOU,m,d)*(SL-UL(VELOCITY,m,d))-UR(RHOU,m,d)*(SR-UR(VELOCITY,m,d)))/(UL(RHO,m)*(SL-UL(VELOCITY,m,d)) -  UR(RHO,m)*(SR-UR(VELOCITY,m,d)));
    }
}

/** Calculates the HLLC intermediate state.
 */
void getStarState(Cell& U, Cell& UStar, Real SK, Real Sstar, ParameterStruct& parameters, Direction_enum d, int m)
{
    static const Real tolerance  = 1E-20;

    Real multiplier = ( std::abs((Sstar-U(VELOCITY,m,d))/Sstar) < tolerance ? 1.0 : (SK-U(VELOCITY,m,d))/(SK-Sstar)   );

    UStar(RHO,m)       = multiplier*U(RHO,m);

    for(int row = 0; row < U.numberOfComponents; row++)
    {
        if(row == d)
        {
            UStar(RHOU,m,row)         = multiplier*U(RHO,m)*Sstar;
        }
        else
        {
            UStar(RHOU,m,row)         = multiplier*U(RHOU,m,row);
        }

        UStar(VELOCITY,m,row)     = UStar(RHOU,m,row)/UStar(RHO,m);
    }


    UStar(TOTAL_E,m)      = multiplier*U(RHO,m)*(U(TOTAL_E,m)/U(RHO,m) + (Sstar-U(VELOCITY,m,d))*(Sstar - (U(SIGMA,m,d,d))/(U(RHO,m)*(SK-U(VELOCITY,m,d))) ));

    return;
}

/*void getSigmaStar(Cell& UKStar, Real Sstar, Real SLT, Real SRT, Cell& UL, Cell& UR, Cell& ULStar, Cell& URStar, Direction_enum d, ParameterStruct& parameters)
{
    for(int i=0;i<UL.numberOfComponents;i++)
    {
        UKStar(SIGMA,0,i,d) = (ULStar(RHO)*(Sstar-SLT)*URStar(RHO)*(Sstar-SRT)*(UL(VELOCITY,0,i) -UR(VELOCITY,0,i)) + ULStar(RHO)*(Sstar-SLT)*UR(SIGMA,0,i,d) - URStar(RHO)*(Sstar-SRT)*UL(SIGMA,0,i,d))/(ULStar(RHO)*(Sstar-SLT)-URStar(RHO)*(Sstar-SRT));
    }

    return;
}*/

/** Finds the (inner) intermediate states for the HLLD solver.
 */
/*void getStarStarState(Cell& UL, Cell& UR, Cell& ULStar, Cell& URStar, Cell& UStarStar, Real SL, Real SR, Real SLT, Real SRT, Real Sstar, ParameterStruct& parameters, Direction_enum d, Cell& UK, Cell& UKStar, Real SKT)
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
*/
/** Calculates the HLLC fluxes.
 */
/*void calc_5Wave_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& ULStarbox, BoxAccessCellArray& URStarbox, BoxAccessCellArray& UStarStarbox, ParameterStruct& parameters, Direction_enum d, const Real *dx, const Real *prob_lo)
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
                    UL.parent->checking = 1;
                    UR.parent->checking = 1;

                    UL.parent->checkLimits(UL.accessPattern.allVariables);
                    UR.parent->checkLimits(UL.accessPattern.allVariables);

                    UL.parent->conservativeToPrimitive();
                    UR.parent->conservativeToPrimitive();

                    UL.parent->checking = 0;
                    UR.parent->checking = 0;

                    SR = std::max(std::abs(UL(VELOCITY,0,d))+UL(SOUNDSPEED),std::abs(UR(VELOCITY,0,d))+UR(SOUNDSPEED));
                    SL = -SR;

                    Sstar = getSstar(UL,UR,SL,SR,i,j,k,d);

                    if(std::isnan(SL) || std::isnan(SR))
                    {
                        amrex::Abort("Nan in SL SR wavespeeds before flux");
                    }

                    if(std::isnan(Sstar))
                    {
                        if(std::isnan(-UR(SIGMA,0,d,d)+UL(SIGMA,0,d,d)+UL(RHOU,0,d)*(SL-UL(VELOCITY,0,d))-UR(RHOU,0,d)*(SR-UR(VELOCITY,0,d))))
                        {
                            if(std::isnan(-UR(SIGMA,0,d,d)+UL(SIGMA,0,d,d)))
                            {
                                if(std::isnan(UR(P)) || std::isnan(UL(P)))
                                {
                                    amrex::Abort("Nan in Sstar wavespeed before flux: num, sigma, p");
                                }
                                else
                                {
                                    if( (UR(ALPHA,0) < UR(ALPHA,1)) && (UL(ALPHA,0) < UL(ALPHA,1)))
                                    {
                                        for(int row = 0; row<3;row++)
                                        {
                                            for(int col = 0; col<3;col++)
                                            {
                                                UR(SIGMA,0,row,col) = -UR(P)*delta<Real>(row,col);
                                                UL(SIGMA,0,row,col) = -UL(P)*delta<Real>(row,col);
                                            }
                                        }

                                        Sstar = getSstar(UL,UR,SL,SR,i,j,k,d);

                                     }
                                    else
                                    {
                                        amrex::Abort("Nan in Sstar wavespeed before flux: num, sigma, not p, in solid");
                                    }
                                }

                            }
                            else if(std::isnan(UL(RHOU,0,d)*(SL-UL(VELOCITY,0,d))-UR(RHOU,0,d)*(SR-UR(VELOCITY,0,d))))
                            {
                                amrex::Abort("Nan in Sstar wavespeed before flux: num, RHOU or vel");
                            }
                            else
                            {
                                Sstar = 0.0;
                            }
                        }
                        else if(std::isnan((UL(RHO)*(SL-UL(VELOCITY,0,d)) -  UR(RHO)*(SR-UR(VELOCITY,0,d)))))
                        {
                            amrex::Abort("Nan in Sstar wavespeed before flux: den");
                        }
                        else
                        {
                            Sstar = 0.0;
                        }
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
}*/

/** Calculates the HLLC fluxes.
 */
void calc_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, ParameterStruct& parameters, Direction_enum d, const Real *dx, const Real *prob_lo)
{
    const auto lo = lbound(ULbox.box);
    const auto hi = ubound(ULbox.box);

    Vector<Real> SR(parameters.numberOfMaterials);
    Vector<Real> SL(parameters.numberOfMaterials);
    Vector<Real> Sstar(parameters.numberOfMaterials);

    IntVect extra(AMREX_D_DECL(0,0,0));

    extra[d]=1;

    for 		   (int k = lo.z; k <= hi.z+extra[z]; ++k)
    {
        for 	   (int j = lo.y; j <= hi.y+extra[y]; ++j)
        {
            for    (int i = lo.x; i <= hi.x+extra[x]; ++i)
            {
                Cell UL(URbox,i-extra[x],j-extra[y],k-extra[z],fluid);
                Cell UR(ULbox,i,j,k,fluid);
                Cell UStar(UStarbox,i,j,k,fluid);

                for(int m = 0; m < parameters.numberOfMaterials; m++)
                {
                    SR[m] = std::max(std::abs(UL(VELOCITY,m,d))+UL(SOUNDSPEED,m),std::abs(UR(VELOCITY,m,d))+UR(SOUNDSPEED,m));
                    SL[m] = -SR[m];

                    Sstar[m] = getSstar(UL,UR,SL[m],SR[m],i,j,k,d,m);

                    fluxbox(i,j,k,USTAR,m) = Sstar[m];

                    if(SL[m]>=0.0)
                    {
                        for(auto n : ULbox.accessPattern.material_conservativeVariables[m])
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL,d);
                        }
                    }
                    else if(Sstar[m]>=0.0)
                    {
                        getStarState(UL,UStar,SL[m],Sstar[m],parameters,d,m);

                        for(auto n : ULbox.accessPattern.material_conservativeVariables[m])
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL,d)+SL[m]*(UStar(n)-UL(n));
                        }

                    }
                    else if(SR[m]>=0.0)
                    {
                        getStarState(UR,UStar,SR[m],Sstar[m],parameters,d,m);

                        for(auto n : ULbox.accessPattern.material_conservativeVariables[m])
                        {
                            fluxbox(i,j,k,n)	= flux(n,UR,d)+SR[m]*(UStar(n)-UR(n));
                        }
                    }
                    else
                    {
                        for(auto n : ULbox.accessPattern.material_conservativeVariables[m])
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
                    U1box(i,j,k,n) += (dt/dx[d])*(fluxbox(i,j,k,n) - fluxbox.right(d,i,j,k,n));
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



