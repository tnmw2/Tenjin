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

    if(parameters.materialInfo[m].plastic)
    {
        UStar(RHOEPSILON,m) = multiplier*U(RHOEPSILON,m);
    }

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

        if(U.accessPattern.materialInfo[m].phase == solid)
        {
            for(int col=0;col<U.numberOfComponents; col++)
            {
                if(row == d)
                {
                    UStar(V_TENSOR,m,row,col) = U(V_TENSOR,m,row,col);
                }
                else
                {
                    UStar(V_TENSOR,m,row,col) = multiplier*U(V_TENSOR,m,row,col);
                }
            }
        }
    }

    UStar(TOTAL_E,m)      = multiplier*U(RHO,m)*(U(TOTAL_E,m)/U(RHO,m) + (Sstar-U(VELOCITY,m,d))*(Sstar - (U(SIGMA,m,d,d))/(U(RHO,m)*(SK-U(VELOCITY,m,d))) ));

    return;
}

void getSigmaStar(Cell& UKStar, Real Sstar, Real SLT, Real SRT, Cell& UL, Cell& UR, Cell& ULStar, Cell& URStar, Direction_enum d, ParameterStruct& parameters, int m)
{
    for(int i=0;i<UL.numberOfComponents;i++)
    {
        UKStar(SIGMA,m,i,d) = (ULStar(RHO,m)*(Sstar-SLT)*URStar(RHO,m)*(Sstar-SRT)*(UL(VELOCITY,m,i) -UR(VELOCITY,m,i)) + ULStar(RHO,m)*(Sstar-SLT)*UR(SIGMA,m,i,d) - URStar(RHO,m)*(Sstar-SRT)*UL(SIGMA,m,i,d))/(ULStar(RHO,m)*(Sstar-SLT)-URStar(RHO,m)*(Sstar-SRT));
    }

    return;
}

/** Finds the (inner) intermediate states for the HLLD solver.
 */
void getStarStarState(Cell& UL, Cell& UR, Cell& ULStar, Cell& URStar, Cell& UStarStar, Real SL, Real SR, Real SLT, Real SRT, Real Sstar, ParameterStruct& parameters, Direction_enum d, Cell& UK, Cell& UKStar, Real SKT, int m)
{
    int N = UL.numberOfComponents;

    UStarStar = UKStar;

    getSigmaStar(UKStar,Sstar,SLT,SRT,UL,UR,ULStar,URStar,d,parameters,m);

    for(int row =0;row<N;row++)
    {
        if(row == d)
        {
            continue;
        }
        else
        {
            UStarStar(RHOU,m,row) += (UKStar(SIGMA,m,row,d)-UK(SIGMA,m,row,d))/(Sstar-SKT);

        }

        UStarStar(VELOCITY,m,row) = UStarStar(RHOU,m,row)/UStarStar(RHO,m);

    }

    for(int row =0;row<N;row++)
    {
        if(row == d)
        {
            continue;
        }
        else
        {
            UStarStar(TOTAL_E,m) += (UStarStar(VELOCITY,m,row)*UKStar(SIGMA,m,row,d)-UK(VELOCITY,m,row)*UK(SIGMA,m,row,d))/(Sstar-SKT);

            for(int col =0;col<N;col++)
            {
                UStarStar(V_TENSOR,m,row,col) += (UStarStar(V_TENSOR,m,d,col)*(UStarStar(VELOCITY,m,row)-UK(VELOCITY,m,row)))/(Sstar-SKT);
            }
        }
    }

    return;
}

/** Calculates the HLLC fluxes.
 */
void calc_5Wave_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& ULStarbox, BoxAccessCellArray& URStarbox, BoxAccessCellArray& UStarStarbox, ParameterStruct& parameters, Direction_enum d, const Real *dx, const Real *prob_lo)
{
    const auto lo = lbound(ULbox.box);
    const auto hi = ubound(ULbox.box);

    Real SR     = 0.0;
    Real SL     = 0.0;
    Real SLT    = 0.0;
    Real SRT    = 0.0;
    Real Sstar  = 0.0;

    int extra[3] = {0,0,0};

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

                for(int m = 0; m < parameters.numberOfMaterials; m++)
                {
                    SR = std::max(std::abs(UL(VELOCITY,m,d))+UL(SOUNDSPEED,m),std::abs(UR(VELOCITY,m,d))+UR(SOUNDSPEED,m));
                    SL = -SR;

                    Sstar = getSstar(UL,UR,SL,SR,i,j,k,d,m);


                    if(std::isnan(SL) || std::isnan(SR))
                    {
                        Print() << "material: " << m << std::endl;
                        Print() << SL << " " << Sstar << " " << SR << std::endl;
                        Print() << UL(VELOCITY,m,d) << " " << UR(VELOCITY,m,d) << std::endl;
                        Print() << UL(SOUNDSPEED,m) << " " << UR(SOUNDSPEED,m) << std::endl;

                       amrex::Abort("Nan in SL SR wavespeeds before flux");
                    }

                    if(std::isnan(Sstar))
                    {
                       amrex::Abort("Nan in SStar wavespeeds before flux");
                    }

                    fluxbox(i,j,k,USTAR,m) = Sstar;

                    if(SL>=0.0)
                    {
                        for(auto n : ULbox.accessPattern.material_conservativeVariables[m])
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL,d);
                        }

                        for(int row=0;row<ULbox.numberOfComponents;row++)
                        {
                            for(int col=0;col<ULbox.numberOfComponents;col++)
                            {
                                fluxbox(i,j,k,VSTAR,m,row,col) = UL(V_TENSOR,m,row,col);
                            }
                        }
                    }
                    else if(Sstar>=0.0)
                    {
                        getStarState(UL,ULStar,SL,Sstar,parameters,d,m);

                        SLT = Sstar - ULStar.parent->transverseWaveSpeed(ULStar.parent_i,ULStar.parent_j,ULStar.parent_k,m);

                        if(SLT>=0.0)
                        {
                            for(auto n : ULbox.accessPattern.material_conservativeVariables[m])
                            {
                                fluxbox(i,j,k,n)	= flux(n,UL,d)+SL*(ULStar(n)-UL(n));
                            }

                            for(int row=0;row<ULbox.numberOfComponents;row++)
                            {
                                for(int col=0;col<ULbox.numberOfComponents;col++)
                                {
                                    fluxbox(i,j,k,VSTAR,m,row,col) = ULStar(V_TENSOR,m,row,col);
                                }
                            }
                        }
                        else
                        {
                            getStarState(UR,URStar,SR,Sstar,parameters,d,m);

                            SRT = Sstar + URStar.parent->transverseWaveSpeed(URStar.parent_i,URStar.parent_j,URStar.parent_k,m);

                            getStarStarState(UL,UR,ULStar,URStar,UStarStar,SL,SR,SLT,SRT,Sstar,parameters,d,UL,ULStar,SLT,m);

                            for(auto n : ULbox.accessPattern.material_conservativeVariables[m])
                            {
                                fluxbox(i,j,k,n)	= flux(n,UL,d)+SL*(ULStar(n)-UL(n))+SLT*(UStarStar(n)-ULStar(n));
                            }

                            for(int row=0;row<ULbox.numberOfComponents;row++)
                            {
                                for(int col=0;col<ULbox.numberOfComponents;col++)
                                {
                                    fluxbox(i,j,k,VSTAR,m,row,col) = UStarStar(V_TENSOR,m,row,col);
                                }
                            }
                        }
                    }
                    else if(SR>=0.0)
                    {
                        getStarState(UR,URStar,SR,Sstar,parameters,d,m);

                        SRT = Sstar + URStar.parent->transverseWaveSpeed(URStar.parent_i,URStar.parent_j,URStar.parent_k,m);

                        if(SRT>=0.0)
                        {
                            getStarState(UL,ULStar,SL,Sstar,parameters,d,m);

                            SLT = Sstar - ULStar.parent->transverseWaveSpeed(ULStar.parent_i,ULStar.parent_j,ULStar.parent_k,m);

                            getStarStarState(UL,UR,ULStar,URStar,UStarStar,SL,SR,SLT,SRT,Sstar,parameters,d,UR,URStar,SRT,m);

                            for(auto n : ULbox.accessPattern.material_conservativeVariables[m])
                            {
                                fluxbox(i,j,k,n)	= flux(n,UR,d)+SR*(URStar(n)-UR(n))+SRT*(UStarStar(n)-URStar(n));
                            }

                            for(int row=0;row<ULbox.numberOfComponents;row++)
                            {
                                for(int col=0;col<ULbox.numberOfComponents;col++)
                                {
                                    fluxbox(i,j,k,VSTAR,m,row,col) = UStarStar(V_TENSOR,m,row,col);
                                }
                            }
                        }
                        else
                        {
                            for(auto n : ULbox.accessPattern.material_conservativeVariables[m])
                            {
                                fluxbox(i,j,k,n)	= flux(n,UR,d)+SR*(URStar(n)-UR(n));
                            }

                            for(int row=0;row<ULbox.numberOfComponents;row++)
                            {
                                for(int col=0;col<ULbox.numberOfComponents;col++)
                                {
                                    fluxbox(i,j,k,VSTAR,m,row,col) = URStar(V_TENSOR,m,row,col);
                                }
                            }
                        }
                    }
                    else
                    {
                        getStarState(UR,URStar,SR,Sstar,parameters,d,m);

                        for(auto n : ULbox.accessPattern.material_conservativeVariables[m])
                        {
                            fluxbox(i,j,k,n)	= flux(n,UR,d);
                        }

                        for(int row=0;row<ULbox.numberOfComponents;row++)
                        {
                            for(int col=0;col<ULbox.numberOfComponents;col++)
                            {
                                fluxbox(i,j,k,VSTAR,m,row,col) = UR(V_TENSOR,m,row,col);
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
    }

    return;
}

/** Calculates the HLLC fluxes.
 */
void calc_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, ParameterStruct& parameters, Direction_enum d, const Real *dx, const Real *prob_lo, BoxAccessLevelSet& LS)
{
    const auto lo = lbound(ULbox.box);
    const auto hi = ubound(ULbox.box);

    Vector<Real> SR(parameters.numberOfMaterials);
    Vector<Real> SL(parameters.numberOfMaterials);
    Vector<Real> Sstar(parameters.numberOfMaterials);

    int extra[3] = {0,0,0};

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
                    if(n.var == V_TENSOR)
                    {
                        U1box(i,j,k,n) += (dt/dx[d])*(fluxbox(i,j,k,n) - fluxbox.right(d,i,j,k,n) -(2.0/3.0)*Ubox(i,j,k,n)*(fluxbox(i,j,k,USTAR,n.mat) - fluxbox.right(d,i,j,k,USTAR,n.mat)) + Ubox(i,j,k,VELOCITY,n.mat,n.row)*(fluxbox(i,j,k,VSTAR,n.mat,d,n.col) - fluxbox.right(d,i,j,k,VSTAR,n.mat,d,n.col)));
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
void MUSCLextrapolate(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, Direction_enum d)
{
    Real r;

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    Real grad;

    for    			(auto n : U.accessPattern.conservativeVariables)
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

                    grad = vanLeerlimiter(r)*0.5*(U.right(d,i,j,k,n)-U.left(d,i,j,k,n));


                    UL(i,j,k,n) = U(i,j,k,n) - 0.5*grad;
                    UR(i,j,k,n) = U(i,j,k,n) + 0.5*grad;

                }
            }
        }
    }

    UL.checkLimits(U.accessPattern.primitiveVariables);
    UR.checkLimits(U.accessPattern.primitiveVariables);

    return;
}

void halfTimeStepEvolution(BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& ULboxnew, BoxAccessCellArray& URboxnew, Direction_enum d, ParameterStruct& parameters, const Real* dx, Real dt)
{
    const auto lo = lbound(ULbox.box);
    const auto hi = ubound(ULbox.box);

    Material_type phase = (parameters.SOLID == 1 ? solid : fluid);

    Real fluxL, fluxR;

    for    			(auto n : ULbox.accessPattern.conservativeVariables)
    {
        for 		   (int k = lo.z; k <= hi.z; ++k)
        {
            for 	   (int j = lo.y; j <= hi.y; ++j)
            {
                for    (int i = lo.x; i <= hi.x; ++i)
                {
                    ULbox.conservativeToPrimitive(i,j,k);
                    URbox.conservativeToPrimitive(i,j,k);

                    Cell UL(ULbox,i,j,k,phase);
                    Cell UR(URbox,i,j,k,phase);

                    fluxL = flux(n,UL,d);
                    fluxR = flux(n,UR,d);

                    ULboxnew(i,j,k,n) = ULbox(i,j,k,n) + 0.5*(dt/dx[d])*(fluxL - fluxR);
                    URboxnew(i,j,k,n) = URbox(i,j,k,n) + 0.5*(dt/dx[d])*(fluxL - fluxR);
                }
            }
        }
    }
}
