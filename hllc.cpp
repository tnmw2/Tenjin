#include "simulationheader.h"

/** Calculates the contact wave speed estimate.
 */
double getSstar(Cell& UL, Cell& UR, Real SL, Real SR, int i, int j, int k, Direction_enum d)
{
    if( std::abs(UL(RHO)*(SL-UL(VELOCITY,0,d)) -  UR(RHO)*(SR-UR(VELOCITY,0,d))) < 1E-10 )
    {
        return 0.0;
    }
    else
    {
        return (-UR(SIGMA,0,d,d)+UL(SIGMA,0,d,d)+UL(RHOU,0,d)*(SL-UL(VELOCITY,0,d))-UR(RHOU,0,d)*(SR-UR(VELOCITY,0,d)))/(UL(RHO)*(SL-UL(VELOCITY,0,d)) -  UR(RHO)*(SR-UR(VELOCITY,0,d)));
    }
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


                SR = std::max(std::abs(UL(VELOCITY,0,d))+UL(SOUNDSPEED),std::abs(UR(VELOCITY,0,d))+UR(SOUNDSPEED));
                SL = -SR;

                Sstar = getSstar(UL,UR,SL,SR,i,j,k,d);

                if(std::isnan(SL) || std::isnan(SR) || std::isnan(Sstar))
                {
                    //UL.parent->checking = 1;
                    //UR.parent->checking = 1;

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
                        Vector<Real> err;

                        err.push_back(SL);
                        err.push_back(SR);

                        std::string message = "Nan in SL SR wavespeeds before flux: ";

                        customAbort(err,message);
                    }

                    if(std::isnan(Sstar))
                    {

                        Vector<Real> err;

                        err.push_back(prob_lo[0] + (Real(i)+0.5)*dx[0]);
                        err.push_back(prob_lo[1] + (Real(j)+0.5)*dx[1]);
                        err.push_back(Sstar);
                        err.push_back(UL(SIGMA,0,d,d));
                        err.push_back(UR(SIGMA,0,d,d));
                        err.push_back(UL(P));
                        err.push_back(UR(P));
                        err.push_back(UL(V_TENSOR,0,d,d));
                        err.push_back(UR(V_TENSOR,0,d,d));
                        err.push_back(UL(RHOU,0,d));
                        err.push_back(UR(RHOU,0,d));
                        err.push_back(UL(VELOCITY,0,d));
                        err.push_back(UR(VELOCITY,0,d));
                        err.push_back(UL(RHO));
                        err.push_back(UR(RHO));
                        err.push_back(SL);
                        err.push_back(SR);

                        err.push_back(UL.parent->operator()(UL.parent_i,UL.parent_j,k,P));

                        std::string message = "Nan in Sstar wavespeed before flux at: ";

                        customAbort(err,message);

                    }
                }

                //fluxbox(i,j,k,USTAR) = Sstar;

                if(SL>=0.0)
                {
                    //getStarState(UL,ULStar,SL,Sstar,parameters,d);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        fluxbox(i,j,k,n)	= flux(n,UL,d);

                        /*if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,ULStar,d);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL,d);
                        }*/
                    }


                    /*for(int row=0;row<ULbox.numberOfComponents;row++)
                    {
                        for(int col=0;col<ULbox.numberOfComponents;col++)
                        {
                            fluxbox(i,j,k,VSTAR,0,row,col) = UL(V_TENSOR,0,row,col);
                        }
                    }*/

                    UStarStar = UL;

                }
                else if(Sstar>=0.0)
                {
                    getStarState(UL,ULStar,SL,Sstar,parameters,d);

                    SLT = Sstar - ULStar.parent->transverseWaveSpeed(ULStar.parent_i,ULStar.parent_j,ULStar.parent_k);

                    if(SLT>=0.0)
                    {
                        for(auto n : ULbox.accessPattern.conservativeVariables)
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL,d)+SL*(ULStar(n)-UL(n));

                            /*if(n.var == ALPHA)
                            {
                                 fluxbox(i,j,k,n)	= flux(n,ULStar,d);
                            }
                            else
                            {
                                fluxbox(i,j,k,n)	= flux(n,UL,d)+SL*(ULStar(n)-UL(n));
                            }*/
                        }


                        /*for(int row=0;row<ULbox.numberOfComponents;row++)
                        {
                            for(int col=0;col<ULbox.numberOfComponents;col++)
                            {
                                fluxbox(i,j,k,VSTAR,0,row,col) = ULStar(V_TENSOR,0,row,col);
                            }
                        }*/

                        UStarStar = ULStar;

                    }
                    else
                    {
                        getStarState(UR,URStar,SR,Sstar,parameters,d);

                        SRT = Sstar + URStar.parent->transverseWaveSpeed(URStar.parent_i,URStar.parent_j,URStar.parent_k);

                        getStarStarState(UL,UR,ULStar,URStar,UStarStar,SL,SR,SLT,SRT,Sstar,parameters,d,UL,ULStar,SLT);

                        for(auto n : ULbox.accessPattern.conservativeVariables)
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL,d)+SL*(ULStar(n)-UL(n))+SLT*(UStarStar(n)-ULStar(n));

                            /*if(n.var == ALPHA)
                            {
                                 fluxbox(i,j,k,n)	= flux(n,ULStar,d);
                            }
                            else
                            {
                                fluxbox(i,j,k,n)	= flux(n,UL,d)+SL*(ULStar(n)-UL(n))+SLT*(UStarStar(n)-ULStar(n));
                            }*/
                        }


                        /*for(int row=0;row<ULbox.numberOfComponents;row++)
                        {
                            for(int col=0;col<ULbox.numberOfComponents;col++)
                            {
                                fluxbox(i,j,k,VSTAR,0,row,col) = UStarStar(V_TENSOR,0,row,col);
                            }
                        }*/
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
                            fluxbox(i,j,k,n)	= flux(n,UR,d)+SR*(URStar(n)-UR(n))+SRT*(UStarStar(n)-URStar(n));

                            /*if(n.var == ALPHA)
                            {
                                 fluxbox(i,j,k,n)	= flux(n,URStar,d);
                            }
                            else
                            {
                                fluxbox(i,j,k,n)	= flux(n,UR,d)+SR*(URStar(n)-UR(n))+SRT*(UStarStar(n)-URStar(n));
                            }*/
                        }

                        /*for(int row=0;row<ULbox.numberOfComponents;row++)
                        {
                            for(int col=0;col<ULbox.numberOfComponents;col++)
                            {
                                fluxbox(i,j,k,VSTAR,0,row,col) = UStarStar(V_TENSOR,0,row,col);
                            }
                        }*/
                    }
                    else
                    {
                        for(auto n : ULbox.accessPattern.conservativeVariables)
                        {
                            fluxbox(i,j,k,n)	= flux(n,UR,d)+SR*(URStar(n)-UR(n));

                            /*if(n.var == ALPHA)
                            {
                                 fluxbox(i,j,k,n)	= flux(n,URStar,d);
                            }
                            else
                            {
                                fluxbox(i,j,k,n)	= flux(n,UR,d)+SR*(URStar(n)-UR(n));
                            }*/
                        }

                        /*for(int row=0;row<ULbox.numberOfComponents;row++)
                        {
                            for(int col=0;col<ULbox.numberOfComponents;col++)
                            {
                                fluxbox(i,j,k,VSTAR,0,row,col) = URStar(V_TENSOR,0,row,col);
                            }
                        }*/

                        UStarStar = URStar;
                    }
                }
                else
                {
                    getStarState(UR,URStar,SR,Sstar,parameters,d);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        fluxbox(i,j,k,n)	= flux(n,UR,d);

                        /*if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,URStar,d);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UR,d);
                        }*/
                    }

                    /*for(int row=0;row<ULbox.numberOfComponents;row++)
                    {
                        for(int col=0;col<ULbox.numberOfComponents;col++)
                        {
                            fluxbox(i,j,k,VSTAR,0,row,col) = UR(V_TENSOR,0,row,col);
                        }
                    }*/

                    UStarStar = UR;

                }

                UStarStarbox.conservativeToPrimitive(i,j,k);

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
void calc_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, ParameterStruct& parameters, Direction_enum d, const Real *dx, const Real *prob_lo)
{
    const auto lo = lbound(ULbox.box);
    const auto hi = ubound(ULbox.box);

    Real SR;
    Real SL;
    Real Sstar;

    int extra[3] = {0,0,0};

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

                //fluxbox(i,j,k,USTAR) = Sstar;

                if(SL>=0.0)
                {
                    //getStarState(UL,UStar,SL,Sstar,parameters,d);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        fluxbox(i,j,k,n)	= flux(n,UL,d);

                        /*if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,UStar,d);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL,d);
                        }*/
                    }

                    UStar = UL;
                }
                else if(Sstar>=0.0)
                {
                    getStarState(UL,UStar,SL,Sstar,parameters,d);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        fluxbox(i,j,k,n)	= flux(n,UL,d)+SL*(UStar(n)-UL(n));

                        /*if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,UStar,d);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UL,d)+SL*(UStar(n)-UL(n));
                        }*/
                    }

                }
                else if(SR>=0.0)
                {
                    getStarState(UR,UStar,SR,Sstar,parameters,d);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        fluxbox(i,j,k,n)	= flux(n,UR,d)+SR*(UStar(n)-UR(n));

                        /*if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,UStar,d);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UR,d)+SR*(UStar(n)-UR(n));
                        }*/
                    }
                }
                else
                {
                    //getStarState(UR,UStar,SR,Sstar,parameters,d);

                    for(auto n : ULbox.accessPattern.conservativeVariables)
                    {
                        fluxbox(i,j,k,n)	= flux(n,UR,d);

                        /*if(n.var == ALPHA)
                        {
                             fluxbox(i,j,k,n)	= flux(n,UStar,d);
                        }
                        else
                        {
                            fluxbox(i,j,k,n)	= flux(n,UR,d);
                        }*/
                    }

                    UStar = UR;
                }

                //UStarbox.conservativeToPrimitive(i,j,k);
                //UStarbox.getSoundSpeed(i,j,k);
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
                    if(n.var == ALPHA || n.var == V_TENSOR)
                    {
                        continue;
                        //U1box(i,j,k,n) += (dt/dx[d])*(fluxbox(i,j,k,n) - fluxbox.right(d,i,j,k,n) -Ubox(i,j,k,n)*(fluxbox(i,j,k,USTAR) - fluxbox.right(d,i,j,k,USTAR)));
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
void MUSCLextrapolate(const BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, Direction_enum d)
{
    Real r,grad;

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    int extra[3] = {0,0,0};

    extra[d]=1;

    for    			(auto n : U.accessPattern.primitiveVariables)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    r = (U(i,j,k,n)-U(i-extra[0],j-extra[1],k-extra[2],n))/(U(i+extra[0],j+extra[1],k+extra[2],n)-U(i,j,k,n));

                    if(std::isinf(r) || std::isnan(r))
                    {
                        r = 0.0;
                    }

                    grad = vanLeerlimiter(r)*0.5*(U(i+extra[0],j+extra[1],k+extra[2],n)-U(i-extra[0],j-extra[1],k-extra[2],n));


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




