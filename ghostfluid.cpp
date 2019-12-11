#include "simulationheader.h"

void getStarState(Cell& U, Cell& UStar, Real SK, Real Sstar, ParameterStruct& parameters, Direction_enum d, int m);

void fastSweep(BoxAccessCellArray& U, BoxAccessLevelSet& LS, const Real* dx, int xsense, int ysense, int sign)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    Real phix,phiy;
    Real phinew;
    Real discriminant;

    Real nx,ny;

    int xlo,xhi,ylo,yhi;

    if(xsense > 0)
    {
        xlo = lo.x;
        xhi = hi.x;
    }
    else
    {
        xlo = hi.x;
        xhi = lo.x;
    }

    if(ysense > 0)
    {
        ylo = lo.y;
        yhi = hi.y;
    }
    else
    {
        ylo = hi.y;
        yhi = lo.y;
    }

    Real dx2 = dx[0]*dx[0];
    Real dy2 = dx[1]*dx[1];

    int m;

    if(sign > 0)
    {
        m = 0;
    }
    else
    {
        m = 1;
    }

    int xsign,ysign;


    for             (auto n : U.accessPattern.material_primitiveVariables[m])
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = ylo; LS.customComparator(j,yhi,ysense) ; LS.customChanger(j,ysense))
            {
                for (int i = xlo; LS.customComparator(i,xhi,xsense) ; LS.customChanger(i,xsense))
                {
                    if( (sgn<Real,int>(LS(i,j,k,0)) != sgn<int,int>(sign)) && !LS.cellIsNextToAnInterface(i,j,k,0)  &&  LS.cellIsNearInterface(i,j,k,dx))
                    {
                        phix = ( std::abs(LS(i-1,j  ,k,0)) < std::abs(LS(i+1,j  ,k,0)) ? U(i-1,j  ,k,n) : U(i+1,j  ,k,n));
                        phiy = ( std::abs(LS(i  ,j-1,k,0)) < std::abs(LS(i  ,j+1,k,0)) ? U(i  ,j-1,k,n) : U(i  ,j+1,k,n));

                        xsign = ( std::abs(LS(i-1,j  ,k,0)) < std::abs(LS(i+1,j  ,k,0)) ? -1 : 1 );
                        ysign = ( std::abs(LS(i  ,j-1,k,0)) < std::abs(LS(i  ,j+1,k,0)) ? -1 : 1 );

                        U(i,j,k,n) = (std::abs(LS(i+xsign,j,k,0)) < std::abs(LS(i,j+ysign,k,0)) ? U(i+xsign,j,k,n) : U(i,j+ysign,k,n));


                        //U(i,j,k,n) = std::min(phix,phiy);

                       /* //discriminant = (phix/dx2+phiy/dy2)*(phix/dx2+phiy/dy2)-(1.0/dx2+1.0/dy2)*(phix*phix/dx2+phiy*phiy/dy2);

                        //Print() << discriminant << " " << -1.0/(dx2*dy2)*(phix-phiy)*(phix-phiy) << std::endl;

                        discriminant = -1.0/(dx2*dy2)*(phix-phiy)*(phix-phiy);



                        if(discriminant < 0.0)
                        {
                            if(std::abs(phix) < std::abs(phiy))
                            {
                                phinew = phix;
                            }
                            else
                            {
                                discriminant = 0.0;

                                phinew = phiy;
                            }
                        }
                        else
                        {
                            phinew = ((phix/dx2+phiy/dy2)+sqrt(discriminant))/(1.0/dx2+1.0/dy2);
                        }

                        //if(phinew < U(i,j,k,n))
                        //{
                            U(i,j,k,n) = phinew;
                        //}*/

                    }
                }
            }
        }
    }

    return;
}

void sweep(MultiFab& S_new, CellArray& U, LevelSet& LS, const Real* dx, Geometry& geom, Vector<BCRec>& bc, int dir, int sense, int sign)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.validbox();

        BoxAccessCellArray baca(mfi,box,U);
        BoxAccessLevelSet  bals(mfi,box,LS);

        fastSweep(baca,bals,dx,dir,sense,sign);
    }

    U.data.FillBoundary(geom.periodicity());
    FillDomainBoundary(U.data, geom, bc);
}

Real getSstarFromDifferentMaterials(Cell& UL, Cell& UR, Real SL, Real SR, int i, int j, int k, Direction_enum d, Vector<int>& m)
{
    if( std::abs(UL(RHO,m[0])*(SL-UL(VELOCITY,m[0],d)) -  UR(RHO,m[1])*(SR-UR(VELOCITY,m[1],d))) < 1E-10 )
    {
        return 0.0;
    }
    else
    {
        return (-UR(SIGMA,m[1],d,d)+UL(SIGMA,m[0],d,d)+UL(RHOU,m[0],d)*(SL-UL(VELOCITY,m[0],d))-UR(RHOU,m[1],d)*(SR-UR(VELOCITY,m[1],d)))/(UL(RHO,m[0])*(SL-UL(VELOCITY,m[0],d)) -  UR(RHO,m[1])*(SR-UR(VELOCITY,m[1],d)));
    }
}

void originalGFM(Cell& U, Cell& UL, Cell& UR, Vector<int>& m)
{
    U(P,m[0]) = UR(P,m[1]);
    U(P,m[1]) = UL(P,m[0]);

    U(VELOCITY,m[0],x) = UR(VELOCITY,m[1],x);
    U(VELOCITY,m[1],x) = UL(VELOCITY,m[0],x);
    U(VELOCITY,m[0],y) = UR(VELOCITY,m[1],y);
    U(VELOCITY,m[1],y) = UL(VELOCITY,m[0],y);
    U(VELOCITY,m[0],z) = UR(VELOCITY,m[1],z);
    U(VELOCITY,m[1],z) = UL(VELOCITY,m[0],z);

    U(RHO,m[0]) = UL(RHO,m[0])*std::pow(U(P,m[0])/UL(P,m[0]), 1.0/1.4);
    U(RHO,m[1]) = UR(RHO,m[1])*std::pow(U(P,m[1])/UR(P,m[1]), 1.0/1.4);

    for(int row = 0; row < U.numberOfComponents; row++)
    {
        for(int col = 0; col < U.numberOfComponents; col++)
        {
            if(U.accessPattern.materialInfo[m[0]].phase == solid)
            {
                U(V_TENSOR,m[0],row,col) = UL(V_TENSOR,m[0],row,col);
            }
            if(U.accessPattern.materialInfo[m[1]].phase == solid)
            {
                U(V_TENSOR,m[1],row,col) = UR(V_TENSOR,m[1],row,col);
            }
        }
    }
}

void getBothStarStatesAndSetNewValue(int i, int j, int k, BoxAccessCellArray& Ubox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, Vector<int>& probe_m, ParameterStruct& parameters)
{
    Cell U      (Ubox,i,j,k,fluid);
    Cell UL     (ULbox,i,j,k,fluid);
    Cell UR     (URbox,i,j,k,fluid);
    Cell UStar  (UStarbox,i,j,k,fluid);


    Real SR    = std::max(std::abs(UL(VELOCITY,probe_m[0],x))+UL(SOUNDSPEED,probe_m[0]),std::abs(UR(VELOCITY,probe_m[1],x))+UR(SOUNDSPEED,probe_m[1]));
    Real SL    = -SR;

    Real Sstar = getSstarFromDifferentMaterials(UL,UR,SL,SR,i,j,k,x,probe_m);

    getStarState(UL,UStar,SL,Sstar,parameters,x,probe_m[0]);
    getStarState(UR,UStar,SR,Sstar,parameters,x,probe_m[1]);

    U = UStar;

    Ubox.conservativeToPrimitive(i,j,k);

    //originalGFM(U,UL,UR,probe_m);


}

void boxGhostFluidValues(BoxAccessCellArray& U, BoxAccessCellArray& U1, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& Ustar, BoxAccessLevelSet& LS0, const Real* dx, const Real* prob_lo, ParameterStruct& parameters)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    Real nx,ny;
    Real interface_x, interface_y;
    Real current_x,   current_y,  current_z;

    Vector<Real> probe_x(2);
    Vector<Real> probe_y(2);
    Vector<int>  probe_m(2);

    Vector<int>  corner_x(2);
    Vector<int>  corner_y(2);
    Vector<Real> xdiff(2);
    Vector<Real> ydiff(2);
    Vector<Real> B(4);

    Vector<int> probe_int(3);

    int m;

    for 		(int k = lo.z; k <= hi.z; ++k)
    {
                current_z = prob_lo[2] + (Real(k)+0.5)*dx[2];

        for 	(int j = lo.y; j <= hi.y; ++j)
        {
                current_y = prob_lo[1] + (Real(j)+0.5)*dx[1];

            for (int i = lo.x; i <= hi.x; ++i)
            {

                current_x = prob_lo[0] + (Real(i)+0.5)*dx[0];

                if(LS0.cellIsNextToAnInterface(i,j,k,0))
                {
                    LS0.calculateNormal(i,j,k,0,dx,nx,ny);
                    LS0.calculateInterpolationPoint(i,j,k,0,dx,nx,ny,current_x,current_y,interface_x,interface_y);
                    LS0.calculateProbes(i,j,k,0,dx,nx,ny,interface_x,interface_y,probe_x,probe_y);

                    m = LS0.whatMaterialIsValid(i,j,k);

                    U.realPositionToCell(probe_int[0],probe_int[1],probe_int[0],probe_x[0],probe_y[0],probe_x[0],dx,prob_lo);
                    probe_m[0] = LS0.whatMaterialIsValid(probe_int[0],probe_int[1],k);

                    U.realPositionToCell(probe_int[0],probe_int[1],probe_int[0],probe_x[1],probe_y[1],probe_x[1],dx,prob_lo);
                    probe_m[1] = LS0.whatMaterialIsValid(probe_int[0],probe_int[1],k);

                    if(probe_m[0] == probe_m[1])
                    {
                        if((current_x-probe_x[0])*(current_x-probe_x[0])+(current_y-probe_y[0])*(current_y-probe_y[0]) <  (current_x-probe_x[1])*(current_x-probe_x[1])+(current_y-probe_y[1])*(current_y-probe_y[1]) )
                        {
                            probe_m[0] = m;
                            probe_m[1] = ( m==0 ? 1 : 0);
                        }
                        else
                        {
                            probe_m[1] = m;
                            probe_m[0] = ( m==0 ? 1 : 0);
                        }

                        if(probe_m[0] == probe_m[1])
                        {
                            Abort("Error in probe extrapolation: materials are the same");
                        }
                    }

                    UL.bilinearInterpolation(U,i,j,k,dx,prob_lo,probe_x[0],probe_y[0],probe_m[0], corner_x, corner_y, xdiff, ydiff, B);
                    UR.bilinearInterpolation(U,i,j,k,dx,prob_lo,probe_x[1],probe_y[1],probe_m[1], corner_x, corner_y, xdiff, ydiff, B);

                    UL.rotateFrameSoXPointsAlongNormal(i,j,k,nx,ny,probe_m[0]);
                    UR.rotateFrameSoXPointsAlongNormal(i,j,k,nx,ny,probe_m[1]);

                    UL.primitiveToConservative(i,j,k);
                    UR.primitiveToConservative(i,j,k);

                    UL.getSoundSpeed(i,j,k);
                    UR.getSoundSpeed(i,j,k);

                    getBothStarStatesAndSetNewValue(i,j,k,U1,UL,UR,Ustar,probe_m,parameters);

                    U1.rotateFrameBack(i,j,k,nx,ny);

                    U1.primitiveToConservative(i,j,k);
                }
                else
                {
                    if(LS0.cellIsNearInterface(i,j,k,dx))
                    {
                        if(LS0(i,j,k,0) > 0.0)
                        {
                            for(auto n : U.accessPattern.material_primitiveVariables[1])
                            {
                                U1(i,j,k,n) = 1E20;
                            }
                        }
                        else
                        {
                            for(auto n : U.accessPattern.material_primitiveVariables[0])
                            {
                                U1(i,j,k,n) = 1E20;
                            }
                        }
                    }
                }
            }
        }
    }

    return;
}

void setGhostFluidValues(MultiFab& S_new,CellArray& U, CellArray& U1, CellArray& UL, CellArray& UR, CellArray& Ustar, LevelSet& LS0, const Real* dx, const Real* prob_lo, ParameterStruct& parameters, Geometry& geom, Vector<BCRec>& bc)
{


    U1 = U;


#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray Ubox (mfi,bx,U );
        BoxAccessCellArray U1box(mfi,bx,U1);
        BoxAccessCellArray ULbox(mfi,bx,UL);
        BoxAccessCellArray URbox(mfi,bx,UR);
        BoxAccessCellArray Ustarbox(mfi,bx,Ustar);
        BoxAccessLevelSet  bals (mfi,bx,LS0);

        boxGhostFluidValues(Ubox,U1box,ULbox,URbox,Ustarbox,bals,dx,prob_lo,parameters);
    }

    U1.data.FillBoundary(geom.periodicity());
    FillDomainBoundary(U1.data, geom, bc);


    int forward   =  1;
    int backward  = -1;

    int positive  =  1;
    int negative  = -1;

    for(int it = 0; it < 10; it++)
    {
        int sweepingDone = 1;

        sweep(S_new, U1, LS0, dx,geom,bc,forward,  forward,   positive);
        sweep(S_new, U1, LS0, dx,geom,bc,backward, forward,   positive);
        sweep(S_new, U1, LS0, dx,geom,bc,backward, backward,  positive);
        sweep(S_new, U1, LS0, dx,geom,bc,forward,  backward,  positive);


        sweep(S_new, U1, LS0, dx,geom,bc,forward,  forward,   negative);
        sweep(S_new, U1, LS0, dx,geom,bc,backward, forward,   negative);
        sweep(S_new, U1, LS0, dx,geom,bc,backward, backward,  negative);
        sweep(S_new, U1, LS0, dx,geom,bc,forward,  backward,  negative);

        for(auto n : U1.accessPattern.primitiveVariables)
        {
            if( U1.data.max(U.getArrayPosition(n)) > 1E19)
            {
                sweepingDone *= 0;
            }
        }

        //if(it == 0)
         //   break;

        //break;

        if(sweepingDone)
        {
            break;
        }
    }

    U1.primitiveToConservative();

    return;
}


