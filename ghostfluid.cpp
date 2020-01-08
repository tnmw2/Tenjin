#include "simulationheader.h"

Real fK(const Real p, Cell& W, int m)
{
    Real adiabaticIndex = W.accessPattern.materialInfo[m].EOS->adiabaticIndex;

    Real ans;

    Real A = 2.0/((adiabaticIndex+1.0)*W(RHO,m));

    Real B = (adiabaticIndex-1.0)/(adiabaticIndex+1.0)*W(P,m);

    if(p > W(P,m))
    { //shock
        ans=(p-W(P,m))*std::sqrt(A/(p+B));
    }
    else
    { //rarefaction
        ans=2.0*W(SOUNDSPEED,m)/(adiabaticIndex-1.0)*(std::pow(p/(W(P,m)),(adiabaticIndex-1.0)/(2.0*adiabaticIndex))-1.0);
    }

    return ans;
}

Real f (const Real p, Cell& WL, Cell& WR, Real du, Vector<int>& m)
{
    Real ans;

    ans = fK(p,WL,m[0])+fK(p,WR,m[1])+du;

    return ans;

}

void exactSolver_bisection(Cell& UL, Cell& UR, ParameterStruct& parameters, Cell& UStar, Vector<int>& m)
{
    Real du = UR(VELOCITY,m[1],x)-UL(VELOCITY,m[0],x);
    Real p;
    Real pA = 0.0; // std::min(UL(P,m[0]),UR(P,m[1]));
    Real pB = 100.0*std::max(UL(P,m[0]),UR(P,m[1]));

    Real tolerance = 1E-10;

    do
    {
        p=0.5*(pA+pB);


        if(sgn<Real,int>(f(pB,UL,UR,du,m)) == sgn<Real,int>(f(pA,UL,UR,du,m)))
        {
            Abort("Sign error in rGFM exact solver");
        }

        if(sgn<Real,int>(f(p,UL,UR,du,m)) == sgn<Real,int>(f(pA,UL,UR,du,m)))
        {
            pA=p;
        }
        else
        {
            pB=p;
        }

    }
    while(std::abs(pA-pB)/p > tolerance);

    UStar(P,m[0]) = p;
    UStar(P,m[1]) = p;

    UStar(VELOCITY,m[0],x) = 0.5*(UL(VELOCITY,m[0],x)+UR(VELOCITY,m[1],x)+fK(p,UR,m[1])-fK(p,UL,m[0]));
    UStar(VELOCITY,m[1],x) = UStar(VELOCITY,m[0],x);

    UStar(VELOCITY,m[0],y) = UL(VELOCITY,m[0],y);
    UStar(VELOCITY,m[1],y) = UR(VELOCITY,m[1],y);

    UStar(VELOCITY,m[0],z) = UL(VELOCITY,m[0],z);
    UStar(VELOCITY,m[1],z) = UR(VELOCITY,m[1],z);

    return;
}

void exactSolver_FindIntermediateState(Cell& UL, Cell&  UR, ParameterStruct& parameters, Cell& UStar, Vector<int>& m )
{
    Real adiabaticIndexL = UL.accessPattern.materialInfo[m[0]].EOS->adiabaticIndex;
    Real adiabaticIndexR = UR.accessPattern.materialInfo[m[1]].EOS->adiabaticIndex;


    exactSolver_bisection(UL,UR,parameters,UStar,m);


    if(UL(P,m[0]) < UStar(P,m[0]))
    {
        //shock -> uses Rankine Hugoniot
        UStar(RHO,m[0]) = UL(RHO,m[0])*(((adiabaticIndexL-1.0)/(adiabaticIndexL+1.0)+UStar(P,m[0])/UL(P,m[0])) / ((adiabaticIndexL-1.0)/(adiabaticIndexL+1.0)*UStar(P,m[0])/UL(P,m[0])+1.0));
    }
    else
    {
        //rarefaction -> uses constant entropy extrapolation
        UStar(RHO,m[0]) = UL(RHO,m[0])*std::pow((UStar(P,m[0]))/(UL(P,m[0])),1.0/adiabaticIndexL);
    }
    if(UR(P,m[1]) < UStar(P,m[1]))
    {
        //shock -> uses Rankine Hugoniot
        UStar(RHO,m[1]) = UR(RHO,m[1])*(((adiabaticIndexR-1.0)/(adiabaticIndexR+1.0)+UStar(P,m[1])/UR(P,m[1])) / ((adiabaticIndexR-1.0)/(adiabaticIndexR+1.0)*UStar(P,m[1])/UR(P,m[1])+1.0));
    }
    else
    {
        //rarefaction -> uses constant entropy extrapolation
        UStar(RHO,m[1]) = UR(RHO,m[1])*std::pow((UStar(P,m[1]))/(UR(P,m[1])),1.0/adiabaticIndexR);
    }

    return;

}

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
    //FillDomainBoundary(U.data, geom, bc);
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

void getSigmaStarFromDifferentMaterials(Cell& UKStar, Real Sstar, Real SLT, Real SRT, Cell& UL, Cell& UR, Cell& ULStar, Cell& URStar, Direction_enum d, ParameterStruct& parameters, int matL, int matR, int matK)
{
    /*for(int i=0;i<UL.numberOfComponents;i++)
    {
        UKStar(SIGMA,matK,i,d) = (ULStar(RHO,matL)*(Sstar-SLT)*URStar(RHO,matR)*(Sstar-SRT)*(UL(VELOCITY,matL,i) -UR(VELOCITY,matR,i)) + ULStar(RHO,matL)*(Sstar-SLT)*UR(SIGMA,matR,i,d) - URStar(RHO,matR)*(Sstar-SRT)*UL(SIGMA,matL,i,d))/(ULStar(RHO,matL)*(Sstar-SLT)-URStar(RHO,matR)*(Sstar-SRT));
    }*/

    for(int i=0;i<UL.numberOfComponents;i++)
    {
        if(i == d)
        {
            if(matK == matL)
            {
                UKStar(SIGMA,matK,i,d) = UL(SIGMA,matK,i,d);
            }
            else
            {
                UKStar(SIGMA,matK,i,d) = UR(SIGMA,matK,i,d);
            }
        }
        else
        {
            UKStar(SIGMA,matK,i,d) = 0.0;
        }
    }


    return;
}

void getStarStarStateFromDifferentMaterials(Cell& UL, Cell& UR, Cell& UStar, Cell& UStarStar, Real SL, Real SR, Real SLT, Real SRT, Real Sstar, ParameterStruct& parameters, Direction_enum d, Cell& UK, Cell& UKStar, Real SKT, int matL, int matR, int matK)
{
    int N = UL.numberOfComponents;

    //UStarStar = UKStar;

    UStarStar.setMaterial(UStar,matK);

    if(UL.accessPattern.materialInfo[matK].phase == solid)
    {
        getSigmaStarFromDifferentMaterials(UKStar,Sstar,SLT,SRT,UL,UR,UStar,UStar,d,parameters,matL,matR,matK);

        for(int row =0;row<N;row++)
        {
            if(row == d)
            {
                continue;
            }
            else
            {
                UStarStar(RHOU,matK,row) += (UKStar(SIGMA,matK,row,d)-UK(SIGMA,matK,row,d))/(Sstar-SKT);
            }

            UStarStar(VELOCITY,matK,row) = UStarStar(RHOU,matK,row)/UStarStar(RHO,matK);
        }


        for(int row =0;row<N;row++)
        {
            if(row == d)
            {
                continue;
            }
            else
            {
                UStarStar(TOTAL_E,matK) += (UStarStar(VELOCITY,matK,row)*UKStar(SIGMA,matK,row,d)-UK(VELOCITY,matK,row)*UK(SIGMA,matK,row,d))/(Sstar-SKT);

                for(int col =0;col<N;col++)
                {
                    UStarStar(V_TENSOR,matK,row,col) += (UStarStar(V_TENSOR,matK,d,col)*(UStarStar(VELOCITY,matK,row)-UK(VELOCITY,matK,row)))/(Sstar-SKT);
                }
            }
        }
    }

    return;
}

void originalGFM(Cell& Uold, Cell& U, Cell& UL, Cell& UR, Vector<int>& m, int validMat)
{
    /*U = Uold;

    if(validMat == m[0])
    {
        U(P,m[1]) = UL(P,m[0]);

        U(VELOCITY,m[1],x) = UL(VELOCITY,m[0],x);
        U(VELOCITY,m[1],y) = UL(VELOCITY,m[0],y);
        U(VELOCITY,m[1],z) = UL(VELOCITY,m[0],z);

        U(RHO,m[1]) = UR(RHO,m[1])*std::pow(U(P,m[1])/UR(P,m[1]), 1.0/1.4);
    }
    else
    {
        U(P,m[0]) = UR(P,m[1]);

        U(VELOCITY,m[0],x) = UR(VELOCITY,m[1],x);
        U(VELOCITY,m[0],y) = UR(VELOCITY,m[1],y);
        U(VELOCITY,m[0],z) = UR(VELOCITY,m[1],z);

        U(RHO,m[0]) = UL(RHO,m[0])*std::pow(U(P,m[0])/UL(P,m[0]), 1.0/1.4);
    }*/

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
}

void getBothStarStatesAndSetNewValue_5Wave(int i, int j, int k, BoxAccessCellArray& Ubox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, BoxAccessCellArray& UStarStarbox, Vector<int>& probe_m, ParameterStruct& parameters)
{

    Cell U          (Ubox,        i,j,k,solid);
    Cell UL         (ULbox,       i,j,k,solid);
    Cell UR         (URbox,       i,j,k,solid);
    Cell UStar      (UStarbox,    i,j,k,solid);
    Cell UStarStar  (UStarStarbox,i,j,k,solid);


    Real SR    = std::max(std::abs(UL(VELOCITY,probe_m[0],x))+UL(SOUNDSPEED,probe_m[0]),std::abs(UR(VELOCITY,probe_m[1],x))+UR(SOUNDSPEED,probe_m[1]));
    Real SL    = -SR;

    Real Sstar = getSstarFromDifferentMaterials(UL,UR,SL,SR,i,j,k,x,probe_m);

    getStarState(UL,UStar,SL,Sstar,parameters,x,probe_m[0]);
    getStarState(UR,UStar,SR,Sstar,parameters,x,probe_m[1]);

    Real SLT = Sstar - UStar.parent->transverseWaveSpeed(UStar.parent_i,UStar.parent_j,UStar.parent_k,probe_m[0]);
    Real SRT = Sstar + UStar.parent->transverseWaveSpeed(UStar.parent_i,UStar.parent_j,UStar.parent_k,probe_m[1]);

    getStarStarStateFromDifferentMaterials(UL,UR,UStar,UStarStar,SL,SR,SLT,SRT,Sstar,parameters,x,UL,UStar,SLT,probe_m[0],probe_m[1],probe_m[0]);
    getStarStarStateFromDifferentMaterials(UL,UR,UStar,UStarStar,SL,SR,SLT,SRT,Sstar,parameters,x,UR,UStar,SRT,probe_m[0],probe_m[1],probe_m[1]);

    U = UStarStar;

    Ubox.conservativeToPrimitive(i,j,k);

}

void getBothStarStatesAndSetNewValue_3Wave(int i, int j, int k, BoxAccessCellArray& Uoldbox, BoxAccessCellArray& Ubox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, BoxAccessCellArray& UStarStarbox, Vector<int>& probe_m, ParameterStruct& parameters, int validMat)
{

    Cell Uold       (Uoldbox,     i,j,k,fluid);
    Cell U          (Ubox,        i,j,k,fluid);
    Cell UL         (ULbox,       i,j,k,fluid);
    Cell UR         (URbox,       i,j,k,fluid);
    Cell UStar      (UStarbox,    i,j,k,fluid);

    Real SR    = std::max(std::abs(UL(VELOCITY,probe_m[0],x))+UL(SOUNDSPEED,probe_m[0]),std::abs(UR(VELOCITY,probe_m[1],x))+UR(SOUNDSPEED,probe_m[1]));
    Real SL    = -SR;

    Real Sstar = getSstarFromDifferentMaterials(UL,UR,SL,SR,i,j,k,x,probe_m);

    getStarState(UL,UStar,SL,Sstar,parameters,x,probe_m[0]);
    getStarState(UR,UStar,SR,Sstar,parameters,x,probe_m[1]);

    U = UStar;

    Ubox.conservativeToPrimitive(i,j,k);

    //originalGFM(Uold,U,UL,UR,probe_m,validMat);

    //exactSolver_FindIntermediateState(UL,UR,parameters,U,probe_m);

}

void boxGhostFluidValues(BoxAccessCellArray& U, BoxAccessCellArray& U1, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& Ustar, BoxAccessCellArray& Ustarstar, BoxAccessLevelSet& LS0, const Real* dx, const Real* prob_lo, ParameterStruct& parameters)
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
                            Vector<Real> err;
                            err.push_back(current_x);
                            err.push_back(current_y);
                            err.push_back(nx);
                            err.push_back(ny);
                            err.push_back(interface_x);
                            err.push_back(interface_y);
                            err.push_back(probe_x[0]);
                            err.push_back(probe_y[0]);
                            err.push_back(probe_x[1]);
                            err.push_back(probe_y[1]);
                            err.push_back(probe_m[0]);
                            err.push_back(probe_m[1]);

                            std::string mess = "Error in probe extrapolation: materials are the same";

                            customAbort(err,mess);

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

                    if(parameters.SOLID)
                    {
                        getBothStarStatesAndSetNewValue_5Wave(i,j,k,U1,UL,UR,Ustar,Ustarstar,probe_m,parameters);
                    }
                    else
                    {
                        getBothStarStatesAndSetNewValue_3Wave(i,j,k,U,U1,UL,UR,Ustar,Ustarstar,probe_m,parameters,m);
                    }

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

void boxGhostFluidValues_1D(BoxAccessCellArray& U, BoxAccessCellArray& U1, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& Ustar, BoxAccessCellArray& Ustarstar, BoxAccessLevelSet& LS0, const Real* dx, const Real* prob_lo, ParameterStruct& parameters)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    Real current_x,   current_y,  current_z;



    Vector<int>  probe_m(2);

    for 		(int k = lo.z; k <= hi.z; ++k)
    {
                current_z = prob_lo[2] + (Real(k)+0.5)*dx[2];

        for 	(int j = lo.y; j <= hi.y; ++j)
        {
                current_y = prob_lo[1] + (Real(j)+0.5)*dx[1];

            for (int i = lo.x; i <= hi.x; ++i)
            {

                current_x = prob_lo[0] + (Real(i)+0.5)*dx[0];

                if(sgn<Real,int>(LS0(i,j,k,0)) != sgn<Real,int>(LS0(i+1,j,k,0)))
                {
                    probe_m[0] = LS0.whatMaterialIsValid(i,j,k);
                    probe_m[1] = LS0.whatMaterialIsValid(i+1,j,k);

                    U1(i  ,j,k,P,probe_m[1]) = U(i  ,j,k,P,probe_m[0]);
                    U1(i+1,j,k,P,probe_m[0]) = U(i+1,j,k,P,probe_m[1]);

                    U1(i  ,j,k,VELOCITY,probe_m[1],x) = U(i  ,j,k,VELOCITY,probe_m[0],x);
                    U1(i+1,j,k,VELOCITY,probe_m[0],x) = U(i+1,j,k,VELOCITY,probe_m[1],x);
                    U1(i  ,j,k,VELOCITY,probe_m[1],y) = U(i  ,j,k,VELOCITY,probe_m[0],y);
                    U1(i+1,j,k,VELOCITY,probe_m[0],y) = U(i+1,j,k,VELOCITY,probe_m[1],y);
                    U1(i  ,j,k,VELOCITY,probe_m[1],z) = U(i  ,j,k,VELOCITY,probe_m[0],z);
                    U1(i+1,j,k,VELOCITY,probe_m[0],z) = U(i+1,j,k,VELOCITY,probe_m[1],z);

                    U1(i  ,j,k,RHO,probe_m[1]) = U1(i+1,j,k,RHO,probe_m[1])*std::pow(U1(i  ,j,k,P,probe_m[1])/U1(i+1,j,k,P,probe_m[1]), 1.0/1.4);
                    U1(i+1,j,k,RHO,probe_m[0]) = U1(i  ,j,k,RHO,probe_m[0])*std::pow(U1(i+1,j,k,P,probe_m[0])/U1(i  ,j,k,P,probe_m[0]), 1.0/1.4);

                    U1.primitiveToConservative(i,j,k);
                }
                else if(sgn<Real,int>(LS0(i,j,k,0)) != sgn<Real,int>(LS0(i-1,j,k,0)))
                {
                    continue;
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

void setGhostFluidValues(MultiFab& S_new,CellArray& U, CellArray& U1, CellArray& UL, CellArray& UR, CellArray& Ustar, CellArray& Ustarstar, LevelSet& LS0, const Real* dx, const Real* prob_lo, ParameterStruct& parameters, Geometry& geom, Vector<BCRec>& bc)
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
        BoxAccessCellArray Ustarstarbox(mfi,bx,Ustarstar);
        BoxAccessLevelSet  bals (mfi,bx,LS0);

        boxGhostFluidValues(Ubox,U1box,ULbox,URbox,Ustarbox,Ustarstarbox,bals,dx,prob_lo,parameters);
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

    FillDomainBoundary(U1.data, geom, bc);

    U1.primitiveToConservative();

    return;
}


