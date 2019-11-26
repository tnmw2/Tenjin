#include "simulationheader.h"

void getStarState(Cell& U, Cell& UStar, Real SK, Real Sstar, ParameterStruct& parameters, Direction_enum d, int m);

double getSstarFromDifferentMaterials(Cell& UL, Cell& UR, Real SL, Real SR, int i, int j, int k, Direction_enum d, Vector<int>& m)
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

    UStarbox.conservativeToPrimitive(i,j,k);

    U = UStar;
}

void boxGhostFluidValues(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& Ustar, BoxAccessLevelSet& LS0, const Real* dx, const Real* prob_lo, ParameterStruct& parameters)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    Real nx,ny;
    Real interface_x, interface_y;
    Real current_x,   current_y,  current_z;

    Vector<Real> probe_x(2);
    Vector<Real> probe_y(2);
    Vector<int>  probe_m(2);

    Vector< std::pair<int,int> > coordList;


    int m;

    for 		(int k = lo.z; k <= hi.z; ++k)
    {
                current_z = prob_lo[2] + (Real(k)+0.5)*dx[2];

        for 	(int j = lo.y; j <= hi.y; ++j)
        {
                current_y = prob_lo[1] + (Real(j)+0.5)*dx[1];

            for (int i = lo.x; i <= hi.x; ++i)
            {
                coordList.resize(1);

                current_x = prob_lo[0] + (Real(i)+0.5)*dx[0];

                if(LS0.cellIsNextToAnInterface(i,j,k,0))
                {
                    coordList.push_back(std::make_pair(i,j));

                    //Print()<<"Point  : \t " << i << " " << j << " " << k << std::endl;
                    //Print()<<"Physic : \t " << current_x << " " <<current_y << std::endl;

                    LS0.calculateNormal(i,j,k,0,dx,nx,ny);

                    //Print()<<"Normal :\t " << nx << " " << ny << std::endl;

                    LS0.calculateInterpolationPoint(i,j,k,0,dx,nx,ny,current_x,current_y,interface_x,interface_y);

                    //Print()<<"Surface:\t " << interface_x << " " << interface_y << std::endl;

                    LS0.calculateProbes(i,j,k,0,dx,nx,ny,interface_x,interface_y,probe_x,probe_y);

                    //Print()<<"Probe1 :\t " << probe_x[0] << " " << probe_y[0] << std::endl;
                    //Print()<<"Probe2 :\t " << probe_x[1] << " " << probe_y[1] << std::endl;



                    if(LS0(i,j,k,0) > 0.0)
                    {
                        m = 0;
                    }
                    else
                    {
                        m = 1;
                    }

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

                    UL.bilinearInterpolation(U,i,j,k,dx,prob_lo,probe_x[0],probe_y[0],probe_m[0]);
                    UR.bilinearInterpolation(U,i,j,k,dx,prob_lo,probe_x[1],probe_y[1],probe_m[1]);

                    //Print()<<"MombefL:\t " << UL(i,j,k,RHOU,m,0) << " " << UL(i,j,k,RHOU,m,1) << std::endl;
                    //Print()<<"MombefR:\t " << UR(i,j,k,RHOU,m,0) << " " << UR(i,j,k,RHOU,m,1) << std::endl;

                    UL.rotateFrameSoXPointsAlongNormal(i,j,k,nx,ny,probe_m[0]);
                    UR.rotateFrameSoXPointsAlongNormal(i,j,k,nx,ny,probe_m[1]);

                    //Print()<<"MomrotL:\t " << UL(i,j,k,RHOU,m,0) << " " << UL(i,j,k,RHOU,m,1) << std::endl;
                    //Print()<<"MomrotR:\t " << UR(i,j,k,RHOU,m,0) << " " << UR(i,j,k,RHOU,m,1) << std::endl;

                    UL.primitiveToConservative(i,j,k);
                    UR.primitiveToConservative(i,j,k);

                    UL.getSoundSpeed(i,j,k);
                    UR.getSoundSpeed(i,j,k);

                    getBothStarStatesAndSetNewValue(i,j,k,U,UL,UR,Ustar,probe_m,parameters);

                    U.rotateFrameBack(i,j,k,nx,ny);

                    U.primitiveToConservative(i,j,k);



                }
            }
        }
    }

    return;
}

void setGhostFluidValues(MultiFab& S_new,CellArray& U, CellArray& UL, CellArray& UR, CellArray& Ustar, LevelSet& LS0, const Real* dx, const Real* prob_lo, ParameterStruct& parameters)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        BoxAccessCellArray Ubox (mfi,bx,U );
        BoxAccessCellArray ULbox(mfi,bx,UL);
        BoxAccessCellArray URbox(mfi,bx,UR);
        BoxAccessCellArray Ustarbox(mfi,bx,Ustar);
        BoxAccessLevelSet  bals (mfi,bx,LS0);

        boxGhostFluidValues(Ubox,ULbox,URbox,Ustarbox,bals,dx,prob_lo,parameters);
    }

    return;
}


