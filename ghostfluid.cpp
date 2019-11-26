#include "simulationheader.h"

double getSstar(Cell& UL, Cell& UR, Real SL, Real SR, int i, int j, int k, Direction_enum d, int m);

/*void getStarState(BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, int m)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    Real SL,SR,Sstar;

    for 		(int k = lo.z; k <= hi.z; ++k)
    {
        for 	(int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                Cell UL(URbox,i,j,k,fluid);
                Cell UR(ULbox,i,j,k,fluid);
                Cell UStar(UStarbox,i,j,k,fluid);


                SR = std::max(std::abs(UL(VELOCITY,m,d))+UL(SOUNDSPEED,m),std::abs(UR(VELOCITY,m,d))+UR(SOUNDSPEED,m));
                SL = -SR;

                Sstar = getSstar(UL,UR,SL,SR,i,j,k,d,m);

                if(SL>=0.0)
                {
                    UStar = UL;
                }
                else if(Sstar[m]>=0.0)
                {
                    getStarState(UL,UStar,SL,Sstar,m);
                }
                else if(SR[m]>=0.0)
                {
                    getStarState(UR,UStar,SR,Sstar,m);
                }
                else
                {
                    UStar = UR;
                }
            }
        }
    }
}*/

void boxGhostFluidValues(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& ULstar, BoxAccessCellArray& URstar, BoxAccessLevelSet& LS0, const Real* dx, const Real* prob_lo)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    Real nx,ny;
    Real interface_x, interface_y;
    Real current_x,   current_y,  current_z;

    Vector<Real> probe_x(2);
    Vector<Real> probe_y(2);

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


                    UL.bilinearInterpolation(U,i,j,k,dx,prob_lo,probe_x[0],probe_y[0],m);
                    UR.bilinearInterpolation(U,i,j,k,dx,prob_lo,probe_x[1],probe_y[1],m);

                    //Print()<<"MombefL:\t " << UL(i,j,k,RHOU,m,0) << " " << UL(i,j,k,RHOU,m,1) << std::endl;
                    //Print()<<"MombefR:\t " << UR(i,j,k,RHOU,m,0) << " " << UR(i,j,k,RHOU,m,1) << std::endl;

                    UL.rotateFrameSoXPointsAlongNormal(i,j,k,nx,ny,m);
                    UR.rotateFrameSoXPointsAlongNormal(i,j,k,nx,ny,m);

                    //Print()<<"MomrotL:\t " << UL(i,j,k,RHOU,m,0) << " " << UL(i,j,k,RHOU,m,1) << std::endl;
                    //Print()<<"MomrotR:\t " << UR(i,j,k,RHOU,m,0) << " " << UR(i,j,k,RHOU,m,1) << std::endl;

                }
            }
        }
    }

    //UL.conservativeToPrimitive();
    //UR.conservativeToPrimitive();

    //UL.getSoundSpeed();
    //UR.getSoundSpeed();

   // getStarState(UL,UR,Ustar);

    //Ustar.conservativeToPrimitive();

    //propagateValues...

    return;
}

void setGhostFluidValues(MultiFab& S_new,CellArray& U, CellArray& UL, CellArray& UR, CellArray& ULstar, CellArray &URstar, LevelSet& LS0, const Real* dx, const Real* prob_lo)
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
        BoxAccessCellArray ULstarbox(mfi,bx,ULstar);
        BoxAccessCellArray URstarbox(mfi,bx,URstar);
        BoxAccessLevelSet  bals (mfi,bx,LS0);

        boxGhostFluidValues(Ubox,ULbox,URbox,ULstarbox,URstarbox,bals,dx,prob_lo);
    }

    return;
}


