#include "simulationheader.h"

THINCArray::THINCArray(BoxArray& ba, DistributionMapping& dm, const int Nghost, ParameterStruct &parameters) : data(ba,dm,parameters.numberOfMaterials+1,Nghost){}

BoxAccessTHINCArray::BoxAccessTHINCArray(MFIter& mfi, const Box &bx, THINCArray &U) : box{bx}, iab{U.data[mfi]}{}

int& BoxAccessTHINCArray::mixedCellFlag(int i, int j, int k)
{
    return (iab.array())(i, j, k, 0);
}

int& BoxAccessTHINCArray::TBVFlag(int i, int j, int k, int m)
{
    return (iab.array())(i, j, k, 1+m);
}

/** Determines whether THINC should be performed based on local values of the volume fraction.
 */
bool mixedCell(Real C, BoxAccessCellArray& U, Direction_enum d, int m, int i, int j, int k)
{
    double delta = 1E-4;

    return ((C > delta) && (C < (0.999999-delta)) && ((U.right(d,i,j,k,ALPHA,m) - U(i,j,k,ALPHA,m))*(U(i,j,k,ALPHA,m) - U.left(d,i,j,k,ALPHA,m)) >= 0.0));
}

/** Estimates the local interface normal using Young's method in 2D.
 */
Real youngsInterfaceConstruction(Real nx, Real ny, BoxAccessCellArray& U, ParameterStruct& parameters, Direction_enum d,const Real* dx, int m, int i, int j, int k)
{
    Real norm;

    nx = -1.0/(8.0*dx[0])*(U.neighbour(1,1,0,i,j,k,ALPHA,m)+2.0*U.neighbour(1,0,0,i,j,k,ALPHA,m)+U.neighbour(1,-1,0,i,j,k,ALPHA,m)-U.neighbour(-1,1,0,i,j,k,ALPHA,m) -2.0*U.neighbour(-1,0,0,i,j,k,ALPHA,m)-U.neighbour(-1,-1,0,i,j,k,ALPHA,m));
    ny = -1.0/(8.0*dx[1])*(U.neighbour(1,1,0,i,j,k,ALPHA,m)+2.0*U.neighbour(0,1,0,i,j,k,ALPHA,m)+U.neighbour(-1,1,0,i,j,k,ALPHA,m)-U.neighbour(-1,-1,0,i,j,k,ALPHA,m)-2.0*U.neighbour(0,-1,0,i,j,k,ALPHA,m)-U.neighbour(1,-1,0,i,j,k,ALPHA,m));

    norm = std::sqrt(nx*nx+ny*ny);

    if(norm < 1E-20)
    {
        return 0;
    }
    else
    {
        if(d==x)
        {
            return std::abs(nx/norm);
        }
        else
        {
            return std::abs(ny/norm);
        }
    }

}

/** Calculates the total boundary variation of any variable for a reconstructed cell. Compares against both MUSCL and THINC.
 */
double TBV(BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& UTHINC_L, BoxAccessCellArray& UTHINC_R, BoxAccessCellArray& UPL, BoxAccessCellArray& UPR, Direction_enum d, int i, int j, int k, MaterialSpecifier m)
{
    using namespace std;

    return std::min(std::min(abs(UR.left(d,i,j,k,m)-UPL(i,j,k,m))+abs(UPR(i,j,k,m)-UL.right(d,i,j,k,m)),abs(UTHINC_R.left(d,i,j,k,m)-UPL(i,j,k,m))+abs(UPR(i,j,k,m)-UTHINC_L.right(d,i,j,k,m))),std::min(abs(UR.left(d,i,j,k,m)-UPL(i,j,k,m))+abs(UPR(i,j,k,m)-UTHINC_L.right(d,i,j,k,m)),abs(UTHINC_R.left(d,i,j,k,m)-UPL(i,j,k,m))+abs(UPR(i,j,k,m)-UL.right(d,i,j,k,m))));

}

/** Performs the BVD-THINC update from Deng for the volume fraction.
 */
void BoxAccessTHINCArray::THINCreconstruction(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& UTHINC_L, BoxAccessCellArray& UTHINC_R, ParameterStruct& parameters,const Real* dx, Direction_enum d)
{

    Real beta;
    Real beta0 = parameters.THINCbeta;
    Real coshBeta;
    Real tanhBeta;
    Real epsilon = 1E-20;
    Real normalisedVectorComponent;

    Real nx;
    Real ny;

    Real min, max, theta;
    Real A, B, C ;

    Real TBVMUSCL, TBVTHINC;
    Real TBVMUSCLRHO, TBVTHINCRHO;

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    for    		   (int k = lo.z; k <= hi.z; ++k)
    {
        for        (int j = lo.y; j <= hi.y; ++j)
        {
            for    (int i = lo.x; i <= hi.x; ++i)
            {
                mixedCellFlag(i,j,k) = 0;

                for(int m = 0;    m < U.numberOfMaterials;m++)
                {
                    min         = std::min(U.left(d,i,j,k,ALPHA,m),U.right(d,i,j,k,ALPHA,m));
                    max         = std::max(U.left(d,i,j,k,ALPHA,m),U.right(d,i,j,k,ALPHA,m))-min;
                    theta       = sgn<Real,Real>(U.right(d,i,j,k,ALPHA,m)-U.left(d,i,j,k,ALPHA,m));

                    if(max < epsilon)
                    {
                        C = 0.0;
                    }
                    else
                    {
                        C = (U(i,j,k,ALPHA,m)-min+epsilon)/(max+epsilon);
                    }


                    if(mixedCell(C,U,d,m,i,j,k))
                    {
                        mixedCellFlag(i,j,k) = 1;

                        if(AMREX_SPACEDIM == 2)
                        {
                            normalisedVectorComponent = youngsInterfaceConstruction(nx,ny,U,parameters,d,dx,m,i,j,k);
                        }
                        else if(AMREX_SPACEDIM == 1)
                        {
                            normalisedVectorComponent = 1.0;
                        }
                        else
                        {
                            Print() << "Haven't implemented THINC in 3D yet" << std::endl;
                        }

                        beta = normalisedVectorComponent*beta0 + 0.01;

                        coshBeta = cosh(beta);
                        tanhBeta = tanh(beta);

                        B = std::exp(theta*beta*((2.0*C)-1.0));
                        A = ( (B/coshBeta)-1.0 )/tanhBeta;


                        UTHINC_R(i,j,k,ALPHA,m) = min + (max/2.0)*(1.0+ theta*((tanhBeta+A)/(1.0+A*tanhBeta)) );
                        UTHINC_L(i,j,k,ALPHA,m) = min + (max/2.0)*(1.0+ theta*A);

                        min         = std::min(U.left(d,i,j,k,ALPHARHO,m),U.right(d,i,j,k,ALPHARHO,m));
                        max         = std::max(U.left(d,i,j,k,ALPHARHO,m),U.right(d,i,j,k,ALPHARHO,m))-min;

                        theta       = sgn<Real,Real>(U.right(d,i,j,k,ALPHARHO,m)-U.left(d,i,j,k,ALPHARHO,m));

                         C = (U(i,j,k,ALPHARHO,m)-min+epsilon)/(max+epsilon);
                         B = std::exp(theta*beta*((2.0*C)-1.0));
                         A = ( (B/coshBeta)-1.0 )/tanhBeta;

                        UTHINC_R(i,j,k,ALPHARHO,m) = min + (max/2.0)*(1.0+ theta*((tanhBeta+A)/(1.0+A*tanhBeta)) );
                        UTHINC_L(i,j,k,ALPHARHO,m) = min + (max/2.0)*(1.0+ theta*A);


                    }
                }
            }
        }
    }

    for    		   (int k = lo.z; k <= hi.z; ++k)
    {
        for        (int j = lo.y; j <= hi.y; ++j)
        {
            for    (int i = lo.x; i <= hi.x; ++i)
            {
                for(int m = 0;    m < U.numberOfMaterials;m++)
                {
                    TBVMUSCL = TBV(UL,UR,UTHINC_L,UTHINC_R,UL,      UR,      d,i,j,k,MaterialSpecifier(ALPHA,m));
                    TBVTHINC = TBV(UL,UR,UTHINC_L,UTHINC_R,UTHINC_L,UTHINC_R,d,i,j,k,MaterialSpecifier(ALPHA,m));

                    TBVMUSCLRHO = TBV(UL,UR,UTHINC_L,UTHINC_R,UL,      UR,      d,i,j,k,MaterialSpecifier(ALPHARHO,m));
                    TBVTHINCRHO = TBV(UL,UR,UTHINC_L,UTHINC_R,UTHINC_L,UTHINC_R,d,i,j,k,MaterialSpecifier(ALPHARHO,m));


                    if( mixedCellFlag(i,j,k) && TBVTHINC < TBVMUSCL && TBVTHINCRHO < TBVMUSCLRHO)
                    {
                       TBVFlag(i,j,k,m) = 1;

                    }
                    else
                    {
                       TBVFlag(i,j,k,m) = 0;
                    }
                }
            }
        }
    }

    for    		   (int k = lo.z; k <= hi.z; ++k)
    {
        for        (int j = lo.y; j <= hi.y; ++j)
        {
            for    (int i = lo.x; i <= hi.x; ++i)
            {
                for(int m = 0;    m < U.numberOfMaterials;m++)
                {
                    if(TBVFlag(i,j,k,m))
                    {
                        UL(i,j,k,ALPHA,m)	= UTHINC_L(i,j,k,ALPHA,m);
                        UR(i,j,k,ALPHA,m)	= UTHINC_R(i,j,k,ALPHA,m);

                        /*UL(i,j,k,ALPHARHO,m)	= UTHINC_L(i,j,k,ALPHARHO,m);
                        UR(i,j,k,ALPHARHO,m)	= UTHINC_R(i,j,k,ALPHARHO,m);

                        UL(i,j,k,RHO_K,m)	= UL(i,j,k,ALPHARHO,m)/UL(i,j,k,ALPHA,m);
                        UR(i,j,k,RHO_K,m)	= UR(i,j,k,ALPHARHO,m)/UR(i,j,k,ALPHA,m);*/
                    }
                }
            }
        }
    }


    return;

}



/** Performs the BVD-THINC update from Deng for the volume fraction.
 */
void BoxAccessTHINCArray::cautiousTHINCreconstruction(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& UTHINC_L, BoxAccessCellArray& UTHINC_R, ParameterStruct& parameters,const Real* dx, Direction_enum d)
{

    Real beta;
    Real beta0 = parameters.THINCbeta;
    Real coshBeta;
    Real tanhBeta;
    Real epsilon = 1E-20;
    Real normalisedVectorComponent;

    Real nx;
    Real ny;

    Real min, max, theta;
    Real A, B, C ;

    Real TBVMUSCL, TBVTHINC;
    Real TBVMUSCLRHO, TBVTHINCRHO;

    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    for    		   (int k = lo.z; k <= hi.z; ++k)
    {
        for        (int j = lo.y; j <= hi.y; ++j)
        {
            for    (int i = lo.x; i <= hi.x; ++i)
            {
                mixedCellFlag(i,j,k) = 0;

                for(int m = 0;    m < U.numberOfMaterials;m++)
                {
                    min         = std::min(U.left(d,i,j,k,ALPHA,m),U.right(d,i,j,k,ALPHA,m));
                    max         = std::max(U.left(d,i,j,k,ALPHA,m),U.right(d,i,j,k,ALPHA,m))-min;
                    theta       = sgn<Real,Real>(U.right(d,i,j,k,ALPHA,m)-U.left(d,i,j,k,ALPHA,m));

                    if(max < epsilon)
                    {
                        C = 0.0;
                    }
                    else
                    {
                        C = (U(i,j,k,ALPHA,m)-min+epsilon)/(max+epsilon);
                    }

                    if(mixedCell(C,U,d,m,i,j,k))
                    {
                        mixedCellFlag(i,j,k) = 1;
                    }
                }
            }
        }
    }

    IntVect extra(AMREX_D_DECL(0,0,0));

    extra[d]=1;

    for    		   (int k = lo.z; k <= hi.z; ++k)
    {
        for        (int j = lo.y; j <= hi.y; ++j)
        {
            for    (int i = lo.x; i <= hi.x; ++i)
            {
                for(int m = 0;    m < U.numberOfMaterials;m++)
                {
                    min         = std::min(U.left(d,i,j,k,ALPHA,m),U.right(d,i,j,k,ALPHA,m));
                    max         = std::max(U.left(d,i,j,k,ALPHA,m),U.right(d,i,j,k,ALPHA,m))-min;
                    theta       = sgn<Real,Real>(U.right(d,i,j,k,ALPHA,m)-U.left(d,i,j,k,ALPHA,m));

                    if(max < epsilon)
                    {
                        C = 0.0;
                    }
                    else
                    {
                        C = (U(i,j,k,ALPHA,m)-min+epsilon)/(max+epsilon);
                    }

                    if(mixedCellFlag(i,j,k) && mixedCellFlag(i-extra[d],j-extra[d],k-extra[d]) && mixedCellFlag(i+extra[d],j+extra[d],k+extra[d]))
                    {
                        mixedCellFlag(i,j,k) = 2;

                        if(AMREX_SPACEDIM == 2)
                        {
                            normalisedVectorComponent = youngsInterfaceConstruction(nx,ny,U,parameters,d,dx,m,i,j,k);
                        }
                        else if(AMREX_SPACEDIM == 1)
                        {
                            normalisedVectorComponent = 1.0;
                        }
                        else
                        {
                            Print() << "Haven't implemented THINC in 3D yet" << std::endl;
                        }

                        beta = normalisedVectorComponent*beta0 + 0.01;

                        coshBeta = cosh(beta);
                        tanhBeta = tanh(beta);

                        B = std::exp(theta*beta*((2.0*C)-1.0));
                        A = ( (B/coshBeta)-1.0 )/tanhBeta;


                        UTHINC_R(i,j,k,ALPHA,m) = min + (max/2.0)*(1.0+ theta*((tanhBeta+A)/(1.0+A*tanhBeta)) );
                        UTHINC_L(i,j,k,ALPHA,m) = min + (max/2.0)*(1.0+ theta*A);

                       /* min         = std::min(U.left(d,i,j,k,ALPHARHO,m),U.right(d,i,j,k,ALPHARHO,m));
                        max         = std::max(U.left(d,i,j,k,ALPHARHO,m),U.right(d,i,j,k,ALPHARHO,m))-min;

                        theta       = sgn<Real,Real>(U.right(d,i,j,k,ALPHARHO,m)-U.left(d,i,j,k,ALPHARHO,m));

                         C = (U(i,j,k,ALPHARHO,m)-min+epsilon)/(max+epsilon);
                         B = std::exp(theta*beta*((2.0*C)-1.0));
                         A = ( (B/coshBeta)-1.0 )/tanhBeta;

                        UTHINC_R(i,j,k,ALPHARHO,m) = min + (max/2.0)*(1.0+ theta*((tanhBeta+A)/(1.0+A*tanhBeta)) );
                        UTHINC_L(i,j,k,ALPHARHO,m) = min + (max/2.0)*(1.0+ theta*A);*/


                    }
                }
            }
        }
    }

    for    		   (int k = lo.z; k <= hi.z; ++k)
    {
        for        (int j = lo.y; j <= hi.y; ++j)
        {
            for    (int i = lo.x; i <= hi.x; ++i)
            {
                for(int m = 0;    m < U.numberOfMaterials;m++)
                {
                    TBVMUSCL = TBV(UL,UR,UTHINC_L,UTHINC_R,UL,      UR,      d,i,j,k,MaterialSpecifier(ALPHA,m));
                    TBVTHINC = TBV(UL,UR,UTHINC_L,UTHINC_R,UTHINC_L,UTHINC_R,d,i,j,k,MaterialSpecifier(ALPHA,m));

                    //TBVMUSCLRHO = TBV(UL,UR,UTHINC_L,UTHINC_R,UL,      UR,      d,i,j,k,MaterialSpecifier(ALPHARHO,m));
                    //TBVTHINCRHO = TBV(UL,UR,UTHINC_L,UTHINC_R,UTHINC_L,UTHINC_R,d,i,j,k,MaterialSpecifier(ALPHARHO,m));


                    if( mixedCellFlag(i,j,k) == 2 && TBVTHINC < TBVMUSCL)// && TBVTHINCRHO < TBVMUSCLRHO)
                    {
                       TBVFlag(i,j,k,m) = 1;
                    }
                    else
                    {
                       TBVFlag(i,j,k,m) = 0;
                    }
                }
            }
        }
    }

    for    		   (int k = lo.z; k <= hi.z; ++k)
    {
        for        (int j = lo.y; j <= hi.y; ++j)
        {
            for    (int i = lo.x; i <= hi.x; ++i)
            {
                for(int m = 0;    m < U.numberOfMaterials;m++)
                {
                    if(TBVFlag(i,j,k,m))
                    {
                        UL(i,j,k,ALPHA,m)	= UTHINC_L(i,j,k,ALPHA,m);
                        UR(i,j,k,ALPHA,m)	= UTHINC_R(i,j,k,ALPHA,m);

                        /*UL(i,j,k,ALPHARHO,m)	= UTHINC_L(i,j,k,ALPHARHO,m);
                        UR(i,j,k,ALPHARHO,m)	= UTHINC_R(i,j,k,ALPHARHO,m);

                        UL(i,j,k,RHO_K,m)	= UL(i,j,k,ALPHARHO,m)/UL(i,j,k,ALPHA,m);
                        UR(i,j,k,RHO_K,m)	= UR(i,j,k,ALPHARHO,m)/UR(i,j,k,ALPHA,m);*/
                    }
                }
            }
        }
    }


    return;

}
