#include "simulationheader.h"

void THINCArray::addTHINCvariable(Variable var, int materialNumber, int rowNumber, int colNumber)
{
    if(materialNumber*rowNumber*colNumber > 0)
    {
        for(int m = 0; m < materialNumber ;m++)
        {
            for(int row = 0; row < rowNumber ; row++)
            {
                for(int col = 0; col < colNumber ; col++)
                {
                    THINCvariables.push_back(MaterialSpecifier(var,m,row,col));
                }
            }
        }
    }
}

THINCArray::THINCArray(BoxArray& ba, DistributionMapping& dm, const int Nghost, ParameterStruct &parameters) : data(ba,dm,parameters.numberOfMaterials+1,Nghost)
{
    addTHINCvariable(ALPHA,parameters.numberOfMaterials);
    //addTHINCvariable(ALPHARHO,parameters.numberOfMaterials);
    //addTHINCvariable(P);
    //addTHINCvariable(VELOCITY,0,3);
    //addTHINCvariable(V_TENSOR,0,3,3);

    //udaykumar groove had just alpha i think?

}

BoxAccessTHINCArray::BoxAccessTHINCArray(MFIter& mfi, const Box &bx, THINCArray &U) : box{bx}, iab{U.data[mfi]}, THINCvariables{U.THINCvariables}{}

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
Real youngsInterfaceConstruction(Real& nx, Real& ny, BoxAccessCellArray& U, ParameterStruct& parameters, Direction_enum d,const Real* dx, int m, int i, int j, int k)
{
    Real norm;

    nx = -1.0/(8.0*dx[0])*(U.neighbour(1,1,0,i,j,k,ALPHA,m)+2.0*U.neighbour(1,0,0,i,j,k,ALPHA,m)+U.neighbour(1,-1,0,i,j,k,ALPHA,m)-U.neighbour(-1,1,0,i,j,k,ALPHA,m) -2.0*U.neighbour(-1,0,0,i,j,k,ALPHA,m)-U.neighbour(-1,-1,0,i,j,k,ALPHA,m));
    ny = -1.0/(8.0*dx[1])*(U.neighbour(1,1,0,i,j,k,ALPHA,m)+2.0*U.neighbour(0,1,0,i,j,k,ALPHA,m)+U.neighbour(-1,1,0,i,j,k,ALPHA,m)-U.neighbour(-1,-1,0,i,j,k,ALPHA,m)-2.0*U.neighbour(0,-1,0,i,j,k,ALPHA,m)-U.neighbour(1,-1,0,i,j,k,ALPHA,m));

    norm = std::sqrt(nx*nx+ny*ny);

    nx = nx/norm;
    ny = ny/norm;

    /*if(norm < 1E-20)
    {
        return 0;
    }
    else
    {*/
        if(d==x)
        {
            return std::abs(nx);
        }
        else
        {
            return std::abs(ny);
        }
    //}

}

/** Calculates the total boundary variation of any variable for a reconstructed cell. Compares against both MUSCL and THINC.
 */
double TBV(BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& UTHINC_L, BoxAccessCellArray& UTHINC_R, BoxAccessCellArray& UPL, BoxAccessCellArray& UPR, Direction_enum d, int i, int j, int k, MaterialSpecifier m)
{
    using namespace std;

    return std::min(std::min(abs(UR.left(d,i,j,k,m)-UPL(i,j,k,m))+abs(UPR(i,j,k,m)-UL.right(d,i,j,k,m)),abs(UTHINC_R.left(d,i,j,k,m)-UPL(i,j,k,m))+abs(UPR(i,j,k,m)-UTHINC_L.right(d,i,j,k,m))),std::min(abs(UR.left(d,i,j,k,m)-UPL(i,j,k,m))+abs(UPR(i,j,k,m)-UTHINC_L.right(d,i,j,k,m)),abs(UTHINC_R.left(d,i,j,k,m)-UPL(i,j,k,m))+abs(UPR(i,j,k,m)-UL.right(d,i,j,k,m))));

}

void BoxAccessTHINCArray::variableTHINCreconstruction(MaterialSpecifier& n, BoxAccessCellArray& U, BoxAccessCellArray& UTHINC_L, BoxAccessCellArray& UTHINC_R, Vector<Real>& beta, Vector<Real>& coshBeta, Vector<Real>& tanhBeta, Real epsilon, int i, int j ,int k, Direction_enum d)
{
    Real min         = std::min(U.left(d,i,j,k,n),U.right(d,i,j,k,n));
    Real max         = std::max(U.left(d,i,j,k,n),U.right(d,i,j,k,n))-min;

    Real theta       = sgn<Real,Real>(U.right(d,i,j,k,n)-U.left(d,i,j,k,n));

    Real C = (U(i,j,k,n)-min+epsilon)/(max+epsilon);
    Real B = std::exp(theta*beta[n.mat]*((2.0*C)-1.0));
    Real A = ( (B/coshBeta[n.mat])-1.0 )/tanhBeta[n.mat];

    UTHINC_R(i,j,k,n) = min + (max/2.0)*(1.0+ theta*((tanhBeta[n.mat]+A)/(1.0+A*tanhBeta[n.mat])) );
    UTHINC_L(i,j,k,n) = min + (max/2.0)*(1.0+ theta*A);
}

/** Performs the BVD-THINC update from Deng for the volume fraction.
 */
void BoxAccessTHINCArray::THINCreconstruction(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& UTHINC_L, BoxAccessCellArray& UTHINC_R, ParameterStruct& parameters,const Real* dx, Direction_enum d)
{
    Real beta0 = parameters.THINCbeta;
    Vector<Real> beta(U.numberOfMaterials);
    Vector<Real> coshBeta(U.numberOfMaterials);
    Vector<Real> tanhBeta(U.numberOfMaterials);

    Real epsilon = 1E-20;
    Vector<Real> normalisedVectorComponent(U.numberOfMaterials);

    Vector<Real> nx(U.numberOfMaterials);
    Vector<Real> ny(U.numberOfMaterials);


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


                    if(mixedCell(C,U,d,m,i,j,k) || mixedCellFlag(i,j,k))
                    {

                        mixedCellFlag(i,j,k) = 1;

                    }
                }

                if(mixedCellFlag(i,j,k))
                {
                    for(int m = 0; m < U.numberOfMaterials;m++)
                    {
                        if(AMREX_SPACEDIM == 2)
                        {
                            normalisedVectorComponent[m] = youngsInterfaceConstruction(nx[m],ny[m],U,parameters,d,dx,m,i,j,k);
                        }
                        else if(AMREX_SPACEDIM == 1)
                        {
                            normalisedVectorComponent[m] = 1.0;
                        }
                        else
                        {
                            Print() << "Haven't implemented THINC in 3D yet" << std::endl;
                        }

                        beta[m] = normalisedVectorComponent[m]*beta0 + 0.01;

                        coshBeta[m] = cosh(beta[m]);
                        tanhBeta[m] = tanh(beta[m]);


                    }

                    for(auto n : THINCvariables)
                    {
                        variableTHINCreconstruction(n,U,UTHINC_L,UTHINC_R,beta,coshBeta,tanhBeta,epsilon,i,j,k,d);
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
                if(mixedCellFlag(i,j,k))
                {
                    for(int m = 0; m < U.numberOfMaterials;m++)
                    {
                        TBVFlag(i,j,k,m) = 1;
                    }

                    for(auto n : THINCvariables)
                    {
                        TBVMUSCL = TBV(UL,UR,UTHINC_L,UTHINC_R,UL,      UR,      d,i,j,k,n);
                        TBVTHINC = TBV(UL,UR,UTHINC_L,UTHINC_R,UTHINC_L,UTHINC_R,d,i,j,k,n);

                        if(TBVTHINC < TBVMUSCL)
                        {
                            UL(i,j,k,n) = UTHINC_L(i,j,k,n);
                            UR(i,j,k,n) = UTHINC_R(i,j,k,n);
                        }
                        else
                        {
                           TBVFlag(i,j,k,0) = 0;
                        }
                    }

                    /*if(TBVFlag(i,j,k,0))
                    {
                        for(auto n : THINCvariables)
                        {
                            UL(i,j,k,n) = UTHINC_L(i,j,k,n);
                            UR(i,j,k,n) = UTHINC_R(i,j,k,n);
                        }
                    }*/

                    /*for(int m = 0; m < U.numberOfMaterials; m++)
                    {
                        UL(i,j,k,RHO_K,m) = UL(i,j,k,ALPHARHO,m)/UL(i,j,k,ALPHA,m);
                        UR(i,j,k,RHO_K,m) = UR(i,j,k,ALPHARHO,m)/UR(i,j,k,ALPHA,m);
                    }*/
                }
            }
        }
    }



    return;

}
