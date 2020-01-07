#include "levelset.h"
#include "tensor.h"

BoxAccessLevelSet::BoxAccessLevelSet(MFIter& mfi, const Box& bx, LevelSet &U) : box{bx}, fab{U.data[mfi]}, NLevelSets{U.NLevelSets}{}

void  BoxAccessLevelSet::initialise(const Real* dx, const Real* prob_lo)
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Real x,y,z;

    for             (int n = 0; n < NLevelSets; n++)
    {

        for 		(int k = lo.z; k <= hi.z; ++k)
        {
                    z = prob_lo[2] + (Real(k)+0.5)*dx[2];

            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                    y = prob_lo[1] + (Real(j)+0.5)*dx[1];

                for (int i = lo.x; i <= hi.x; ++i)
                {
                    x = prob_lo[0] + (Real(i)+0.5)*dx[0];

                    //(*this)(i,j,k,n) = (0.2 -sqrt(std::abs((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))));// ( (0.2 -sqrt(std::abs((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)))) < 0.0 ? -1E20 : 1E20);

                    //(*this)(i,j,k,n) = (1.0-(x+y))/sqrt(2.0);

                    //(*this)(i,j,k,n) = 0.5-x;

                    (*this)(i,j,k,n) = 0.005-x;

                }
            }
        }
    }
}

Real& BoxAccessLevelSet::operator()(int i, int j, int k, int n)
{
    return (fab.array())(i, j, k, n);
}

Real BoxAccessLevelSet::levelSetDerivative(Real v, const Real* dx, int dir, int i, int j, int k, int n)
{
    if(v > 0.0)
    {
        return D1(dir,i,j,k,n,-1,dx);
    }
    else
    {
        return D1(dir,i,j,k,n,1,dx);
    }
}

Real BoxAccessLevelSet::D1(int dir, int i, int j, int k, int n, int sign, const Real* dx)
{
    int extra[3] = {0,0,0};

    extra[dir]=1;

    if(sign>0)
    {
       return ((*this)(i+extra[0],j+extra[1],k+extra[2],n)-(*this)(i,j,k,n))/dx[dir];
    }
    else
    {
       return ((*this)(i,j,k,n)-(*this)(i-extra[0],j-extra[1],k-extra[2],n))/dx[dir];
    }
}

void BoxAccessLevelSet::resetLevelSet()
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for             (int n = 0; n < NLevelSets; n++)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    if(cellIsNextToAnInterface(i,j,k,n))
                    {
                        continue;
                    }
                    else
                    {
                        (*this)(i,j,k,n) = sgn<Real,Real>((*this)(i,j,k,n))*1E20;
                    }
                }
            }
        }
    }
}

bool BoxAccessLevelSet::customComparator(int i, int lim, int sense)
{
    if(sense > 0)
    {
        return i<=lim;
    }
    else
    {
        return i>=lim;
    }
}

void BoxAccessLevelSet::customChanger(int& i, int sense)
{
    if(sense > 0)
    {
        i++;
    }
    else
    {
        i--;
    }

    return;
}

void BoxAccessLevelSet::fastSweep(const Real* dx, int xsense, int ysense, int sign)
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Real phix,phiy;
    Real phinew;
    Real discriminant;

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


    for             (int n = 0; n < NLevelSets; n++)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = ylo; customComparator(j,yhi,ysense) ; customChanger(j,ysense))
            {
                for (int i = xlo; customComparator(i,xhi,xsense) ; customChanger(i,xsense))
                {
                    if(sign == 1)
                    {
                        if((sgn<Real,int>((*this)(i,j,k,n)) > 0) && !cellIsNextToAnInterface(i,j,k,n))
                        {

                            phix = std::min((*this)(i+1,j  ,k  ,n),(*this)(i-1,j  ,k  ,n));
                            phiy = std::min((*this)(i  ,j+1,k  ,n),(*this)(i  ,j-1,k  ,n));

                            discriminant = (phix/dx2+phiy/dy2)*(phix/dx2+phiy/dy2)-(1.0/dx2+1.0/dy2)*( (phix*phix/dx2+phiy*phiy/dy2)  -1.0);


                            if(discriminant < 0.0)
                            {
                                if(std::abs(phix) < std::abs(phiy))
                                {
                                    discriminant = (1.0/dx2);

                                    phinew = phix+sqrt(discriminant)*dx2;
                                }
                                else
                                {
                                    discriminant = (1.0/dy2);

                                    phinew = phiy+sqrt(discriminant)*dy2;
                                }
                            }
                            else
                            {
                                phinew = ((phix/dx2+phiy/dy2)+sqrt(discriminant))/(1.0/dx2+1.0/dy2);
                            }

                            if(phinew < (*this)(i,j,k,n))
                            {
                                (*this)(i,j,k,n) = phinew;
                            }
                        }
                    }
                    else if(sign == -1)
                    {
                        if((sgn<Real,int>((*this)(i,j,k,n)) < 0) && !cellIsNextToAnInterface(i,j,k,n))
                        {

                            phix = std::max((*this)(i+1,j  ,k  ,n),(*this)(i-1,j  ,k  ,n));
                            phiy = std::max((*this)(i  ,j+1,k  ,n),(*this)(i  ,j-1,k  ,n));

                            discriminant = (phix/dx2+phiy/dy2)*(phix/dx2+phiy/dy2)-(1.0/dx2+1.0/dy2)*( (phix*phix/dx2+phiy*phiy/dy2)  -1.0);


                            if(discriminant < 0.0)
                            {
                                if(std::abs(phix) < std::abs(phiy))
                                {
                                    discriminant = (1.0/dx2);

                                    phinew = phix-sqrt(discriminant)*dx2;
                                }
                                else
                                {
                                    discriminant = (1.0/dy2);

                                    phinew = phiy-sqrt(discriminant)*dy2;
                                }
                            }
                            else
                            {
                                phinew = ((phix/dx2+phiy/dy2)-sqrt(discriminant))/(1.0/dx2+1.0/dy2);
                            }

                            if(phinew > (*this)(i,j,k,n))
                            {
                                (*this)(i,j,k,n) = phinew;
                            }
                        }
                    }
                }
            }
        }
    }

    return;
}

bool BoxAccessLevelSet::cellIsNextToAnInterface(int i, int j, int k, int n, int limiter)
{
    Vector<int> pm {-1,1};


    if(limiter<0)
    {
        for(auto dif : pm)
        {
            if( (sgn<Real,int>((*this)(i,j,k,n)) != sgn<Real,int>((*this)(i+dif,j,k,n))) )
            {
                return true;
            }
            if( (sgn<Real,int>((*this)(i,j,k,n)) != sgn<Real,int>((*this)(i,j+dif,k,n))) )
            {
                return true;
            }
        }
    }
    else
    {
        if(limiter == x)
        {
            for(auto dif : pm)
            {
                if( (sgn<Real,int>((*this)(i,j,k,n)) != sgn<Real,int>((*this)(i+dif,j,k,n))) )
                {
                    return true;
                }
            }
        }
        if(limiter == y)
        {
            for(auto dif : pm)
            {
                if( (sgn<Real,int>((*this)(i,j,k,n)) != sgn<Real,int>((*this)(i,j+dif,k,n))) )
                {
                    return true;
                }
            }
        }
    }


    return false;
}

void BoxAccessLevelSet::calculateNormal(int i , int j , int k, int n, const Real* dx, Real& nx, Real& ny)
{
    nx = ((*this)(i+1,j,   k,n)-(*this)(i-1,j,  k,n))/(2.0*dx[0]);
    ny = ((*this)(i   ,j+1,k,n)-(*this)(i,  j-1,k,n))/(2.0*dx[1]);

    if(std::isnan(nx) || std::isnan(ny))
    {
        Print() << i << " " << j << std::endl;
        Print() << nx << " " << ny << std::endl;
        Abort("Nan in normal calculation");
    }

    Real norm = sqrt(nx*nx+ny*ny);

    nx = nx/norm;
    ny = ny/norm;

    if((std::abs(nx) > 1.0) || (std::abs(ny) > 1.0) )
    {
        Print() << i << " " << j << std::endl;
        Print() << nx << " " << ny << std::endl;
        Abort("Normal is not normalised");
    }

    return;
}

void BoxAccessLevelSet::calculateInterpolationPoint(int i , int j , int k, int n, const Real* dx, Real& nx, Real& ny, Real& cx, Real& cy, Real& ix, Real& iy)
{
    ix = cx - (*this)(i,j,k,n)*nx;
    iy = cy - (*this)(i,j,k,n)*ny;

    return;
}

void BoxAccessLevelSet::calculateProbes(int i , int j , int k, int n, const Real* dx, Real& nx, Real& ny, Real& ix, Real& iy, Vector<Real>& px, Vector<Real>& py)
{
    static const Real probe_length = 1.5;

    px[0] = ix - probe_length*dx[0]*nx;
    px[1] = ix + probe_length*dx[0]*nx;

    py[0] = iy - probe_length*dx[1]*ny;
    py[1] = iy + probe_length*dx[1]*ny;

    if( (std::abs(px[0]-ix) > 2.0*probe_length*dx[0]) || (std::abs(px[1]-ix) > 2.0*probe_length*dx[0]))
    {
        Abort("Error in probe calculation");
    }
    if( (std::abs(py[0]-iy) > 2.0*probe_length*dx[1]) || (std::abs(py[1]-iy) > 2.0*probe_length*dx[1]))
    {
        Abort("Error in probe calculation");
    }

}

bool BoxAccessLevelSet::cellIsValid(int i, int j, int k, int m)
{
    if(m == 0)
    {
        return ((*this)(i,j,k,0) > 0.0);
    }
    else
    {
        return ((*this)(i,j,k,0) < 0.0);
    }
}

int BoxAccessLevelSet::whatMaterialIsValid(int i, int j, int k)
{
    if((*this)(i,j,k,0) > 0.0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

bool BoxAccessLevelSet::cellIsNearInterface(int i, int j, int k, const Real* dx)
{
    return (std::abs((*this)(i,j,k,0)) < 10.0*dx[0]);
}

