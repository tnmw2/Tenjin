#include "levelset.h"
#include "tensor.h"

BoxAccessLevelSet::BoxAccessLevelSet(MFIter& mfi, const Box& bx, LevelSet &U) : box{bx}, fab{U.data[mfi]}, NLevelSets{U.NLevelSets}{}

void  BoxAccessLevelSet::initialise(InitialStruct& initial, const Real* dx, const Real* prob_lo)
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Real x,y,z;

    /*Real chamfer = 0.0;//1E-3;
    Real length  = 2.347E-2;
    Real radius  = initial.interface;*/

    Real bubbleCentre_x  = 0.7;
    Real radius          = 0.2;
    Real shock           = initial.interface;



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

                    //(*this)(i,j,k,n) = x-0.005;

                    //(*this)(i,j,k,n) = 0.005-x;

                    /*******************************
                     * Rod Impact
                     ******************************/

                    /*if(x<radius)
                    {
                        if(y< length-chamfer)
                        {
                            (*this)(i,j,k,n) = 0.5*dx[0];
                        }
                        else if(y<length)
                        {
                            if( (x< radius-chamfer) || (x - (radius-chamfer))*(x - (radius-chamfer))+(y - (length-chamfer))*(y - (length-chamfer)) < chamfer*chamfer     )
                            {
                                (*this)(i,j,k,n) = 0.5*dx[0];
                            }
                            else
                            {
                                (*this)(i,j,k,n) = -0.5*dx[0];
                            }
                        }
                        else
                        {
                            (*this)(i,j,k,n) = -0.5*dx[0];
                        }
                    }
                    else
                    {
                        (*this)(i,j,k,n) = -0.5*dx[0];
                    }*/


                    /*******************************
                     * Bubble
                     ******************************/


                    if(x<shock)
                    {
                        (*this)(i,j,k,n) = -0.5*dx[0];
                    }
                    else
                    {
                        if((x-bubbleCentre_x)*(x-bubbleCentre_x)+y*y < radius*radius) //  if((x > bubbleCentre_x-radius) && (x<bubbleCentre_x+radius) && (y < radius)) //
                        {
                            (*this)(i,j,k,n) = 0.5*dx[0];
                        }
                        else
                        {
                            (*this)(i,j,k,n) = -0.5*dx[0];
                        }
                    }
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
    /*********************************
     * First order spatial derivative
     ********************************/
    /*
    if(v > 0.0)
    {
        return D1(dir,i,j,k,n,-1,dx);
    }
    else
    {
        return D1(dir,i,j,k,n,1,dx);
    }*/

    /**********************************
     * Third order ENO spatial derivative
     ***********************************/
    /*int extra[3]       = {0,0,0};
    int extra_k[3]     = {0,0,0};
    int extra_kstar[3] = {0,0,0};

    Real c     = 0.0;
    Real cstar = 0.0;

    extra[dir] = 1;

    if(v > 0.0)
    {
        extra_k[dir]=-1;
    }

    if(std::abs(D2(dir,i+extra_k[0],j+extra_k[1],k+extra_k[2],n,0,dx)) <= std::abs(D2(dir,i+extra_k[0]+extra[0],j+extra_k[1]+extra[1],k+extra_k[2]+extra[2],n,0,dx)))
    {
        c = D2(dir,i+extra_k[0],j+extra_k[1],k+extra_k[2],n,0,dx);
        extra_kstar[dir] = extra_k[dir]-1;
    }
    else
    {
        c = D2(dir,i+extra_k[0]+extra[0],j+extra_k[1]+extra[1],k+extra_k[2]+extra[2],n,0,dx);
        extra_kstar[dir] = extra_k[dir];
    }

    if(std::abs(D3(dir,i+extra_kstar[0],j+extra_kstar[1],k+extra_kstar[2],n,-1,dx)) <= std::abs(D3(dir,i+extra_kstar[0],j+extra_kstar[1],k+extra_kstar[2],n,1,dx)))
    {
        cstar = D3(dir,i+extra_kstar[0],j+extra_kstar[1],k+extra_kstar[2],n,-1,dx);
    }
    else
    {
        cstar = D3(dir,i+extra_kstar[0],j+extra_kstar[1],k+extra_kstar[2],n, 1,dx);
    }

    return  D1(dir,i+extra_k[0],j+extra_k[1],k+extra_k[2],n,1,dx);// + c*(2.0*((Real)(-extra_k[dir]))-1.0)*dx[dir] + cstar*(3.0*((Real)(-extra_kstar[dir]))*((Real)(-extra_kstar[dir]))-6.0*((Real)(-extra_kstar[dir]))+2.0)*dx[dir]*dx[dir];
    */

    /**********************************
     * Fifth order WENO spatial derivative
     ***********************************/
    Real v1,v2,v3,v4,v5;
    Real S[3];
    Real a[3];
    Real w[3];

    Real epsilon = 1E-6;

    int extra[3]       = {0,0,0};

    extra[dir] = 1;

    if(v > 0.0)
    {
        v1 = D1(dir,i-extra[0]*2,j-extra[1]*2,k-extra[2]*2,n,-1,dx);
        v2 = D1(dir,i-extra[0]*1,j-extra[1]*1,k-extra[2]*1,n,-1,dx);
        v3 = D1(dir,i-extra[0]*0,j-extra[1]*0,k-extra[2]*0,n,-1,dx);
        v4 = D1(dir,i+extra[0]*1,j+extra[1]*1,k+extra[2]*1,n,-1,dx);
        v5 = D1(dir,i+extra[0]*2,j+extra[1]*2,k+extra[2]*2,n,-1,dx);
    }
    else
    {
        v1 = D1(dir,i+extra[0]* 3,j+extra[1]* 3,k+extra[2]* 3,n,-1,dx);
        v2 = D1(dir,i+extra[0]* 2,j+extra[1]* 2,k+extra[2]* 2,n,-1,dx);
        v3 = D1(dir,i+extra[0]* 1,j+extra[1]* 1,k+extra[2]* 1,n,-1,dx);
        v4 = D1(dir,i+extra[0]* 0,j+extra[1]* 0,k+extra[2]* 0,n,-1,dx);
        v5 = D1(dir,i+extra[0]*-1,j+extra[1]*-1,k+extra[2]*-1,n,-1,dx);
    }

    S[0] = (13.0/12.0)*(v1-2.0*v2+v3)*(v1-2.0*v2+v3)+0.25*(v1-4.0*v2+3.0*v3)*(v1-4.0*v2+3.0*v3);
    S[1] = (13.0/12.0)*(v2-2.0*v3+v4)*(v2-2.0*v3+v4)+0.25*(v2-v4)*(v2-v4);
    S[2] = (13.0/12.0)*(v3-2.0*v4+v5)*(v1-2.0*v2+v3)+0.25*(3.0*v3-4.0*v4+v5)*(3.0*v3-4.0*v4+v5);

    a[0] = (1.0/10.0)*(1.0/(epsilon+S[0]))*(1.0/(epsilon+S[0]));
    a[1] = (6.0/10.0)*(1.0/(epsilon+S[1]))*(1.0/(epsilon+S[1]));
    a[2] = (3.0/10.0)*(1.0/(epsilon+S[2]))*(1.0/(epsilon+S[2]));

    for(int c = 0; c<3; c++)
    {
        w[c] = a[c]/(a[0]+a[1]+a[2]);
    }

    return w[0]*((1.0/3.0)*v1-(7.0/6.0)*v2+(11.0/6.0)*v3)+w[1]*(-(1.0/6.0)*v2+(5.0/6.0)*v3+(1.0/3.0)*v4)+w[2]*((1.0/3.0)*v3+(5.0/6.0)*v4-(1.0/6.0)*v5);

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

Real BoxAccessLevelSet::D2(int dir, int i, int j, int k, int n, int sign, const Real* dx)
{
    return (D1(dir,i,j,k,n,1,dx) - D1(dir,i,j,k,n,-1,dx))/(2.0*dx[dir]);
}

Real BoxAccessLevelSet::D3(int dir, int i, int j, int k, int n, int sign, const Real* dx)
{
    int extra[3] = {0,0,0};

    extra[dir]=1;

    if(sign>0)
    {
       return (D2(dir,i+extra[0],j+extra[1],k+extra[2],n,0,dx)-D2(dir,i,j,k,n,0,dx))/(3.0*dx[dir]);
    }
    else
    {
       return (D2(dir,i,j,k,n,0,dx) - D2(dir,i-extra[0],j-extra[1],k-extra[2],n,0,dx))/(3.0*dx[dir]);
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
    Vector<int> pm {-1,0,1};


    if(limiter<0)
    {
        for(auto xdif : pm)
        {
            for(auto ydif : pm)
            {
                if(xdif == 0 && ydif == 0)
                {
                    continue;
                }
                else if( (sgn<Real,int>((*this)(i,j,k,n)) != sgn<Real,int>((*this)(i+xdif,j+ydif,k,n))) )
                {
                    return true;
                }
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

