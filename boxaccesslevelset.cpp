#include "levelset.h"

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



                    (*this)(i,j,k,n) = 0.5- y;
                }
            }
        }
    }
}

Real& BoxAccessLevelSet::operator()(int i, int j, int k, int n)
{
    return (fab.array())(i, j, k, n);
}

void  BoxAccessLevelSet::advanceLevelSet(BoxAccessCellArray& U, BoxAccessLevelSet& LS, Real dt, const Real* dx)
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
                    (*this)(i,j,k,n) = LS(i,j,k,n);

                    for(int row = 0; row < AMREX_SPACEDIM; row++)
                    {
                        (*this)(i,j,k,n) -= dt*(U(i,j,k,VELOCITY,0,row)*levelSetDerivative(LS,U(i,j,k,VELOCITY,0,row),dx,row,i,j,k,n));
                    }
                }
            }
        }
    }
}

Real BoxAccessLevelSet::levelSetDerivative(BoxAccessLevelSet& LS, Real v, const Real* dx, int dir, int i, int j, int k, int n)
{
    if(v > 0.0)
    {
        return D1(LS,dir,i,j,k,n,-1,dx);
    }
    else
    {
        return D1(LS,dir,i,j,k,n,1,dx);
    }
}

Real BoxAccessLevelSet::D1(BoxAccessLevelSet& LS, int dir, int i, int j, int k, int n, int sign, const Real* dx)
{
    IntVect extra(AMREX_D_DECL(0,0,0));

    for(int n = 0; n < 3 ;n++)
    {
        extra[n] = 0;
    }

    extra[dir]=1;

    if(sign>0)
    {
       return (LS(i+extra[0],j+extra[1],k+extra[2],n)-LS(i,j,k,n))/dx[dir];
    }
    else
    {
       return (LS(i,j,k,n)-LS(i-extra[0],j-extra[1],k-extra[2],n))/dx[dir];
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
                    if(cellIsNextToAnInterface(i,j,k,n))//((sgn<Real,int>((*this)(i,j,k,n)) != sgn<Real,int>((*this)(i,j-1,k,n))) ||  (sgn<Real,int>((*this)(i,j,k,n)) != sgn<Real,int>((*this)(i,j+1,k,n)))) //
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

void BoxAccessLevelSet::fastSweep(const Real* dx, int dir, int sense, int sign)
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);


    IntVect extra(AMREX_D_DECL(0,0,0));

    for(int n = 0; n < 3 ;n++)
    {
        extra[n] = 0;
    }

    extra[dir]=1;

    Real phix,phiy;
    Real phinew;
    Real discriminant;

    Real dx2 = dx[0]*dx[0];
    Real dy2 = dx[1]*dx[1];

    if(sense > 0)
    {
        for             (int n = 0; n < NLevelSets; n++)
        {
            for 		(int k = lo.z; k <= hi.z; ++k)
            {
                for 	(int j = lo.y; j <= hi.y; ++j)
                {
                    for (int i = lo.x; i <= hi.x; ++i)
                    {
                        if(sign == 1)
                        {
                            if((sgn<Real,int>((*this)(i,j,k,n)) > 0) && (sgn<Real,int>((*this)(i+extra[0],j+extra[1],k+extra[2],n)) > 0) && (sgn<Real,int>((*this)(i-extra[0],j-extra[1],k-extra[2],n)) > 0))
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
                            if((sgn<Real,int>((*this)(i,j,k,n)) < 0) && (sgn<Real,int>((*this)(i+extra[0],j+extra[1],k+extra[2],n)) < 0) && (sgn<Real,int>((*this)(i-extra[0],j-extra[1],k-extra[2],n)) < 0))
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
    }
    else
    {
        for             (int n = 0; n < NLevelSets; n++)
        {
            for 		(int k = hi.z; k >= lo.z; --k)
            {
                for 	(int j = hi.y; j >= lo.y; --j)
                {
                    for (int i = hi.x; i >= lo.x; --i)
                    {
                        if(sign == 1)
                        {
                            if((sgn<Real,int>((*this)(i,j,k,n)) > 0) && (sgn<Real,int>((*this)(i+extra[0],j+extra[1],k+extra[2],n)) > 0) && (sgn<Real,int>((*this)(i-extra[0],j-extra[1],k-extra[2],n)) > 0))
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
                            if((sgn<Real,int>((*this)(i,j,k,n)) < 0) && (sgn<Real,int>((*this)(i+extra[0],j+extra[1],k+extra[2],n)) < 0) && (sgn<Real,int>((*this)(i-extra[0],j-extra[1],k-extra[2],n)) < 0))
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
    }

    return;
}

bool BoxAccessLevelSet::cellIsNextToAnInterface(int i, int j, int k, int n)
{
    Vector<int> pm {-1,0,1};

    //for(auto nk : pm)
    //{
        for(auto nj : pm)
        {
            for(auto ni : pm)
            {
                if(ni == 0 && nj == 0)// && nk == 0)
                {
                    continue;
                }
                else if( (sgn<Real,int>((*this)(i,j,k,n)) != sgn<Real,int>((*this)(i+ni,j+nj,k,n))) )
                {
                    return true;
                }
            }
        }
    //}

    return false;
}
