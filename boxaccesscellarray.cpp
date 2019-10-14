#include "cellarray.h"
#include "simulationheader.h"

BoxAccessCellArray::BoxAccessCellArray(const Box& bx, FArrayBox& fb,  CellArray& U) : box{bx}, fab{fb}, accessPattern{U.accessPattern}, numberOfMaterials{U.numberOfMaterials}{}

BoxAccessCellArray::BoxAccessCellArray(MFIter& mfi, const Box& bx, CellArray &U) : box{bx}, fab{U.data[mfi]}, accessPattern{U.accessPattern}, numberOfMaterials{U.numberOfMaterials}{}

Real& BoxAccessCellArray::operator()(int i, int j, int k, Variable var, int mat, int row, int col)
{
    MaterialSpecifier temp(var,mat,row,col);
    return (*this)(i,j,k,temp);
}

Real& BoxAccessCellArray::operator()(int i, int j, int k, int var, int mat, int row, int col)
{
    MaterialSpecifier temp((Variable)var,mat,row,col);
    return (*this)(i,j,k,temp);
}

Real& BoxAccessCellArray::operator()(int i, int j, int k, MaterialSpecifier& m)
{
    switch(m.var)
    {
    case ALPHA:             return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    case ALPHARHO:          return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    case RHO_K:             return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    case RHO:               return (fab.array())(i, j, k, (accessPattern[m.var]));
        break;
    case RHOU:              return (fab.array())(i, j, k, (accessPattern[m.var]+m.row));
        break;
    case TOTAL_E:           return (fab.array())(i, j, k, (accessPattern[m.var]));
        break;
    case VELOCITY:          return (fab.array())(i, j, k, (accessPattern[m.var]+m.row));
        break;
    case P:                 return (fab.array())(i, j, k, (accessPattern[m.var]));
        break;
    case SOUNDSPEED:        return (fab.array())(i, j, k, (accessPattern[m.var]));
        break;
    case USTAR:             return (fab.array())(i, j, k, (accessPattern[m.var]));
        break;
    case RHO_MIX:           return (fab.array())(i, j, k, (accessPattern[m.var]+accessPattern.materialInfo[m.mat].mixtureIndex+m.row));
        break;
    case LAMBDA:            return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    case ALPHARHOLAMBDA:    return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    default: Print() << "Incorrect Access variable " << m.var << std::endl;
        exit(1);
    }
}

void  BoxAccessCellArray::conservativeToPrimitive()
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Real kineticEnergy;

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                (*this)(i,j,k,RHO) = 0;

                for(int m = 0; m < numberOfMaterials ; m++)
                {
                    (*this)(i,j,k,RHO_K,m)	  = (*this)(i,j,k,ALPHARHO,m)/(*this)(i,j,k,ALPHA,m);
                    (*this)(i,j,k,RHO)       += (*this)(i,j,k,ALPHARHO,m);

                    if(accessPattern.materialInfo[m].mixture)
                    {
                        (*this)(i,j,k,LAMBDA,m)	  = (*this)(i,j,k,ALPHARHOLAMBDA,m)/(*this)(i,j,k,ALPHARHO,m);

                        accessPattern.materialInfo[m].EOS->rootFind((*this),i,j,k,m);
                    }
                }

                kineticEnergy = 0.0;

                for(int row = 0; row < AMREX_SPACEDIM ; row++)
                {
                    (*this)(i,j,k,VELOCITY,0,row) = (*this)(i,j,k,RHOU,0,row)/(*this)(i,j,k,RHO);

                    kineticEnergy += 0.5*(*this)(i,j,k,RHO)*(*this)(i,j,k,VELOCITY,0,row)*(*this)(i,j,k,VELOCITY,0,row);
                }

                (*this)(i,j,k,P) = ((*this)(i,j,k,TOTAL_E)-kineticEnergy)/(getEffectiveInverseGruneisen(i,j,k));
            }
        }
    }
}

void  BoxAccessCellArray::primitiveToConservative()
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Real kineticEnergy;

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                (*this)(i,j,k,RHO) = 0;

                for(int m = 0; m < numberOfMaterials ; m++)
                {
                    (*this)(i,j,k,ALPHARHO,m)	 = (*this)(i,j,k,ALPHA,m)*(*this)(i,j,k,RHO_K,m);
                    (*this)(i,j,k,RHO)          += (*this)(i,j,k,ALPHARHO,m);

                    if(accessPattern.materialInfo[m].mixture)
                    {
                        (*this)(i,j,k,ALPHARHOLAMBDA,m)  = (*this)(i,j,k,LAMBDA,m)*(*this)(i,j,k,ALPHA,m)*(*this)(i,j,k,RHO_K,m);
                    }
                }

                kineticEnergy = 0.0;

                for(int row = 0; row < AMREX_SPACEDIM ; row++)
                {
                    (*this)(i,j,k,RHOU,0,row) = (*this)(i,j,k,VELOCITY,0,row)*(*this)(i,j,k,RHO);

                    kineticEnergy += 0.5*(*this)(i,j,k,RHO)*(*this)(i,j,k,VELOCITY,0,row)*(*this)(i,j,k,VELOCITY,0,row);
                }


                (*this)(i,j,k,TOTAL_E)   	= (*this)(i,j,k,P)*getEffectiveInverseGruneisen(i,j,k) + kineticEnergy;
            }
        }
    }
}

Real BoxAccessCellArray::getEffectiveInverseGruneisen(int i, int j, int k)
{
    Real tot = 0.0;

    for(int m=0;m<numberOfMaterials;m++)
    {
        tot += (accessPattern.materialInfo[m].EOS->inverseGruneisen((*this),i,j,k,m));// (*this)(i,j,k,ALPHA,m)/(0.4);//  (accessPattern.materialInfo[m].EOS->inverseGruneisen((*this),i,j,k,m)); // (*this)(i,j,k,ALPHA,m)/(accessPattern.materialInfo[m].EOS->GruneisenGamma);
    }

    return tot;
}

void BoxAccessCellArray::getSoundSpeed()
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                Real a = 0.0;

                for(int m=0; m<numberOfMaterials;m++)
                {
                    a     += accessPattern.materialInfo[m].EOS->getSoundSpeedContribution((*this),i,j,k,m);// (((*this)(i,j,k,P)*(parameters.adiabaticIndex[m])/(*this)(i,j,k,RHO_K,m))*(*this)(i,j,k,ALPHARHO,m)/(parameters.adiabaticIndex[m]-1.0))/(*this)(i,j,k,RHO);
                }

                (*this)(i,j,k,SOUNDSPEED) = std::sqrt(a/getEffectiveInverseGruneisen(i,j,k));


            }
        }
    }
}

Real& BoxAccessCellArray::left(Direction_enum d, int i, int j, int k, MaterialSpecifier& m)
{
    switch(d)
    {
    case x: return (*this)(i-1,j,k,m);
        break;
    case y: return (*this)(i,j-1,k,m);
        break;
    case z: return (*this)(i,j,k-1,m);
        break;
    default: Print() << "Bad Direction in Left function" << std::endl; exit(1);
    }
}

Real& BoxAccessCellArray::right(Direction_enum d, int i, int j, int k, MaterialSpecifier& m)
{
    switch(d)
    {
    case x: return (*this)(i+1,j,k,m);
        break;
    case y: return (*this)(i,j+1,k,m);
        break;
    case z: return (*this)(i,j,k+1,m);
        break;
    default: Print() << "Bad Direction in Right function" << std::endl; exit(1);
    }
}

Real& BoxAccessCellArray::left(Direction_enum d, int i, int j, int k, Variable var, int mat, int row, int col)
{
    MaterialSpecifier temp(var,mat,row,col);
    return left(d,i,j,k,temp);
}

Real& BoxAccessCellArray::right(Direction_enum d, int i, int j, int k, Variable var, int mat, int row, int col)
{
    MaterialSpecifier temp(var,mat,row,col);
    return right(d,i,j,k,temp);
}
