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
    case LAMBDA:            return (fab.array())(i, j, k, (accessPattern[m.var]+accessPattern.materialInfo[m.mat].mixtureIndex));
        break;
    case ALPHARHOLAMBDA:    return (fab.array())(i, j, k, (accessPattern[m.var]+accessPattern.materialInfo[m.mat].mixtureIndex));
        break;
    case SIGMA:             return (fab.array())(i, j, k, (accessPattern[m.var]+m.row*numberOfComponents+m.col));
        break;
    case V_TENSOR:          return (fab.array())(i, j, k, (accessPattern[m.var]+m.row*numberOfComponents+m.col));
        break;
    case VSTAR:             return (fab.array())(i, j, k, (accessPattern[m.var]+m.row*numberOfComponents+m.col));
        break;
    case DEVH:              return (fab.array())(i, j, k, (accessPattern[m.var]+m.row*numberOfComponents+m.col));
        break;
    case HJ2:               return (fab.array())(i, j, k, (accessPattern[m.var]));
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
                }

                getHenckyJ2(i,j,k);

                kineticEnergy = 0.0;

                for(int row = 0; row < numberOfComponents ; row++)
                {
                    (*this)(i,j,k,VELOCITY,0,row) = (*this)(i,j,k,RHOU,0,row)/(*this)(i,j,k,RHO);

                    kineticEnergy += 0.5*(*this)(i,j,k,RHO)*(*this)(i,j,k,VELOCITY,0,row)*(*this)(i,j,k,VELOCITY,0,row);
                }

                for(int m = 0; m < numberOfMaterials ; m++)
                {
                    if(accessPattern.materialInfo[m].mixture)
                    {
                        (*this)(i,j,k,LAMBDA,m)	  = (*this)(i,j,k,ALPHARHOLAMBDA,m)/(*this)(i,j,k,ALPHARHO,m);


                        accessPattern.materialInfo[m].EOS->rootFind((*this),i,j,k,m,kineticEnergy);
                    }
                }

                (*this)(i,j,k,P) = ((*this)(i,j,k,TOTAL_E)-kineticEnergy - getEffectiveNonThermalInternalEnergy(i,j,k)+ getEffectiveNonThermalPressure(i,j,k))/(getEffectiveInverseGruneisen(i,j,k));

                stressTensor(i,j,k);
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

                getHenckyJ2(i,j,k);

                kineticEnergy = 0.0;

                for(int row = 0; row < numberOfComponents ; row++)
                {
                    (*this)(i,j,k,RHOU,0,row) = (*this)(i,j,k,VELOCITY,0,row)*(*this)(i,j,k,RHO);

                    kineticEnergy += 0.5*(*this)(i,j,k,RHO)*(*this)(i,j,k,VELOCITY,0,row)*(*this)(i,j,k,VELOCITY,0,row);
                }

                (*this)(i,j,k,TOTAL_E) = (*this)(i,j,k,P)*getEffectiveInverseGruneisen(i,j,k) + getEffectiveNonThermalInternalEnergy(i,j,k) - getEffectiveNonThermalPressure(i,j,k) + kineticEnergy;

                stressTensor(i,j,k);
            }
        }
    }
}

void BoxAccessCellArray::stressTensor(int i, int j, int k)
{
        Real TotalShearModulus = 0.0;

        int checkSolid = 0;

        for(int m=0;m<numberOfMaterials;m++)
        {
            if(accessPattern.materialInfo[m].phase == solid)
            {
                checkSolid = 1;
            }

            TotalShearModulus += accessPattern.materialInfo[m].EOS->componentShearModulus((*this),i,j,k,m)*(accessPattern.materialInfo[m].EOS->inverseGruneisen((*this),i,j,k,m));
        }

        TotalShearModulus = TotalShearModulus/getEffectiveInverseGruneisen(i,j,k);

        if(checkSolid)
        {
            for(int row=0;row<numberOfComponents;row++)
            {
                for(int col=0;col<numberOfComponents;col++)
                {
                    (*this)(i,j,k,SIGMA,0,row,col)=-(*this)(i,j,k,P)*delta<Real>(row,col) + 2.0*TotalShearModulus*(*this)(i,j,k,DEVH,0,row,col);
                }
            }
        }
        else
        {
            for(int row=0;row<numberOfComponents;row++)
            {
                for(int col=0;col<numberOfComponents;col++)
                {
                    (*this)(i,j,k,SIGMA,0,row,col)=-(*this)(i,j,k,P)*delta<Real>(row,col);
                }
            }
        }

        return;
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

Real BoxAccessCellArray::getEffectiveNonThermalInternalEnergy(int i, int j, int k)
{
    Real tot = 0.0;

    for(int m=0;m<numberOfMaterials;m++)
    {
        tot += (accessPattern.materialInfo[m].EOS->coldCompressionInternalEnergy((*this),i,j,k,m)+accessPattern.materialInfo[m].EOS->shearInternalEnergy((*this),i,j,k,m));
    }

    return tot;
}

Real BoxAccessCellArray::getEffectiveNonThermalPressure(int i, int j, int k)
{
    Real tot = 0.0;

    for(int m=0;m<numberOfMaterials;m++)
    {
        tot += accessPattern.materialInfo[m].EOS->coldCompressionPressure((*this),i,j,k,m)+accessPattern.materialInfo[m].EOS->shearPressure((*this),i,j,k,m);
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
                Real a      = 0.0;
                Real xiTot  = 0.0;

                for(int m=0; m<numberOfMaterials;m++)
                {
                    a     += accessPattern.materialInfo[m].EOS->xi((*this),i,j,k,m)*accessPattern.materialInfo[m].EOS->getSoundSpeedContribution((*this),i,j,k,m)*(*this)(i,j,k,ALPHARHO,m)/(*this)(i,j,k,RHO);
                    xiTot += accessPattern.materialInfo[m].EOS->xi((*this),i,j,k,m)*(*this)(i,j,k,ALPHA,m);
                }

                (*this)(i,j,k,SOUNDSPEED) = std::sqrt(a/xiTot);


            }
        }
    }
}

Real BoxAccessCellArray::transverseWaveSpeed(int i, int j, int k)
{
    Real b      = 0.0;

    for(int m=0;m<numberOfMaterials;m++)
    {
        if(accessPattern.materialInfo[m].phase == solid)
        {
            b       += accessPattern.materialInfo[m].EOS->inverseGruneisen((*this),i,j,k,m)*accessPattern.materialInfo[m].EOS->componentShearModulus((*this),i,j,k,m)/((*this)(i,j,k,RHO_K,m));
        }
    }

    return sqrt(b/getEffectiveInverseGruneisen(i,j,k));
}

void BoxAccessCellArray::getHenckyJ2(int i, int j, int k)
{
    /*Real totalAlpha = 0.0;

    for(int m=0; m<numberOfMaterials;m++)
    {
        if(accessPattern.materialInfo[m].phase == solid)
        {
            totalAlpha += (*this)(i,j,k,ALPHA,m);
        }
    }*/

    for(int m=0; m<numberOfMaterials;m++)
    {
        if(accessPattern.materialInfo[m].phase == solid)
        {
            double tempdevH[numberOfComponents*numberOfComponents];

            getDeviatoricHenckyStrain(i,j,k);

            amrexToArray(i,j,k,DEVH,0,tempdevH);

            double product[numberOfComponents*numberOfComponents];

            (*this)(i,j,k,HJ2) = trace(squareMatrixMultiplyTranspose(tempdevH,tempdevH,product,numberOfComponents)); //totalAlpha

            return;
        }
    }

    return;
}

void BoxAccessCellArray::getDeviatoricHenckyStrain(int i, int j, int k)
{

    double temp[numberOfComponents*numberOfComponents];
    double tempV[numberOfComponents*numberOfComponents];

    amrexToArray(i,j,k,V_TENSOR,0,tempV);

    squareMatrixMultiplyTranspose(tempV,tempV,temp);

    invert(temp,tempV);

    for(int row=0;row<numberOfComponents;row++)
    {
        for(int col=0;col<numberOfComponents;col++)
        {
            (*this)(i,j,k,DEVH,0,row,col) = 0.5*0.5*(temp[row*numberOfComponents+col]-tempV[row*numberOfComponents+col]);
        }
    }

    return;

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

Real& BoxAccessCellArray::neighbour(int di, int dj, int dk, int i, int j, int k, Variable var, int mat, int row, int col)
{
    MaterialSpecifier temp(var,mat,row,col);

    return neighbour(di,dj,dk,i,j,k,temp);
}

Real& BoxAccessCellArray::neighbour(int di, int dj, int dk, int i, int j, int k, MaterialSpecifier& m)
{
    return (*this)(i+di,j+dj,k+dk,m);
}

void BoxAccessCellArray::amrexToArray(int i, int j, int k, Variable var, int m, double* copy, int nx, int ny)
{
    for(int row = 0; row<nx ;row++)
    {
        for(int col = 0; col<ny ; col++)
        {
            copy[row*nx+col] = (*this)(i,j,k,var,m,row,col);
        }
    }
}

void BoxAccessCellArray::normaliseV()
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    double temp[numberOfComponents*numberOfComponents];

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                amrexToArray(i,j,k,V_TENSOR,0,temp);

                Real norm = std::pow(det(temp,numberOfComponents),-1.0/3.0);

                for(int row=0;row<numberOfComponents;row++)
                {
                    for(int col=0;col<numberOfComponents;col++)
                    {
                        (*this)(i,j,k,V_TENSOR,0,row,col) *= norm;
                    }
                }
            }
        }
    }

    return;
}

void BoxAccessCellArray::cleanUpV()
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    double tempV[numberOfComponents*numberOfComponents];
    double temp [numberOfComponents*numberOfComponents];
    Real norm;
    Real totalSolid;
    Real totalAlpha;

    for    		   (int k = lo.z; k <= hi.z; ++k)
    {
        for        (int j = lo.y; j <= hi.y; ++j)
        {
            for    (int i = lo.x; i <= hi.x; ++i)
            {
                totalAlpha = 0.0;
                totalSolid = 0.0;

                for(int m = 0; m < numberOfMaterials; m++)
                {
                    totalAlpha += (*this)(i,j,k,ALPHA,m);

                    if(accessPattern.materialInfo[m].phase == solid)
                    {
                        totalSolid += (*this)(i,j,k,ALPHA,m);
                    }
                }

                //Print() << totalAlpha << std::endl;

                amrexToArray(i,j,k,V_TENSOR,0,tempV);

                squareMatrixMultiplyTranspose(tempV,tempV,temp);

                matrixSquareRoot(tempV,temp);

                norm = std::pow(det(tempV),-1.0/3.0);

                for(int row=0;row<numberOfComponents;row++)
                {
                    for(int col=0;col<numberOfComponents;col++)
                    {
                        (*this)(i,j,k,V_TENSOR,0,row,col) = ((tempV[row*numberOfComponents+col]*norm*totalSolid)+delta<Real>(row,col)*(totalAlpha-totalSolid))/totalAlpha;
                    }
                }
            }
        }
    }

    return;
}
