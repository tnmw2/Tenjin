#include "simulationheader.h"

BoxAccessCellArray::BoxAccessCellArray(const Box& bx, FArrayBox& fb,  CellArray& U) : box{bx}, fab{fb}, accessPattern{U.accessPattern}, numberOfMaterials{U.numberOfMaterials}{}

BoxAccessCellArray::BoxAccessCellArray(MFIter& mfi, const Box& bx, CellArray &U) : box{bx}, fab{U.data[mfi]}, accessPattern{U.accessPattern}, numberOfMaterials{U.numberOfMaterials}{}

Real& BoxAccessCellArray::operator()(int i, int j, int k, Variable var, int mat, int row, int col)
{
    MaterialSpecifier temp(var,mat,row,col);
    return (*this)(i,j,k,temp);
}

const Real BoxAccessCellArray::operator()(int i, int j, int k, Variable var, int mat, int row, int col) const
{
    MaterialSpecifier temp(var,mat,row,col);
    return (*this)(i,j,k,temp);
}

Real& BoxAccessCellArray::operator()(int i, int j, int k, int var, int mat, int row, int col)
{
    MaterialSpecifier temp((Variable)var,mat,row,col);
    return (*this)(i,j,k,temp);
}

const Real BoxAccessCellArray::operator()(int i, int j, int k, MaterialSpecifier& m) const
{
    return (fab.array())(i, j, k, (accessPattern[m.var]+accessPattern.numberOfMaterialsForVariable[m.var]*m.mat+accessPattern.numberOfRowsForVariable[m.var]*m.row+m.col));
}

Real& BoxAccessCellArray::operator()(int i, int j, int k, MaterialSpecifier& m)
{
    return (fab.array())(i, j, k, (accessPattern[m.var]+accessPattern.numberOfMaterialsForVariable[m.var]*m.mat+accessPattern.numberOfRowsForVariable[m.var]*m.row+m.col));

}

void  BoxAccessCellArray::conservativeToPrimitive(int i, int j, int k)
{
    Real kineticEnergy = 0.0;

    (*this)(i,j,k,RHO) = 0.0;

    for(int m = 0; m < numberOfMaterials ; m++)
    {
        (*this)(i,j,k,RHO_K,m)	  = (*this)(i,j,k,ALPHARHO,m)/(*this)(i,j,k,ALPHA,m);
        (*this)(i,j,k,RHO)       += (*this)(i,j,k,ALPHARHO,m);

        if(accessPattern.materialInfo[m].plastic)
        {
            (*this)(i,j,k,EPSILON,m) = (*this)(i,j,k,ALPHARHOEPSILON,m)/(*this)(i,j,k,ALPHARHO,m);
        }

    }

    getHenckyJ2(i,j,k);

    for(int row = 0; row < numberOfComponents ; row++)
    {
        (*this)(i,j,k,VELOCITY,0,row) = (*this)(i,j,k,RHOU,0,row)/(*this)(i,j,k,RHO);

        (*this)(i,j,k,VELOCITY,0,row) = (std::abs((*this)(i,j,k,VELOCITY,0,row)) < 1E-10 ? 0.0 :(*this)(i,j,k,VELOCITY,0,row));

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

void  BoxAccessCellArray::primitiveToConservative(int i, int j, int k)
{

    Real kineticEnergy = 0.0;

    (*this)(i,j,k,RHO) = 0.0;

    for(int m = 0; m < numberOfMaterials ; m++)
    {
        (*this)(i,j,k,ALPHARHO,m)	 = (*this)(i,j,k,ALPHA,m)*(*this)(i,j,k,RHO_K,m);
        (*this)(i,j,k,RHO)          += (*this)(i,j,k,ALPHARHO,m);

        if(accessPattern.materialInfo[m].mixture)
        {
            (*this)(i,j,k,ALPHARHOLAMBDA,m)  = (*this)(i,j,k,LAMBDA,m)*(*this)(i,j,k,ALPHA,m)*(*this)(i,j,k,RHO_K,m);

            accessPattern.materialInfo[m].EOS->defineMixtureDensities((*this),i,j,k,m);
        }

        if(accessPattern.materialInfo[m].plastic)
        {
            (*this)(i,j,k,ALPHARHOEPSILON,m) = (*this)(i,j,k,EPSILON,m)*(*this)(i,j,k,ALPHARHO,m);
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

void  BoxAccessCellArray::conservativeToPrimitive()
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    checkLimits(accessPattern.conservativeVariables);

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                conservativeToPrimitive(i,j,k);
            }
        }
    }

    checkLimits(accessPattern.primitiveVariables);
}

void  BoxAccessCellArray::primitiveToConservative()
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    checkLimits(accessPattern.primitiveVariables);

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                primitiveToConservative(i,j,k);
            }
        }
    }

    checkLimits(accessPattern.conservativeVariables);

}

void BoxAccessCellArray::stressTensor(int i, int j, int k)
{
    if(accessPattern.parameters.SOLID)
    {
        Real TotalShearModulus = 0.0;

        for(int m=0;m<numberOfMaterials;m++)
        {
            TotalShearModulus += accessPattern.materialInfo[m].EOS->componentShearModulus((*this),i,j,k,m)*(accessPattern.materialInfo[m].EOS->inverseGruneisen((*this),i,j,k,m));
        }

        TotalShearModulus = TotalShearModulus/getEffectiveInverseGruneisen(i,j,k);

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
                getSoundSpeed(i,j,k);
            }
        }
    }
}

void BoxAccessCellArray::getSoundSpeed(int i, int j, int k)
{

    Real a      = 0.0;
    Real xiTot  = 0.0;

    for(int m=0; m<numberOfMaterials;m++)
    {
        a     += std::max(0.0,accessPattern.materialInfo[m].EOS->xi((*this),i,j,k,m)*accessPattern.materialInfo[m].EOS->getSoundSpeedContribution((*this),i,j,k,m)*(*this)(i,j,k,ALPHARHO,m)/(*this)(i,j,k,RHO));
        xiTot += accessPattern.materialInfo[m].EOS->xi((*this),i,j,k,m)*(*this)(i,j,k,ALPHA,m);
    }

    if(a<=0.0)
    {
        a = 1E-10;// soundSpeedTolerance
    }

    (*this)(i,j,k,SOUNDSPEED) = std::sqrt(a/xiTot);

}

Real BoxAccessCellArray::transverseWaveSpeed(int i, int j, int k)
{
    Real b      = 0.0;

    getHenckyJ2(i,j,k);

    for(int m=0;m<numberOfMaterials;m++)
    {
        if(accessPattern.materialInfo[m].phase == solid)
        {
            b       += accessPattern.materialInfo[m].EOS->inverseGruneisen((*this),i,j,k,m)*accessPattern.materialInfo[m].EOS->componentShearModulus((*this),i,j,k,m)/((*this)(i,j,k,RHO_K,m));
        }
    }

    if(std::isnan(b) || b < 0.0)
    {
        b = 0.0;
    }

    Real temp = sqrt(b/getEffectiveInverseGruneisen(i,j,k));

    if(std::isnan(temp))
    {
        temp = 0.0;
    }

    return temp;

    //return sqrt(b/getEffectiveInverseGruneisen(i,j,k));
}

void BoxAccessCellArray::getHenckyJ2(int i, int j, int k)
{
    for(int m=0; m<numberOfMaterials;m++)
    {
        if(accessPattern.materialInfo[m].phase == solid)
        {
            double tempdevH[numberOfComponents*numberOfComponents];
            double tempdevHT[numberOfComponents*numberOfComponents];
            double product[numberOfComponents*numberOfComponents];


            getDeviatoricHenckyStrain(i,j,k);

            amrexToArray(i,j,k,DEVH,0,tempdevH);

            matrixCopy(tempdevH,tempdevHT);

            (*this)(i,j,k,HJ2) = trace(squareMatrixMultiplyTranspose(tempdevH,tempdevHT,product,numberOfComponents)); //totalAlpha

            return;
        }
    }

    return;
}

void BoxAccessCellArray::getDeviatoricHenckyStrain(int i, int j, int k)
{

    double temp[numberOfComponents*numberOfComponents];
    double tempV[numberOfComponents*numberOfComponents];
    double tempV1[numberOfComponents*numberOfComponents];

    amrexToArray(i,j,k,V_TENSOR,0,tempV);

    matrixCopy(tempV,tempV1);

    squareMatrixMultiplyTranspose(tempV,tempV1,temp);

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

void BoxAccessCellArray::amrexToArray(int i, int j, int k, Variable var, int m, double* copy, int nx, int ny) const
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

                /*if(totalSolid < 0.0001)
                {
                    for(int row=0;row<numberOfComponents;row++)
                    {
                        for(int col=0;col<numberOfComponents;col++)
                        {
                            (*this)(i,j,k,V_TENSOR,0,row,col) = delta<Real>(row,col);
                        }
                    }
                }
                else
                {*/
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
                //}

                getHenckyJ2(i,j,k);
            }
        }
    }

    return;
}

void BoxAccessCellArray::cleanUpAlpha(const Real* dx, const Real* prob_lo)
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Real totalAlpha;
    int Nan;

    for 		(int k = lo.z; k <= hi.z; ++k)
    {
        for 	(int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {

                totalAlpha = 0.0;
                Nan = 0;

                for(int m=0;m<numberOfMaterials;m++)
                {
                    if(std::isnan((*this)(i,j,k,ALPHA,m)))
                    {
                        (*this)(i,j,k,ALPHA,m) = 0.0;

                        if(Nan == 0)
                        {
                            Nan = m;
                        }
                        else
                        {

                            Vector<Real> err;

                            err.push_back(prob_lo[0] + (Real(i)+0.5)*dx[0]);
                            err.push_back(prob_lo[1] + (Real(j)+0.5)*dx[1]);

                            for(int m=0;m<numberOfMaterials;m++)
                            {
                                err.push_back((*this)(i,j,k,ALPHA,m));
                            }

                            std::string message = "Error in scaling volume fractions, too many Nans. at : ";

                            customAbort(err,message);
                        }
                    }
                    else if((*this)(i,j,k,ALPHA,m)<0.0)
                    {
                        (*this)(i,j,k,ALPHA,m)= 1E-6;
                    }
                    else if((*this)(i,j,k,ALPHA,m)>1.0)
                    {
                        (*this)(i,j,k,ALPHA,m) = 1.0;
                    }

                    totalAlpha += (*this)(i,j,k,ALPHA,m);
                }

                if(Nan>0)
                {
                    (*this)(i,j,k,ALPHA,Nan) = 1.0 - totalAlpha;

                    totalAlpha += (*this)(i,j,k,ALPHA,Nan);
                }

                if(totalAlpha <= 0.0)
                {
                    Abort("Error in scaling volume fractions because total is less than zero.");
                }

                for(int m=0;m<numberOfMaterials;m++)
                {
                    (*this)(i,j,k,ALPHA,m) *= 1.0/totalAlpha;
                }
            }
        }
    }

    return;
}

bool BoxAccessCellArray::check(MaterialSpecifier& m)
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                if(std::isnan((*this)(i,j,k,m)))
                {
                    return true;
                }
            }
        }
    }

    return false;
}

bool BoxAccessCellArray::contains_nan()
{
    bool checker = false;

    for(auto n : accessPattern.conservativeVariables)
    {
        if(check(n))
        {
            Print() << "Nan in variable: " << accessPattern.variableNames[accessPattern[n.var]] << std::endl;
            checker = true;
        }
    }
    for(auto n : accessPattern.primitiveVariables)
    {
        if(check(n))
        {
            Print() << "Nan in variable: " << accessPattern.variableNames[accessPattern[n.var]] << std::endl;
            checker = true;
        }
    }

    return checker;
}

bool BoxAccessCellArray::cellIsMostlyFluid(int i, int j, int k)
{
    double TotalSolid = 0.0;

    for(int m=0;m<numberOfMaterials;m++)
    {
        if(accessPattern.materialInfo[m].phase == solid )
        {
            TotalSolid += (*this)(i,j,k,ALPHA,m);
        }
    }

    return TotalSolid <= 0.5 ;
}

void BoxAccessCellArray::checkAndAmendVariable(MaterialSpecifier& m, int i, int j, int k)
{
    if(std::isnan((*this)(i,j,k,m)))
    {
        smoothFromNeighbours(m,i,j,k);
    }
    else if((*this)(i,j,k,m) < accessPattern.limits[m.var].first)
    {
        //smoothFromNeighbours(m,i,j,k);

        (*this)(i,j,k,m) = accessPattern.limits[m.var].first;

        return;
    }
    else if((*this)(i,j,k,m) > accessPattern.limits[m.var].second)
    {
        //smoothFromNeighbours(m,i,j,k);

        (*this)(i,j,k,m) = accessPattern.limits[m.var].second;

        return;
    }

    return;
}

void BoxAccessCellArray::checkLimits(Vector<MaterialSpecifier>& list)
{
    if(!checking)
    {
        return;
    }

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for(auto n : list)
    {
        for    		(int k = lo.z; k <= hi.z; ++k)
        {
            for     (int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    checkAndAmendVariable(n,i,j,k);
                }
            }
        }
    }
}

void BoxAccessCellArray::checkLimits(MaterialSpecifier& m)
{
    if(!checking)
    {
        return;
    }

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                checkAndAmendVariable(m,i,j,k);
            }
        }
    }
}

void BoxAccessCellArray::smoothFromNeighbours(MaterialSpecifier& m, int i, int j, int k)
{
    Real sum     = 0.0;
    int  counter = 0;

    Vector<int> pm {-1,0,1};

    for(auto nk : pm)
    {
        for(auto nj : pm)
        {
            for(auto ni : pm)
            {
                if(ni == 0 && nj == 0 && nk == 0)
                {
                    continue;
                }
                else if(!(std::isnan(neighbour(ni,nj,nk,i,j,k,m))))
                {
                    sum += neighbour(ni,nj,nk,i,j,k,m);

                    counter++;
                }
            }
        }
    }

    if(counter > 0)
    {
        (*this)(i,j,k,m) = (1.0/counter)*(sum);
    }
    else
    {
        amrex::Abort("Couldn't smooth from Neighbours");
    }
}
