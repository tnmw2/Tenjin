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
    return (fab.array())(i, j, k, (accessPattern[m.var]+accessPattern.numberOfMaterialsForVariable[m.var]*m.mat+accessPattern.numberOfRowsForVariable[m.var]*m.row+m.col));


    /*switch(m.var)
    {
    case RHO:               return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    case RHOU:              return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat*numberOfComponents+m.row));
        break;
    case TOTAL_E:           return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    case VELOCITY:          return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat*numberOfComponents+m.row));
        break;
    case P:                 return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    case SOUNDSPEED:        return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    case USTAR:             return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    case SIGMA:             return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat*numberOfComponents*numberOfComponents+m.row*numberOfComponents+m.col));
        break;
    case V_TENSOR:          return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat*numberOfComponents*numberOfComponents+m.row*numberOfComponents+m.col));
        break;
    case DEVH:              return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat*numberOfComponents*numberOfComponents+m.row*numberOfComponents+m.col));
        break;
    case VSTAR:             return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat*numberOfComponents*numberOfComponents+m.row*numberOfComponents+m.col));
        break;
    case HJ2:               return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    case EPSILON:           return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    case ALPHARHOEPSILON:   return (fab.array())(i, j, k, (accessPattern[m.var]+m.mat));
        break;
    default: Print() << "Incorrect Access variable " << m.var << " in boxaccesscellarray: " << accessPattern.variableNames[getArrayPosition(m)] << std::endl;
        exit(1);
    }*/
}

int BoxAccessCellArray::getArrayPosition(MaterialSpecifier& m)
{
   // Note: If you change this, change the corresponding function for cellarray.

    return (accessPattern[m.var]+accessPattern.numberOfMaterialsForVariable[m.var]*m.mat+accessPattern.numberOfRowsForVariable[m.var]*m.row+m.col);


    /*switch(m.var)
    {
        case RHO:               return (accessPattern[m.var]+m.mat);
            break;
        case RHOU:              return (accessPattern[m.var]+m.mat*numberOfComponents+m.row);
            break;
        case TOTAL_E:           return (accessPattern[m.var]+m.mat);
            break;
        case VELOCITY:          return (accessPattern[m.var]+m.mat*numberOfComponents+m.row);
            break;
        case P:                 return (accessPattern[m.var]+m.mat);
            break;
        case SOUNDSPEED:        return (accessPattern[m.var]+m.mat);
            break;
        case USTAR:             return (accessPattern[m.var]+m.mat);
            break;
        case SIGMA:             return (accessPattern[m.var]+m.mat*numberOfComponents*numberOfComponents+m.row*numberOfComponents+m.col);
            break;
        case V_TENSOR:          return (accessPattern[m.var]+m.mat*numberOfComponents*numberOfComponents+m.row*numberOfComponents+m.col);
            break;
        case DEVH:              return (accessPattern[m.var]+m.mat*numberOfComponents*numberOfComponents+m.row*numberOfComponents+m.col);
            break;
        case VSTAR:             return (accessPattern[m.var]+m.mat*numberOfComponents*numberOfComponents+m.row*numberOfComponents+m.col);
            break;
        case HJ2:               return (accessPattern[m.var]+m.mat);
            break;
        case EPSILON:           return (accessPattern[m.var]+m.mat);
            break;
        case ALPHARHOEPSILON:   return (accessPattern[m.var]+m.mat);
            break;
    default: Print() << "Incorrect Access variable " << m.var << " in boxaccesscellarray: " << accessPattern.variableNames[getArrayPosition(m)] << std::endl;
        exit(1);
    }*/
}

void  BoxAccessCellArray::conservativeToPrimitive()
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    checkLimits(accessPattern.conservativeVariables);

    Real kineticEnergy;

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                (*this).conservativeToPrimitive(i,j,k);
            }
        }
    }

    checkLimits(accessPattern.primitiveVariables);

}

void  BoxAccessCellArray::primitiveToConservative()
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Real kineticEnergy;

    checkLimits(accessPattern.primitiveVariables);

    for    		(int k = lo.z; k <= hi.z; ++k)
    {
        for     (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                (*this).primitiveToConservative(i,j,k);
            }
        }
    }

    checkLimits(accessPattern.conservativeVariables);

}

void  BoxAccessCellArray::conservativeToPrimitive(int i, int j, int k)
{
    Real kineticEnergy;

    for(int m = 0; m < numberOfMaterials ; m++)
    {
        if(accessPattern.materialInfo[m].plastic)
        {
            (*this)(i,j,k,EPSILON,m) = (*this)(i,j,k,RHOEPSILON,m)/(*this)(i,j,k,RHO,m);
        }

        kineticEnergy = 0.0;

        for(int row = 0; row < numberOfComponents ; row++)
        {
            (*this)(i,j,k,VELOCITY,m,row) = (*this)(i,j,k,RHOU,m,row)/(*this)(i,j,k,RHO,m);

            kineticEnergy += 0.5*(*this)(i,j,k,RHO,m)*(*this)(i,j,k,VELOCITY,m,row)*(*this)(i,j,k,VELOCITY,m,row);
        }

        getHenckyJ2(i,j,k,m);

        (*this)(i,j,k,P,m) = ((*this)(i,j,k,TOTAL_E,m)-kineticEnergy - getEffectiveNonThermalInternalEnergy(i,j,k,m)+ getEffectiveNonThermalPressure(i,j,k,m))/(getEffectiveInverseGruneisen(i,j,k,m));

        stressTensor(i,j,k,m);
    }

    return;
}

void  BoxAccessCellArray::primitiveToConservative(int i, int j, int k)
{
    Real kineticEnergy;

    for(int m = 0; m < numberOfMaterials ; m++)
    {
        if(accessPattern.materialInfo[m].plastic)
        {
            (*this)(i,j,k,RHOEPSILON,m) = (*this)(i,j,k,EPSILON,m)*(*this)(i,j,k,RHO,m);
        }

        kineticEnergy = 0.0;

        for(int row = 0; row < numberOfComponents ; row++)
        {
            (*this)(i,j,k,RHOU,m,row) = (*this)(i,j,k,VELOCITY,m,row)*(*this)(i,j,k,RHO,m);

            kineticEnergy += 0.5*(*this)(i,j,k,RHO,m)*(*this)(i,j,k,VELOCITY,m,row)*(*this)(i,j,k,VELOCITY,m,row);
        }

        getHenckyJ2(i,j,k,m);

        (*this)(i,j,k,TOTAL_E,m) = (*this)(i,j,k,P,m)*getEffectiveInverseGruneisen(i,j,k,m) + getEffectiveNonThermalInternalEnergy(i,j,k,m) - getEffectiveNonThermalPressure(i,j,k,m) + kineticEnergy;

        stressTensor(i,j,k,m);
    }
}

void BoxAccessCellArray::stressTensor(int i, int j, int k, int m)
{
    if(accessPattern.materialInfo[m].phase == solid)
    {
        Real shearModulus = accessPattern.materialInfo[m].EOS->componentShearModulus((*this),i,j,k,m);

        for(int row=0;row<numberOfComponents;row++)
        {
            for(int col=0;col<numberOfComponents;col++)
            {
                (*this)(i,j,k,SIGMA,m,row,col)=-(*this)(i,j,k,P,m)*delta<Real>(row,col) + 2.0*shearModulus*(*this)(i,j,k,DEVH,m,row,col);
            }
        }
    }
    else
    {
        for(int row=0;row<numberOfComponents;row++)
        {
            for(int col=0;col<numberOfComponents;col++)
            {
                (*this)(i,j,k,SIGMA,m,row,col)=-(*this)(i,j,k,P,m)*delta<Real>(row,col);
            }
        }
    }

    return;
}

Real BoxAccessCellArray::getEffectiveInverseGruneisen(int i, int j, int k, int m)
{
    return (accessPattern.materialInfo[m].EOS->inverseGruneisen((*this),i,j,k,m));
}

Real BoxAccessCellArray::getEffectiveNonThermalInternalEnergy(int i, int j, int k, int m)
{
    return (accessPattern.materialInfo[m].EOS->coldCompressionInternalEnergy((*this),i,j,k,m)+accessPattern.materialInfo[m].EOS->shearInternalEnergy((*this),i,j,k,m));
}

Real BoxAccessCellArray::getEffectiveNonThermalPressure(int i, int j, int k, int m)
{
    return accessPattern.materialInfo[m].EOS->coldCompressionPressure((*this),i,j,k,m)+accessPattern.materialInfo[m].EOS->shearPressure((*this),i,j,k,m);
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
                for(int m=0; m<numberOfMaterials;m++)
                {
                    (*this)(i,j,k,SOUNDSPEED,m) = sqrt(std::max(1E-10,accessPattern.materialInfo[m].EOS->getSoundSpeedContribution((*this),i,j,k,m)));
                }
            }
        }
    }
}

void BoxAccessCellArray::getSoundSpeed(int i, int j, int k)
{
    for(int m=0; m<numberOfMaterials;m++)
    {
        (*this)(i,j,k,SOUNDSPEED,m) = sqrt(std::max(1E-10,accessPattern.materialInfo[m].EOS->getSoundSpeedContribution((*this),i,j,k,m)));
    }
}

Real BoxAccessCellArray::transverseWaveSpeed(int i, int j, int k, int m)
{
    if(accessPattern.materialInfo[m].phase == solid)
    {
        getHenckyJ2(i,j,k,m);

        Real b  = accessPattern.materialInfo[m].EOS->componentShearModulus((*this),i,j,k,m)/((*this)(i,j,k,RHO,m));

        if(std::isnan(b) || b < 0.0)
        {
            b = 0.0;
        }

        Real temp = sqrt(b);

        if(std::isnan(temp))
        {
            temp = 0.0;
        }

        return temp;
    }
    else
    {
        return 0.0;
    }
}

void BoxAccessCellArray::getHenckyJ2(int i, int j, int k, int m)
{
    if(accessPattern.materialInfo[m].phase == solid)
    {
        double tempdevH[numberOfComponents*numberOfComponents];
        double tempdevHT[numberOfComponents*numberOfComponents];
        double product[numberOfComponents*numberOfComponents];

        getDeviatoricHenckyStrain(i,j,k,m);

        amrexToArray(i,j,k,DEVH,m,tempdevH);

        matrixCopy(tempdevH,tempdevHT);

        (*this)(i,j,k,HJ2,m) = trace(squareMatrixMultiplyTranspose(tempdevH,tempdevHT,product,numberOfComponents));

        return;
    }

    return;
}

void BoxAccessCellArray::getDeviatoricHenckyStrain(int i, int j, int k, int m)
{

    double temp[numberOfComponents*numberOfComponents];
    double tempV[numberOfComponents*numberOfComponents];
    double tempV1[numberOfComponents*numberOfComponents];

    amrexToArray(i,j,k,V_TENSOR,m,tempV);

    matrixCopy(tempV,tempV1);

    squareMatrixMultiplyTranspose(tempV,tempV1,temp);

    invert(temp,tempV);

    for(int row=0;row<numberOfComponents;row++)
    {
        for(int col=0;col<numberOfComponents;col++)
        {
            (*this)(i,j,k,DEVH,m,row,col) = 0.5*0.5*(temp[row*numberOfComponents+col]-tempV[row*numberOfComponents+col]);
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

    for             (int m = 0;    m <  numberOfMaterials; m++)
    {
        if(accessPattern.materialInfo[m].phase == solid)
        {
            for    		(int k = lo.z; k <= hi.z; ++k)
            {
                for     (int j = lo.y; j <= hi.y; ++j)
                {
                    for (int i = lo.x; i <= hi.x; ++i)
                    {
                        amrexToArray(i,j,k,V_TENSOR,m,temp);

                        Real norm = std::pow(det(temp,numberOfComponents),-1.0/3.0);

                        for(int row=0;row<numberOfComponents;row++)
                        {
                            for(int col=0;col<numberOfComponents;col++)
                            {
                                (*this)(i,j,k,V_TENSOR,m,row,col) *= norm;
                            }
                        }
                    }
                }
            }
        }
    }

    return;
}

void BoxAccessCellArray::cleanUpV()
{
    normaliseV();

    return;
}

void BoxAccessCellArray::cleanUpAlpha()
{
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

void BoxAccessCellArray::bilinearInterpolation(BoxAccessCellArray& U, int i, int j, int k, const Real* dx, const Real* prob_lo, Real probe_x, Real probe_y, int m,Vector<int>& corner_x, Vector<int>& corner_y, Vector<Real>& xdiff, Vector<Real>& ydiff, Vector<Real>& B)
{
    //Print() << "HERE" << std::endl;

    int incell_x,incell_y;

    if((probe_x-prob_lo[0]) > 0.0)
    {
       incell_x = (int)((probe_x-prob_lo[0])/dx[0]);

       if(probe_x > (prob_lo[0] + (Real(incell_x)+0.5)*dx[0]))
       {
           corner_x[0] = incell_x;
           corner_x[1] = corner_x[0]+1;
       }
       else
       {
           corner_x[0] = incell_x-1;
           corner_x[1] = corner_x[0]+1;
       }
    }
    else
    {
       incell_x = (int)((probe_x-prob_lo[0])/dx[0]) -1;

       if(probe_x > (prob_lo[0] + (Real(incell_x)+0.5)*dx[0]))
       {
           corner_x[0] = incell_x;
           corner_x[1] = corner_x[0]+1;
       }
       else
       {
           corner_x[0] = incell_x-1;
           corner_x[1] = corner_x[0]+1;
       }
    }

    if((probe_y-prob_lo[1]) > 0.0)
    {
       incell_y = (int)((probe_y-prob_lo[1])/dx[1]);

       if(probe_y > (prob_lo[1] + (Real(incell_y)+0.5)*dx[1]))
       {
           corner_y[0] = incell_y;
           corner_y[1] = corner_y[0]+1;
       }
       else
       {
           corner_y[0] = incell_y-1;
           corner_y[1] = corner_y[0]+1;
       }
    }
    else
    {
       incell_y = (int)((probe_y-prob_lo[1])/dx[1]) -1;

       if(probe_y > (prob_lo[1] + (Real(incell_y)+0.5)*dx[1]))
       {
           corner_y[0] = incell_y;
           corner_y[1] = corner_y[0]+1;
       }
       else
       {
           corner_y[0] = incell_y-1;
           corner_y[1] = corner_y[0]+1;
       }
    }

    //Print()<<"corners_x:\t " << corner_x[0] << " " << corner_x[1] << std::endl;
    //Print()<<"corners_y:\t " << corner_y[0] << " " << corner_y[1] << std::endl;

    xdiff[0] = prob_lo[0] + (Real(corner_x[1])+0.5)*dx[0] - probe_x;
    xdiff[1] = prob_lo[0] + (Real(corner_x[0])+0.5)*dx[0] - probe_x;

    ydiff[0] = prob_lo[1] + (Real(corner_y[1])+0.5)*dx[1] - probe_y;
    ydiff[1] = prob_lo[1] + (Real(corner_y[0])+0.5)*dx[1] - probe_y;

    //Print()<<"xdiff:\t " << xdiff[0] << " " << xdiff[1] << std::endl;
    //Print()<<"ydiff:\t " << ydiff[0] << " " << ydiff[1] << std::endl;

    if((sgn<Real,int>(xdiff[0]) == sgn<Real,int>(xdiff[1])) || (sgn<Real,int>(ydiff[0]) == sgn<Real,int>(ydiff[1])) )
    {
        Print()<<"xdiff:\t "  << xdiff[0] << " " << xdiff[1] << std::endl;
        Print()<<"ydiff:\t "  << ydiff[0] << " " << ydiff[1] << std::endl;

        Print()<<"corners_x:\t "  << corner_x[0] << " " << corner_x[1] << std::endl;
        Print()<<"corners_y:\t "  << corner_y[0] << " " << corner_y[1] << std::endl;

        Print()<<"incell:\t "  << incell_x << " " << incell_y << std::endl;
        Print()<<"probe :\t "  << probe_x  << " " << probe_y  << std::endl;

        Print() << dx[0] << std::endl;



        //Print() << dx[1] << std::endl;
        //
        //Print() << "y " << probe_y << " " << prob_lo[1] + (Real(corner_y[1])+0.5)*dx[1] << " " << prob_lo[1] + (Real(corner_y[0])+0.5)*dx[1] << std::endl;
        //Print() <<  remainder((probe_y-prob_lo[1])/dx[1],1.0) << std::endl;
        //
        //Print() << probe_y << " " << prob_lo[1] << " " << dx[1] << " " <<  (probe_y-prob_lo[1])/dx[1] << " so  " <<  (int)((probe_y-prob_lo[1])/dx[1]) << " " << corner_y[1]<< std::endl;
        Abort("Sign error in interpolation");
    }


    for(auto n : accessPattern.material_primitiveVariables[m])
    {
        //Print() << "B" << std::endl;

        for(int row = 0; row< 2;row++)
        {
            for(int col = 0; col < 2; col++)
            {
                if(row == col)
                {
                    B[row*2+col] = U(corner_x[row],corner_y[col],k,n);
                }
                else
                {
                    B[row*2+col] = -U(corner_x[row],corner_y[col],k,n);
                }

                //Print() << B[row*2+col]<< " ";
            }

            //Print() << std::endl;
        }

        (*this)(i,j,k,n) = 0.0;

        for(int row = 0; row< 2;row++)
        {
            for(int col = 0; col < 2; col++)
            {
                (*this)(i,j,k,n) +=  B[row*2+col]*xdiff[row]*ydiff[col];
            }
        }

        (*this)(i,j,k,n) *= 1.0/(dx[0]*dx[1]);


        //Print() <<accessPattern.variableNames[accessPattern[n.var]] << " " << (*this)(i,j,k,n) << std::endl;
    }

}

void BoxAccessCellArray::rotateFrameSoXPointsAlongNormal(int i, int j, int k, Real nx, Real ny, int m)
{
    //Just a rotation matrix (nx,ny //-ny,nx)

    Real normalVel  = (*this)(i,j,k,VELOCITY,m,x)* nx +(*this)(i,j,k,VELOCITY,m,y)*ny;
    Real tangentVel = (*this)(i,j,k,VELOCITY,m,x)*-ny +(*this)(i,j,k,VELOCITY,m,y)*nx;

    (*this)(i,j,k,VELOCITY,m,x) = normalVel;
    (*this)(i,j,k,VELOCITY,m,y) = tangentVel;

    if(accessPattern.materialInfo[m].phase == solid)
    {
        Real rotationMatrix[numberOfComponents*numberOfComponents] = {nx, ny, 0.0, -ny, nx, 0.0, 0.0, 0.0, 1.0};
        Real Vmat          [numberOfComponents*numberOfComponents];
        Real prod          [numberOfComponents*numberOfComponents];

        amrexToArray(i,j,k,V_TENSOR,m,Vmat);

        squareMatrixMultiply(rotationMatrix,Vmat,prod);

        for(int row = 0; row< numberOfComponents; row++)
        {
            for(int col = 0; col < numberOfComponents; col++)
            {
                (*this)(i,j,k,V_TENSOR,m,row,col) = prod[row*numberOfComponents+col];
            }
        }
    }

    return;
}

void BoxAccessCellArray::rotateFrameBack(int i, int j, int k, Real nx, Real ny)
{
    //Just a rotation matrix (nx,-ny //ny,nx)

    for(int m = 0; m<numberOfMaterials;m++)
    {
        Real normalVel  = (*this)(i,j,k,VELOCITY,m,0)* nx +(*this)(i,j,k,VELOCITY,m,1)*-ny;
        Real tangentVel = (*this)(i,j,k,VELOCITY,m,0)* ny +(*this)(i,j,k,VELOCITY,m,1)* nx;

        (*this)(i,j,k,VELOCITY,m,0) = normalVel;
        (*this)(i,j,k,VELOCITY,m,1) = tangentVel;

        if(accessPattern.materialInfo[m].phase == solid)
        {
            Real rotationMatrix[numberOfComponents*numberOfComponents] = {nx, -ny, 0.0, ny, nx, 0.0, 0.0, 0.0, 1.0};
            Real Vmat          [numberOfComponents*numberOfComponents];
            Real prod          [numberOfComponents*numberOfComponents];

            amrexToArray(i,j,k,V_TENSOR,m,Vmat);

            squareMatrixMultiply(rotationMatrix,Vmat,prod);

            for(int row = 0; row< numberOfComponents; row++)
            {
                for(int col = 0; col < numberOfComponents; col++)
                {
                    (*this)(i,j,k,V_TENSOR,m,row,col) = prod[row*numberOfComponents+col];
                }
            }
        }
    }

    return;
}

void BoxAccessCellArray::realPositionToCell(int& i, int& j, int& k, Real x, Real y, Real z, const Real* dx, const Real* prob_lo)
{

    if((x-prob_lo[0]) > 0.0)
    {
       i = (int)((x-prob_lo[0])/dx[0]);
    }
    else
    {
       j = (int)((x-prob_lo[0])/dx[0])-1;
    }

    if((y-prob_lo[1]) > 0.0)
    {
       j = (int)((y-prob_lo[1])/dx[1]);
    }
    else
    {
       j = (int)((y-prob_lo[1])/dx[1])-1;
    }

    if(AMREX_SPACEDIM == 3)
    {
        if((z-prob_lo[2]) > 0.0)
        {
           k = (int)((z-prob_lo[2])/dx[2]);
        }
        else
        {
           k = (int)((z-prob_lo[2])/dx[2])-1;
        }
    }

}
