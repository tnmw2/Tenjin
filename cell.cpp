#include "cell.h"

Cell::Cell(BoxAccessCellArray& U, int i, int j, int k)
{
    rhoU.resize(numberOfComponents);
    u.resize(numberOfComponents);
    sigma.resize(numberOfComponents*numberOfComponents);

    materials.resize(U.numberOfMaterials);

    rho     =   (&U(i,j,k,RHO));
    E       =   (&U(i,j,k,TOTAL_E));
    p       =   (&U(i,j,k,P));
    a       =   (&U(i,j,k,SOUNDSPEED));

    for(int row = 0; row < U.numberOfComponents ; row++)
    {
        u[row]      = (&U(i,j,k,VELOCITY,0,row));
        rhoU[row]   = (&U(i,j,k,RHOU,0,row));

        for(int col = 0; col < U.numberOfComponents ; col++)
        {
            sigma[row*numberOfComponents+col] = (&U(i,j,k,SIGMA,0,row,col));
        }
    }

    for(int m=0;m<U.numberOfMaterials;m++)
    {
        materials[m].alpha      = (&U(i,j,k,ALPHA,m));
        materials[m].alphaRho   = (&U(i,j,k,ALPHARHO,m));
        materials[m].rho        = (&U(i,j,k,RHO_K,m));

        if(U.accessPattern.materialInfo[m].mixture)
        {
            materials[m].lambda         = (&U(i,j,k,LAMBDA,m));
            materials[m].alphaRhoLambda = (&U(i,j,k,ALPHARHOLAMBDA,m));
        }
    }

    uStar     =   (&U(i,j,k,USTAR));

}

Real& Cell::operator()(Variable var, int mat, int row, int col)
{
    return (*this)(MaterialSpecifier(var,mat,row,col));
}

Real& Cell::operator()(MaterialSpecifier m)
{
    switch(m.var)
    {
    case ALPHA:             return *(materials[m.mat].alpha);
        break;
    case ALPHARHO:          return *(materials[m.mat].alphaRho);
        break;
    case RHO_K:             return *(materials[m.mat].rho);
        break;
    case RHO:               return *rho;
        break;
    case RHOU:              return *(rhoU[m.row]);
        break;
    case TOTAL_E:           return *E;
        break;
    case VELOCITY:          return *(u[m.row]);
        break;
    case P:                 return *p;
        break;
    case SOUNDSPEED:        return *a;
        break;
    case USTAR:             return *uStar;
        break;
    case LAMBDA:            return *(materials[m.mat].lambda);
        break;
    case ALPHARHOLAMBDA:    return *(materials[m.mat].alphaRhoLambda);
        break;
    case SIGMA:             return *(sigma[m.row*numberOfComponents+m.col]);
    default: Print() << "Incorrect cell variable" << std::endl;
        exit(1);
    }
}
