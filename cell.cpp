#include "cell.h"

Cell::Cell(BoxAccessCellArray& U, int i, int j, int k)
{
    rhoU.resize(numberOfComponents);
    u.resize(numberOfComponents);

    materials.resize(U.numberOfMaterials);

    rho     =   (&U(i,j,k,RHO));
    rhoU[0] =   (&U(i,j,k,RHOU));
    E       =   (&U(i,j,k,TOTAL_E));
    u[0]    =   (&U(i,j,k,VELOCITY));
    p       =   (&U(i,j,k,P));
    a       =   (&U(i,j,k,SOUNDSPEED));

    for(int row = 0; row < numberOfComponents ; row++)
    {
        u[row]      = (&U(i,j,k,VELOCITY,0,row));
        rhoU[row]   = (&U(i,j,k,RHOU,0,row));
    }

    for(int m=0;m<U.numberOfMaterials;m++)
    {
        materials[m].alpha      = (&U(i,j,k,ALPHA,m));
        materials[m].alphaRho   = (&U(i,j,k,ALPHARHO,m));
        materials[m].rho        = (&U(i,j,k,RHO_K,m));
    }
}

//Cell::Cell(Array4<Real> const& arr, int i, int j, int k) : rho(&arr(i,j,k,RHO)), rhoU(&arr(i,j,k,RHOU)), E(&arr(i,j,k,TOTAL_E)), u(&arr(i,j,k,VELOCITY)), p(&arr(i,j,k,P)), a(&arr(i,j,k,SOUNDSPEED)){}

Real& Cell::operator()(Variable var, int mat, int row, int col)
{
    return (*this)(MaterialSpecifier(var,mat,row,col));
}

Real& Cell::operator()(MaterialSpecifier m)
{
    switch(m.var)
    {
    case ALPHA:     return *(materials[m.mat].alpha);
        break;
    case ALPHARHO:  return *(materials[m.mat].alphaRho);
        break;
    case RHO_K:     return *(materials[m.mat].rho);
        break;
    case RHO:       return *rho;
        break;
    case RHOU:      return *(rhoU[m.row]);
        break;
    case TOTAL_E:   return *E;
        break;
    case VELOCITY:  return *(u[m.row]);
        break;
    case P:         return *p;
        break;
    case SOUNDSPEED:return *a;
        break;
    default: Print() << "Incorrect cell variable" << std::endl;
        exit(1);
    }
}
