#include "cell.h"

Cell::Cell(BoxAccessCellArray& U, int i, int j, int k, Material_type _phase) : phase(_phase), accessPattern(U.accessPattern), parent_i(i), parent_j(j), parent_k(k), parent(&U)
{
    rhoU.resize(numberOfComponents);
    u.resize(numberOfComponents);
    sigma.resize(numberOfComponents*numberOfComponents);

    if(phase == solid)
    {
        V.resize(numberOfComponents*numberOfComponents);
        VStar.resize(numberOfComponents*numberOfComponents);
    }

    materials.resize(U.numberOfMaterials);

    for(auto n : accessPattern.cellVariables)
    {
        assignPointer(U,i,j,k,n);
    }
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
        break;
    case V_TENSOR:          return *(V[m.row*numberOfComponents+m.col]);
        break;
    case VSTAR:             return *(VStar[m.row*numberOfComponents+m.col]);
        break;
    case EPSILON:           return *(materials[m.mat].epsilon);
        break;
    case ALPHARHOEPSILON:   return *(materials[m.mat].alphaRhoEpsilon);
        break;
    default: Print() << "Incorrect cell variable" << std::endl;
        exit(1);
    }
}

void Cell::assignPointer(BoxAccessCellArray& U, int i, int j, int k, MaterialSpecifier m)
{
    switch(m.var)
    {
    case ALPHA:           materials[m.mat].alpha                    = &U(i,j,k,m);
        break;
    case ALPHARHO:        (materials[m.mat].alphaRho)               = &U(i,j,k,m);
        break;
    case RHO_K:           (materials[m.mat].rho)                    = &U(i,j,k,m);
        break;
    case RHO:             rho                                       = &U(i,j,k,m);
        break;
    case RHOU:            (rhoU[m.row])                             = &U(i,j,k,m);
        break;
    case TOTAL_E:         E                                         = &U(i,j,k,m);
        break;
    case VELOCITY:        (u[m.row])                                = &U(i,j,k,m);
        break;
    case P:               p                                         = &U(i,j,k,m);
        break;
    case SOUNDSPEED:      a                                         = &U(i,j,k,m);
        break;
    case USTAR:           uStar                                     = &U(i,j,k,m);
        break;
    case LAMBDA:          (materials[m.mat].lambda)                 = &U(i,j,k,m);
        break;
    case ALPHARHOLAMBDA:  (materials[m.mat].alphaRhoLambda)         = &U(i,j,k,m);
        break;
    case SIGMA:           (sigma[m.row*numberOfComponents+m.col])   = &U(i,j,k,m);
        break;
    case V_TENSOR:        (V[m.row*numberOfComponents+m.col])       = &U(i,j,k,m);
        break;
    case VSTAR:           (VStar[m.row*numberOfComponents+m.col])   = &U(i,j,k,m);
        break;
    case EPSILON:         (materials[m.mat].epsilon)                = &U(i,j,k,m);
        break;
    case ALPHARHOEPSILON: (materials[m.mat].alphaRhoEpsilon)        = &U(i,j,k,m);
        break;
    default: Print() << "Incorrect cell variable" << std::endl;
        exit(1);
    }
}

void Cell::operator= (Cell& U)
{
    for(auto n : accessPattern.cellVariables)
    {
        (*this)(n) = U(n);
    }
}

bool Cell::check(MaterialSpecifier& n)
{
    return std::isnan((*this)(n));
}

bool Cell::contains_nan()
{

    for(auto n : accessPattern.cellVariables)
    {
        if(check(n))
        {
            Print() << accessPattern.variableNames[accessPattern[n.var]] << std::endl;
            return true;
        }
    }

    return false;

}
