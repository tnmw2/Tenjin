#include "cell.h"

Cell::Cell(BoxAccessCellArray& U, int i, int j, int k, Material_type _phase) : accessPattern(U.accessPattern), parent_i(i), parent_j(j), parent_k(k), parent(&U)
{
    materials.resize(U.numberOfMaterials);

    for(int n = 0; n< U.numberOfMaterials; n++)
    {
        materials[n].allocateSpace(this,n);

        materials[n].phase = _phase;
    }

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
    case RHO:               return *(materials[m.mat].rho);
        break;
    case RHOU:              return *(materials[m.mat].rhoU[m.row]);
        break;
    case TOTAL_E:           return *(materials[m.mat].E);
        break;
    case VELOCITY:          return *(materials[m.mat].u[m.row]);
        break;
    case P:                 return *(materials[m.mat].p);
        break;
    case SOUNDSPEED:        return *(materials[m.mat].a);
        break;
    case USTAR:             return *(materials[m.mat].uStar);
        break;
    case SIGMA:             return *(materials[m.mat].sigma[m.row*numberOfComponents+m.col]);
        break;
    case V_TENSOR:          return *(materials[m.mat].V[m.row*numberOfComponents+m.col]);
        break;
    case VSTAR:             return *(materials[m.mat].VStar[m.row*numberOfComponents+m.col]);
        break;
    case EPSILON:           return *(materials[m.mat].epsilon);
        break;
    case RHOEPSILON:        return *(materials[m.mat].rhoEpsilon);
        break;
    default: Print() << "Incorrect cell variable" << std::endl;
        exit(1);
    }
}

void Cell::assignPointer(BoxAccessCellArray& U, int i, int j, int k, MaterialSpecifier m)
{
    switch(m.var)
    {
    case RHO:             (materials[m.mat].rho)                                     = &U(i,j,k,m);
        break;
    case RHOU:            (materials[m.mat].rhoU[m.row])                             = &U(i,j,k,m);
        break;
    case TOTAL_E:         (materials[m.mat].E)                                       = &U(i,j,k,m);
        break;
    case VELOCITY:        (materials[m.mat].u[m.row])                                = &U(i,j,k,m);
        break;
    case P:               (materials[m.mat].p)                                       = &U(i,j,k,m);
        break;
    case SOUNDSPEED:      (materials[m.mat].a)                                       = &U(i,j,k,m);
        break;
    case USTAR:           (materials[m.mat].uStar)                                   = &U(i,j,k,m);
        break;
    case SIGMA:           (materials[m.mat].sigma[m.row*numberOfComponents+m.col])   = &U(i,j,k,m);
        break;
    case V_TENSOR:        (materials[m.mat].V[m.row*numberOfComponents+m.col])       = &U(i,j,k,m);
        break;
    case VSTAR:           (materials[m.mat].VStar[m.row*numberOfComponents+m.col])   = &U(i,j,k,m);
        break;
    case EPSILON:         (materials[m.mat].epsilon)                                 = &U(i,j,k,m);
        break;
    case RHOEPSILON:      (materials[m.mat].rhoEpsilon)                              = &U(i,j,k,m);
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

void Material::allocateSpace(Cell* ptr, int n)
{
    parent          = ptr;
    materialNumber  = n;

    rhoU.resize(3);
    u.resize(3);
    sigma.resize(9);
    V.resize(9);
    VStar.resize(9);

}

void Material::operator= (Material& M)
{
    for(auto n : parent->accessPattern.cellVariables)
    {
        if(n.mat == M.materialNumber)
        {
            parent->operator()(n) = M.parent->operator()(n);
        }
    }
}
