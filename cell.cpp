#include "cell.h"

Cell::Cell(BoxAccessCellArray& U, int i, int j, int k, Material_type _phase) : phase(_phase), accessPattern(U.accessPattern), parent_i(i), parent_j(j), parent_k(k), parent(&U){}

Real& Cell::operator()(Variable var, int mat, int row, int col)
{
    return (*this)(MaterialSpecifier(var,mat,row,col));
}

Real& Cell::operator()(MaterialSpecifier m)
{
    return parent->operator ()(parent_i,parent_j,parent_k,m);
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
