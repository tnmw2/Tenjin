#include "cellarray.h"

BoxAccessCellArray::BoxAccessCellArray(const Box& bx, FArrayBox& fb, Array4<Real> const& prop_arr, CellArray& U) : box{bx}, fab{fb}, arr{prop_arr}, accessPattern{U.accessPattern}{}

Real& BoxAccessCellArray::operator()(int i, int j, int k, Variable var, int mat, int row, int col)
{
    return (*this)(MaterialSpecifier(var,mat,row,col),i,j,k);
}

Real& BoxAccessCellArray::operator()(int i, int j, int k, int var, int mat, int row, int col)
{
    return (*this)(MaterialSpecifier((Variable)var,mat,row,col),i,j,k);
}

Real& BoxAccessCellArray::operator()(Variable var, int mat, int row, int col, int i, int j, int k)
{
    return (*this)(MaterialSpecifier(var,mat,row,col),i,j,k);
}

Real& BoxAccessCellArray::operator()(MaterialSpecifier& m, int i, int j, int k)
{
    switch(m.var)
    {
    case RHO:       return arr(i,j,k,accessPattern[m.var]);
        break;
    case RHOU:      return arr(i,j,k,accessPattern[m.var]+m.row);
        break;
    case TOTAL_E:   return arr(i,j,k,accessPattern[m.var]);
        break;
    case VELOCITY:  return arr(i,j,k,accessPattern[m.var]+m.row);
        break;
    case P:         return arr(i,j,k,accessPattern[m.var]);
        break;
    case SOUNDSPEED:return arr(i,j,k,accessPattern[m.var]);
        break;
    default: Print() << "Incorrect Access variable" << std::endl;
        exit(1);
    }
}

Real& BoxAccessCellArray::operator()(MaterialSpecifier m, int i, int j, int k)
{
    switch(m.var)
    {
    case RHO:       return arr(i,j,k,accessPattern[m.var]);
        break;
    case RHOU:      return arr(i,j,k,accessPattern[m.var]+m.row);
        break;
    case TOTAL_E:   return arr(i,j,k,accessPattern[m.var]);
        break;
    case VELOCITY:  return arr(i,j,k,accessPattern[m.var]+m.row);
        break;
    case P:         return arr(i,j,k,accessPattern[m.var]);
        break;
    case SOUNDSPEED:return arr(i,j,k,accessPattern[m.var]);
        break;
    default: Print() << "Incorrect Access variable" << std::endl;
        exit(1);
    }
}
