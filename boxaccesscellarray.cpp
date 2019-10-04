#include "cellarray.h"

BoxAccessCellArray::BoxAccessCellArray(const Box& bx, FArrayBox& fb, Array4<Real> const& prop_arr) : box{bx}, fab{fb}, arr{prop_arr}{}

