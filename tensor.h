#ifndef TENSOR_H
#define TENSOR_H

#include "cblas.h"
//#include "cblas_f77.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <algorithm>

double intpow(double, int);

double* squareMatrixMultiply(double*,double*,double*,int=3);
double* squareMatrixTransposeMultiply(double*,double*,double*,int=3);
double* squareMatrixMultiplyTranspose(double*,double*,double*,int=3);
double* transpose(double*,double*,int=3);
double trace(double*, int=3);
double I2(double*, int=3 );
double det(double*, int=3 );
double* invert(double* ,double*,int=3 );
void find3x3EigenValues(double*, double*);
void singularValueDecomposition(double*, double *, double*, double*, int=3);
double epsilon(int,int,int);
double sphericalpart(double*, int=3);
double* deviator(double*, double*, int=3);
double* fingerTensor(double*, double*,int=3);
void getStretchTensor(double*,double*,int =3);
void matrixLogarithm(double*, double*, int=3);
void matrixExponential(double*, double*, int, int=3);
void matrixSquareRoot(double*, double*, int=3);
void matrixCopy(double*, double*, int N=3);
double Jnorm(double*, int =3);
double Heaviside(double);


/** Dirac delta function.
 */
template<typename T> T delta(int i, int j)
{
    return (T)(i==j);
}

/** Sign function.
 */
template <typename I,typename O> O sgn(I val)
{
    return O((I(0) < val) - (val < I(0)));
}

/** Heaviside function.
 */
template <typename In,typename Out> Out Heaviside(In x)
{
    return (Out)(x > In(0));
}


#endif // TENSOR_H

