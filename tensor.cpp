#include "cblas.h"
#include "cblas_f77.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <algorithm>

#include "tensor.h"

extern "C"
{
    //LAPACK

    //Eigenvalues
    extern int  dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
    //SVD
    //extern void dgesdd_(char, int, int, double*, int, double*, double*, int, double*, int, double*, int, int*, int);
    extern void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
    extern void dgesdd_( char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*, int*);
}

//BLAS matrix Multiply
void cblas_dgemm(CBLAS_LAYOUT, CBLAS_TRANSPOSE, CBLAS_TRANSPOSE, int, int, int, double, const double*, int, const double*, int, double, double*, int);
void cblas_dcopy(const int N,const double* X, const int incX, double* Y,const int incY);

void matrixPrint(double* U,int N =3)
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            std::cout<<U[i*N+j] << "    ";
        }

        std::cout<<std::endl;
    }

    std::cout<<std::endl;

    return;
}

/** Calculates the singular value decomposition of the matrix \f$ A = LSR^T \f$.
 * All matrices are assumed to be NxN square matrices. We call the BLAS functions
 * dcopy and dgesdd. Note that \f$S\f$ is a vector not a matrix.
 * \param [in] A Input matrix
 * \param [out] S Singular value vector
 * \param [out] L Orthogonal matrix
 * \param [out] RT Orthogonal matrix
 */
void singularValueDecomposition(double* S, double* A, double* L, double* RT, int N)
{
    int INFO=0;
    int LWORK =4*N*N+6*N+N;
    //int LWORK =5*N;

    int IWORK[8*N];
    double WORK[LWORK];
    char all = 'A';
    int one=1;

    double Acopy[N*N];

    cblas_dcopy(N*N,A,one,Acopy,one);

    dgesdd_(&all,&N,&N,Acopy,&N,S,L,&N,RT,&N,WORK,&LWORK,IWORK,&INFO);

    //dgesvd_(&all,&all,&N,&N,Acopy,&N,S,L,&N,RT,&N,WORK,&LWORK,&INFO);

    return;
}


/**Multiplies two square matrices together.
 * All matrices are assumed to be NxN square matrices. We call the BLAS function
 * dgemm.
 * \param [in] M1 First matrix
 * \param [in] M2 Second matrix
 * \param [out] product Product matrix
 */
double* squareMatrixMultiply(double* M1,double*M2,double*product,int N)
{
    cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,N,N,N,1.0,M1,N,M2,N,0.0,product,N);

    return product;
}

/**Multiplies two square matrices together, transposing the first.
 * All matrices are assumed to be NxN square matrices. We call the BLAS function
 * dgemm.
 * \param [in] M1 First matrix
 * \param [in] M2 Second matrix
 * \param [out] product Product matrix
 */
double* squareMatrixTransposeMultiply(double* M1,double*M2,double*product,int N)
{
    cblas_dgemm(CblasRowMajor,CblasTrans, CblasNoTrans,N,N,N,1.0,M1,N,M2,N,0.0,product,N);

    return product;
}


/**Multiplies two square matrices together, transposing the second.
 * All matrices are assumed to be NxN square matrices. We call the BLAS function
 * dgemm.
 * \param [in] M1 First matrix
 * \param [in] M2 Second matrix
 * \param [out] product Product matrix
 */
double* squareMatrixMultiplyTranspose(double* M1,double*M2,double*product,int N)
{
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,N,N,N,1.0,M1,N,M2,N,0.0,product,N);

    return product;
}

/**Transposes an NxN matrix.
 * \param [in] M Input matrix
 * \param [out] MT Transpose matrix
 */
double* transpose(double* M, double*MT, int N)
{
    int i,j;
    //#pragma omp parallel for schedule(static) private(i,j) shared(MT,M,N)
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            MT[i*N+j]=M[j*N+i];
        }
    }

    return MT;
}

/**Finds the trace of a square matrix.
 * \param [in] M Input matrix
 * \param [out] tr The trace = \f$ M_{ii} \f$
 */
double trace(double*M, int N)
{
    double tr =0.0;

    for(int i=0;i<N;i++)
    {
        tr+=M[i*N+i];
    }

    return tr;
}

/**Finds the second matrix invariant of a square matrix.
 * \param [in] M Input matrix
 */
double I2(double*M, int N)
{
    double Msquared [N*N];
    double tempM[N*N];

    matrixCopy(M,tempM);

    double Tr = trace(M,N);

    return 0.5*((Tr*Tr)-trace(squareMatrixMultiply(M,tempM,Msquared,N),N));
}

/**Finds \f$ x^y \f$.
 * \param [in] x Base
 * \param [in] y Exponent
 */
double intpow(double x, int y)
{
    double z = 1.0;

    for(;y>0;y--)
    {
        z*=x;
    }

    return z;
}

/**Finds the inverse of a square matrix.
 * \param [in] M1 Input matrix
 * \param [out] Minv The inverse matrix
 */
double* invert(double* M1,double*Minv,int N)
{
    double cofactor[N*N];
    double cofactorT[N*N];
    double submatrix[(N-1)*(N-1)];
    int counter;

    double D = det(M1,N);

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            counter=0;

            for(int k=0;k<N;k++)
            {
                for(int m=0;m<N;m++)
                {
                    if(k==i || m==j)
                    {
                        continue;
                    }
                    else
                    {
                        submatrix[counter]=M1[k*N+m];
                        counter++;
                    }
                }
            }
            //cofactor[i*N+j]=std::pow(-1.0,i+j)*det(submatrix,N-1);
            cofactor[i*N+j]  =intpow(-1.0,i+j)*det(submatrix,N-1);
        }
    }

    transpose(cofactor,cofactorT,N);

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            Minv[i*N+j]=cofactorT[i*N+j]/D;
        }
    }

    return Minv;
}

/** Deprecated function to find the Eigenvalues of a 3x3 matrix.
 * \param [in] M input matrix.
 * \param [out] lambda vector of eigenvalues.
 */
void find3x3EigenValues(double* M, double* lambda)
{
    int n=3;
    char Nchar='N';
    double vl[n];
    double vr[n];
    int one=1;
    int lwork=6*n;
    double* eigImag=new double[n];
    double* work   =new double[lwork];

    double Mcopy[n*n];

    for(int i=0;i<n;i++)
    {
        vl[i]=0.0;
        vr[i]=0.0;
    }

    cblas_dcopy(n*n,M,one,Mcopy,one);


    int info;

    //LAPACK Eigenvalue subroutine

    dgeev_(&Nchar,&Nchar,&n,Mcopy,&n,lambda,eigImag,vl,&one,vr,&one, work,&lwork,&info);

    if (info!=0)
    {
        std::cout << "Error: dgeev returned error code " << info << std::endl;
        return;
    }

    for(int i=0;i<n;i++)
    {
        if(eigImag[i]!=0.0)
        {
            std::cout << eigImag[i] << std::endl;
            throw "Imaginary eigenvalues";
        }
    }


    delete [] eigImag;
    delete [] work;

    return;

}


/** Calculates the determinant of a matrix.
 *  Recursively calls this function until a 2x2 matrix is reached.
 * \param [in] M input matrix.
 * \param [out] determinant The determinant.
 */
double det(double*M, int N)
{
    double determinant=0.0;

    if(N==2)
    {
        return M[0]*M[N+1]-M[1]*M[N];
    }
    else
    {
        for(int x=0;x<N;x++)
        {
            double cofactor[(N-1)*(N-1)];

            int counter=0;

            for(int i=0;i<N;i++)
            {
                for(int j=0;j<N;j++)
                {
                    if(i==0 || j==x)
                    {
                        continue;
                    }
                    else
                    {
                        cofactor[counter]=M[i*N+j];
                        counter++;
                    }
                }
            }

            determinant+=intpow(-1.0,x%2)*M[x]*det(cofactor,N-1);
        }

        return determinant;
    }
}

/** The Levi-Civita tensor.
 */
double epsilon(int i, int j, int k)
{
    if(i==j || i==k || j==k )
    {
        return 0.0;
    }
    else if((i==1 && j==2 && k==3) || (i==3 && j==1 && k==2) || (i==2 && j==3 && k==1) )
    {
        return 1.0;
    }
    else
    {
        return -1.0;
    }
}

/** Finds the spherical part of a tensor.
 * \param [in] M input matrix
 */
double sphericalpart(double* M, int N)
{
    return trace(M,N)/((double)N);
}

/** Finds the deviatoric part of a tensor.
 * \param [in] M input matrix
 * \param [out] dev deviatoric matrix
 */
double* deviator(double* dev, double* M, int N)
{
    double p = sphericalpart(M,N);

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            dev[i*N+j] = M[i*N+j]-p*delta<double>(i,j);
        }
    }

    return dev;
}

/** Deprecated function to find the Finger tensor from the inverse deformation tensor.
 */
double* fingerTensor(double* g, double* G, int N)
{
    squareMatrixTransposeMultiply(g,g,G,N);
    return G;
}

/** Find the stretch tensor using a polar decomposition.
 * Goes via a singular value decomposition.
 * \param [in] F input tensor
 * \param [out] V Stretch tensor
 */
void getStretchTensor(double* V, double* F, int N)
{
    double Svec [N];
    double S	[N*N];
    double L 	[N*N];
    double RT 	[N*N];
    double temp [N*N];

    singularValueDecomposition(Svec,F,RT,L,N); // !!! Single value decomposition is wrong! had to swap L and RT!


    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            S[i*N+j] = Svec[i]*delta<double>(i,j);
        }
    }

    squareMatrixMultiply(L,S,temp,N);
    squareMatrixMultiplyTranspose(temp,L,V,N);

    return;
}

/** Find the logarithm of a matrix using a singular value decomposition.
 * \param [in] A input matrix
 * \param [out] Log Matrix logarithm
 */
void matrixLogarithm(double* Log, double* A, int N)
{
    double Svec [N];
    double S	[N*N];
    double L 	[N*N];
    double RT 	[N*N];
    double temp [N*N];

    //#pragma omp critical
    singularValueDecomposition(Svec,A,RT,L,N); // !!! Single value decomposition is wrong! had to swap L and RT!


    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            S[i*N+j] = std::log(Svec[i])*delta<double>(i,j);
        }
    }

    squareMatrixMultiply(L,S,temp,N);
    squareMatrixMultiply(temp,RT,Log,N);

    return;
}

/** Find the square root of a matrix using a singular value decomposition.
 * \param [in] A input matrix
 * \param [out] Sqrt Matrix square root
 */
void matrixSquareRoot(double* Sqrt, double* A, int N)
{
    double Svec [N];
    double S	[N*N];
    double L 	[N*N];
    double RT 	[N*N];
    double temp [N*N];

    singularValueDecomposition(Svec,A,RT,L,N); // !!! Single value decomposition is wrong! had to swap L and RT!


    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            S[i*N+j] = std::sqrt(Svec[i])*delta<double>(i,j);
        }
    }

    squareMatrixMultiply(L,S,temp,N);
    squareMatrixMultiply(temp,RT,Sqrt,N);

    return;
}

/** Find the exponential of a matrix using Taylor expansion.
 * \param [in] A input matrix
 * \param [in] order Expansion order
 * \param [out] exp Matrix exponential
 */
void matrixExponential(double* exp, double* A, int order, int N)
{
    double temp1[N*N];
    double temp2[N*N];

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            exp[i*N+j] = delta<double>(i,j);
            temp1[i*N+j] = delta<double>(i,j);
        }
    }

    int factorial = 1;

    for(int n=1;n<order;n++)
    {
        factorial *= n;

        squareMatrixMultiply(temp1,A,temp2,N);

        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                exp[i*N+j] += temp2[i*N+j]/factorial;
                temp1[i*N+j] = temp2[i*N+j];
            }
        }

    }

    return;
}

/** Produce a copy of a matrix using BLAS.
 * \param [in] M input matrix
 * \param [out] Mcopy Copy
 */
void matrixCopy(double* M, double* Mcopy, int N)
{
    int one = 1;

    cblas_dcopy(N*N,M,one,Mcopy,one);

    return;
}

/** Find the Frobenius norm of a matrix
 * \param [in] M input matrix
 * \param [out] Mcopy Copy
 */
double Jnorm(double* M, int N)
{
    double temp[N*N];

    matrixCopy(M,temp,N);

    double product[N*N];

    squareMatrixMultiplyTranspose(M,temp,product,N);

    return sqrt(trace(product));
}

/** Heaviside function.
 */
double Heaviside(double x)
{
    return (double)(x>=0.0);
}
