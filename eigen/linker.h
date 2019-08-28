using namespace std;
#include <complex>

#ifndef LINKER_H
#define LINKER_H


extern "C"
{
    //EIGENSOLVERMXM (M, A, rows, cols, nev, ncv, w, v, neigh, info)
    //This function solves a generalized eigenproblem using LAPACK functions for
    //solving systems and a custom (compacted) matrix-vector product
    void eigensolvermxm_(complex<double> *m, complex<double> *a,
                         int *rows, int *cols, int *nev, int *ncv,
                         complex<double> *w, complex<double> *v,
                         int *neigh, int *info, int* numIter,
                         int* numOpS, int* numOpBX, int* numReOR);

    //EIGENSOLVERMXMAXB (M, A, rows, cols, nev, ncv, w, v, neigh, info)
    //This method solves a generalized eigenproblem using a custom (compacted)
    //matrix-vector product and CG for solving systems
    void eigensolvermxmaxb_(complex<double> *m, complex<double> *a,
                            int *rows, int *cols, int *nev, int *ncv,
                            complex<double> *w, complex<double> *v,
                            int *neigh, int *info);

    //EIGENSOLVERQ (M, A, n, nev, ncv, w, v, info)
    //This function solves a generalized eigenproblem using LAPACK functions for
    //matrix-vector product and solving systems
    void eigensolverq_(complex<double> *m, complex<double> *a,
                       int *n, int *nev, int *ncv,
                       complex<double> *w, complex<double> *v,
                       int *info);

    //ZPOTRF( UPLO, N, A, LDA, INFO )
    //computes the Cholesky factorization of a complex Hermitian positive definite matrix A.
    void zpotrf_(char *uplo, int *n, complex<double> *A, int *lda, int* info);

    //ZPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
    //solves a system of linear equations A*X = B with a Hermitian
    //positive definite matrix A using the Cholesky factorization
    //A = U**H * U or A = L * L**H computed by ZPOTRF.
    void zpotrs_(char *uplo, int *aRows, int *bCols, complex<double> *A,
                 int *aCols, complex<double> *B, int *bRows, int* info);

    //ZGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO)
    //computes the singular value decomposition (SVD) of a complex
    //M-by-N matrix A, optionally computing the left and/or right singular
    //vectors. The SVD is written
    //
    //      A = U * SIGMA * conjugate-transpose(V)
    //
    //where SIGMA is an M-by-N matrix which is zero except for its
    //min(m,n) diagonal elements, U is an M-by-M unitary matrix, and
    //V is an N-by-N unitary matrix.  The diagonal elements of SIGMA
    //are the singular values of A; they are real and non-negative, and
    //are returned in descending order.  The first min(m,n) columns of
    //U and V are the left and right singular vectors of A.
    //Note that the routine returns V**H, not V.
    void zgesvd_(char *jobu, char *jobvt, int *aRows, int *aCols,
                 complex<double> *A, int *aLda, double *S,
                 complex<double> *U, int *uRows, complex<double> *VT,
                 int *vtRows, complex<double> *work, int *lengthWork,
                 double *rwork, int *info);

    //DOUBLE PRECISION DZNRM2( N, X, INCX )
    //DZNRM2 returns the euclidean norm of a vector via the function
    //name, so that DZNRM2 := sqrt( conjg( x' )*x )
    double dznrm2_(int *n, complex<double> *x, int *incx);

    //ZGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK, INFO )
    //This function estimates the reciprocal of the condition number of a general
    //complex matrix A, in either the 1-norm or the infinity-norm, using
    //the LU factorization computed by ZGETRF.
    //An estimate is obtained for norm(inv(A)), and the reciprocal of the
    //condition number is computed as RCOND = 1 / ( norm(A) * norm(inv(A)) ).
    void zgecon_(char *norm, int *row, complex<double> *a, int *lda, double *anorm,
                 double *rcond, complex<double> *work, double *rwork, int *info);

    //ZGETRF( M, N, A, LDA, IPIV, INFO )
    //ZGETRF computes an LU factorization of a general M-by-N matrix A
    //using partial pivoting with row interchanges.
    void zgetrf_(int *rows, int *cols, complex<double> *A, int *lda,
                 int *ipiv, int *info);

    //DOUBLE PRECISION ZLANGE( NORM, M, N, A, LDA, WORK )
    //ZLANGE  returns the value of the one norm,  or the Frobenius norm, or
    //the  infinity norm,  or the  element of  largest absolute value  of a
    //complex matrix A.
    //ZLANGE returns the value
    //ZLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
    //         (
    //         ( norm1(A),         NORM = '1', 'O' or 'o'
    //         (
    //         ( normI(A),         NORM = 'I' or 'i'
    //         (
    //         ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
    //where  norm1  denotes the  one norm of a matrix (maximum column sum),
    //normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    //normF  denotes the  Frobenius norm of a matrix (square root of sum of
    //squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
    double zlange_(char *norm, int *rows, int *cols, complex<double> *a,
                 int *lda, complex<double> *work);
}
#endif // LINKER_H
