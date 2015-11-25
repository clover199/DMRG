/* This is a calling function of the ZGEMM in the BLAS lib. 
  It performs one of the matrix-matrix operations

    C = alpha * op(A) * op(B) + beta * C,

  where op(X) is one of

    op(X) = X ,  op(X) = X**T  or  op(X) = X**H

  alpha and beta are scalars, and A, B and C are matrices, 
  with op(A) an M-by-K matrix, op(B) a K-by-N matrix and 
  C an M-by-N matrix.

  Note: Fortran stores matrices columns after columns, which is 
  different with C or C++. It's like transforming all the matrices, 
  so the actual formula I use is

  c**T = alpha * op(b)**T * op(a)**T + beta * c**T
   n-m             n-k         k-m             n-m
   M-N             M-K         K-N             M-N
*/

#include <complex>

using namespace std;
 
extern "C" void zgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, 
                       complex<double> *alpha, complex<double> *A, int *LDA, 
                       complex<double> *B, int *LDB, complex<double> *beta, 
                       complex<double> *C, int *LDC);
				
void tgemm(char transa, char transb, int m, int n, int K,
           complex<double> *a, complex<double> *b, complex<double> *C,
           complex<double> alpha, complex<double> beta)
{
  char TRANSA = transb;
  char TRANSB = transa;
  int M = n;
  int N = m;
  int LDA = K;
  int LDB = N;
  complex<double> *A = b;
  if(TRANSA == 'N' or TRANSA == 'n') LDA = M;
  complex<double> *B = a;
  if(TRANSB == 'N' or TRANSB == 'n') LDB = K;
  int LDC = M;
  zgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}
