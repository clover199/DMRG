/* This is a calling function of the DGEMM in the MKL lib. 
  It performs one of the matrix-matrix operations

    C = alpha * op(A) * op(B) + beta * C,

  where op(X) is one of

    op(X) = X   or  op(X) = X**T,

  alpha and beta are scalars, and A, B and C are matrices, 
  with op(A) an M-by-K matrix, op(B) a K-by-N matrix and 
  C an M-by-N matrix.
*/

#include <iostream>
#include "mkl.h"

using namespace std;

void zgemm(char transa, char transb, int m, int n, int k,
           double *a, double *b, double *c, double alpha, double beta)
{
  CBLAS_TRANSPOSE TRANSA = CblasTrans;
  int LDA = m;
  if(transa == 'N' or transa == 'n')
  {
    TRANSA = CblasNoTrans;
    LDA = k;
  }
  CBLAS_TRANSPOSE TRANSB = CblasTrans;
  int LDB = k;
  if(transb == 'N' or transb == 'n')
  {
    TRANSB = CblasNoTrans;
    LDB = n;
  }
  cblas_dgemm(CblasRowMajor, TRANSA, TRANSB, 
              m, n, k, alpha, a, LDA, b, LDB, beta, c, n);
}
