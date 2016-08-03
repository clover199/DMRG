/* This is a calling function of the DGEMM in the BLAS lib. 
  It performs one of the matrix-matrix operations

    C = alpha * op(A) * op(B) + beta * C,

  where op(X) is one of

    op(X) = X   or  op(X) = X**T,

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

using namespace std;

extern "C" void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
                       double *alpha, double *A, int *LDA, double *B, int *LDB,
                       double *beta, double *C, int *LDC);

void dgemm(char transa, char transb, int m, int n, int K,
           double *a, double *b, double *C, double alpha, double beta)
{
  if(m>0 and n>0 and K>0)
  {
    char TRANSA = transb;
    char TRANSB = transa;
    int M = n;
    int N = m;
    int LDA = K;
    int LDB = N;
    double *A = b;
    if(TRANSA == 'N' or TRANSA == 'n') LDA = M;
    double *B = a;
    if(TRANSB == 'N' or TRANSB == 'n') LDB = K;
    int LDC = M;
    dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
  }
}
