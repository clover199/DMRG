/* This is a calling function of the ZGESVD in the LAPACK lib. 
  It computes the singular value decomposition (SVD) of a 
  M-by-N matrix A, optionally computing the left and/or right 
  singular vectors.

  A = U Sigma V

  where Sigma is an M-by-N matrix which is zero except for its 
  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, 
  and V is an N-by-N orthogonal matrix. The diagonal elements of 
  SIGMA are the singular values of A; they are real and non-negative, 
  and are returned in descending order. The first min(m,n) columns 
  of U and rows of V are the left and right singular vectors of A.

  Note1: I changed the VT in the LAPACK Documentation into V.
  Note2: Fortran stores matrices columns after columns, which is 
  different with C or C++. It's like transforming all the matrices, 
  so the actual formula I use is

  a**T = v**T Sigma u**T
   n-m    n-n  n-m   m-m
   M-N    M-M  M-N   N-N
*/

#include <complex>
#include <iostream>

using namespace std;

extern "C" void zgesvd_(char *JOBU, char *JOBVT, int *M, int *N,
                        complex<double> *A, int *LDA, double *s,
                        complex<double> *U, int *LDU, complex<double> *VT, int *LDVT, 
                        complex<double> *WORK, int *LWORK, double *RWORK, int *INFO);


void tgesvd(complex<double> *a, int m, int n, double *S,
            complex<double> *u, complex<double> *v)
{
  char JOBU = 'A'; /* 'A':  all M columns of U are returned in array U;
                      'S':  the first min(m,n) columns of U (the left singular vectors)
                            are returned in the array U;
                      'O':  the first min(m,n) columns of U (the left singular vectors)
                            are overwritten on the array A;
                      'N':  no columns of U (no left singular vectors) are computed. */
  char JOBV = 'A'; /* 'A':  all N rows of V are returned in the array V;
                      'S':  the first min(m,n) rows of V (the right singular vectors)
                            are returned in the array VT;
                      'O':  the first min(m,n) rows of V (the right singular vectors)
                            are overwritten on the array A;
                      'N':  no rows of V (no right singular vectors) are computed. */
  int M = n;
  int N = m;
  complex<double> *A = new complex<double> [M*N];
  for(int i=0;i<M*N;i++) A[i] = a[i];
  int LDA = M;
  complex<double> *U = v; /* If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
                             if JOBU = 'S', U contains the first min(M,N) columns of U
                                            (the left singular vectors, stored columnwise);
                             if JOBU = 'N' or 'O', U is not referenced. */
  int LDU = M; /* if JOBU = 'S' or 'A', LDU >= M. */
  complex<double> *V = u; /* If JOBVT = 'A', V contains the N-by-N orthogonal matrix V;
                             if JOBVT = 'S', V contains the first min(M,N) rows of V
                                             (the right singular vectors, stored rowwise);
                             if JOBVT = 'N' or 'O', VT is not referenced. */
  int LDV = N; /* if JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N). */
  complex<double> work;
  int LWORK = -1; /* The dimension of the array WORK.
                     LWORK >= MAX(1,2*MIN(M,N)+MAX(M,N)).
                     For good performance, LWORK should generally be larger.

                     If LWORK = -1, then a workspace query is assumed; the routine only
                                    calculates the optimal size of the WORK array, returns
                                    this value as the first entry of the WORK array, and no
                                    error message related to LWORK is issued by XERBLA. */
  int min = M<N ? M : N;
  double *RWORK = new double [5*min]; /* dimension (5*min(M,N))
                                         On exit, if INFO > 0, RWORK(1:MIN(M,N)-1) contains
                                         the unconverged superdiagonal elements of an upper
                                         bidiagonal matrix B whose diagonal is in S (not
                                         necessarily sorted). B satisfies A = U * B * VT, so
                                         it has the same singular values as A, and singular
                                         vectors related by U and VT. */
  int INFO; /* = 0:  successful exit.
               < 0:  if INFO = -i, the i-th argument had an illegal value.
               > 0:  if DBDSQR did not converge, INFO specifies how many superdiagonals of
                     an intermediate bidiagonal form B did not converge to zero. See the
                     description of WORK above for details.*/

  /* workspace query */
  zgesvd_(&JOBU, &JOBV, &M, &N, A, &LDA, S, U, &LDU, V, &LDV, &work, &LWORK, RWORK, &INFO);
  if(INFO) cout << "Error in zgesvd: first INFO = " << INFO << endl;
  LWORK = work.real();

  int max = min+M+N > 5*min ? min+M+N : 5*min;
  if(max>LWORK) cout << "Warning in zgesvd: small LWORK. " << max << " > " << LWORK << endl;
	
  complex<double> *WORK = new complex<double> [LWORK];
  zgesvd_(&JOBU, &JOBV, &M, &N, A, &LDA, S, U, &LDU, V, &LDV, WORK, &LWORK, RWORK, &INFO);
  if(INFO) cout << "Error in zgesvd: second INFO = " << INFO << endl;

  delete A;
  delete RWORK;
  delete WORK;
}
