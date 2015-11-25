/* This is a calling function of the DSYEV in the LAPACK lib. 
  It computes all eigenvalues and, optionally, eigenvectors of a
  real symmetric matrix A.

  I defined two functions

  void dsyev(double *in, int n, double *val)
  void dsyev(double *in, int n, double *val, double *vec)

  vec is only the ground state vector.
  Note: Fortran stores matrices columns after columns, which is 
  different with C or C++. It's like transforming all the matrices.
*/

#include <iostream>
#include <cmath>
using namespace std;

extern "C" void dsyev_(char *JOBZ, char *UPLO, int *n, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO);

void dsyev(double *a, int N, double *W)
{
  char JOBZ = 'N'; /* Compute eigenvalues only */
  char UPLO = 'U'; /* Upper triangle of A is stored. */
  double *A = new double[N*N];

  /* This part transforms the input matrix into the required form and tests whether
     the input matrix is symmetric. */
  for(int i=0;i<N*N;i++) A[i] = 0;
  double max = 0;
  for(int i=0;i<N;i++) for(int j=i;j<N;j++) 
  {
    A[j*N+i] = a[i*N+j];
    double temp = abs(a[i*N+j]-a[j*N+i]);
    if(max<temp) max = temp;
  }
  if(max>1e-10) cerr << "Error in dsyev: not symmetric. The maximum difference is " << max << endl;
  
  int LDA = N;
  double work;
  int LWORK = -1; /* LWORK >= max(1,3*N-1).
                     For optimal efficiency, LWORK >= (NB+2)*N, where NB is the
                     blocksize for DSYTRD returned by ILAENV.

                     If LWORK = -1, then a workspace query is assumed; the routine only
                     calculates the optimal size of the WORK array, returns this value
                     as the first entry of the WORK array, and no error message related
                     to LWORK is issued by XERBLA. */
  int INFO; /* = 0:  successful exit
               < 0:  if INFO = -i, the i-th argument had an illegal value
               > 0:  if INFO = i, the algorithm failed to converge; i off-diagonal elements
                     of an intermediate tridiagonal form did not converge to zero. */

  dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, &work, &LWORK, &INFO); /* workspace query */
  if(INFO) cerr << "Error in dsyev: first INFO = " << INFO << endl;

  LWORK = work;
  double *WORK = new double [LWORK];
  dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);
  if(INFO) cerr << "Error in dsyev: second INFO = " << INFO << endl;

  delete A;
  delete WORK;
}

void dsyev(double *a, int N, double *W, double *vec)
{
  char JOBZ = 'V'; /* Compute eigenvalues and eigenvectors. */
  char UPLO = 'U'; /* Upper triangle of A is stored. */
  double *A = new double[N*N];  /* INFO = 0, A returns the orthonormal eigenvectors
                                   of the matrix A. */
  
  /* This part transforms the input matrix into the required form and tests
     whether the input matrix is symmetric. */
  for(int i=0;i<N*N;i++) A[i] = 0;
  double max = 0;
  for(int i=0;i<N;i++) for(int j=i;j<N;j++) 
  {
    A[j*N+i] = a[i*N+j];
    double temp = abs(a[i*N+j]-a[j*N+i]);
    if(max<temp) max = temp;
  }
  if(max>1e-10) cerr << "Error in dsyev: not symmetrix. The maximum difference is " << max << endl;
  
  int LDA = N;
  double work;
  int LWORK = -1; /* LWORK >= max(1,3*N-1).
                     For optimal efficiency, LWORK >= (NB+2)*N, where NB is the
                     blocksize for DSYTRD returned by ILAENV.

                     If LWORK = -1, then a workspace query is assumed; the routine only
                     calculates the optimal size of the WORK array, returns this value
                     as the first entry of the WORK array, and no error message related
                     to LWORK is issued by XERBLA. */
  int INFO; /* = 0:  successful exit
               < 0:  if INFO = -i, the i-th argument had an illegal value
               > 0:  if INFO = i, the algorithm failed to converge; i off-diagonal elements
                     of an intermediate tridiagonal form did not converge to zero. */

  dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, &work, &LWORK, &INFO); /* workspace query */
  if(INFO) cerr << "Error in dsyev: first INFO = " << INFO << endl;

  LWORK = work;
  double *WORK = new double [LWORK];
  dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);
  if(INFO) cerr << "Error in dsyev: second INFO = " << INFO << endl;
  
  /* Transform the eigenstate with lowest eigenvalue into vec. */
  for(int i=0;i<N;i++) vec[i] = A[i];
  delete A;
  delete WORK;
}

