/* This is a calling function of the ZHEEV in the LAPACK lib. 
  It computes all eigenvalues and, optionally, eigenvectors of a
  hermitian matrix A.

  I defined two functions

  void dsyev(complex<double> *in, int n, double *val)
  void dsyev(complex<double> *in, int n, double *val, complex<double> *vec)

  vec is only the ground state vector.
  Note: Fortran stores matrices columns after columns, which is 
  different with C or C++. It's like transforming all the matrices.
*/

#include <iostream>
#include <complex>
#include <cmath>

using namespace std;

extern "C" void zheev_(char *JOBZ, char *UPLO, int *N, complex<double> *A, int *LDA, 
double *W, complex<double> *WORK, int *LWORK, double *RWORK, int *INFO);

void zheev(complex<double> *a, int N, double *W)
{
  char JOBZ='N'; /* Compute eigenvalues only */
  char UPLO='U'; /* Upper triangle of A is stored. */
  complex<double> *A = new complex<double> [N*N];

  /* This part transforms the input matrix into the required form and tests whether
     the input matrix is hermitian. */
  for(int i=0;i<N*N;i++) A[i] = 0;
  double max = 0;
  for(int i=0;i<N;i++) for(int j=i;j<N;j++) 
  {
    A[j*N+i] = a[i*N+j];
    double temp = abs(a[i*N+j]-conj(a[j*N+i]));
    if(max<temp) max = temp;
  }
  if(max>1e-10) cout << "Error in zheev: not hermitian. The maximum difference is "
                     << max << endl;

  int LDA = N;
  complex<double> work;
  int LWORK=-1; /* LWORK >= max(1,2*N-1).
                   For optimal efficiency, LWORK >= (NB+1)*N, where NB is the blocksize
                   for ZHETRD returned by ILAENV.

                   If LWORK = -1, then a workspace query is assumed; the routine only
                   calculates the optimal size of the WORK array, returns this value as
                   the first entry of the WORK array, and no error message related to LWORK
                   is issued by XERBLA. */
  double *RWORK = new double [3*N-2];
  int INFO; /* = 0:  successful exit
               < 0:  if INFO = -i, the i-th argument had an illegal value
               > 0:  if INFO = i, the algorithm failed to converge; i off-diagonal elements
                     of an intermediate tridiagonal form did not converge to zero. */

  zheev_(&JOBZ, &UPLO, &N, A, &LDA, W, &work, &LWORK, RWORK, &INFO); /* workspace query */
  if(INFO) cout << "Error in zheev: first INFO = " << INFO << endl;

  LWORK=work.real();
  complex<double> *WORK = new complex<double> [LWORK];
  zheev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, RWORK, &INFO);
  if(INFO) cout << "Error in zheev: second INFO = " << INFO << endl;

  delete A;
  delete RWORK;
  delete WORK;
}

void zheev(complex<double> *a, int N, double *W, complex<double> *vec)
{
  char JOBZ='V'; /* Compute eigenvalues and eigenvectors. */
  char UPLO='U'; /* Upper triangle of A is stored. */
  complex<double> *A = new complex<double> [N*N];

  /* This part transforms the input matrix into the required form and tests whether
     the input matrix is hermitian. */
  for(int i=0;i<N*N;i++) A[i] = 0;
  double max = 0;
  for(int i=0;i<N;i++) for(int j=i;j<N;j++) 
  {
    A[j*N+i] = a[i*N+j];
    double temp = abs(a[i*N+j]-conj(a[j*N+i]));
    if(max<temp) max = temp;
  }
  if(max>1e-10) cout << "Error in zheev: not hermitian. The maximum difference is " << max << endl;

  int LDA = N;
  complex<double> work;
  int LWORK=-1; /* LWORK >= max(1,2*N-1).
                   For optimal efficiency, LWORK >= (NB+1)*N, 
                   where NB is the blocksize for ZHETRD returned by ILAENV.

                   If LWORK = -1, then a workspace query is assumed; the routine only
                   calculates the optimal size of the WORK array, returns this value as
                   the first entry of the WORK array, and no error message related to
                   LWORK is issued by XERBLA. */
  double *RWORK = new double [3*N-2];
  int INFO; /* = 0:  successful exit
               < 0:  if INFO = -i, the i-th argument had an illegal value
               > 0:  if INFO = i, the algorithm failed to converge; i off-diagonal elements
                     of an intermediate tridiagonal form did not converge to zero. */

  zheev_(&JOBZ, &UPLO, &N, A, &LDA, W, &work, &LWORK, RWORK, &INFO); /* workspace query */
  if(INFO) cout << "Error in zheev: first INFO = " << INFO << endl;

  LWORK=work.real();
  complex<double> *WORK = new complex<double> [LWORK];
  zheev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, RWORK, &INFO);
  if(INFO) cout << "Error in zheev: second INFO = " << INFO << endl;

  /* Transform the eigenstate with lowest eigenvalue into vec. */
  for(int i=0;i<N;i++) vec[i] = A[i];

  delete A;
  delete RWORK;
  delete WORK;

}
