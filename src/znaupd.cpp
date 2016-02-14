/*
  In this header file, I have defined a simplified function call to
  the ARPACK solver routine for a complex, asymmetric eigenvalue
  problem Av = lv.  The looping procedure and final extraction of
  eigenvalues and vectors is handled automatically.  Most of the
  parameters to the FORTRAN functions are hidden from the user, since
  most of them are determined from user input anyway.
  
  The remaining parameters to the function calls are as follows:
  
    znaupd(int n, int nev, complex<double> *Evals)
    znaupd(int n, int nev, complex<double> *Evals, complex<double> *Evecs)

    n: the order of the square matrix A
    nev: the number of eigenvalues to be found, starting at the
         bottom.  Note that the highest eigenvalues, or some
	 other choices, can be found.  For now, the choice of
	 the lowest nev eigenvalues is hard-coded.
    Evals: a one-dimensional array of length nev to hold the
           eigenvalues.
    Evecs: a two-dimensional array of size nev by n to hold the
           eigenvectors.  If this argument is not provided, the
	   eigenvectors are not calculated.  Note that the
	   vectors are stored as columns of Evecs, so that the
	   elements of vector i are the values Evecs[j*nev+i].

  The function is overdefined, so that you can call it in either
  fashion and the appropriate version will be run.

  To use these function calls, there must be a function
  defined in the calling program of the form

    av(int n, complex<double> *in, complex<double> *out)

  where the function av finds out = A.in, and n is the order of the
  matrix A.  This function must be defined before the statement that
  includes this header file, as it needs to know about this function.
  It is used in the looping procedure.

  Note that 0 < nev < n-1.

  Scot Shaw
  30 August 1999 

  I made some changes to the original file. 
  Evecs is changed to be a one-dimensional array of size n to hold
  the eigenvector with lowest eigenvalue.
  -- Ye Zhuang
*/

#include <iostream>
#include <complex>
#include <cmath>

using namespace std;

extern "C" void znaupd_(int *ido, char *bmat, int *n, char *which,
			int *nev, double *tol, complex<double> *resid,
			int *ncv, complex<double> *v, int *ldv,
			int *iparam, int *ipntr, complex<double> *workd,
			complex<double> *workl, int *lworkl,
			double *rwork, int *info);

extern "C" void zneupd_(int *rvec, char *All, int *select,
			complex<double> *d, complex<double> *z, int *ldz,
			double *sigma, complex<double> *workev, char *bmat,
			int *n, char *which, int *nev, double *tol,
			complex<double> *resid, int *ncv,
			complex<double> *v, int *ldv, int *iparam,
			int *ipntr, complex<double> *workd,
			complex<double> *workl, int *lworkl,
			double *rwork, int *ierr);

void znaupd(int n, int nev, double *Evals, complex<double> *Evecs, 
            void av(int n, complex<double> *in, complex<double> *out))
{
  int ido = 0; /* Initialization of the reverse communication parameter. */
  char bmat = 'I'; /* Specifies that the right hand side matrix should be the identity matrix;
                      this makes the problem a standard eigenvalue problem.
                      Setting bmat = "G" would have us solve the problem Av = lBv (this would
                      involve using some other programs from BLAS, however). */
  char which[3] = "SR"; /* LM: largest magnitude
			   SM: smallest magnitude
			   LR: largest real component
			   SR: smallest real compoent
			   LI: largest imaginary component
			   SI: smallest imaginary component */
  double tol = 1e-10; /* Sets the tolerance; tol<=0 specifies machine precision */
  complex<double> *resid = new complex<double>[n];
  int ncv = 6*nev; /* The largest number of basis vectors that will  be used in the
                      Implicitly Restarted Arnoldi Process.  Work per major iteration
                      is proportional to N*NCV*NCV. */
  if (ncv>n) ncv = n;
  int ldv = n;
  complex<double> *v = new complex<double>[ldv*ncv];
  int iparam[11]; /* An array used to pass information to the routines about their
                     functional modes. */
  iparam[0] = 1;   // Specifies the shift strategy (1->exact)
  iparam[2] = 3*n; // Maximum number of iterations
  iparam[6] = 1;   /* Sets the mode of dsaupd.
		      1 is exact shifting,
		      2 is user-supplied shifts,
		      3 is shift-invert mode,
		      4 is buckling mode,
		      5 is Cayley mode. */
  int ipntr[14]; /* Indicates the locations in the work array workd
			  where the input and output vectors in the
			  callback routine are located. */
  complex<double> *workd = new complex<double>[3*n];
  int lworkl = 5*ncv*(ncv+1); /* Length of the workl array */
  complex<double> *workl = new complex<double>[lworkl];
  double *rwork = new double[ncv];
  int info = 0; /* Passes convergence information out of the iteration routine. */
  int rvec = 1; /* =0: Specifies that eigenvectors should NOT be calculated
                   =1: Specifies that eigenvectors should be calculated */
  char howmny = 'A';
  int *select = new int[ncv];
  complex<double> *d = new complex<double>[ncv]; /* This vector will return the eigenvalues
                                                    from the second routine, dseupd. */
  double sigma;
  complex<double> *workev = new complex<double>[3*ncv];
  int ierr;

  /* Here we enter the main loop where the calculations are performed. The communication
     parameter ido tells us when the desired tolerance is reached, and at that point we
     exit and extract the solutions. */

 do {
   znaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, 
	   &ncv, v, &ldv, iparam, ipntr, workd, workl,
	   &lworkl, rwork, &info);
    
   if ((ido==1)||(ido==-1)) av(n, workd+ipntr[0]-1, workd+ipntr[1]-1);
  } while ((ido==1)||(ido==-1));

  /* From those results, the eigenvalues and vectors are extracted. */

  if (info<0)
  {
    cerr << "Error with znaupd, info = " << info << "\n";
    cerr << "Check documentation in dsaupd\n\n";
  }
  else
  {
    zneupd_(&rvec, &howmny, select, d, v, &ldv, &sigma, workev,
	    &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv,
	    iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

    if (ierr!=0)
    {
      cerr << "Error with zneupd, info = " << ierr << "\n";
      cerr << "Check the documentation of zneupd.\n\n";
    }
    else if (info==1)
    {
      cerr << "Maximum number of iterations reached.\n\n";
    }
    else if (info==3)
    {
      cerr << "No shifts could be applied during implicit\n";
      cerr << "Arnoldi update, try increasing NCV.\n\n";
    }

    /* Before exiting, we copy the solution information over to the arrays of
       the calling program */
    complex<double> temp;
    for(int i=0; i<nev; i++) for (int j=i; j<nev; j++)
      if (d[j].real() < d[i].real())
      {
        temp = d[j];
        d[j] = d[i];
        d[i] = temp;
        for (int k=0; k<n; k++)
        {
          temp = v[i*n+k];
          v[i*n+k] = v[j*n+k];
          v[j*n+k] = temp;
        }
      }
   
    for(int i=0;i<nev;i++) if(abs(d[i].imag())>1e-6)
      cerr << "Error in znaupd: eigenvalues not real. " << d[i] << endl;
    for (int i=0; i<nev; i++) Evals[i] = d[i].real();
    for (int i=0; i<n; i++) Evecs[i] = v[i];

    delete resid;
    delete v;
    delete workd;
    delete workl;
    delete select;
    delete d;
  }
}
