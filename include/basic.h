#ifndef _BASIC_
#define _BASIC_

#include "val_for_lanczos.h"


// dgemm
void zgemm(char transa, char transb, int m, int n, int K,
           double *a, double *b, double *C, double alpha, double beta);

// dgesvd
void zgesvd(double *a, int m, int n, double *S, double *u, double *v);

// dsaupd
void znaupd(int n, int nev, double *Evals, double *Evecs, lanczos<double>& pass_val);

//dsyev
void zheev(double *a, int N, double *W);
void zheev(double *a, int N, double *W, double *vec);


using std::complex;
#include <complex>

void zgemm(char transa, char transb, int m, int n, int K,
           complex<double> *a, complex<double> *b, complex<double> *C,
           complex<double> alpha, complex<double> beta);

void zgesvd(complex<double> *a, int m, int n, double *S,
            complex<double> *u, std::complex<double> *v);

void znaupd(int n, int nev, double *Evals, complex<double> *Evecs,
            lanczos<complex<double> >& pass_val);

void zheev(complex<double> *a, int N, double *W);
void zheev(complex<double> *a, int N, double *W, complex<double> *vec);

#endif
