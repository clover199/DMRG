#ifndef _BASIC_
#define _BASIC_

void tgemm(char transa, char transb, int m, int n, int K,
           double *a, double *b, double *C, double alpha, double beta);

void tgesvd(double *a, int m, int n, double *S, double *u, double *v);

void dsaupd(int n, int nev, double *Evals, double *Evecs, void av(int n, double *in, double *out));

void dsyev(double *a, int N, double *W);
void dsyev(double *a, int N, double *W, double *vec);

#include <complex>

void tgemm(char transa, char transb, int m, int n, int K,
           std::complex<double> *a, std::complex<double> *b, std::complex<double> *C,
           std::complex<double> alpha, std::complex<double> beta);

void tgesvd(std::complex<double> *a, int m, int n, double *S,
            std::complex<double> *u, std::complex<double> *v);

void znaupd(int n, int nev, double *Evals, std::complex<double> *Evecs,
            void av(int n, complex<double> *in, complex<double> *out));

void zheev(std::complex<double> *a, int N, double *W);
void zheev(std::complex<double> *a, int N, double *W, std::complex<double> *vec);

#endif
