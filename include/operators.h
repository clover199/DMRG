#ifndef _MY_OPERATORS_
#define _MY_OPERATORS_

#include "qtensor.h"

// ***********
// *    3    *
// *    |    *
// *  0-O-2  *
// *    |    *
// *    1    *
// ***********

qtensor<double> number();

#ifdef FERMION
qtensor<double> fermion_c();

qtensor<double> fermion_c_up();

qtensor<double> fermion_c_down();

qtensor<double> fermion_pair();

qtensor<double> fermion_n_up();

qtensor<double> fermion_n_down();

qtensor<double> fermion_n_pair();

qtensor<double> fermion_parity(int n=2);

qtensor<double> pair_parity();
#endif

// ******************************************************************
// Ising model (without symmetry): H Sx_{i} + J Sz_{i} Sz_{i+1}

qtensor<double> H_Ising(double h, double j);

// ******************************************************************
// 2-Potts model / Ising model: H Sz_{i} + J Sx_{i} Sx_{i+1}

qtensor<double> H_Potts2(double h, double j);

// ******************************************************************
// 3-state chiral Potts model (for the paper):
//   p=0, t=0: three state Potts model
//   - f e^{-ip} T^d_{j} - Je^{-it} S^d_{j} S_{j+1}

qtensor< complex<double> > H_Potts3(double f, double j, double p=0, double t=0);

// ******************************************************************
// XYZ Chain: Jx Sx_{i} Sx_{i+1} + Jy Sy_{i} Sy_{i+1} + Jz Sz_{i} Sz_{i+1}
//   spin 1/2

qtensor<complex<double> > H_XYZ(double Jx, double Jy, double Jz);

// ******************************************************************
// Kitaev chain: - i u a_{j} b_{j} - i v b_{j] a_{j+1}
//   u=2U, v=T=P:
//   T c^d_{i} c_{i+1} - P c_{i} c_{i+1} + h.c. - 2U c^d_{i] c_{i}

#ifdef FERMION
qtensor<double> H_Kitaev(double t, double p, double u);
qtensor<complex<double> > H_Kitaev(complex<double> t, complex<double> p, double u);
#endif

#endif
