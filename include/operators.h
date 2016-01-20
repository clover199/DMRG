#ifndef _MY_OPERATORS_
#define _MY_OPERATORS_

#include "qtensor.h"

// *******
// *  1  *
// *  |  *
// *  O  *
// *  |  *
// *  0  *
// *******

qtensor<double> sigma_x();  // no symmetry

qtensor< complex<double> > sigma_y();  // no symmetry

qtensor<double> sigma_z();  // no symmetry

qtensor<double> fermion_c();

qtensor<double> fermion_c_up();

qtensor<double> fermion_c_down();

// ***********
// *    3    *
// *    |    *
// *  0-O-2  *
// *    |    *
// *    1    *
// ***********

// ******************************************************************
// Ising model: J Sz_{i} Sz_{i+1} + H Sz{i}

qtensor<double> H_Ising(double j=1., double h=1.);

qtensor<double> H_Ising_ledge(double j=1., double h=1.);

qtensor<double> H_Ising_redge(double j=1., double h=1.);

// ******************************************************************
// spinless fermion Hamiltonian:
//   U=0: p-wave superconductor
//   P=0: Hubbard model
// T c^d_{i} c_{i+1} + P c_{i} c_{i+1} + h.c. + U c^d_{i] c_{i}
// Where T is the hopping parameter and P is the pairing parameter

qtensor<double> H_spinless_fermion(double t=1., double p=1., double u=1.);

#endif
