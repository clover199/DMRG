#ifndef _OPERATORS_
#define _OPERATORS_

#include "tensor.h"


tensor fermion_c();

tensor fermion_c_up();

tensor fermion_c_down();

// *************************************
// index: alpha phy phy alpha

// spinless fermion Hamiltonian:
// T c^d_{i} c_{i+1} + P c_{i} c_{i+1} + h.c. + U c^d_{i] c_{i}
// Where T is the hopping parameter and P is the pairing parameter
tensor H_spinless(double t, double p, double u);

#endif
