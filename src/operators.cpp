
#include "global.h"
#include "qtensor.h"


// 0 1  without symmetry
// 1 0
qtensor<double> sigma_x()
{
  tensor<double> x(2,2);
  x.update(1, 0,1);
  x.update(1, 1,0);
  qtensor<double> ret(1,1);
  ret.update(x, 0,0);
  return ret;
}


// 0 -i  without symmetry
// i  0
qtensor< complex<double> > sigma_y()
{
  tensor< complex<double> > x(2,2);
  x.update(complex<double>(0,-1), 0,1);
  x.update(complex<double>(0,1), 1,0);
  qtensor< complex<double> > ret(1,1);
  ret.update(x, 0,0);
  return ret;
}


// 1  0  without symmetry
// 0 -1
qtensor<double> sigma_z()
{
  tensor<double> x(2,2);
  x.update(1, 0,0);
  x.update(-1, 1,1);
  qtensor<double> ret(1,1);
  ret.update(x, 0,0);
  return ret;
}


// spinless fermion annihilation operator c
// basis: 0 1.  symmetry: 0 1.
// 0 1
// 0 0
qtensor<double> fermion_c()
{
  tensor<double> ld(1,1);
  tensor<double> ru(1,1); ru.update(1, 0,0);
  qtensor<double> ret(2,2);
  ret.update(ld, 1,0);
  ret.update(ru, 0,1);
  return ret;
}


// spinfull fermion annihilation operator c_up
// basis: 0 2 u d.  symmetry: 0 1.
// 0 0 1 0
// 0 0 0 0
// 0 0 0 0
// 0 1 0 0
qtensor<double> fermion_c_up()
{
  tensor<double> ld(2,2); ld.update(1, 1,1);
  tensor<double> ru(2,2); ru.update(1, 0,0);
  qtensor<double> ret(2,2);
  ret.update(ru, 0,1);
  ret.update(ld, 1,0);
  return ret;
}


// spinfull fermion annihilation operator c_down
// basis: 0 2 u d.  symmetry: 0 1.
// 0 0 0 1
// 0 0 0 0
// 0 - 0 0
// 0 0 0 0
qtensor<double> fermion_c_down()
{
  tensor<double> ld(2,2); ld.update(-1, 0,1);
  tensor<double> ru(2,2); ru.update(1, 0,1);
  qtensor<double> ret(2,2);
  ret.update(ru, 0,1);
  ret.update(ld, 1,0);
  return ret;
}

// ***********
// *    3    *
// *    |    *
// *  0-O-2  *
// *    |    *
// *    1    *
// ***********

// ******************************************************************
// Ising model: J Sz_{i} Sz_{i+1} + H Sz{i}
//
//  MPO: I     0     0      Sz: 1  0
//       Sz    0     0          0 -1
//      H*Sz  J*Sz   I

qtensor<double> H_Ising(double j, double h)
{
  tensor<double> x(3,2,3,2);
  x.update( 1, 0,0,0,0);
  x.update( 1, 0,1,0,1);
  x.update( 1, 1,0,0,0);
  x.update(-1, 1,1,0,1);
  x.update( h, 2,0,0,0);
  x.update(-h, 2,1,0,1);
  x.update( j, 2,0,1,0);
  x.update(-j, 2,1,1,1);
  x.update( 1, 2,0,2,0);
  x.update(-1, 2,1,2,1);

  qtensor<double> ret(1,1,1,1);
  ret.update(x, 0,0,0,0);
  return ret;
}

// ******************************************************************
// spinless fermion Hamiltonian:
// T c^d_{i} c_{i+1} + P c_{i} c_{i+1} + h.c. + U c^d_{i] c_{i}
//
//   MPO: I     0     0     0     0     0         c: 0 1     c^d: 0 0
//        c     0     0     0     0     0            0 0          1 0
//        c     0     0     0     0     0
//       c^d    0     0     0     0     0
//       c^d    0     0     0     0     0
//       U*n  T*c^d  P*c -P*c^d -T*c    I
//
// Where T is the hopping parameter and P is the pairing parameter

qtensor<double> H_spinless_fermion(double t, double p, double u)
{
  tensor<double> lu(6,1,6,1);
  lu.update(1, 0,0,0,0);
  lu.update(1, 5,0,5,0);
  tensor<double> rd(6,1,6,1);
  rd.update(1, 0,0,0,0);
  rd.update(u, 5,0,0,0);
  rd.update(1, 5,0,5,0);
  tensor<double> ld(6,1,6,1); // c^d non-zero
  ld.update( 1, 3,0,0,0);
  ld.update( 1, 4,0,0,0);
  ld.update( t, 5,0,1,0);
  ld.update(-p, 5,0,3,0);
  tensor<double> ru(6,1,6,1); // c non-zero
  ru.update( 1, 1,0,0,0);
  ru.update( 1, 2,0,0,0);
  ru.update( p, 5,0,2,0);
  ru.update(-t, 5,0,4,0);

  qtensor<double> ret(1,2,1,2);
  ret.update(lu, 0,0,0,0);
  ret.update(ru, 0,0,0,1);
  ret.update(ld, 0,1,0,0);
  ret.update(rd, 0,1,0,1);
  return ret;
}

