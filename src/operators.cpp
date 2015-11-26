
#include "global.h"
#include "tensor.h"

// spinless fermion annihilation operator c
// basis: 0 1.  symmetry: 0 1.
// 0 1
// 0 0
tensor fermion_c()
{
  dtensor<double> ld(1,1);
  dtensor<double> ru(1,1); ru.update(1, 0,0);
  qtensor< dtensor<double> > ret(2,2);
  ret.update(ld, 1,0);
  ret.update(ru, 0,1);
  return tensor(ret);
}


// spinfull fermion annihilation operator c_up
// basis: 0 2 u d.  symmetry: 0 1.
// 0 0 1 0
// 0 0 0 0
// 0 0 0 0
// 0 1 0 0
tensor fermion_c_up()
{
  dtensor<double> ld(2,2); ld.update(1, 1,1);
  dtensor<double> ru(2,2); ru.update(1, 0,0);
  qtensor< dtensor<double> > ret(2,2);
  ret.update(ld, 1,0);
  ret.update(ru, 0,1);
  return tensor(ret);
}


// spinfull fermion annihilation operator c_down
// basis: 0 2 u d.  symmetry: 0 1.
// 0 0 0 1
// 0 0 0 0
// 0 - 0 0
// 0 0 0 0
tensor fermion_c_down()
{
  dtensor<double> ld(2,2); ld.update(-1, 0,1);
  dtensor<double> ru(2,2); ru.update(1, 0,1);
  qtensor< dtensor<double> > ret(2,2);
  ret.update(ld, 1,0);
  ret.update(ru, 0,1);
  return tensor(ret);
}

// ***********************************************
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
// index: a p p a
tensor H_spinless(double t, double p, double u)
{
  dtensor<double> lu(6,1,1,6);
  lu.update(1, 0,0,0,0);
  lu.update(1, 5,0,0,5);
  dtensor<double> rd(6,1,1,6);
  rd.update(1, 0,0,0,0);
  rd.update(u, 5,0,0,0);
  rd.update(1, 5,0,0,5);
  dtensor<double> ld(6,1,1,6); // c^d non-zero
  ld.update( 1, 3,0,0,0);
  ld.update( 1, 4,0,0,0);
  ld.update( t, 5,0,0,1);
  ld.update(-p, 5,0,0,3);
  dtensor<double> ru(6,1,1,6); // c non-zero
  ru.update( 1, 1,0,0,0);
  ru.update( 1, 2,0,0,0);
  ru.update( p, 5,0,0,2);
  ru.update(-t, 5,0,0,4);
  qtensor< dtensor<double> > ret(1,2,2,1);
  ret.update(lu, 0,0,0,0);
  ret.update(rd, 0,1,1,0);
  ret.update(ld, 0,1,0,0);
  ret.update(ru, 0,0,1,0);
  return tensor(ret);
}

