
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

#ifdef FERMION
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
  ret.add_sign(-1, 1);
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
  ret.add_sign(-1, 1);
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
  ret.add_sign(-1, 1);
  return ret;
}
#endif

// ***********************
// *  2       3       2  *
// *  |       |       |  *
// *  O-1   0-O-2   1-O  *
// *  |       |       |  *
// *  0       1       0  *
// ***********************

// ******************************************************************
// Ising model: J Sz_{i} Sz_{i+1} + H Sx{i} (no symmetry considered)
//
//  MPO: I     0     0      Sz: 1  0     Sx: 0 1
//       Sz    0     0          0 -1         1 0
//      H*Sx  J*Sz   I

qtensor<double> H_Ising(double j, double h)
{
  tensor<double> x(3,2,3,2);
  x.update( 1, 0,0,0,0);
  x.update( 1, 0,1,0,1);
  x.update( 1, 1,0,0,0);
  x.update(-1, 1,1,0,1);
  x.update( h, 2,0,0,1);
  x.update( h, 2,1,0,0);
  x.update( j, 2,0,1,0);
  x.update(-j, 2,1,1,1);
  x.update( 1, 2,0,2,0);
  x.update( 1, 2,1,2,1);

  qtensor<double> ret(1,1,1,1);
  ret.update(x, 0,0,0,0);
  return ret;
}

//  MPO: H*Sx  J*Sz  I
qtensor<double> H_Ising_ledge(double j, double h)
{
  tensor<double> x(2,3,2);
  x.update( h, 0,0,1);
  x.update( h, 1,0,0);
  x.update( j, 0,1,0);
  x.update(-j, 1,1,1);
  x.update( 1, 0,2,0);
  x.update( 1, 1,2,1);

  qtensor<double> ret(1,1,1);
  ret.update(x, 0,0,0);
  return ret;
}

//  MPO: I  Sz  H*Sx
qtensor<double> H_Ising_redge(double j, double h)
{
  tensor<double> x(2,3,2);
  x.update( 1, 0,0,0);
  x.update( 1, 1,0,1);
  x.update( 1, 0,1,0);
  x.update(-1, 1,1,1);
  x.update( h, 0,2,1);
  x.update( h, 1,2,0);

  qtensor<double> ret(1,1,1);
  ret.update(x, 0,0,0);
  return ret;
}

// ******************************************************************
// 2-Potts model / Ising model: (with symmetry)
//  J Sx_{i} Sx_{i+1} + H Sz{i} 
//
//  MPO: I     0     0      Sz: 1  0     Sx: 0 1
//       Sx    0     0          0 -1         1 0
//      H*Sz  J*Sx   I

qtensor<double> H_2Potts(double j, double h)
{
  tensor<double> l00(3,1,3,1);  // I Sz
  l00.update( 1, 0,0,0,0);
  l00.update( h, 2,0,0,0);
  l00.update( 1, 2,0,2,0);
  tensor<double> l11(3,1,3,1);
  l11.update( 1, 0,0,0,0);
  l11.update(-h, 2,0,0,0);
  l11.update( 1, 2,0,2,0);
  tensor<double> l01(3,1,3,1);  // Sx
  l01.update( 1, 1,0,0,0);
  l01.update( j, 2,0,1,0);
  tensor<double> l10(3,1,3,1);
  l10.update( 1, 1,0,0,0);
  l10.update( j, 2,0,1,0);

  qtensor<double> ret(1,2,1,2);
  ret.update(l00, 0,0,0,0);
  ret.update(l01, 0,0,0,1);
  ret.update(l10, 0,1,0,0);
  ret.update(l11, 0,1,0,1);
  return ret;
}

qtensor<double> H_2Potts_ledge(double j, double h)
{
  tensor<double> l00(1,3,1);  // I Sz
  l00.update( h, 0,0,0);
  l00.update( 1, 0,2,0);
  tensor<double> l11(1,3,1);
  l11.update(-h, 0,0,0);
  l11.update( 1, 0,2,0);
  tensor<double> l01(1,3,1);  // Sx
  l01.update( j, 0,1,0);
  tensor<double> l10(1,3,1);
  l10.update( j, 0,1,0);

  qtensor<double> ret(2,1,2);
  ret.update(l00, 0,0,0);
  ret.update(l01, 0,0,1);
  ret.update(l10, 1,0,0);
  ret.update(l11, 1,0,1);
  return ret;
}

qtensor<double> H_2Potts_redge(double j, double h)
{
  tensor<double> l00(1,3,1);  // I Sz
  l00.update( 1, 0,0,0);
  l00.update( h, 0,2,0);
  tensor<double> l11(1,3,1);
  l11.update( 1, 0,0,0);
  l11.update(-h, 0,2,0);
  tensor<double> l01(1,3,1);  // Sx
  l01.update( 1, 0,1,0);
  tensor<double> l10(1,3,1);
  l10.update( 1, 0,1,0);

  qtensor<double> ret(2,1,2);
  ret.update(l00, 0,0,0);
  ret.update(l01, 0,0,1);
  ret.update(l10, 1,0,0);
  ret.update(l11, 1,0,1);
  return ret;
}

// ******************************************************************
// Three state Potts model (same as the paper):
// -fe^{ip} T_{j} - fe^{-ip} T^d_{j} - Je^{it} S_{j} S^d_{j+1} - Je^{-it} S^d_{j} S_{j+1}
//
//   MPO: I         0          0        0            T: 1  0  0     S: 0 1 0
//        S         0          0        0               0  w  0        0 0 1
//       S^d        0          0        0               0  0 w^2       1 0 0
//        h  -Je{-it}*S^d  -Je{it}*c    I
//
//    h = -f e^{ip} T_{j} - f e^{-ip} T^d_{j}

qtensor< complex<double> > H_Potts(double f, double j, double p, double t)
{
  tensor< complex<double> > l00(4,1,4,1);
  complex<double> val = -2*f*cos(p);
  l00.update(  1, 0,0,0,0);
  l00.update(val, 3,0,0,0);
  l00.update(  1, 3,0,3,0);
  tensor< complex<double> > l11(4,1,4,1);
  val = -2*f*cos(p+2*PI/3);
  l11.update(  1, 0,0,0,0);
  l11.update(val, 3,0,0,0);
  l11.update(  1, 3,0,3,0);
  tensor< complex<double> > l22(4,1,4,1);
  val = -2*f*cos(p+4*PI/3);
  l22.update(  1, 0,0,0,0);
  l22.update(val, 3,0,0,0);
  l22.update(  1, 3,0,3,0);
  val = -j*complex<double>(cos(t), sin(t));
  tensor< complex<double> > l01(4,1,4,1);  // S
  l01.update(  1, 1,0,0,0);
  l01.update(val, 3,0,2,0);
  tensor< complex<double> > l12(4,1,4,1);
  l12.update(  1, 1,0,0,0);
  l12.update(val, 3,0,2,0);
  tensor< complex<double> > l20(4,1,4,1);
  l20.update(  1, 1,0,0,0);
  l20.update(val, 3,0,2,0);
  val = -j*complex<double>(cos(t), -sin(t));
  tensor< complex<double> > l10(4,1,4,1);  // S^d
  l10.update(  1, 2,0,0,0);
  l10.update(val, 3,0,1,0);
  tensor< complex<double> > l21(4,1,4,1);
  l21.update(  1, 2,0,0,0);
  l21.update(val, 3,0,1,0);
  tensor< complex<double> > l02(4,1,4,1);
  l02.update(  1, 2,0,0,0);
  l02.update(val, 3,0,1,0);

  qtensor< complex<double> > ret(1,3,1,3);
  ret.update(l00, 0,0,0,0);
  ret.update(l01, 0,0,0,1);
  ret.update(l02, 0,0,0,2);
  ret.update(l10, 0,1,0,0);
  ret.update(l11, 0,1,0,1);
  ret.update(l12, 0,1,0,2);
  ret.update(l20, 0,2,0,0);
  ret.update(l21, 0,2,0,1);
  ret.update(l22, 0,2,0,2);
  return ret;
}

qtensor< complex<double> > H_Potts_ledge(double f, double j, double p, double t)
{
  tensor< complex<double> > l00(1,4,1);
  complex<double> val = -2*f*cos(p);
  l00.update(val, 0,0,0);
  l00.update(  1, 0,3,0);
  tensor< complex<double> > l11(1,4,1);
  val = -2*f*cos(p+2*PI/3);
  l11.update(val, 0,0,0);
  l11.update(  1, 0,3,0);
  tensor< complex<double> > l22(1,4,1);
  val = -2*f*cos(p+4*PI/3);
  l22.update(val, 0,0,0);
  l22.update(  1, 0,3,0);
  val = -j*complex<double>(cos(t), sin(t));
  tensor< complex<double> > l01(1,4,1);  // S
  l01.update(val, 0,2,0);
  tensor< complex<double> > l12(1,4,1);
  l12.update(val, 0,2,0);
  tensor< complex<double> > l20(1,4,1);
  l20.update(val, 0,2,0);
  val = -j*complex<double>(cos(t), -sin(t));
  tensor< complex<double> > l10(1,4,1);  // S^d
  l10.update(val, 0,1,0);
  tensor< complex<double> > l21(1,4,1);
  l21.update(val, 0,1,0);
  tensor< complex<double> > l02(1,4,1);
  l02.update(val, 0,1,0);

  qtensor< complex<double> > ret(3,1,3);
  ret.update(l00, 0,0,0);
  ret.update(l01, 0,0,1);
  ret.update(l02, 0,0,2);
  ret.update(l10, 1,0,0);
  ret.update(l11, 1,0,1);
  ret.update(l12, 1,0,2);
  ret.update(l20, 2,0,0);
  ret.update(l21, 2,0,1);
  ret.update(l22, 2,0,2);
  return ret;
}

qtensor< complex<double> > H_Potts_redge(double f, double j, double p, double t)
{
  tensor< complex<double> > l00(1,4,1);
  complex<double> val = -2*f*cos(p);
  l00.update(  1, 0,0,0);
  l00.update(val, 0,3,0);
  tensor< complex<double> > l11(1,4,1);
  val = -2*f*cos(p+2*PI/3);
  l11.update(  1, 0,0,0);
  l11.update(val, 0,3,0);
  tensor< complex<double> > l22(1,4,1);
  val = -2*f*cos(p+4*PI/3);
  l22.update(  1, 0,0,0);
  l22.update(val, 0,3,0);
  tensor< complex<double> > l01(1,4,1);  // S
  l01.update(  1, 0,1,0);
  tensor< complex<double> > l12(1,4,1);
  l12.update(  1, 0,1,0);
  tensor< complex<double> > l20(1,4,1);
  l20.update(  1, 0,1,0);
  tensor< complex<double> > l10(1,4,1);  // S^d
  l10.update(  1, 0,2,0);
  tensor< complex<double> > l21(1,4,1);
  l21.update(  1, 0,2,0);
  tensor< complex<double> > l02(1,4,1);
  l02.update(  1, 0,2,0);

  qtensor< complex<double> > ret(3,1,3);
  ret.update(l00, 0,0,0);
  ret.update(l01, 0,0,1);
  ret.update(l02, 0,0,2);
  ret.update(l10, 1,0,0);
  ret.update(l11, 1,0,1);
  ret.update(l12, 1,0,2);
  ret.update(l20, 2,0,0);
  ret.update(l21, 2,0,1);
  ret.update(l22, 2,0,2);
  return ret;
}

#ifdef FERMION
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
  ret.add_sign(0, -1, 0, 1);
  return ret;
}

//  MPO: U*n  T*c^d  P*c -P*c^d -T*c  I
qtensor<double> H_spinless_fermion_ledge(double t, double p, double u)
{
  tensor<double> lu(1,6,1);
  lu.update(1, 0,5,0);
  tensor<double> rd(1,6,1);
  rd.update(u, 0,0,0);
  rd.update(1, 0,5,0);
  tensor<double> ld(1,6,1); // c^d non-zero
  ld.update( t, 0,1,0);
  ld.update(-p, 0,3,0);
  tensor<double> ru(1,6,1); // c non-zero
  ru.update( p, 0,2,0);
  ru.update(-t, 0,4,0);

  qtensor<double> ret(2,1,2);
  ret.update(lu, 0,0,0);
  ret.update(ru, 0,0,1);
  ret.update(ld, 1,0,0);
  ret.update(rd, 1,0,1);
  ret.add_sign(-1, 0, 1);
  return ret;
}

//   MPO: I  c  c  c^d  c^d  U*n 
qtensor<double> H_spinless_fermion_redge(double t, double p, double u)
{
  tensor<double> lu(1,6,1);
  lu.update(1, 0,0,0);
  tensor<double> rd(1,6,1);
  rd.update(1, 0,0,0);
  rd.update(u, 0,5,0);
  tensor<double> ld(1,6,1); // c^d non-zero
  ld.update( 1, 0,3,0);
  ld.update( 1, 0,4,0);
  tensor<double> ru(1,6,1); // c non-zero
  ru.update( 1, 0,1,0);
  ru.update( 1, 0,2,0);

  qtensor<double> ret(2,1,2);
  ret.update(lu, 0,0,0);
  ret.update(ru, 0,0,1);
  ret.update(ld, 1,0,0);
  ret.update(rd, 1,0,1);
  ret.add_sign(-1, 0, 1);
  return ret;
}
#endif