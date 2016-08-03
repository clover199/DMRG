
#include "global.h"
#include "qtensor.h"


// basis: 0 1
// 0 0
// 0 1
qtensor<double> number() // number operator
{
  tensor<double> t00(1,1,1,1);
  tensor<double> t01(1,1,1,1);
  tensor<double> t10(1,1,1,1);
  tensor<double> t11(1,1,1,1);
  t11.update(1, 0,0,0,0);
  qtensor<double> ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
#ifdef FERMION
  ret.add_sign(0, -1, 0, 1);
#endif
  return ret;
}


#ifdef FERMION
// basis: 0 1
// 0 1
// 0 0
qtensor<double> fermion_c() // spinless fermion annihilation operator
{
  tensor<double> t00(1,1,1,1);
  tensor<double> t01(1,1,1,1);
  tensor<double> t10(1,1,1,1);
  tensor<double> t11(1,1,1,1);
  t01.update(1, 0,0,0,0);
  qtensor<double> ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  ret.add_sign(0, -1, 0, 1);
  return ret;
}


// basis: 0 2 u d
// 0 0 1 0
// 0 0 0 0
// 0 0 0 0
// 0 1 0 0
qtensor<double> fermion_c_up() // spin up fermion annihilation operator
{
  tensor<double> t00(1,2,1,2);
  tensor<double> t01(1,2,1,2);
  tensor<double> t10(1,2,1,2);
  tensor<double> t11(1,2,1,2);
  t01.update(1, 0,0,0,0);
  t10.update(1, 0,1,0,1);
  qtensor<double> ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  ret.add_sign(0, -1, 0, 1);
  return ret;
}


// basis: 0 2 u d
// 0 0 0 1
// 0 0 0 0
// 0 - 0 0
// 0 0 0 0
qtensor<double> fermion_c_down() // spin down fermion annihilation operator
{
  tensor<double> t00(1,2,1,2);
  tensor<double> t01(1,2,1,2);
  tensor<double> t10(1,2,1,2);
  tensor<double> t11(1,2,1,2);
  t01.update( 1, 0,0,0,1);
  t10.update(-1, 0,0,0,1);
  qtensor<double> ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  ret.add_sign(0, -1, 0, 1);
  return ret;
}


// basis: 0 2 u d
// 0 - 0 0
// 0 0 0 0
// 0 0 0 0
// 0 0 0 0
qtensor<double> fermion_pair() // cooper pair annihilation operator
{
  tensor<double> t00(1,2,1,2);
  tensor<double> t01(1,2,1,2);
  tensor<double> t10(1,2,1,2);
  tensor<double> t11(1,2,1,2);
  t00.update(-1, 0,0,0,1);
  qtensor<double> ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  ret.add_sign(0, -1, 0, 1);
  return ret;
}


// basis: 0 2 u d
// 0 0 0 0
// 0 1 0 0
// 0 0 1 0
// 0 0 0 0
qtensor<double> fermion_n_up() // spin up number operator
{
  tensor<double> t00(1,2,1,2);
  tensor<double> t01(1,2,1,2);
  tensor<double> t10(1,2,1,2);
  tensor<double> t11(1,2,1,2);
  t00.update(1, 0,1,0,1);
  t11.update(1, 0,0,0,0);
  qtensor<double> ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  ret.add_sign(0, -1, 0, 1);
  return ret;
}

// basis: 0 2 u d
// 0 0 0 0
// 0 1 0 0
// 0 0 0 0
// 0 0 0 1
qtensor<double> fermion_n_down() // spin down number operator
{
  tensor<double> t00(1,2,1,2);
  tensor<double> t01(1,2,1,2);
  tensor<double> t10(1,2,1,2);
  tensor<double> t11(1,2,1,2);
  t00.update(1, 0,1,0,1);
  t11.update(1, 0,1,0,1);
  qtensor<double> ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  ret.add_sign(0, -1, 0, 1);
  return ret;
}

// basis: 0 2 u d
// 0 0 0 0
// 0 1 0 0
// 0 0 0 0
// 0 0 0 0
qtensor<double> fermion_n_pair() // cooper pair number operator
{
  tensor<double> t00(1,2,1,2);
  tensor<double> t01(1,2,1,2);
  tensor<double> t10(1,2,1,2);
  tensor<double> t11(1,2,1,2);
  t00.update(1, 0,1,0,1);
  qtensor<double> ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  ret.add_sign(0, -1, 0, 1);
  return ret;
}

qtensor<double> fermion_parity(int n) // fermion parity operator
{
  n = n/2;
  tensor<double> t00(1,n,1,n);
  tensor<double> t01(1,n,1,n);
  tensor<double> t10(1,n,1,n);
  tensor<double> t11(1,n,1,n);
  for(int i=0;i<n;i++) t00.update( 1, 0,i,0,i);
  for(int i=0;i<n;i++) t11.update(-1, 0,i,0,i);
  qtensor<double> ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  ret.add_sign(0, -1, 0, 1);
  return ret;
}

// basis: 0 2 u d
// 1 0 0 0
// 0 - 0 0
// 0 0 0 0
// 0 0 0 0
qtensor<double> pair_parity() // cooper-pair parity operator
{
  tensor<double> t00(1,2,1,2);
  tensor<double> t01(1,2,1,2);
  tensor<double> t10(1,2,1,2);
  tensor<double> t11(1,2,1,2);
  t00.update( 1, 0,0,0,0);
  t00.update(-1, 0,1,0,1);
  qtensor<double> ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  ret.add_sign(0, -1, 0, 1);
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
// Ising model (without symmetry): H Sx_{i} + J Sz_{i} Sz_{i+1}
//
//  MPO: I     0     0      Sz: 1  0     Sx: 0 1
//       Sz    0     0          0 -1         1 0
//      H*Sx  J*Sz   I

qtensor<double> H_Ising(double h, double j)
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

// ******************************************************************
// 2-Potts model / Ising model: H Sz_{i} + J Sx_{i} Sx_{i+1}
//
//  MPO: I     0     0      Sz: 1  0     Sx: 0 1
//       Sx    0     0          0 -1         1 0
//      H*Sz  J*Sx   I

qtensor<double> H_Potts2(double h, double j)
{ 
  tensor<double> t00(3,1,3,1);  // I Sz
  t00.update( 1, 0,0,0,0);
  t00.update( h, 2,0,0,0);
  t00.update( 1, 2,0,2,0);
  tensor<double> t11(3,1,3,1);
  t11.update( 1, 0,0,0,0);
  t11.update(-h, 2,0,0,0);
  t11.update( 1, 2,0,2,0);
  tensor<double> t01(3,1,3,1);  // Sx
  t01.update( 1, 1,0,0,0);
  t01.update( j, 2,0,1,0);
  tensor<double> t10(3,1,3,1);
  t10.update( 1, 1,0,0,0);
  t10.update( j, 2,0,1,0);

  qtensor<double> ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  return ret;
}

// ******************************************************************
// 3-state chiral Potts model (for the paper):
//   -fe^{ip} T_{j} - fe^{-ip} T^d_{j} - Je^{it} S_{j} S^d_{j+1} - Je^{-it} S^d_{j} S_{j+1}
//
//   MPO: I         0          0        0            T: 1  0  0     S: 0 1 0
//        S         0          0        0               0  w  0        0 0 1
//       S^d        0          0        0               0  0 w^2       1 0 0
//        h  -Je{-it}*S^d  -Je{it}*S    I
//
//    h = -f e^{ip} T_{j} - f e^{-ip} T^d_{j}

qtensor< complex<double> > H_Potts3(double f, double j, double p, double t)
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

// ******************************************************************
// XYZ Chain: Jx Sx_{i} Sx_{i+1} + Jy Sy_{i} Sy_{i+1} + Jz Sz_{i} Sz_{i+1}
//
//  MPO:  I     0      0      0      0   Sx: 0 1     Sy: 0 -i    Sz: 1  0
//       Sx     0      0      0      0       1 0         i  0        0 -1
//       Sy     0      0      0      0
//       Sz     0      0      0      0
//        0   Jx*Sx  Jy*Sy  Jz*Sz    I

qtensor<complex<double> > H_XYZ(double Jx, double Jy, double Jz)
{
  complex<double> i = complex<double>(0,1);
  tensor<complex<double> > t00(5,1,5,1);  // I, Sz
  t00.update(    1, 0,0,0,0);
  t00.update(    1, 3,0,0,0);
  t00.update(   Jz, 4,0,3,0);
  t00.update(    1, 4,0,4,0);
  tensor<complex<double> > t01(5,1,5,1);  // Sx, Sy
  t01.update(    1, 1,0,0,0);
  t01.update(   -i, 2,0,0,0);
  t01.update(   Jx, 4,0,1,0);
  t01.update(-i*Jy, 4,0,2,0);
  tensor<complex<double> > t10(5,1,5,1);  // Sx, Sy
  t10.update(    1, 1,0,0,0);
  t10.update(    i, 2,0,0,0);
  t10.update(   Jx, 4,0,1,0);
  t10.update( i*Jy, 4,0,2,0);
  tensor<complex<double> > t11(5,1,5,1);  // I, Sz
  t11.update(    1, 0,0,0,0);
  t11.update(   -1, 3,0,0,0);
  t11.update(  -Jz, 4,0,3,0);
  t11.update(    1, 4,0,4,0);

  qtensor<complex<double> > ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  return ret;
}


#ifdef FERMION
// ******************************************************************
// Kitaev chain: T c^d_{i} c_{i+1} - P c_{i} c_{i+1} + h.c. - 2U c^d_{i] c_{i}
//
//   MPO: I        0          0         0         c: 0 1     c^d: 0 0
//        c        0          0         0            0 0          1 0
//       c^d       0          0         0
//      -2U*n   T*c^d-P*c  P*c^d-T*c    I
//
// Where T is the hopping parameter and P is the pairing parameter


qtensor<double> H_Kitaev(double t, double p, double u)
{
  tensor<double> t00(4,1,4,1); // I
  t00.update(   1, 0,0,0,0);
  t00.update(   1, 3,0,3,0);
  tensor<double> t11(4,1,4,1); // I, n
  t11.update(   1, 0,0,0,0);
  t11.update(-2*u, 3,0,0,0);
  t11.update(   1, 3,0,3,0);
  tensor<double> t01(4,1,4,1); // c
  t01.update(   1, 1,0,0,0);
  t01.update(  -p, 3,0,1,0);
  t01.update(  -t, 3,0,2,0);
  tensor<double> t10(4,1,4,1); // c^d
  t10.update(   1, 2,0,0,0);
  t10.update(   t, 3,0,1,0);
  t10.update(   p, 3,0,2,0);

  qtensor<double> ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  ret.add_sign(0, -1, 0, 1);
  return ret;
}

qtensor<complex<double> > H_Kitaev(complex<double> t, complex<double> p, double u)
{
// T c^d_{i} c_{i+1} - P c_{i} c_{i+1}
//
//   MPO: I         0          0         0         c: 0 1     c^d: 0 0
//        c         0          0         0            0 0          1 0
//       c^d        0          0         0
//      -2U*n   T*c^d-P*c Pz*c^d-Tz*c    I

  complex<double> tz = conj(t);
  complex<double> pz = conj(p);
  tensor<complex<double> > t00(4,1,4,1); // I
  t00.update(   1, 0,0,0,0);
  t00.update(   1, 3,0,3,0);
  tensor<complex<double> > t11(4,1,4,1); // I, n
  t11.update(   1, 0,0,0,0);
  t11.update(-2*u, 3,0,0,0);
  t11.update(   1, 3,0,3,0);
  tensor<complex<double> > t01(4,1,4,1); // c
  t01.update(   1, 1,0,0,0);
  t01.update(  -p, 3,0,1,0);
  t01.update( -tz, 3,0,2,0);
  tensor<complex<double> > t10(4,1,4,1); // c^d
  p = conj(p);
  t10.update(   1, 2,0,0,0);
  t10.update(   t, 3,0,1,0);
  t10.update(  pz, 3,0,2,0);

  qtensor<complex<double> > ret(1,2,1,2);
  ret.update(t00, 0,0,0,0);
  ret.update(t01, 0,0,0,1);
  ret.update(t10, 0,1,0,0);
  ret.update(t11, 0,1,0,1);
  ret.add_sign(0, -1, 0, 1);
  return ret;
}

#endif