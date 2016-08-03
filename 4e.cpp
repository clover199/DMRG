/* This is the main function to calculate the 4e superconductor (paper)
   with I-breaking hopping terms (the 2nd phase diagram).
   Periodic boundary condition (PBC) can be applied to perturbation only. */


#include "global.h"
#include "functions.h"

qtensor<complex<double> > H_4e1(double t1, double t3, double f0,
                                double f1, double p1, double frac=0);

qtensor<complex<double> > H_4e2(double t1, double t3, double f0,
                                double f1, double p1, double frac=0);

qtensor<complex<double> > c_up1();
qtensor<complex<double> > c_down1();
qtensor<complex<double> > c_up2();
qtensor<complex<double> > c_down2();
qtensor<complex<double> > c_pair1();
qtensor<complex<double> > c_pair2();
qtensor<complex<double> > c_pairpair();
qtensor<complex<double> > c_pair_dagger1();
qtensor<complex<double> > c_pair_dagger2();

int main(int argc, char *argv[])
{

  if(SYMMETRY-8)
  {
    cerr << "Define the SYMMETRY in global.h as 8." << endl;
    return 0;
  }

  time_t start, end;
  time(&start);

  // ********** initialize the system **********
  vector<string> para_name; // name of the parameters
  //                          0    1    2    3
  para_name = set_para_name("t1","t3","f0","f1","frac");

  int sites = 4;
  int cutoff = 10;
  int sweep = -1; 
  int symmetry_sector = 0;
  
  vector<double> para;
  string filename = "";
  para = set_para_val(argc, argv, sites, cutoff, sweep, symmetry_sector,
                      para_name, filename);
  if(sites%4!=0) cerr << "Error: invalid size l=" << sites << endl;
  sites /= 2;
  
  mpo<complex<double> > my_mpo(sites);
  for(int i=0;i<sites;i++)
  {
    if(i%2==0) my_mpo[i] = H_4e1(para[0], para[1], para[2], para[3], para[3], para[4]*i/sites);
    else my_mpo[i] = H_4e2(para[0], para[1], para[2], para[3], para[3], para[4]*i/sites);
  }
  my_mpo[0] = H_4e1(para[0], para[1], para[2], para[3], para[3], 0).left();
  my_mpo[sites-1] = H_4e2(para[0], para[1], para[2], para[3], para[3], para[4]*(sites-1)/sites).right();
  mps<complex<double> > my_mps(my_mpo);

#ifdef PBC
  cout << "Using periodic boundary condition\n" << endl;
  my_mps.add_edge(H_4e1(0, -para[1], 0, para[3], para[3], para[4]*(sites-1)/sites).right(),
                  H_4e2(0, -para[1], 0, para[3], para[3], para[4]*(sites-1)/sites).left() );
#endif

  cout.precision(10);
  dmrg(my_mps, my_mpo, cutoff, sweep, symmetry_sector, filename);
  my_mps.prep_calc();

  cout << "\n********** Calculate expectation values **********\n\n";
  ofstream data;
  filename = "data"+ filename;
  data.open(filename.c_str());

  qtensor<complex<double> > u1 = c_up1();
  qtensor<complex<double> > u1_d = u1.exchange(1,3);
  u1_d.conjugate();
  qtensor<complex<double> > d1 = c_down1();
  qtensor<complex<double> > d1_d = d1.exchange(1,3);
  d1_d.conjugate();
  qtensor<complex<double> > u2 = c_up2();
  qtensor<complex<double> > u2_d = u2.exchange(1,3);
  u2_d.conjugate();
  qtensor<complex<double> > d2 = c_down2();
  qtensor<complex<double> > d2_d = d2.exchange(1,3);
  d2_d.conjugate();
  qtensor<complex<double> > p1 = c_pair1();
  qtensor<complex<double> > p1_d = c_pair_dagger1();
  qtensor<complex<double> > p2 = c_pair2();
  qtensor<complex<double> > p2_d = c_pair_dagger2();
  qtensor<complex<double> > pp = c_pairpair();
  qtensor<complex<double> > id = p1.id();
  complex<double> result;
  
  // ********** correlation function ***********
  my_mpo.resize(sites); 

  for(int i=0;i<sites;i++) my_mpo[i] = id;
  my_mpo.edge(p1_d.left(), p2.right());
  result = calc(my_mps, my_mpo);
  cout << "L1 R1: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

  my_mpo.edge(p1_d.left(), p1.right());
  result = calc(my_mps, my_mpo);
  cout << "L1 R2: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

  my_mpo.edge(p2_d.left(), p2.right());
  result = calc(my_mps, my_mpo);
  cout << "L2 R1: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

  my_mpo.edge(p2_d.left(), p1.right());
  result = calc(my_mps, my_mpo);
  cout << "L2 R2: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

  my_mpo.edge(u1_d.left(), u2.right());
  result = calc(my_mps, my_mpo);
  cout << "up: L1 R1: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

  my_mpo.edge(u1_d.left(), u1.right());
  result = calc(my_mps, my_mpo);
  cout << "up: L1 R2: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

  my_mpo.edge(u2_d.left(), u2.right());
  result = calc(my_mps, my_mpo);
  cout << "up: L2 R1: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

  my_mpo.edge(u2_d.left(), u1.right());
  result = calc(my_mps, my_mpo);
  cout << "up: L2 R2: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

  my_mpo.edge(d1_d.left(), d2.right());
  result = calc(my_mps, my_mpo);
  cout << "down: L1 R1: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

  my_mpo.edge(d1_d.left(), d1.right());
  result = calc(my_mps, my_mpo);
  cout << "down: L1 R2: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

  my_mpo.edge(d2_d.left(), d2.right());
  result = calc(my_mps, my_mpo);
  cout << "down: L2 R1: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

  my_mpo.edge(d2_d.left(), d1.right());
  result = calc(my_mps, my_mpo);
  cout << "down: L2 R2: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

  result = calc(my_mps);
  cout << "Normalization: " << result << endl;
  data << sites << "\t" << result.real() << "\t" << result.imag() << endl;

for(int p=1;p<sites-1;p++)
{
  for(int i=0;i<sites;i++) my_mpo[i] = id;
  my_mpo.edge(id.left(), id.right());
  my_mpo[p] = p1;
  result = calc(my_mps, my_mpo);
  cout << "Cooper-pair: " << 2*p+1 << " \t " << result << endl;
  data << 2*p+1 << "\t" << result.real() << "\t" << result.imag() << endl;
  my_mpo[p] = p2;
  result = calc(my_mps, my_mpo);
  cout << "Cooper-pair: " << 2*p+2 << " \t " << result << endl;
  data << 2*p+2 << "\t" << result.real() << "\t" << result.imag() << endl;
}

for(int p=1;p<sites-1;p++)
{
  for(int i=0;i<sites;i++) my_mpo[i] = id;
  my_mpo.edge(id.left(), id.right());
  my_mpo[p] = pp;
  result = calc(my_mps, my_mpo);
  cout << "pair-pair: " << 2*p+1.5 << " \t " << result << endl;
  data << 2*p+1.5 << "\t" << result.real() << "\t" << result.imag() << endl;
if(p<sites-2)
{
  my_mpo[p] = p2;
  my_mpo[p+1] = p1;
  result = calc(my_mps, my_mpo);
  cout << "pair-pair: " << 2*p+2.5 << " \t " << result << endl;
  data << 2*p+2.5 << "\t" << result.real() << "\t" << result.imag() << endl;
}
}

  time(&end);
  cout << "Time used: " << difftime(end,start) << endl;
}

namespace gg
{
  int lcu[32][5] = { // left c_up
	  {0,4,0,0, 1},
	  {0,4,1,1, 1},
	  {0,5,0,0, 1},
	  {0,5,1,1,-1},
	  {1,4,0,0, 1},
	  {1,4,1,1,-1},
	  {1,5,0,0, 1},
	  {1,5,1,1, 1},
	  {2,6,0,0, 1},
	  {2,6,1,1,-1},
	  {2,7,0,0, 1},
	  {2,7,1,1, 1},
	  {3,6,0,0, 1},
	  {3,6,1,1,-1},
	  {3,7,0,0, 1},
	  {3,7,1,1, 1},
	  {4,0,0,0, 1},
	  {4,0,1,1, 1},
	  {4,1,0,0,-1},
	  {4,1,1,1, 1},
	  {5,0,0,0,-1},
	  {5,0,1,1, 1},
	  {5,1,0,0, 1},
	  {5,1,1,1, 1},
	  {6,2,0,0, 1},
	  {6,2,1,1,-1},
	  {6,3,0,0,-1},
	  {6,3,1,1, 1},
	  {7,2,0,0,-1},
	  {7,2,1,1,-1},
	  {7,3,0,0, 1},
	  {7,3,1,1, 1} };
 
  int lcd[32][5] = { // left c_down
	  {0,2,0,1, 1},
	  {0,2,1,0, 1},
	  {0,3,0,1, 1},
	  {0,3,1,0,-1},
	  {1,2,0,1, 1},
	  {1,2,1,0, 1},
	  {1,3,0,1, 1},
	  {1,3,1,0,-1},
	  {2,0,0,1, 1},
	  {2,0,1,0, 1},
	  {2,1,0,1,-1},
	  {2,1,1,0,-1},
	  {3,0,0,1, 1},
	  {3,0,1,0,-1},
	  {3,1,0,1,-1},
	  {3,1,1,0, 1},
	  {4,6,0,1, 1},
	  {4,6,1,0,-1},
	  {4,7,0,1,-1},
	  {4,7,1,0, 1},
	  {5,6,0,1, 1},
	  {5,6,1,0,-1},
	  {5,7,0,1,-1},
	  {5,7,1,0, 1},
	  {6,4,0,1,-1},
	  {6,4,1,0, 1},
	  {6,5,0,1, 1},
	  {6,5,1,0,-1},
	  {7,4,0,1,-1},
	  {7,4,1,0, 1},
	  {7,5,0,1, 1},
	  {7,5,1,0,-1} };

  int rcu[32][5] = { // right c_up
	  {0,4,1,0,-1},
	  {0,4,0,1, 1},
	  {0,5,1,0, 1},
	  {0,5,0,1, 1},
	  {1,4,1,0, 1},
	  {1,4,0,1, 1},
	  {1,5,1,0,-1},
	  {1,5,0,1, 1},
	  {2,6,1,0, 1},
	  {2,6,0,1, 1},
	  {2,7,1,0,-1},
	  {2,7,0,1, 1},
	  {3,6,1,0, 1},
	  {3,6,0,1, 1},
	  {3,7,1,0,-1},
	  {3,7,0,1, 1},
	  {4,0,1,0, 1},
	  {4,0,0,1,-1},
	  {4,1,1,0,-1},
	  {4,1,0,1,-1},
	  {5,0,1,0,-1},
	  {5,0,0,1,-1},
	  {5,1,1,0, 1},
	  {5,1,0,1,-1},
	  {6,2,1,0, 1},
	  {6,2,0,1, 1},
	  {6,3,1,0,-1},
	  {6,3,0,1,-1},
	  {7,2,1,0,-1},
	  {7,2,0,1, 1},
	  {7,3,1,0, 1},
	  {7,3,0,1,-1} };
 
  int rcd[32][5] = { // right c_down
	  {0,2,0,0, 1},
	  {0,2,1,1,-1},
	  {0,3,0,0, 1},
	  {0,3,1,1, 1},
	  {1,2,0,0, 1},
	  {1,2,1,1,-1},
	  {1,3,0,0, 1},
	  {1,3,1,1, 1},
	  {2,0,0,0, 1},
	  {2,0,1,1,-1},
	  {2,1,0,0,-1},
	  {2,1,1,1, 1},
	  {3,0,0,0,-1},
	  {3,0,1,1,-1},
	  {3,1,0,0, 1},
	  {3,1,1,1, 1},
	  {4,6,0,0,-1},
	  {4,6,1,1,-1},
	  {4,7,0,0,-1},
	  {4,7,1,1,-1},
	  {5,6,0,0,-1},
	  {5,6,1,1,-1},
	  {5,7,0,0,-1},
	  {5,7,1,1,-1},
	  {6,4,0,0,-1},
	  {6,4,1,1,-1},
	  {6,5,0,0, 1},
	  {6,5,1,1, 1},
	  {7,4,0,0, 1},
	  {7,4,1,1, 1},
	  {7,5,0,0,-1},
	  {7,5,1,1,-1} };
	  
  int lp[16][5] = { // left c_up c_down
	  {0,6,0,1, 1},
	  {0,7,0,1,-1},
	  {1,6,0,1, 1},
	  {1,7,0,1,-1},
	  {2,4,0,1,-1},
	  {2,5,0,1, 1},
	  {3,4,0,1,-1},
	  {3,5,0,1, 1},
	  {4,2,1,0, 1},
	  {4,3,1,0,-1},
	  {5,2,1,0, 1},
	  {5,3,1,0,-1},
	  {6,0,1,0,-1},
	  {6,1,1,0, 1},
	  {7,0,1,0,-1},
	  {7,1,1,0, 1} };

  int rp[16][5] = { // right c_up c_down
	  {0,6,0,1,-1},
	  {0,7,0,1,-1},
	  {1,6,0,1,-1},
	  {1,7,0,1,-1},
	  {2,4,1,0,-1},
	  {2,5,1,0, 1},
	  {3,4,1,0,-1},
	  {3,5,1,0, 1},
	  {4,2,0,1, 1},
	  {4,3,0,1,-1},
	  {5,2,0,1, 1},
	  {5,3,0,1,-1},
	  {6,0,1,0, 1},
	  {6,1,1,0,-1},
	  {7,0,1,0,-1},
	  {7,1,1,0, 1} };

  int pp[4][5] = { // cccc 
	  {0,0,0,0, 1},
	  {0,1,0,0,-1},
	  {1,0,0,0, 1},
	  {1,1,0,0,-1} };
	  
  int ln[8][5] = { // left n_up + n_down -1
	  {0,1,0,0,-1},
	  {1,0,0,0,-1},
	  {2,3,0,0,-1},
	  {3,2,0,0,-1},
	  {4,5,1,1,-1},
	  {5,4,1,1,-1},
	  {6,7,1,1,-1},
	  {7,6,1,1,-1} };

  int rn[8][5] = { // right n_up + n_down -1
	  {0,1,0,0,-1},
	  {1,0,0,0,-1},
	  {2,3,1,1,-1},
	  {3,2,1,1,-1},
	  {4,5,0,0,-1},
	  {5,4,0,0,-1},
	  {6,7,1,1, 1},
	  {7,6,1,1, 1} };
}
	  
using namespace gg;

qtensor<complex<double> > H_4e1(double t1, double t3, double f0,
                                double f1, double p1, double frac)
{
  complex<double> i = complex<double> (0,1);
  
  vector< vector< tensor<complex<double> > > > t;
  t.resize(8);
  for(int k=0;k<8;k++)
  {
    t[k].resize(8);
    for(int j=0;j<8;j++) t[k][j] = tensor<complex<double> >(9,2,9,2);
  }
  
  // I: (0,0) (8,8)
  for(int k=0;k<8;k++) 
  {
    t[k][k].update( 1, 0,0,0,0);
    t[k][k].update( 1, 0,1,0,1);
    t[k][k].update( 1, 8,0,8,0);
    t[k][k].update( 1, 8,1,8,1);
  }
  // 2i t1 c^d c - f0(2n-1)(2n-1): (8,0)
  t[0][0].update(-2*f0, 8,0,0,0);
  t[0][0].update( 2*f0, 8,1,0,1);
  t[1][1].update(-2*f0, 8,0,0,0);
  t[1][1].update( 2*f0, 8,1,0,1);
  t[6][6].update( 2*f0, 8,0,0,0);
  t[6][6].update(-2*f0, 8,1,0,1);
  t[7][7].update( 2*f0, 8,0,0,0);
  t[7][7].update(-2*f0, 8,1,0,1);
  t[2][2].update(-2.*i*t1, 8,0,0,1);
  t[2][2].update( 2.*i*t1, 8,1,0,0);
  t[3][3].update(-2.*i*t1, 8,0,0,1);
  t[3][3].update( 2.*i*t1, 8,1,0,0);
  t[4][4].update( 2.*i*t1, 8,0,0,1);
  t[4][4].update(-2.*i*t1, 8,1,0,0);
  t[5][5].update( 2.*i*t1, 8,0,0,1);
  t[5][5].update(-2.*i*t1, 8,1,0,0);
  t[6][6].update( 4.*i*t1, 8,0,0,1);
  t[6][6].update(-4.*i*t1, 8,1,0,0); 
  
  for(int k=0;k<32;k++) 
  {
    // c_u: left (1,0) right (8,2)
    t[ lcu[k][0] ][ lcu[k][1] ].update(0.5*lcu[k][4], 1, lcu[k][2], 0, lcu[k][3]);
    t[ rcu[k][0] ][ rcu[k][1] ].update(i*(t3*rcu[k][4]), 8, rcu[k][2], 2, rcu[k][3]);
    // c_u^d: left (2,0) right (8,1)
    t[ lcu[k][1] ][ lcu[k][0] ].update(0.5*lcu[k][4], 2, lcu[k][3], 0, lcu[k][2]);
    t[ rcu[k][1] ][ rcu[k][0] ].update(i*(t3*rcu[k][4]), 8, rcu[k][3], 1, rcu[k][2]);
  
    // c_d: left (3,0) right (8,4)
    t[ lcd[k][0] ][ lcd[k][1] ].update(0.5*lcd[k][4], 3, lcd[k][2], 0, lcd[k][3]);
    t[ rcd[k][0] ][ rcd[k][1] ].update(i*(t3*rcd[k][4]), 8, rcd[k][2], 4, rcd[k][3]);
    // c_d^d: left (4,0) right (8,3)
    t[ lcd[k][1] ][ lcd[k][0] ].update(0.5*lcd[k][4], 4, lcd[k][3], 0, lcd[k][2]);
    t[ rcd[k][1] ][ rcd[k][0] ].update(i*(t3*rcd[k][4]), 8, rcd[k][3], 3, rcd[k][2]);
  }
  
  complex<double> pc = p1*complex<double>( cos(2*PI*frac) , sin(2*PI*frac) );
  complex<double> pd = p1*complex<double>( cos(2*PI*frac) , -sin(2*PI*frac) );
  for(int k=0;k<16;k++) 
  {
    // c_u c_d: left (5,0) right (8,5)
    t[ lp[k][0] ][ lp[k][1] ].update(0.5*lp[k][4], 5, lp[k][2], 0, lp[k][3]);
    t[ rp[k][0] ][ rp[k][1] ].update(pc*(4.0*rp[k][4]), 8, rp[k][2], 5, rp[k][3]);
    // c_d^d: left (6,0) right (8,6)
    t[ lp[k][1] ][ lp[k][0] ].update(0.5*lp[k][4], 6, lp[k][3], 0, lp[k][2]);
    t[ rp[k][1] ][ rp[k][0] ].update(pd*(4.0*rp[k][4]), 8, rp[k][3], 6, rp[k][2]);
  }
  for(int k=0;k<8;k++) 
  {
    // n_u + n_d -1: left (7,0) right (8,7)
    t[ ln[k][0] ][ ln[k][1] ].update(ln[k][4], 7, ln[k][2], 0, ln[k][3]);
    t[ rn[k][0] ][ rn[k][1] ].update(-4*f1*rn[k][4], 8, rn[k][2], 7, rn[k][3]);
  }
  
  qtensor<complex<double> > ret(1,8,1,8);
  ret.u_dim(9,2,9,2);
  for(int k=0;k<8;k++) for(int j=0;j<8;j++) 
    ret.update(t[j][k], 0,j,0,k);
  ret.add_sign(0, -1, 0, 1);
  ret = ret.simplify();
  return ret;
}


qtensor<complex<double> > H_4e2(double t1, double t3, double f0,
                                double f1, double p1, double frac)
{
  complex<double> i = complex<double> (0,1);
  
  vector< vector< tensor<complex<double> > > > t;
  t.resize(8);
  for(int k=0;k<8;k++)
  {
    t[k].resize(8);
    for(int j=0;j<8;j++) t[k][j] = tensor<complex<double> >(9,2,9,2);
  }
  
  // I: (0,0) (8,8)
  for(int k=0;k<8;k++) 
  {
    t[k][k].update( 1, 0,0,0,0);
    t[k][k].update( 1, 0,1,0,1);
    t[k][k].update( 1, 8,0,8,0);
    t[k][k].update( 1, 8,1,8,1);
  }
  // 2i t1 c^d c - f0(2n-1)(2n-1): (8,0)
  t[0][0].update(-2*f0, 8,0,0,0);
  t[0][0].update( 2*f0, 8,1,0,1);
  t[1][1].update(-2*f0, 8,0,0,0);
  t[1][1].update( 2*f0, 8,1,0,1);
  t[6][6].update( 2*f0, 8,0,0,0);
  t[6][6].update(-2*f0, 8,1,0,1);
  t[7][7].update( 2*f0, 8,0,0,0);
  t[7][7].update(-2*f0, 8,1,0,1);
  t[2][2].update(-2.*i*t1, 8,0,0,1);
  t[2][2].update( 2.*i*t1, 8,1,0,0);
  t[3][3].update(-2.*i*t1, 8,0,0,1);
  t[3][3].update( 2.*i*t1, 8,1,0,0);
  t[4][4].update( 2.*i*t1, 8,0,0,1);
  t[4][4].update(-2.*i*t1, 8,1,0,0);
  t[5][5].update( 2.*i*t1, 8,0,0,1);
  t[5][5].update(-2.*i*t1, 8,1,0,0);
  t[6][6].update( 4.*i*t1, 8,0,0,1);
  t[6][6].update(-4.*i*t1, 8,1,0,0); 
 
  for(int k=0;k<32;k++) 
  {
    // c_u: right (1,0) left (8,2)
    t[ rcu[k][0] ][ rcu[k][1] ].update(0.5*rcu[k][4], 1, rcu[k][2], 0, rcu[k][3]);
    t[ lcu[k][0] ][ lcu[k][1] ].update(i*(t3*lcu[k][4]), 8, lcu[k][2], 2, lcu[k][3]);
    // c_u^d: right (2,0) left (8,1)
    t[ rcu[k][1] ][ rcu[k][0] ].update(0.5*rcu[k][4], 2, rcu[k][3], 0, rcu[k][2]);
    t[ lcu[k][1] ][ lcu[k][0] ].update(i*(t3*lcu[k][4]), 8, lcu[k][3], 1, lcu[k][2]);
  
    // c_d: right (3,0) left (8,4)
    t[ rcd[k][0] ][ rcd[k][1] ].update(0.5*rcd[k][4], 3, rcd[k][2], 0, rcd[k][3]);
    t[ lcd[k][0] ][ lcd[k][1] ].update(i*(t3*lcd[k][4]), 8, lcd[k][2], 4, lcd[k][3]);
    // c_d^d: right (4,0) left (8,3)
    t[ rcd[k][1] ][ rcd[k][0] ].update(0.5*rcd[k][4], 4, rcd[k][3], 0, rcd[k][2]);
    t[ lcd[k][1] ][ lcd[k][0] ].update(i*(t3*lcd[k][4]), 8, lcd[k][3], 3, lcd[k][2]);
  }
  
  complex<double> pc = p1*complex<double>( cos(2*PI*frac) , sin(2*PI*frac) );
  complex<double> pd = p1*complex<double>( cos(2*PI*frac) , -sin(2*PI*frac) );
  for(int k=0;k<16;k++) 
  {
    // c_u c_d: left (5,0) right (8,5)
    t[ lp[k][0] ][ lp[k][1] ].update(0.5*lp[k][4], 5, lp[k][2], 0, lp[k][3]);
    t[ rp[k][0] ][ rp[k][1] ].update(pc*(4.0*rp[k][4]), 8, rp[k][2], 5, rp[k][3]);
    // c_d^d: left (6,0) right (8,6)
    t[ lp[k][1] ][ lp[k][0] ].update(0.5*lp[k][4], 6, lp[k][3], 0, lp[k][2]);
    t[ rp[k][1] ][ rp[k][0] ].update(pd*(4.0*rp[k][4]), 8, rp[k][3], 6, rp[k][2]);
  }
  for(int k=0;k<8;k++) 
  {
    // n_u + n_d -1: left (7,0) right (8,7)
    t[ ln[k][0] ][ ln[k][1] ].update(ln[k][4], 7, ln[k][2], 0, ln[k][3]);
    t[ rn[k][0] ][ rn[k][1] ].update(-4*f1*rn[k][4], 8, rn[k][2], 7, rn[k][3]);
  }
  
  qtensor<complex<double> > ret(1,8,1,8);
  ret.u_dim(9,2,9,2);
  for(int k=0;k<8;k++) for(int j=0;j<8;j++) 
    ret.update(t[j][k], 0,j,0,k);
  ret.add_sign(0, -1, 0, 1);
  ret = ret.simplify();
  return ret;
}


qtensor<complex<double> > c_up1()
{
  vector< vector< tensor<complex<double> > > > t;
  t.resize(8);
  for(int k=0;k<8;k++)
  {
    t[k].resize(8);
    for(int j=0;j<8;j++) t[k][j] = tensor<complex<double> >(1,2,1,2);
  }

  for(int k=0;k<32;k++) 
    t[ lcu[k][0] ][ lcu[k][1] ].update(0.5*lcu[k][4], 0, lcu[k][2], 0, lcu[k][3]);

  qtensor<complex<double> > ret(1,8,1,8);
  ret.u_dim(1,2,1,2);
  for(int k=0;k<8;k++) for(int j=0;j<8;j++) 
    ret.update(t[j][k], 0,j,0,k);
  ret.add_sign(0, -1, 0, 1);
  ret = ret.simplify();
  return ret;
}

qtensor<complex<double> > c_down1()
{
  vector< vector< tensor<complex<double> > > > t;
  t.resize(8);
  for(int k=0;k<8;k++)
  {
    t[k].resize(8);
    for(int j=0;j<8;j++) t[k][j] = tensor<complex<double> >(1,2,1,2);
  }

  for(int k=0;k<32;k++) 
    t[ lcd[k][0] ][ lcd[k][1] ].update(0.5*lcd[k][4], 0, lcd[k][2], 0, lcd[k][3]);

  qtensor<complex<double> > ret(1,8,1,8);
  ret.u_dim(1,2,1,2);
  for(int k=0;k<8;k++) for(int j=0;j<8;j++) 
    ret.update(t[j][k], 0,j,0,k);
  ret.add_sign(0, -1, 0, 1);
  ret = ret.simplify();
  return ret;
}

qtensor<complex<double> > c_up2()
{
  vector< vector< tensor<complex<double> > > > t;
  t.resize(8);
  for(int k=0;k<8;k++)
  {
    t[k].resize(8);
    for(int j=0;j<8;j++) t[k][j] = tensor<complex<double> >(1,2,1,2);
  }

  for(int k=0;k<32;k++) 
    t[ rcu[k][0] ][ rcu[k][1] ].update(0.5*rcu[k][4], 0, rcu[k][2], 0, rcu[k][3]);

  qtensor<complex<double> > ret(1,8,1,8);
  ret.u_dim(1,2,1,2);
  for(int k=0;k<8;k++) for(int j=0;j<8;j++) 
    ret.update(t[j][k], 0,j,0,k);
  ret.add_sign(0, -1, 0, 1);
  ret = ret.simplify();
  return ret;
}

qtensor<complex<double> > c_down2()
{
  vector< vector< tensor<complex<double> > > > t;
  t.resize(8);
  for(int k=0;k<8;k++)
  {
    t[k].resize(8);
    for(int j=0;j<8;j++) t[k][j] = tensor<complex<double> >(1,2,1,2);
  }

  for(int k=0;k<32;k++) 
    t[ rcd[k][0] ][ rcd[k][1] ].update(0.5*rcd[k][4], 0, rcd[k][2], 0, rcd[k][3]);

  qtensor<complex<double> > ret(1,8,1,8);
  ret.u_dim(1,2,1,2);
  for(int k=0;k<8;k++) for(int j=0;j<8;j++) 
    ret.update(t[j][k], 0,j,0,k);
  ret.add_sign(0, -1, 0, 1);
  ret = ret.simplify();
  return ret;
}


qtensor<complex<double> > c_pair1()
{
  vector< vector< tensor<complex<double> > > > t;
  t.resize(8);
  for(int k=0;k<8;k++)
  {
    t[k].resize(8);
    for(int j=0;j<8;j++) t[k][j] = tensor<complex<double> >(1,2,1,2);
  }

  for(int k=0;k<16;k++) 
    t[ lp[k][0] ][ lp[k][1] ].update(0.5*lp[k][4], 0, lp[k][2], 0, lp[k][3]);

  qtensor<complex<double> > ret(1,8,1,8);
  ret.u_dim(1,2,1,2);
  for(int k=0;k<8;k++) for(int j=0;j<8;j++) 
    ret.update(t[j][k], 0,j,0,k);
  ret.add_sign(0, -1, 0, 1);
  ret = ret.simplify();
  return ret;
}


qtensor<complex<double> > c_pair2()
{
  vector< vector< tensor<complex<double> > > > t;
  t.resize(8);
  for(int k=0;k<8;k++)
  {
    t[k].resize(8);
    for(int j=0;j<8;j++) t[k][j] = tensor<complex<double> >(1,2,1,2);
  }

  for(int k=0;k<16;k++) 
    t[ rp[k][0] ][ rp[k][1] ].update(0.5*rp[k][4], 0, rp[k][2], 0, rp[k][3]);

  qtensor<complex<double> > ret(1,8,1,8);
  ret.u_dim(1,2,1,2);
  for(int k=0;k<8;k++) for(int j=0;j<8;j++) 
    ret.update(t[j][k], 0,j,0,k);
  ret.add_sign(0, -1, 0, 1);
  ret = ret.simplify();
  return ret;
}


qtensor<complex<double> > c_pairpair()
{
  vector< vector< tensor<complex<double> > > > t;
  t.resize(8);
  for(int k=0;k<8;k++)
  {
    t[k].resize(8);
    for(int j=0;j<8;j++) t[k][j] = tensor<complex<double> >(1,2,1,2);
  }

  for(int k=0;k<4;k++) 
    t[ pp[k][0] ][ pp[k][1] ].update(0.5*pp[k][4], 0, pp[k][2], 0, pp[k][3]);

  qtensor<complex<double> > ret(1,8,1,8);
  ret.u_dim(1,2,1,2);
  for(int k=0;k<8;k++) for(int j=0;j<8;j++) 
    ret.update(t[j][k], 0,j,0,k);
  ret.add_sign(0, -1, 0, 1);
  ret = ret.simplify();
  return ret;
}


qtensor<complex<double> > c_pair_dagger1()
{
  vector< vector< tensor<complex<double> > > > t;
  t.resize(8);
  for(int k=0;k<8;k++)
  {
    t[k].resize(8);
    for(int j=0;j<8;j++) t[k][j] = tensor<complex<double> >(1,2,1,2);
  }

  for(int k=0;k<16;k++) 
    t[ lp[k][1] ][ lp[k][0] ].update(-0.5*lp[k][4], 0, lp[k][3], 0, lp[k][2]);

  qtensor<complex<double> > ret(1,8,1,8);
  ret.u_dim(1,2,1,2);
  for(int k=0;k<8;k++) for(int j=0;j<8;j++) 
    ret.update(t[j][k], 0,j,0,k);
  ret.add_sign(0, -1, 0, 1);
  ret = ret.simplify();
  return ret;
}


qtensor<complex<double> > c_pair_dagger2()
{
  vector< vector< tensor<complex<double> > > > t;
  t.resize(8);
  for(int k=0;k<8;k++)
  {
    t[k].resize(8);
    for(int j=0;j<8;j++) t[k][j] = tensor<complex<double> >(1,2,1,2);
  }

  for(int k=0;k<16;k++) 
    t[ rp[k][1] ][ rp[k][0] ].update(-0.5*rp[k][4], 0, rp[k][3], 0, rp[k][2]);

  qtensor<complex<double> > ret(1,8,1,8);
  ret.u_dim(1,2,1,2);
  for(int k=0;k<8;k++) for(int j=0;j<8;j++) 
    ret.update(t[j][k], 0,j,0,k);
  ret.add_sign(0, -1, 0, 1);
  ret = ret.simplify();
  return ret;
}
