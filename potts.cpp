/* This is the main function for the 3-Potts model. */


#include "global.h"
#include "functions.h"
#include "operators.h"

int main(int argc, char *argv[])
{
  if(SYMMETRY-3)
  {
    cerr << "Define the SYMMETRY in global.h as 3." << endl;
    return 0;
  }
#ifdef FERMION
  cerr << "Do not define FERMION in global.h" << endl;
  return 0;
#endif

  time_t start, end;
  time(&start);
  
  // ********** initialize the system **********

  vector<string> para_name; // name of the parameters
  para_name = set_para_name("f","j","p","t");

  int sites = 4;
  int cutoff = 10;
  int sweep = -1;
  int symmetry_sector = 0;
  
  /* We assume the parameters are set by 'argv' in the form of (for example):
    L=4  m=10  n=-1  t=1  delta=1  u=1 ...
    The order doesn't matter. If the value is not set, we use default values. */

  vector<double> para;
  string filename = "";
  para = set_para_val(argc, argv, sites, cutoff, sweep, symmetry_sector,
                      para_name, filename);

  mpo< complex<double> > my_mpo(sites, H_Potts3(para[0], para[1], para[2], para[3]));
  my_mpo.edge();
  mps< complex<double> > my_mps(my_mpo);
#ifdef PBC
  cout << "Using periodic boundary condition\n" << endl;
  my_mps.add_edge(H_Potts3(0, para[1], para[2], para[3]).left(),
                  H_Potts3(0, para[1], para[2], para[3]).right());
#else
  cout << "Using open boundary condition\n" << endl;
#endif

  dmrg(my_mps, my_mpo, cutoff, sweep, symmetry_sector, filename);
  
  time(&end);
  cout << "Time used: " << difftime(end,start) << endl;
}
