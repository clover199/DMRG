
#include "global.h"
#include "functions.h"
#include "operators.h"

#define FERMION
//#define SYMMETRY 3

extern string filename = "";
extern bool print = false;
extern int symmetry_sector = -1;

int main(int argc, char *argv[])
{
  time_t start, end;
  time(&start);
  
  // ********** initialize the system **********

  vector<string> para_name; // name of the parameters
  para_name = set_para_name("j","h");

  int sites = 4;
  int cutoff = 10;
  int sweep = -1;
  
  /* We assume the parameters are set by 'argv' in the form of (for example):
    L=4  m=10  n=-1  t=1  delta=1  u=1 ...
    The order doesn't matter. If the value is not set, we use default values. */
  vector<double> para( para_name.size(), 1.0 );
  para = set_para_val(argc, argv, sites, cutoff, sweep, symmetry_sector, para_name);

  mpo<double> my_mpo(sites, H_Ising(para[0], para[1]));
  my_mpo.edge(H_Ising_ledge(para[0], para[1]), H_Ising_redge(para[0], para[1]));
  mps<double> my_mps(my_mpo);

  dmrg(my_mps, my_mpo, cutoff, sweep);
  
  time(&end);
  cout << "Time used: " << difftime(end,start) << endl;
}
