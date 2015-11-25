
#include "global.h"
#include "functions.h"
#include "operators.h"
#include "dmrg.h"
#include "tensor.h"

extern string filename = "";
extern bool print = false;


int main(int argc, char *argv[])
{
  bool write = false;
  time_t start, end;
  time(&start);
  
  // ********** initialize the system **********

  vector<string> para_name; // name of the parameters
  para_name = set_para_name("t","delta","u");

  int sites = 4;
  int cutoff = 10;
  int sweep = -1;
  
  // We assume the parameters are set by 'argv' in the form of (for example):
  // L=4  m=10  n=-1  t=1  delta=1  u=1 ...
  // The order doesn't matter. If the value is not set, we use the defalut values.
  vector<double> para( para_name.size(), 1.0 );
  para = set_para_val(argc, argv, sites, cutoff, sweep, para_name);
  
  mps my_mps(sites);
  tensor ham;
  ham = H_spinless(para[0], para[1], para[2]);
  mpo my_mpo(sites, ham);

  dmrg(my_mps, my_mpo, cutoff, sweep);
 
  time(&end);
  cout << "Time used: " << difftime(end,start) << endl;
}
