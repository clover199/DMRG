
#include "global.h"
#include "functions.h"
#include "operators.h"

int main(int argc, char *argv[])
{
  time_t start, end;
  time(&start);
  
  // ********** initialize the system **********

  vector<string> para_name; // name of the parameters
  para_name = set_para_name("j","h");

  int sites = 6;
  int cutoff = 10;
  int sweep = -1; 
  int symmetry_sector = 0;
  
  vector<double> para( para_name.size(), 1.0 );
  string filename = "";
  para = set_para_val(argc, argv, sites, cutoff, sweep, symmetry_sector,
                      para_name, filename);

  mpo<double> my_mpo(sites, H_2Potts(para[0], para[1]));
  my_mpo.edge(H_2Potts_ledge(para[0], para[1]),
              H_2Potts_redge(para[0], para[1]));
  mps<double> my_mps(my_mpo);

  cout.precision(10);
  dmrg(my_mps, my_mpo, cutoff, sweep, symmetry_sector, filename);
  
  time(&end);
  cout << "Time used: " << difftime(end,start) << endl;
}
