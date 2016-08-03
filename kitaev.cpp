/* This is the main function to calculate the Kitaev Chain. */

#include "global.h"
#include "functions.h"
#include "operators.h"

int main(int argc, char *argv[])
{
  if(SYMMETRY-2)
  {
    cerr << "Define the SYMMETRY in global.h as 2." << endl;
    return 0;
  }
#ifndef FERMION
  cerr << "Define FERMION in global.h" << endl;
  return 0;
#endif

  time_t start, end;
  time(&start);
  
  // ********** initialize the system **********

  vector<string> para_name; // name of the parameters
  para_name = set_para_name("t","p","u");

  int sites = 4;
  int cutoff = 10;
  int sweep = -1; 
  int symmetry_sector = 0;
  
  vector<double> para;
  string filename = "";
  para = set_para_val(argc, argv, sites, cutoff, sweep, symmetry_sector,
                      para_name, filename);

  mpo<double> my_mpo(sites);
  for(int i=0;i<sites;i++) my_mpo[i] = H_Kitaev(para[0], para[1], para[2]);
  my_mpo.edge(H_Kitaev(para[0], para[1], para[2]).left(),
              H_Kitaev(para[0], para[1], para[2]).right());
  mps<double> my_mps(my_mpo);

#ifdef PBC
  cout << "Using periodic boundary condition\n" << endl;
  my_mps.add_edge(H_Kitaev(para[0], -para[1], 0.0).left(),
                  H_Kitaev(para[0], -para[1], 0.0).right());
#else
  cout << "Using open boundary condition\n" << endl;
#endif

  cout.precision(10);
  dmrg(my_mps, my_mpo, cutoff, sweep, symmetry_sector, filename);
  my_mps.prep_calc();
  cout << "Normalization: " << calc(my_mps) << endl;

  time(&end);
  cout << "Time used: " << difftime(end,start) << endl;
}

