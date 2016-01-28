/* This is the main function for the spinless fermion

Make sure to define FERMION and SYMMETRY is defined as 2

*/

#include "global.h"
#include "functions.h"
#include "operators.h"

int main(int argc, char *argv[])
{
  time_t start, end;
  time(&start);
  
  // ********** initialize the system **********

  vector<string> para_name; // name of the parameters
  para_name = set_para_name("t","p","u");

  int sites = 4;
  int cutoff = 10;
  int sweep = -1; 
  int symmetry_sector = 0;
  
  vector<double> para( para_name.size(), 1.0 );
  string filename = "";
  para = set_para_val(argc, argv, sites, cutoff, sweep, symmetry_sector,
                      para_name, filename);

  mpo<double> my_mpo(sites, H_spinless_fermion(para[0], para[1], para[2]));
  my_mpo.edge(H_spinless_fermion_ledge(para[0], para[1], para[2]),
              H_spinless_fermion_redge(para[0], para[1], para[2]));
  mps<double> my_mps(my_mpo);

  cout.precision(10);
  dmrg(my_mps, my_mpo, cutoff, sweep, symmetry_sector, filename);
  
  time(&end);
  cout << "Time used: " << difftime(end,start) << endl;
}
