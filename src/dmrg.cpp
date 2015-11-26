
#include "global.h"
#include "mps.h"
#include "mpo.h"

namespace{
  vector<tensor> store_mps; // used to store_mps the multiplied MPS
  ofstream data_energy;
  ofstream data_singular;
}


template <typename T>
void av2 (int n, T *in, T *out)
{

}


void update_two(int l, int r, int cutoff, mps& my_mps, mpo& my_mpo)
{
  
//  znaupd(dim, NEV, val, vec, av2);

  tensor vec;
  tensor U, V;
  vector<double> S;
  vec.svd(U, S, V);

}


void update_one(int l, int cutoff, mps& my_mps, mpo& my_mpo)
{
//  znaupd(dim, NEV, val, vec, av1);
}


void dmrg(mps my_mps, mpo my_mpo, int cutoff, int sweep)
{
  int L;
  L = my_mps.initialize(my_mpo);
  my_mps.create_store(store_mps);
  
  string name = "energy"+filename;
  data_energy.open(name.c_str());
  name = "singular"+filename;
  data_singular.open(name.c_str());

  int pre_cutoff = 10;
  int pre_sweep = cutoff/50;
  
  print = false;
  for(int i=0;i<L/2;i++)
    update_two(i, L-i-1, pre_cutoff, my_mps, my_mpo);
  for(int i=L/2;i<L;i++)
    update_one(i, pre_cutoff, my_mps, my_mpo);
/*
  for(int t=0;t<pre_sweep;t++)
  {
    pre_cutoff = 50*t;
	for(int i=L;i>0;i--)
      update_one();
	for(int i=0;i<L;i++)
      update_one();
  }
  
  for(int t=0;t<sweep;t++)
  {
	for(int i=L;i>0;i--)
      update_one();
	for(int i=0;i<L;i++)
      update_one();
  }
*/
  data_energy.close();
  data_singular.close();
}

