
#include "global.h"
#include "mps.h"
#include "mpo.h"
#include "useful.h"
#include "basic.h"

// The notation of the index for the MPS and MPO
//
// U,V:  0-O-2
//         |
//         1
//
//  --O-2     3     2-O--
//    |       |       |
//  --O-1   0-O-2   1-O--
//    |       |       |
//  --O-0     1     0-O--
//    l               r

namespace{
  ofstream data_energy;
  ofstream data_singular;
  qtensor<double> lenv;  // the left environment
  qtensor<double> renv;  // the right environment
  double * store;
  vector< vector<int> > lmap;
  vector< vector<int> > rmap;
}


void av2 (int n, double *in, double *out)
{
  renv.contract(store, in, rmap, 3, 'N', 'T');
  lenv.contract(out, store, lmap, 2, 'N', 'T');
}


void update_two(int l, int r, int cutoff, mps<double>& my_mps, mpo<double>& my_mpo)
{
  int dim = 1;
  dim *= my_mps(l).dimension(2);  // left environment
  dim *= my_mpo[l+1].dimension(3);  // the left point
  dim *= my_mpo[r-1].dimension(3);  // the right point
  dim *= my_mps(r).dimension(2);  // right environment
  
  lenv.contract(my_mps(l), 1, my_mpo[l+1], 0);
  qtensor<double> ope = my_mpo[r-1].exchange(0,2);
  renv.contract(my_mps(r), 1, ope, 0);
  
  vector< vector<int> > fullmap;
  vector<int> index;
  index.resize(4);
  index[0] = my_mps(l).index(2);
  index[1] = my_mpo[l+1].index(3);
  index[2] = my_mpo[r-1].index(3);
  index[3] = my_mps(r).index(2);
  generate_map(fullmap, index);
    
  double *val = new double [NEV];
  double *vecs = new double [dim];
  store = new double [dim*my_mpo[r-1].dimension(2)];
  dsaupd(dim, NEV, val, vecs, av2);
  delete store;
  
  qtensor<double> vec(vecs, fullmap, fullmap);
  qtensor<double> U, V, S;
  vec.svd(U, S, V, 2, cutoff);
  my_mps.update_state(l, U);
  my_mps.update(S);
  my_mps.update_state(r, V);
  
  ope.contract(my_mpo[l], U, 'N', 'N');
  my_mps(l+1).contract(U, ope, 'T', 'N');
  ope.contract(my_mpo[r], V, 'N', 'T');
  my_mps(r-1).contract(V, ope, 'N', 'N');
}


void move2right(int l, int cutoff, mps<double>& my_mps, mpo<double>& my_mpo)
{

}


void move2left(int r, int cutoff, mps<double>& my_mps, mpo<double>& my_mpo)
{

}


void dmrg(mps<double> my_mps, mpo<double> my_mpo, int cutoff, int sweep)
{
  int L = my_mps.size();
  
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
    move2right(i, pre_cutoff, my_mps, my_mpo);

  data_energy.close();
  data_singular.close();
}

