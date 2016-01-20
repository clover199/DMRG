
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
  int r_num;
  int l_num;
}


void av2 (int n, double *in, double *out)
{
  renv.contract(store, in, rmap, 'N', 'T', r_num);
  lenv.contract(out, store, lmap, 'N', 'T', l_num);
}


void check_hermitian(int n)
{
  double* in = new double [n];
  double* out = new double [n];
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
    {
      double temp;
      for(int s=0;s<n;s++) in[s] = 0;
      in[j] = 1;
      av2(n, in, out);
      temp = out[i];
//      cout <<  out[i] << "\t";
      for(int s=0;s<n;s++) in[s] = 0;
      in[i] = 1;
      av2(n, in, out);
      if(abs(temp-out[j])>TOL) cout << "Not hermitian!\n";
    }
//    cout << endl;
  }
}


void two_sites(int size, int cutoff, mps<double>& my_mps, mpo<double>& my_mpo)
{
  int l = 0;
  int r = size-1;
  lenv = my_mpo[l];
  renv = my_mpo[r];
  int d = 1;
  d *= lenv.dimension(2);  // the left point
  d *= renv.dimension(2);  // the right point

  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, lenv, 1, renv, 1);
  double *val;
  double *vecs = new double [d];
  if(d<NEV*2)
  {
    qtensor<double> tt;
    tt.contract(lenv, 1, renv, 1);
    tt = tt.exchange(2,3);
    tt = tt.combine(2,3);
    tt = tt.combine(0,1);
    val = new double [d];
    tt.eig(val, vecs);
    cout << "Energy:\n";
    for(int i=0;i<d;i++) cout << "  " << val[i] << endl;
    data_energy << l << "\t" << r;
    if(NEV<d) for(int i=0;i<NEV;i++) data_energy << "\t" << val[i];
    else
    {
      for(int i=0;i<d;i++) data_energy << "\t" << val[i];
      for(int i=d;i<NEV;i++) data_energy << "\t" << 0;
    }
    data_energy << endl;
  }
  else
  {
    r_num = 1;
    l_num = 2;
    renv.get_map(rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', r_num, 0);
    lenv.get_map(lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', l_num, 0);

    val = new double [NEV];
    int d_store = 1;
    for(int i=0;i<dim_ret.size();i++)
    {
      int temp = 0;
      for(int j=0;j<dim_ret[i].size();j++) temp += dim_ret[i][j];
      d_store *= temp;
    }
    store = new double [d_store];
    check_hermitian(d);
    dsaupd(d, NEV, val, vecs, av2);
    delete store;
    cout << "Energy:\n";
    for(int i=0;i<NEV;i++) cout << "  " << val[i] << endl;
    data_energy << l << "\t" << r;
    for(int i=0;i<NEV;i++) data_energy << "\t" << val[i];
    data_energy << endl;
  }

  qtensor<double> vec(vecs, dim, sym);
  qtensor<double> U, V, S;
  vector<double> s;
  s = vec.svd(U, S, V, 1, cutoff);
  double entropy = 0;
  for(int i=0;i<s.size();i++) entropy += -s[i]*s[i]*log(s[i]*s[i]);
  cout << "Singular value:\n";
  for(int i=0;i<s.size();i++) cout << "  " << s[i] << endl;
  cout << "Entropy: " << entropy << endl;
  data_singular << l << "\t" << r << "\t" << entropy;
  int max = 20;
  if(s.size()>max) for(int i=0;i<max;i++) data_singular << "\t" << s[i];
  else
  {
    for(int i=0;i<s.size();i++) data_singular << "\t" << s[i];
    for(int i=s.size();i<max;i++) data_singular << "\t" << 0;
  }
  data_singular << endl;

  my_mps[l] = U;
  my_mps.center() = S;
  my_mps[r] = V;
  vec.contract(lenv, U, 'N', 'N');
  my_mps(l).contract(U, vec, 'T', 'N');
  vec.contract(renv, V, 'N', 'T');
  my_mps(r).contract(V, vec, 'N', 'N');
}

void update_two(int l, int r, int cutoff, mps<double>& my_mps, mpo<double>& my_mpo)
{
  int d = 1;
  d *= my_mps(l-1).dimension(2);  // left environment
  d *= my_mpo[l].dimension(3);  // the left point
  d *= my_mpo[r].dimension(3);  // the right point
  d *= my_mps(r+1).dimension(2);  // right environment

  lenv.contract(my_mps(l-1), 1, my_mpo[l], 0);
  lenv = lenv.exchange(3,4);
  qtensor<double> ope = my_mpo[r].exchange(0,2);
  renv.contract(my_mps(r+1), 1, ope, 0);
  renv = renv.exchange(0,1);
  
  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, lenv, 2, renv, 2);
  r_num = 2;
  l_num = 3;
  renv.get_map(rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', r_num, 0);
  lenv.get_map(lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', l_num, 0);

  double *val = new double [NEV];
  double *vecs = new double [d];
  int d_store = 1;
  for(int i=0;i<dim_ret.size();i++)
  {
    int temp = 0;
    for(int j=0;j<dim_ret[i].size();j++) temp += dim_ret[i][j];
    d_store *= temp;
  }
  store = new double [d_store];
//  check_hermitian(d);
  cout << "Energy:\n";
  dsaupd(d, NEV, val, vecs, av2);
  delete store;
  for(int i=0;i<NEV;i++) cout << "  " << val[i] << endl;
  data_energy << l << "\t" << r;
  for(int i=0;i<NEV;i++) data_energy << "\t" << val[i];
  data_energy << endl;

  qtensor<double> vec(vecs, dim, sym);
  qtensor<double> U, V, S;
  vector<double> s;
  s = vec.svd(U, S, V, 2, cutoff);
  double entropy = 0;
  for(int i=0;i<s.size();i++) entropy += -s[i]*s[i]*log(s[i]*s[i]);
  cout << "Singular value:\n";
  for(int i=0;i<s.size();i++) cout << "  " << s[i] << endl;
  cout << "Entropy: " << entropy << endl;
  data_singular << l << "\t" << r << "\t" << entropy;
  int max = 20;
  if(s.size()>max) for(int i=0;i<max;i++) data_singular << "\t" << s[i];
  else
  {
    for(int i=0;i<s.size();i++) data_singular << "\t" << s[i];
    for(int i=s.size();i<max;i++) data_singular << "\t" << 0;
  }
  data_singular << endl;

  my_mps[l] = U;
  my_mps.center() = S;
  my_mps[r] = V;
  vec.contract(lenv, U, 'N', 'N', 2);
  my_mps(l).contract(U, vec, 'T', 'N', 2);
  vec.contract(renv, V, 'N', 'T', 2);
  my_mps(r).contract(V, vec, 'N', 'N', 2);
}


void move2right(int l, int cutoff, mps<double>& my_mps, mpo<double>& my_mpo)
{
  int d = 1;
  d *= my_mps(l-1).dimension(2);  // left environment
  d *= my_mpo[l].dimension(3);  // the left point
  d *= my_mps(l+1).dimension(2);  // right environment

  lenv.contract(my_mps(l-1), 1, my_mpo[l], 0);
  lenv = lenv.exchange(3,4);
  renv = my_mps(l+1);
  
  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, lenv, 2, renv, 1);
  r_num = 1;
  l_num = 3;
  renv.get_map(rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', r_num, 0);
  lenv.get_map(lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', l_num, 0);

  double *val = new double [NEV];
  double *vecs = new double [d];
  int d_store = 1;
  for(int i=0;i<dim_ret.size();i++)
  {
    int temp = 0;
    for(int j=0;j<dim_ret[i].size();j++) temp += dim_ret[i][j];
    d_store *= temp;
  }
  store = new double [d_store];
//  check_hermitian(d);
  cout << "Energy:\n";
  dsaupd(d, NEV, val, vecs, av2);
  delete store;
  for(int i=0;i<NEV;i++) cout << "  " << val[i] << endl;
  data_energy << l << "\t" << l+1;
  for(int i=0;i<NEV;i++) data_energy << "\t" << val[i];
  data_energy << endl;

  qtensor<double> vec(vecs, dim, sym);
  qtensor<double> U, V, S;
  vector<double> s;
  s = vec.svd(U, S, V, 2, cutoff);
  double entropy = 0;
  for(int i=0;i<s.size();i++) entropy += -s[i]*s[i]*log(s[i]*s[i]);
  cout << "Singular value:\n";
  for(int i=0;i<s.size();i++) cout << "  " << s[i] << endl;
  cout << "Entropy: " << entropy << endl;
  data_singular << l << "\t" << l+1 << "\t" << entropy;
  int max = 20;
  if(s.size()>max) for(int i=0;i<max;i++) data_singular << "\t" << s[i];
  else
  {
    for(int i=0;i<s.size();i++) data_singular << "\t" << s[i];
    for(int i=s.size();i<max;i++) data_singular << "\t" << 0;
  }
  data_singular << endl;

  my_mps[l] = U;
  my_mps.center() = S;
  vec.contract(lenv, U, 'N', 'N', 2);
  my_mps(l).contract(U, vec, 'T', 'N', 2);
}


void move2left(int r, int cutoff, mps<double>& my_mps, mpo<double>& my_mpo)
{
  int d = 1;
  d *= my_mps(r-1).dimension(2);  // left environment
  d *= my_mpo[r].dimension(3);  // the right point
  d *= my_mps(r+1).dimension(2);  // right environment

  lenv = my_mps(r-1);
  qtensor<double> ope = my_mpo[r].exchange(0,2);
  renv.contract(my_mps(r+1), 1, ope, 0);
  renv = renv.exchange(0,1);
  
  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, lenv, 1, renv, 2);
  r_num = 2;
  l_num = 2;
  renv.get_map(rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', r_num, 0);
  lenv.get_map(lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', l_num, 0);

  double *val = new double [NEV];
  double *vecs = new double [d];
  int d_store = 1;
  for(int i=0;i<dim_ret.size();i++)
  {
    int temp = 0;
    for(int j=0;j<dim_ret[i].size();j++) temp += dim_ret[i][j];
    d_store *= temp;
  }
  store = new double [d_store];
//  check_hermitian(d);
  cout << "Energy:\n";
  dsaupd(d, NEV, val, vecs, av2);
  delete store;
  for(int i=0;i<NEV;i++) cout << "  " << val[i] << endl;
  data_energy << r-1 << "\t" << r;
  for(int i=0;i<NEV;i++) data_energy << "\t" << val[i];
  data_energy << endl;

  qtensor<double> vec(vecs, dim, sym);
  qtensor<double> U, V, S;
  vector<double> s;
  s = vec.svd(U, S, V, 1, cutoff);
  double entropy = 0;
  for(int i=0;i<s.size();i++) entropy += -s[i]*s[i]*log(s[i]*s[i]);
  cout << "Singular value:\n";
  for(int i=0;i<s.size();i++) cout << "  " << s[i] << endl;
  cout << "Entropy: " << entropy << endl;
  data_singular << r-1 << "\t" << r << "\t" << entropy;
  int max = 20;
  if(s.size()>max) for(int i=0;i<max;i++) data_singular << "\t" << s[i];
  else
  {
    for(int i=0;i<s.size();i++) data_singular << "\t" << s[i];
    for(int i=s.size();i<max;i++) data_singular << "\t" << 0;
  }
  data_singular << endl;

  my_mps.center() = S;
  my_mps[r] = V;
  vec.contract(renv, V, 'N', 'T', 2);
  my_mps(r).contract(V, vec, 'N', 'N', 2);
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
  cout << "********** starting l=" << 0 << " r=" << L-1 << " **********\n";
  two_sites(L, pre_cutoff, my_mps, my_mpo);
  for(int i=1;i<L/2;i++)
  {
    cout << "********** starting l=" << i << " r=" << L-1-i << " **********\n";
    update_two(i, L-i-1, pre_cutoff, my_mps, my_mpo);
  }
  for(int i=L/2;i<L-1;i++)
  {
    cout << "********** starting l=" << i << " r=" << i+1 << " **********\n";
    move2right(i, pre_cutoff, my_mps, my_mpo);
  }
  for(int i=L-2;i>0;i--)
  {
    cout << "********** sweep l=" << i << " r=" << i+1 << " **********\n";
    move2left(i, pre_cutoff, my_mps, my_mpo);
  }
  data_energy.close();
  data_singular.close();
}

