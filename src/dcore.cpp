
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
  if(n>9) for(int i=0;i<n;i++)
  {
    for(int j=i+1;j<n;j++)
    {
      double temp;
      for(int s=0;s<n;s++) in[s] = 0;
      in[j] = 1;
      av2(n, in, out);
      temp = out[i];
      for(int s=0;s<n;s++) in[s] = 0;
      in[i] = 1;
      av2(n, in, out);
      if(abs(temp-out[j])>TOL)
        cout << "Not hermitian! " << i << " " << j << " " 
             << temp << " " << out[j] << endl;
      ;
    }
  }
  else
  {
    cout << "Print matrix:" << endl;
    for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
      {
        double temp;
        for(int s=0;s<n;s++) in[s] = 0;
        in[j] = 1;
        av2(n, in, out);
        temp = out[i];
        if(abs(out[i])>TOL) cout <<  out[i] << "\t";
        else cout << 0 << "\t";
      }
      cout << endl;
    }
  }
  delete in;
  delete out;
}


void print_energy(int l ,int r, double* val, int n=NEV)
{
  cout << "Energy:\n";
  for(int i=0;i<n;i++) cout << "  " << val[i] << endl;
  data_energy << l << "\t" << r;
  if(NEV<=n) for(int i=0;i<NEV;i++) data_energy << "\t" << val[i];
  else
  {
    for(int i=0;i<n;i++) data_energy << "\t" << val[i];
    for(int i=n;i<NEV;i++) data_energy << "\t" << 0;
  }
  data_energy << endl;
}


void print_singular(int l, int r, const vector<double>& s)
{
  cout << "Singular value:\n";
  for(int i=0;i<s.size();i++) cout << "  " << s[i] << endl;
  double entropy = 0;
  for(int i=0;i<s.size();i++) entropy += -s[i]*s[i]*log(s[i]*s[i]+TOL);
  cout << "Entropy: " << entropy << endl;
  data_singular << l << "\t" << r << "\t" << entropy;
  if(s.size()>=NSI) for(int i=0;i<NSI;i++) data_singular << "\t" << s[i];
  else
  {
    for(int i=0;i<s.size();i++) data_singular << "\t" << s[i];
    for(int i=s.size();i<NSI;i++) data_singular << "\t" << 0;
  }
  data_singular << endl;
}


void two_sites(int size, int cutoff,
               mps<double>& my_mps, mpo<double>& my_mpo,
               int sector)
{
  int l = 0;
  int r = size-1;
  lenv = my_mpo[l];
  renv = my_mpo[r];

  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, lenv, 1, renv, 1, sector);
  int d = get_dimension(dim, sym);
  double *val = new double [NEV*10];
  double *vecs = new double [d];
  if(d<NEV*10)
  {
    qtensor<double> whole;
    whole.contract(lenv, 1, renv, 1);
    whole = whole.exchange(2,3);
    whole = whole.combine(2,3);
    whole = whole.combine(0,1);
    whole = whole.simplify();
    whole.eig(val, vecs, sector);
    print_energy(l, r, val, d);
  }
  else
  {
    r_num = 1;
    l_num = 2;
    renv.get_map(rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', r_num, 0);
    lenv.get_map(lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', l_num, 0, true);
    int d_store = get_dimension(dim_ret, sym_ret);
    store = new double [d_store];
//    check_hermitian(d);
    cout << "Using Lanczos, d=" << d << endl;
    dsaupd(d, 1, val, vecs, av2);
    delete store;
    print_energy(l, r, val);
  }

  qtensor<double> vec(vecs, dim, sym);
  qtensor<double> S, U, V;
  vector<double> s;
  s = vec.svd(U, S, V, 1, cutoff);
  print_singular(l, r, s);

  my_mps[l] = U;
  my_mps.center() = S;
  my_mps[r] = V;
  vec.contract(lenv, U, 'N', 'N');
  my_mps(l).contract(U, vec, 'T', 'N');
  vec.contract(renv, V, 'N', 'T');
  my_mps(r).contract(V, vec, 'N', 'N');
}


void update_two(int l, int r, int cutoff,
                mps<double>& my_mps, mpo<double>& my_mpo,
               int sector)
{
  lenv.contract(my_mps(l-1), 1, my_mpo[l], 0);
  lenv = lenv.exchange(3,4);
  qtensor<double> ope = my_mpo[r].exchange(0,2);
  renv.contract(my_mps(r+1), 1, ope, 0);
  renv = renv.exchange(0,1);
  
  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, lenv, 2, renv, 2, sector);
  int d = get_dimension(dim, sym);
  double *val = new double [NEV*10];
  double *vecs = new double [d];
  if(d<NEV*10)
  {  // check energy using ED
    qtensor<double> whole, tlenv, trenv;
    tlenv = lenv.combine(3,4);
    tlenv = tlenv.combine(0,1);
    trenv = renv.combine(3,4);
    trenv = trenv.combine(0,1);
    whole.contract(tlenv, 1, trenv, 1);
    whole = whole.simplify();
    whole = whole.exchange(2,3);
    whole = whole.combine(2,3);
    whole = whole.combine(0,1);
    whole.eig(val, vecs, sector);
    print_energy(l, r, val, d);
  }
  else
  {
    r_num = 2;
    l_num = 3;
    renv.get_map(rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', r_num, 0);
    lenv.get_map(lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', l_num, 0, true);
    int d_store = get_dimension(dim_ret, sym_ret);
    store = new double [d_store];
//    check_hermitian(d);
    cout << "Using Lanczos, d=" << d << endl;
    dsaupd(d, NEV, val, vecs, av2);
    delete store;
    print_energy(l, r, val);
  }

  qtensor<double> vec(vecs, dim, sym);
  qtensor<double> U, V;
  qtensor<double> S;
  vector<double> s;
  s = vec.svd(U, S, V, 2, cutoff);
  print_singular(l, r, s);

  my_mps[l] = U;
  my_mps.center() = S;
  my_mps[r] = V;
  vec.contract(lenv, U, 'N', 'N', 2);
  my_mps(l).contract(U, vec, 'T', 'N', 2);
  vec.contract(renv, V, 'N', 'T', 2);
  my_mps(r).contract(V, vec, 'N', 'N', 2);
}


void move2right(int l, int r, int cutoff,
                mps<double>& my_mps, mpo<double>& my_mpo,
               int sector)
{
  lenv.contract(my_mps(l-1), 1, my_mpo[l], 0);
  lenv = lenv.exchange(3,4);
  renv = my_mps(r);
  
  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, lenv, 2, renv, 1, sector);
  int d = get_dimension(dim, sym);
  double *val = new double [NEV*10];
  double *vecs = new double [d];
  if(d<NEV*10)
  {  // check energy using ED
    qtensor<double> whole, tlenv;
    tlenv = lenv.combine(3,4);
    tlenv = tlenv.combine(0,1);
    whole.contract(tlenv, 1, renv, 1);
    whole = whole.simplify();
    whole = whole.exchange(2,3);
    whole = whole.combine(2,3);
    whole = whole.combine(0,1);
    whole.eig(val, vecs, sector);
    print_energy(l, l, val, d);
  }
  else
  {
    r_num = 1;
    l_num = 3;
    renv.get_map(rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', r_num, 0);
    lenv.get_map(lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', l_num, 0, true);
    int d_store = get_dimension(dim_ret, sym_ret);
    store = new double [d_store];
  //  check_hermitian(d);
    cout << "Using Lanczos, d=" << d << endl;
    dsaupd(d, NEV, val, vecs, av2);
    delete store;
    print_energy(l, l, val);
  }

  qtensor<double> vec(vecs, dim, sym);
  qtensor<double> U, V;
  qtensor<double> S;
  vector<double> s;
  s = vec.svd(U, S, V, 2, cutoff);
  print_singular(l, l, s);

  my_mps[l] = U;
  my_mps.center() = S;
//  my_mps[r] = V;
  vec.contract(lenv, U, 'N', 'N', 2);
  my_mps(l).contract(U, vec, 'T', 'N', 2);
//  vec.contract(renv, V, 'N', 'T');
//  my_mps(r).contract(V, vec, 'N', 'N');
}


void move2left(int l, int r, int cutoff,
               mps<double>& my_mps, mpo<double>& my_mpo,
               int sector)
{
  lenv = my_mps(l);
  qtensor<double> ope = my_mpo[r].exchange(0,2);
  renv.contract(my_mps(r+1), 1, ope, 0);
  renv = renv.exchange(0,1);
  
  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, lenv, 1, renv, 2, sector);
  int d = get_dimension(dim, sym);
  double *val = new double [NEV*10];
  double *vecs = new double [d];
  if(d<NEV*10)
  {  // check energy using ED
    qtensor<double> whole, trenv;
    trenv = renv.combine(3,4);
    trenv = trenv.combine(0,1);
    whole.contract(lenv, 1, trenv, 1);
    whole = whole.simplify();
    whole = whole.exchange(2,3);
    whole = whole.combine(2,3);
    whole = whole.combine(0,1);
    val = new double [d];
    whole.eig(val, vecs, sector);
    print_energy(r, r, val, d);
  }
  else
  {
    r_num = 2;
    l_num = 2;
    renv.get_map(rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', r_num, 0);
    lenv.get_map(lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', l_num, 0, true);
    int d_store = get_dimension(dim_ret, sym_ret);
    store = new double [d_store];
//    check_hermitian(d);
    cout << "Using Lanczos, d=" << d << endl;
    dsaupd(d, NEV, val, vecs, av2);
    delete store;
    print_energy(r, r, val);
  }

  qtensor<double> vec(vecs, dim, sym);
  qtensor<double> U, V;
  qtensor<double> S;
  vector<double> s;
  s = vec.svd(U, S, V, 1, cutoff);
  print_singular(r, r, s);

//  my_mps[l] = U;
  my_mps.center() = S;
  my_mps[r] = V;
//  vec.contract(lenv, U, 'N', 'N');
//  my_mps(l).contract(U, vec, 'T', 'N');
  vec.contract(renv, V, 'N', 'T', 2);
  my_mps(r).contract(V, vec, 'N', 'N', 2);
}


void dmrg(mps<double>& my_mps, mpo<double>& my_mpo,
          int cutoff, int sweep, int sector, const string& filename)
{
  int L = my_mps.size();

  string name = "energy"+filename;
  data_energy.open(name.c_str());
  name = "singular"+filename;
  data_singular.open(name.c_str());

  int pre_cutoff = 10;
  int pre_sweep = cutoff/50;
  
  cout << "********** starting l=" << 0 << " r=" << L-1 << " **********\n";
  two_sites(L, pre_cutoff, my_mps, my_mpo, sector);
  for(int i=1;i<L/2;i++)
  {
    cout << ">>>>>>>>>>> starting l=" << i << " r=" << L-1-i << " <<<<<<<<<<\n";
    update_two(i, L-i-1, pre_cutoff, my_mps, my_mpo, sector);
  }
  for(int i=L/2;i<L-1;i++)
  {
    cout << ">>>>>>>>>> starting l=" << i << " r=" << i+1 << " >>>>>>>>>>\n";
    move2right(i, i+1, pre_cutoff, my_mps, my_mpo, sector);
  }
  for(int i=L-2;i>0;i--)
  {
    cout << "<<<<<<<<<< sweep=1 l=" << i << " r=" << i+1 << " <<<<<<<<<<\n";
    move2left(i-1, i, cutoff, my_mps, my_mpo, sector);
  }
  for(int i=1;i<L-1;i++)
  {
    cout << ">>>>>>>>>> sweep=1 l=" << i << " r=" << i+1 << " >>>>>>>>>>\n";
    move2right(i, i+1, cutoff, my_mps, my_mpo, sector);
  }
  data_energy.close();
  data_singular.close();
}

