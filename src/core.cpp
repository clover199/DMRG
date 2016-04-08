
#include "global.h"
#include "mps.h"
#include "mpo.h"
#include "useful.h"
#include "basic.h"
#include "val_for_lanczos.h"

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

template <typename T>
void av(T *in, T *out, lanczos<T>& pass_val)
{
  pass_val.renv.contract(pass_val.store, in, pass_val.rmap, 'N', 'T', pass_val.r_num);
  pass_val.lenv.contract(out, pass_val.store, pass_val.lmap, 'N', 'T', pass_val.l_num);
#ifdef PBC
  pass_val.redge.contract(pass_val.store, in, pass_val.remap, 'N', 'T', pass_val.r_num);
  pass_val.ledge.contract(out, pass_val.store, pass_val.lemap, 'N', 'T', pass_val.l_num);
#endif
}

template void av(double *in, double *out, lanczos<double>& pass_val);
template void av(complex<double> *in, complex<double> *out, lanczos<complex<double> >& pass_val);


double my_conj(double x) { return x; }
complex<double> my_conj(complex<double> x) { return std::conj(x); }

template <typename T>
void check_hermitian(int n, lanczos<T>& pass_val)
{
  T* in = new T [n];
  T* out = new T [n];
  if(n>9) for(int i=0;i<n;i++)
  {
    for(int j=i+1;j<n;j++)
    {
      T temp;
      for(int s=0;s<n;s++) in[s] = 0;
      in[j] = 1;
      av(in, out, pass_val);
      temp = out[i];
      for(int s=0;s<n;s++) in[s] = 0;
      in[i] = 1;
      av(in, out, pass_val);
      if(abs(my_conj(temp)-out[j])>TOL)
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
        T temp;
        for(int s=0;s<n;s++) in[s] = 0;
        in[j] = 1;
        av(in, out, pass_val);
        temp = out[i];
        if(abs(out[i])>TOL) cout <<  out[i] << "\t";
        else { temp = 0; cout << temp << "\t";}
      }
      cout << endl;
    }
  }
  delete in;
  delete out;
}

template void check_hermitian(int n, lanczos<double>& pass_val);
template void check_hermitian(int n, lanczos<complex<double> >& pass_val);

template <typename T>
double two_sites(int l, int r, int cutoff, mps<T>& my_mps, mpo<T>& my_mpo,
               int sector, ofstream& data_energy, ofstream& data_singular)
{
  lanczos<T> pass_val;

  pass_val.lenv = my_mpo[l];
  pass_val.renv = my_mpo[r];
#ifdef PBC
  pass_val.ledge = my_mps.edge(l);
  pass_val.redge = my_mps.edge(r);
#endif

  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, pass_val.lenv, 1, pass_val.renv, 1, sector);
  int d = get_dimension(dim, sym);
  double *val = new double [NEV*10];
  T *vecs = new T [d];
  if(d<NEV*10)
  {
    qtensor<T> whole, edge, all;
    whole.contract(pass_val.lenv, 1, pass_val.renv, 1);
    whole = whole.simplify();
#ifdef PBC
    edge.contract(pass_val.ledge, 1, pass_val.redge, 1);
    edge = edge.simplify();
    all.plus(edge, whole);
    whole = all.simplify();
#endif
    whole = whole.exchange(2,3);
    whole = whole.combine(2,3);
    whole = whole.combine(0,1);
    cout << "Using ED, d=" << d << endl;
//    whole.print_matrix();
    whole.eig(val, vecs, sector);
    print_energy(data_energy, l, r, val, d);
  }
  else
  {
    pass_val.r_num = 1;
    pass_val.l_num = 2;
    pass_val.renv.get_map(pass_val.rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', pass_val.r_num, 0);
    int d_store = get_dimension(dim_ret, sym_ret);
    pass_val.lenv.get_map(pass_val.lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', pass_val.l_num, 0, true);
#ifdef PBC
    pass_val.redge.get_map(pass_val.remap, dim_ret, sym_ret, dim, sym, 'N', 'T', pass_val.r_num, 0);
    pass_val.ledge.get_map(pass_val.lemap, dim, sym, dim_ret, sym_ret, 'N', 'T', pass_val.l_num, 1, true);
#endif
    pass_val.store = new T [d_store];
//    check_hermitian(d, pass_val);
    cout << "Using Lanczos, d=" << d << endl;
    znaupd(d, NEV, val, vecs, pass_val);
    delete pass_val.store;
    print_energy(data_energy, l, r, val);
  }

  qtensor<T> vec(vecs, dim, sym);
  qtensor<T> U, V;
  qtensor<double> S;
  vector<double> s;
  s = vec.svd(U, S, V, 1, cutoff);
  print_singular(data_singular, l, r, s, cutoff);

  my_mps[l] = U;
  my_mps.center(l+1) = S;
  my_mps[r] = V;

  vec.contract(pass_val.lenv, U, 'N', 'N');
  U.conjugate();
  my_mps(l).contract(U, vec, 'T', 'N');
  vec.contract(pass_val.renv, V, 'N', 'T');
  V.conjugate();
  my_mps(r).contract(V, vec, 'N', 'N');
#ifdef PBC
  U.conjugate();
  vec.contract(pass_val.ledge, U, 'N', 'N');
  U.conjugate();
  my_mps.edge(l).contract(U, vec, 'T', 'N');
  V.conjugate();
  vec.contract(pass_val.redge, V, 'N', 'T');
  V.conjugate();
  my_mps.edge(r).contract(V, vec, 'N', 'N');
#endif
  return val[0];
}

template double two_sites(int l, int r, int cutoff,
               mps<double>& my_mps, mpo<double>& my_mpo,
               int sector, ofstream& data_energy, ofstream& data_singular);
template double two_sites(int l, int r, int cutoff,
               mps<complex<double> >& my_mps, mpo<complex<double> >& my_mpo,
               int sector, ofstream& data_energy, ofstream& data_singular);


template <typename T>
double update_two(int l, int r, int cutoff, mps<T>& my_mps, mpo<T>& my_mpo,
               int sector, ofstream& data_energy, ofstream& data_singular)
{
  lanczos<T> pass_val;

  pass_val.lenv.contract(my_mps(l-1), 1, my_mpo[l], 0);
  pass_val.lenv = pass_val.lenv.exchange(3,4);
  qtensor<T> ope = my_mpo[r].exchange(0,2);
  pass_val.renv.contract(my_mps(r+1), 1, ope, 0, true);
  pass_val.renv = pass_val.renv.exchange(0,1);
#ifdef PBC
  ope = my_mps.edge(l-1).id(my_mpo[l]);
  pass_val.ledge.contract(my_mps.edge(l-1), 1, ope, 0);
  pass_val.ledge = pass_val.ledge.exchange(3,4);
  pass_val.redge.contract(my_mps.edge(r+1), 1, ope, 0, true);
  pass_val.redge = pass_val.redge.exchange(0,1);
#endif
  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, pass_val.lenv, 2, pass_val.renv, 2, sector);
  int d = get_dimension(dim, sym);
  double *val = new double [NEV];
  T *vecs = new T [d];
  if(0)
  {  // check energy using ED
    val = new double [d];
    qtensor<T> whole, tlenv, trenv;
    tlenv = pass_val.lenv.combine(3,4);
    tlenv = tlenv.combine(0,1);
    trenv = pass_val.renv.combine(3,4);
    trenv = trenv.combine(0,1);
    whole.contract(tlenv, 1, trenv, 1);
    whole = whole.simplify();
#ifdef PBC
    qtensor<T> tledge, tredge, edge, all;
    tledge = pass_val.ledge.combine(3,4);
    tledge = tledge.combine(0,1);
    tredge = pass_val.redge.combine(3,4);
    tredge = tredge.combine(0,1);
    edge.contract(tledge, 1, tredge, 1);
    edge = edge.simplify();
    all.plus(whole, edge);
    whole = all.simplify();
#endif
    whole = whole.exchange(2,3);
    whole = whole.combine(2,3);
    whole = whole.combine(0,1);
    cout << "Using ED, d=" << d << endl;
//    whole.print_matrix();
    whole.eig(val, vecs, sector);
    print_energy(data_energy, l, r, val, NEV);
  }
//  else
  {
    pass_val.r_num = 2;
    pass_val.l_num = 3;
    pass_val.renv.get_map(pass_val.rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', pass_val.r_num, 0);
    int d_store = get_dimension(dim_ret, sym_ret);
    pass_val.lenv.get_map(pass_val.lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', pass_val.l_num, 0, true);
#ifdef PBC
    pass_val.redge.get_map(pass_val.remap, dim_ret, sym_ret, dim, sym, 'N', 'T', pass_val.r_num, 0);
    pass_val.ledge.get_map(pass_val.lemap, dim, sym, dim_ret, sym_ret, 'N', 'T', pass_val.l_num, 1, true);
#endif
    pass_val.store = new T [d_store];
//    check_hermitian(d, pass_val);
    cout << "Using Lanczos, d=" << d << endl;
    znaupd(d, NEV, val, vecs, pass_val);
    delete pass_val.store;
    print_energy(data_energy, l, r, val);
  }

  qtensor<T> vec(vecs, dim, sym);
  qtensor<T> U, V;
  qtensor<double> S;
  vector<double> s;
  s = vec.svd(U, S, V, 2, cutoff);
  print_singular(data_singular, l, r, s, cutoff);

  my_mps[l] = U;
  my_mps.center(l+1) = S;
  my_mps[r] = V;
  vec.contract(pass_val.lenv, U, 'N', 'N', 2);
  U.conjugate();
  my_mps(l).contract(U, vec, 'T', 'N', 2);
  vec.contract(pass_val.renv, V, 'N', 'T', 2);
  V.conjugate();
  my_mps(r).contract(V, vec, 'N', 'N', 2);
#ifdef PBC
  U.conjugate();
  vec.contract(pass_val.ledge, U, 'N', 'N', 2);
  U.conjugate();
  my_mps.edge(l).contract(U, vec, 'T', 'N', 2);
  V.conjugate();
  vec.contract(pass_val.redge, V, 'N', 'T', 2);
  V.conjugate();
  my_mps.edge(r).contract(V, vec, 'N', 'N', 2);
#endif
  return val[0];
}

template double update_two(int l, int r, int cutoff,
                mps<double>& my_mps, mpo<double>& my_mpo,
                int sector, ofstream& data_energy, ofstream& data_singular);
template double update_two(int l, int r, int cutoff,
                mps<complex<double> >& my_mps, mpo<complex<double> >& my_mpo,
                int sector, ofstream& data_energy, ofstream& data_singular);


template <typename T>
double move2right(int l, int r, int cutoff, mps<T>& my_mps, mpo<T>& my_mpo,
                int sector, ofstream& data_energy, ofstream& data_singular)
{
  lanczos<T> pass_val;

  pass_val.lenv.contract(my_mps(l-1), 1, my_mpo[l], 0);
  pass_val.lenv = pass_val.lenv.exchange(3,4);
  pass_val.renv = my_mps(r);
#ifdef PBC
  qtensor<T> ope = my_mps.edge(l-1).id(my_mpo[l]);
  pass_val.ledge.contract(my_mps.edge(l-1), 1, ope, 0);
  pass_val.ledge = pass_val.ledge.exchange(3,4);
  pass_val.redge = my_mps.edge(r);
#endif

  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, pass_val.lenv, 2, pass_val.renv, 1, sector);
  int d = get_dimension(dim, sym);
  double *val = new double [NEV];
  T *vecs = new T [d];
  if(d<NEV)
  {  // check energy using ED
    qtensor<T> whole, tlenv;
    tlenv = pass_val.lenv.combine(3,4);
    tlenv = tlenv.combine(0,1);
    whole.contract(tlenv, 1, pass_val.renv, 1);
    whole = whole.simplify();
#ifdef PBC
    qtensor<T> tledge, edge, all;
    tledge = pass_val.ledge.combine(3,4);
    tledge = tledge.combine(0,1);
    edge.contract(tledge, 1, pass_val.redge, 1);
    edge = edge.simplify();
    all.plus(whole, edge);
    whole = all.simplify();
#endif
    whole = whole.exchange(2,3);
    whole = whole.combine(2,3);
    whole = whole.combine(0,1);
    cout << "Using ED, d=" << d << endl;
    whole.eig(val, vecs, sector);
    print_energy(data_energy, l, r, val, d);
  }
  else
  {
    pass_val.r_num = 1;
    pass_val.l_num = 3;
    pass_val.renv.get_map(pass_val.rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', pass_val.r_num, 0);
    int d_store = get_dimension(dim_ret, sym_ret);
    pass_val.lenv.get_map(pass_val.lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', pass_val.l_num, 0, true);
#ifdef PBC
    pass_val.redge.get_map(pass_val.remap, dim_ret, sym_ret, dim, sym, 'N', 'T', pass_val.r_num, 0);
    pass_val.ledge.get_map(pass_val.lemap, dim, sym, dim_ret, sym_ret, 'N', 'T', pass_val.l_num, 1, true);
#endif
    pass_val.store = new T [d_store];
  //  check_hermitian(d, pass_val);
    cout << "Using Lanczos, d=" << d << endl;
    znaupd(d, NEV, val, vecs, pass_val);
    delete pass_val.store;
    print_energy(data_energy, l, r, val);
  }

  qtensor<T> vec(vecs, dim, sym);
  qtensor<T> U, V;
  qtensor<double> S;
  vector<double> s;
  s = vec.svd(U, S, V, 2, cutoff);
  print_singular(data_singular, l, r, s, cutoff);

  my_mps[l] = U;
//  my_mps[r] = V;
  my_mps.center(l+1) = S;
  vec.contract(pass_val.lenv, U, 'N', 'N', 2);
  U.conjugate();
  my_mps(l).contract(U, vec, 'T', 'N', 2);
#ifdef PBC
  U.conjugate();
  vec.contract(pass_val.ledge, U, 'N', 'N', 2);
  U.conjugate();
  my_mps.edge(l).contract(U, vec, 'T', 'N', 2);
#endif
  return val[0];
}

template double move2right(int l, int r, int cutoff,
                mps<double>& my_mps, mpo<double>& my_mpo,
                int sector, ofstream& data_energy, ofstream& data_singular);
template double move2right(int l, int r, int cutoff,
                mps<complex<double> >& my_mps, mpo<complex<double> >& my_mpo,
                int sector, ofstream& data_energy, ofstream& data_singular);


template <typename T>
double move2left(int l, int r, int cutoff, mps<T>& my_mps, mpo<T>& my_mpo,
               int sector, ofstream& data_energy, ofstream& data_singular)
{
  lanczos<T> pass_val;

  pass_val.lenv = my_mps(l);
  qtensor<T> ope = my_mpo[r].exchange(0,2);
  pass_val.renv.contract(my_mps(r+1), 1, ope, 0, true);
  pass_val.renv = pass_val.renv.exchange(0,1);
#ifdef PBC
  pass_val.ledge = my_mps.edge(l);
  ope = my_mps.edge(r+1).id(my_mpo[r]);
  pass_val.redge.contract(my_mps.edge(r+1), 1, ope, 0, true);
  pass_val.redge = pass_val.redge.exchange(0,1);
#endif

  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, pass_val.lenv, 1, pass_val.renv, 2, sector);
  int d = get_dimension(dim, sym);
  double *val = new double [NEV];
  T *vecs = new T [d];
  if(d<NEV)
  {  // check energy using ED
    qtensor<T> whole, trenv;
    trenv = pass_val.renv.combine(3,4);
    trenv = trenv.combine(0,1);
    whole.contract(pass_val.lenv, 1, trenv, 1);
    whole = whole.simplify();
#ifdef PBC
    qtensor<T> tredge, edge, all;
    tredge = pass_val.redge.combine(3,4);
    tredge = tredge.combine(0,1);
    edge.contract(pass_val.ledge, 1, tredge, 1);
    edge = edge.simplify();
    all.plus(whole, edge);
    whole = all.simplify();
#endif
    whole = whole.exchange(2,3);
    whole = whole.combine(2,3);
    whole = whole.combine(0,1);
    val = new double [d];
    cout << "Using ED, d=" << d << endl;
    whole.eig(val, vecs, sector);
    print_energy(data_energy, l, r, val, d);
  }
  else
  {
    pass_val.r_num = 2;
    pass_val.l_num = 2;
    pass_val.renv.get_map(pass_val.rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', pass_val.r_num, 0);
    int d_store = get_dimension(dim_ret, sym_ret);
    pass_val.lenv.get_map(pass_val.lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', pass_val.l_num, 0, true);
#ifdef PBC
    pass_val.redge.get_map(pass_val.remap, dim_ret, sym_ret, dim, sym, 'N', 'T', pass_val.r_num, 0);
    pass_val.ledge.get_map(pass_val.lemap, dim, sym, dim_ret, sym_ret, 'N', 'T', pass_val.l_num, 1, true);
#endif
    pass_val.store = new T [d_store];
//    check_hermitian(d, pass_val);
    cout << "Using Lanczos, d=" << d << endl;
    znaupd(d, NEV, val, vecs, pass_val);
    delete pass_val.store;
    print_energy(data_energy, l, r, val);
  }

  qtensor<T> vec(vecs, dim, sym);
  qtensor<T> U, V;
  qtensor<double> S;
  vector<double> s;
  s = vec.svd(U, S, V, 1, cutoff);
  print_singular(data_singular, l, r, s, cutoff);

  my_mps.center(l+1) = S;
  my_mps[r] = V;
  vec.contract(pass_val.renv, V, 'N', 'T', 2);
  V.conjugate();
  my_mps(r).contract(V, vec, 'N', 'N', 2);
#ifdef PBC
  V.conjugate();
  vec.contract(pass_val.redge, V, 'N', 'T', 2);
  V.conjugate();
  my_mps.edge(r).contract(V, vec, 'N', 'N', 2);
#endif
  return val[0];
}

template double move2left(int l, int r, int cutoff,
               mps<double>& my_mps, mpo<double>& my_mpo,
               int sector, ofstream& data_energy, ofstream& data_singular);
template double move2left(int l, int r, int cutoff,
               mps<complex<double> >& my_mps, mpo<complex<double> >& my_mpo,
               int sector, ofstream& data_energy, ofstream& data_singular);


template <typename T>
void grow2right(int l, mps<T>& my_mps, mpo<T>& my_mpo)
{
  if(l==0)  my_mps(l) = my_mpo[l];
  else
  {
    qtensor<T> temp;
    temp.contract(my_mps(l-1), 1, my_mpo[l], 0);
    temp = temp.exchange(3,4);
    temp = temp.combine(3,4);
    my_mps(l) = temp.combine(0,1);
#ifdef PBC
    qtensor<T> ope = my_mps.edge(l-1).id(my_mpo[l]);
    temp.contract(my_mps.edge(l-1), 1, ope, 0);
    temp = temp.exchange(3,4);
    temp = temp.combine(3,4);
    my_mps.edge(l) = temp.combine(0,1);
#endif
    my_mps.clear(l-1);
  }
}

template void grow2right(int l, mps<double>& my_mps, mpo<double>& my_mpo);
template void grow2right(int l, mps<complex<double> >& my_mps, mpo<complex<double> >& my_mpo);


template <typename T>
void grow2left(int r, mps<T>& my_mps, mpo<T>& my_mpo)
{
  if(r==(my_mps.size()-1)) my_mps(r) = my_mpo[r];
  else
  {
    qtensor<T> temp, ope;
    ope = my_mpo[r].exchange(0,2);
    temp.contract(my_mps(r+1), 1, ope, 0, true);
    temp = temp.exchange(0,1);
    temp = temp.combine(3,4);
    my_mps(r) = temp.combine(0,1);

#ifdef PBC
    ope = my_mps.edge(r+1).id(my_mpo[r]);
    temp.contract(my_mps.edge(r+1), 1, ope, 0, true);
    temp = temp.exchange(0,1);
    temp = temp.combine(3,4);
    my_mps.edge(r) = temp.combine(0,1);
#endif
    my_mps.clear(r+1);
  }
}

template void grow2left(int r, mps<double>& my_mps, mpo<double>& my_mpo);
template void grow2left(int r, mps<complex<double> >& my_mps, mpo<complex<double> >& my_mpo);


template <typename T>
void two_sites(int l, int r, mps<T>& my_mps, int sector, 
               ofstream& data_energy, ofstream& data_singular)
{
  lanczos<T> pass_val;

  pass_val.lenv = my_mps(l);
  pass_val.renv = my_mps(r);
#ifdef PBC
  pass_val.ledge = my_mps.edge(l);
  pass_val.redge = my_mps.edge(r);
#endif

  vector< vector<int> > dim, sym, dim_ret, sym_ret;
  generate_dim_sym(dim, sym, pass_val.lenv, 1, pass_val.renv, 1, sector);
  int d = get_dimension(dim, sym);
  double *val = new double [NEV*10];
  T *vecs = new T [d];
  if(d<2000)
  {
    qtensor<T> whole, edge, all;
    whole.contract(pass_val.lenv, 1, pass_val.renv, 1);
    whole = whole.simplify();
#ifdef PBC
    edge.contract(pass_val.ledge, 1, pass_val.redge, 1);
    edge = edge.simplify();
    all.plus(edge, whole);
    whole = all.simplify();
#endif
    whole = whole.exchange(2,3);
    whole = whole.combine(2,3);
    whole = whole.combine(0,1);
    cout << "Using ED, d=" << d << endl;
//    whole.print_matrix();
    whole.eig(val, vecs, sector);
    print_energy(data_energy, l, r, val, d);
  }
  else
  {
    pass_val.r_num = 1;
    pass_val.l_num = 2;
    pass_val.renv.get_map(pass_val.rmap, dim_ret, sym_ret, dim, sym, 'N', 'T', pass_val.r_num, 0);
    int d_store = get_dimension(dim_ret, sym_ret);
    pass_val.lenv.get_map(pass_val.lmap, dim, sym, dim_ret, sym_ret, 'N', 'T', pass_val.l_num, 0, true);
#ifdef PBC
    pass_val.redge.get_map(pass_val.remap, dim_ret, sym_ret, dim, sym, 'N', 'T', pass_val.r_num, 0);
    pass_val.ledge.get_map(pass_val.lemap, dim, sym, dim_ret, sym_ret, 'N', 'T', pass_val.l_num, 1, true);
#endif
    pass_val.store = new T [d_store];
//    check_hermitian(d, pass_val);
    cout << "Using Lanczos, d=" << d << endl;
    znaupd(d, NEV, val, vecs, pass_val);
    delete pass_val.store;
    print_energy(data_energy, l, r, val);
  }

  qtensor<T> vec(vecs, dim, sym);
  qtensor<T> U, V;
  qtensor<double> S;
  vector<double> s;
  s = vec.svd(U, S, V, 1, 0);
  print_singular(data_singular, l, r, s, s.size());
}

template void two_sites(int l, int r, mps<double>& my_mps, int sector,
               ofstream& data_energy, ofstream& data_singular);
template void two_sites(int l, int r, mps<complex<double> >& my_mps, int sector,
               ofstream& data_energy, ofstream& data_singular);

