
#include "global.h"
#include "tensor_quantum.h"
/*
template <typename T>
void update_one(vector<tensor<T> >& mps, vector<tensor>& mpo,
                int l, int r, int cutoff, vector<double> & data_ret)
{
  matrix U, V;
  int dim_ope = ta_ope.size();

  int left = mps[l].dimension();
  int mid = mps[l-1].dimension();
  int right = mps[r].dimension();
  complex<double> *vec = new complex<double> [left*mid*right];
  energy(mps[l], mps[l+1], mps[r], mpo[l], mpo[l+1], mpo[r], vec, data_ret);
  
  grow(s_ope[l], ta_ope, para, s_ope[l+1][0]);
  for(int i=1;i<dim_ope;i++)
  {
    s_ope[l+1][i].direct(s_ope[l][i].row(), tb_ope[i]);
    e_ope[l+1][i].direct(e_ope[l][i], tb_ope[i].row());
  }
  svd(vec, s_ope[l+1][0].row(), s_ope[r][0].row(), U, V, cutoff, data_ret);
  for(int i=0;i<dim_ope;i++) s_ope[l+1][i].transform(U); 
  for(int i=1;i<dim_ope;i++) e_ope[l+1][i].transform(U); 

  mps = U;
  delete vec;
}

void update_two(vector< vector<matrix> > & s_ope, vector< vector<matrix> >  & e_ope,
               vector<matrix> & ta_ope, vector<matrix> & tb_ope, matrix & mps,
               int l, int r, int cutoff, complex<double> *para, vector<double> & data_ret)
{
  int left = s_ope[l][0].row();
  int mid = ta_ope[0].row();
  int right = s_ope[r][0].row();
  complex<double> *vec = new complex<double> [left*mid*right];
  matrix U, V;
  int dim_ope = ta_ope.size();

  energy(e_ope[l], s_ope[l], ta_ope, tb_ope, s_ope[r], e_ope[r], para, vec, data_ret);
  grow(tb_ope, s_ope[r], para, s_ope[r-1][0]);
  for(int i=1;i<dim_ope;i++)
  {
    s_ope[r-1][i].direct(ta_ope[i], s_ope[r][i].row());
    e_ope[r-1][i].direct(ta_ope[i].row(), e_ope[r][i]);
  }
  svd(vec, s_ope[l][0].row(), s_ope[r-1][0].row(), U, V, cutoff, data_ret);
  for(int i=0;i<dim_ope;i++) s_ope[r-1][i].transform(V);
  for(int i=1;i<dim_ope;i++) e_ope[r-1][i].transform(V);

  mps = V;
  delete vec;
}

void move2left(vector< vector<matrix> > & s_ope, vector< vector<matrix> >  & e_ope,
               vector<matrix> & ta_ope, vector<matrix> & tb_ope, matrix & mps,
               int l, int r, int cutoff, complex<double> *para, vector<double> & data_ret)
{
  int left = s_ope[l][0].row();
  int mid = ta_ope[0].row();
  int right = s_ope[r][0].row();
  complex<double> *vec = new complex<double> [left*mid*right];
  matrix U, V;
  int dim_ope = ta_ope.size();

  energy(e_ope[l], s_ope[l], ta_ope, tb_ope, s_ope[r], e_ope[r], para, vec, data_ret);
  grow(tb_ope, s_ope[r], para, s_ope[r-1][0]);
  for(int i=1;i<dim_ope;i++)
  {
    s_ope[r-1][i].direct(ta_ope[i], s_ope[r][i].row());
    e_ope[r-1][i].direct(ta_ope[i].row(), e_ope[r][i]);
  }
  svd(vec, s_ope[l][0].row(), s_ope[r-1][0].row(), U, V, cutoff, data_ret);
  for(int i=0;i<dim_ope;i++) s_ope[r-1][i].transform(V);
  for(int i=1;i<dim_ope;i++) e_ope[r-1][i].transform(V);

  mps = V;
  delete vec;
}
*/