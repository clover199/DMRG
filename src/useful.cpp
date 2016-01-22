#include "global.h"
#include "qtensor.h"

#ifdef SYMMETRY
  int add(int a, int b) {return (a+b)%SYMMETRY;}
#elif defined FERMION
  int add(int a, int b) {return (a+b)%2;}
#else
  int add(int a, int b) {return a+b;}
#endif


void print_matrix(const vector< vector<int> >& map)
{
  cout << "Print matrix:" << endl;
  for(int i=0;i<map.size();i++)
  {
    for(int j=0;j<map[i].size();j++) cout << map[i][j] << " ";
    cout << endl;
  }
  cout << endl;
}


// generate the full array of all possible symmetries (sorted) given index 
// i.e. index=[4,1,3], then
// map=[ [0,0,0]
//       [0,0,1]
//       [0,0,2]
//       [1,0,0]
//       [1,0,1]
//       [1,0,2]
//       [2,0,0]
//       [2,0,1]
//       [2,0,2]
//       [3,0,0]
//       [3,0,1]
//       [3,0,2] ]
void generate_map(vector< vector<int> >& map, const vector<int>& index)
{
  int dim = 1;  // dimension of map
  for(int i=0;i<index.size();i++) dim *= index[i];
  map.resize(dim);
  for(int i=0;i<dim;i++) map[i].resize(index.size());
  int count = 1;  // the periodicity of the symmetry for the index i
  for(int i=0;i<index.size();i++)
  {
    for(int p=0;p<count;p++)
      for(int s=0;s<index[i];s++)
        for(int k=0;k<dim/count/index[i];k++)
          map[k+s*dim/count/index[i]+p*dim/count][i] = s;
    count *= index[i];
  }
}


int get_dimension(const vector< vector<int> >& dim, const vector< vector<int> >& sym)
{
  int d = 0;
  for(int i=0;i<sym.size();i++) 
  {
    int temp = 1;
    for(int j=0;j<dim.size();j++) temp *= dim[j][sym[i][j]];
    d += temp;
  }
  return d;
}


template <typename T>
void generate_dim_sym(vector< vector<int> >& dim, vector< vector<int> >& sym,
                      const qtensor<T>& lenv, int lnum,
                      const qtensor<T>& renv, int rnum)
{
  dim.resize(lnum+rnum);
  for(int i=0;i<lnum;i++)
  {
    int temp = lenv.index()-lnum+i;
    dim[i].resize(lenv.index(temp));
    for(int j=0;j<dim[i].size();j++) dim[i][j] = lenv.dimension(temp,j);
  }
  for(int i=0;i<rnum;i++)
  {
    int temp = renv.index()-rnum+i;
    dim[lnum+i].resize(renv.index(temp));
    for(int j=0;j<dim[lnum+i].size();j++) dim[lnum+i][j] = renv.dimension(temp,j);
  }

  vector<int> index(dim.size(), 0);  // used as index here
  for(int i=0;i<dim.size();i++) index[i] = dim[i].size();
  vector< vector<int> > sym_ret;
  generate_map(sym_ret, index);
  sym.clear();
  if(symmetry_sector!=-1) for(int i=0;i<sym_ret.size();i++)
  {
    int sum = 0;
    for(int j=0;j<dim.size();j++) sum = add(sum, sym_ret[i][j]);
    if(sum==symmetry_sector) sym.push_back(sym_ret[i]);
  }
  else
  {
    sym.resize(1);
    sym[0] = vector<int> (dim.size(), 0);
  }
}


template void generate_dim_sym(vector< vector<int> >& dim,
                               vector< vector<int> >& sym,
                               const qtensor<double>& lenv, int lnum,
                               const qtensor<double>& renv, int rnum);
template void generate_dim_sym(vector< vector<int> >& dim,
                               vector< vector<int> >& sym,
                               const qtensor<complex<double> >& lenv, int lnum,
                               const qtensor<complex<double> >& renv, int rnum);

