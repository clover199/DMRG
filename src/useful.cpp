#include "global.h"
#include "qtensor.h"

#ifdef S4E

int ss(int a, int b)
{
// 000 001 010 011 100 101 110 111
//  0   1   2   3   4   5   6   7
  int p[8][8] = {{0,1,2,3,4,5,6,7},
                 {1,0,3,2,5,4,7,6},
                 {2,3,0,1,6,7,4,5},
                 {3,2,1,0,7,6,5,4},
                 {4,5,6,7,0,1,2,3},
                 {5,4,7,6,1,0,3,2},
                 {6,7,4,5,2,3,0,1},
                 {7,6,5,4,3,2,1,0}};
  if(b<0) return p[a][-b];
  else return p[a][b];
}

int ff(int a)
{
  int f[8] = {0,0,1,1,1,1,0,0};
  return f[a];
}

#else
#ifdef SYMMETRY
  int ss(int a, int b) {return (a+b+SYMMETRY)%SYMMETRY;}

  int ff(int a) {return a%2;}

#else
  int ss(int a, int b) {return a+b;}

  int ff(int a) {return a%2;}

#endif
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


void print_energy(ofstream& data, int l ,int r, double* val, int n)
{
  cout.precision(ACC);
  data.precision(ACC);
  data.setf(std::ios::fixed);
  data.setf(std::ios::showpoint);

  cout << "Energy:\n";
  if(PRINT_TO_SCREEN) for(int i=0;i<n;i++) cout << "  " << val[i] << endl;
  else cout << val[0] << endl;
  data << l << "\t" << r;
  if(NEV<=n) for(int i=0;i<NEV;i++) data << "\t" << val[i];
  else
  {
    for(int i=0;i<n;i++) data << "\t" << val[i];
    for(int i=n;i<NEV;i++) data << "\t" << 0;
  }
  data << endl;
}


void print_singular(ofstream& data, int l, int r, const vector<double>& s, int cutoff)
{
  cout.precision(ACC);
  data.precision(ACC);
  data.setf(std::ios::fixed);
  data.setf(std::ios::showpoint);

  if(PRINT_TO_SCREEN)
  {
    cout << "Singular value:\n";
    for(int i=0;i<s.size();i++) if(s[i]>TOL) cout << "  " << s[i] << endl;
  }
  double entropy = 0;
  for(int i=0;i<s.size();i++) entropy += -s[i]*s[i]*log(s[i]*s[i]+TOL);
  cout << "Entropy: " << entropy << endl;
  data << l << "\t" << r << "\t" << entropy;
  if(s.size()>=NSI) for(int i=0;i<NSI;i++) data << "\t" << s[i];
  else
  {
    for(int i=0;i<s.size();i++) data << "\t" << s[i];
    for(int i=s.size();i<NSI;i++) data << "\t" << 0;
  }
  double error = 1;
  cutoff = s.size() < cutoff ? s.size() : cutoff;
  for(int i=0;i<cutoff;i++) error -= s[i]*s[i];
  cout << "Truncation error: " << error << endl;
  data << "\t" << error << endl;
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
  for(int i=0;i<index.size();i++)  // loop over each index
  {
    for(int p=0;p<count;p++)  // loop over each period. for index i, we have p(i) periods
      for(int s=0;s<index[i];s++)  // loop over index within periods
        for(int k=0;k<dim/count/index[i];k++)
          map[k+s*(dim/count/index[i])+p*(dim/count)][i] = s;
    count *= index[i];  // p(i+1) = p(i) * index[i]
  }
}


int get_dimension(const vector< vector<int> >& dim,
                  const vector< vector<int> >& sym)
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
                      const qtensor<T>& renv, int rnum,
                      int sector)
{
  vector< vector<int> > ldim = lenv.dimension_all();
  vector< vector<int> > rdim = renv.dimension_all();
  dim.resize(lnum+rnum);
  for(int i=0;i<lnum;i++) dim[i] = ldim[ldim.size()-lnum+i];
  for(int i=0;i<rnum;i++) dim[lnum+i] = rdim[rdim.size()-rnum+i];

  vector<int> index(lnum, 0);
  for(int i=0;i<lnum;i++) index[i] = dim[i].size();
  generate_map(ldim, index);  // generate the sym for the left part
  index.resize(rnum, 0);
  for(int i=0;i<rnum;i++) index[i] = dim[lnum+i].size();
  generate_map(rdim, index);  // generate the sym for the right part

  sym.clear();
  index.resize(lnum+rnum);  // used to store the symmetry sector
  if(sector!=-1) for(int s=0;s<SYMMETRY;s++)
  {
    int t = ss(sector, -s);
    for(int l=0;l<ldim.size();l++) for(int r=0;r<rdim.size();r++)
    {
      int lsym = 0;
      for(int j=0;j<lnum;j++) lsym = ss(lsym, ldim[l][j]);
      int rsym = 0;
      for(int j=0;j<rnum;j++) rsym = ss(rsym, rdim[r][j]);
      if(lsym==s and rsym==t)
      {
        for(int i=0;i<lnum;i++) index[i] = ldim[l][i];
        for(int i=0;i<rnum;i++) index[lnum+i] = rdim[r][i];
        sym.push_back(index);
      }
    }
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
                               const qtensor<double>& renv, int rnum,
                               int sector);
template void generate_dim_sym(vector< vector<int> >& dim,
                               vector< vector<int> >& sym,
                               const qtensor<complex<double> >& lenv, int lnum,
                               const qtensor<complex<double> >& renv, int rnum,
                               int sector);

