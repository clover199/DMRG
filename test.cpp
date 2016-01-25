
#include "global.h"
#include "functions.h"
#include "operators.h"
#include "useful.h"

extern string filename = "";
extern bool print = false;
extern int symmetry_sector = 0;

int main()
{
  std::srand(std::time(0));
 
  qtensor<double> A(2,2,2);
  tensor<double> A0(2,2,2);
  A0.rand();
  A.update(A0, 0,0,0);
  A0.rand();
  A.update(A0, 0,0,1);
  A0.rand();
  A.update(A0, 1,1,0);
  A0.rand();
  A.update(A0, 1,0,1);
//  cout << "\ntensor A:\n"; A.print_matrix(); cout << endl;

  qtensor<double> B(2,2,2);
  A0.rand();
  B.update(A0, 0,0,1);
  A0.rand();
  B.update(A0, 0,1,0);
  A0.rand();
  B.update(A0, 1,0,0);
  A0.rand();
  B.update(A0, 1,0,1);

  char transa = 'T';
  char transb = 'N';
  int num = 2;
  qtensor<double> D;
  D.contract(A, B, transa, transb, num);
//  cout << "\ntensor D:\n"; D.print_matrix(); cout << endl;

  double *in = new double [B.dimension()];
  double *out = new double [D.dimension()];
  vector< vector<int> > dim;
  vector< vector<int> > sym;
  int d = B.val(in, dim, sym);
  vector< vector<int> > dim_ret;
  vector< vector<int> > sym_ret;
  vector< vector<int> > map;
  A.get_map(map, dim_ret, sym_ret, dim, sym, transa, transb, num, 0, 0);
  A.contract(out, in, map, transa, transb, num);
  qtensor<double> E(out, dim_ret, sym_ret);
//  cout << "\ntensor E:\n"; E.print_matrix();
  A.plus(D,E,-1);
  A.print_all();
}
