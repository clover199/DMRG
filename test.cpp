
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
  cout << "\ntensor A:\n";
  qtensor<double> A(3,3,1,1);
  tensor<double> A0(3,2,4,1);
  A0.rand();
  A.update(A0, 0,1,0,0);
  A0 = tensor<double> (2,1,4,1);
  A0.rand();
  A.update(A0, 1,0,0,0);
  A0 = tensor<double> (3,2,4,1);
  A0.rand();
  A.update(A0, 2,2,0,0);
  A0 = tensor<double> (2,2,4,1);
  A0.rand();
//  A.update(A0, 1,1,0,0);
  A.print();
  cout << endl;

  qtensor<double> U,S,V,B;
  A.svd(U,S,V,2,10);
  cout<<"S:\n";S.print_all();cout<<endl;
  B = U * S * V;
  S.plus(A,B,-1);
  S.print_all(1);
}