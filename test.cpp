
#include "global.h"
#include "tensor.h"

int main()
{
  srand(time(0));

  int a1=3, a2=4, a3=5;
  stensor< complex<double> > A(a1,a2,a3);
  for(int i=0;i<10;i++) A.update(5.0*rand()/RAND_MAX-10, rand()%a1, rand()%a2, rand()%a3);
  A.sort();
//  cout << "A" << endl; A.print();

  int b1=5, b2=5;
  stensor< complex<double> > B(b1,b2);
  for(int i=0;i<20;i++) B.update(5.0*rand()/RAND_MAX-10, rand()%b1, rand()%b2);
  B.sort();
//  cout << "B" << endl; B.print();  

  int c1=5;
  stensor< complex<double> > C(c1);
  for(int i=0;i<10;i++) C.update(5.0*rand()/RAND_MAX-10, rand()%c1);
  C.sort();
//  cout << "C" << endl; C.print();

  int d=4;
  stensor< complex<double> > D(d);
  for(int i=0;i<d;i++) D.update(1, i);
  D.diag(3);
//  cout << "D" << endl; D.print();

/*
  tensor< double > S;
  A.svd(B,S,C,2);
  S.print();
  B= B.transpose(2);
  C = C.transpose();
  D.contract(B,A,2);
  B.contract(D,C);
  B.print();*/
}