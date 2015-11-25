#ifndef _MATRIX_PRODUCT_STATE_
#define _MATRIX_PRODUCT_STATE_

#include "global.h"
#include "tensor_quantum.h"
#include "mpo.h"

class mps
{
	
private:

  int cut_;
  vector< qtensor< dtensor<double> > > dstate_;
  vector< qtensor< dtensor< complex<double> > > > zstate_;
  vector<double> S_;
  
public:

  // set the number of sites
  mps(int sites);

  // initialize the structure of the MPS according to MPO
  // return the size of the system
  int initialize(mpo& in);
  
  void create_store(vector<tensor>& store);
  
  qtensor< dtensor<double> > operator[](int n);
  
  void move(int n);

};

#endif
