#ifndef _MATRIX_PRODUCT_STATE_
#define _MATRIX_PRODUCT_STATE_

#include "global.h"
#include "tensor.h"
#include "mpo.h"

class mps
{
private:
  int cut_;
  vector<tensor> state_;
  vector<double> singular_;
  
public:

  // set the number of sites
  mps(int sites);

  // initialize the structure of the MPS according to MPO
  // return the size of the system
  int initialize(mpo& in);
  
  void create_store(vector<tensor>& store);
  
  tensor operator[](int n) { return state_[n]; }
  
  vector<double> middle() { return singular_; }
  
  void move(int n);

};

#endif
