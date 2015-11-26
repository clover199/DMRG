#ifndef _MATRIX_PRODUCT_OPERATOR_
#define _MATRIX_PRODUCT_OPERATOR_

#include "global.h"
#include "tensor.h"

class mpo
{
  friend class mps;
private:

  int sites_;
  vector<tensor> A_;
  
public:

  mpo(int sites, tensor ham);


};

#endif
