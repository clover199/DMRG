#ifndef _MATRIX_PRODUCT_OPERATOR_
#define _MATRIX_PRODUCT_OPERATOR_

#include "global.h"
#include "tensor_quantum.h"
#include "tensor.h"

class mpo
{
	
private:

  int sites_;
  vector< qtensor< stensor<double> > > A_;
  
public:

  mpo(int sites, tensor ham);


};

#endif
