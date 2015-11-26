#ifndef _TENSOR_ALL_TYPES_
#define _TENSOR_ALL_TYPES_

#include "global.h"
#include "tensor_quantum.h"

// the tensor class that includes all types of qtensor

class tensor
{
private:
  bool dense_;
  bool complex_;
  qtensor< dtensor<double> > dd_;
  qtensor< stensor<double> > ds_;
  qtensor< dtensor< complex<double> > > zd_;
  qtensor< stensor< complex<double> > > zs_;
	
public:
  tensor();
  tensor(qtensor< dtensor<double> >& in);
  tensor(qtensor< stensor<double> >& in);
  tensor(qtensor< dtensor< complex<double> > >& in);
  tensor(qtensor< stensor< complex<double> > >& in);
  
  void back(qtensor< dtensor<double> >& out);
  void back(qtensor< stensor<double> >& out);
  void back(qtensor< dtensor< complex<double> > >& out);
  void back(qtensor< stensor< complex<double> > >& out);

  void svd(tensor& U, vector<double>& S, tensor& V);

};
#endif
