
#include "global.h"
#include "tensor.h"

tensor::tensor() {dense_ = false; complex_ = false;}

tensor::tensor(qtensor< dtensor<double> >& in) { dd_ = in; dense_ = true; complex_ = false;}

tensor::tensor(qtensor< stensor<double> >& in) { ds_ = in; dense_ = false; complex_ = false;}

tensor::tensor(qtensor< dtensor< complex<double> > >& in)
{ zd_ = in; dense_ = true; complex_ = true;}

tensor::tensor(qtensor< stensor< complex<double> > >& in)
{ zs_ = in; dense_ = false; complex_ = true;}


void tensor::back(qtensor< dtensor<double> >& out)
{
  if( (dense_==true) and (complex_==false) )
    out = dd_;
  else cerr << "Error in tensor::back: no dense-double tensor available\n";
}

void tensor::back(qtensor< stensor<double> >& out)
{
  if( (dense_==false) and (complex_==false) )
    out = ds_;
  else cerr << "Error in tensor::back: no sparse-double tensor available\n";
}

void tensor::back(qtensor< dtensor< complex<double> > >& out)
{
  if( (dense_==true) and (complex_==true) )
    out = zd_;
  else cerr << "Error in tensor::back: no dense-complex tensor available\n";
}

void tensor::back(qtensor< stensor< complex<double> > >& out)
{
  if( (dense_==false) and (complex_==true) )
    out = zs_;
  else cerr << "Error in tensor::back: no sparse-complex tensor available\n";
}

