
#include "mps.h"

mps::mps(int sites)
{
  
}


int mps::initialize(mpo& in)
{

}


void mps::create_store(vector<tensor>& store)
{

}

qtensor< dtensor<double> > mps::operator[](int n)
{
  return dstate_[n];
}

void mps::move(int n)
{
  
}
