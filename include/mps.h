#ifndef _MATRIX_PRODUCT_STATE_
#define _MATRIX_PRODUCT_STATE_

#include "global.h"
#include "qtensor.h"
#include "mpo.h"

template <typename T>
class mps
{
private:

  int cut_;

  vector< qtensor<T> > state_;

  vector< qtensor<T> > store_;

  qtensor<double> singular_;
  
public:

  mps(const mpo<T>& in) {
    cut_ = 0;
    state_.resize(in.size()); 
    store_.resize(in.size()); 
  }
  
  void print() const {
    cout << "Print all MPS state" << endl;
    for(int i=0;i<state_.size();i++)
    {
      cout << i << ": ";
      state_[i].print();
      cout << endl;
    }
    cout << "Print all stored tensor" << endl;
    for(int i=0;i<store_.size();i++)
    {
      cout << i << ": ";
      store_[i].print();
      cout << endl;
    }
  }

  int size() const { return state_.size(); }
  
  qtensor<T>& operator[](int n) { return state_[n]; }
  
  qtensor<T>& operator()(int n) { return store_[n]; }
  
  qtensor<double>& center() { return singular_; }
};

#endif
