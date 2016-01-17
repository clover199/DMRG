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
  
  int size() const { return state_.size(); }
  
  qtensor<T>& operator[](int n) { return state_[n]; }
  
  qtensor<T>& operator()(int n) { return store_[n]; }
  
  qtensor<double>& middle() { return singular_; }
  
  void update_state(int n, qtensor<T>& in) { state_[n] = in; }
  
  void update_store(int n, qtensor<T>& in) { state_[n] = in; }
  
  void update(qtensor<double>& in) { singular_ = in; }

};

#endif
