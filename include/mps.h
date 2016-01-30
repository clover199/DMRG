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
#ifdef PBC
  vector< qtensor<T> > edge_;
#endif
  qtensor<double> singular_;
  
public:

  mps(const mpo<T>& in) {
    cut_ = 0;
    state_.resize(in.size()); 
    store_.resize(in.size()); 
#ifdef PBC
    edge_.resize(in.size());
#endif
  }
  
  void print() const {
    cout << "Print all MPS state" << endl;
    for(int i=0;i<state_.size();i++)
    {
      cout << i << ": ";cout<<endl;
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
#ifdef PBC
    cout << "Print all edge tensor" << endl;
    for(int i=0;i<edge_.size();i++)
    {
      cout << i << ": ";
      edge_[i].print();
      cout << endl;
    }
#endif
  }
  
#ifdef PBC
  void add_edge(const qtensor<T>& left, const qtensor<T>& right) {
    edge_[0] = left;
    edge_[edge_.size()-1] = right;
  }

  qtensor<T>& edge(int n) { return edge_[n];}
#endif

  int size() const { return state_.size(); }
  
  qtensor<T>& operator[](int n) { return state_[n]; }
  
  qtensor<T>& operator()(int n) { return store_[n]; }
  
  qtensor<double>& center() { return singular_; }
};

#endif
