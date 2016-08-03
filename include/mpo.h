#ifndef _MATRIX_PRODUCT_OPERATOR_
#define _MATRIX_PRODUCT_OPERATOR_

#include "global.h"
#include "qtensor.h"

template <typename T>
class mpo
{
private:

  int sites_;
  
  // Store only the different operators. The operators for sites are
  // O0 O1 O2 ... O0 O1 O2 ... O0 O1 O2 ... ... start from site 0
  vector< qtensor<T> > operators_;
  
  qtensor<T> l_edge_;  // overlap O0
  
  qtensor<T> r_edge_;  // overlap the operator for the last site
  
public:

  mpo(int sites, const qtensor<T>& ham) {
    sites_ = sites;
    operators_.resize(1);
    operators_[0] = ham;
 }
 
  mpo(int sites) {
    sites_ = sites;
    operators_.resize(sites);
  }

  void update(const qtensor<T>& ham) {
    operators_.push_back(ham);
  }

  int size() const { return sites_; }  // return the system size

  void resize(int n) { operators_.resize(n); }

  int num() const { return operators_.size(); }  // return the number of stored MPO

  void edge() {
    l_edge_ = operators_[0].left();
    r_edge_ = operators_[(sites_-1)%operators_.size()].right();
  }

  void edge(const qtensor<T>& left, const qtensor<T>& right) {
    l_edge_ = left;
    r_edge_ = right;
  }

  qtensor<T>& operator[](int n) {
    if(n==0) return l_edge_;
      else if(n==sites_-1) return r_edge_;
    else return operators_[n%operators_.size()];
  }
  
};

#endif
