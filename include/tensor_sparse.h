#ifndef _TENSOR_SPARSE_
#define _TENSOR_SPARSE_

#include "global.h"
#include "tensor_dense.h"

// Attention: the constructor function and functions:
// 'exchange', 'transpose', 'cutoff'
// can only accept at most 8 indexes

// the sparse tensor, T = 'double', 'complex<double>'

template <typename T> 
class stensor
{
  template <typename>  friend class dtensor;
  template <typename>  friend class stensor;

private:
  vector<int> index_;
  vector< vector<int> > cor_;
  vector< T > val_;

  // array to number
  int a2n(const vector<int> in) const;

  // number to array
  vector<int> n2a(int num) const;
	
public:
  // claim space with all elements as zero
  stensor(int i0=1, int i1=0, int i2=0, int i3=0, 
          int i4=0, int i5=0, int i6=0, int i7=0);

  void update(T val,
              int i0, int i1=-1, int i2=-1, int i3=-1, 
              int i4=-1, int i5=-1, int i6=-1, int i7=-1);

  // from dense to sparse
  stensor& convertfrom(const dtensor<T>& in);

  // from sparse to dense
  dtensor<T>& convertto(dtensor<T>& in) const;

  void sort();
  
  int rank() const {return index_.size();}

  int dimension() const;
  int dimension(int i) const;

  void print() const;

  // exchange the two indexes at positions a, b (start from 0)
  // for a two-tensor and a=0, b=1, this is just matrix transpose.
  stensor exchange(int a=0, int b=1) const;

  // conjugate
  stensor conjugate() const;

  // transpose as a matrix with the first num indexes as left index
  // and the rest indexes as right index
  stensor transpose(int num=1) const;
  
  // retain 'cutoff' number of index 'index'
  stensor<T> cut(int index, int cutoff) const;
  
  // result = alpha * A + beta * B
  stensor& plus(const stensor& A, const stensor& B, T alpha=1, T beta=1);

  // create a num-tensor with only diagonal elements non-zero
  stensor& diag(int num=1);
};
#endif
