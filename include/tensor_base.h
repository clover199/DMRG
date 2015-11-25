#ifndef _TENSOR_DENSE_
#define _TENSOR_DENSE_

#include "global.h"

// Attention: the constructor function and functions:
// 'exchange', 'transpose', 'cutoff'
// can only accept at most 8 indexes

// the dense tensor, T = 'double', 'complex<double>'

template <typename T> 
class btensor
{
  template <typename>  friend class btensor;
  template <typename>  friend class stensor;

private:
  vector<int> index_;
  vector< T > val_;
  
  inline T* begin() { return &*val_.begin(); }
	
public:
  // claim space with all elements as zero
  btensor(int i0=1, int i1=0, int i2=0, int i3=0, 
         int i4=0, int i5=0, int i6=0, int i7=0);

  void update(T val,
              int i0, int i1=-1, int i2=-1, int i3=-1, 
              int i4=-1, int i5=-1, int i6=-1, int i7=-1);
  
  int rank() const {return index_.size();}
 
  int dimension() const;
  int dimension(int i) const;

  void print() const;

  // exchange the two indexes at positions a, b (start from 0)
  // for a two-tensor and a=0, b=1, this is just matrix transpose.
  btensor exchange(int a=0, int b=1) const;

  // conjugate
  btensor conjugate() const;

  // transpose as a matrix with the first num indexes as left index
  // and the rest indexes as right index
  btensor transpose(int num=1) const;
  
  // retain 'cut' number of index 'index'
  btensor<T> cut(int index, int cutoff) const;
  
  // result = alpha * A + beta * B
  btensor& plus(const btensor& A, const btensor& B, T alpha=1, T beta=1);
  
  // contract the last/first num indexes for the left/right tensor
  btensor& contract(btensor& A, btensor& B, int num=1);
  
  // multiply the tensor by a tensor (in the form of array) wiht index l,m,r
  // return a tensor also in the form of array
  // i.e. update the ii index of the input tensor
  T* contract(T* in, int ii, int left, int mid, int right);

  // combine the index of tensor A, B
  btensor& combine(btensor& A, btensor& B);
  
  // create a num-tensor with only diagonal elements non-zero
  btensor& diag(int num=1);
  
  // SVD between first num and the rest indexes
  // A = U * S * V
  // returns the number of singular values (length of S)
  int svd(btensor& U, btensor<double>& S, btensor& V, int num=1);
};
#endif
